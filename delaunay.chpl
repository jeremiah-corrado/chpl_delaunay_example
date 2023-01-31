module test {
    use Delaunay;
    use Subprocess;

    config const numPoints =  6;

    proc main() {
        // get a point cloud
        const xyPoints = randomXY(numPoints);
        // const xyPoints = try! readXY("./input/points_in.txt");

        for (p, i) in zip(xyPoints, xyPoints.domain) {
            writef("%i \t (%1.5r, %1.5r)\n", i, p[0], p[1]);
        }

        // create triangulation
        const tris = computeDelaunayTri(xyPoints);

        // create image
        try! writeTriData("./results/", xyPoints, tris);
        try! Subprocess.spawn(["python3", "plotTri.py", xyPoints.size:string]);
    }

    module Delaunay {
        private use List;
        private use IO;
        private use IO.FormattedIO;
        private use Sort;
        import super.V2D.v2D;

        proc computeDelaunayTri(coordinates: [?d] 2*real) : [] tri
            where d.rank == 1
        {
            // get a sorted array of points from the input points
            var points = sortedPoints(coordinates);

            // define a starting simplex
            var s = new simplex(points[0]!, points[2]!, points[1]!),
                simplices = new list([s,]);

            /*
                cwt := currentWorkingTriangle (end of the simplex list)
                while there are unclaimed points in 'points'
                    pick the edge on 'cwt' closest to the cloud's un-triangulated center of mass
                    with those two points, create a new simplex with any point that meets the Delaunay condition

            */
            var numClaimed = 3;
            label triangulating while numClaimed < d.size {
                const com = openCenterOfMass(points),
                      workingEdge = simplices.last().closestEdgeTo(com);

                writef("nc: %i \t com: (%1.5r %1.5r) \n", numClaimed, com.x, com.y);

                // writeln("nc: ", numClaimed, " com: (", com, " we: --[", workingEdge, "]--");

                for p in points {
                    if !simplices.last().inCircle(p!) {
                        writeln("\t adding: ", workingEdge[0].id, ", ", workingEdge[1].id, ", ", p!.id);

                        simplices.append(new simplex(
                            workingEdge[0],
                            workingEdge[1],
                            p!
                        ));
                        numClaimed += 1;
                        continue triangulating;
                    }
                }

                writeln(simplices);
                halt("found nothing...");
            }

            const tris = [s in simplices] tri.fromSimplex(s);

            // delete the temporary Points
            for p in points do delete p;

            return tris;
        }

        proc sortedPoints(coordinates: [?d] 2*real) : [d] unmanaged Point? {
            var allPoints : [d] unmanaged Point? =
                [i in d] new unmanaged Point(coordinates[i][0], coordinates[i][1], i);

            // define comparator
            record magComp { }; proc magComp.key(a): real { return a!.mag(); }
            sort(allPoints, comparator = new magComp());

            return allPoints;
        }

        record tri {
            var ids: 3*int;

            proc init(a: int, b: int, c: int) {
                this.ids = (a, b, c);
            }

            proc type fromSimplex(s: simplex) {
                return new tri(
                    s.points[0].id,
                    s.points[1].id,
                    s.points[2].id
                );
            }

            iter edges() : 2*int {
                yield (ids[0], ids[1]);
                yield (ids[1], ids[2]);
                yield (ids[2], ids[0]);
            }
        }

        class Point {
            var vec: v2D;
            var id: int = -1;
            var numEdges: int;
            var onHull: bool;

            proc init() {
                this.vec = new v2D(max(real), max(real));
                this.id = -1;
                this.numEdges = 0;
                this.onHull = false;
            }

            proc init(x: real, y: real, id: int) {
                this.vec = new v2D(x, y);
                this.id = id;
                this.numEdges = 0;
                this.onHull = false;
            }

            inline proc dist(other: Point): real {
                return this.vec.distTo(other.vec);
            }

            inline proc mag(): real {
                return this.vec.mag;
            }

            inline proc hasEdge(): bool {
                return this.numEdges > 0;
            }

            inline proc matComponent(const ref p: Point): 3*real {
                return this.vec.matComponent(p.vec);
            }
        }

        record simplex {
            var points: 3*unmanaged Point;

            proc init(in p0: Point, in p1: Point, in p2: Point) {
                // put the points in a clockwise ordering by:
                const center = p0.vec + p1.vec + p2.vec;

                const det = (p0.vec.x - center.x) * (p1.vec.x - center.x) -
                            (p0.vec.y - center.y) * (p1.vec.y - center.y);

                if det < 0 {
                    this.points = (p0, p1, p2);
                } else {
                    this.points = (p0, p2, p1);
                }

                p0.numEdges += 2;
                p1.numEdges += 2;
                p2.numEdges += 2;
            }

            proc this(i: int): borrowed Point {
                return this.points(i).borrow();
            }

            iter edges(): 2*borrowed Point {
                yield (this[0], this[1]);
                yield (this[1], this[2]);
                yield (this[2], this[0]);
            }

            proc inCircle(p: Point): bool {
                if p.id == this[0].id || p.id == this[1].id || p.id == this[2].id then return true;

                // https://en.wikipedia.org/wiki/Delaunay_triangulation#Algorithms
                // (assumes clockwise ordering of points)
                const m = [
                    this[0].matComponent(p),
                    this[1].matComponent(p),
                    this[2].matComponent(p),
                ];

                const det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

                return det > 0;
            }

            proc closestEdgeTo(v: v2D): 2*Point {
                var minDist = max(real), minIdx = -1;
                for ((p1, p2), i) in zip(this.edges(), 0..) {
                    const dist = v.distToLine(p1.vec, p2.vec);

                    if dist < minDist {
                        minIdx = i;
                        minDist = dist;
                    }
                }

                select minIdx {
                    when 0 do return (this.points[0], this.points[1]);
                    when 1 do return (this.points[1], this.points[2]);
                    otherwise do return (this.points[2], this.points[0]);
                }
            }
        }

        proc openCenterOfMass(points: [?d] unmanaged Point?): v2D {
            var sum = new v2D(),
                count = 0;

            forall p in points with (+ reduce sum, + reduce count) {
                if !p!.hasEdge() {
                    sum += p!.vec;
                    count += 1;
                }
            }

            return sum / count;
        }
    }

    module V2D {
        record v2D {
            var x: real;
            var y: real;

            proc init() {
                this.x = 0.0;
                this.y = 0.0;
            }

            proc init(x: real, y: real) {
                this.x = x;
                this.y = y;
            }

            operator +=(const ref other: v2D): v2D {
                this.x += other.x;
                this.y += other.y;
            }

            operator +(const ref a: v2D, const ref b: v2D): v2D {
                return new v2D(a.x + b.x, a.y + b.y);
            }

            operator -=(const ref other: v2D): v2D {
                this.x -= other.x;
                this.y -= other.y;
            }

            operator -(const ref a: v2D, const ref b: v2D): v2D {
                return new v2D(a.x - b.x, a.y - b.y);
            }

            operator /(const ref v: v2D, s: real): v2D {
                return new v2D(v.x / s, v.y / s);
            }

            inline proc mag: real {
                return sqrt(this.x**2 + this.y**2);
            }

            inline proc distTo(const ref other: v2D): real {
                return sqrt((this.x - other.x)**2 + (this.y - other.y)**2);
            }

            inline proc matComponent(const ref other: v2D): 3*real {
                return (this.x - other.x, this.y - other.y, (this.x**2 - other.x**2) + (this.y**2 - other.y**2));
            }

            proc distToLine(const ref a: v2D, const ref b: v2D): real {
                return abs(
                    (b.x-a.x)*(a.y-this.y) -
                    (a.x-this.x)*(b.y-a.y)
                ) / a.distTo(b);
            }
        }
    }

    proc readXY(path) throws {
        use IO;

        var pointLines = openreader(path).readAll(string);
        const nPoints = pointLines.count("\n");

        var points : [0..<nPoints] 2*real;
        for (p, i) in zip(pointLines.split("\n", ignoreEmpty=true), 0..) {
            const (x, _, y) = p.partition(" ");
            points[i] = (x:real, y:real);
        }
        return points;
    }

    proc randomXY(n: int) {
        use Random;

        var d = {0..<n},
            xRand : [d] real,
            yRand : [d] real;
        fillRandom(xRand);
        fillRandom(yRand);

        return [i in d] (xRand[i], yRand[i]);
    }

    proc writeTriData(path, points, tris) throws {
        use IO;
        use Set;

        var pw = openwriter(path + "points.txt");
        for p in points do pw.writef("%1.8r %1.8r\n", (...p));

        var ew = openwriter(path + "edges.txt"),
            edgeSet = new set(2*int);

        for t in tris do edgeSet.add(t.edges());
        for e in edgeSet do ew.writeln(e[0], " ", e[1]);
    }
}
