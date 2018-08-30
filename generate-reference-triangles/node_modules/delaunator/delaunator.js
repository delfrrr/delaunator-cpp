(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? module.exports = factory() :
    typeof define === 'function' && define.amd ? define(factory) :
    (global.Delaunator = factory());
}(this, (function () { 'use strict';

    var Delaunator = function Delaunator(coords) {
        var this$1 = this;

        var minX = Infinity;
        var minY = Infinity;
        var maxX = -Infinity;
        var maxY = -Infinity;

        var n = coords.length >> 1;
        var ids = this.ids = new Uint32Array(n);

        if (n > 0 && typeof coords[0] !== 'number') { throw new Error('Expected coords to contain numbers.'); }

        this.coords = coords;

        for (var i = 0; i < n; i++) {
            var x = coords[2 * i];
            var y = coords[2 * i + 1];
            if (x < minX) { minX = x; }
            if (y < minY) { minY = y; }
            if (x > maxX) { maxX = x; }
            if (y > maxY) { maxY = y; }
            ids[i] = i;
        }

        var cx = (minX + maxX) / 2;
        var cy = (minY + maxY) / 2;

        var minDist = Infinity;
        var i0, i1, i2;

        // pick a seed point close to the centroid
        for (var i$1 = 0; i$1 < n; i$1++) {
            var d = dist(cx, cy, coords[2 * i$1], coords[2 * i$1 + 1]);
            if (d < minDist) {
                i0 = i$1;
                minDist = d;
            }
        }

        minDist = Infinity;

        // find the point closest to the seed
        for (var i$2 = 0; i$2 < n; i$2++) {
            if (i$2 === i0) { continue; }
            var d$1 = dist(coords[2 * i0], coords[2 * i0 + 1], coords[2 * i$2], coords[2 * i$2 + 1]);
            if (d$1 < minDist && d$1 > 0) {
                i1 = i$2;
                minDist = d$1;
            }
        }

        var minRadius = Infinity;

        // find the third point which forms the smallest circumcircle with the first two
        for (var i$3 = 0; i$3 < n; i$3++) {
            if (i$3 === i0 || i$3 === i1) { continue; }

            var r = circumradius(
                coords[2 * i0], coords[2 * i0 + 1],
                coords[2 * i1], coords[2 * i1 + 1],
                coords[2 * i$3], coords[2 * i$3 + 1]);

            if (r < minRadius) {
                i2 = i$3;
                minRadius = r;
            }
        }

        if (minRadius === Infinity) {
            throw new Error('No Delaunay triangulation exists for this input.');
        }

        // swap the order of the seed points for counter-clockwise orientation
        if (area(coords[2 * i0], coords[2 * i0 + 1],
            coords[2 * i1], coords[2 * i1 + 1],
            coords[2 * i2], coords[2 * i2 + 1]) < 0) {

            var tmp = i1;
            i1 = i2;
            i2 = tmp;
        }

        var i0x = coords[2 * i0];
        var i0y = coords[2 * i0 + 1];
        var i1x = coords[2 * i1];
        var i1y = coords[2 * i1 + 1];
        var i2x = coords[2 * i2];
        var i2y = coords[2 * i2 + 1];

        var center = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);
        this._cx = center.x;
        this._cy = center.y;

        // sort the points by distance from the seed triangle circumcenter
        quicksort(ids, coords, 0, ids.length - 1, center.x, center.y);

        // initialize a hash table for storing edges of the advancing convex hull
        this._hashSize = Math.ceil(Math.sqrt(n));
        this._hash = [];
        for (var i$4 = 0; i$4 < this._hashSize; i$4++) { this$1._hash[i$4] = null; }

        // initialize a circular doubly-linked list that will hold an advancing convex hull
        var e = this.hull = insertNode(coords, i0);
        this._hashEdge(e);
        e.t = 0;
        e = insertNode(coords, i1, e);
        this._hashEdge(e);
        e.t = 1;
        e = insertNode(coords, i2, e);
        this._hashEdge(e);
        e.t = 2;

        var maxTriangles = 2 * n - 5;
        var triangles = this.triangles = new Uint32Array(maxTriangles * 3);
        var halfedges = this.halfedges = new Int32Array(maxTriangles * 3);

        this.trianglesLen = 0;

        this._addTriangle(i0, i1, i2, -1, -1, -1);

        for (var k = 0, xp = (void 0), yp = (void 0); k < ids.length; k++) {
            var i$5 = ids[k];
            var x$1 = coords[2 * i$5];
            var y$1 = coords[2 * i$5 + 1];

            // skip duplicate points
            if (x$1 === xp && y$1 === yp) { continue; }
            xp = x$1;
            yp = y$1;

            // skip seed triangle points
            if ((x$1 === i0x && y$1 === i0y) ||
                (x$1 === i1x && y$1 === i1y) ||
                (x$1 === i2x && y$1 === i2y)) { continue; }

            // find a visible edge on the convex hull using edge hash
            var startKey = this$1._hashKey(x$1, y$1);
            var key = startKey;
            var start = (void 0);
            do {
                start = this$1._hash[key];
                key = (key + 1) % this$1._hashSize;
            } while ((!start || start.removed) && key !== startKey);

            e = start;
            while (area(x$1, y$1, e.x, e.y, e.next.x, e.next.y) >= 0) {
                e = e.next;
                if (e === start) {
                    throw new Error('Something is wrong with the input points.');
                }
            }

            var walkBack = e === start;

            // add the first triangle from the point
            var t = this$1._addTriangle(e.i, i$5, e.next.i, -1, -1, e.t);

            e.t = t; // keep track of boundary triangles on the hull
            e = insertNode(coords, i$5, e);

            // recursively flip triangles from the point until they satisfy the Delaunay condition
            e.t = this$1._legalize(t + 2);

            // walk forward through the hull, adding more triangles and flipping recursively
            var q = e.next;
            while (area(x$1, y$1, q.x, q.y, q.next.x, q.next.y) < 0) {
                t = this$1._addTriangle(q.i, i$5, q.next.i, q.prev.t, -1, q.t);
                q.prev.t = this$1._legalize(t + 2);
                this$1.hull = removeNode(q);
                q = q.next;
            }

            if (walkBack) {
                // walk backward from the other side, adding more triangles and flipping
                q = e.prev;
                while (area(x$1, y$1, q.prev.x, q.prev.y, q.x, q.y) < 0) {
                    t = this$1._addTriangle(q.prev.i, i$5, q.i, -1, q.t, q.prev.t);
                    this$1._legalize(t + 2);
                    q.prev.t = t;
                    this$1.hull = removeNode(q);
                    q = q.prev;
                }
            }

            // save the two new edges in the hash table
            this$1._hashEdge(e);
            this$1._hashEdge(e.prev);
        }

        // trim typed triangle mesh arrays
        this.triangles = triangles.subarray(0, this.trianglesLen);
        this.halfedges = halfedges.subarray(0, this.trianglesLen);
    };

    Delaunator.from = function from (points, getX, getY) {
        if (!getX) { getX = defaultGetX; }
        if (!getY) { getY = defaultGetY; }

        var n = points.length;
        var coords = new Float64Array(n * 2);

        for (var i = 0; i < n; i++) {
            var p = points[i];
            coords[2 * i] = getX(p);
            coords[2 * i + 1] = getY(p);
        }

        return new Delaunator(coords);
    };

    Delaunator.prototype._hashEdge = function _hashEdge (e) {
        this._hash[this._hashKey(e.x, e.y)] = e;
    };

    Delaunator.prototype._hashKey = function _hashKey (x, y) {
        var dx = x - this._cx;
        var dy = y - this._cy;
        // use pseudo-angle: a measure that monotonically increases
        // with real angle, but doesn't require expensive trigonometry
        var p = 1 - dx / (Math.abs(dx) + Math.abs(dy));
        return Math.floor((2 + (dy < 0 ? -p : p)) / 4 * this._hashSize);
    };

    Delaunator.prototype._legalize = function _legalize (a) {
        var ref = this;
            var triangles = ref.triangles;
            var coords = ref.coords;
            var halfedges = ref.halfedges;

        var b = halfedges[a];

        /* if the pair of triangles doesn't satisfy the Delaunay condition
         * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
         * then do the same check/flip recursively for the new pair of triangles
         *
         *       pl                pl
         *      /||\              /  \
         *   al/ || \bl        al/\a
         *    /  ||  \          /  \
         *   /  a||b  \flip/___ar___\
         * p0\   ||   /p1   =>   p0\---bl---/p1
         *    \  ||  /          \  /
         *   ar\ || /br         b\/br
         *      \||/              \  /
         *       pr                pr
         */
        var a0 = a - a % 3;
        var b0 = b - b % 3;

        var al = a0 + (a + 1) % 3;
        var ar = a0 + (a + 2) % 3;
        var bl = b0 + (b + 2) % 3;

        var p0 = triangles[ar];
        var pr = triangles[a];
        var pl = triangles[al];
        var p1 = triangles[bl];

        var illegal = inCircle(
            coords[2 * p0], coords[2 * p0 + 1],
            coords[2 * pr], coords[2 * pr + 1],
            coords[2 * pl], coords[2 * pl + 1],
            coords[2 * p1], coords[2 * p1 + 1]);

        if (illegal) {
            triangles[a] = p1;
            triangles[b] = p0;

            var hbl = halfedges[bl];

            // edge swapped on the other side of the hull (rare); fix the halfedge reference
            if (hbl === -1) {
                var e = this.hull;
                do {
                    if (e.t === bl) {
                        e.t = a;
                        break;
                    }
                    e = e.next;
                } while (e !== this.hull);
            }
            this._link(a, hbl);
            this._link(b, halfedges[ar]);
            this._link(ar, bl);

            var br = b0 + (b + 1) % 3;

            this._legalize(a);
            return this._legalize(br);
        }

        return ar;
    };

    Delaunator.prototype._link = function _link (a, b) {
        this.halfedges[a] = b;
        if (b !== -1) { this.halfedges[b] = a; }
    };

    // add a new triangle given vertex indices and adjacent half-edge ids
    Delaunator.prototype._addTriangle = function _addTriangle (i0, i1, i2, a, b, c) {
        var t = this.trianglesLen;

        this.triangles[t] = i0;
        this.triangles[t + 1] = i1;
        this.triangles[t + 2] = i2;

        this._link(t, a);
        this._link(t + 1, b);
        this._link(t + 2, c);

        this.trianglesLen += 3;

        return t;
    };

    function dist(ax, ay, bx, by) {
        var dx = ax - bx;
        var dy = ay - by;
        return dx * dx + dy * dy;
    }

    function area(px, py, qx, qy, rx, ry) {
        return (qy - py) * (rx - qx) - (qx - px) * (ry - qy);
    }

    function inCircle(ax, ay, bx, by, cx, cy, px, py) {
        var dx = ax - px;
        var dy = ay - py;
        var ex = bx - px;
        var ey = by - py;
        var fx = cx - px;
        var fy = cy - py;

        var ap = dx * dx + dy * dy;
        var bp = ex * ex + ey * ey;
        var cp = fx * fx + fy * fy;

        return dx * (ey * cp - bp * fy) -
               dy * (ex * cp - bp * fx) +
               ap * (ex * fy - ey * fx) < 0;
    }

    function circumradius(ax, ay, bx, by, cx, cy) {
        var dx = bx - ax;
        var dy = by - ay;
        var ex = cx - ax;
        var ey = cy - ay;

        var bl = dx * dx + dy * dy;
        var cl = ex * ex + ey * ey;
        var d = dx * ey - dy * ex;

        var x = (ey * bl - dy * cl) * 0.5 / d;
        var y = (dx * cl - ex * bl) * 0.5 / d;

        return bl && cl && d && (x * x + y * y) || Infinity;
    }

    function circumcenter(ax, ay, bx, by, cx, cy) {
        var dx = bx - ax;
        var dy = by - ay;
        var ex = cx - ax;
        var ey = cy - ay;

        var bl = dx * dx + dy * dy;
        var cl = ex * ex + ey * ey;
        var d = dx * ey - dy * ex;

        var x = ax + (ey * bl - dy * cl) * 0.5 / d;
        var y = ay + (dx * cl - ex * bl) * 0.5 / d;

        return {x: x, y: y};
    }

    // create a new node in a doubly linked list
    function insertNode(coords, i, prev) {
        var node = {
            i: i,
            x: coords[2 * i],
            y: coords[2 * i + 1],
            t: 0,
            prev: null,
            next: null,
            removed: false
        };

        if (!prev) {
            node.prev = node;
            node.next = node;

        } else {
            node.next = prev.next;
            node.prev = prev;
            prev.next.prev = node;
            prev.next = node;
        }
        return node;
    }

    function removeNode(node) {
        node.prev.next = node.next;
        node.next.prev = node.prev;
        node.removed = true;
        return node.prev;
    }

    function quicksort(ids, coords, left, right, cx, cy) {
        var i, j, temp;

        if (right - left <= 20) {
            for (i = left + 1; i <= right; i++) {
                temp = ids[i];
                j = i - 1;
                while (j >= left && compare(coords, ids[j], temp, cx, cy) > 0) { ids[j + 1] = ids[j--]; }
                ids[j + 1] = temp;
            }
        } else {
            var median = (left + right) >> 1;
            i = left + 1;
            j = right;
            swap(ids, median, i);
            if (compare(coords, ids[left], ids[right], cx, cy) > 0) { swap(ids, left, right); }
            if (compare(coords, ids[i], ids[right], cx, cy) > 0) { swap(ids, i, right); }
            if (compare(coords, ids[left], ids[i], cx, cy) > 0) { swap(ids, left, i); }

            temp = ids[i];
            while (true) {
                do { i++; } while (compare(coords, ids[i], temp, cx, cy) < 0);
                do { j--; } while (compare(coords, ids[j], temp, cx, cy) > 0);
                if (j < i) { break; }
                swap(ids, i, j);
            }
            ids[left + 1] = ids[j];
            ids[j] = temp;

            if (right - i + 1 >= j - left) {
                quicksort(ids, coords, i, right, cx, cy);
                quicksort(ids, coords, left, j - 1, cx, cy);
            } else {
                quicksort(ids, coords, left, j - 1, cx, cy);
                quicksort(ids, coords, i, right, cx, cy);
            }
        }
    }

    function compare(coords, i, j, cx, cy) {
        var d1 = dist(coords[2 * i], coords[2 * i + 1], cx, cy);
        var d2 = dist(coords[2 * j], coords[2 * j + 1], cx, cy);
        return (d1 - d2) || (coords[2 * i] - coords[2 * j]) || (coords[2 * i + 1] - coords[2 * j + 1]);
    }

    function swap(arr, i, j) {
        var tmp = arr[i];
        arr[i] = arr[j];
        arr[j] = tmp;
    }

    function defaultGetX(p) {
        return p[0];
    }
    function defaultGetY(p) {
        return p[1];
    }

    return Delaunator;

})));
