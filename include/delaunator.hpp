#pragma once

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace delaunator {

constexpr double EPSILON = std::numeric_limits<double>::epsilon();
constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

//@see https://stackoverflow.com/questions/33333363/built-in-mod-vs-custom-mod-function-improve-the-performance-of-modulus-op/33333636#33333636
inline size_t fast_mod(const size_t i, const size_t c) {
    return i >= c ? i % c : i;
}

template <typename PointType>
double get_x_value(PointType const& pt) {
    return pt.x;
}

template <typename PointType>
double get_y_value(PointType const& pt) {
    return pt.y;
}

struct Point {
    double x;
    double y;
};

inline double circumradius(const double ax,
                           const double ay,
                           const double bx,
                           const double by,
                           const double cx,
                           const double cy) {
    const double dx = bx - ax;
    const double dy = by - ay;
    const double ex = cx - ax;
    const double ey = cy - ay;

    const double bl = dx * dx + dy * dy;
    const double cl = ex * ex + ey * ey;
    const double d = dx * ey - dy * ex;

    const double x = (ey * bl - dy * cl) * 0.5 / d;
    const double y = (dx * cl - ex * bl) * 0.5 / d;

    if ((bl > 0.0 || bl < 0.0) && (cl > 0.0 || cl < 0.0) && (d > 0.0 || d < 0.0)) {
        return x * x + y * y;
    } else {
        return std::numeric_limits<double>::max();
    }
}

struct HullPoint;

struct HullPoint {
    std::size_t index;
    std::size_t triangle_index;
    double x;
    double y;
    HullPoint * next;
    HullPoint * prev;

    HullPoint(std::size_t index_, double x_, double y_) :
            index(index_),
            triangle_index(0),
            x(x_),
            y(y_),
            next(nullptr),
            prev(nullptr) {}
};

inline bool orient(HullPoint * pt_p,
                   HullPoint * pt_q,
                   HullPoint * pt_r) {
    return (pt_q->y - pt_p->y) * (pt_r->x - pt_q->x) - (pt_q->x - pt_p->x) * (pt_r->y - pt_q->y) < 0.0;
}

inline bool in_circle(HullPoint const& a,
                      HullPoint const& b,
                      HullPoint const& c,
                      HullPoint const& p) {
    const double dx = a.x - p.x;
    const double dy = a.y - p.y;
    const double ex = b.x - p.x;
    const double ey = b.y - p.y;
    const double fx = c.x - p.x;
    const double fy = c.y - p.y;

    const double ap = dx * dx + dy * dy;
    const double bp = ex * ex + ey * ey;
    const double cp = fx * fx + fy * fy;

    return (dx * (ey * cp - bp * fy) -
            dy * (ex * cp - bp * fx) +
            ap * (ex * fy - ey * fx)) < 0.0;
}

inline bool check_pts_equal(HullPoint * pt1, HullPoint * pt2) {
    return std::fabs(pt1->x - pt2->x) <= EPSILON &&
           std::fabs(pt1->y - pt2->y) <= EPSILON;
}

struct DelaunatorMap {

    std::vector<HullPoint*> values;
    Point center;

    DelaunatorMap(std::size_t size,
                  Point const& center_):
            values(size, nullptr),
            center(center_) {}

    static double pseudo_angle(const double dx, const double dy) {
        const double p = dx / (std::abs(dx) + std::abs(dy));
        return (dy > 0.0 ? 3.0 - p : 1.0 + p) / 4.0; // [0..1)
    }

    std::size_t hash(const double x, const double y) const {
        const double dx = x - center.x;
        const double dy = y - center.y;
        return fast_mod(
            static_cast<std::size_t>(std::llround(std::floor(pseudo_angle(dx, dy) * static_cast<double>(values.size())))),
            values.size());
    }

    void add(HullPoint * pt) {
        values[hash(pt->x, pt->y)] = pt;
    }
    
    HullPoint * find(const double x, const double y) {
        auto start = values.begin() + static_cast<std::vector<HullPoint*>::difference_type>(hash(x, y));
        
        for (auto iter = start; iter != values.end(); ++iter) {
            if (*iter != nullptr && (*iter)->next != nullptr) {
                return *iter;
            }
        }
        for (auto iter = values.begin(); iter != start; ++iter) {
            if (*iter != nullptr && (*iter)->next != nullptr) {
                return *iter;
            }
        }
        // Should not reach here under normal operation
        return nullptr;
    }
};

struct DelaunatorHull {

    std::vector<HullPoint> points;
    std::vector<HullPoint*> sorted_points;
    std::size_t start_index;

    template <typename PointContainer>
    DelaunatorHull(PointContainer const& coords):
            points(),
            start_index(0) {
        points.reserve(coords.size());
        sorted_points.reserve(coords.size());
        for (std::size_t i = 0; i < coords.size(); ++i) {
            points.push_back(HullPoint(i, get_x_value(coords[i]), get_y_value(coords[i])));
        }
        for (auto & pt : points) {
            sorted_points.push_back(&pt);   
        }
    }

    static double dist(const double ax,
                       const double ay,
                       const double bx,
                       const double by) {
        const double dx = ax - bx;
        const double dy = ay - by;
        return dx * dx + dy * dy;
    }

    Point find_circumcenter(HullPoint * pt_a,
                            HullPoint * pt_b,
                            HullPoint * pt_c) {
        const double dx = pt_b->x - pt_a->x;
        const double dy = pt_b->y - pt_a->y;
        const double ex = pt_c->x - pt_a->x;
        const double ey = pt_c->y - pt_a->y;

        const double bl = dx * dx + dy * dy;
        const double cl = ex * ex + ey * ey;
        const double d = dx * ey - dy * ex;

        const double x = pt_a->x + (ey * bl - dy * cl) * 0.5 / d;
        const double y = pt_a->y + (dx * cl - ex * bl) * 0.5 / d;

        return {x, y};
    }

    void sort(Point const& center) {
        std::sort(sorted_points.begin(), sorted_points.end(), [center](HullPoint * a, HullPoint * b) {
            const double d1 = dist(a->x, a->y, center.x, center.y);
            const double d2 = dist(b->x, b->y, center.x, center.y);
            const double diff1 = d1 - d2;
            const double diff2 = a->x - b->x;
            const double diff3 = a->y - b->y;
            if (diff1 > 0.0 || diff1 < 0.0) {
                return diff1 < 0;
            } else if (diff2 > 0.0 || diff2 < 0.0) {
                return diff2 < 0;
            } else {
                return diff3 < 0;
            }
        });
    }

    Point find_centroid() {
        double max_x = std::numeric_limits<double>::min();
        double max_y = std::numeric_limits<double>::min();
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        for (auto const& pt : points) {
            min_x = std::min(pt.x, min_x);
            min_y = std::min(pt.y, min_y);
            max_x = std::max(pt.x, max_x);
            max_y = std::max(pt.y, max_y);
        }
        return {(static_cast<double>(min_x) + static_cast<double>(max_x)) / 2.0,
                (static_cast<double>(min_y) + static_cast<double>(max_y)) / 2.0 }; 
    }

    HullPoint * find_point_nearest_centroid(Point const& centroid) {

        HullPoint * closest = nullptr;
        double min_dist = std::numeric_limits<double>::max();
        
        // pick a seed point close to the centroid
        for (auto & pt: points) {
            const double d = dist(centroid.x, centroid.y, pt.x, pt.y);
            if (d < min_dist) {
                closest = &pt;
                min_dist = d;
            }
        }
        return closest;
    }

    HullPoint * find_next_closest_point(HullPoint * centroid_pt) {
        HullPoint * closest = nullptr;
        double min_dist = std::numeric_limits<double>::max();
        for (auto & pt : points) {
            HullPoint * pt_ptr = &pt;
            if (centroid_pt == pt_ptr) continue;
            const double d = dist(centroid_pt->x, centroid_pt->y, pt.x, pt.y);
            if (d < min_dist && d > 0.0) {
                closest = pt_ptr;
                min_dist = d;
            }
        }
        return closest;
    }

    HullPoint * find_closest_circumradius_point(HullPoint * pt1,
                                                HullPoint * pt2) {
        HullPoint * closest = nullptr;
        double min_radius = std::numeric_limits<double>::max();
        for (auto & pt : points) {
            const double r = circumradius(pt1->x, pt1->y, pt2->x, pt2->y, pt.x, pt.y);
            if (r < min_radius) {
                closest = &pt;
                min_radius = r;
            }
        }
        if (!(min_radius < std::numeric_limits<double>::max())) {
            throw std::runtime_error("not triangulation");
        }
        return closest;
    }

};

struct DelaunatorResult {

    std::vector<std::size_t> triangles;
    std::vector<std::size_t> halfedges;
    std::size_t hull_start;
    bool success;

    DelaunatorResult(): 
            triangles(), 
            halfedges(),
            hull_start(0), 
            success(false) {}
    
    DelaunatorResult(std::size_t size): 
            triangles(),
            halfedges(),
            hull_start(0),
            success(false) {
        std::size_t max_triangles = size < 3 ? 1 : 2 * size - 5;
        triangles.reserve(max_triangles * 3);
        halfedges.reserve(max_triangles * 3);
    }

    void link(std::size_t a, std::size_t b) {
        std::size_t s = halfedges.size();
        if (a == s) {
            halfedges.push_back(b);
        } else if (a < s) {
            halfedges[a] = b;
        } else {
            throw std::runtime_error("Cannot link edge");
        }
        if (b != INVALID_INDEX) {
            std::size_t s2 = halfedges.size();
            if (b == s2) {
                halfedges.push_back(a);
            } else if (b < s2) {
                halfedges[b] = a;
            } else {
                throw std::runtime_error("Cannot link edge");
            }
        }
    }

    std::size_t add_triangle(std::size_t i0,
                             std::size_t i1,
                             std::size_t i2,
                             std::size_t a,
                             std::size_t b,
                             std::size_t c) {
        std::size_t t = triangles.size();
        triangles.push_back(i0);
        triangles.push_back(i1);
        triangles.push_back(i2);
        link(t, a);
        link(t + 1, b);
        link(t + 2, c);
        return t;
    }

    std::size_t legalize(DelaunatorHull & hull, std::size_t a) {
        const std::size_t b = halfedges[a];

        /* if the pair of triangles doesn't satisfy the Delaunay condition
        * (p1 is inside the circumcircle of [p0, pl, pr]), flip them,
        * then do the same check/flip recursively for the new pair of triangles
        *
        *           pl                    pl
        *          /||\                  /  \
        *       al/ || \bl            al/    \a
        *        /  ||  \              /      \
        *       /  a||b  \    flip    /___ar___\
        *     p0\   ||   /p1   =>   p0\---bl---/p1
        *        \  ||  /              \      /
        *       ar\ || /br             b\    /br
        *          \||/                  \  /
        *           pr                    pr
        */
        const std::size_t a0 = a - a % 3;
        const std::size_t b0 = b - b % 3;

        const std::size_t al = a0 + (a + 1) % 3;
        const std::size_t ar = a0 + (a + 2) % 3;
        const std::size_t bl = b0 + (b + 2) % 3;

        const std::size_t p0 = triangles[ar];
        const std::size_t pr = triangles[a];
        const std::size_t pl = triangles[al];
        const std::size_t p1 = triangles[bl];

        if (b == INVALID_INDEX) {
            return ar;
        }

        bool illegal = in_circle(hull.points[p0],
                                 hull.points[pr],
                                 hull.points[pl],
                                 hull.points[p1]);
        
        if (illegal) {
            triangles[a] = p1;
            triangles[b] = p0;

            std::size_t hbl = halfedges[bl];

            // edge swapped on the other side of the hull (rare); fix the halfedge reference
            if (hbl == INVALID_INDEX) {
                HullPoint * start = &(hull.points[hull_start]);
                HullPoint * e = start;
                do {
                    if (e->triangle_index == bl) {
                        e->triangle_index = a;
                        break;
                    }
                    e = e->next;
                } while (e != start);
            }
            link(a, hbl);
            link(b, halfedges[ar]);
            link(ar, bl);

            std::size_t br = b0 + (b + 1) % 3;

            legalize(hull, a);
            return legalize(hull, br);
        }
        return ar;
    }
};

template <typename PointType>
DelaunatorResult delaunator(std::vector<PointType> const& coords) {
 
    std::size_t point_size = coords.size();
    
    if (point_size < 3) {
        // Not enough points to triangulate
        return DelaunatorResult();
    }
    
    DelaunatorResult result(point_size);

    DelaunatorHull hull(coords);

    Point centroid = hull.find_centroid();
    HullPoint * pt0 = hull.find_point_nearest_centroid(centroid);
    HullPoint * pt1 = hull.find_next_closest_point(pt0);
    HullPoint * pt2 = hull.find_closest_circumradius_point(pt0, pt1);
    
    if (orient(pt0, pt1, pt2)) {
        std::swap(pt1, pt2);
    }

    Point center = hull.find_circumcenter(pt0, pt1, pt2);

    hull.sort(center);
    
    std::size_t hash_size = static_cast<std::size_t>(std::llround(std::ceil(std::sqrt(point_size))));
    DelaunatorMap hash(hash_size, center);
    
    // Setup linking of first points now
    pt0->next = pt2->prev = pt1;
    pt1->next = pt0->prev = pt2;
    pt2->next = pt1->prev = pt0;

    pt0->triangle_index = 0;
    pt1->triangle_index = 1;
    pt2->triangle_index = 2;

    result.hull_start = pt0->index;
    
    hash.add(pt0);
    hash.add(pt1);
    hash.add(pt2);

    result.add_triangle(pt0->index, pt1->index, pt2->index, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
    HullPoint * last_pt = nullptr;

    for (auto pt : hull.sorted_points) {

        // skip near-duplicate points
        if (last_pt && check_pts_equal(pt, last_pt)) continue;
        last_pt = pt;

        // skip seed triangle points
        if (check_pts_equal(pt, pt0) ||
            check_pts_equal(pt, pt1) ||
            check_pts_equal(pt, pt2)) continue;

        // find a visible edge on the convex hull using edge hash
        auto start = hash.find(pt->x, pt->y);

        if (!start) {
            throw std::runtime_error("No values in hash");
        }
        start = start->prev;
        
        HullPoint * e = start;
        HullPoint * q = e->next;

        while (!orient(pt, e, q)) {
            e = q;
            if (e == start) {
                e = nullptr;
                break;
            }
            q = e->next;
        }

        if (!e) continue; // likely a near-duplicate point; skip it

        // add the first triangle from the point
        std::size_t t = result.add_triangle(e->index,
                                            pt->index,
                                            e->next->index,
                                            INVALID_INDEX,
                                            INVALID_INDEX,
                                            e->triangle_index);
        
        pt->triangle_index = result.legalize(hull, t + 2);
        e->triangle_index = t;

        // walk forward through the hull, adding more triangles and flipping recursively

        HullPoint * n = e->next;
        q = n->next;

        while (orient(pt, n, q)) {
            t = result.add_triangle(n->index, pt->index, q->index, pt->triangle_index, INVALID_INDEX, n->triangle_index);
            pt->triangle_index = result.legalize(hull, t + 2);
            n->next = nullptr; // mark as removed
            n = q;
            q = n->next;
        }

        // walk backward from the other side, adding more triangles and flipping
        if (e == start) {
            q = e->prev;
            while (orient(pt, q, e)) {
                t = result.add_triangle(q->index, pt->index, e->index, INVALID_INDEX, e->triangle_index, q->triangle_index);
                result.legalize(hull, t + 2);
                q->triangle_index = t;
                e->next = nullptr; // mark as removed
                e = q;
                q = e->prev;
            }
        }

        // update the hull indices
        pt->prev = e;
        result.hull_start = e->index;
        n->prev = pt;
        e->next = pt;
        pt->next = n;

        hash.add(pt);
        hash.add(e);
    }
    result.success = true;
    return result;
}
/*
// Kahan and Babuska summation, Neumaier variant; accumulates less FP error
inline double sum(const std::vector<double>& x) {
    double sum = x[0];
    double err = 0.0;

    for (size_t i = 1; i < x.size(); i++) {
        const double k = x[i];
        const double m = sum + k;
        err += std::fabs(sum) >= std::fabs(k) ? sum - m + k : k - m + sum;
        sum = m;
    }
    return sum + err;
}

double Delaunator::get_hull_area() {
    std::vector<double> hull_area;
    std::size_t e = hull_start;
    do {
        hull_area.push_back((coords[2 * e] - coords[2 * hull_prev[e]]) * (coords[2 * e + 1] + coords[2 * hull_prev[e] + 1]));
        e = hull_next[e];
    } while (e != hull_start);
    return sum(hull_area);
}
*/
} //namespace delaunator
