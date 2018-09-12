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

inline double dist(
    const double ax,
    const double ay,
    const double bx,
    const double by) {
    const double dx = ax - bx;
    const double dy = ay - by;
    return dx * dx + dy * dy;
}

inline double circumradius(
    const double ax,
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

inline double area(
    const double px,
    const double py,
    const double qx,
    const double qy,
    const double rx,
    const double ry) {
    return (qy - py) * (rx - qx) - (qx - px) * (ry - qy);
}

inline std::pair<double, double> circumcenter(
    const double ax,
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

    const double x = ax + (ey * bl - dy * cl) * 0.5 / d;
    const double y = ay + (dx * cl - ex * bl) * 0.5 / d;

    return std::make_pair(x, y);
}

inline double compare(
    std::vector<double> const& coords,
    std::size_t i,
    std::size_t j,
    double cx,
    double cy) {
    const double d1 = dist(coords[2 * i], coords[2 * i + 1], cx, cy);
    const double d2 = dist(coords[2 * j], coords[2 * j + 1], cx, cy);
    const double diff1 = d1 - d2;
    const double diff2 = coords[2 * i] - coords[2 * j];
    const double diff3 = coords[2 * i + 1] - coords[2 * j + 1];

    if (diff1 > 0.0 || diff1 < 0.0) {
        return diff1;
    } else if (diff2 > 0.0 || diff2 < 0.0) {
        return diff2;
    } else {
        return diff3;
    }
}

struct sort_to_center {

    std::vector<double> const& coords;
    double cx;
    double cy;

    bool operator()(std::size_t i, std::size_t j) {
        return compare(coords, i, j, cx, cy) < 0;
    }
};

inline bool in_circle(
    double ax,
    double ay,
    double bx,
    double by,
    double cx,
    double cy,
    double px,
    double py) {
    const double dx = ax - px;
    const double dy = ay - py;
    const double ex = bx - px;
    const double ey = by - py;
    const double fx = cx - px;
    const double fy = cy - py;

    const double ap = dx * dx + dy * dy;
    const double bp = ex * ex + ey * ey;
    const double cp = fx * fx + fy * fy;

    return (dx * (ey * cp - bp * fy) -
            dy * (ex * cp - bp * fx) +
            ap * (ex * fy - ey * fx)) < 0.0;
}

inline bool check_pts_equal(double x1, double y1, double x2, double y2) {
    return std::fabs(x1 - x2) < std::numeric_limits<double>::epsilon() &&
           std::fabs(y1 - y2) < std::numeric_limits<double>::epsilon();
}

constexpr std::size_t INVALID_INDEX = std::numeric_limits<std::size_t>::max();

struct DelaunatorPoint {
    std::size_t i;
    double x;
    double y;
    std::size_t t;
    std::size_t prev;
    std::size_t next;
    bool removed;
};

class Delaunator {

public:
    std::vector<double> const& coords;
    std::vector<std::size_t> triangles;
    std::vector<std::size_t> halfedges;

    Delaunator(std::vector<double> const& in_coords);

private:
    std::vector<std::size_t> m_hash;
    std::vector<DelaunatorPoint> m_hull;
    double m_center_x;
    double m_center_y;
    std::size_t m_hash_size;

    std::size_t remove_node(std::size_t node);
    std::size_t legalize(std::size_t a);
    std::size_t insert_node(std::size_t i);
    std::size_t insert_node(std::size_t i, std::size_t prev);
    std::size_t hash_key(double x, double y);
    void hash_edge(std::size_t e);
    std::size_t add_triangle(
        std::size_t i0,
        std::size_t i1,
        std::size_t i2,
        std::size_t a,
        std::size_t b,
        std::size_t c);
    void link(std::size_t a, std::size_t b);
};

Delaunator::Delaunator(std::vector<double> const& in_coords)
    : coords(in_coords),
      triangles(),
      halfedges(),
      m_hash(),
      m_hull(),
      m_center_x(),
      m_center_y(),
      m_hash_size() {
    std::size_t n = coords.size() >> 1;

    double max_x = std::numeric_limits<double>::min();
    double max_y = std::numeric_limits<double>::min();
    double min_x = std::numeric_limits<double>::max();
    double min_y = std::numeric_limits<double>::max();
    std::vector<std::size_t> ids;
    ids.reserve(n);

    for (std::size_t i = 0; i < n; i++) {
        const double x = coords[2 * i];
        const double y = coords[2 * i + 1];

        if (x < min_x) min_x = x;
        if (y < min_y) min_y = y;
        if (x > max_x) max_x = x;
        if (y > max_y) max_y = y;

        ids.push_back(i);
    }
    const double cx = (min_x + max_x) / 2;
    const double cy = (min_y + max_y) / 2;
    double min_dist = std::numeric_limits<double>::max();

    std::size_t i0 = INVALID_INDEX;
    std::size_t i1 = INVALID_INDEX;
    std::size_t i2 = INVALID_INDEX;

    // pick a seed point close to the centroid
    for (std::size_t i = 0; i < n; i++) {
        const double d = dist(cx, cy, coords[2 * i], coords[2 * i + 1]);
        if (d < min_dist) {
            i0 = i;
            min_dist = d;
        }
    }

    min_dist = std::numeric_limits<double>::max();

    // find the point closest to the seed
    for (std::size_t i = 0; i < n; i++) {
        if (i == i0) continue;
        const double d = dist(coords[2 * i0], coords[2 * i0 + 1], coords[2 * i], coords[2 * i + 1]);
        if (d < min_dist && d > 0.0) {
            i1 = i;
            min_dist = d;
        }
    }

    double min_radius = std::numeric_limits<double>::max();

    // find the third point which forms the smallest circumcircle with the first two
    for (std::size_t i = 0; i < n; i++) {
        if (i == i0 || i == i1) continue;

        const double r = circumradius(
            coords[2 * i0], coords[2 * i0 + 1], coords[2 * i1], coords[2 * i1 + 1], coords[2 * i], coords[2 * i + 1]);

        if (r < min_radius) {
            i2 = i;
            min_radius = r;
        }
    }

    if (!(min_radius < std::numeric_limits<double>::max())) {
        throw std::runtime_error("not triangulation");
        ;
    }

    bool coord_area = area(
                          coords[2 * i0], coords[2 * i0 + 1], coords[2 * i1], coords[2 * i1 + 1], coords[2 * i2], coords[2 * i2 + 1]) < 0.0;
    if (coord_area) {
        std::swap(i1, i2);
    }

    const double i0x = coords[2 * i0];
    const double i0y = coords[2 * i0 + 1];
    const double i1x = coords[2 * i1];
    const double i1y = coords[2 * i1 + 1];
    const double i2x = coords[2 * i2];
    const double i2y = coords[2 * i2 + 1];

    std::tie(m_center_x, m_center_y) = circumcenter(i0x, i0y, i1x, i1y, i2x, i2y);

    // sort the points by distance from the seed triangle circumcenter
    // cerr << ids << endl;
    std::sort(ids.begin(), ids.end(), sort_to_center{ coords, m_center_x, m_center_y });
    // quicksort(ids, coords, 0, n - 1, m_center_x, m_center_y);
    // cerr << ids << endl;

    m_hash_size = static_cast<std::size_t>(std::llround(std::ceil(std::sqrt(n))));
    m_hash.reserve(m_hash_size);
    for (std::size_t i = 0; i < m_hash_size; i++) {
        m_hash.push_back(INVALID_INDEX);
    }

    m_hull.reserve(coords.size());

    std::size_t e = insert_node(i0);
    hash_edge(e);
    m_hull[e].t = 0;

    e = insert_node(i1, e);
    hash_edge(e);
    m_hull[e].t = 1;

    e = insert_node(i2, e);
    hash_edge(e);
    m_hull[e].t = 2;

    std::size_t max_triangles = n < 3 ? 1 : 2 * n - 5;
    triangles.reserve(max_triangles * 3);
    halfedges.reserve(max_triangles * 3);
    add_triangle(i0, i1, i2, INVALID_INDEX, INVALID_INDEX, INVALID_INDEX);
    double xp = std::numeric_limits<double>::quiet_NaN();
    double yp = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t k = 0; k < n; k++) {
        const std::size_t i = ids[k];
        const double x = coords[2 * i];
        const double y = coords[2 * i + 1];
        if (check_pts_equal(x, y, xp, yp)) continue;
        xp = x;
        yp = y;
        if (
            check_pts_equal(x, y, i0x, i0y) ||
            check_pts_equal(x, y, i1x, i1y) ||
            check_pts_equal(x, y, i2x, i2y)) continue;

        const std::size_t start_key = hash_key(x, y);
        std::size_t key = start_key;
        std::size_t start = INVALID_INDEX;
        do {
            start = m_hash[key];
            key = (key + 1) % m_hash_size;
        } while (
            (start == INVALID_INDEX || m_hull[start].removed) &&
            (key != start_key));

        e = start;

        while (
            area(
                x, y, m_hull[e].x, m_hull[e].y, m_hull[m_hull[e].next].x, m_hull[m_hull[e].next].y) >= 0.0) {
            e = m_hull[e].next;

            if (e == start) {
                throw std::runtime_error("Something is wrong with the input points.");
            }
        }

        const bool walk_back = e == start;

        // add the first triangle from the point
        std::size_t t = add_triangle(
            m_hull[e].i,
            i,
            m_hull[m_hull[e].next].i,
            INVALID_INDEX,
            INVALID_INDEX,
            m_hull[e].t);

        m_hull[e].t = t; // keep track of boundary triangles on the hull
        e = insert_node(i, e);

        // recursively flip triangles from the point until they satisfy the Delaunay condition
        m_hull[e].t = legalize(t + 2);
        if (m_hull[m_hull[m_hull[e].prev].prev].t == halfedges[t + 1]) {
            m_hull[m_hull[m_hull[e].prev].prev].t = t + 2;
        }
        // walk forward through the hull, adding more triangles and flipping recursively
        std::size_t q = m_hull[e].next;
        while (
            area(
                x, y, m_hull[q].x, m_hull[q].y, m_hull[m_hull[q].next].x, m_hull[m_hull[q].next].y) < 0.0) {
            t = add_triangle(
                m_hull[q].i, i, m_hull[m_hull[q].next].i, m_hull[m_hull[q].prev].t, INVALID_INDEX, m_hull[q].t);
            m_hull[m_hull[q].prev].t = legalize(t + 2);
            remove_node(q);
            q = m_hull[q].next;
        }
        if (walk_back) {
            // walk backward from the other side, adding more triangles and flipping
            q = m_hull[e].prev;
            while (
                area(
                    x, y, m_hull[m_hull[q].prev].x, m_hull[m_hull[q].prev].y, m_hull[q].x, m_hull[q].y) < 0.0) {
                t = add_triangle(
                    m_hull[m_hull[q].prev].i, i, m_hull[q].i, INVALID_INDEX, m_hull[q].t, m_hull[m_hull[q].prev].t);
                legalize(t + 2);
                m_hull[m_hull[q].prev].t = t;
                remove_node(q);
                q = m_hull[q].prev;
            }
        }
        hash_edge(e);
        hash_edge(m_hull[e].prev);
    }
}

std::size_t Delaunator::remove_node(std::size_t node) {
    m_hull[m_hull[node].prev].next = m_hull[node].next;
    m_hull[m_hull[node].next].prev = m_hull[node].prev;
    m_hull[node].removed = true;
    return m_hull[node].prev;
}

std::size_t Delaunator::legalize(std::size_t a) {
    std::size_t b = halfedges[a];

    std::size_t a0 = a - a % 3;
    std::size_t b0 = b - b % 3;

    std::size_t al = a0 + (a + 1) % 3;
    std::size_t ar = a0 + (a + 2) % 3;
    std::size_t bl = b0 + (b + 2) % 3;

    std::size_t p0 = triangles[ar];
    std::size_t pr = triangles[a];
    std::size_t pl = triangles[al];
    std::size_t p1 = triangles[bl];

    const bool illegal = in_circle(
        coords[2 * p0], coords[2 * p0 + 1], coords[2 * pr], coords[2 * pr + 1], coords[2 * pl], coords[2 * pl + 1], coords[2 * p1], coords[2 * p1 + 1]);

    if (illegal) {
        triangles[a] = p1;
        triangles[b] = p0;
        link(a, halfedges[bl]);
        link(b, halfedges[ar]);
        link(ar, bl);

        std::size_t br = b0 + (b + 1) % 3;

        legalize(a);
        return legalize(br);
    }
    return ar;
}

std::size_t Delaunator::insert_node(std::size_t i) {
    std::size_t node = m_hull.size();
    DelaunatorPoint p = {
        i,
        coords[2 * i],
        coords[2 * i + 1],
        0,
        node,
        node,
        false
    };
    m_hull.push_back(p);
    return node;
}

std::size_t Delaunator::insert_node(std::size_t i, std::size_t prev) {
    std::size_t node = insert_node(i);
    m_hull[node].next = m_hull[prev].next;
    m_hull[node].prev = prev;
    m_hull[m_hull[node].next].prev = node;
    m_hull[prev].next = node;
    return node;
}

std::size_t Delaunator::hash_key(double x, double y) {
    const double dx = x - m_center_x;
    const double dy = y - m_center_y;
    // use pseudo-angle: a measure that monotonically increases
    // with real angle, but doesn't require expensive trigonometry
    const double p = 1.0 - dx / (std::abs(dx) + std::abs(dy));
    return static_cast<std::size_t>(std::llround(std::floor(
        (2.0 + (dy < 0.0 ? -p : p)) / 4.0 * static_cast<double>(m_hash_size) //TODO:is this conversion save?
        )));
}

void Delaunator::hash_edge(std::size_t e) {
    m_hash[hash_key(m_hull[e].x, m_hull[e].y)] = e;
}

std::size_t Delaunator::add_triangle(
    std::size_t i0,
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

void Delaunator::link(std::size_t a, std::size_t b) {
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

} //namespace delaunator
