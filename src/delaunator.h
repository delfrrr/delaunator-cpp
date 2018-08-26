#pragma once
#include <vector>
#include <exception>
#include <deque>
#include <memory>

struct DelaunatorPoint {
    size_t i;
    double x;
    double y;
    long int t;
    long int prev;
    long int next;
    bool removed;
};

class Delaunator{
    public:
        Delaunator(const std::vector<double> &in_coords);
        std::vector<unsigned long int> triangles;
        std::vector<long int> halfedges;
        // size_t triangles_len;
    private:
        double m_center_x;
        double m_center_y;
        size_t m_hash_size;
        std::vector<int> m_hash;
        std::vector<double> m_coords;
        std::vector<DelaunatorPoint> m_hull;
        size_t insert_node(size_t i);
        size_t insert_node(size_t i, size_t prev);
        size_t hash_key(double x, double y);
        void hash_edge(size_t e);
        size_t add_triangle(
            size_t i0, size_t i1, size_t i2,
            size_t a, size_t b, size_t c
        );
        void link(size_t a, size_t b);
        size_t legalize(size_t a);
        size_t remove_node(size_t node);
};
