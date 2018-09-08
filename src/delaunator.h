#pragma once
#include <vector>
#include <exception>
#include <deque>
#include <memory>

struct DelaunatorPoint {
    int64_t i;
    double x;
    double y;
    int64_t t;
    int64_t prev;
    int64_t next;
    bool removed;
};

class Delaunator{
    public:
        Delaunator(std::vector<double>& in_coords);
    
        std::vector<uint64_t>  triangles;
        std::vector<int64_t>           halfedges;
        std::vector<double>             coords;
    
    private:
    
        Delaunator();
    
        double                          m_center_x;
        double                          m_center_y;
        int64_t                        m_hash_size;
        int64_t                        m_hull_index;
    
        std::vector<int>                m_hash;
        std::vector<DelaunatorPoint>    m_hl;
    
        double                          m_epilon;
 
    private:
    
        double      pseudo_angle(const double dx, const double dy);
        int64_t    insert_node(int64_t i);
        int64_t    insert_node(int64_t i, int64_t prev);
        int64_t    hash_key(double x, double y);
        void        hash_edge(int64_t e);
    
        int64_t    add_triangle(
            int64_t i0, int64_t i1, int64_t i2,
            int64_t a, int64_t b, int64_t c
        );
    
        void        link(int64_t a, int64_t b);
        int64_t    legalize(int64_t a, int64_t &e);
        int64_t    remove_node(int64_t node);
    
};
