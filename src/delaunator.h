#pragma once
#include <vector>
#include <exception>
#include <deque>

struct DelaunatorPoint {
    long int i;
    double x;
    double y;
    long int t;
};

class Delaunator{
    public:
        Delaunator(const std::vector<double> &coords);
    private:
        double m_center_x;
        double m_center_y;
        double m_hash_size;
        std::vector<int> m_hash;
        std::deque<DelaunatorPoint> m_hull;
        long int hash_key(double x, double y);
        void hash_edge(std::deque<DelaunatorPoint>::iterator e);
};
