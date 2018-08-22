#pragma once
#include <vector>
#include <exception>

class Delaunator{
    public:
        Delaunator(const std::vector<double> &coords);
    private:
        double m_center_x;
        double m_center_y;
        double m_hash_size;
        std::vector<int> m_hash;
};
