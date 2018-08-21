#pragma once
#include <vector>
#include <exception>

using namespace std;

class Delaunator{
    public:
        Delaunator(const vector<double> &coords);
        struct No_triangulation : public exception {
            const char * what () const throw () {
                return "No delaunay triangulation";
            };
        };
    private:
        double m_center_x;
        double m_center_y;
        double m_hash_size;
        vector<int> m_hash;
};
