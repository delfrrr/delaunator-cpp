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
        double center_x;
        double center_y;
};
