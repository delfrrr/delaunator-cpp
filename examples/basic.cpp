#include <delaunator.hpp>
#include <cstdio>

// Define your point type
struct my_point {
    double x;
    double y;
};

int main() {
    /* x0, y0, x1, y1, ... */
    std::vector<my_point> coords = {{-1, 1}, {1, 1}, {1, -1}, {-1, -1}};

    //triangulation happens here
    auto d = delaunator::delaunator<my_point>(coords);

    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        printf(
            "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
            coords[d.triangles[i]].x,        //tx0
            coords[d.triangles[i]].y,        //ty0
            coords[d.triangles[i + 1]].x,    //tx1
            coords[d.triangles[i + 1]].y,    //ty1
            coords[d.triangles[i + 2]].x,    //tx2
            coords[d.triangles[i + 2]].y     //ty2
        );
    }
}
