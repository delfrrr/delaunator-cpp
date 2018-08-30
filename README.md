# delaunator-cpp

A really fast C++ library for
[Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of 2D points.

delaunator-cpp is a C++ port from https://github.com/mapbox/delaunator a JavaScript implementation of very fast 2D Delaunay algorithm.

## Features

* Probably the fastest C++ open source 2D Delaunay implementation
* Roughly 3 times faster then JS version.
* Example showing triangulation of GeoJson points

## Usage

```CPP
#include "delaunator.h"
#include  <cstdio>
using namespace std;
//...
int main(int, char* argv[]) {
    //...
    const vector<double> coords = {/* x0, y0, x1, y1, ... */};

    //triangulation happens here
    //note moving points to constructor
    Delaunator delaunator(move(coords));

    for(long int i = 0; i < delaunator.triangles.size(); i+=3) {
        printf(
            "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
            delaunator.coords[2 * delaunator.triangles[i]],         //tx0
            delaunator.coords[2 * delaunator.triangles[i] + 1],     //ty0
            delaunator.coords[2 * delaunator.triangles[i + 1]],     //tx1
            delaunator.coords[2 * delaunator.triangles[i + 1] + 1], //ty1
            delaunator.coords[2 * delaunator.triangles[i + 2]],     //tx2
            delaunator.coords[2 * delaunator.triangles[i + 2] + 1], //ty2
        )
    }
}
```

For full example see `src/triangulate.cpp`

## TODO

* Benchmarks
* Unit tests
