# delaunator-cpp

A really fast C++ library for
[Delaunay triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation) of 2D points.

delaunator-cpp is a C++ port from https://github.com/mapbox/delaunator a JavaScript implementation of very fast 2D Delaunay algorithm.

[![Build Status](https://travis-ci.org/delfrrr/delaunator-cpp.svg?branch=master)](https://travis-ci.org/delfrrr/delaunator-cpp)
[![badge](https://mapbox.s3.amazonaws.com/cpp-assets/hpp-skel-badge_blue.svg)](https://github.com/mapbox/hpp-skel)

## Features

* Probably the fastest C++ open source 2D Delaunay implementation
* Roughly 6 times faster then JS version `delaunator@v2.0.3` ([needs confirmation](https://github.com/delfrrr/delaunator-cpp/pull/8))
* Example showing triangulation of GeoJson points

## Usage

`examples/basic.cpp`

```CPP
#include <delaunator.hpp>
#include <cstdio>

int main() {
    /* x0, y0, x1, y1, ... */
    std::vector<double> coords = {-1, 1, 1, 1, 1, -1, -1, -1};

    //triangulation happens here
    delaunator::Delaunator d(coords);

    for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
        printf(
            "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
            d.coords[2 * d.triangles[i]],        //tx0
            d.coords[2 * d.triangles[i] + 1],    //ty0
            d.coords[2 * d.triangles[i + 1]],    //tx1
            d.coords[2 * d.triangles[i + 1] + 1],//ty1
            d.coords[2 * d.triangles[i + 2]],    //tx2
            d.coords[2 * d.triangles[i + 2] + 1] //ty2
        );
    }
}
```

[See more examples here](./examples)
