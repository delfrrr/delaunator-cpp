#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <string>
#include <vector>

#include "delaunator.h"
#include "json-helpers.h"

// test data
#include "vertex_points.h"
#include "triangle_indices.h"

using namespace std;

TEST_CASE("regular grid")
{
    std::vector<double> points;
    
    int w = 100;
    int h = 100;
    
    for(int r = 0; r < h; r++)
    {
        for(int c = 0; c < w; c++)
        {
            points.push_back(c);
            points.push_back(r);
        }
    }
    
     Delaunator dn(points);
    
     std::cout << "number of triangles: " << dn.triangles.size()/3 << std::endl;
     std::cout << "number of vertices: " << points.size()/2 << std::endl;
     CHECK(dn.triangles.size()/3 == (w - 1) * (h - 1) * 2);
}

TEST_CASE("triangles match JS version ouput", "[Delaunator]")
{
    std::vector<double> vertex_list_big;
    unittests::init_vertex_list(vertex_list_big);
    
    Delaunator dn(vertex_list_big);
    
    std::vector<int> tri_big;
    unittests::init_triangle_list(tri_big);
    
    CHECK(dn.triangles.size() == tri_big.size());
    for(int i = 0; i < dn.triangles.size(); i++)
    {

        if(dn.triangles[i] != tri_big[i])
        {
            std::cout << "faile at i: " << i << std::endl;
        }
        
        CHECK(dn.triangles[i] == tri_big[i]);
    }
    
    CHECK(true);
}

