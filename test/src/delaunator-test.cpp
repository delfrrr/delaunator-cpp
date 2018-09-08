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

