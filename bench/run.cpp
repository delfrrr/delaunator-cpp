#include <benchmark/benchmark.h>
#include <string>
#include <delaunator.hpp>
#include "../test/test_utils.hpp"

static void BM_simple(benchmark::State& state) // NOLINT google-runtime-references
{
    std::string points_str = test_utils::read_file("./test/test-files/osm-nodes-45331-epsg-3857.geojson");
    std::vector<double> coords = test_utils::get_geo_json_points(points_str);
    
    while (state.KeepRunning())
    {
        delaunator::Delaunator delaunator(coords);
    }
}

int main(int argc, char* argv[])
{
    benchmark::RegisterBenchmark("BM_simple", BM_simple);      // NOLINT clang-analyzer-cplusplus.NewDeleteLeaks
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
