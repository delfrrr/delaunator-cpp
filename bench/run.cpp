#include <benchmark/benchmark.h>
#include <string>
#include <delaunator.hpp>
#include "../examples/utils.hpp"

namespace {
    void BM_simple(benchmark::State& state)
    {
        std::string points_str = utils::read_file("./test/test-files/osm-nodes-45331-epsg-3857.geojson");
        std::vector<double> coords = utils::get_geo_json_points(points_str);

        while (state.KeepRunning())
        {
            delaunator::Delaunator delaunator(coords);
        }
    }
}

int main(int argc, char* argv[])
{
    benchmark::RegisterBenchmark("BM_simple", BM_simple);
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    return 0;
}
