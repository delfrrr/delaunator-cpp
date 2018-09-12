#include <benchmark/benchmark.h>
#include <string>
#include <delaunator.hpp>
#include "../examples/utils.hpp"

namespace {
    void BM_45K_geojson_nodes(benchmark::State& state)
    {
        std::string points_str = utils::read_file("./test/test-files/osm-nodes-45331-epsg-3857.geojson");
        std::vector<double> coords = utils::get_geo_json_points(points_str);

        while (state.KeepRunning())
        {
            delaunator::Delaunator delaunator(coords);
        }
    }
}

BENCHMARK(BM_45K_geojson_nodes)->Unit(benchmark::kMillisecond);

BENCHMARK_MAIN();
