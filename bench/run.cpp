#include "../examples/utils.hpp"
#include <benchmark/benchmark.h>
#include <delaunator.hpp>
#include <random>
#include <string>
#include <vector>

std::vector<utils::point> generate_uniform(std::size_t n) {
    std::vector<utils::point> coords;
    coords.reserve(n);
    std::srand(350);
    double norm = static_cast<double>(RAND_MAX) / 1e3;
    for (size_t i = 0; i < n; i++) {
        coords.push_back({(static_cast<double>(std::rand()) / norm), (static_cast<double>(std::rand()) / norm)});
    }
    return coords;
}

void BM_45K_geojson_nodes(benchmark::State& state) {
    std::string points_str = utils::read_file("./test/test-files/osm-nodes-45331-epsg-3857.geojson");
    auto coords = utils::get_geo_json_points(points_str);

    while (state.KeepRunning()) {
        auto result = delaunator::delaunator(coords);
        benchmark::DoNotOptimize(result);
    }   
}

void BM_uniform(benchmark::State& state) {
    std::vector<utils::point> coords = generate_uniform(static_cast<std::size_t>(state.range(0)));
    while (state.KeepRunning()) {
        auto result = delaunator::delaunator(coords);
        benchmark::DoNotOptimize(result);
    }
    state.SetComplexityN(state.range(0));
}

BENCHMARK(BM_45K_geojson_nodes)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_uniform)->Arg(2000)->Arg(100000)->Arg(200000)->Arg(500000)->Arg(1000000)->Unit(benchmark::kMillisecond);

#if BENCHMARK_BIG_O
BENCHMARK(BM_uniform)->RangeMultiplier(2)->Range(1 << 12, 1 << 22)->Unit(benchmark::kMillisecond)->Complexity();
#endif
#if BENCHMARK_10M
BENCHMARK(BM_uniform)->Arg(1000000 * 10)->Unit(benchmark::kMillisecond);
#endif
#if BENCHMARK_100M
BENCHMARK(BM_uniform)->Arg(1000000 * 100)->Unit(benchmark::kMillisecond);
#endif

BENCHMARK_MAIN()
