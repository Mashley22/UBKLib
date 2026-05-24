#include <thread>
#include <mutex>
#include <vector>
#include <atomic>
#include <filesystem>

import UBKLib;

#define DEFAULT_NUM_FIELD_LINES_TO_TRACE 1000
#define DEFAULT_K_VAL 100
#define DEFAULT_THREAD_COUNT 8

#define NUM_FIELD_LINES_ARG_IDX 1
#define THREAD_COUNT_ARG_IDX 2

using T = double;

constexpr std::size_t MIN_FIELD_LINE_POINT_COUNT = 50;

struct Point {
  T x; 
  T y;
  T B;
};

struct constKEquatorialSurface {
  T k_val;
  std::vector<Point> points;
};

static_assert(sizeof(Point) == 3 * sizeof(T));

static std::vector<constKEquatorialSurface> m_constKEquatorialSurfaces;
static unsigned long long num = DEFAULT_NUM_FIELD_LINES_TO_TRACE;
static unsigned long thread_count = DEFAULT_THREAD_COUNT;

struct Field {
  ubk::Igrf13<T> igrf13;
  ubk::Ts89<T> ts89;

  [[nodiscard]] ubk::Vector3<ubk::nanoTesla<T>>
  getField(ubk::Vector3<T> pos) const {
    return ts89.getField(pos) + igrf13.getField(pos);
  }
};

struct SharedResults {
  std::mutex mtx;
  std::size_t validLines = 0;
};

void worker(SharedResults& results, std::atomic<bool>& done){
  constexpr ubk::FieldLineParams<T> params = {
    .innerLim = 1.05,
    .outterLim = 20.0,
    .maxStepDotField = 0.01,
    .failRatio = 5,
    .maxStepSize = 0.005,
    .maxStepCount = 50000,
  };
  
  constexpr ubk::Time time{
    .year = 2010,
    .month = 2,
    .day = 15,
    .hours = 0,
    .minutes = 0,
    .seconds = 0
  };
  Field field;
  field.igrf13.setTime(time);
  field.ts89.dipole_tilt() = field.igrf13.dipole_tilt();
  field.ts89.setTime(time);
  field.ts89.iop() = 0;

  ubk::FieldLineGenerator<T, Field, params> generator;
  generator.assignModel(field);

  ubk::UniformEquatorGenerator<T> rng(1.0, 16.0);

  std::vector<Point> k_points;
  k_points.reserve(m_constKEquatorialSurfaces.size());
  
  while(!done.load()) {
    k_points.resize(0);
    auto seed = rng.gen();
    ubk::FieldLine<T, Field, params> fieldLine;

    try {
      fieldLine = generator.generateFieldLine(seed);
    } catch(std::runtime_error& e) {
      continue;
    }

    if (fieldLine.points().back().loc.amp() > 1.1 &&
        fieldLine.points().front().loc.amp() > 1.1) {
      continue;
    }

    if (fieldLine.points().size() < MIN_FIELD_LINE_POINT_COUNT) {
      continue;
    }
    
    try {
      calculateLongitudinalInvariants(fieldLine);
    } catch(ubk::BifurcatingFieldLine) {
      continue;
    } catch(ubk::NoMinimaFound) {
      continue;
    }
    
    for (std::size_t i = 0; i < m_constKEquatorialSurfaces.size(); i++) {
      if (fieldLine.maxLongitudinalInvariant() < m_constKEquatorialSurfaces[i].k_val) break;
      
      auto points = fieldLine.getPointsWithK(m_constKEquatorialSurfaces[i].k_val);
      k_points.push_back({
        .x = seed.x,
        .y = seed.y,
        .B = points[1].magneticIntensity
      });
    }

    if (k_points.empty()) continue;

    std::lock_guard lock(results.mtx);

    results.validLines++;
    if (results.validLines >= num) {
      done.store(true);
    }

    for (std::size_t i = 0; i < k_points.size(); i++) {
      m_constKEquatorialSurfaces[i].points.push_back(k_points[i]);
    }
    
  }
}

int main(int argc, char** argv) {
  SharedResults results;
  std::atomic<bool> done{false};
  std::vector<std::thread> threads;

  if (argc >= THREAD_COUNT_ARG_IDX + 1) {
    num = std::stoull(argv[NUM_FIELD_LINES_ARG_IDX]);
    thread_count = std::stoul(argv[THREAD_COUNT_ARG_IDX]);

    for (int i = THREAD_COUNT_ARG_IDX + 1; i < argc; i++) {
      m_constKEquatorialSurfaces.push_back({.k_val = std::stod(argv[i]), .points = {}});
    }
  }

  if (m_constKEquatorialSurfaces.empty()) {
    m_constKEquatorialSurfaces.push_back({.k_val = DEFAULT_K_VAL, .points = {}});
  }
  
  for (std::size_t t = 0; t < thread_count; t++) {
    threads.emplace_back(worker, std::ref(results), std::ref(done));
  }
  
  for (auto& t : threads) {
    t.join();
  }

  std::filesystem::create_directories("data");

  for (std::size_t i = 0; i < m_constKEquatorialSurfaces.size(); i++) {
    
    std::string filename{};
    if (argc > THREAD_COUNT_ARG_IDX + 1) {
      filename = "data/" + std::string(argv[i + THREAD_COUNT_ARG_IDX]) + ".bin"; 
    }
    else {
      filename = "data/data.bin";
    }

    std::ofstream out(filename, std::ios::binary);

    out.write(
      reinterpret_cast<const char*>(m_constKEquatorialSurfaces[i].points.data()),
      static_cast<std::streamsize>(m_constKEquatorialSurfaces[i].points.size() * sizeof(Point))
    );

  }

  return 0;
}
