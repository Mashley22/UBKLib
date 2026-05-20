#include <iostream>
#include <thread>
#include <mutex>
#include <vector>
#include <atomic>

import UBKLib;

#define DEFAULT_NUM_FIELD_LINES_TO_TRACE 1000
#define DEFAULT_K_VAL 100
#define DEFAULT_THREAD_COUNT 8

#define K_VAL_ARG_IDX 1
#define NUM_FIELD_LINES_ARG_IDX 2
#define THREAD_COUNT_ARG_IDX 3

using T = double;

constexpr std::size_t MIN_FIELD_LINE_POINT_COUNT = 50;

static T k_val = DEFAULT_K_VAL;
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
  std::size_t totalValidPointsTraced = 0;
  std::size_t bifurcatingFieldLines = 0;
  std::size_t shortFieldLines = 0;
  std::size_t validLines = 0;
};

void worker(SharedResults& results, std::atomic<bool>& done){
  constexpr ubk::FieldLineParams<T> params = {
    .innerLim = 1.05,
    .outterLim = 15.0,
    .maxStepDotField = 0.01,
    .failRatio = 2,
    .maxStepSize = 0.01,
    .maxStepCount = 10000,
  };
  
  constexpr ubk::Time time{};
  Field field;
  field.igrf13.setTime(time);
  field.ts89.dipole_tilt() = field.igrf13.dipole_tilt();
  field.ts89.setTime(time);

  ubk::FieldLineGenerator<T, Field, params> generator;
  generator.assignModel(field);

  ubk::UniformEquatorGenerator<T> rng(1.0, 15.0);
  
  while(!done.load()) {
    auto seed = rng.gen();
    ubk::FieldLine<T, Field, params> fieldLine;

    try {
      fieldLine = generator.generateFieldLine(seed);
    } catch(std::runtime_error& e) {
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
    
    if (fieldLine.maxLongitudinalInvariant() < k_val) continue;
    
    auto points = fieldLine.getPointsWithK(k_val);

    {
      std::lock_guard lock(results.mtx);
      std::cout << seed.x << ',' << seed.y << ',' << points[1].magneticIntensity << '\n';
      results.totalValidPointsTraced += fieldLine.points().size();
      results.validLines++;
      if (results.validLines > num) {
        done.store(true);
      }
    }
    
  }
}

int main(int argc, char** argv) {
  SharedResults results;
  std::atomic<bool> done{false};
  std::vector<std::thread> threads;

  if (argc == 4) {
    k_val = std::stod(argv[K_VAL_ARG_IDX]);
    num = std::stoull(argv[NUM_FIELD_LINES_ARG_IDX]);
    thread_count = std::stoul(argv[THREAD_COUNT_ARG_IDX]);
  }
  
  for (std::size_t t = 0; t < thread_count; t++) {
    threads.emplace_back(worker, std::ref(results), std::ref(done));
  }
  
  for (auto& t : threads) {
    t.join();
  }
  
  std::cout << "\nValid field lines traced: " << results.validLines << '\n';
  std::cout << "Total valid points traced: " << results.totalValidPointsTraced << '\n';
  std::cout << "Bifurcating field lines: " << results.bifurcatingFieldLines << '\n';
  std::cout << "Short field lines: " << results.shortFieldLines << '\n';

  return 0;
}
