#include <iostream>

import UBKLib;

using T = double;

#define DEFAULT_NUM_FIELD_LINES_TO_TRACE 100
#define DEFAULT_K_VAL 100

#define NUM_FIELD_LINES_ARG_IDX 2
#define K_VAL_ARG_IDX 1

static unsigned long long num_field_lines = DEFAULT_NUM_FIELD_LINES_TO_TRACE;
static T k_val = DEFAULT_K_VAL;

constexpr std::size_t MIN_FIELD_LINE_POINT_COUNT = 50;

constexpr ubk::FieldLineParams<T> params = {
  .innerLim = 1.05,
  .outterLim = 15.0,
  .maxStepDotField = 0.01,
  .failRatio = 2,
  .maxStepSize = 0.01,
  .maxStepCount = 10000,
};

struct Field {
  ubk::Igrf13<T> igrf13;
  ubk::Ts89<T> ts89;

  [[nodiscard]] ubk::Vector3<ubk::nanoTesla<T>>
  getField(ubk::Vector3<T> pos) const {
    return ts89.getField(pos) + igrf13.getField(pos);
  }

};

int main(int argc, char** argv) {
  if (argc == 3) {
    k_val = std::stod(argv[K_VAL_ARG_IDX]);
    num_field_lines = std::stoull(argv[NUM_FIELD_LINES_ARG_IDX]);
  }

  ubk::Time time{};
  Field field;
  field.igrf13.setTime(time);
  field.ts89.dipole_tilt() = field.igrf13.dipole_tilt();
  field.ts89.setTime(time);

  ubk::FieldLineGenerator<T, Field, params> generator;
  generator.assignModel(field);

  ubk::UniformEquatorGenerator<T> rng(1.0, 15.0);
  
  std::size_t totalValidPointsTraced = 0;
  std::size_t bifurcatingFieldLines = 0;
  std::size_t shortFieldLines = 0;

  std::size_t i = 0;
  while(i < num_field_lines) {
    auto seed = rng.gen();
    ubk::FieldLine<T, Field, params> fieldLine;

    try {
      fieldLine = generator.generateFieldLine(seed);
    } catch(std::runtime_error) {
      continue;
    }

    if (fieldLine.points().size() < MIN_FIELD_LINE_POINT_COUNT) {
      shortFieldLines++;
      continue;
    }
    
    try {
      calculateLongitudinalInvariants(fieldLine);
    } catch(ubk::BifurcatingFieldLine) {
      bifurcatingFieldLines++;
      continue;
    } catch(ubk::NoMinimaFound) {
      continue;
    }
    
    std::array<ubk::FieldLine<T, Field, params>::PointInfo, 2> points;
  
    if (fieldLine.maxLongitudinalInvariant() < k_val) continue;
    points = fieldLine.getPointsWithK(k_val);

    std::cout << vec3ToStr(points[0].loc) << '\n';
    std::cout << vec3ToStr(points[1].loc) << '\n';
    totalValidPointsTraced += fieldLine.points().size();
    i++;
  }
  
  std::cout << "Valid field lines traced: " << num_field_lines << '\n';
  std::cout << "Total valid points traced: " << totalValidPointsTraced << '\n';
  std::cout << "Bifurcating field lines: " << bifurcatingFieldLines << '\n';
  std::cout << "Short field lines: " << shortFieldLines;

  return 0;
}
