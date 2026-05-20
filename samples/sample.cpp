#include <iostream>

import UBKLib;

#define NUM_FIELD_LINES_TO_TRACE 100
#define K_VAL 100

using T = double;

constexpr std::size_t MIN_FIELD_LINE_POINT_COUNT = 50;

struct Field {
  ubk::Igrf13<T> igrf13;
  ubk::Ts89<T> ts89;

  [[nodiscard]] ubk::Vector3<ubk::nanoTesla<T>>
  getField(ubk::Vector3<T> pos) const {
    return ts89.getField(pos) + igrf13.getField(pos);
  }

};

int main() {
  constexpr ubk::FieldLineParams<T> params = {
    .innerLim = 1.05,
    .outterLim = 15.0,
    .maxStepDotField = 0.01,
    .failRatio = 2,
    .maxStepSize = 0.01,
    .maxStepCount = 10000,
  };
  
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
  while(i < NUM_FIELD_LINES_TO_TRACE) {
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
  
    if (fieldLine.maxLongitudinalInvariant() < K_VAL) continue;
    points = fieldLine.getPointsWithK(K_VAL);

    std::cout << vec3ToStr(points[0].loc) << '\n';
    std::cout << vec3ToStr(points[1].loc) << '\n';
    totalValidPointsTraced += fieldLine.points().size();
    i++;
  }
  
  std::cout << "Valid field lines traced: " << NUM_FIELD_LINES_TO_TRACE << '\n';
  std::cout << "Total valid points traced: " << totalValidPointsTraced << '\n';
  std::cout << "Bifurcating field lines: " << bifurcatingFieldLines << '\n';
  std::cout << "Short field lines: " << shortFieldLines;

  return 0;
}
