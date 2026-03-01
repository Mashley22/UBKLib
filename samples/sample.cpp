#include <iostream>

import UBKLib;

#define NUM_FIELD_LINES_TO_TRACE 100
#define K_VAL 100

using T = double;

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
  field.igrf13.init(time);
  field.ts89.dipole_tilt() = field.igrf13.dipole_tilt();
  field.ts89.setTime(time);

  ubk::FieldLineGenerator<T, Field, params> generator;
  generator.assignModel(field);

  ubk::UniformEquatorGenerator<T> rng(1.0, 15.0);
  
  std::size_t totalValidPointsTraced = 0;
  std::size_t bifercatingFieldLines = 0;
  std::size_t shortFieldLines = 0;

  std::size_t i = 0;
  while(i < NUM_FIELD_LINES_TO_TRACE) {
    auto seed = rng.gen();
    ubk::FieldLine<T, Field, params> fieldLine;

    try {
      fieldLine = generator.generateFieldLine(seed);
    } catch(ubk::BifercatingFieldLine& e) {
      (void)e;
      bifercatingFieldLines++;
      continue;
    }

    if (fieldLine.points()[0].loc.ampSquared() < 1.1 ||
        fieldLine.points().back().loc.ampSquared() < 1.1) {
      shortFieldLines++;
      continue;
    }

    calculateLongitudinalInvariants(fieldLine);
    
    std::array<ubk::FieldLine<T, Field, params>::UBKInfos, 2> points;
  
    if (fieldLine.maxLongitudinalInvariant() < K_VAL) continue;
    points = fieldLine.getPointsWithK(K_VAL);

    std::cout << vec3ToStr(points[0].loc) << '\n';
    std::cout << vec3ToStr(points[1].loc) << '\n';
    totalValidPointsTraced += fieldLine.points().size();
    i++;
  }
  
  std::cout << "Valid field lines traced: " << NUM_FIELD_LINES_TO_TRACE << '\n';
  std::cout << "Total valid points traced: " << totalValidPointsTraced << '\n';
  std::cout << "Bifercating field lines: " << bifercatingFieldLines << '\n';
  std::cout << "Short field lines: " << shortFieldLines;

  return 0;
}
