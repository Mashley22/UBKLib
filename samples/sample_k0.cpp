#include <iostream>

import UBKLib;

#define NUM_FIELD_LINES_TO_TRACE 100

using T = double;

constexpr ubk::FieldLineParams<T> params = {
  .innerLim = 1.05,
  .outterLim = 15.0,
  .maxStepDotField = 0.01,
  .failRatio = 2,
  .maxStepSize = 0.01,
  .maxStepCount = 10000,
};

int main() {

  ubk::FieldLineGenerator<T, ubk::Dipole<T>, params> generator;
  ubk::UniformEquatorGenerator<T> rng(1.0, 15.0);

  std::size_t totalValidPointsTraced = 0;

  std::size_t i = 0;
  while(i < NUM_FIELD_LINES_TO_TRACE) {
    auto seed = rng.gen();
    ubk::FieldLine<T, ubk::Dipole<T>, params> fieldLine; 

    try {
      fieldLine = generator.generateFieldLine(seed);
    } catch(ubk::BifercatingFieldLine& e) {
      (void)e;
      continue;
    }

    if (fieldLine.points()[0].loc.ampSquared() < 1.1 ||
        fieldLine.points().back().loc.ampSquared() < 1.1) {
      continue;
    }

    auto point = fieldLine.getMinima();

    std::cout << vec3ToStr(point.loc) << '\n';
    totalValidPointsTraced += fieldLine.points().size();
    i++;
  }

  std::cout << "Total valid points traced: " << totalValidPointsTraced;

  return 0;
}
