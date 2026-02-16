#include <iostream>
#include <iomanip>

import UBKLib;

using T = double;

int main() {
  constexpr ubk::FieldLineParams<T> params = {
    .innerLim = 1.05,
    .outterLim = 15.0,
    .maxStepDotField = 0.01,
    .failRatio = 2,
    .maxStepSize = 0.01,
    .maxStepCount = 10000,
  };

  ubk::FieldLineGenerator<T, ubk::Dipole<T>, params> generator;
  ubk::UniformEquatorGenerator<T> rng(1.0, 15.0);
  
  for (int i = 0; i < 100; i++) {
    auto seed = rng.gen();
    std::cout << std::setprecision(14);
    ubk::FieldLine<T, ubk::Dipole<T>, params> fieldLine = generator.generateFieldLine(seed);
    if (fieldLine.points().size() < 1000 ||
        fieldLine.points()[0].loc.ampSquared() < 1.1 ||
        fieldLine.points().back().loc.ampSquared() < 1.1) {
      throw;
    }
    calculateLongitudinalInvariants(fieldLine);
    auto point = fieldLine.getMinima();
    std::cout << vec3ToStr(point.loc) << '\n';
  }

  return 0;
}
