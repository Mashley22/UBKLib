#include <iostream>

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
  
  for (int i = 1; i < 100; i++) {
    auto seed = rng.gen();
    ubk::FieldLine<T, ubk::Dipole<T>, params> fieldLine = generator.generateFieldLine(seed);
    if (fieldLine.points().size() < 1000 ||
        fieldLine.points()[0].loc.ampSquared() < 1.1 ||
        fieldLine.points().back().loc.ampSquared() < 1.1) {
      throw;
    }
    calculateLongitudinalInvariants(fieldLine);
    auto points = fieldLine.getPointsWithK(5.0);
    (void) points;
    std::cout << seed.x << "  " << seed.y << "    " << seed.z << '\n';
    std::cout << fieldLine.points().size() << '\n';
    std::cout << fieldLine.maxLongitudinalInvariant() << '\n';
  }

  return 0;
}
