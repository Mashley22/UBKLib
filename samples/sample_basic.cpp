import UBKLib;

#include <iostream>

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
  
  for (int i = 1; i < 100; i++) {
    ubk::FieldLine<T, ubk::Dipole<T>, params> fieldLine = generator.generateFieldLine({2.0, 1.0 / i, 1.0 / i});
    if (fieldLine.points().size() < 1000 ||
        fieldLine.points()[0].loc.ampSquared() < 1.1 ||
        fieldLine.points().back().loc.ampSquared() < 1.1) {
      throw;
    }
    calculateLongitudinalInvariants(fieldLine);
    std::cout << fieldLine.points().size() << '\n';
  }

  return 0;
}
