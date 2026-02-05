import UBKLib;

#include <iostream>

int main() {
  constexpr ubk::FieldLineParams<double> params = {
    .innerLim = 1.05,
    .outterLim = 15,
    .maxStepDotField = 0.01,
    .failRatio = 2,
    .maxStepCount = 10000,
  };
  
  ubk::FieldLineGenerator<double, ubk::Dipole<double>, params> generator;
  
  for (int i = 1; i < 10000; i++) {
    ubk::FieldLine<double, ubk::Dipole<double>, params> fieldLine = generator.generateFieldLine({2.0, 1.0 / i, 1.0 / i});
    if (fieldLine.points().size() < 1000 ||
        fieldLine.points()[0].ampSquared() < 1.1 ||
        fieldLine.points().back().ampSquared() < 1.1) {
      throw;
    }
    std::cout << fieldLine.points().size() << '\n';
  }

  return 0;
}
