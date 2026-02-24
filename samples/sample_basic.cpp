#include <iostream>

import UBKLib;

using T = double;

struct Field {
  ubk::Ts89<T> ts89;
  ubk::Dipole<T> dipole;

  [[nodiscard]] ubk::Vector3<ubk::nanoTesla<T>>
  getField(ubk::Vector3<T> pos) const {
    return ts89.getField(pos) + dipole.getField(pos);
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
  
  ubk::FieldLineGenerator<T, Field, params> generator;
  ubk::FieldLineGenerator<T, ubk::Dipole<T>, params> generator2;
  
  for (int i = 0; i < 1; i++) {
    ubk::Vector3<ubk::Re<T>> seed = {2.0, 0.0, 0.0};
    ubk::FieldLine<T, Field, params> fieldLine = generator.generateFieldLine(seed);
    ubk::FieldLine<T, ubk::Dipole<T>, params> fieldLine2 = generator2.generateFieldLine(seed);
    if (fieldLine.points().size() < 1000 ||
        fieldLine.points()[0].loc.ampSquared() < 1.1 ||
        fieldLine.points().back().loc.ampSquared() < 1.1) {
      throw;
    }
    for (const auto& point : fieldLine.points()) {
      std::cout << ubk::vec3ToStr(point.loc) << '\n';
    }
    for (const auto& point : fieldLine2.points()) {
      std::cout << ubk::vec3ToStr(-1 * point.loc) << '\n';
    }
  }

  return 0;
}
