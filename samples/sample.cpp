#include <iostream>

import UBKLib;

using T = double;

struct Field {
  ubk::Dipole<T> dipole;

  [[nodiscard]] ubk::Vector3<ubk::nanoTesla<T>>
  getField(ubk::Vector3<T> pos) const {
    return dipole.getField(pos);
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
  ubk::UniformEquatorGenerator<T> rng(1.0, 15.0);
  
  for (int i = 1; i < 100; i++) {
    auto seed = rng.gen();
    ubk::FieldLine<T, Field, params> fieldLine;
    //std::cout << "generating field line" << std::endl;
    try {
      fieldLine = generator.generateFieldLine(seed);
    } catch(ubk::BifercatingFieldLine& e) {
      (void)e;
      continue;
    }

    if (fieldLine.points().size() < 1000 ||
        fieldLine.points()[0].loc.ampSquared() < 1.1 ||
        fieldLine.points().back().loc.ampSquared() < 1.1) {
      continue;
    }

    //std::cout << "calculating invariant" << std::endl;
    calculateLongitudinalInvariants(fieldLine);
    
    std::array<ubk::FieldLine<T, Field, params>::UBKInfos, 2> points;
  
    if (fieldLine.maxLongitudinalInvariant() < 100) { continue; }
    //std::cout << "getting k points" << std::endl;
    points = fieldLine.getPointsWithK(100);

    std::cout << vec3ToStr(points[0].loc) << '\n';
    std::cout << vec3ToStr(points[1].loc) << '\n';
  }

  return 0;
}
