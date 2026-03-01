#include <iostream>

import UBKLib;

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
  ubk::Vector3<ubk::Re<T>> seed = rng.gen();

  std::cout << "Tracing field line with seed point: " << ubk::vec3ToStr(seed) << '\n';
  
  ubk::FieldLine<T, ubk::Dipole<T>, params> fieldLine;

  try {
    fieldLine = generator.generateFieldLine(seed);
  } catch(ubk::BifercatingFieldLine& e) {
    (void)e;
    std::cout << "Bifercating field line!";
    return 0;
  }

  if (fieldLine.points()[0].loc.ampSquared() < 1.1 ||
      fieldLine.points().back().loc.ampSquared() < 1.1) {
    std::cout << "Field line did not trace the whole length!";
    return 0;
  }

  for (const auto& point : fieldLine.points()) {
    std::cout << ubk::vec3ToStr(point.loc) << '\n';
  }
  
  std::cout << "Total valid points traced: " << fieldLine.points().size();

  return 0;
}
