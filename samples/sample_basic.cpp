#include <iostream>

import UBKLib;

using T = double;

constexpr std::size_t MIN_FIELD_LINE_POINT_COUNT = 1000;

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
  } catch(std::runtime_error) {
    std::cout << "Couldn't even take one step!";
    return 0;
  }

  try {
    calculateLongitudinalInvariants(fieldLine);
  } catch(ubk::BifurcatingFieldLine) {
    std::cout << "Bifurcating field line!";
    return 0;
  }

  if (fieldLine.points().size() < MIN_FIELD_LINE_POINT_COUNT) {
    std::cout << "Field line did not trace the whole length!";
    return 0;
  }

  for (const auto& point : fieldLine.points()) {
    std::cout << ubk::vec3ToStr(point.loc) << '\n';
  }
  
  std::cout << "Total valid points traced: " << fieldLine.points().size();

  return 0;
}
