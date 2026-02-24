#include <catch2/catch_test_macros.hpp>

import UBKLib;

#define MIN_POINTS 1000
#define MAX_TERMINATION_DIST 1.1

namespace ubk {

namespace {

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
[[nodiscard]] constexpr bool
singleMinima(const FieldLine<T, FieldModel, Params>& fieldLine) {
  T minima = fieldLine.points().front().longitudinalInvariant;
  bool foundMinima = false;
  for (const auto& point : fieldLine.points()) {
    if (point.longitudinalInvariant > minima) {
      foundMinima = true;
    }
    else {
      if (foundMinima) {
        return false;
      }
      minima = point.longitudinalInvariant;
    }
  }

  return true;
}

}

TEST_CASE( "test tracer with dipole field", "[tracer][dipole]" ) {
  constexpr FieldLineParams<double> params = {
    .innerLim = 1.05,
    .outterLim = 15,
    .maxStepDotField = 0.01,
    .failRatio = 1.5,
    .maxStepCount = 10000,
  };
  
  FieldLineGenerator<double, Dipole<double>, params> generator;

  SECTION( "Generated field line is somewhat valid" ) {
    FieldLine<double, Dipole<double>, params> fieldLine = generator.generateFieldLine({2.0, 0.0, 0.0});
    REQUIRE(fieldLine.points().size() > 100);
    REQUIRE(fieldLine.points().front().loc.amp() <= MAX_TERMINATION_DIST);
    REQUIRE(fieldLine.points().back().loc.amp() <= MAX_TERMINATION_DIST);

    for (const auto& points : fieldLine.points()) {
      REQUIRE(points.loc.y == 0);
    }
  }

  SECTION( "K calculations are working" ) {
    FieldLine<double, Dipole<double>, params> fieldLine = generator.generateFieldLine({2.0, 0.0, 0.0});
    calculateLongitudinalInvariants(fieldLine);
    REQUIRE(fieldLine.points().back().longitudinalInvariant == fieldLine.points().front().longitudinalInvariant);
    REQUIRE(singleMinima(fieldLine));
  }

  SECTION( "K = 0 point for a dipole should have z = 0" ) {
    FieldLine<double, Dipole<double>, params> fieldLine = generator.generateFieldLine({2.0, 0.0, 0.0});
    calculateLongitudinalInvariants(fieldLine);
    auto minima = fieldLine.getMinima();
    REQUIRE(minima.loc.z == 0);
  }
}

}
