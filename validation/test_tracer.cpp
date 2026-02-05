import UBKLib;

#include <catch2/catch_test_macros.hpp>

namespace ubk {

TEST_CASE( "test tracer with dipole field", "[tracer][dipole]" ) {
  constexpr FieldLineParams<double> params = {
    .innerLim = 1.05,
    .outterLim = 15,
    .maxStepDotField = 0.01,
    .failRatio = 1.5,
    .maxStepCount = 1000,
  };
  
  FieldLineGenerator<double, Dipole<double>, params> generator;

  FieldLine<double, Dipole<double>, params> fieldLine = generator.generateFieldLine({2.0, 0.0, 0.0});
  REQUIRE(fieldLine.points().size() > 1000);
}

}
