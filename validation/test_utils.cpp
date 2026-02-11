#include <catch2/catch_all.hpp>

import UBKLib.utils;

namespace ubk {

TEST_CASE( "utils_math_vector3", "[utils][math][Vector3]" ) { // not testing both doubles and floats

  constexpr Vector3<double> vec1 = {1.0, 0.1, 0.009};
  constexpr Vector3<double> vec2 = {3.0, 1.256, 90898.1};
  
  SECTION( "vector-vector addition" ) {
    constexpr Vector3<double> vec1add2 = {vec1.x + vec2.x, vec1.y + vec2.y, vec1.z + vec2.z};
    STATIC_REQUIRE(vec1 + vec2 == vec1add2);
    STATIC_REQUIRE(vec2 + vec1 == vec1add2);
  }

  SECTION( "vector-vector subtraction" ) {
    constexpr Vector3<double> vec1sub2 = {vec1.x - vec2.x, vec1.y - vec2.y, vec1.z - vec2.z};
    STATIC_REQUIRE(vec1 - vec2 == vec1sub2);
  }

  SECTION( "vector-scalar multiplication" ) {
    constexpr double scalar = 5.2;
    STATIC_REQUIRE(vec1 * scalar == Vector3{vec1.x * scalar, vec1.y * scalar, vec1.z * scalar});
    STATIC_REQUIRE(vec2 * scalar == Vector3{vec2.x * scalar, vec2.y * scalar, vec2.z * scalar});

    STATIC_REQUIRE(scalar * vec1 == Vector3{vec1.x * scalar, vec1.y * scalar, vec1.z * scalar});
    STATIC_REQUIRE(scalar * vec2 == Vector3{vec2.x * scalar, vec2.y * scalar, vec2.z * scalar});
  }

  SECTION( "vector-scalar division" ) {
    constexpr double scalar = 5.2;
    STATIC_REQUIRE(vec1 / scalar == Vector3{vec1.x / scalar, vec1.y / scalar, vec1.z / scalar});
    STATIC_REQUIRE(vec2 / scalar == Vector3{vec2.x / scalar, vec2.y / scalar, vec2.z / scalar});
  }
}

}
