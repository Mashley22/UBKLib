module;
#include <concepts>
#include <cmath>

#include <UBK/macros.hpp>

export module ubk.utils:math;

import :check;

export namespace ubk {

template<std::floating_point T>
using Re = T;

template<std::floating_point T>
using kV = T;

template<std::floating_point T>
using microTesla = T;

template<std::floating_point T>
struct Vector3 {
  T x, y, z;
  
  constexpr
  Vector3(double x = 0, double y = 0, double z = 0) UBK_NOEXCEPT
    : x(x), y(y), z(z) {}

  [[nodiscard]] constexpr 
  Vector3 operator+(const Vector3& other) const UBK_NOEXCEPT {
    return {
      .x = x + other.x,
      .y = y + other.y,
      .z = z + other.z
    };
  }

  [[nodiscard]] constexpr 
  Vector3 operator-(const Vector3& other) const UBK_NOEXCEPT {
    return {
      .x = x - other.x,
      .y = y - other.y,
      .z = z - other.z
    };
  }

  [[nodiscard]] constexpr
  Vector3 operator*(const T scalar) const UBK_NOEXCEPT {
    return {
      .x = scalar * x,
      .y = scalar * y,
      .z = scalar * z,
    };
  }

  [[nodiscard]] constexpr
  Vector3 operator/(const T scalar) const UBK_NOEXCEPT {
    return {
      .x = scalar / x,
      .y = scalar / y,
      .z = scalar / z,
    };
  }

  [[nodiscard]] constexpr T
  dot(const Vector3& other) const UBK_NOEXCEPT {
    return x * other.x + y * other.y + z * other.z;
  }

  [[nodiscard]] constexpr
  Vector3& operator+= (const Vector3& other) UBK_NOEXCEPT {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  [[nodiscard]] constexpr
  Vector3& operator-= (const Vector3& other) UBK_NOEXCEPT {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  [[nodiscard]] constexpr
  Vector3& operator*= (const Vector3& other) UBK_NOEXCEPT {
    x *= other.x;
    y *= other.y;
    z *= other.z;
    return *this;
  }

  [[nodiscard]] constexpr
  Vector3& operator/= (const Vector3& other) UBK_NOEXCEPT {
    x /= other.x;
    y /= other.y;
    z /= other.z;
    return *this;
  }

  [[nodiscard]] constexpr T
  ampSquared(void) const UBK_NOEXCEPT {
    return x * x + y * y + z * z;
  }

  [[nodiscard]] T
  amp(void) const UBK_NOEXCEPT {
    return std::sqrt(ampSquared()); // in c++ 26 ...
  }

  [[nodiscard]] constexpr Vector3 
  cross(const Vector3& other) const UBK_NOEXCEPT {
    return {
      .x = y * other.z - z * other.y,  
      .y = z * other.x - x * other.z,  
      .z = x * other.y - y * other.x   
    };
  }

  [[nodiscard]] Vector3
  normalised(void) const UBK_NOEXCEPT {
    return Vector3(*this/amp());
  }
};

template<std::floating_point T, typename dist_t = Re<T>>
class SphericalPolar {
public:
  T theta;
  dist_t r;
  T phi;

  [[nodiscard]] constexpr T
  ampSquared(void) const UBK_NOEXCEPT {
    return r * r;
  }

  [[nodiscard]] constexpr T
  amp(void) const UBK_NOEXCEPT {
    return r;
  };

  [[nodiscard]] operator Vector3<T>() const UBK_NOEXCEPT {
    return Vector3<T>({
      .x = r * std::sin(phi) * std::cos(theta),
      .y = r * std::sin(phi) * std::sin(theta),
      .z = r * std::sin(phi)  
    });
  }
};

};
