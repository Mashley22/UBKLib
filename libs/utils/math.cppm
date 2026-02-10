module;

#include <concepts>
#include <cmath>

#include <UBK/macros.hpp>

#include <xoshiro.h>

export module UBKLib.utils:math;

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
  T x{0}, y{0}, z{0};
  
  [[nodiscard]] constexpr Vector3 
  operator+(const Vector3& other) const UBK_NOEXCEPT {
    return {
      .x = x + other.x,
      .y = y + other.y,
      .z = z + other.z
    };
  }

  [[nodiscard]] constexpr Vector3
  operator-(const Vector3& other) const UBK_NOEXCEPT {
    return {
      .x = x - other.x,
      .y = y - other.y,
      .z = z - other.z
    };
  }

  [[nodiscard]] constexpr Vector3 
  operator*(const T scalar) const UBK_NOEXCEPT {
    return {
      .x = scalar * x,
      .y = scalar * y,
      .z = scalar * z,
    };
  }

  [[nodiscard]] friend constexpr Vector3
  operator*(const T scalar, Vector3<T> vec) {
    return vec * scalar;
  }

  [[nodiscard]] constexpr Vector3 
  operator/(const T scalar) const UBK_NOEXCEPT {
    return {
      .x = x / scalar,
      .y = y / scalar,
      .z = z / scalar,
    };
  }

  [[nodiscard]] constexpr T
  dot(const Vector3& other) const UBK_NOEXCEPT {
    return x * other.x + y * other.y + z * other.z;
  }

  [[nodiscard]] constexpr Vector3&
  operator+= (const Vector3& other) UBK_NOEXCEPT {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }

  [[nodiscard]] constexpr Vector3&
  operator-= (const Vector3& other) UBK_NOEXCEPT {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }

  [[nodiscard]] constexpr Vector3&
  operator*= (T scalar) UBK_NOEXCEPT {
    x *= scalar;
    y *= scalar;
    z *= scalar;
    return *this;
  }

  [[nodiscard]] constexpr Vector3&
  operator/= (T scalar) UBK_NOEXCEPT {
    x /= scalar;
    y /= scalar;
    z /= scalar;
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

  bool operator==(const Vector3&) const UBK_NOEXCEPT = default;
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
      .x = r * std::sin(theta) * std::cos(phi),
      .y = r * std::sin(theta) * std::sin(phi),
      .z = r * std::cos(theta)  
    });
  }
};

template<std::floating_point T, 
  xso::Distribution Distx,
  xso::Distribution Disty,
  xso::Distribution Distz, 
  class Generator = xso::rng>
[[nodiscard]] Vector3<T>
genRndVector3(Distx& distx, Disty& disty, Distz& distz) {
  thread_local Generator gen;
  return { .x = static_cast<T>(distx(gen)),
           .y = static_cast<T>(disty(gen)),
           .z = static_cast<T>(distz(gen)) };
}

template<std::floating_point T,
  xso::Distribution Distx,
  xso::Distribution Disty,
  xso::Distribution Distz, 
  class Generator = xso::rng>
class RandomVector3Generator {
public:
  RandomVector3Generator();

  RandomVector3Generator(const Distx& distx, const Disty& disty, const Distz& distz) UBK_NOEXCEPT
    : m_distx(distx), m_disty(disty), m_distz(distz) {}

  RandomVector3Generator(const Distx& dist) UBK_NOEXCEPT
    : m_distx(dist), m_disty(dist), m_distz(dist) {}
 
  [[nodiscard]] Vector3<T>
  gen(void) {
    return { .x = static_cast<T>(m_distx(m_gen)),
             .y = static_cast<T>(m_disty(m_gen)),
             .z = static_cast<T>(m_distz(m_gen)) };
  }
  
  [[nodiscard]] const Distx&
  distx(void) const { return m_distx; }

  [[nodiscard]] const Distx&
  disty(void) const { return m_disty; }

  [[nodiscard]] const Distx&
  distz(void) const { return m_distz; }

  [[nodiscard]] Distx&
  distx(void) { return m_distx; }

  [[nodiscard]] Distx&
  disty(void) { return m_disty; }

  [[nodiscard]] Distx&
  distz(void) { return m_distz; }

private:
  static thread_local inline Generator m_gen{};
  Distx m_distx;
  Disty m_disty;
  Distz m_distz;
};

};
