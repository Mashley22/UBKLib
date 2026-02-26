module;

#include <array>
#include <concepts>
#include <cmath>

#include <UBK/macros.hpp>

#include <xoshiro.h>

export module UBKLib.utils:math;

import :check;

export namespace ubk {

[[nodiscard]] constexpr std::size_t
factorial(std::size_t val) UBK_NOEXCEPT {
  std::size_t retVal = 1;
  for (std::size_t i = 2; i <= val; i++) {
    retVal *= i;
  }
  return retVal;
}

template<std::floating_point T>
[[nodiscard]] constexpr T
sin(T val) {
  if (std::is_constant_evaluated()) {
    T val2 = val * val;
    return val * (1 - val2 / factorial(3));
  }
  return std::sin(val);
}

template<std::floating_point T>
[[nodiscard]] constexpr T
cos(T val) {
  if (std::is_constant_evaluated()) {
    T val2 = val * val;
    return 1 - val2 / factorial(2);
  }
  return std::cos(val);
}

template<std::floating_point T>
[[nodiscard]] T
acos(T val) {
  return std::acos(val);
}

template<std::floating_point T>
[[nodiscard]] constexpr T
tan(T val) {
  return std::tan(val);
}

template<std::floating_point T>
[[nodiscard]] constexpr T
pow(T val, std::size_t exponent) {
  T retVal = 1;
  for (std::size_t i = 0; i < exponent; i++) {
    retVal *= val; 
  }
  return retVal;
}

template<std::floating_point T>
using Re = T;

template<std::floating_point T>
using kV = T;

template<std::floating_point T>
using nanoTesla = T;

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

  template<std::floating_point T_arr = T>
  [[nodiscard]] constexpr static Vector3<T>
  fromArr(const std::array<T_arr, 3> arr) UBK_NOEXCEPT {
    return {
      .x = static_cast<T>(arr[0]),
      .y = static_cast<T>(arr[1]),
      .z = static_cast<T>(arr[2]),
    };
  }

  template<std::floating_point T_arr = T>
  [[nodiscard]] constexpr std::array<T_arr, 3>
  toArr(void) const UBK_NOEXCEPT {
    std::array<T_arr, 3> arr;
    arr[0] = static_cast<T_arr>(x);
    arr[1] = static_cast<T_arr>(y);
    arr[2] = static_cast<T_arr>(z);
    return arr;
  }
};

/**
 *@brief pure struct to stll allow the initializer lists
 *
 */
template<std::floating_point T, typename dist_t = Re<T>>
struct SphericalPolar_t {
  T theta;
  T phi;
  dist_t r;

  [[nodiscard]] constexpr T
  ampSquared(void) const UBK_NOEXCEPT {
    return r * r;
  }

  [[nodiscard]] constexpr T
  amp(void) const UBK_NOEXCEPT {
    return r;
  };

  [[nodiscard]] explicit operator Vector3<dist_t>() const UBK_NOEXCEPT {
    return Vector3<dist_t>({
      .x = r * sin(theta) * cos(phi),
      .y = r * sin(theta) * sin(phi),
      .z = r * cos(theta)  
    });
  }
};

template<std::floating_point T, typename dist_t = Re<T>>
class SphericalPolar : public SphericalPolar_t<T, dist_t> {
public:
  using Base = SphericalPolar_t<T, dist_t>;

  SphericalPolar(const SphericalPolar_t<T, dist_t> base) UBK_NOEXCEPT : SphericalPolar_t<T, dist_t>(base) {}

  SphericalPolar(const Vector3<dist_t> cartesian) UBK_NOEXCEPT {

    check(cartesian.ampSquared() != 0);

    Base::r = cartesian.amp();
    
    Base::theta = acos(cartesian.z / Base::r);
    
    if (cartesian.x == 0 && cartesian.y == 0) {
      Base::phi = 0;
    }
    else {
      Base::phi = acos(cartesian.x / (Base::r * sin(Base::theta)));
    }
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
