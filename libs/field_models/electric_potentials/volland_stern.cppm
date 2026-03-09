module;

#include <concepts>
#include <UBK/macros.hpp>

export module UBKLib.field_models:electric_potentials.volland_stern;

import :traits;
import UBKLib.utils;

export namespace ubk {

template<std::floating_point T>
class VollandStern : public TimeDependentField<T> {
public:
  using Base = TimeDependentField<T>;

  VollandStern(void) UBK_NOEXCEPT = default;

  VollandStern(int kp, const Time& time)
    : Base(time), m_kp(kp) {}

  static constexpr kV<T> SURFACE_POTENTIAL = static_cast<T>(94.2);

  [[nodiscard]] static constexpr T
  E0(T kp) {
    return 0.045 / pow(1 - 0.159 * kp + 0.0093 * pow(kp, 2), 3);
  }
  
  // from Maynard and Chen (1975), UBKikuity
  [[nodiscard]] kV<T>
  getField(T mlt, T theta) const {
    T sinTheta_2 = pow(sin(theta), 2);
    return -1 * SURFACE_POTENTIAL * sinTheta_2 - E0() * sin((mlt - 12) * PI<T> / 12) / pow(sinTheta_2, 2);
  }

  [[nodiscard]] constexpr T
  E0(void) const {
    return E0(m_kp);
  }

  [[nodiscard]] kV<T>
  getField(const Vector3<Re<T>> pos) const {
    return getField(mltFromGeo(pos, Base::time_es()), SphericalPolar(pos).theta);
  }

private:
  T m_kp{0};
};

static_assert(ElectricPotentialModel<VollandStern<double>, double>);
static_assert(ElectricPotentialModel<VollandStern<float>, float>);

static_assert(!ElectricPotentialModel<VollandStern<float>, double>);
static_assert(!ElectricPotentialModel<VollandStern<double>, float>);

}
