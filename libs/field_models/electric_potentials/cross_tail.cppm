module;

#include <concepts>
#include <numbers>

export module UBKLib.field_models:electric_potentials.cross_tail;

import :traits;
import UBKLib.utils;

export namespace ubk {

template<std::floating_point T>
class CrossTailPotential : public TimeDependentField<T> {
public:
  using Base = TimeDependentField<T>;
  static constexpr kV<T> SURFACE_POTENTIAL = static_cast<T>(92.0); 
  static constexpr T TAIL_FIELD = static_cast<T>(10.0); 

  [[nodiscard]] kV<T>
  getField(T mlt, T theta) const {
    T sinTheta_2 = sin(theta);

    return -1 * SURFACE_POTENTIAL * sinTheta_2 - sin((mlt - 12) * std::numbers::pi / 12) * TAIL_FIELD / sinTheta_2;
  }

  [[nodiscard]] kV<T>
  getField(const Vector3<Re<T>> pos) const {
    return getField(mltFromGeo(pos, Base::time_es()), SphericalPolar(pos).theta);
  }

  T m_time{0};
};

static_assert(ElectricPotentialModel<CrossTailPotential<double>, double>);
static_assert(ElectricPotentialModel<CrossTailPotential<float>, float>);

static_assert(!ElectricPotentialModel<CrossTailPotential<float>, double>);
static_assert(!ElectricPotentialModel<CrossTailPotential<double>, float>);

}
