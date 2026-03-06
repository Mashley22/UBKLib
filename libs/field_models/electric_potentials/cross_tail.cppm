module;

#include <concepts>

export module UBKLib.field_models:electric_potentials.cross_tail;

import :traits;
import UBKLib.utils;

export namespace ubk {

template<std::floating_point T>
class CrossTailPotential {
public:

static constexpr kV<T> SURFACE_POTENTIAL = static_cast<T>(92.0); 
static constexpr T TAIL_FIELD = static_cast<T>(10.0); 

  [[nodiscard]] kV<T>
  getField(T mlt, T theta) const {
    T sinTheta_2 = sin(theta);

    return -1 * SURFACE_POTENTIAL * sinTheta_2 - sin((mlt - 12) * PI<T> / 12) * TAIL_FIELD / sinTheta_2;
  }

  [[nodiscard]] kV<T>
  getField(const Vector3<T> pos) const {
    return pos();
  }

  T m_time{0};
};

static_assert(ElectricPotentialModel<CrossTailPotential<double>, double>);
static_assert(ElectricPotentialModel<CrossTailPotential<float>, float>);

static_assert(!ElectricPotentialModel<CrossTailPotential<float>, double>);
static_assert(!ElectricPotentialModel<CrossTailPotential<double>, float>);

}
