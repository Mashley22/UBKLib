module;

#include <concepts>
#include <cmath>

export module ubk.field_models:electric_potentials.cross_tail;

import :traits;
import ubk.utils;

namespace ubk {

template<std::floating_point T>
class CrossTailPotential {
public:

static constexpr kV<T> ROTATION_VOLTAGE = 92.0; 
static constexpr kV<T> TAIL_POTENTIAL   = 50.0; 
static constexpr Re<T> TAIL_WIDTH       = 30.0; 

  [[nodiscard]] kV<T>
  getField(Vector3<T> pos) const {
      
    T r = std::sqrt(pos.x * pos.x + pos.y * pos.y);

    T v_rot = -ROTATION_VOLTAGE / r;

    T E_field = TAIL_POTENTIAL / TAIL_WIDTH; 
    T v_conv  = -E_field * pos.y;

    return v_rot + v_conv;
  }

};

static_assert(ElectricPotentialModel<CrossTailPotential<double>, double>);
static_assert(ElectricPotentialModel<CrossTailPotential<float>, float>);

static_assert(!ElectricPotentialModel<CrossTailPotential<float>, double>);
static_assert(!ElectricPotentialModel<CrossTailPotential<double>, float>);

}
