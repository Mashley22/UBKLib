module;

#include <concepts>
#include <cmath>

export module ubk.field_models:electric_potentials.cross_tail;

import ubk.utils;

namespace ubk {

template<std::floating_point T>
class CrossTailPotential {
public:

static constexpr T ROTATION_VOLTAGE_kV = 92.0; 
static constexpr T TAIL_POTENTIAL_kV   = 50.0; 
static constexpr T TAIL_WIDTH_Re       = 30.0; 

  [[nodiscard]] T 
  getField(Vector3<T> pos) const {
      
      
      T r = std::sqrt(pos.x * pos.x + pos.y * pos.y);

      T v_rot = -ROTATION_VOLTAGE_kV / r;

      T E_field = TAIL_POTENTIAL_kV / TAIL_WIDTH_Re; 
      T v_conv  = -E_field * pos.y;

      return v_rot + v_conv;
  }

};

}
