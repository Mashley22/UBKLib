module;

#include <concepts>

export module UBKLib.field_models:electric_potentials.volland_stern;

import UBKLib.utils;

export namespace ubk {

template<std::floating_point T>
class VollandStern {
public:
  static constexpr kV<T> SURFACE_POTENTIAL = static_cast<T>(94.2);
  static constexpr kV<T> TAIL_FIELD = static_cast<T>(10);
  
  // from Maynard and Chen (1975), UBKikuity
  [[nodiscard]] kV<T>
  getField(T mlt, T theta) const {
    T E0 = 0.045 / pow(1 - 0.159 * m_kp + 0.0093 * pow(m_kp, 2), 3);
    
    T sinTheta_2 = pow(sin(theta), 2);
    return -1 * SURFACE_POTENTIAL * sinTheta_2 - E0 * sin((mlt - 12) * PI<T> / 12) / pow(sinTheta_2, 2);
  }

private:
  T m_kp{0};

};

}
