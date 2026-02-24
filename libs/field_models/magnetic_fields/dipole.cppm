module;

#include <concepts>

export module UBKLib.field_models:magnetic_fields.dipole;

import :traits;
import UBKLib.utils;

export namespace ubk {

template<std::floating_point T>
class Dipole {
public:
  static constexpr nanoTesla<T> B_0 = static_cast<T>(31200);

  [[nodiscard]] Vector3<nanoTesla<T>>
  getField(Vector3<T> pos) const {
    T r_5 = pos.amp() * pos.ampSquared() * pos.ampSquared();
    return Vector3<nanoTesla<T>> {
      .x = -3 * B_0 * pos.x * pos.z / r_5,
      .y = -3 * B_0 * pos.y * pos.z / r_5,
      .z = B_0 * (pos.ampSquared() - 3 * pos.z * pos.z) / r_5
    };
  }
};

static_assert(MagneticFieldModel<Dipole<double>, double>);
static_assert(MagneticFieldModel<Dipole<float>, float>);

static_assert(!MagneticFieldModel<Dipole<float>, double>);
static_assert(!MagneticFieldModel<Dipole<double>, float>);

}
