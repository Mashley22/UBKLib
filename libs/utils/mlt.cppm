module;

#include <array>
#include <concepts>
#include <cmath>
#include <numbers>

#include <cxform.h>

export module UBKLib.utils:mlt;

import :math;
import :units;

export namespace ubk {

template<std::floating_point T>
[[nodiscard]] T
mltFromGsm(const Vector3<Re<T>> gsm) {
  return std::fmod((std::atan2(gsm.y, gsm.x) * 12 / std::numbers::pi) + 12, 24);
}

template<std::floating_point T>
[[nodiscard]] T
mltFromGeo(const Vector3<Re<T>> geo, double t) {
  Vector3<Re<T>> gsm;
  {
    std::array<Re<T>, 3> tempGsm;
    cxform2(GEO, GSM, t, geo.toArr().data(), tempGsm.data());
    gsm = Vector3<Re<T>>::fromArr(tempGsm);
  }
  return mltFromGsm(gsm);
}

}
