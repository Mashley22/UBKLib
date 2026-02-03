module;

#include <concepts>
#include <vector>
#include <span>
#include <cmath>

#include <UBK/macros.hpp>

export module ubk.tracer;

import ubk.utils;
import ubk.field_models;

export namespace ubk {

template<std::floating_point T>
struct PointsPair {
  Vector3<Re<T>> begin, end;

  Vector3<Re<T>> diff(void) const UBK_NOEXCEPT {
    return begin - end;
  }
};

template<std::floating_point T>
using MirrorPoints = PointsPair<T>;

// here to make field model stuff simpler
template<std::floating_point T, class FieldModel>
requires MagneticFieldModel<FieldModel, T>
[[nodiscard]] T
integrationStep(PointsPair<T> points, const FieldModel& field, microTesla<T> mirrorPointMagneticIntensity) {
  microTesla<T> localMagneticIntensity = field.getField(points.begin).amp();
  T distSquared = points.diff().ampSquared();
  return std::sqrt((localMagneticIntensity - mirrorPointMagneticIntensity) * distSquared);
}

template<std::floating_point T>
class FieldLine {

  [[nodiscard]] std::span<const Vector3<Re<T>>>
  points(void) const UBK_NOEXCEPT {
    return m_points;
  }
  
  [[nodiscard]] MirrorPoints<T>
  mirrorPoints(void) const UBK_NOEXCEPT {
    return {.begin = m_points.front(), .end = m_points.back()};
  }
  
  template<class FieldModel>
  requires MagneticFieldModel<FieldModel, T>
  [[nodiscard]] T
  maxModifiedLongitudinalInvariant(const FieldModel& field) {
    T retVal = 0;

    for (std::size_t i = 0; i < m_points.size() - 1; i++) {
      retVal += integrationStep({.begin = m_points[i], .end = m_points[i + 1]},
                                field);
    }

    return retVal;
  }

private:
  microTesla<T> m_mirrorPointMagneticIntensity;
  std::vector<Vector3<Re<T>>> m_points;

  template<class FieldModel>
  requires MagneticFieldModel<FieldModel, T>
  [[nodiscard]] T
  integrationStep_(PointsPair<T> points, const FieldModel& field) {
    microTesla<T> localMagneticIntensity = field.getField(points.begin).amp();
    T distSquared = points.diff().ampSquared();
    return std::sqrt((localMagneticIntensity - m_mirrorPointMagneticIntensity) * distSquared);
  }
};

}
