module;

#include <concepts>
#include <optional>
#include <vector>
#include <span>
#include <cmath>

#include <UBK/macros.hpp>

export module UBKLib.tracer;

import UBKLib.utils;
import UBKLib.field_models;

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

template<std::floating_point T>
struct FieldLineParams {
  T innerLim = 1.05;
  T outterLim = 15;
  T maxStepDotField = 0.01;
  T failRatio = 1.5;
  std::size_t maxTries = 100;
};

// here to make field model stuff simpler
template<std::floating_point T, class FieldModel>
requires MagneticFieldModel<FieldModel, T>
[[nodiscard]] T
integrationStep(PointsPair<T> points, const FieldModel& field, microTesla<T> mirrorPointMagneticIntensity) {
  microTesla<T> localMagneticIntensity = field.getField(points.begin).amp();
  T distSquared = points.diff().ampSquared();
  return std::sqrt((localMagneticIntensity - mirrorPointMagneticIntensity) * distSquared);
}

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
class FieldLineGenerator;

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
class FieldLine {
public:

  friend class FieldLineGenerator<T, FieldModel, Params>;

  [[nodiscard]] std::span<const Vector3<Re<T>>>
  points(void) const UBK_NOEXCEPT {
    return m_points;
  }
  
  [[nodiscard]] MirrorPoints<T>
  mirrorPoints(void) const UBK_NOEXCEPT {
    return {.begin = m_points.front(), .end = m_points.back()};
  }
  
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

  [[nodiscard]] T
  integrationStep_(PointsPair<T> points, const FieldModel& field) {
    microTesla<T> localMagneticIntensity = field.getField(points.begin).amp();
    T distSquared = points.diff().ampSquared();
    return std::sqrt((localMagneticIntensity - m_mirrorPointMagneticIntensity) * distSquared);
  }
};

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
class FieldLineGenerator {
public:
  [[nodiscard]] FieldLine<T, FieldModel, Params>
  generateFieldLine(Vector3<Re<T>> startPoint) {
    clearAll_();
    fill_<FillDirection::FORWARD>(startPoint);
    fill_<FillDirection::BACKWARD>(startPoint);
    
    FieldLine<T, FieldModel, Params> fieldLine;

    fieldLine.m_mirrorPointMagneticIntensity = m_fieldModel.getField(m_backward.back()).amp();

    for (auto it = m_backward.rbegin(); it != m_backward.rend(); it++) {
      fieldLine.m_points.push_back(*it); //could pop it
    }
    fieldLine.m_points.push_back(startPoint);

    for (auto it = m_forward.begin(); it != m_forward.end(); it++) {
      fieldLine.m_points.push_back(*it); //could pop it
    }

    return fieldLine;
  }

  void
  assignModel(const FieldModel& fieldModel) {
    m_fieldModel = fieldModel;
  }

private:
  std::vector<Vector3<Re<T>>> m_forward;
  std::vector<Vector3<Re<T>>> m_backward;
  FieldModel m_fieldModel;

  enum class FillDirection {
    FORWARD,
    BACKWARD
  };
  
  [[nodiscard]] bool
  validStep_(Vector3<Re<T>> loc, Vector3<Re<T>> step, Vector3<microTesla<T>> field) {
    T distSquared = loc.ampSquared();
    return (step.dot(field) < Params.maxStepDotField) &&
           (Params.innerLim * Params.innerLim < distSquared) &&
           (Params.outterLim * Params.outterLim > distSquared);
  }
  
  template<FillDirection direc>
  [[nodiscard]] std::optional<Vector3<Re<T>>>
  takeStep_(Vector3<Re<T>> loc) {
    Vector3<microTesla<T>> field = m_fieldModel.getField(loc);
    Vector3<Re<T>> step = field.normalised();
    
    if constexpr (direc == FillDirection::BACKWARD) {
      field = -1 * field;
      step = -1 * step;
    }

    for (std::size_t i = 0; i < Params.maxTries; i++) {
      Vector3<Re<T>> newLoc = step + loc;
      if (validStep_(newLoc, step, field)) {
        return newLoc;
      }
      step = step / Params.failRatio;
    }
    return std::nullopt;
  }
    
  template<FillDirection direc>
  [[nodiscard]]
  std::vector<Vector3<Re<T>>>& buf_(void) UBK_NOEXCEPT {
    if constexpr (direc == FillDirection::FORWARD) {
      return m_forward;
    }
    else if (direc == FillDirection::BACKWARD) {
      return m_backward;
    }
    else {
      throw;
    }
  }
  
  template<FillDirection direc>
  void 
  fill_(Vector3<Re<T>> starting) {

    std::optional<Vector3<Re<T>>> nextLoc = takeStep_<direc>(starting);

    if (!nextLoc.has_value()) {
      return;
      // throw;
    }
    
    buf_<direc>().push_back(nextLoc.value());
    while(true) {
      nextLoc = takeStep_<direc>(nextLoc.value());
      if (!nextLoc.has_value()) {
        return;
      }
      buf_<direc>().push_back(nextLoc.value());
    }
  }

  void
  clearAll_(void) {
    m_forward.resize(0);
    m_backward.resize(0);
  }
};

}
