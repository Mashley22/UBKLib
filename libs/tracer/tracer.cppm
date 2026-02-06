module;

#include <concepts>
#include <optional>
#include <vector>
#include <span>
#include <cmath>
#include <stdexcept>

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
  T failRatio = 2;
  T maxStepSize = 0.01;
  std::size_t maxStepCount = 10000;

  constexpr T minStepSize(void) const UBK_NOEXCEPT {
    if constexpr (std::is_same_v<T, double>) {
      return outterLim * 1e-15;
    }
    if constexpr (std::is_same_v<T, double>) {
      return outterLim * 1e-7;
    }
  }
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

  struct FullPointInfo {
    Vector3<Re<T>> loc{};
    Vector3<microTesla<T>> magneticField{};
    T magneticIntensity{};
    T longitudinalInvariant{};
  };

  friend class FieldLineGenerator<T, FieldModel, Params>;

  [[nodiscard]] constexpr std::span<const FullPointInfo>
  points(void) const UBK_NOEXCEPT {
    return m_points;
  }
  
private:
  std::vector<FullPointInfo> m_points;
};

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
class FieldLineGenerator {
public:
  using FieldLinePoint = FieldLine<T, FieldModel, Params>::FullPointInfo;

  [[nodiscard]] FieldLine<T, FieldModel, Params>
  generateFieldLine(Vector3<Re<T>> startPoint) {
    clearAll_();
    fill_<FillDirection::FORWARD>(startPoint);
    fill_<FillDirection::BACKWARD>(startPoint);
    
    FieldLine<T, FieldModel, Params> fieldLine;

    if (m_backward.back().magneticIntensity < m_forward.back().magneticIntensity) {
      trimForward_();
    }
    else {
      trimBackward_();
    }


    for (auto it = m_backward.rbegin(); it != m_backward.rend(); it++) {
      fieldLine.m_points.push_back(*it); //could pop it
    }

    {
      FieldLinePoint point = {
        .loc = startPoint,
        .magneticField = m_fieldModel.getField(startPoint)
      };
      point.magneticIntensity = point.magneticField.amp();

      fieldLine.m_points.push_back(point);
    }

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
  std::vector<FieldLinePoint> m_forward;
  std::vector<FieldLinePoint> m_backward;
  FieldModel m_fieldModel;

  enum class FillDirection {
    FORWARD,
    BACKWARD
  };
  
  [[nodiscard]] bool
  validStep_(Vector3<Re<T>> loc) {
    T distSquared = loc.ampSquared();
    return (Params.innerLim * Params.innerLim < distSquared) &&
           (Params.outterLim * Params.outterLim > distSquared);
  }
  
  template<FillDirection direc>
  [[nodiscard]] std::optional<Vector3<Re<T>>>
  takeStep_(Vector3<Re<T>> loc, Vector3<microTesla<T>> field, microTesla<T> fieldIntensity) {
    Vector3<Re<T>> step = field / fieldIntensity * std::min(Params.maxStepSize, Params.maxStepDotField / fieldIntensity);
    
    if constexpr (direc == FillDirection::BACKWARD) {
      field = -1 * field;
      step = -1 * step;
    }

    while(true) {
      Vector3<Re<T>> newLoc = step + loc;
      if (validStep_(newLoc)) {
        return newLoc;
      }
      else if (step.ampSquared() < (Params.minStepSize() * Params.minStepSize())) {
        return std::nullopt;
      }
      step = step / Params.failRatio;
    }
  }
    
  template<FillDirection direc>
  [[nodiscard]]
  std::vector<FieldLinePoint>&
  buf_(void) UBK_NOEXCEPT {
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
    
    FieldLinePoint point = {
      .loc = starting,
      .magneticField = m_fieldModel.getField(point.loc),
    };
    point.magneticIntensity = point.magneticField.amp();

    std::optional<Vector3<Re<T>>> nextLoc = takeStep_<direc>(point.loc, point.magneticField, point.magneticIntensity);

    if (!nextLoc.has_value()) {
      std::runtime_error("Couldn't even take one step!");
    }
    
    point.loc = nextLoc.value();
    point.magneticField = m_fieldModel.getField(point.loc);
    point.magneticIntensity = point.magneticField.amp();
    buf_<direc>().push_back(point);

    for (std::size_t i = 0; i < Params.maxStepCount; i++) {
      nextLoc = takeStep_<direc>(point.loc, point.magneticField, point.magneticIntensity);
      if (!nextLoc.has_value()) {
        return;
      }

      point.loc = nextLoc.value();
      point.magneticField = m_fieldModel.getField(point.loc);
      point.magneticIntensity = point.magneticField.amp();
      buf_<direc>().push_back(point);
    }
  }

  void
  clearAll_(void) {
    m_forward.resize(0);
    m_backward.resize(0);
  }

  void
  trimForward_(void) {
    microTesla<T> targetIntensity = m_backward.back().magneticIntensity;
    FieldLinePoint prevPoint = m_forward.back();
    m_forward.pop_back();

    while(targetIntensity < m_forward.back().magneticIntensity) {
      prevPoint = m_forward.back();
      m_forward.pop_back();
    }

    FieldLinePoint res = {
      .loc = (targetIntensity - m_forward.back().magneticIntensity) /
        (prevPoint.magneticIntensity - m_forward.back().magneticIntensity) *
        (m_forward.back().loc - prevPoint.loc) +
        m_forward.back().loc,
      .magneticIntensity = targetIntensity
    };
    m_forward.push_back(res);
  }

  void
  trimBackward_(void) {
    microTesla<T> targetIntensity = m_forward.back().magneticIntensity;
    FieldLinePoint prevPoint = m_backward.back();
    m_backward.pop_back();

    while(targetIntensity < m_backward.back().magneticIntensity) {
      prevPoint = m_backward.back();
      m_backward.pop_back();
    }

    FieldLinePoint res = {
      .loc = (targetIntensity - m_backward.back().magneticIntensity) /
        (prevPoint.magneticIntensity - m_backward.back().magneticIntensity) *
        (m_backward.back().loc - prevPoint.loc) +
        m_backward.back().loc,
      .magneticIntensity = targetIntensity
    };
    m_backward.push_back(res);
  }
};

}
