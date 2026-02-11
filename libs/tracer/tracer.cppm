module;

#include <concepts>
#include <optional>
#include <vector>
#include <span>
#include <cmath>
#include <stdexcept>
#include <cstdint>
#include <utility>

#include <UBK/macros.hpp>

export module UBKLib.tracer;

import UBKLib.utils;
import UBKLib.field_models;

template<std::floating_point T>
[[nodiscard]] static constexpr bool // anonymous namespace/ even putting this in the ubk namespace doesnt work....
oppositeSigns(T a, T b) {
  return (a * b) < 0;
}

export namespace ubk {

class BifercatingFieldLine {};

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
    if constexpr (std::is_same_v<T, float>) {
      return outterLim * 1e-7f;
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
  
  struct UBKInfos {
    Vector3<Re<T>> loc{};
    Vector3<microTesla<T>> magneticField{};
    T magneticIntensity{};
    T electricPotential{};

  };

  struct FullPointInfo {
    Vector3<Re<T>> loc{};
    Vector3<microTesla<T>> magneticField{};
    T magneticIntensity{};
    T longitudinalInvariant{};

    operator UBKInfos() const UBK_NOEXCEPT {
      return {
        .loc = loc,
        .magneticField = magneticField,
        .magneticIntensity = magneticIntensity
      };
    };
  };

  friend class FieldLineGenerator<T, FieldModel, Params>;

  [[nodiscard]] constexpr std::span<const FullPointInfo>
  points(void) const UBK_NOEXCEPT {
    return m_points;
  }

  [[nodiscard]] constexpr std::span<FullPointInfo>
  points(void) UBK_NOEXCEPT {
    return m_points;
  }

  [[nodiscard]] constexpr const T
  maxLongitudinalInvariant(void) const UBK_NOEXCEPT {
    return m_points.front().longitudinalInvariant; // Note that this should ALWAYS be subsituteable for back()
  }

  [[nodiscard]] constexpr std::array<UBKInfos, 2>
  getPointsWithK(T k_val) {  
    check(k_val != 0);
    std::array<UBKInfos, 2> retVal;
    std::size_t which = 0;
    for (std::size_t i = 0; i < m_points.size() - 1; i++) {
      if (k_val == m_points[i].longitudinalInvariant) {
        retVal[which] = m_points[i];
        if (which == 1) {
          return retVal;
        }
        which++;
      }
      
      T deltaK_1 = k_val - m_points[i].longitudinalInvariant;
      T deltaK_2 = k_val - m_points[i + 1].longitudinalInvariant;
      
      if (oppositeSigns<T>(deltaK_1, deltaK_2)) {
        Vector3<Re<T>> deltaLoc = m_points[i + 1].loc - m_points[i].loc;
        Vector3<microTesla<T>> deltaField = m_points[i + 1].magneticField - m_points[i].magneticField;
        T deltaK = m_points[i + 1].longitudinalInvariant - m_points[i].longitudinalInvariant;

        T scale = deltaK_1 / deltaK;

        retVal[which].loc = deltaLoc * scale + m_points[i].loc;
        retVal[which].magneticField = deltaField * scale + m_points[i].magneticField;
        retVal[which].magneticIntensity = retVal[which].magneticField.amp();

        retVal[which] = m_points[i];
        if (which == 1) {
          return retVal;
        }
        which++;
      }
    }
    std::unreachable();
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
    T h = 1 / fieldIntensity * std::min(Params.maxStepSize, Params.maxStepDotField / fieldIntensity);
    
    while(true) {
      Vector3<Re<T>> step = static_cast<T>(0.5) * h * (field + m_fieldModel.getField(loc + h * field));
      if constexpr (direc == FillDirection::BACKWARD) {
        step = -1 * step;
      }
      Vector3<Re<T>> newLoc = step + loc;
      if (validStep_(newLoc)) {
        return newLoc;
      }
      else if (step.ampSquared() < (Params.minStepSize() * Params.minStepSize())) {
        return std::nullopt;
      }
      h = h / Params.failRatio;
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

    microTesla<T> minMagneticIntensity = point.magneticIntensity;
    bool foundMinima = false;

    std::optional<Vector3<Re<T>>> nextLoc = takeStep_<direc>(point.loc, point.magneticField, point.magneticIntensity);

    if (!nextLoc.has_value()) {
      std::runtime_error("Couldn't even take one step!");
    }

    auto checkNotBifercating = [&](microTesla<T> newIntensity) {
      if (newIntensity > minMagneticIntensity) {
        foundMinima = true;
      }
      else {
        if (foundMinima) {
          throw BifercatingFieldLine{};
        }
        minMagneticIntensity = newIntensity;
      }
    };

    auto assignPoint = [&]() {
      point.loc = nextLoc.value();
      point.magneticField = m_fieldModel.getField(point.loc);
      point.magneticIntensity = point.magneticField.amp();
      
      checkNotBifercating(point.magneticIntensity);

      buf_<direc>().push_back(point);
    };

    assignPoint();
    
    for (std::size_t i = 0; i < Params.maxStepCount; i++) {
      nextLoc = takeStep_<direc>(point.loc, point.magneticField, point.magneticIntensity);
      if (!nextLoc.has_value()) {
        return;
      }
      
      assignPoint();
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

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
void
calculateLongitudinalInvariants(FieldLine<T, FieldModel, Params>& fieldLine) {

  auto longitudinalInvariant = [&](std::size_t startIdx, int direc) {
    T K = 0;
    std::size_t idx = startIdx;
    auto nextIdx = [&](){ return static_cast<std::size_t>(static_cast<std::int64_t>(idx) + direc); };

    while (fieldLine.points()[startIdx].magneticIntensity >
      fieldLine.points()[nextIdx()].magneticIntensity) {

      K += std::sqrt((fieldLine.points()[startIdx].magneticIntensity -
                     fieldLine.points()[nextIdx()].magneticIntensity) *
                     (fieldLine.points()[nextIdx()].loc -
                     fieldLine.points()[idx].loc).ampSquared());

      idx = nextIdx();

      if (idx == 0 || idx == fieldLine.points().size()) {
        break;
      }
    }

    return K;
  };
  
  fieldLine.points().front().longitudinalInvariant = longitudinalInvariant(0, 1);
  fieldLine.points().back().longitudinalInvariant = longitudinalInvariant(fieldLine.points().size() - 1, -1);
  
  // for a bit of certainty here
  fieldLine.points().front().longitudinalInvariant = (fieldLine.points().front().longitudinalInvariant +
    fieldLine.points().front().longitudinalInvariant) / 2;
  fieldLine.points().front().longitudinalInvariant = fieldLine.points().back().longitudinalInvariant;

  for (std::size_t i = 1; i < fieldLine.points().size() - 1; i++) {
    bool forward = fieldLine.points()[i].magneticIntensity > fieldLine.points()[i + 1].magneticIntensity;
    bool backward = fieldLine.points()[i].magneticIntensity > fieldLine.points()[i - 1].magneticIntensity;

    if (forward && backward) {
      std::terminate();
    }
    else if (!forward && !backward) {
      fieldLine.points()[i].longitudinalInvariant = 0;
    }
    else {
      if (forward) {
        fieldLine.points()[i].longitudinalInvariant = longitudinalInvariant(i, 1);
      }
      if (backward) {
        fieldLine.points()[i].longitudinalInvariant = longitudinalInvariant(i, -1);
      }
    }
  }
}

}
