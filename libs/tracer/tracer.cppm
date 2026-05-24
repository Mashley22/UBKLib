module;

#include <concepts>
#include <optional>
#include <vector>
#include <span>
#include <cmath>
#include <stdexcept>
#include <cstdint>
#include <utility>
#include <fstream>

#include <iostream>

#include <UBK/macros.hpp>

export module UBKLib.tracer;

import UBKLib.utils;
import UBKLib.field_models;

template<std::floating_point T>
[[nodiscard]] static constexpr bool
oppositeSigns(T a, T b) {
  return (a * b) <= 0;
}

export namespace ubk {

class BifurcatingFieldLine {};

class NoMinimaFound {};

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

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
class FieldLineGenerator;

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
class FieldLine {
public:

  struct PointInfo {
    Vector3<Re<T>> loc{};
    nanoTesla<T> magneticIntensity{};
    T longitudinalInvariant{};
  };

  friend class FieldLineGenerator<T, FieldModel, Params>;

  [[nodiscard]] constexpr std::span<const PointInfo>
  points(void) const UBK_NOEXCEPT {
    return m_points;
  }

  [[nodiscard]] constexpr std::span<PointInfo>
  points(void) UBK_NOEXCEPT {
    return m_points;
  }

  [[nodiscard]] constexpr const T
  maxLongitudinalInvariant(void) const UBK_NOEXCEPT {
    return m_points.front().longitudinalInvariant; // Note that this should ALWAYS be subsituteable for back()
  }
  
  /**
   *@brief get the two mirror points with a value of k larger equal to the value requested
   *
   *@note undefined if k_val > max_k along the field line
   *
  */
  [[nodiscard]] constexpr std::array<PointInfo, 2>
  getPointsWithK(T k_val) {  
    check(k_val != 0);
    std::array<PointInfo, 2> retVal;
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
        T deltaK = m_points[i + 1].longitudinalInvariant - m_points[i].longitudinalInvariant;
        nanoTesla<T> deltaMagneticIntensity = m_points[i + 1].magneticIntensity - m_points[i].magneticIntensity;

        if (deltaK == 0) {
          retVal[which] = m_points[i];
        }
        else {
          T scale = deltaK_1 / deltaK;

          retVal[which].loc = deltaLoc * scale + m_points[i].loc;
          retVal[which].magneticIntensity = deltaMagneticIntensity * scale + m_points[i].magneticIntensity;
        }

        if (which == 1) {
          return retVal;
        }
        which++;
      }
    }
    std::unreachable();
  }

  [[nodiscard]] constexpr PointInfo
  getMinima(void) {
    nanoTesla<T> minIntensity = m_points[0].magneticIntensity;
    for (std::size_t i = 1; i < m_points.size() - 1; i++) {

      if (minIntensity < m_points[i].magneticIntensity) {
        
        // Note that we can ignore the determinant = 0 check;
        auto invert3x3 = [](const T m[3][3]) {
          std::array<std::array<T, 3>, 3> inv;
          T det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
                  m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                  m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

          T invDet = 1.0 / det;

          inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * invDet;
          inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet;
          inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet;

          inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet;
          inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet;
          inv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invDet;

          inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * invDet;
          inv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invDet;
          inv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invDet;

          return inv;
        };

        auto getQuadraticCoefficients = [&](T s_forward, T s_backward, const std::array<T, 3>& B_vals) {
          std::array<T, 3> coeffs;
          T forwardMat[3][3] = {
            {pow(s_backward, 2), s_backward, 1},
            {0, 0, 1},
            {pow(s_forward, 2), s_forward, 1}
          };

          std::array<std::array<T, 3>, 3> inverseMat = invert3x3(forwardMat);

          for (std::size_t i = 0; i < 3; i++) {
            coeffs[i] = 0;
            for (std::size_t j = 0; j < 3; j++) {
              coeffs[i] += inverseMat[i][j] * B_vals[i];
            }
          }

          return coeffs;

        };

        T s_backward = -1 * (m_points[i].loc - m_points[i - 1].loc).amp();
        T s_forward = (m_points[i].loc - m_points[i + 1].loc).amp();

        std::array<T, 3> quadraticCoeffs =
          getQuadraticCoefficients(s_forward, s_backward, 
                                   std::array<T, 3>{
                                   m_points[i - 1].magneticIntensity,
                                   m_points[i].magneticIntensity,
                                   m_points[i + 1].magneticIntensity
                                   });

        T s_min = -1 * quadraticCoeffs[1] / (2 * quadraticCoeffs[0]);
        
        Vector3<Re<T>> secondPointLoc;
        if (s_min < 0) {
          secondPointLoc = m_points[i - 1].loc;
        }
        else {
          secondPointLoc = m_points[i + 1].loc;
        }

        Vector3<Re<T>> deltaLoc = secondPointLoc - m_points[i].loc;
        Vector3<Re<T>> loc_min = m_points[i].loc - deltaLoc * s_min / deltaLoc.amp();

        return {
          .loc = loc_min,
        };
      }
      else {
        minIntensity = m_points[i].magneticIntensity;
      }
    }
    throw NoMinimaFound{};
  }
  
private:
  std::vector<PointInfo> m_points;
};

template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
class FieldLineGenerator {
public:
  using FieldLinePoint = FieldLine<T, FieldModel, Params>::PointInfo;

  [[nodiscard]] FieldLine<T, FieldModel, Params>
  generateFieldLine(Vector3<Re<T>> startPoint) {
    clearAll_();

    trace_<FillDirection::FORWARD>(startPoint);
    trace_<FillDirection::BACKWARD>(startPoint);

    FieldLine<T, FieldModel, Params> fieldLine;
    std::size_t backIdx{0};
    std::size_t frontIdx{0};

    if (m_backward.back().magneticIntensity < m_forward.back().magneticIntensity) {
      backIdx = trimForward_();
    }
    else {
      frontIdx = trimBackward_();
    }
    
    if (!m_backward.empty()) {
      for (auto it = m_backward.rbegin();
      it != std::next(m_backward.rend(), -1 * static_cast<std::ptrdiff_t>(backIdx));
      it++) {
        fieldLine.m_points.push_back(*it); //could pop it
      }
    }

    {
      FieldLinePoint point = {
        .loc = startPoint,
        .magneticIntensity = m_fieldModel.getField(startPoint).amp()
      };

      fieldLine.m_points.push_back(point);
    }

    for (auto it = std::next(m_forward.begin(), static_cast<std::ptrdiff_t>(frontIdx));
    it != m_forward.end();
    it++) {
      fieldLine.m_points.push_back(*it); //could pop it
    }

    return fieldLine;
  }

  void
  assignModel(const FieldModel& fieldModel) {
    m_fieldModel = fieldModel;
  }

  [[nodiscard]] const FieldModel& 
  fieldModel(void) const UBK_NOEXCEPT {
    return m_fieldModel;
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
  takeStep_(Vector3<Re<T>> loc, Vector3<nanoTesla<T>> field, Re<T> prevStepSize = Params.maxStepSize) {
    T h = prevStepSize;
    field = field.normalised();
    if constexpr (direc == FillDirection::BACKWARD) {
      h = -1 * h;
    }
    
    while(true) {
      Vector3<Re<T>> step = h * (field + m_fieldModel.getField(loc + h * field)).normalised();
      Vector3<Re<T>> newLoc = step + loc;
      if (validStep_(newLoc)) {
        return step;
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
  trace_(Vector3<Re<T>> starting) {
    auto field = m_fieldModel.getField(starting);
    FieldLinePoint point = {
      .loc = starting,
      .magneticIntensity = field.amp(),
    };

    std::optional<Vector3<Re<T>>> nextStep = takeStep_<direc>(point.loc, field);

    if (!nextStep.has_value()) {
      throw std::runtime_error("Couldn't even take one step!");
    }

    auto assignPoint = [&](std::optional<Vector3<Re<T>>> step) {
      point.loc = step.value() + point.loc;
      field = m_fieldModel.getField(point.loc);
      point.magneticIntensity = field.amp();
      
      buf_<direc>().push_back(point);
    };

    assignPoint(nextStep);
    
    for (std::size_t i = 0; i < Params.maxStepCount; i++) {
      nextStep = takeStep_<direc>(point.loc, field, nextStep.value().amp());
      if (!nextStep.has_value()) {
        return;
      }
      
      assignPoint(nextStep);
    }
  }

  void
  clearAll_(void) {
    m_forward.resize(0);
    m_backward.resize(0);
  }
    
  /**
   *@brief pops off elements of the m_forward vector until the end points have the same
   *       magnetic intensity, may start trimming from the front of m_backward
   *
   *@returns how many elements off the front of the m_backward to trim off
  */
  [[nodiscard]] std::size_t
  trimForward_(void) {
    nanoTesla<T> targetIntensity = m_backward.back().magneticIntensity;
    FieldLinePoint prevPoint = m_forward.back();
    m_forward.pop_back();

    while(targetIntensity < m_forward.back().magneticIntensity) {
      prevPoint = m_forward.back();
      m_forward.pop_back();
      if (m_forward.empty()) {
        break;
      }
    }
  
    if (m_forward.empty()) {
      std::size_t backIdx = 0;
      while(m_backward[backIdx].magneticIntensity < targetIntensity) {
        backIdx++;
      }
      return backIdx;
    }
    else {
      FieldLinePoint res = {
        .loc = (targetIntensity - m_forward.back().magneticIntensity) /
          (prevPoint.magneticIntensity - m_forward.back().magneticIntensity) *
          (m_forward.back().loc - prevPoint.loc) +
          m_forward.back().loc,
        .magneticIntensity = targetIntensity
      };
      m_forward.push_back(res);
      return 0;
    }
  }
  
  /**
   *@brief 
  */
  [[nodiscard]] std::size_t
  trimBackward_(void) {
    nanoTesla<T> targetIntensity = m_forward.back().magneticIntensity;
    FieldLinePoint prevPoint = m_backward.back();
    m_backward.pop_back();

    while(targetIntensity < m_backward.back().magneticIntensity) {
      prevPoint = m_backward.back();
      m_backward.pop_back();
      if (m_backward.empty()) {
        break;
      }
    }

    if (m_backward.empty()) {
      std::size_t frontIdx = 0;
      while(m_forward[frontIdx].magneticIntensity < targetIntensity) {
        frontIdx++;
      }

      return frontIdx;
    }
    else {
      FieldLinePoint res = {
        .loc = (targetIntensity - m_backward.back().magneticIntensity) /
          (prevPoint.magneticIntensity - m_backward.back().magneticIntensity) *
          (m_backward.back().loc - prevPoint.loc) +
          m_backward.back().loc,
        .magneticIntensity = targetIntensity
      };
      m_backward.push_back(res);
      return 0;
    }
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

      if (idx == 1 || idx == fieldLine.points().size() - 1) {
        break;
      }
    }

    return K;
  };

  {
    T K = 0;
    constexpr std::size_t startIdx = 0;
    std::size_t idx = 0;
    bool minimaFound = false;
    T minIntensity = fieldLine.points()[startIdx].magneticIntensity;

    auto nextIdx = [&](){ return static_cast<std::size_t>(static_cast<std::int64_t>(idx) + 1); };

    while (fieldLine.points()[startIdx].magneticIntensity >
      fieldLine.points()[nextIdx()].magneticIntensity) {
      
      if (fieldLine.points()[nextIdx()].magneticIntensity < minIntensity) {
        if (minimaFound) { throw BifurcatingFieldLine{}; }
        minimaFound = true;
      }
      else {
        minIntensity = fieldLine.points()[nextIdx()].magneticIntensity;
      }

      K += std::sqrt((fieldLine.points()[startIdx].magneticIntensity -
                     fieldLine.points()[nextIdx()].magneticIntensity) *
                     (fieldLine.points()[nextIdx()].loc -
                     fieldLine.points()[idx].loc).ampSquared());

      idx = nextIdx();

      if (idx == 1 || idx == fieldLine.points().size() - 1) {
        break;
      }
    }
      
    if (!minimaFound) {
      throw NoMinimaFound{};
    }
    fieldLine.points().front().longitudinalInvariant = K;
  }

  fieldLine.points().back().longitudinalInvariant = longitudinalInvariant(fieldLine.points().size() - 1, -1);
  
  // for a bit of certainty here
  fieldLine.points().front().longitudinalInvariant = (fieldLine.points().front().longitudinalInvariant +
    fieldLine.points().back().longitudinalInvariant) / 2;
  fieldLine.points().front().longitudinalInvariant = fieldLine.points().back().longitudinalInvariant;

  for (std::size_t i = 1; i < fieldLine.points().size() - 1; i++) {
    auto dB_forward = fieldLine.points()[i + 1].magneticIntensity - fieldLine.points()[i].magneticIntensity;
    auto dB_backward = fieldLine.points()[i].magneticIntensity - fieldLine.points()[i - 1].magneticIntensity;

    bool turningPointRegion = oppositeSigns(dB_backward, dB_forward);

    if (!turningPointRegion) {
      if (dB_forward > 0) {
        fieldLine.points()[i].longitudinalInvariant = longitudinalInvariant(i, -1);
      }
      else {
        fieldLine.points()[i].longitudinalInvariant = longitudinalInvariant(i, 1);
      }
    }
    else {
      fieldLine.points()[i].longitudinalInvariant = 0;
    }
  }
}


template<std::floating_point T, class FieldModel, FieldLineParams<T> Params>
requires MagneticFieldModel<FieldModel, T>
void
saveFieldLine(const FieldLine<T, FieldModel, Params>& fieldLine, const std::string& filename) {
  std::ofstream outFile(filename, std::ios::binary);

  if (outFile.is_open()) {
    outFile.write(reinterpret_cast<const char*>(fieldLine.points().data()), static_cast<std::streamsize>(fieldLine.points().size() * sizeof(typename FieldLine<T, FieldModel, Params>::PointInfo)));
  }
}

}
