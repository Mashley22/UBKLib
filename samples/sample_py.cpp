#include <pybind11/pybind11.h>

#include <array>
#include <numeric>
#include <optional>

#include <cxform.h>

import UBKLib;

using T = double;

namespace py = pybind11;

namespace {

constexpr ubk::FieldLineParams<T> params = {
  .innerLim = 1.05,
  .outterLim = 15.0,
  .maxStepDotField = 1,
  .failRatio = 2,
  .maxStepSize = 0.01,
  .maxStepCount = 1000,
};

constexpr ubk::Time simTime {
  .year = 1960, 
  .month = 1,
  .day = 15,
};

struct Field {
  ubk::Igrf13<T> igrf13;
  ubk::Ts89<T> ts89;

  [[nodiscard]] ubk::Vector3<ubk::nanoTesla<T>>
  getField(ubk::Vector3<T> pos) const {
    return ts89.getField(pos) + igrf13.getField(pos);
  }

};

ubk::FieldLineGenerator<T, Field, params> generator;

static std::optional<ubk::Vector3<ubk::Re<T>>> 
gsmToGeo(ubk::Vector3<ubk::Re<T>> gsm) {
  std::array<T, 3> geo;
  if (cxform2(GSM, GEO, static_cast<double>(ubk::timeToEs(simTime)), gsm.toArr().data(), geo.data()) != 0) {
    return std::nullopt;
  }
  return std::make_optional(ubk::Vector3<ubk::Re<T>>::fromArr(geo));
}

}

T
calculateB(T x, T y, T k) {
  Field field;
  field.igrf13.setTime(simTime);
  field.ts89.dipole_tilt() = field.igrf13.dipole_tilt();
  field.ts89.setTime(simTime);

  generator.assignModel(field);

  ubk::Vector3<ubk::Re<T>> gsm{x, y, 0};
  auto geo = gsmToGeo(gsm);

  if (!geo.has_value()) {
    return std::numeric_limits<T>::quiet_NaN();
  }
  
  ubk::FieldLine<T, Field, params> fieldLine;
  try {
    fieldLine = generator.generateFieldLine(geo.value());
  } catch(...) {
    return std::numeric_limits<T>::quiet_NaN();
  }
 
  try {
    ubk::calculateLongitudinalInvariants(fieldLine);
  } catch(...) {
    return std::numeric_limits<T>::quiet_NaN();
  }

  auto val = fieldLine.getPointsWithK(k);

  return val[1].magneticIntensity;
}

PYBIND11_MODULE(UBKLibpp, m) {
  m.def("calculateB", &calculateB, 
    "Calculates B for a K value on a field line with a given equatorial point");
}
