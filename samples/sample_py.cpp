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
  .year = 2010, 
  .month = 2,
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
smToGeo(ubk::Vector3<ubk::Re<T>> sm) {
  std::array<T, 3> geo;
  if (cxform2(SM, GEO, static_cast<double>(ubk::timeToEs(simTime)), sm.toArr().data(), geo.data()) != 0) {
    return std::nullopt;
  }
  return std::make_optional(ubk::Vector3<ubk::Re<T>>::fromArr(geo));
}

}

py::list
calculateB(T x, T y, py::list k_vals) {
  Field field;
  field.igrf13.setTime(simTime);
  field.ts89.dipole_tilt() = field.igrf13.dipole_tilt();
  field.ts89.setTime(simTime);

  generator.assignModel(field);

  ubk::Vector3<ubk::Re<T>> sm{x, y, 0};
  auto geo = smToGeo(sm);

  py::list retVal;

  const auto& invalid_field_line = [&]() {
    for (std::size_t i = 0; i < k_vals.size(); i++) {
      retVal.append(std::numeric_limits<T>::quiet_NaN());
    }
    return retVal;
  };

  if (!geo.has_value()) {
    return invalid_field_line();
  }
  
  ubk::FieldLine<T, Field, params> fieldLine;
  try {
    fieldLine = generator.generateFieldLine(geo.value());
  } catch(...) {
    return invalid_field_line();
  }
 
  try {
    ubk::calculateLongitudinalInvariants(fieldLine);
  } catch(...) {
    return invalid_field_line();
  }

  for (const auto& k : k_vals) {
    T k_value = k.cast<T>();
    if (k_value > fieldLine.maxLongitudinalInvariant()) {
      retVal.append(std::numeric_limits<T>::quiet_NaN());
    }
    else {
      retVal.append(fieldLine.getPointsWithK(k.cast<T>())[1].magneticIntensity);
    }
  }
  
  return retVal;
}

PYBIND11_MODULE(UBKLibpp, m) {
  m.def("calculateB", &calculateB, 
    "Calculates B for a K value(s) on a field line with a given equatorial point");
}
