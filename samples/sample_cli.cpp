#include <print>
#include <string>
#include <vector>
#include <optional>

#include <cxform.h>

import UBKLib;

#define STARTPOINT_X_ARG_IDX 1
#define STARTPOINT_Y_ARG_IDX 2
#define STARTPOINT_Z_ARG_IDX 3

using T = double;

std::vector<T> k_vals;

constexpr ubk::FieldLineParams<T> params = {
  .innerLim = 1.05,
  .outterLim = 15.0,
  .maxStepDotField = 0.01,
  .failRatio = 2,
  .maxStepSize = 0.01,
  .maxStepCount = 10000,
};

ubk::Time simTime {
  .year = 1960, 
  .month = 1,
  .day = 15,
};

struct Field {
  ubk::Dipole<T> dipole;

  [[nodiscard]] ubk::Vector3<ubk::nanoTesla<T>>
  getField(ubk::Vector3<T> pos) const {
    return dipole.getField(pos);
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

int main(int argc, char** argv) {
  if (argc < (STARTPOINT_Y_ARG_IDX + 1)) {
    std::println("Please provide atleast the x, y (gsm) components of the start point for field line tracing!");
    return 1;
  }
  
  ubk::Vector3<ubk::Re<T>> startPoint_gsm = {
    .x = std::stod(argv[STARTPOINT_X_ARG_IDX]),
    .y = std::stod(argv[STARTPOINT_Y_ARG_IDX]),
  };
  
  if (argc == STARTPOINT_Y_ARG_IDX + 1) {
    startPoint_gsm.z = 0;
    k_vals.push_back(0);
  }
  else {
    startPoint_gsm.z = std::stod(argv[STARTPOINT_Z_ARG_IDX]);
  }

  ubk::Vector3<ubk::Re<T>> startPoint;
 
  { 
    auto val = gsmToGeo(startPoint_gsm);
    if (!val.has_value()) {
      std::println("Error in coordinate transform of gsm to geo of start point");
      return 2;
    }
    startPoint = val.value();
  }
  
  for (int i = STARTPOINT_Z_ARG_IDX + 1; i < argc; i++) {
    T val = std::stod(argv[i]);
    if (val < 0) {
      std::println("Supplied K = {}. K must be greater than 0!", val);
      return 2; 
    }
    k_vals.push_back(val);
  }

  if (k_vals.empty()) {
    k_vals.push_back(0);
  }

  auto fieldLine = generator.generateFieldLine(startPoint);
  try { 
    ubk::calculateLongitudinalInvariants(fieldLine);
  } catch(ubk::BifercatingFieldLine& e) {
    (void)e;
    std::println("Bifercating field line!");
    return 3;
  }

  for (const T k : k_vals) {
    if (k != 0) {
      std::println("{}", fieldLine.getPointsWithK(k)[0].magneticIntensity);
    }
    else {
      std::println("{}", generator.fieldModel().getField(fieldLine.getMinima().loc).amp());
    }
  }
}
