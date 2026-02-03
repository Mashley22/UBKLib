module;

#include <concepts>

export module UBKLib.field_models:traits;

import UBKLib.utils;

export namespace ubk {

template<typename Field, typename inFP_t = double, typename outFP_t = inFP_t>
concept VectorField =
  std::floating_point<inFP_t> &&
  std::floating_point<outFP_t> &&
  std::is_default_constructible_v<Field> &&
  requires(const Field model, Vector3<inFP_t> pos) {
    { model.getField(pos) } -> std::convertible_to<Vector3<outFP_t>>;
  };

template<typename Field, typename inFP_t = double, typename outFP_t = inFP_t>
concept ScalarField = 
  std::floating_point<inFP_t> &&
  std::floating_point<outFP_t> &&
  std::is_default_constructible_v<Field> &&
  requires(const Field model, Vector3<inFP_t> pos) {
    { model.getField(pos) } -> std::convertible_to<outFP_t>;
  };

template<typename Field, typename inFP_t = double, typename outFP_t = inFP_t>
concept PotentialField = ScalarField<Field, inFP_t, outFP_t>;

/**
 *@brief the pos is in Re
*/
template<typename Field, typename fp_t = double>
concept MagneticFieldModel = VectorField<Field, Re<fp_t>, microTesla<fp_t>>;

template<typename Field, typename fp_t = double>
concept ElectricPotentialModel = PotentialField<Field, Re<fp_t>, kV<fp_t>>;

}
