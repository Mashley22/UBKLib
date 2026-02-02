module;

#include <concepts>

export module ubk.field_models:traits;

import ubk.utils;

export namespace ubk {

template<typename Field, typename fp_t = double>
concept VectorField =
  std::floating_point<fp_t> &&
  std::is_default_constructible_v<Field> &&
  requires(const Field model, Vector3<fp_t> pos) {
    { model.getField(pos) } -> std::convertible_to<Vector3<fp_t>>;
  };

template<typename Field, typename fp_t = double>
concept ScalarField = 
  std::floating_point<fp_t> &&
  std::is_default_constructible_v<Field> &&
  requires(const Field model, Vector3<fp_t> pos) {
    { model.getField(pos) } -> std::convertible_to<fp_t>;
  };

template<typename Field, typename fp_t = double>
concept PotentialField = ScalarField<Field, fp_t>;

/**
 *@brief the pos is in Re
*/
template<typename Field, typename fp_t = double>
concept MagneticFieldModel = VectorField<Field, fp_t>;

template<typename Field, typename fp_t = double>
concept ElectricPotentialModel = PotentialField<Field, fp_t>;

}
