module;

#include <concepts>

#include <UBK/macros.hpp>

export module UBKLib.field_models:traits;

import UBKLib.utils;

export namespace ubk {

template<std::floating_point T>
struct TimeDependentField {
public:
  TimeDependentField(void) UBK_NOEXCEPT = default;
  TimeDependentField(const Time& time) 
    : m_time(time), m_time_es(static_cast<T>(timeToEs(time))) {}

  [[nodiscard]] T
  time_es(void) const UBK_NOEXCEPT {
    return m_time_es;
  }

  [[nodiscard]] const Time&
  time(void) const UBK_NOEXCEPT {
    return m_time;
  }

  void
  setTime(const Time& newTime) UBK_NOEXCEPT {
    m_time = newTime;
    m_time_es = static_cast<T>(timeToEs(m_time));
  }
  
private:
  Time m_time;
  T m_time_es;
};

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
concept MagneticFieldModel = VectorField<Field, Re<fp_t>, nanoTesla<fp_t>>;

template<typename Field, typename fp_t = double>
concept ElectricPotentialModel = 
  ScalarField<Field, Re<fp_t>, kV<fp_t>> &&
  requires(const Field model, fp_t mlt, fp_t colatitude) {
    { model.getField(mlt, colatitude) } -> std::convertible_to<kV<fp_t>>;
  };

}
