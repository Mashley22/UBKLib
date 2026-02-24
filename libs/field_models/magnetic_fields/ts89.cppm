module;

#include <array>
#include <concepts>
#include <cmath>
#include <span>

#include <cxform.h>

#include <UBK/macros.hpp>

export module UBKLib.field_models:magnetic_fields.ts89;

import :traits;
import UBKLib.utils;

template<std::floating_point T>
[[nodiscard]] static ubk::Vector3<ubk::Re<T>>
m_geoToGsm(ubk::Vector3<ubk::Re<T>> val, const double time) UBK_NOEXCEPT {
  std::array<T, 3> out;
  cxform2(GEO, GSM, time, val.toArr().data(), out.data());

  return ubk::Vector3<ubk::Re<T>>::template fromArr<T>(out);
}

template<std::floating_point T>
[[nodiscard]] static ubk::Vector3<ubk::Re<T>>
m_gsmToGeo(ubk::Vector3<ubk::Re<T>> val, const double time) UBK_NOEXCEPT {
  std::array<T, 3> out;
  cxform2(GSM, GEO, time, val.toArr().data(), out.data());

  return ubk::Vector3<ubk::Re<T>>::template fromArr<T>(out);
}

export namespace ubk {



// NOTE NEED TO BE CAREFUL OF FLOATING POINT MATH STUFF, SHOULD BE FINE BUT YOU NEVER KNOW
// For ease of comparing to original will provide both big fuck off param tables as well as
// proper named variables, following the geopack python implementation
// C ordering has been used, although its constexpr so really should not be affecting performance
// may be some llm hallucination here but oh well
template<std::floating_point T>
class Ts89 { 
public:

  static constexpr T params[7][30] = {
    {
      -116.53, -10719, 42.375, 59.753, -11363, 1.7844, 30.268,
      -0.035372, -0.066832, 0.016456, -1.3024, 0.0016529, 0.0020293, 20.289,
      -0.025203, 224.91, -9234.8, 22.788, 7.8813, 1.8362, -0.27228,
      8.8184, 2.8714, 14.468, 32.177, 0.01, 0, 7.0459,
      4, 20
    },
    {
      -55.553, -13198, 60.647, 61.072, -16064, 2.2534, 34.407,
      -0.038887, -0.094571, 0.027154, -1.3901, 0.001346, 0.0013238, 23.005,
      -0.030565, 55.047, -3875.7, 20.178, 7.9693, 1.4575, 0.89471,
      9.4039, 3.5215, 14.474, 36.555, 0.01, 0, 7.0787,
      4, 20
    },
    {
      -101.34, -13480, 111.35, 12.386, -24699, 2.6459, 38.948,
      -0.03408, -0.12404, 0.029702, -1.4052, 0.0012103, 0.0016381, 24.49,
      -0.037705, -298.32, 4400.9, 18.692, 7.9064, 1.3047, 2.4541,
      9.7012, 7.1624, 14.288, 33.822, 0.01, 0, 6.7442,
      4, 20
    },
    {
      -181.69, -12320, 173.79, -96.664, -39051, 3.2633, 44.968,
      -0.046377, -0.16686, 0.048298, -1.5473, 0.0010277, 0.0031632, 27.341,
      -0.050655, -514.1, 12482, 16.257, 8.5834, 1.0194, 3.6148,
      8.6042, 5.5057, 13.778, 32.373, 0.01, 0, 7.3195,
      4, 20
    },
    {
      -436.54, -9001, 323.66, -410.08, -50340, 3.9932, 58.524,
      -0.038519, -0.26822, 0.074528, -1.4268, -0.0010985, 0.0096613, 27.557,
      -0.056522, -867.03, 20652, 14.101, 8.3501, 0.72996, 3.8149,
      9.2908, 6.4674, 13.729, 28.353, 0.01, 0, 7.4237,
      4, 20
    },
    {
      -707.77, -4471.9, 432.81, -435.51, -60400, 4.6229, 68.178,
      -0.088245, -0.21002, 0.11846, -2.6711, 0.0022305, 0.01091, 27.547,
      -0.05408, -424.23, 1100.2, 13.954, 7.5337, 0.89714, 3.7813,
      8.2945, 5.174, 14.213, 25.237, 0.01, 0, 7.0037,
      4, 20
    },
    {
      -1190.4, 2749.9, 742.56, -1110.3, -77193, 7.6727, 102.05,
      -0.096015, -0.74507, 0.11214, -1.3614, 0.0015157, 0.022283, 23.164,
      -0.074146, -2219.1, 48253, 12.714, 7.6777, 0.57138, 2.9633,
      9.3909, 9.7263, 11.123, 21.558, 0.01, 0, 4.4518,
      4, 20
    }
  };

  static constexpr T A02 = 25.0;
  static constexpr T XLW2 = 170;
  static constexpr T YN = 30;
  static constexpr T RPI = 0.318309890;
  static constexpr T RT = 30;

  static constexpr T XD = 0;
  static constexpr T XLD2 = 40;

  static constexpr T SXC = 4;
  static constexpr T XLWC2 = 50;

  [[nodiscard]] constexpr T 
  DYC_closureCurrentCharacteristicScale(void) const UBK_NOEXCEPT {
    return params[m_iop][29]; // 20
  }

  [[nodiscard]] constexpr T 
  GAM_rateOfTailSheetThickening(void) const UBK_NOEXCEPT {
    return params[m_iop][28]; // 4
  }

  [[nodiscard]] constexpr T 
  SX_discVectorPotentialScaleXOffset(void) const UBK_NOEXCEPT {
    return params[m_iop][27];
  }

  [[nodiscard]] constexpr T
  Q_dyXDependence(void) const UBK_NOEXCEPT {
    return params[m_iop][26]; // 0
  }

  [[nodiscard]] constexpr T
  DEL_tailCurrentRateOfThickeningAlongY(void) const UBK_NOEXCEPT {
    return params[m_iop][25]; // 0.01
  }
  
  [[nodiscard]] constexpr T
  P_dyCharacteristicWScale(void) const UBK_NOEXCEPT {
    return params[m_iop][24];
  }

  [[nodiscard]] constexpr T
  AT_tailCurrentCharacteristicRadius(void) const UBK_NOEXCEPT {
    return params[m_iop][23];
  }


  [[nodiscard]] constexpr T
  G_degreeOfTransverseBending(void) const UBK_NOEXCEPT {
    return params[m_iop][22];
  }

  [[nodiscard]] constexpr T
  RC_hingingDistance(void) const UBK_NOEXCEPT {
    return params[m_iop][21];
  }

  [[nodiscard]] constexpr T
  DD_ringCurrentNightToDaySideThickeningRate(void) const UBK_NOEXCEPT {
    return params[m_iop][20];
  }

  [[nodiscard]] constexpr T
  D0_tailCurrentSheetHalfThickness(void) const UBK_NOEXCEPT {
    return params[m_iop][19];
  }

  [[nodiscard]] constexpr T
  ADR_ringCurrentCharacteristicRadius(void) const UBK_NOEXCEPT {
    return params[m_iop][18];
  }

  [[nodiscard]] constexpr T
  DX_chapmanFerraroCharacteristicScale(void) const UBK_NOEXCEPT {
    return params[m_iop][17];
  }

  struct ZsCurrentSheetShapeInfo {
    T shape;
    T x_derivative;
    T y_derivative;
  };

  [[nodiscard]] constexpr ZsCurrentSheetShapeInfo
  zsTailCurrentSheetShape(const Vector3<T> coords,
                          const Vector3<T> coords_sm) const {

    T x_rc = coords_sm.x + RC_hingingDistance();
    T x_rc_16 = pow(x_rc, 2) + 16;

    T sx_rc = std::sqrt(x_rc_16);
    T zs1 = 0.5 * tan(m_dipole_tilt) * (x_rc - sx_rc);

    T y410 = pow(coords.y, 4) + static_cast<T>(1e4);
    T sy4 = sin(m_dipole_tilt) / y410;

    T gsy4 = G_degreeOfTransverseBending() * sy4;
    T d2zsgy = -1 * sy4 / y410 * 4e4 * pow(coords.y, 3);
    
    return {
      .shape = -1 * zs1 / sx_rc,
      .x_derivative = zs1 - gsy4 * pow(coords.y, 4),
      .y_derivative = G_degreeOfTransverseBending() * d2zsgy};
  }

  struct RingCurrentInfo {
    T zr;
    Vector3<T> der4;
  };

  [[nodiscard]] RingCurrentInfo 
  ringCurrent(const Vector3<T> coords, 
              const Vector3<T> coords_sm,
              const ZsCurrentSheetShapeInfo currentSheetShape) const {

    T dsqt = std::sqrt(pow(coords_sm.x, 2) + A02);
    T fa0 = 0.5 * (1 + coords_sm.x / dsqt);
    T ddr = D0_tailCurrentSheetHalfThickness() + DD_ringCurrentNightToDaySideThickeningRate() * fa0;
    T dfa0 = (A02 / pow(dsqt, 3)) / 2;
    T zr = coords_sm.z - currentSheetShape.shape;
    T tr = std::sqrt(pow(zr, 2) + pow(ddr, 2));
    T rtr = 1 / tr;
    T ro2 = pow(coords_sm.x, 2) + pow(coords.y, 2);
    T adrt = ADR_ringCurrentCharacteristicRadius() + tr;
    T fk = 1 / (pow(adrt, 2) + ro2);
    T fc = pow(fk, 2) * std::sqrt(fk);
    T facxy = 3.0 * adrt * fc * rtr;
    T xzr = coords_sm.x * zr;
    T yzr = coords.y * zr;
    T dbxdp = facxy * xzr;
    T xzyz = coords_sm.x * currentSheetShape.x_derivative + coords.y * currentSheetShape.y_derivative;
    T faq = zr * xzyz - ddr * DD_ringCurrentNightToDaySideThickeningRate() * dfa0 * coords_sm.x;
    T dbzdp = fc * (2 * pow(adrt, 2) - ro2) + facxy * faq;

    
    return {
      .zr = zr,
      .der4 = {
        .x = dbxdp * cos(m_dipole_tilt) + dbzdp * sin(m_dipole_tilt),
        .y = facxy * yzr,
        .z = dbzdp * cos(m_dipole_tilt) - dbxdp * sin(m_dipole_tilt)
      }
    }; 
  }
  
  struct TailCurrentSheetInfo {
    Vector3<T> der0;
    Vector3<T> der1;
    Vector3<T> der15;
    Vector3<T> der16;

    TailCurrentSheetInfo(Vector3<T> der0_val, Vector3<T> der1_val, T m_dipole_tilt) 
      : der0(der0_val), der1(der1_val), der15(der0_val * pow(m_dipole_tilt, 2)), der16(der1_val * pow(m_dipole_tilt, 2)) {}
  };

  [[nodiscard]] TailCurrentSheetInfo
  tailCurrentSheet(const Vector3<T> coords,
                   const Vector3<T> coords_sm,
                   const ZsCurrentSheetShapeInfo currentSheetShape,
                   const RingCurrentInfo ringCurrentInfo) const {

    T dely2 = DEL_tailCurrentRateOfThickeningAlongY() * pow(coords.y, 2);
    T d = D0_tailCurrentSheetHalfThickness() + dely2;
    
    T adsl = 0;

    if (std::abs(GAM_rateOfTailSheetThickening()) >= 1e-6) {
      T xxd = coords_sm.x - XD;
      T rqd = 1 / (pow(xxd, 2) + XLD2);
      T rqds = std::sqrt(rqd);
      T h = 0.5 * (1 + xxd * rqds);
      T hs = -0.5 * XLD2 * rqd * rqds;
      T gamh = GAM_rateOfTailSheetThickening() * h;
      d += gamh;
      T xghs = coords_sm.x * GAM_rateOfTailSheetThickening() * hs;
      adsl = -d * xghs;
    }

    T t = std::sqrt(pow(ringCurrentInfo.zr, 2) + pow(d, 2));
    T xsmx = coords_sm.x - SX_discVectorPotentialScaleXOffset();
    T rdsq = std::sqrt(1 / (pow(xsmx, 2) + XLW2));
    T v = (1 - xsmx * rdsq) / 2;
    T dvx = (XLW2 / 2) * pow(rdsq, 3);
    T om = std::sqrt(std::sqrt(pow(coords_sm.x, 2) + 16) - coords_sm.x);
    T oms =- om / (om * om + coords_sm.x) / 2;
    T rdy = 1 / (P_dyCharacteristicWScale() + Q_dyXDependence() * om);
    T omsv = oms * v;
    T fy = 1 / (1 + pow(coords.y * rdy, 2));
    T w = v * fy;
    T yfy1 = 2 * fy * pow(coords.y * rdy, 2);
    T fypr = yfy1 * rdy;
    T fydy = fypr * fy;
    T dwx = dvx * fy + fydy * Q_dyXDependence() * omsv;
    T ydwy =- v * yfy1 * fy;
    T ddy = 2 * DEL_tailCurrentRateOfThickeningAlongY() * coords.y;
    T att = AT_tailCurrentCharacteristicRadius() + t;
    T s1 = std::sqrt(pow(att, 2) + pow(coords_sm.x, 2) + pow(coords.y, 2));
    T f5 = 1 / s1;
    T f7 = 1 / (s1 + att);
    T f1 = f5 * f7;
    T f9 = att * pow(f5 , 3);
    T fs = ringCurrentInfo.zr *
           coords_sm.x * currentSheetShape.x_derivative + coords.y * currentSheetShape.y_derivative -
           d * coords.y * ddy + adsl;
    T xdwx = coords_sm.x * dwx + ydwy;
    T rtt = 1 / t;
    T wt = w * rtt;
    T brrz1 = wt * f1;
    T brrz2 = wt * pow(f5, 3);
    T dbxc1 = brrz1 * coords_sm.x * ringCurrentInfo.zr;
    T dbxc2 = brrz2 * coords_sm.x * ringCurrentInfo.zr;

    T wtfs = wt * fs;
    T dbzc1 = w * f5 + xdwx * f7 + wtfs * f1;
    T dbzc2 = w * f9 + xdwx * f1 + wtfs * pow(f5, 3);
    
    return {
      { .x = dbxc1 * cos(m_dipole_tilt) + dbzc1 * sin(m_dipole_tilt),
        .y = brrz1 * coords.y * ringCurrentInfo.zr,
        .z = dbzc1 * cos(m_dipole_tilt) - dbxc1 * sin(m_dipole_tilt)}, 

      { .x = dbxc2 * cos(m_dipole_tilt) + dbzc2 * sin(m_dipole_tilt),
        .y = brrz2 * ringCurrentInfo.zr,
        .z = dbzc2 * cos(m_dipole_tilt) - dbxc2 * sin(m_dipole_tilt)},
      
      m_dipole_tilt
    };
  }

  [[nodiscard]] Vector3<nanoTesla<T>>
  tailCurrentSheetField(const Vector3<T> coords,
                        const Vector3<T> coords_sm,
                        const ZsCurrentSheetShapeInfo currentSheetShape,
                        const RingCurrentInfo ringCurrentInfo) const {

    auto tailCurrentSheetInfo = tailCurrentSheet(coords, coords_sm, currentSheetShape, ringCurrentInfo);

    return tailCurrentSheetInfo.der0 * m_params()[0] +
           tailCurrentSheetInfo.der1 * m_params()[1] +
           tailCurrentSheetInfo.der15 * m_params()[15] +
           tailCurrentSheetInfo.der16 * m_params()[16];
  }

  [[nodiscard]] constexpr Vector3<nanoTesla<T>>
  closureCurrentField(const Vector3<T> coords) const UBK_NOEXCEPT {

    T zpl = coords.z + RT;
    T zmn = coords.z - RT;
    T rho_2 = pow(coords.x, 2) + pow(coords.y, 2);
    T spl = std::sqrt(pow(zpl, 2));
    T smn = std::sqrt(pow(zmn, 2) + rho_2);
    T xsxc = coords.x - SXC;
    T rqc2 = 1 / (pow(xsxc, 2) + XLWC2);
    T rqc = std::sqrt(rqc2);
    T fyc = 1 / (1 + pow(coords.y / DYC_closureCurrentCharacteristicScale(), 2));
    T wc =  (1 - xsxc * rqc) * fyc / 2;
    T dwcx = XLWC2 * rqc2 * rqc * fyc / 2;
    T dwcy = -2 * wc * fyc * coords.y / pow(DYC_closureCurrentCharacteristicScale(), 2);
    T szrp = 1 / (spl + zpl);
    T szrm = 1 / (smn - zmn);
    T xywc = coords.x * dwcx + coords.y * dwcy;
    T wcsp = wc / spl;
    T wcsm = wc / smn;
    T fxyp = wcsp * szrp;
    T fxym = wcsm * szrm;
    T fxpl = coords.x * fxyp;
    T fxmn = -1 * coords.x * fxym;
    T fypl = coords.y * fxyp;
    T fymn = -1 * coords.y * fxym;
    T fzpl = wcsp + xywc * szrp;
    T fzmn = wcsm + xywc * szrm;

    Vector3<T> der2 = {
      .x = fxpl + fxmn,
      .y = fypl + fymn,
      .z = fzpl + fzmn
    };
    Vector3<T> der3 = {
      .x = (fxpl - fxmn) * sin(m_dipole_tilt),
      .y = (fypl - fymn) * sin(m_dipole_tilt),
      .z = (fzpl - fzmn) * sin(m_dipole_tilt)
    };

    return m_params()[2] * der2 + m_params()[3] * der3;
  }

  [[nodiscard]] Vector3<nanoTesla<T>>
  chapmanFerraroField(const Vector3<T> coords) const UBK_NOEXCEPT {

    T ex = std::exp(coords.x / DX_chapmanFerraroCharacteristicScale());
    T ec = ex * cos(m_dipole_tilt);
    T es = ex * sin(m_dipole_tilt);
    T ecz = ec * coords.z;
    T esz = es * coords.z;
    T eszy2 = esz * pow(coords.y, 2);
    T eszz2 = esz * pow(coords.z, 2);
    T ecz2 = ecz * coords.z;
    T esy = es * coords.y;

    return { 
      .x = m_params()[5] * ecz +
                  m_params()[6] * es +
                  m_params()[7] * esy * coords.y +
                  m_params()[8] * esz * coords.z,

      .y = m_params()[9] * ecz * coords.y +
           m_params()[10] * esy +
           m_params()[11] * esy * pow(coords.y, 2) + 
           m_params()[12] * esy * pow(coords.z, 2), 

      .z = m_params()[13] * ec +
           m_params()[14] * ec * pow(coords.y, 2) +
           (m_params()[5] / DX_chapmanFerraroCharacteristicScale() + m_params()[9]) * ecz2 / (-2) +
           -1 * (m_params()[6] / DX_chapmanFerraroCharacteristicScale() + m_params()[10]) * esz + 
           -1 * (m_params()[7] / DX_chapmanFerraroCharacteristicScale() + 3 * m_params()[11]) * eszy2 + 
           (-1.0 / 3.0) * (m_params()[8] / DX_chapmanFerraroCharacteristicScale() + m_params()[12]) * eszz2
    };
  }

  [[nodiscard]] Vector3<nanoTesla<T>>
  getField(Vector3<T> coords) const UBK_NOEXCEPT {
    coords = geoToGsm(coords);

    Vector3<T> coords_sm = {
      .x = coords.x * cos(m_dipole_tilt) - coords.z * sin(m_dipole_tilt),
      .z = coords.x * sin(m_dipole_tilt) + coords.z * cos(m_dipole_tilt)
    };
    
    ZsCurrentSheetShapeInfo zs = zsTailCurrentSheetShape(coords, coords_sm);
    RingCurrentInfo ringCurrentInfo = ringCurrent(coords, coords_sm, zs);

    Vector3<T> ringCurrentField = ringCurrentInfo.der4 * m_params()[4];

    return gsmToGeo(ringCurrentField + 
           tailCurrentSheetField(coords, coords_sm, zs, ringCurrentInfo) +
           chapmanFerraroField(coords) +
           closureCurrentField(coords));
  }

private:
  int m_iop{0};
  double m_time{0};
  T m_dipole_tilt{0};

  [[nodiscard]] constexpr std::span<const T>
  m_params(void) const UBK_NOEXCEPT{
    return {params[m_iop], 30};
  }

  [[nodiscard]] Vector3<T>
  geoToGsm(Vector3<T> val) const UBK_NOEXCEPT {
    return m_geoToGsm(val, m_time);
  }

  [[nodiscard]] Vector3<T>
  gsmToGeo(Vector3<T> val) const UBK_NOEXCEPT {
    return m_gsmToGeo(val, m_time);
  }

};

static_assert(MagneticFieldModel<Ts89<double>, double>);
static_assert(MagneticFieldModel<Ts89<float>, float>);

static_assert(!MagneticFieldModel<Ts89<float>, double>);
static_assert(!MagneticFieldModel<Ts89<double>, float>);

}
