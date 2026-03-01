module;

#include <concepts>
#include <cstdint>

#include <UBK/macros.hpp>

#include <cxform.h>

export module UBKLib.utils:units;

export namespace ubk {

template<std::floating_point T>
using Re = T;

template<std::floating_point T>
using kV = T;

template<std::floating_point T>
using nanoTesla = T;

/**
 *@brief a gregorian calender date and UT time
 */
struct Time {
  std::size_t year{1950};
  std::size_t month{1};
  std::size_t day{1};
  std::size_t hours{0};
  std::size_t minutes{0};
  std::size_t seconds{0};

};

[[nodiscard]] std::int64_t
timeToEs(const Time& time) UBK_NOEXCEPT {
  return date2es(
    static_cast<int>(time.year),
    static_cast<int>(time.month),
    static_cast<int>(time.day),
    static_cast<int>(time.hours),
    static_cast<int>(time.minutes),
    static_cast<int>(time.seconds)
  );
}

}
