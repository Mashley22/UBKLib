module;

#include <concepts>
#include <string>

export module UBKLib.string_helpers;

import UBKLib.utils;

export namespace ubk {

template<std::floating_point T, std::size_t separation = 4>
[[nodiscard]] std::string
vec3ToStr(const Vector3<T>& vector) {
  std::string str;
  str.append(std::to_string(vector.x));
  str.append(5, ' ');
  str.append(std::to_string(vector.y));
  str.append(5, ' ');
  str.append(std::to_string(vector.z));
  return str;
}

}
