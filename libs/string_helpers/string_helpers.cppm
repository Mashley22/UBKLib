module;

#include <concepts>
#include <string>
#include <string_view>

export module UBKLib.string_helpers;

import UBKLib.utils;

export namespace ubk {

template<std::floating_point T>
[[nodiscard]] std::string
vec3ToStr(const Vector3<T>& vector, std::string_view seperator = ",") {
  std::string str;
  str.append(std::to_string(vector.x));
  str.append(seperator);
  str.append(std::to_string(vector.y));
  str.append(seperator);
  str.append(std::to_string(vector.z));
  return str;
}

}
