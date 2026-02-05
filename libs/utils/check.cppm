module;

#include<source_location>
#include <string>
#include <exception>

export module UBKLib.utils:check;

export namespace ubk {

struct CheckFail {
  std::string msg;
  std::source_location loc;
};

constexpr void
check(bool expr,
      std::string&& msg = {},
      std::source_location loc = std::source_location::current()) {
  if (!expr) {
    #ifndef UBKLIB_UNIT_TEST
    (void)msg;
    (void)loc;
    std::terminate();
    #else 
    throw CheckFail{.msg = msg, .loc = loc};
    #endif 
  }
}

}
