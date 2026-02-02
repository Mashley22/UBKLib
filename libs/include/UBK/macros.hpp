#ifndef UBK_MACROS_HPP
#define UBK_MACROS_HPP

#define UBK_CHECK(...) \
  UBK_CHECK_IMPL(__VA_ARGS__, \
    UBK_CHECK_EXPR_MSG, \
    UBK_CHECK_EXPR_ONLY)(__VA_ARGS__)

#define UBK_CHECK_ASSUME(...) \
  UBK_CHECK_ASSUME_IMPL(__VA_ARGS__, \
    UBK_CHECK_ASSUME_EXPR_MSG, \
    UBK_CHECK_ASSUME_EXPR_ONLY)(__VA_ARGS__)

#define UBK_CHECK_IMPL(_1, _2, NAME, ...) NAME

#define UBK_CHECK_EXPR_ONLY(expr) ubk::check(expr)
#define UBK_CHECK_EXPR_MSG(expr, msg) ubk::check(expr, msg)

#define UBK_CHECK_ASSUME_IMPL(_1, _2, NAME, ...) NAME

#define UBK_CHECK_ASSUME_EXPR_ONLY(expr) \
do { \
  UBK_CHECK_EXPR_ONLY(expr) \
  [[asume(expr)]]; \
} while(0)

#define UBK_CHECK_ASSUME_EXPR_MSG(expr, msg) \
do { \
  UBK_CHECK_EXPR_MSG(expr, msg) \
  [[asume(expr)]]; \
} while(0)

#ifndef UBKLIB_CHECK_FAIL_THROWS
#define UBK_NOEXCEPT noexcept
#else
#define UBK_NOEXCEPT
#endif

#endif /* UBK_MACROS_HPP */
