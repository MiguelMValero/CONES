#ifndef __PDM_PRIV_H__
#define __PDM_PRIV_H__

#include <stdio.h>
#include <math.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Suppress warning
 *============================================================================*/

#if defined(__INTEL_COMPILER)
#define PDM_PRAGMA_TO_STR(x) _Pragma(#x)
#define PDM_INTEL_SUPPRESS_WARNING_PUSH _Pragma("warning(push)")
#define PDM_INTEL_SUPPRESS_WARNING(w) PDM_PRAGMA_TO_STR(warning(disable:w))
#define PDM_INTEL_SUPPRESS_WARNING_POP _Pragma("warning(pop)")
#define PDM_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    PDM_INTEL_SUPPRESS_WARNING_PUSH PDM_INTEL_SUPPRESS_WARNING(w)
#else // PDM_INTEL
#define PDM_INTEL_SUPPRESS_WARNING_PUSH
#define PDM_INTEL_SUPPRESS_WARNING(w)
#define PDM_INTEL_SUPPRESS_WARNING_POP
#define PDM_INTEL_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // PDM_INTEL

#if defined(__clang__)
#define PDM_PRAGMA_TO_STR(x) _Pragma(#x)
#define PDM_CLANG_SUPPRESS_WARNING_PUSH _Pragma("clang diagnostic push")
#define PDM_CLANG_SUPPRESS_WARNING(w) PDM_PRAGMA_TO_STR(clang diagnostic ignored w)
#define PDM_CLANG_SUPPRESS_WARNING_POP _Pragma("clang diagnostic pop")
#define PDM_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)                                                \
    PDM_CLANG_SUPPRESS_WARNING_PUSH PDM_CLANG_SUPPRESS_WARNING(w)
#else // PDM_CLANG
#define PDM_CLANG_SUPPRESS_WARNING_PUSH
#define PDM_CLANG_SUPPRESS_WARNING(w)
#define PDM_CLANG_SUPPRESS_WARNING_POP
#define PDM_CLANG_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // PDM_CLANG

#if defined(__GNUC__) && !defined(__INTEL_COMPILER)
#define PDM_PRAGMA_TO_STR(x) _Pragma(#x)
#define PDM_GCC_SUPPRESS_WARNING_PUSH _Pragma("GCC diagnostic push")
#define PDM_GCC_SUPPRESS_WARNING(w) PDM_PRAGMA_TO_STR(GCC diagnostic ignored w)
#define PDM_GCC_SUPPRESS_WARNING_POP _Pragma("GCC diagnostic pop")
#define PDM_GCC_SUPPRESS_WARNING_WITH_PUSH(w)                                                  \
    PDM_GCC_SUPPRESS_WARNING_PUSH PDM_GCC_SUPPRESS_WARNING(w)
#else // PDM_GCC
#define PDM_GCC_SUPPRESS_WARNING_PUSH
#define PDM_GCC_SUPPRESS_WARNING(w)
#define PDM_GCC_SUPPRESS_WARNING_POP
#define PDM_GCC_SUPPRESS_WARNING_WITH_PUSH(w)
#endif // DOCTEST_GCC

PDM_GCC_SUPPRESS_WARNING("-Wcast-qual")
PDM_GCC_SUPPRESS_WARNING("-Wunknown-pragmas")

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Macro definitions
 *============================================================================*/

enum {X, Y, Z};

/**
 * Absolute value
 */

#define PDM_ABS(a)     ((a) <  0  ? -(a) : (a))

/**
 * Minimum value
 */

#define PDM_MIN(a,b)   ((a) < (b) ?  (a) : (b))

/**
 * Maximum value
 */

#define PDM_MAX(a,b)   ((a) > (b) ?  (a) : (b))

/**
 * Sign of value
 */

#define PDM_SIGN(a)  (((a) > 0) ? 1 : (((a) < 0) ? -1 : 0))

/**
 * Dot product
 */

#define PDM_DOT_PRODUCT(vect1, vect2)                                   \
    ((vect1)[X] * (vect2)[X] + (vect1)[Y] * (vect2)[Y] + (vect1)[Z] * (vect2)[Z])

/**
 * Module
 */

#define PDM_MODULE(vect)                                                \
  sqrt((vect)[X] * (vect)[X] + (vect)[Y] * (vect)[Y] + (vect)[Z] * (vect)[Z])

#define PDM_CROSS_PRODUCT(prod_vect, vect1, vect2)  \
    ((prod_vect)[X] = (vect1)[Y] * (vect2)[Z] - (vect2)[Y] * (vect1)[Z], \
     (prod_vect)[Y] = (vect2)[X] * (vect1)[Z] - (vect1)[X] * (vect2)[Z], \
     (prod_vect)[Z] = (vect1)[X] * (vect2)[Y] - (vect2)[X] * (vect1)[Y])

#define PDM_DETERMINANT2X2(vect1, vect2) \
    ((vect1)[X] * (vect2)[Y] - (vect2)[X] * (vect1)[Y] )

#define PDM_DOT_PRODUCT_2D(vect1, vect2) \
    ((vect1)[X] * (vect2)[X] + (vect1)[Y] * (vect2)[Y])

#define PDM_PI 3.1415926535897931

#define PDM_SQRT2     1.415//1.4142135623730951
#define PDM_SQRT2_INV 0.707//0.7071067811865476

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PRIV_H__ */
