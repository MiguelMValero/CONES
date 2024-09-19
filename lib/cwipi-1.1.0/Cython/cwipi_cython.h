#ifndef __CWIPI_CYTHON_H__
#define __CWIPI_CYTHON_H__

#include "cwipi.h"
#include "cwipi_priv.h"
#include "cwp.h"
#include "cwp_priv.h"

// CWIPI
CWP_CLANG_SUPPRESS_WARNING("-Wshadow")
CWP_CLANG_SUPPRESS_WARNING("-Wdeprecated")
CWP_CLANG_SUPPRESS_WARNING("-Wpedantic")
CWP_CLANG_SUPPRESS_WARNING("-W#warnings")
CWP_CLANG_SUPPRESS_WARNING("-Wincompatible-function-pointer-types")
CWP_CLANG_SUPPRESS_WARNING("-Wdeprecated-declarations")

typedef enum  {
  CWIPI_FAKE_ENUM
} cwipi_fake_enum_t;


#endif 
