#ifndef __BFTC_MEM_H__
#define __BFTC_MEM_H__

/*============================================================================
 * Base memory allocation wrappers with optional tracing
 *============================================================================*/

/*
  This file is part of the "Base Functions and Types" library, intended to
  simplify and enhance portability, memory and I/O use for scientific codes.

  Copyright (C) 2004-2006  EDF

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

/*----------------------------------------------------------------------------*/

/*
 * Obtain definitions such as that of size_t through stddef.h (C99 standard)
 * if available (preferred method), or through stdlib.h (which defines
 * malloc() and family and so must define size_t some way) otherwise.
 */

#if defined(__STDC_VERSION__)
#  if (__STDC_VERSION__ == 199901L)
#    include <stddef.h>
#  else
#    include <stdlib.h>
#  endif
#else
#  include <stdlib.h>
#endif

/* BFT headers */

#include "bftc_config.h"
#include "bftc_error.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public macros
 *============================================================================*/

/*
 * Allocate memory for _ni items of type _type.
 *
 * This macro calls bftc_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  --> pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define BFTC_MALLOC(_ptr, _ni, _type) \
_ptr = (_type *) bftc_mem_malloc(_ni, sizeof(_type), \
                                #_ptr, __FILE__, __LINE__)

/*
 * Reallocate memory for _ni items of type _type.
 *
 * This macro calls bftc_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 *   _ni   <-- number of items.
 *   _type <-- element type.
 */

#define BFTC_REALLOC(_ptr, _ni, _type) \
_ptr = (_type *) bftc_mem_realloc(_ptr, _ni, sizeof(_type), \
                                 #_ptr, __FILE__, __LINE__)

/*
 * Free allocated memory.
 *
 * This macro calls bftc_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * parameters:
 *   _ptr  <->  pointer to allocated memory.
 */

#ifdef __cplusplus /* avoid casting from void for C++ */

#define BFTC_FREE(_ptr) \
bftc_mem_free(_ptr, #_ptr, __FILE__, __LINE__), _ptr = NULL

#else

#define BFTC_FREE(_ptr) \
_ptr = bftc_mem_free(_ptr, #_ptr, __FILE__, __LINE__)

#endif /* __cplusplus */

/*
 * Allocate aligned memory for _ni items of type _type.
 *
 * This macro calls bftc_mem_memalign(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * parameters:
 *   _ptr    --> pointer to allocated memory.
 *   _align <-- alignment.
 *   _ni    <-- number of items.
 *   _type  <-- element type.
 */

#define BFTC_MEMALIGN(_ptr, _align, _ni, _type) \
_ptr = (_type *) bftc_mem_memalign(_align, _ni, sizeof(_type), \
                                  #_ptr, __FILE__, __LINE__)

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*
 * Initialize memory handling.
 *
 * This function should be called before any other bftc_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * parameter:
 *   log_file_name <-- name of optional log_file (if NULL, no log).
 */

void
bftc_mem_init(const char  *log_file_name);

/*
 * End memory handling.
 *
 * This function should be called after all other bftc_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 */

void
bftc_mem_end(void);

/*
 * Initialize memory handling.
 *
 * This function should be called before any other bftc_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * parameter:
 *   log_file_name <-- name of optional log_file (if NULL, no log).
 */

/*
 * Indicates if bftc_mem_...() functions are initialized.
 *
 * returns:
 *   1 if bftc_mem_init has been called, 0 otherwise.
 */

int
bftc_mem_initialized(void);

/*
 * Allocate memory for ni items of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the bftc_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * parameters:
 *   ni        <-- number of items.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to allocated memory.
 */

void *
bftc_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num);

/*
 * Reallocate memory for ni items of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the bftc_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, bftc_alloc() called).
 *   ni        <-- number of items.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num   -> line number in calling source file
 *
 * returns:
 *   pointer to allocated memory.
 */

void *
bftc_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num);

/*
 * Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the bftc_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * parameters:
 *   ptr       <-> pointer to previous memory location
 *                 (if NULL, bftc_alloc() called).
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   NULL pointer.
 */

void *
bftc_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num);

/*
 * Allocate aligned memory for ni elements of size bytes.
 *
 * This function calls posix_memalign() if available, but adds tracing
 * capabilities, and automatically calls the bftc_error() errorhandler if
 * it fails to allocate the required memory.
 *
 * The associated function bftc_mem_have_memalign() indicates if this
 * type of allocation may be used on this system.
 *
 * parameters:
 *   alignment <-- alignent.
 *   ni        <-- number of items.
 *   size      <-- element size.
 *   var_name  <-- allocated variable name string.
 *   file_name <-- name of calling source file.
 *   line_num  <-- line number in calling source file.
 *
 * returns:
 *   pointer to allocated memory.
 */

void *
bftc_mem_memalign(size_t       alignment,
                 size_t       ni,
                 size_t       size,
                 const char  *var_name,
                 const char  *file_name,
                 int          line_num);

/*
 * Indicate if a memory aligned allocation variant is available.
 *
 * If no such function is available, bftc_mem_memalign() will always fail.
 *
 * returns:
 *   1 if memory aligned allocation is possible, 0 otherwise.
 */

int
bftc_mem_have_memalign(void);

/*
 * Return current theoretical dynamic memory allocated.
 *
 * returns:
 *   current memory handled through bftc_mem_...() (in kB).
 */

size_t
bftc_mem_size_current(void);

/*
 * Return maximum theoretical dynamic memory allocated.
 *
 * returns:
 *   maximum memory handled through bftc_mem_...() (in kB).
 */

size_t
bftc_mem_size_max(void);

/* Returns the error handler associated with the bftc_mem_...() functions.
 *
 * returns:
 *   pointer to the error handler function.
 */

bftc_error_handler_t *
bftc_mem_error_handler_get(void);

/*
 * Associates an error handler with the bftc_mem_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after bftc_print_flush() is called), and the general error handler used
 * by bftc_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameter:
 *   handler <-- pointer to the error handler function.
 */

void
bftc_mem_error_handler_set(bftc_error_handler_t *handler);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __BFTC_MEM_H__ */
