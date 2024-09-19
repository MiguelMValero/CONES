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

/*-----------------------------------------------------------------------------*/

#include "bftc_config_defs.h"

/*
 * Standard C library headers
 */

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * Optional library and BFT headers
 */

#include "bftc_error.h"
#include "bftc_mem_usage.h"
#include "bftc_mem.h"
#include "bftc_printf.h"
#include "cwp_priv.h"

/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*-------------------------------------------------------------------------------
 * Local type definitions
 *-----------------------------------------------------------------------------*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS

/*
 * Structure defining an allocated memory block (for memory tracing)
 */

struct _bftc_mem_block_t {

  void    *p_bloc;  /* Allocated memory block start adress */
  size_t   size;    /* Allocated memory block length */

};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

/*-------------------------------------------------------------------------------
 * Local macro documentation
 *-----------------------------------------------------------------------------*/

/*! \fn BFTC_MALLOC(_ptr, _ni, _type)
 * \brief Allocate memory for _ni elements of type _type.
 *
 * This macro calls bftc_mem_malloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr  pointer to allocated memory.
 * \param [in]  _ni   number of elements.
 * \param [in]  _type element type.
 */

/*! \fn BFTC_REALLOC(_ptr, _ni, _type)
 * \brief Reallocate memory for _ni elements of type _type.
 *
 * This macro calls bftc_mem_realloc(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 * \param [in]      _ni   number of elements.
 * \param [in]      _type element type.
 */

/*! \fn BFTC_FREE(_ptr)
 * \brief Free allocated memory.
 *
 * This macro calls bftc_mem_free(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * The freed pointer is set to NULL to avoid accidental reuse.
 *
 * \param [in, out] _ptr  pointer to allocated memory.
 */

/*! \fn BFTC_MEMALIGN(_ptr, _align, _ni, _type)
 * \brief Allocate aligned memory for _ni elements of type _type.
 *
 * This macro calls bftc_mem_memalign(), automatically setting the
 * allocated variable name and source file name and line arguments.
 *
 * \param [out] _ptr   pointer to allocated memory.
 * \param [in]  _align alignment.
 * \param [in]  _ni    number of elements.
 * \param [in]  _type  element type.
 */

/*-------------------------------------------------------------------------------
 * Local macro definitions
 *-----------------------------------------------------------------------------*/

/*-----------------------------------------------------------------------------
 * Local function prototypes
 *-----------------------------------------------------------------------------*/

/*
 * Default memory error handler.
 *
 * Memory status is written to stderr, and the general error handler
 * used by bftc_error() is then called (which results in the
 * termination of the current process).
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_bftc_mem_error_handler_default(const char  *file_name,
                               int          line_num,
                               int          sys_error_code,
                               const char  *format,
                               va_list      arg_ptr);

/*-----------------------------------------------------------------------------
 * Local static variable definitions
 *-----------------------------------------------------------------------------*/

static int  _bftc_mem_global_initialized = 0;

static FILE *_bftc_mem_global_file = NULL;

static struct _bftc_mem_block_t  *_bftc_mem_global_block_array = NULL;

static unsigned long  _bftc_mem_global_block_nbr = 0 ;
static unsigned long  _bftc_mem_global_block_max = 512 ;

static size_t  _bftc_mem_global_alloc_cur = 0;
static size_t  _bftc_mem_global_alloc_max = 0;

static size_t  _bftc_mem_global_n_allocs = 0;
static size_t  _bftc_mem_global_n_reallocs = 0;
static size_t  _bftc_mem_global_n_frees = 0;

static bftc_error_handler_t  *_bftc_mem_error_handler
                              = (_bftc_mem_error_handler_default);

/*-----------------------------------------------------------------------------
 * Local function definitions
 *-----------------------------------------------------------------------------*/

/*
 * Given a character string representing a file name, returns
 * pointer to that part of the string corresponding to the base name.
 *
 * parameters:
 *   file_name:      <-- full name of source file.
 *
 * return:
 *   pointer to part of file name corresponding to base name.
 */

static const char *
_bftc_mem_basename(const char  *file_name)
{
  int i;

  if (file_name == NULL)
    return NULL;

  for (i = strlen(file_name) - 1;
       i > 0 && file_name[i] != BFTC_DIR_SEPARATOR;
       i--);

  if (file_name[i] == BFTC_DIR_SEPARATOR)
    i++;

  return (file_name + i);
}

/*
 * Determines values associated with an array representing a
 * long integer
 *
 * parameters:
 *   counter: <-- counter to update.
 *   value:   --> counter values in output unit [out].
 *   unit:    --> counter value unit : ' ', 'k', 'm', 'g', 't', or 'p'
 *                for bytes, Kilobytes, Megabytes, Gigabytes, Terabytes,
 *                or Petabytes [out].
 */

static void
_bftc_mem_size_val(const size_t    counter,
                  unsigned long   value[2],
                  char           *unit)
{
  int i;
  size_t _counter[2];
  const char units[] = {' ', 'k', 'm', 'g', 't', 'p', 'e'};

  for (i = 0, _counter[0] = counter, _counter[1] = 0;
       _counter[0] >= 1024 && i < 6;
       i++) {
    _counter[1] = _counter[0] % 1024;
    _counter[0] /= 1024;
  }

  value[0] = _counter[0];
  value[1] = _counter[1];
  *unit = units[i];
}

/*
 * Memory usage summary.
 */

static void _bftc_mem_summary(FILE  *f)
{
  char unit;
  unsigned long value[2];
  size_t mem_usage;

  if (f == NULL)
    return;

  fprintf(f, "\n\n");
  fprintf(f, _("Memory allocation summary\n"
               "-------------------------\n\n"));

  /* Available memory usage information */

  _bftc_mem_size_val(_bftc_mem_global_alloc_cur, value, &unit);
  fprintf(f,
          _("Theoretical current allocated memory:   %8lu.%lu %cB\n"),
          value[0], value[1], unit);

  _bftc_mem_size_val(_bftc_mem_global_alloc_max, value, &unit);
  fprintf(f,
          _("Theoretical maximum allocated memory:   %8lu.%lu %cB\n"),
          value[0], value[1], unit);

  fprintf(f,
          _("\n"
            "Number of allocations:   %lu\n"
            "          reallocations: %lu\n"
            "          frees:         %lu\n\n"),
          (unsigned long)_bftc_mem_global_n_allocs,
          (unsigned long)_bftc_mem_global_n_reallocs,
          (unsigned long)_bftc_mem_global_n_frees);

  if (bftc_mem_usage_initialized() == 1) {

    /* Maximum measured memory */

    mem_usage = bftc_mem_usage_max_pr_size();
    if (mem_usage > 0) {
      fprintf(f,
              _("Maximum program memory measure:  %8lu kB\n"),
              (unsigned long)mem_usage);
    }

    /* Current measured memory */

    mem_usage = bftc_mem_usage_pr_size();
    if (mem_usage > 0)
      fprintf(f,
              _("Current program memory measure:   %8lu kB\n"),
              (unsigned long)mem_usage);
  }

}

/*
 * Default memory error handler.
 *
 * Memory status is written to stderr (after bftc_print_flush() is called),
 * and the general error handler used by bftc_error() is then called
 * (which results in the termination of the current process).
 *
 * parameters:
 *   file_name:      <-- name of source file from which error handler called.
 *   line_num:       <-- line of source file from which error handler called.
 *   sys_error_code: <-- error code if error in system or libc call, 0 otherwise.
 *   format:         <-- format string, as printf() and family.
 *   arg_ptr:        <-> variable argument list based on format string.
 */

static void
_bftc_mem_error_handler_default(const char  *file_name,
                               int          line_num,
                               int          sys_error_code,
                               const char  *format,
                               va_list      arg_ptr)
{
  bftc_error_handler_t * general_err_handler;

  bftc_printf_flush();

  _bftc_mem_summary(stderr);

  general_err_handler = bftc_error_handler_get();
  general_err_handler(file_name, line_num, sys_error_code, format, arg_ptr);
}

/*
 * Calls the error handler (set by bftc_mem_error_handler_set() or default).
 *
 * With the default error handler, an error message is output to stderr,
 * (after bftc_print_flush() is called), and the general error handler used
 * by bftc_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * parameters:
 *   file_name      <-- name of source file from which failed bftc_mem_...()
 *                      function was called.
 *   line_num       <-- line of source file from which failed bftc_mem_...()
 *                      function was called.
 *   sys_error_code <-- error code if error in system or libc call,
 *                      0 otherwise.
 *   format         <-- format string, as printf() and family.
 *   ...            <-- variable arguments based on format string.
 */

static void
_bftc_mem_error(const char  *file_name,
               int          line_num,
               int          sys_error_code,
               const char  *format,
               ...)
{
  va_list  arg_ptr;

  if (_bftc_mem_global_file != NULL) {
    _bftc_mem_summary(_bftc_mem_global_file);
    fflush(_bftc_mem_global_file);
  }

  va_start(arg_ptr, format);

  _bftc_mem_error_handler(file_name, line_num, sys_error_code, format, arg_ptr);

  va_end(arg_ptr);
}

/*
 * Return the _bftc_mem_block structure corresponding to a given
 * allocated block.
 *
 * parameters:
 *   p_in: <-- allocated block's start adress.
 *
 * returns:
 *   corresponding _bftc_mem_block structure.
 */

static struct _bftc_mem_block_t *
_bftc_mem_block_info(const void *p_get)
{
  struct _bftc_mem_block_t  *pinfo = NULL;
  unsigned long idx;

  if (_bftc_mem_global_block_array != NULL) {

    for (idx = _bftc_mem_global_block_nbr - 1;
         idx > 0 && (_bftc_mem_global_block_array + idx)->p_bloc != p_get;
         idx--);

    if ((_bftc_mem_global_block_array + idx)->p_bloc != p_get)
      _bftc_mem_error(__FILE__, __LINE__, 0,
                     _("Adress [%10p] does not correspond to "
                       "the beginning of an allocated block."),
                     p_get);
    else {
      pinfo = _bftc_mem_global_block_array + idx;
      assert(p_get == pinfo->p_bloc);
    }

  }

  return pinfo;
}

/*
 * Return the size of a given allocated block.
 *
 * parameters:
 *   p_in: <-- allocated block's start adress.
 *
 * returns:
 *   block size.
 */

static size_t
_bftc_mem_block_size(const void  *p_in)
{
  struct _bftc_mem_block_t *pinfo = _bftc_mem_block_info(p_in);

  if (pinfo != NULL)
    return pinfo->size;
  else
    return 0;
}

/*
 * Fill a _bftc_mem_block_t structure for an allocated pointer.
 */

static void
_bftc_mem_block_malloc(void          *p_new,
                      const size_t   size_new)
{
  struct _bftc_mem_block_t *pinfo;

  assert(size_new != 0);

  if (_bftc_mem_global_block_array == NULL)
    return;

  if (_bftc_mem_global_block_nbr >= _bftc_mem_global_block_max) {

    _bftc_mem_global_block_max *= 2;
    _bftc_mem_global_block_array
      = (struct _bftc_mem_block_t *) realloc(_bftc_mem_global_block_array,
                                            sizeof(struct _bftc_mem_block_t)
                                            * _bftc_mem_global_block_max);

    if (_bftc_mem_global_block_array == NULL) {
      _bftc_mem_error(__FILE__, __LINE__, errno,
                     _("Memory allocation failure"));
      return;
    }

  }

  _bftc_mem_global_block_nbr += 1;

  pinfo = _bftc_mem_global_block_array + _bftc_mem_global_block_nbr - 1;

  /* Start adress and size of allocated block */

  pinfo->p_bloc = p_new;
  pinfo->size   = size_new;
}

/*
 * Update a _bftc_mem_block_t structure for an reallocated pointer.
 */

static void
_bftc_mem_block_realloc(const void    *p_old,
                       void          *p_new,
                       size_t         size_new)
{
  struct _bftc_mem_block_t *pinfo;

  assert(size_new != 0);

  pinfo = _bftc_mem_block_info(p_old);

  if (pinfo != NULL) {
    pinfo->p_bloc = p_new;
    pinfo->size   = size_new;
  }
}

/*
 * Free a _bftc_mem_block_t structure for a freed pointer.
 */

static void
_bftc_mem_block_free(const void *p_free)
{
  struct _bftc_mem_block_t *pinfo, *pmove;
  unsigned long idx;

  if (_bftc_mem_global_block_array == NULL)
    return;

  for (idx = _bftc_mem_global_block_nbr - 1;
       idx > 0 && (_bftc_mem_global_block_array + idx)->p_bloc != p_free;
       idx--);

  if ((_bftc_mem_global_block_array + idx)->p_bloc != p_free)
    _bftc_mem_error(__FILE__, __LINE__, 0,
                   _("Adress [%10p] does not correspond to "
                     "the beginning of an allocated block."),
                   p_free);

  else {

    /* We move the contents of the array's final block to the position
       of the freed block, and shorten the array's useful part by one. */

    pinfo = _bftc_mem_global_block_array + idx;
    pmove = _bftc_mem_global_block_array + _bftc_mem_global_block_nbr - 1;
    pinfo->p_bloc = pmove->p_bloc;
    pinfo->size   = pmove->size;

    _bftc_mem_global_block_nbr -= 1;

  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*!
 * \brief Initialize memory handling.
 *
 * This function should be called before any other bftc_mem_...()
 * function. To activate memory allocation logging, a logfile
 * name should be given as an argument. The resulting file will
 * be a regular, local file. If this file cannot be opened for
 * some reason, logging is silently de-activated.
 *
 * \param log_file_name name of optional log_file (if NULL, no log).
 */

void
bftc_mem_init(const char *log_file_name)
{
  size_t alloc_size;

  if (_bftc_mem_global_initialized == 1) {
    _bftc_mem_error(__FILE__, __LINE__, 0,
                   _("bftc_mem_init() has already been called"));
  }
  _bftc_mem_global_initialized = 1;

  alloc_size = sizeof(struct _bftc_mem_block_t) * _bftc_mem_global_block_max;

  _bftc_mem_global_block_array
    =  (struct _bftc_mem_block_t  *) malloc(sizeof(struct _bftc_mem_block_t) * _bftc_mem_global_block_max);

  if (_bftc_mem_global_block_array == NULL) {
    _bftc_mem_error(__FILE__, __LINE__, errno,
                   _("Failure to allocate \"%s\" (%lu bytes)"),
                   "_bftc_mem_global_block_array", (unsigned long)alloc_size);
    return;
  }

  if (log_file_name != NULL) {

    _bftc_mem_global_file = fopen(log_file_name, "w");

    /*
      If the file could not be opened, we do not abort, as it is not
      absolutely necessary. We silently continue.
      (We could warn the user, but this would require either using
      bftc_printf(), which we prefer to keep independent of the bftc_mem_...()
      functions to avoid evental crossed definitions when user-defined, or
      "warning handling" similar to error handling, with a possibility
      of user-defined warning handlers, as we are not sure if the calling
      code uses stderr (especially in a distributed environment). This
      is probably not worth the bother.
    */

    if (_bftc_mem_global_file == NULL)
      fprintf(stderr,
              _("Failure to open memory log file \"%s\"\n"),
              log_file_name);

  }

  /* Log file header */

  if (_bftc_mem_global_file != NULL) {

    fprintf(_bftc_mem_global_file,
            _("       :     FILE NAME              : LINE  :"
              "  POINTER NAME                          : N BYTES   :"
              " (+- N BYTES) : TOTAL BYTES  : [    ADRESS]\n"
              "-------:----------------------------:-------:"
              "----------------------------------------:-----------:"
              "-----------------------------:--------------"));

  }

}

/*!
 * \brief End memory handling.
 *
 * This function should be called after all other bftc_mem_...()
 * functions. In case of memory allocation logging, it
 * writes final information to the log file and closes is.
 */

void bftc_mem_end(void)
{
  if (_bftc_mem_global_initialized == 0) {
    _bftc_mem_error(__FILE__, __LINE__, 0,
                   _("bftc_mem_end() called before bftc_mem_init()"));
  }
  _bftc_mem_global_initialized = 0;

  if (_bftc_mem_global_file != NULL) {

    unsigned long  non_free = 0;
    struct _bftc_mem_block_t  *pinfo;

    /* Memory usage summary */

    _bftc_mem_summary(_bftc_mem_global_file);

    /* List of non-freed pointers */

    if (_bftc_mem_global_block_array != NULL) {

      fprintf(_bftc_mem_global_file, "List of non freed pointers:\n");

      for (pinfo = _bftc_mem_global_block_array;
           pinfo < _bftc_mem_global_block_array + _bftc_mem_global_block_nbr;
           pinfo++) {

        fprintf(_bftc_mem_global_file,"[%10p]\n", pinfo->p_bloc);
        non_free++;

      }

      fprintf(_bftc_mem_global_file,
              _("Number of non freed pointers remaining: %lu\n"),
              non_free);

    }

    fclose(_bftc_mem_global_file);
  }

  /* Reset defaults in case of later initialization */

  if (_bftc_mem_global_block_array != NULL) {
    free(_bftc_mem_global_block_array);
    _bftc_mem_global_block_array = NULL;
  }

  _bftc_mem_global_block_nbr   = 0 ;
  _bftc_mem_global_block_max   = 512 ;

  _bftc_mem_global_alloc_cur = 0;
  _bftc_mem_global_alloc_max = 0;

  _bftc_mem_global_n_allocs = 0;
  _bftc_mem_global_n_reallocs = 0;
  _bftc_mem_global_n_frees = 0;

}

/*!
 * \brief Indicates if bftc_mem_...() functions are initialized.
 *
 * \returns 1 if bftc_mem_init has been called, 0 otherwise.
 */

int
bftc_mem_initialized(void)
{
  return _bftc_mem_global_initialized;
}

/*!
 * \brief Allocate memory for ni elements of size bytes.
 *
 * This function calls malloc(), but adds tracing capabilities, and
 * automatically calls the bftc_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to allocated memory.
 */

void *
bftc_mem_malloc(size_t       ni,
               size_t       size,
               const char  *var_name,
               const char  *file_name,
               int          line_num)
{
  void       *p_loc;
  size_t      alloc_size = ni * size;

  if (ni == 0)
    return NULL;

  /* Allocate memory and check return */

  p_loc = malloc(alloc_size);

  if (p_loc == NULL) {
    _bftc_mem_error(file_name, line_num, errno,
                   _("Failure to allocate \"%s\" (%lu bytes)"),
                   var_name, (unsigned long)alloc_size);
    return NULL;
  }
  else if (_bftc_mem_global_initialized == 0)
    return p_loc;

  /* Memory allocation counting */

  _bftc_mem_global_alloc_cur += alloc_size;

  if (_bftc_mem_global_alloc_max < _bftc_mem_global_alloc_cur)
    _bftc_mem_global_alloc_max = _bftc_mem_global_alloc_cur;

  if (_bftc_mem_global_file != NULL) {
    fprintf(_bftc_mem_global_file, "\n  alloc: %-27s:%6d : %-39s: %9lu",
            _bftc_mem_basename(file_name), line_num,
            var_name, (unsigned long)alloc_size);
    fprintf(_bftc_mem_global_file, " : (+%9lu) : %12lu : [%10p]",
            (unsigned long)alloc_size,
            (unsigned long)_bftc_mem_global_alloc_cur,
            p_loc);
    fflush(_bftc_mem_global_file);
  }

  _bftc_mem_block_malloc(p_loc, alloc_size);

  _bftc_mem_global_n_allocs += 1;

  /* Return pointer to allocated memory */

  return p_loc;
}

/*!
 * \brief Reallocate memory for ni elements of size bytes.
 *
 * This function calls realloc(), but adds tracing capabilities, and
 * automatically calls the bftc_error() errorhandler if it fails to
 * allocate the required memory.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, bftc_alloc() called).
 * \param [in] ni        number of elements.
 * \param [in] size      element size.
 * \param [in] var_name  allocated variable name string.
 * \param [in] file_name name of calling source file.
 * \param [in] line_num  line number in calling source file.
 *
 * \returns pointer to reallocated memory.
 */

void *
bftc_mem_realloc(void        *ptr,
                size_t       ni,
                size_t       size,
                const char  *var_name,
                const char  *file_name,
                int          line_num)
{
  void      *p_loc;

  long size_diff;
  size_t old_size;
  size_t new_size = ni * size;

  /*
    Behave as bftc_malloc() if the previous pointer is equal to NULL.
    Note that the operation will then appear as a first allocation
    ('alloc') in the _bftc_mem_global_file trace file.
  */

  if (ptr == NULL)
    return bftc_mem_malloc(ni,
                          size,
                          var_name,
                          file_name,
                          line_num);

  /* If the old size equals the new size, nothing needs to be done. */

  old_size = _bftc_mem_block_size(ptr);

  if (new_size == old_size)
    return ptr;

  /*
    We may also simply free memory. Note that in this case, the operation
    appears as 'free' in the _bftc_mem_global_file trace file.
  */

  else if (ni == 0)
    return bftc_mem_free(ptr,
                        var_name,
                        file_name,
                        line_num);

  /* In the final case, we have a true reallocation */

  else {

    size_diff = new_size - old_size;

    p_loc = realloc(ptr, new_size);

    if (p_loc == NULL) {
      _bftc_mem_error(file_name, line_num, errno,
                     _("Failure to reallocate \"%s\" (%lu bytes)"),
                     var_name, (unsigned long)new_size);
      return NULL;
    }
    else if (_bftc_mem_global_initialized == 0)
      return p_loc;

    _bftc_mem_global_alloc_cur += size_diff;

    if (size_diff > 0) {
      if (_bftc_mem_global_alloc_max < _bftc_mem_global_alloc_cur)
        _bftc_mem_global_alloc_max = _bftc_mem_global_alloc_cur;
    }

    if (_bftc_mem_global_file != NULL) {
      char sgn = (char) ((size_diff > 0) ? '+' : '-');
      fprintf(_bftc_mem_global_file, "\nrealloc: %-27s:%6d : %-39s: %9lu",
              _bftc_mem_basename(file_name), line_num,
              var_name, (unsigned long)new_size);
      fprintf(_bftc_mem_global_file, " : (%c%9lu) : %12lu : [%10p]",
              sgn,
              (unsigned long) ((size_diff > 0) ? size_diff : -size_diff),
              (unsigned long)_bftc_mem_global_alloc_cur,
              p_loc);
      fflush(_bftc_mem_global_file);
    }

CWP_GCC_SUPPRESS_WARNING_WITH_PUSH("-Wuse-after-free")

    _bftc_mem_block_realloc(ptr, p_loc, new_size);

CWP_GCC_SUPPRESS_WARNING_POP

    _bftc_mem_global_n_reallocs += 1;

    return p_loc;
  }

}

/*!
 * \brief Free allocated memory.
 *
 * This function calls free(), but adds tracing capabilities, and
 * automatically calls the bftc_error() errorhandler if it fails to
 * free the corresponding memory. In case of a NULL pointer argument,
 * the function simply returns.
 *
 * \param [in] ptr       pointer to previous memory location
 *                       (if NULL, bftc_alloc() called).
 * \param [in] var_name  allocated variable name string
 * \param [in] file_name name of calling source file
 * \param [in] line_num  line number in calling source file
 *
 * \returns NULL pointer.
 */

void *
bftc_mem_free(void        *ptr,
             const char  *var_name,
             const char  *file_name,
             int          line_num)
{
  size_t  size_info;

  /* NULL pointer case (non-allocated location) */

  if (ptr == NULL)
    return NULL;

  /* General case (free allocated memory) */

  if (_bftc_mem_global_initialized != 0) {

    size_info = _bftc_mem_block_size(ptr);

    _bftc_mem_global_alloc_cur -= size_info;

    if (_bftc_mem_global_file != NULL) {
      fprintf(_bftc_mem_global_file,"\n   free: %-27s:%6d : %-39s: %9lu",
              _bftc_mem_basename(file_name), line_num,
              var_name, (unsigned long)size_info);
      fprintf(_bftc_mem_global_file, " : (-%9lu) : %12lu : [%10p]",
              (unsigned long)size_info,
              (unsigned long)_bftc_mem_global_alloc_cur,
              ptr);
      fflush(_bftc_mem_global_file);
    }

    _bftc_mem_block_free(ptr);

    _bftc_mem_global_n_frees += 1;
  }

  free(ptr);

  return NULL;
}

/*!
 * \brief Allocate aligned memory for ni elements of size bytes.
 *
 * This function calls posix_memalign() if available, but adds tracing
 * capabilities, and automatically calls the bftc_error() errorhandler if
 * it fails to allocate the required memory.
 *
 * The associated function bftc_mem_have_memalign() indicates if this
 * type of allocation may be used on this system.
 *
 * \param [in]  alignment alignment.
 * \param [in]  ni        number of elements.
 * \param [in]  size      element size.
 * \param [in]  var_name  allocated variable name string.
 * \param [in]  file_name name of calling source file.
 * \param [in]  line_num  line number in calling source file.
 *
 * \returns pointer to allocated memory.
 */

void *
bftc_mem_memalign(size_t       alignment,
                 size_t       ni,
                 size_t       size,
                 const char  *var_name,
                 const char  *file_name,
                 int          line_num)
{
#if defined(HAVE_POSIX_MEMALIGN)

  int         retval;
  void       *p_loc;
  size_t      alloc_size = ni * size;

  if (ni == 0)
    return NULL;

  /* Allocate memory and check return */

  BFTC_CLANG_SUPPRESS_WARNING_WITH_PUSH("-Wimplicit-function-declaration")
  retval = posix_memalign(&p_loc, alignment, alloc_size);
  BFTC_CLANG_SUPPRESS_WARNING_POP

  if (retval != 0) {
    switch (retval) {
    case EINVAL:
      _bftc_mem_error(file_name, line_num, 0,
                     _("Alignment %lu for \"%s\" not a power of 2\n"
                       "or a multiple of sizeof(void *) = %lu"),
                     (unsigned long)alignment, var_name,
                     (unsigned long)(sizeof(void *)));
      break;
    default:
      _bftc_mem_error(file_name, line_num, 0,
                     _("Failure to allocate \"%s\" (%lu bytes)"),
                     var_name, (unsigned long)alloc_size);
    }
    return NULL;
  }
  else if (_bftc_mem_global_initialized == 0)
    return p_loc;

  /* Memory allocation counting */

  _bftc_mem_global_alloc_cur += alloc_size;

  if (_bftc_mem_global_alloc_max < _bftc_mem_global_alloc_cur)
    _bftc_mem_global_alloc_max = _bftc_mem_global_alloc_cur;

  if (_bftc_mem_global_file != NULL) {
    fprintf(_bftc_mem_global_file, "\n  alloc: %-27s:%6d : %-39s: %9lu",
            _bftc_mem_basename(file_name), line_num,
            var_name, (unsigned long)alloc_size);
    fprintf(_bftc_mem_global_file, " : (+%9lu) : %12lu : [%10p]",
            (unsigned long)alloc_size,
            (unsigned long)_bftc_mem_global_alloc_cur,
            p_loc);
    fflush(_bftc_mem_global_file);
  }

  _bftc_mem_block_malloc(p_loc, alloc_size);

  _bftc_mem_global_n_allocs += 1;

  /* Return pointer to allocated memory */

  return p_loc;

#else

  _bftc_mem_error(file_name, line_num, errno,
                 _("No aligned allocation function available on this system"));

  return NULL;

#endif
}

/*!
 * \brief Return current theoretical dynamic memory allocated.
 *
 * \return current memory handled through bftc_mem_...() (in kB).
 */

size_t
bftc_mem_size_current(void)
{
  return (_bftc_mem_global_alloc_cur / 1024);
}

/*!
 * \brief Return maximum theoretical dynamic memory allocated.
 *
 * \return maximum memory handled through bftc_mem_...() (in kB).
 */

size_t
bftc_mem_size_max(void)
{
  return (_bftc_mem_global_alloc_max / 1024);
}

/*!
 * \brief Returns the error handler associated with the bftc_mem_...() functions.
 *
 * \return pointer to the error handler function.
 */

bftc_error_handler_t *
bftc_mem_error_handler_get(void)
{
  return _bftc_mem_error_handler;
}

/*!
 * \brief Associates an error handler with the bftc_mem_...() functions.
 *
 * With the default error handler, an error message is output to stderr,
 * (after bftc_print_flush() is called), and the general error handler used
 * by bftc_error() is then called (which results in the termination of the
 * current process or process group).
 *
 * \param handler pointer to the error handler function [in].
 */

void
bftc_mem_error_handler_set(bftc_error_handler_t *handler)
{
  _bftc_mem_error_handler = handler;
}

/*!
 * \brief Indicate if a memory aligned allocation variant is available.
 *
 * If no such function is available, bftc_mem_memalign() will always fail.
 *
 * \returns 1 if memory aligned allocation is possible, 0 otherwise.
 */

int
bftc_mem_have_memalign(void)
{
#if defined(HAVE_POSIX_MEMALIGN)
  return 1;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
