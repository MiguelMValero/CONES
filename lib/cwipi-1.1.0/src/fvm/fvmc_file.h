#ifndef __FVMC_FILE_H__
#define __FVMC_FILE_H__

/*============================================================================
 * Parallel file I/O
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2007-2008  EDF

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

#include "fvmc_config.h"

#if defined(FVMC_HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*
 * File hints and semantics
 */

#define FVMC_FILE_NO_MPI_IO            (1 << 0)
#define FVMC_FILE_NO_PREDISTRIBUTE     (1 << 1)

/* MPI-IO positioning semantics */

#define FVMC_FILE_EXPLICIT_OFFSETS     (1 << 2)
#define FVMC_FILE_INDIVIDUAL_POINTERS  (1 << 3)

/*============================================================================
 * Type definitions
 *============================================================================*/

/* FVM file descriptor */

typedef struct _fvmc_file_t  fvmc_file_t;

/* FVM file modes */

typedef enum {

  FVMC_FILE_MODE_READ,   /* Read mode */
  FVMC_FILE_MODE_WRITE,  /* Write mode */
  FVMC_FILE_MODE_APPEND  /* Append mode */

} fvmc_file_mode_t;

/* Hints for file management */

typedef unsigned int fvmc_file_hints_t;

/* Offset for FVM file position indicator (int64_t in C99) */

#if defined(FVMC_SIZEOF_LONG_LONG)
typedef long long fvmc_file_off_t;
#else
typedef long fvmc_file_off_t;
#endif

/* Possibilities for the third argument of fvmc_file_seek() */

typedef enum {

  FVMC_FILE_SEEK_SET,   /* Seek from beginning of file */
  FVMC_FILE_SEEK_CUR,   /* Seek from current position */
  FVMC_FILE_SEEK_END    /* Seek from end of file */

} fvmc_file_seek_t;

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Create a file descriptor and open the associated file.
 *
 * By default, data is written or read as native data. This behavior may be
 * modified by fvmc_file_set_swap_endian().
 *
 * parameters:
 *   name  <-- file name
 *   mode  <-- file acces mode: read, write, or append
 *   hints <-- file I/O hints (for MPI and MPI I/O behavior)
 *
 * returns:
 *   pointer to fvmc_file_t file descriptor (NULL in case of failure);
 *   currently, errors are fatal.
 *----------------------------------------------------------------------------*/

#if defined(FVMC_HAVE_MPI)

fvmc_file_t *
fvmc_file_open(const char         *name,
              fvmc_file_mode_t     mode,
              fvmc_file_hints_t    hints,
              MPI_Comm            comm);

#else

fvmc_file_t *
fvmc_file_open(const char         *name,
              fvmc_file_mode_t     mode,
              fvmc_file_hints_t    hints);

#endif

/*----------------------------------------------------------------------------
 * Destroy a file descriptor and close the associated file.
 *
 * parameters:
 *   f <-> file descriptor to destroy
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

fvmc_file_t *
fvmc_file_free(fvmc_file_t  *f);

/*----------------------------------------------------------------------------
 * Return a file's name.
 *
 * parameters:
 *   f <-- fvmc_file_t descriptor
 *
 * returns:
 *   pointer to the file's name.
 *----------------------------------------------------------------------------*/

const char *
fvmc_file_get_name(const fvmc_file_t  *f);

/*----------------------------------------------------------------------------
 * Ensure that data is read or written in big-endian
 * (network standard) format.
 *
 * parameters:
 *   f <-> fvmc_file_t descriptor
 *----------------------------------------------------------------------------*/

void
fvmc_file_set_big_endian(fvmc_file_t  *f);

/*----------------------------------------------------------------------------
 * Return a file's byte-swapping behavior.
 *
 * parameters:
 *   f <-- fvmc_file_t descriptor
 *
 * returns:
 *   0 if file's endianness is the same as the system's, 1 otherwise.
 *----------------------------------------------------------------------------*/

int
fvmc_file_get_swap_endian(const fvmc_file_t  *f);

/*----------------------------------------------------------------------------
 * Set a file's byte-swapping behavior.
 *
 * Using this function assumes one is familiar with a file's coding
 * or structure; use with caution.
 *
 * parameters:
 *   f    <-> fvmc_file_t descriptor
 *   swap --> 1 if bytes must be swapped, 0 otherwise
 *----------------------------------------------------------------------------*/

void
fvmc_file_set_swap_endian(fvmc_file_t  *f,
                         int          swap);

/*----------------------------------------------------------------------------
 * Read data to a buffer, distributing it to all processes associated
 * with a file.
 *
 * parameters:
 *   f    <-- fvmc_file_t descriptor
 *   buf  --> pointer to location receiving data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the number of items (not bytes) sucessfully read; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvmc_file_read_global(fvmc_file_t  *f,
                     void        *buf,
                     size_t       size,
                     size_t       ni);

/*----------------------------------------------------------------------------
 * Write global data to a file.
 *
 * Under MPI, data is only written by the associated communicator's root
 * rank. The buffers on other ranks are ignored, though the file offset
 * is updated (i.e. the call to this function is collective).
 *
 * parameters:
 *   f    <-- fvmc_file_t descriptor
 *   buf  <-- pointer to location containing data
 *   size <-- size of each item of data in bytes
 *   ni   <-- number of items to read
 *
 * returns:
 *   the number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvmc_file_write_global(fvmc_file_t  *f,
                      const void  *buf,
                      size_t       size,
                      size_t       ni);

/*----------------------------------------------------------------------------
 * Read data to a buffer, distributing a contiguous part of it to each
 * process associated with a file.
 *
 * Each process should receive a (possibly empty) block of the data,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * parameters:
 *   f                <-- fvmc_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   stride           <-- number of (interlaced) values per block item
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully read; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvmc_file_read_block(fvmc_file_t  *f,
                    void        *buf,
                    size_t       size,
                    size_t       stride,
                    fvmc_gnum_t   global_num_start,
                    fvmc_gnum_t   global_num_end);

/*----------------------------------------------------------------------------
 * Write data to a file, each associated process providing a contiguous part
 * of this data.
 *
 * Each process should provide a (possibly empty) block of the data,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * This function may require an internal copy of the data to ensure that
 * the buffer contents are not modified, so if the buffer contents are
 * temporary values, to be deleted after writing, using
 * fvmc_file_write_block_buffer() instead may be used to avoid an unneeded
 * memory allocation and copy.
 *
 * parameters:
 *   f                <-- fvmc_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   stride           <-- number of (interlaced) values per block item
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvmc_file_write_block(fvmc_file_t  *f,
                     void        *buf,
                     size_t       size,
                     size_t       stride,
                     fvmc_gnum_t   global_num_start,
                     fvmc_gnum_t   global_num_end);

/*----------------------------------------------------------------------------
 * Write data to a file, each associated process providing a contiguous part
 * of this data.
 *
 * Each process should provide a (possibly empty) block of the data,
 * and we should have:
 *   global_num_start at rank 0 = 1
 *   global_num_start at rank i+1 = global_num_end at rank i.
 * Otherwise, behavior (especially positioning for future reads) is undefined.
 *
 * This function is intended to be used mainly data that is already of
 * copy of original data (such as data that has been redistributed across
 * processors just for the sake of output), or that is to be deleted after
 * writing, so it may modify the values in its input buffer (notably to
 * convert from little-endian to big-endian of vice-versa if necessary).
 *
 * parameters:
 *   f                <-- fvmc_file_t descriptor
 *   buf              --> pointer to location receiving data
 *   size             <-- size of each item of data in bytes
 *   stride           <-- number of (interlaced) values per block item
 *   global_num_start <-- global number of first block item (1 to n numbering)
 *   global_num_end   <-- global number of past-the end block item
 *                        (1 to n numbering)
 *
 * returns:
 *   the (local) number of items (not bytes) sucessfully written; currently,
 *   errors are fatal.
 *----------------------------------------------------------------------------*/

size_t
fvmc_file_write_block_buffer(fvmc_file_t  *f,
                            void        *buf,
                            size_t       size,
                            size_t       stride,
                            fvmc_gnum_t   global_num_start,
                            fvmc_gnum_t   global_num_end);

/*----------------------------------------------------------------------------
 * Update the file pointer according to whence.
 *
 * parameters:
 *   f      <-- fvmc_file_t descriptor.
 *   offset <-- add to position specified to whence to obtain new position,
 *              measured in characters from the beginning of the file.
 *   whence <-- beginning if FVMC_FILE_SEEK_SET, current if FVMC_FILE_SEEK_CUR,
 *               or end-of-file if FVMC_FILE_SEEK_END.
 *
 * returns:
 *   0 upon success, nonzero otherwise; currently, errors are fatal.
 *----------------------------------------------------------------------------*/

int
fvmc_file_seek(fvmc_file_t       *f,
              fvmc_file_off_t    offset,
              fvmc_file_seek_t   whence);

/*----------------------------------------------------------------------------
 * Return the position of the file pointer.
 *
 * When using MPI-IO with individual file pointers, we consider the file
 * pointer to be equal to the highest value of then individual file pointers.
 *
 * parameters:
 *   f <-- fvmc_file_t descriptor
 *
 * returns:
 *   current position of the file pointer
 *----------------------------------------------------------------------------*/

fvmc_file_off_t
fvmc_file_tell(fvmc_file_t  *f);

/*----------------------------------------------------------------------------
 * Get the default semantics for file access.
 *
 * returns:
 *   current default semantics for file access
 *----------------------------------------------------------------------------*/

fvmc_file_hints_t
fvmc_file_get_default_semantics(void);

/*----------------------------------------------------------------------------
 * Set the default semantics for file access.
 *
 * This may fail if semantics given contain incompatible values,
 * such as (FVMC_FILE_EXPLICIT_OFFSETS | FVMC_FILE_INDIVIDUAL_POINTERS),
 * or when setting MPI-IO access semantics when MPI-IO is not available.
 *
 * returns:
 *   0 if the semantics were valid, 1 otherwise.
 *----------------------------------------------------------------------------*/

int
fvmc_file_set_default_semantics(fvmc_file_hints_t  hints);

/*----------------------------------------------------------------------------
 * Dump the metadata of a file structure in human readable form
 *
 * parameters:
 *   f <-- pointer to file
 *----------------------------------------------------------------------------*/

void
fvmc_file_dump(const fvmc_file_t  *f);

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FVMC_FILE_H__ */
