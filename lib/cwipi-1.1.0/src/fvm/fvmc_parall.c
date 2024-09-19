/*============================================================================
 * Base functions for parallelism
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2004-2008  EDF

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_config_defs.h"
#include "fvmc_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_parall.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Local structure definitions
 *============================================================================*/

/*============================================================================
 * Global variables
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

/* Basic communicator info */

static MPI_Comm  _fvmc_mpi_parall_comm = MPI_COMM_NULL;  /* Intra-communicator */
static int       _fvmc_mpi_parall_size = 1;
static int       _fvmc_mpi_parall_rank = 0;

/* MPI Datatypes associated with fvm datatypes */

MPI_Datatype  fvmc_datatype_to_mpi[] = {MPI_DATATYPE_NULL,
                                       MPI_CHAR,
                                       MPI_UNSIGNED_CHAR,
                                       MPI_FLOAT,
                                       MPI_DOUBLE,
                                       MPI_SHORT,          /* FVMC_INT16 */
                                       MPI_INT,            /* FVMC_INT32 */
                                       MPI_LONG_INT,       /* FVMC_INT64 */
                                       MPI_UNSIGNED_SHORT, /* FVMC_UINT16 */
                                       MPI_UNSIGNED,       /* FVMC_UINT32 */
                                       MPI_UNSIGNED_LONG}; /* FVMC_UINT64 */

/* Minimum recommended scatter/gather buffer size */

static size_t _fvmc_parall_min_coll_buf_size = 1024*1024*8;

#endif

/*=============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Return default MPI communicator for FVM library functions.
 *
 * returns:
 *   handle to MPI communicator
 *----------------------------------------------------------------------------*/

MPI_Comm
fvmc_parall_get_mpi_comm(void)
{
  return _fvmc_mpi_parall_comm;
}

/*----------------------------------------------------------------------------
 * Set default MPI communicator for FVM library functions.
 *
 * parameters:
 *   comm <-- handle to MPI communicator
 *----------------------------------------------------------------------------*/

void
fvmc_parall_set_mpi_comm(const MPI_Comm  comm)
{
  int mpi_flag;

  MPI_Initialized(&mpi_flag);

  /* Set communicator */

  _fvmc_mpi_parall_comm = comm;

  if (mpi_flag != 0 && _fvmc_mpi_parall_comm != MPI_COMM_NULL) {
    MPI_Comm_size(_fvmc_mpi_parall_comm, &_fvmc_mpi_parall_size);
    MPI_Comm_rank(_fvmc_mpi_parall_comm, &_fvmc_mpi_parall_rank);
  }
  else {
    _fvmc_mpi_parall_size = 1;
    _fvmc_mpi_parall_rank = 0;
  }

  /* Check (correct) fvmc_datatype_to_mpi values */

  if (mpi_flag != 0)
  {
    int size_short, size_int, size_long, size_long_long;

    MPI_Type_size(MPI_SHORT, &size_short);
    MPI_Type_size(MPI_INT,   &size_int);
    MPI_Type_size(MPI_LONG,  &size_long);

#if defined(MPI_LONG_LONG)
    MPI_Type_size(MPI_LONG_LONG, &size_long_long);
#else
    size_long_long = 0;
#endif

    fvmc_datatype_to_mpi[FVMC_CHAR] = MPI_CHAR;
    fvmc_datatype_to_mpi[FVMC_UCHAR] = MPI_UNSIGNED_CHAR;

    if (size_short == 2) {
      fvmc_datatype_to_mpi[FVMC_INT16] = MPI_SHORT;
      fvmc_datatype_to_mpi[FVMC_UINT16] = MPI_UNSIGNED_SHORT;
    }

    if (size_int == 4) {
      fvmc_datatype_to_mpi[FVMC_INT32] = MPI_INT;
      fvmc_datatype_to_mpi[FVMC_UINT32] = MPI_UNSIGNED;
    }
    else if (size_short == 4) {
      fvmc_datatype_to_mpi[FVMC_INT32] = MPI_SHORT;
      fvmc_datatype_to_mpi[FVMC_UINT32] = MPI_UNSIGNED_SHORT;
    }
    else if (size_long == 4) {
      fvmc_datatype_to_mpi[FVMC_INT32] = MPI_LONG;
      fvmc_datatype_to_mpi[FVMC_UINT32] = MPI_UNSIGNED_LONG;
    }

    if (size_int == 8) {
      fvmc_datatype_to_mpi[FVMC_INT64] = MPI_INT;
      fvmc_datatype_to_mpi[FVMC_UINT64] = MPI_UNSIGNED;
    }
    else if (size_long == 8) {
      fvmc_datatype_to_mpi[FVMC_INT64] = MPI_LONG;
      fvmc_datatype_to_mpi[FVMC_UINT64] = MPI_UNSIGNED_LONG;
    }
#if defined(MPI_LONG_LONG)
    else if (size_long_long == 8) {
      fvmc_datatype_to_mpi[FVMC_INT64] = MPI_LONG_LONG;
#if defined(MPI_UNSIGNED_LONG_LONG)
      fvmc_datatype_to_mpi[FVMC_UINT64] = MPI_UNSIGNED_LONG_LONG;
#else
      fvmc_datatype_to_mpi[FVMC_UINT64] = MPI_LONG_LONG;
#endif
    }
#endif
  }

}

#endif

/*----------------------------------------------------------------------------
 * Return rank of current process among associated program processes.
 *
 * returns:
 *   rank of current process in current communicator, or 0 in scalar mode
 *----------------------------------------------------------------------------*/

int
fvmc_parall_get_rank(void)
{
#if defined(FVMC_HAVE_MPI)
  return _fvmc_mpi_parall_rank;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------
 * Return number of processes associated with the current program.
 *
 * returns:
 *   number of processes in current communicator, or 1 in scalar mode
 *----------------------------------------------------------------------------*/

int
fvmc_parall_get_size(void)
{
#if defined(FVMC_HAVE_MPI)
  return _fvmc_mpi_parall_size;
#else
  return 1;
#endif
}

#if defined(FVMC_HAVE_MPI)

/*----------------------------------------------------------------------------
 * Sum counters on all FVM default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

void
fvmc_parall_counter(fvmc_gnum_t  cpt[],
                   const int   n)
{

  if (_fvmc_mpi_parall_size > 1) {

    int        i;
    fvmc_gnum_t *sum;
    fvmc_gnum_t _sum[64];

    if (n > 64)
      BFTC_MALLOC(sum, n, fvmc_gnum_t);
    else
      sum = _sum;

    MPI_Allreduce(cpt, sum, n, FVMC_MPI_GNUM, MPI_SUM,
                  _fvmc_mpi_parall_comm);

    for (i = 0; i < n ; i++)
      cpt[i] = sum[i];

    if (sum != _sum)
      BFTC_FREE(sum);

  }

}

/*----------------------------------------------------------------------------
 * Maximum values of a counter on all FVM default communicator processes.
 *
 * parameters:
 *   cpt <-> local counter value  input, global counter value output (array)
 *   n   <-- number of counter array values
 *----------------------------------------------------------------------------*/

void
fvmc_parall_counter_max(fvmc_lnum_t  cpt[],
                       const int   n)
{

  if (_fvmc_mpi_parall_size > 1) {

    int        i;
    fvmc_lnum_t *maxval;
    fvmc_lnum_t _maxval[64];

    if (n > 64)
      BFTC_MALLOC(maxval, n, fvmc_lnum_t);
    else
      maxval = _maxval;

    MPI_Allreduce(cpt, maxval, n, FVMC_MPI_LNUM, MPI_MAX,
                  _fvmc_mpi_parall_comm);

    for (i = 0; i < n ; i++)
      cpt[i] = maxval[i];

    if (maxval != _maxval)
      BFTC_FREE(maxval);

  }

}

#endif /* defined(FVMC_HAVE_MPI) */

/*----------------------------------------------------------------------------
 * Return minimum recommended scatter or gather buffer size.
 *
 * This is used by FVM's internal strided and indexed array scatter/gather
 * algorithms, for non MPI-IO Input/output.
 *
 * returns:
 *   minimum recommended gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

size_t
fvmc_parall_get_min_coll_buf_size(void)
{
#if defined(FVMC_HAVE_MPI)
  return _fvmc_parall_min_coll_buf_size;
#else
  return 0;
#endif
}

/*----------------------------------------------------------------------------
 * Define minimum recommended gather buffer size.
 *
 * This is used by FVM's internal strided and indexed array scatter/gather
 * algorithms, for non MPI-IO Input/output.
 *
 * parameters:
 *   minimum recommended gather buffer size (in bytes)
 *----------------------------------------------------------------------------*/

void
fvmc_parall_set_min_coll_buf_size(size_t buffer_size)
{
#if defined(FVMC_HAVE_MPI)
  _fvmc_parall_min_coll_buf_size = buffer_size;
#endif
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
