#ifndef __PDM_MPI_PRIV_H__
#define __PDM_MPI_PRIV_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation csback to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

struct _pdm_mpi_win_shared_t {
  MPI_Win  win;
  void    *ptr;
};

struct _pdm_mpi_win_t {
  MPI_Win  win;
  void    *ptr;
};


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MPI_PRIV_H__ */
