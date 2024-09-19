/*============================================================================
 * Set up communication with coupled codes.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2008  EDF

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
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * BFT library headers
 *----------------------------------------------------------------------------*/

#include <bftc_mem.h>
#include <bftc_printf.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_config_defs.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_coupling.h"

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

/*=============================================================================
 * Local type definitions
 *============================================================================*/

/* Structure used manage information about coupling with MPI_COMM_WORLD */

#if defined(FVMC_HAVE_MPI)

struct _fvmc_coupling_mpi_world_t {

  int    n_apps;       /* Number of distinct applications */
  int    app_id;       /* Id of the local application in the application info */
  int    app_names_l;  /* Length of application names array */

  int   *app_info;     /* For each application, 5 integers: application number,
                          associated root, n_ranks, and indexes in app_names */
  char  *app_names;    /* Array of application type names and instance names */

};

#endif /* defined(FVMC_HAVE_MPI) */

/*============================================================================
 * Local function defintions
 *============================================================================*/

#if defined(FVMC_HAVE_MPI)

#endif /* defined(FVMC_HAVE_MPI) */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Discover other applications in the same MPI_COMM_WORLD.
 *
 * The application communicator app_comm is usually obtained from
 * MPI_COMM_WORLD using MPI_Comm_split, with app_num corresponding to
 * the "color" argument in that function.
 *
 * As this function requires communication between applications, it
 * is a collective function in MPI_COMM_WORLD.
 *
 * parameters:
 *   app_num   <-- application number in MPI_COMM_WORLD (nonnegative).
 *   app_name  <-- name of current application.
 *   case_name <-- name of current case, or NULL.
 *   app_comm  <-- communicator associated with local application.
 *
 * returns:
 *   FVM coupling MPI_COMM_WORLD info structure.
 *----------------------------------------------------------------------------*/

fvmc_coupling_mpi_world_t *
fvmc_coupling_mpi_world_create(int          app_num,
                              const char  *app_type,
                              const char  *app_name,
                              MPI_Comm     app_comm)
{
  int i;
  MPI_Status status;

  int world_rank = -1, app_rank = -1, n_app_ranks = 0;
  int root_marker = 0;
  int info_tag = 1, name_tag = 2;

  int counter[2] = {0, 0};
  int l_rank_info[5] = {-1, -1, -1, 1, 1};
  int *rank_info = NULL;
  char *app_names = NULL;

  fvmc_coupling_mpi_world_t *w = NULL;

  /* Initialization */

  BFTC_MALLOC(w, 1, fvmc_coupling_mpi_world_t);

  w->n_apps = 0;
  w->app_id = -1;
  w->app_names_l = 0;
  w->app_info = NULL;
  w->app_names = NULL;

  /* Initialization */

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (app_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(app_comm, &app_rank);
    MPI_Comm_size(app_comm, &n_app_ranks);
  }
  else {
    app_rank = 0;
    n_app_ranks = 1;
  }

  l_rank_info[0] = app_num;
  l_rank_info[1] = world_rank;
  l_rank_info[2] = n_app_ranks;
  if (app_type != NULL)
    l_rank_info[3] = strlen(app_type) + 1;
  if (app_name != NULL)
    l_rank_info[4] = strlen(app_name) + 1;

  if (app_rank == 0)
    root_marker = 1;

  /* Root rank of MPI_COMM_WORLD counts applications and receives info */

  MPI_Reduce(&root_marker, &(counter[0]), 1, MPI_INT, MPI_SUM, 0,
             MPI_COMM_WORLD);

  /* Root of MPI_COMM_WORLD collects all info */

  if (world_rank == 0) {

    int start = 0;

    BFTC_MALLOC(rank_info, counter[0]*5, int);

    if (app_rank == 0) {
      for (i = 0; i < 5; i++)
        rank_info[i] = l_rank_info[i];
      start = 1;
    }

    /* Use of different tags for info and strings is important
       here as we use MPI_ANY_SOURCE and messages could be mixed */

    for (i = start; i < counter[0]; i++)
      MPI_Recv(rank_info + i*5, 5, MPI_INT, MPI_ANY_SOURCE, info_tag,
               MPI_COMM_WORLD, &status);

    /* Convert rank_info count to index values */

    for (i = 0; i < counter[0]; i++)
      counter[1] += (rank_info[i*5 + 3] + rank_info[i*5 + 4]);

    BFTC_MALLOC(app_names, counter[1], char);
    memset(app_names, 0, counter[1]);

    counter[1] = 0;

    if (app_rank == 0) {
      strcpy(app_names, app_type);
      if (app_name != NULL)
        strcpy(app_names + rank_info[3], app_name);
      else
        app_names[rank_info[3]] = '\0';
      counter[1] += (rank_info[3] + rank_info[4]);
      rank_info[4] = rank_info[3];
      rank_info[3] = 0;
    }

    for (i = start; i < counter[0]; i++) {
      int app_type_size = rank_info[i*5 + 3];
      int app_name_size = rank_info[i*5 + 4];
      int msg_len = app_type_size + app_name_size;
      rank_info[i*5 + 3] = counter[1];
      rank_info[i*5 + 4] = counter[1] + app_type_size;
      MPI_Recv(app_names + counter[1], msg_len, MPI_CHAR, rank_info[i*5 +1],
               name_tag, MPI_COMM_WORLD, &status);
      counter[1] += msg_len;
    }

  }

  /* Other root ranks send info to root */

  else if (app_rank == 0) { /* world_rank != 0 */

    char *sendbuf = NULL;
    int   sendbuf_l = l_rank_info[3] + l_rank_info[4];

    BFTC_MALLOC(sendbuf, sendbuf_l, char);

    if (app_type != NULL)
      strcpy(sendbuf, app_type);
    else
      sendbuf[0] = '\0';
    if (app_name != NULL)
      strcpy(sendbuf + l_rank_info[3], app_name);
    else
      sendbuf[l_rank_info[3]] = '\0';

    MPI_Send(l_rank_info, 5, MPI_INT, 0, info_tag, MPI_COMM_WORLD);
    MPI_Send(sendbuf, sendbuf_l, MPI_CHAR, 0, name_tag, MPI_COMM_WORLD);

    BFTC_FREE(sendbuf);
  }

  /* Now root broadcasts application info */

  MPI_Bcast(counter, 2, MPI_INT, 0, MPI_COMM_WORLD);

  if (world_rank != 0) {
    BFTC_MALLOC(rank_info, counter[0]*5, int);
    BFTC_MALLOC(app_names, counter[1], char);
  }

  MPI_Bcast(rank_info, counter[0]*5, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(app_names, counter[1], MPI_CHAR, 0, MPI_COMM_WORLD);

  /* Set global values */

  w->n_apps = counter[0];
  w->app_names_l = counter[1];
  w->app_info = rank_info;
  w->app_names = app_names;

  for (i = 0; i < w->n_apps && w->app_id < 0; i++) {
    if (w->app_info[i*5] == app_num)
      w->app_id = i;
  }

  return w;
}

/*----------------------------------------------------------------------------
 * Free an FVM coupling MPI_COMM_WORLD info structure.
 *
 * parameters:
 *   w <-> pointer to structure that should be freed.
 *----------------------------------------------------------------------------*/

void
fvmc_coupling_mpi_world_destroy(fvmc_coupling_mpi_world_t **w)
{
  fvmc_coupling_mpi_world_t *_w = *w;

  if (_w != NULL) {
    BFTC_FREE(_w->app_info);
    BFTC_FREE(_w->app_names);
    BFTC_FREE(*w);
  }
}

/*----------------------------------------------------------------------------
 * Return the number of applications in MPI_COMM_WORLD.
 *
 * parameters:
 *   w <-- pointer to FVM coupling MPI_COMM_WORLD info structure.
 *
 * returns:
 *   number of application in MPI_COMM_WORLD.
 *----------------------------------------------------------------------------*/

int
fvmc_coupling_mpi_world_n_apps(const fvmc_coupling_mpi_world_t  *w)
{
  int retval = 0;

  if (w != NULL)
    retval = w->n_apps;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return the id of the local application in MPI_COMM_WORLD.
 *
 * parameters:
 *   w <-- pointer to FVM coupling MPI_COMM_WORLD info structure.
 *
 * returns:
 *   id of the local application in MPI_COMM_WORLD.
 *----------------------------------------------------------------------------*/

int
fvmc_coupling_mpi_world_get_app_id(const fvmc_coupling_mpi_world_t  *w)
{
  int retval = -1;

  if (w != NULL)
    retval = w->app_id;

  return retval;
}

/*----------------------------------------------------------------------------
 * Return application information in MPI_COMM_WORLD.
 *
 * parameters:
 *   w      <-- pointer to FVM coupling MPI_COMM_WORLD info structure.
 *   app_id <-- application id
 *
 * returns:
 *   application information structure.
 *----------------------------------------------------------------------------*/

fvmc_coupling_mpi_world_info_t
fvmc_coupling_mpi_world_get_info(const fvmc_coupling_mpi_world_t  *w,
                                int                              app_id)
{
  fvmc_coupling_mpi_world_info_t  retval;

  retval.app_num = -1;
  retval.root_rank = -1;
  retval.n_ranks = 0;
  retval.app_type = NULL;
  retval.app_name = NULL;

  if (w != NULL) {
    if (app_id >= 0 && app_id < w->n_apps) {
      retval.app_num = w->app_info[app_id*5];
      retval.root_rank = w->app_info[app_id*5 + 1];
      retval.n_ranks = w->app_info[app_id*5 + 2];
      retval.app_type = w->app_names + w->app_info[app_id*5 + 3];
      retval.app_name = w->app_names + w->app_info[app_id*5 + 4];
    }
  }

  return retval;
}

/*----------------------------------------------------------------------------
 * Create an intracommunicator from a local and distant communicator
 * within MPI_COMM_WORLD.
 *
 * parameters:
 *   app_comm      <-- communicator associated with local application
 *   distant_root  <-- rank of distant group leader in MPI_COMM_WORLD
 *   new_comm      --> pointer to new communicator
 *   local_range   --> first and past-the last ranks of local application
 *                     in new communicator
 *   distant_range --> first and past-the last ranks of distant application
 *                     in new communicator
 *----------------------------------------------------------------------------*/

void
fvmc_coupling_mpi_intracomm_create(MPI_Comm   app_comm,
                                  int        distant_root,
                                  MPI_Comm  *new_comm,
                                  int        local_range[2],
                                  int        distant_range[2])
{
  int coupling_tag = ('F'+'V'+'M'+'_'+'C'+'O'+'U'+'P'+'L'+'I'+'N'+'G') % 512;
  int  mpi_flag = 0;
  int  n_dist_ranks = 0;
  int  n_loc_ranks, r_glob, r_loc_max;
  MPI_Comm  intercomm_tmp;
  int  r_coupl, r_coupl_min;
  int  high = 1;

  /* Initialization */

  *new_comm = MPI_COMM_NULL;

  MPI_Initialized(&mpi_flag);

  if (mpi_flag == 0)
    return;

  MPI_Comm_rank(MPI_COMM_WORLD, &r_glob);

  MPI_Allreduce(&r_glob, &r_loc_max, 1, MPI_INT, MPI_MAX, app_comm);

  if (distant_root > r_loc_max)
    high = 0;

  MPI_Comm_size(app_comm, &n_loc_ranks);

  /* Create a reserved communicator */

  MPI_Intercomm_create(app_comm, 0, MPI_COMM_WORLD,
                       distant_root, coupling_tag, &intercomm_tmp);

  MPI_Intercomm_merge(intercomm_tmp, high, new_comm);

  MPI_Comm_free(&intercomm_tmp);

  /* Compute number of distant ranks and first distant rank */

  MPI_Comm_size(*new_comm, &n_dist_ranks);
  n_dist_ranks -= n_loc_ranks;

  /* Check rank in new communicator (should not be necessary with correctly
     set "high" value, but seems to be with Open MPI 1.0.1) */

  MPI_Comm_rank(*new_comm, &r_coupl);
  MPI_Allreduce(&r_coupl, &r_coupl_min, 1, MPI_INT, MPI_MIN, app_comm);
  high = (r_coupl_min == 0) ? 0 : 1;

  /* Deduce the position of the first distant rank in the new communicator */

  if (high == 0) {
    local_range[0] = 0;
    distant_range[0] = n_loc_ranks;
  }
  else {
    local_range[0] = n_dist_ranks;
    distant_range[0] = 0;
  }

  local_range[1] = local_range[0] + n_loc_ranks;
  distant_range[1] = distant_range[0] + n_dist_ranks;
}

/*----------------------------------------------------------------------------
 * Dump printout of an FVM coupling MPI_COMM_WORLD info structure.
 *
 * parameters:
 *   w <-- pointer to FVM coupling MPI_COMM_WORLD info structure.
 *----------------------------------------------------------------------------*/

void
fvmc_coupling_mpi_world_dump(const fvmc_coupling_mpi_world_t  *w)
{
  int i;

  if (w == NULL) {
    bftc_printf("  Coupling MPI_COMM_WORLD info: nil\n");
    return;
  }

  bftc_printf("  Coupling MPI_COMM_WORLD info: %p\n"
             "    number of applications:     %d\n"
             "    local application id:       %d\n"
             "    app_names_size:             %d\n\n",
             w, w->n_apps, w->app_id, w->app_names_l);

  for (i = 0; i < w->n_apps; i++)
    bftc_printf("    Application number:  %d\n"
               "      root_rank:         %d\n"
               "      n_ranks:           %d\n"
               "      app_type:          \"%s\"\n"
               "      app_name:          \"%s\"\n\n",
               w->app_info[i*5], w->app_info[i*5 + 1], w->app_info[i*5 + 2],
               w->app_names + w->app_info[i*5 + 3],
               w->app_names + w->app_info[i*5 + 4]);
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

