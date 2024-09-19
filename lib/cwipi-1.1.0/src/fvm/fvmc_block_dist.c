/*============================================================================
 * Utility functions for block distribution.
 *============================================================================*/

/*
  This file is part of the "Finite Volume Mesh" library, intended to provide
  finite volume mesh and associated fields I/O and manipulation services.

  Copyright (C) 2011       ONERA

  Copyright (C) 2008-2009  EDF

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
#include <bftc_error.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "fvmc_defs.h"
#include "fvmc_config_defs.h"
#include "fvmc_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvmc_block_dist.h"

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

/*============================================================================
 * Local function defintions
 *============================================================================*/

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Compute block size and rank info for use with a block distribution.
 *
 * arguments:
 *   rank_id        <-- id of local rank (ignored in serial mode)
 *   n_ranks        <-- number of associated ranks
 *   min_rank_step  <-- minimum rank step between blocks
 *   min_block_size <-- minimum number of entities per block
 *   n_g_ents       <-- total number of associated entities
 *
 * returns:
 *   block size and range info structure
 *----------------------------------------------------------------------------*/

fvmc_block_dist_info_t
fvmc_block_dist_compute_sizes(int         rank_id,
                             int         n_ranks,
                             int         min_rank_step,
                             fvmc_lnum_t  min_block_size,
                             fvmc_gnum_t  n_g_ents)
{
  int _rank_id = rank_id;
  fvmc_gnum_t _min_block_size = 1;

  fvmc_block_dist_info_t bi;

  /* Special case: only 1 rank */

  if (n_ranks == 1) {

    bi.gnum_range[0] = 1;
    bi.gnum_range[1] = n_g_ents + 1;
    bi.n_ranks = 1;
    bi.rank_step = 1;
    bi.block_size = n_g_ents;

    return bi;
  }

  /* Determine rank stepping if necessary */

  assert(rank_id > -1);

  if (min_block_size > 1)
    _min_block_size = min_block_size;

  bi.block_size = 0;
  bi.n_ranks = n_ranks;
  bi.rank_step = 1;

  while (   n_g_ents/bi.n_ranks < _min_block_size
         && bi.n_ranks > 1
         && bi.rank_step < n_ranks) {
    bi.rank_step *= 2;
    bi.n_ranks = n_ranks / bi.rank_step;
  }
  if (bi.rank_step < min_rank_step)
    bi.rank_step = min_rank_step;
  if (bi.rank_step > n_ranks) {
    bi.rank_step = n_ranks;
    bi.n_ranks = 1;
  }

  if (rank_id % bi.rank_step == 0)
    _rank_id = rank_id/bi.rank_step;          /* non-empty block */
  else
    _rank_id = - (rank_id/bi.rank_step + 1);  /* empty block on this rank */

  /* Now determine block size and local range */

  bi.block_size = n_g_ents / bi.n_ranks;

  if (n_g_ents % bi.n_ranks)
    bi.block_size += 1;

  if (_rank_id > -1) {
    int i;
    for (i = 0; i < 2; i++) {
      bi.gnum_range[i] = (_rank_id+i)*bi.block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }
  else {
    int i;
    for (i = 0; i < 2; i++) {
      bi.gnum_range[i] = (-_rank_id)*bi.block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }

  return bi;
}

/*----------------------------------------------------------------------------
 * Compute block size and rank info for use with a block distribution
 * for a new global number of entities with a given number of active
 * ranks.
 *
 * arguments:
 *   rank_id        <-- id of local rank (ignored in serial mode)
 *   n_ranks        <-- number of associated ranks
 *   n_block_ranks  <-- number of ranks associated with a block
 *   n_g_ents       <-- total number of associated entities
 *
 * returns:
 *   block size and range info structure
 *----------------------------------------------------------------------------*/

fvmc_block_dist_info_t
fvmc_block_dist_compute_sizes_nr(int         rank_id,
                                int         n_ranks,
                                int         n_block_ranks,
                                fvmc_gnum_t  n_g_ents)
{
  int _rank_id = rank_id;

  fvmc_block_dist_info_t bi;

  /* Special case: only 1 rank */

  if (n_ranks == 1) {

    bi.gnum_range[0] = 1;
    bi.gnum_range[1] = n_g_ents + 1;
    bi.n_ranks = 1;
    bi.rank_step = 1;
    bi.block_size = n_g_ents;

    return bi;
  }

  /* Determine rank stepping if necessary */

  assert(rank_id > -1);

  bi.block_size = 0;
  bi.n_ranks = n_block_ranks;
  bi.rank_step = n_ranks / n_block_ranks;

  if (n_block_ranks < 1 || bi.rank_step > n_ranks) {
    bi.rank_step = n_ranks;
    bi.n_ranks = 1;
  }
  else if (bi.rank_step < 1) {
    bi.rank_step = 1;
    bi.n_ranks = n_ranks;
  }

  if (rank_id % bi.rank_step == 0)
    _rank_id = rank_id/bi.rank_step;          /* non-empty block */
  else
    _rank_id = - (rank_id/bi.rank_step + 1);  /* empty block on this rank */

  /* Now determine block size and local range */

  bi.block_size = n_g_ents / bi.n_ranks;

  if (n_g_ents % bi.n_ranks)
    bi.block_size += 1;

  if (_rank_id > -1) {
    int i;
    for (i = 0; i < 2; i++) {
      bi.gnum_range[i] = (_rank_id+i)*bi.block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }
  else {
    int i;
    for (i = 0; i < 2; i++) {
      bi.gnum_range[i] = (-_rank_id)*bi.block_size + 1;
      if (bi.gnum_range[i] > n_g_ents + 1)
        bi.gnum_range[i] = n_g_ents + 1;
    }
  }

  return bi;
}

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */

