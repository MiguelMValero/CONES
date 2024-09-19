/*
  This file is part of the CWIPI library.

  Copyright (C) 2011  ONERA

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <unistd.h>

#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -h              this message.\n\n");

  exit(exit_code);
}

/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   verbose   <--   Output logging
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                   argc,
  char                **argv,
  int                  *verbose
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int
main(int argc, char *argv[]) {

  int verbose  = 0;

  _read_args(argc,
             argv,
             &verbose);

  // Initialize MPI
  MPI_Init(&argc, &argv);
  int i_rank;
  int n_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &n_rank);

  // Initialize CWIPI
  int n_code = 1;
  const char  **code_name      = malloc(sizeof(char *) * n_code);
  CWP_Status_t  is_active_rank = CWP_STATUS_ON;
  MPI_Comm     *intra_comm     = malloc(sizeof(MPI_Comm) * n_code);

  int I_am_code1 = 0;
  int I_am_code2 = 0;

  if (i_rank == 0) {
    code_name[0] = "code1";
    I_am_code1   = 1;
  }

  if (i_rank == 1) {
    code_name[0] = "code2";
    I_am_code2   = 1;
  }
  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);

  // Dump properties
  if (I_am_code1) {
    // CWP_Properties_dump();
  }

  if (I_am_code2) {
    CWP_Properties_dump();
  }

  // Add parameters
  int toto = 5;
  const char *tata = "Bonjour code 1!";
  if (I_am_code1) {
    CWP_Param_lock("code1");
    CWP_Param_add("code1", "toto", CWP_INT, &toto);
    CWP_Param_add("code1", "tata", CWP_CHAR, &tata);
    CWP_Param_unlock("code1");
  }

  toto = 12;
  if (I_am_code2) {
    CWP_Param_lock("code2");
    CWP_Param_add("code2", "toto", CWP_INT, &toto);
    CWP_Param_unlock("code2");
  }

  // Get number of parameters
  int n_param = CWP_Param_n_get("code1", CWP_INT);
  if (verbose) {
    log_trace("n_param code1 : %d\n", n_param);
  }

  // Get parameter
  int get_toto = 0;
  if (I_am_code1) {
    CWP_Param_get("code1", "toto", CWP_INT, &get_toto);
    if (verbose) {
      log_trace("toto code1 : %d\n", get_toto);
    }
  }

  char *get_tata = NULL;
  CWP_Param_get("code1", "tata", CWP_CHAR, &get_tata);
  if (verbose) {
    log_trace("tata code1 : %s\n", get_tata);
  }

  // Get parameter list
  char **param_list = NULL;
  CWP_Param_list_get("code1", CWP_INT, &n_param, &param_list);
  if (verbose) {
    for (int i = 0; i < n_param; i++) {
      log_trace("param_list[%d] : %s\n", i, param_list[i]);
    }
  }

  // Code information
  int          n_codes   = CWP_Codes_nb_get();
  const char **code_list = CWP_Codes_list_get();
  if (verbose) {
    for (int i = 0; i < n_codes; i++) {
      log_trace("code_list[%d] : %s\n", i, code_list[i]);
    }
  }

  int          n_loc_codes   = CWP_Loc_codes_nb_get();
  const char **loc_code_list = CWP_Loc_codes_list_get();
  if (verbose) {
    for (int i = 0; i < n_loc_codes; i++) {
      log_trace("loc_code_list[%d] : %s\n", i, loc_code_list[i]);
    }
  }

  // Wait for toto
  int titi;
  while (1) {
    int value = CWP_Param_is("code2", "toto", CWP_INT);
    if (value == 1) {
      CWP_Param_get("code2", "toto", CWP_INT, &titi);
      break;
    }
  }

  // Parameter reduce
  int param_value = 0;
  CWP_Param_reduce(CWP_OP_SUM,
                   "toto",
                   CWP_INT,
                   &param_value,
                   2,
                   code_list);

  if (verbose) {
    log_trace("param_value code1 : %d\n", param_value);
  }

  // Finalize CWIPI
  CWP_Finalize();

  // free
  free(code_name);
  free(intra_comm);
  free(code_list);
  free(loc_code_list);
  for (int i = 0; i < n_param; i++) {
    free(param_list[i]);
  }
  free(param_list);
  free(get_tata);

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}

