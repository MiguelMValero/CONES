/*
  This file is part of the CWIPI library.

  Copyright (C) 2022-2023  ONERA

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

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "client_server/server.h"
#include <pdm_error.h>
#include <pdm_io.h>
#include <pdm_mpi.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include "pdm_logging.h"
#include "pdm_printf.h"

#include "cwp_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Util functions
 *============================================================================*/

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -p     Begin and end port number for the server sockets port range.\n\n"
     "  -cn    Code identifier. \n\n"
     "  -h     This message.\n\n");

  exit(exit_code);
}

static void
_read_args
(
 int            argc,
 char         **argv,
 char         **config,     // filename for server ip adresses + ports
 int           *port_begin, // begin of port range
 int           *port_end,   // end of port range
 char         **code_name   // code name
)
{
  int i = 1;

  /* Parse and check command line */

  while (i < argc) {

    if (strcmp(argv[i], "-h") == 0)
      _usage(EXIT_SUCCESS);

    else if (strcmp(argv[i], "-c") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *config = argv[i];
      }
    }

    else if (strcmp(argv[i], "-p") == 0) {

      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *port_begin = atoi(argv[i]);
      }

      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *port_end = atoi(argv[i]);
      }
    }

    else if (strcmp(argv[i], "-cn") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *code_name = argv[i];
      }
    }

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}

/*=============================================================================
 * Main
 *============================================================================*/

int
main
(
 int argc,
 char *argv[]
)
{
  // default
  char *config     = NULL;
  int   port_begin = 49100;
  int   port_end   = 49150;
  char* code_name  = NULL;

  _read_args(argc,
             argv,
             &config,
             &port_begin,
             &port_end,
             &code_name);

  if (code_name == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Server must be launched with a non NULL code identifier.\n");
  }

  if (config == NULL) {
    config = (char *) "cwp_config_srv.txt";
  }

  // mpi
  int i_rank;
  int n_rank;

  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &i_rank);
  MPI_Comm_size(comm, &n_rank);

  // create intracomm
  size_t  s_code_name = strlen(code_name);
  int *s_recv = malloc(sizeof(int) * n_rank);

  MPI_Allgather(&s_code_name,
                1,
                MPI_INT,
                s_recv,
                1,
                MPI_INT,
                comm);

  int *idx_recv = malloc(sizeof(int) * n_rank);
  int s_total   = s_recv[0];
  idx_recv[0]   = 0;

  for (int i = 1; i < n_rank; i++) {
    s_total    += s_recv[i-1];
    idx_recv[i] = idx_recv[i-1] + s_recv[i-1];
  }

  char *code_names = malloc(sizeof(char) * s_total);
  MPI_Allgatherv(code_name,
                 s_code_name,
                 MPI_CHAR,
                 code_names,
   (const int *) s_recv,
   (const int *) idx_recv,
                 MPI_CHAR,
                 comm);

  // post
  char **post_code_names = malloc(sizeof(char *) * n_rank);
  for (int i = 0; i < n_rank; i++) {
    post_code_names[i] = malloc(sizeof(char) * (s_recv[i] + 1));
  }

  for (int i = 0; i < n_rank; i++) {
    int idx = idx_recv[i];
    post_code_names[i][s_recv[i]] = '\0';
    memcpy(post_code_names[i], &code_names[idx], s_recv[i]);
  }

  // set ids
  int n_code_name    = 0;
  int *code_ids      = malloc(sizeof(int) * n_rank);
  int *code_name_idx = malloc(sizeof(int) * n_rank);

  for (int i = 0; i < n_rank; i++) {
    code_ids[i] = -1;

    for (int j = 0; j < n_code_name; j++) {
      int jj = code_name_idx[j];

      if (strcmp(post_code_names[jj], post_code_names[i]) == 0) {
        code_ids[i] = code_ids[jj];
      }
    }
    if (code_ids[i] == -1) {
      code_name_idx[n_code_name++] = i;
      code_ids[i] = i;
    }
  }
  for (int i = 0; i < n_rank; i++) {
    free(post_code_names[i]);
  }
  free(post_code_names);

  MPI_Barrier(comm);

  MPI_Comm intra_comm;
  MPI_Comm_split(comm, code_ids[i_rank], i_rank, &intra_comm);

  // free
  free(code_name_idx);
  free(code_ids);
  free(code_names);
  free(idx_recv);
  free(s_recv);

  int i_intra_rank;
  int n_intra_rank;
  MPI_Comm_rank(intra_comm, &i_intra_rank);
  MPI_Comm_size(intra_comm, &n_intra_rank);

  // port choice
  MPI_Comm comm_node;
  int i_rank_node;

  // shared comm split
  MPI_Comm_split_type(intra_comm, MPI_COMM_TYPE_SHARED, i_intra_rank, MPI_INFO_NULL, &comm_node);
  MPI_Comm_rank(comm_node, &i_rank_node);

  uint16_t server_port = (uint16_t) (port_begin + i_rank_node);

  // create server
  p_server svr = malloc(sizeof(t_server));
  if (CWP_server_create(comm, server_port, 0, svr) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Server creation failed\n");
    return -1;
  }

  // write config file
  // --> retreive host_name size
  int  irank_host_name_size     = strlen(svr->host_name);
  int *all_jrank_host_name_size = malloc(sizeof(int) * n_intra_rank);

  MPI_Allgather(&irank_host_name_size,
                1,
                MPI_INT,
                all_jrank_host_name_size,
                1,
                MPI_INT,
                intra_comm);

  int max_host_name_size = -1;
  for (int i = 0; i < n_intra_rank; i++) {
    if (all_jrank_host_name_size[i] > max_host_name_size) {
      max_host_name_size = all_jrank_host_name_size[i];
    }
  }

  free(all_jrank_host_name_size);

  // --> create format and data string
  char format[99];
  sprintf(format,"%s%d.%ds/%s9.9d\n", "%", max_host_name_size, max_host_name_size, "%");

  char data[max_host_name_size + 12];
  sprintf(data, format, svr->host_name, svr->port);

  // --> open
  PDM_io_file_t *write = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE, // PDM_IO_KIND_MPIIO_EO,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              PDM_MPI_mpi_2_pdm_mpi_comm(&intra_comm),
              -1.,
              &write,
              &ierr);

  // --> global write of header
  char  buf[99];
  sprintf(buf, "FORMAT hostname/port\nN %10.10ld\nSIZE %5.5ld\n", (long) n_intra_rank, strlen(data)); // header_size = 45 char

  size_t s_buf =  strlen(buf);
  PDM_io_global_write(write,
        (PDM_l_num_t) sizeof(char),
        (PDM_l_num_t) s_buf,
                      buf);

  // --> block write of data
  int one = 1;
  PDM_g_num_t i_rank_gnum = (PDM_g_num_t) (i_intra_rank+1);

  PDM_io_par_interlaced_write(write,
                              PDM_STRIDE_CST_INTERLACED,
              (PDM_l_num_t *) &one,
                (PDM_l_num_t) strlen(data),
                (PDM_l_num_t) one,
                              &i_rank_gnum,
               (const void *) data);

  // --> close
  PDM_io_close(write);
  PDM_io_free(write);

  // verbose
  MPI_Barrier(comm);

  if (i_rank == 0) {
    printf("----------------------------------------------------------------------------\n");
    printf("All servers listening and cwipi config file created. You may connect clients\n");
    printf("----------------------------------------------------------------------------\n");
  }

  // accept
  CWP_server_run(svr);

  // shutdown server
  CWP_server_kill(svr);

  // free
  free(svr);

  // mpi finalize
  MPI_Comm_free(&intra_comm);
  MPI_Comm_free(&comm_node);
  MPI_Finalize();

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
