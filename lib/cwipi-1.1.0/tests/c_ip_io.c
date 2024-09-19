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
#include <string.h>
#include <stdint.h>
#include <math.h>

#include <pdm_error.h>
#include <pdm_io.h>
#include <pdm_mpi.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include "pdm_logging.h"
#include "pdm_printf.h"

#include "cwp_priv.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*
 * Function definitions
 */

static void
_usage(int exit_code)
{
  PDM_printf
    ("\n"
     "  Usage: \n\n"
     "  -c     Filename of the server configuration file.\n\n"
     "  -p     Begin and end port number for the server sockets port range.\n\n"
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
 int           *port_end    // end of port range
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

    else
      _usage(EXIT_FAILURE);
    i++;

  }

}

/*
 * Main
 */

int main(int argc, char *argv[])
{

  // default
  char *config     = NULL;
  int   port_begin = 1024;
  int   port_end   = 49151;

  _read_args(argc,
             argv,
             &config,
             &port_begin,
             &port_end);

  if (config == NULL) {
    config = (char *) "cwp_config_srv.txt";
  }

  // mpi
  int i_rank;
  int n_rank;

  PDM_MPI_Init(&argc, &argv);
  PDM_MPI_Comm comm = PDM_MPI_COMM_WORLD;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // port choice (test numa)
  PDM_MPI_Comm comm_node;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_node); // PDM_MPI_SPLIT_SHARED

  int i_rank_node;
  PDM_MPI_Comm_rank(comm_node, &i_rank_node);

  int port = port_begin + i_rank_node;

  // retreive host_name
  char *host_name = malloc(99);
  gethostname(host_name, 99);

  // determine max host_name size
  int  irank_host_name_size     = strlen(host_name);
  int *all_jrank_host_name_size = malloc(sizeof(int) * n_rank);

  PDM_MPI_Allgather(&irank_host_name_size,
                    1,
                    PDM_MPI_INT,
                    all_jrank_host_name_size,
                    1,
                    PDM_MPI_INT,
                    comm);

  int max_host_name_size = -1;
  for (int i = 0; i < n_rank; i++) {
    if (all_jrank_host_name_size[i] > max_host_name_size) {
      max_host_name_size = all_jrank_host_name_size[i];
    }
  }

  // create string: host_name/port\n "%?.?s/9.9d\n"
  char format[99];
  sprintf(format,"%s%d.%ds/%s9.9d\n", "%", max_host_name_size, max_host_name_size, "%");

  log_trace("%s", format);

  char data[max_host_name_size + 12];
  sprintf(data, format, host_name, port);

  // write with pdm_io
  // --> open

  PDM_io_file_t *unite = NULL;
  PDM_l_num_t    ierr;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE, // PDM_IO_KIND_MPIIO_EO,
              PDM_IO_MOD_WRITE,
              PDM_IO_NATIVE,
              comm,
              -1.,
              &unite,
              &ierr);

  // --> global write: header and offset

  char  buf[99];
  sprintf(buf, "FORMAT hostname/port\nSIZE %5.5ld\n", strlen(data)); // header_size = 30 char

  size_t s_buf =  strlen(buf);
  PDM_io_global_write(unite,
        (PDM_l_num_t) sizeof(char),
        (PDM_l_num_t) s_buf,
                      buf);

  // --> par_block_write

  log_trace("%s", data);

  int one = 1;
  PDM_g_num_t debut_bloc = 0; // i_rank * strlen(data) + strlen(buf)
  PDM_UNUSED(debut_bloc);
  PDM_g_num_t i_rank_gnum = (PDM_g_num_t) (i_rank+1);

  PDM_io_par_interlaced_write(unite,
                         PDM_STRIDE_CST_INTERLACED,
         (PDM_l_num_t *) &one, // n_composantes
           (PDM_l_num_t) strlen(data),
           (PDM_l_num_t) one,  // n_donnees
                        &i_rank_gnum, // debut_bloc,
          (const void *) data);

  // --> close

  PDM_io_close(unite);
  PDM_io_free(unite);

  // read with pdm_io
  // --> open

  PDM_io_file_t *read = NULL;
  PDM_l_num_t    ierr_read;

  PDM_io_open(config,
              PDM_IO_FMT_BIN,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              PDM_IO_KIND_MPI_SIMPLE,
              PDM_IO_MOD_READ,
              PDM_IO_NATIVE,
              comm,
              -1.,
              &read,
              &ierr_read);

  // --> global read of header

  char *buffer = malloc(99);
  for (int i = 0; i < 99; i++) {
    buffer[i] = '\0';
  }

  PDM_io_global_read(read,
                     32 * sizeof(char),
                     1,
                     buffer);

  log_trace("%s", buffer);

  char div[] = "\n";
  char *readbuff = strtok(buffer, div);

  log_trace("%s\n", readbuff);
  readbuff = strtok(NULL, div);

  log_trace("%s\n", readbuff);

  char div1[] = " ";
  char *l2 = strtok(readbuff, div1);
  log_trace("%s\n", l2);
  l2 = strtok(NULL, div1);

  log_trace("%s\n", l2);

  int size = atoi(l2);

  log_trace("size: %ld (real), %d (read)", strlen(data), size);

  // --> read data (hostname/port);
  PDM_g_num_t debut_bloc_read = i_rank * size + 32;
  PDM_UNUSED(debut_bloc_read);

  char *read_data = malloc(size+1);

  for (int i = 0; i < size+1; i++) {
    read_data[i] = '\0';
  }

  PDM_io_par_interlaced_read(read,
                        PDM_STRIDE_CST_INTERLACED,
         (PDM_l_num_t *) &one, // n_composantes
           (PDM_l_num_t) size,
           (PDM_l_num_t) one,  // n_donnees
                         &i_rank_gnum,
                         read_data);

  // retreive hostname and port seperatly
  char divider[] = "/";
  char *read_str = strtok(read_data, divider);
  log_trace("%s\n", read_str);
  read_str = strtok(NULL, divider);
  int read_port = atoi(read_str);
  log_trace("read_port = %d\n", read_port);

  // --> close

  PDM_io_close(read);
  PDM_io_free(read);

  // free
  free(all_jrank_host_name_size);

  PDM_MPI_Finalize();

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
