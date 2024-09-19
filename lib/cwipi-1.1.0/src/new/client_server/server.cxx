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

/*
  This file is inspired from OpenPALM.
  OpenPALM is a free software under the GNU Lesser General Public License.
  See: https://www.cerfacs.fr/globc/PALM_WEB/
*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <arpa/inet.h>
#include <errno.h>
#include <iostream>
#include <map>
#include <cassert>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "server.h"
#include "transfer.h"
#include "client.h"
#include "message.h"
#include "struct.hxx"

#include "cwp_priv.h"

#include <pdm_error.h>
#include <pdm_logging.h>
#include <pdm_mpi.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* debug */

static int svr_debug = 0;

/* file struct definition */

static t_server_mpi svr_mpi;
static t_cwp        svr_cwp;

/*=============================================================================
 * Private function interfaces
 *============================================================================*/

// --> wrapper

static void read_name(char **name,  p_server svr) {
  int name_size;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
  *name = (char *) realloc(*name, name_size);
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) *name, name_size);
}

static void write_name(char * name, p_server svr) {
  int name_size = strlen(name) + 1; // +1 for "\0"
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, (void*) name, name_size);
}

// --> n_vtx

static int n_vtx_get(CWP_Block_t block_type) {
  int n_vtx = 0;

  switch(block_type) {

  case CWP_BLOCK_NODE: {
    n_vtx = 1;
    } break;

  case CWP_BLOCK_EDGE2: {
    n_vtx = 2;
    } break;

  case CWP_BLOCK_FACE_TRIA3: {
    n_vtx = 3;
    } break;

  case CWP_BLOCK_FACE_QUAD4: {
    n_vtx = 4;
    } break;

  case CWP_BLOCK_FACE_POLY: {
    // TO DO: what value ?
    PDM_error(__FILE__, __LINE__, 0, "Number of nodes unknown for CWP_BLOCK_FACE_POLY\n");
    return -1;
    } break;

  case CWP_BLOCK_CELL_TETRA4: {
    n_vtx = 4;
    } break;

  case CWP_BLOCK_CELL_HEXA8: {
    n_vtx = 8;
    } break;

  case CWP_BLOCK_CELL_PRISM6: {
    n_vtx = 6;
    } break;

  case CWP_BLOCK_CELL_PYRAM5: {
    n_vtx = 5;
    } break;

  case CWP_BLOCK_CELL_POLY: {
    // TO DO: what value ?
    PDM_error(__FILE__, __LINE__, 0, "Number of nodes unknown for CWP_BLOCK_CELL_POLY\n");
    return -1;
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown block type %d\n", block_type);
    return -1;
  }

  return n_vtx;
}

// --> ho_n_vtx

static int ho_n_vtx_get(CWP_Block_t block_type, int order) {
switch (block_type) {
  case CWP_BLOCK_NODE :
    return 1;
    break;
  case CWP_BLOCK_EDGE2 :
    return 2;
    break;
  case CWP_BLOCK_FACE_TRIA3 :
    return 3;
    break;
  case CWP_BLOCK_FACE_QUAD4 :
    return 4;
    break;
  case CWP_BLOCK_CELL_TETRA4 :
    return 4;
    break;
  case CWP_BLOCK_CELL_PYRAM5 :
    return 5;
    break;
  case CWP_BLOCK_CELL_PRISM6 :
    return 6;
    break;
  case CWP_BLOCK_CELL_HEXA8 :
    return 8;
    break;

  case CWP_BLOCK_EDGEHO :
    return order + 1;
    break;
  case CWP_BLOCK_FACE_TRIAHO :
    return (order + 1) * (order + 2) / 2;
    break;
  case CWP_BLOCK_FACE_QUADHO :
    return (order + 1) * (order + 1);
    break;
  case CWP_BLOCK_CELL_TETRAHO :
    return (order + 1) * (order + 2) * (order + 3) / 6;
    break;
  case CWP_BLOCK_CELL_PYRAMHO :
    return (order + 1) * (order + 2) * (2*order + 3) / 6;
    break;
  case CWP_BLOCK_CELL_PRISMHO :
    return (order + 1) * (order + 1) * (order + 2) / 2;
    break;
  case CWP_BLOCK_CELL_HEXAHO :
    return (order + 1) * (order + 1) * (order + 1);
    break;
  default :
    PDM_error (__FILE__, __LINE__, 0, "Unknown order for Poly2D and Poly3D (block_type %d)\n", block_type);
  }
  return -1;
}

/*=============================================================================
 * Server CWIPI function interfaces
 *============================================================================*/

void
CWP_server_Init
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.global_comm);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INIT);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  int            n_code;
  char         **code_names     = NULL;
  CWP_Status_t   is_active_rank = CWP_STATUS_OFF;

  // receive data
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &n_code, sizeof(int));

  if (n_code > 1) {
    PDM_error(__FILE__, __LINE__, 0, "CWIPI client-server not implemented yet for n_code > 1\n");
  }

  code_names = (char **) malloc(sizeof(char *) * n_code);
  for (int i = 0; i < n_code; i++) {
    code_names[i] = (char *) malloc(sizeof(char));
    read_name(&code_names[i], svr);
  }
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &is_active_rank, n_code * sizeof(CWP_Status_t));

  // send status msg
  MPI_Barrier(svr_mpi.global_comm);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INIT);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch CWP_Init
  svr_mpi.intra_comms = (MPI_Comm *) malloc(sizeof(MPI_Comm) * n_code);
  CWP_Init(svr_mpi.global_comm,
           n_code,
           (const char **) code_names,
           is_active_rank,
           svr_mpi.intra_comms);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INIT);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  for (int i = 0; i < n_code; i++) {
    free(code_names[i]);
  }
  free(code_names);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Finalize
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FINALIZE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  if (svr_mpi.intra_comms != NULL) free(svr_mpi.intra_comms);

  if (svr_cwp.code_names != NULL) {
    for (int i = 0; i < svr_cwp.n_code_names; i++) {
      if (svr_cwp.code_names[i] != NULL) free((void *) svr_cwp.code_names[i]);
    }
    free(svr_cwp.code_names);
  }

  if (svr_cwp.loc_code_names != NULL) {
    for (int i = 0; i < svr_cwp.n_loc_code_names; i++) {
      if (svr_cwp.loc_code_names[i] != NULL) free((void *) svr_cwp.loc_code_names[i]);
    }
    free(svr_cwp.loc_code_names);
  }

  if (!svr_cwp.char_param_value.empty()) {
    for (const auto& x : svr_cwp.char_param_value) {
      if (x.second != NULL) {
        free((void *) x.second);
      }
    }
    svr_cwp.char_param_value.clear();
  }

  // launch
  CWP_Finalize();

  // send status msg
  MPI_Barrier(svr_mpi.global_comm);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FINALIZE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_lock
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_LOCK);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // receive code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_LOCK);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Param_lock((const char *) code_name);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_LOCK);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_unlock
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_UNLOCK);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // receive code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_UNLOCK);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Param_unlock((const char *) code_name);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_UNLOCK);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_add
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_ADD);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read data type
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  switch (data_type) {

  case CWP_DOUBLE: {
    double double_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &double_initial_value, sizeof(double));

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_ADD);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &double_initial_value);
    } break;

  case CWP_INT: {
    double int_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &int_initial_value, sizeof(int));

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_ADD);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &int_initial_value);
    } break;

  case CWP_CHAR: {
    int name_size;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
    char *char_initial_value = (char *) malloc(name_size);
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) char_initial_value, name_size);
    std::string s(param_name);
    svr_cwp.char_param_value.insert(std::make_pair(s, char_initial_value));

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_ADD);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Param_add(local_code_name,
                  param_name,
                  data_type,
                  &svr_cwp.char_param_value[s]);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_ADD);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_get
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read value data_type
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send value
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  switch (data_type) {

  case CWP_DOUBLE: {
    double double_value;

    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &double_value);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_GET);
      message.flag = CWP_SVR_LCH_END;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &double_value, sizeof(double));
    } break;

  case CWP_INT: {
    int int_value;
    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &int_value);

      // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_GET);
      message.flag = CWP_SVR_LCH_END;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &int_value, sizeof(int));
    } break;

  case CWP_CHAR: {
    char *char_value = NULL;
    CWP_Param_get(local_code_name,
                  param_name,
                  data_type,
                  &char_value);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_GET);
      message.flag = CWP_SVR_LCH_END;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    write_name(char_value, svr);
    free(char_value);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  // free
  free(local_code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_set
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  switch (data_type) {

  case CWP_DOUBLE: {
    double double_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &double_initial_value, sizeof(double));

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &double_initial_value);
    } break;

  case CWP_INT: {
    double int_initial_value;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &int_initial_value, sizeof(int));

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &int_initial_value);
    } break;

  case CWP_CHAR: {
    int name_size;
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) &name_size, sizeof(int));
    std::string s(param_name);
    free(svr_cwp.char_param_value[s]);
    svr_cwp.char_param_value[s] = (char *) malloc(name_size);
    CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, (void*) svr_cwp.char_param_value[s], name_size);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Param_set(local_code_name,
                  param_name,
                  data_type,
                  &svr_cwp.char_param_value[s]);
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown CWP_Type_t %i\n", data_type);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_del
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_DEL);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_DEL);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Param_del(local_code_name,
                param_name,
                data_type);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_DEL);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  if (data_type == CWP_CHAR) {
    std::string s(param_name);
    if (svr_cwp.char_param_value[s] != NULL) free(svr_cwp.char_param_value[s]);
    svr_cwp.char_param_value.erase(s);
  }

  free(local_code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_n_get
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_N_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_N_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int n_param = CWP_Param_n_get(code_name,
                                data_type);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_N_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send n_param
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_param, sizeof(int));

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_list_get
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_LIST_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // read data type
  CWP_Type_t data_type = (CWP_Type_t ) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_LIST_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  svr_cwp.n_param_names = -1;
  svr_cwp.param_names = NULL;
  CWP_Param_list_get(code_name,
                     data_type,
                     &svr_cwp.n_param_names,
                     &svr_cwp.param_names);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_LIST_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send nParam
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &svr_cwp.n_param_names, sizeof(int));

  // send paramNames
  for (int i = 0; i < svr_cwp.n_param_names; i++) {
    write_name(svr_cwp.param_names[i], svr);
  }

  // free
  free(code_name);

  if (svr_cwp.param_names != NULL) {
    for (int i = 0; i < svr_cwp.n_param_names; i++) {
      if (svr_cwp.param_names[i] != NULL) free(svr_cwp.param_names[i]);
    }
    free(svr_cwp.param_names);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_is
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_IS);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read initial value
  CWP_Type_t data_type = (CWP_Type_t) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_IS);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int is_param = CWP_Param_is(code_name,
                              param_name,
                              data_type);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_IS);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send is_param
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &is_param, sizeof(int));

  // free
  free(code_name);
  free(param_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Param_reduce
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_REDUCE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read operation
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  CWP_Op_t op = (CWP_Op_t) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &op, sizeof(CWP_Op_t));

  // read param name
  char *param_name = (char *) malloc(sizeof(char));
  read_name(&param_name, svr);

  // read data type
  CWP_Type_t data_type = (CWP_Type_t) -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &data_type, sizeof(CWP_Type_t));

  // read number of codes
  int nCode = -1;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &nCode, sizeof(int));

  // read code names
  char **code_names = (char **) malloc(sizeof(char *) * nCode);
  for (int i = 0; i < nCode; i++) {
    code_names[i] = (char *) malloc(sizeof(char));
    read_name(&code_names[i], svr);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PARAM_REDUCE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send res
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  switch (data_type) {

  case CWP_DOUBLE: {
    double res = -1.;
    CWP_Param_reduce(op,
                     param_name,
                     data_type,
                     (void *) &res,
                     nCode,
     (const char **) code_names);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_REDUCE);
      message.flag = CWP_SVR_LCH_END;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &res, sizeof(double));
    } break;

  case CWP_INT: {
    int res = -1;
    CWP_Param_reduce(op,
                     param_name,
                     data_type,
                     (void *) &res,
                     nCode,
     (const char **) code_names);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_PARAM_REDUCE);
      message.flag = CWP_SVR_LCH_END;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &res, sizeof(int));
    } break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Unknown CWP_Type_t %i or impossible on CWP_CHAR\n", data_type);
  }

  // free
  free(param_name);
  for (int i = 0; i < nCode; i++) {
    free(code_names[i]);
  }
  free(code_names);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Cpl_create
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_CREATE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read coupled code name
  char *coupled_code_name = (char *) malloc(sizeof(char));
  read_name(&coupled_code_name, svr);

  // read entities dimension
  CWP_Interface_t entities_dim;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &entities_dim, sizeof(CWP_Interface_t));

  // read communication type
  CWP_Comm_t comm_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &comm_type, sizeof(CWP_Comm_t));

  // read spatial interpolation
  CWP_Spatial_interp_t spatial_interp;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &spatial_interp, sizeof(CWP_Spatial_interp_t));

  // read number of partitions
  int n_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_part, sizeof(int));

  // read displacement type
  CWP_Dynamic_mesh_t displacement;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &displacement, sizeof(CWP_Dynamic_mesh_t));

  // read time exchange type
  CWP_Time_exch_t recv_freq_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &recv_freq_type, sizeof(CWP_Time_exch_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_CREATE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Cpl_create(local_code_name,
                 cpl_id,
                 coupled_code_name,
                 entities_dim,
                 comm_type,
                 spatial_interp,
                 n_part,
                 displacement,
                 recv_freq_type);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_CREATE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // create occurence in map
  std::string s(cpl_id);
  t_coupling coupling = t_coupling();
  svr_cwp.coupling.insert(std::make_pair(s, coupling));
  svr_cwp.coupling[s].n_part = n_part;
  svr_cwp.coupling[s].mesh_dynamic = displacement;

  // free
  free(local_code_name);
  free(cpl_id);
  free(coupled_code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Cpl_barrier
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_BARRIER);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_BARRIER);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Cpl_barrier(local_code_name,
                  cpl_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_BARRIER);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
}

void
CWP_server_Cpl_del
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_DEL);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_DEL);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Cpl_del(local_code_name,
              cpl_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_DEL);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // remove from map
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].ijk_grid != NULL) free(svr_cwp.coupling[s].ijk_grid);

  if (svr_cwp.coupling[s].mesh_dynamic == CWP_DYNAMIC_MESH_DEFORMABLE) {
    if (svr_cwp.coupling[s].n_user_vtx != NULL) free(svr_cwp.coupling[s].n_user_vtx);
    if (svr_cwp.coupling[s].n_vtx != NULL) free(svr_cwp.coupling[s].n_vtx);
  }

  // n_part
  if (svr_cwp.coupling[s].vtx_coord != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].vtx_coord[i_part] != NULL) free(svr_cwp.coupling[s].vtx_coord[i_part]);
    }
    free(svr_cwp.coupling[s].vtx_coord);
  }

  if (svr_cwp.coupling[s].vtx_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].vtx_global_num[i_part] != NULL) free(svr_cwp.coupling[s].vtx_global_num[i_part]);
    }
    free(svr_cwp.coupling[s].vtx_global_num);
  }

  if (svr_cwp.coupling[s].usr_coord != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].usr_coord[i_part] != NULL) free(svr_cwp.coupling[s].usr_coord[i_part]);
    }
    free(svr_cwp.coupling[s].usr_coord);
  }

  if (svr_cwp.coupling[s].usr_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].usr_global_num[i_part] != NULL) free(svr_cwp.coupling[s].usr_global_num[i_part]);
    }
    free(svr_cwp.coupling[s].usr_global_num);
  }

  if (svr_cwp.coupling[s].face_edge_idx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_edge_idx[i_part] != NULL) free(svr_cwp.coupling[s].face_edge_idx[i_part]);
    }
    free(svr_cwp.coupling[s].face_edge_idx);
  }

  if (svr_cwp.coupling[s].face_edge != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_edge[i_part] != NULL) free(svr_cwp.coupling[s].face_edge[i_part]);
    }
    free(svr_cwp.coupling[s].face_edge);
  }

  if (svr_cwp.coupling[s].edge_vtx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].edge_vtx[i_part] != NULL) free(svr_cwp.coupling[s].edge_vtx[i_part]);
    }
    free(svr_cwp.coupling[s].edge_vtx);
  }

  if (svr_cwp.coupling[s].face_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_global_num[i_part] != NULL) free(svr_cwp.coupling[s].face_global_num[i_part]);
    }
    free(svr_cwp.coupling[s].face_global_num);
  }

  if (svr_cwp.coupling[s].cell_face_idx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].cell_face_idx[i_part] != NULL) free(svr_cwp.coupling[s].cell_face_idx[i_part]);
    }
    free(svr_cwp.coupling[s].cell_face_idx);
  }

  if (svr_cwp.coupling[s].cell_face != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].cell_face[i_part] != NULL) free(svr_cwp.coupling[s].cell_face[i_part]);
    }
    free(svr_cwp.coupling[s].cell_face);
  }

  if (svr_cwp.coupling[s].face_vtx_idx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_vtx_idx[i_part] != NULL) free(svr_cwp.coupling[s].face_vtx_idx[i_part]);
    }
    free(svr_cwp.coupling[s].face_vtx_idx);
  }

  if (svr_cwp.coupling[s].face_vtx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_vtx[i_part] != NULL) free(svr_cwp.coupling[s].face_vtx[i_part]);
    }
    free(svr_cwp.coupling[s].face_vtx);
  }

  if (svr_cwp.coupling[s].cell_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].cell_global_num[i_part] != NULL) free(svr_cwp.coupling[s].cell_global_num[i_part]);
    }
    free(svr_cwp.coupling[s].cell_global_num);
  }

  // Block
  if (!svr_cwp.coupling[s].block.empty()) {
    std::map<int, t_block>::iterator it_b = svr_cwp.coupling[s].block.begin();
    while (it_b != svr_cwp.coupling[s].block.end()) {

      // n_part
      if ((it_b->second).connec_faces_idx != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_faces_idx[i_part] != NULL) free((it_b->second).connec_faces_idx[i_part]);
        }
        free((it_b->second).connec_faces_idx);
      }

      if ((it_b->second).connec_faces != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_faces[i_part] != NULL) free((it_b->second).connec_faces[i_part]);
        }
        free((it_b->second).connec_faces);
      }

      if ((it_b->second).connec_cells_idx != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_cells_idx[i_part] != NULL) free((it_b->second).connec_cells_idx[i_part]);
        }
        free((it_b->second).connec_cells_idx);
      }

      if ((it_b->second).connec_cells != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_cells[i_part] != NULL) free((it_b->second).connec_cells[i_part]);
        }
        free((it_b->second).connec_cells);
      }

      if ((it_b->second).cell_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).cell_global_num[i_part] != NULL) free((it_b->second).cell_global_num[i_part]);
        }
        free((it_b->second).cell_global_num);
      }

      if ((it_b->second).connec_idx != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_idx[i_part] != NULL) free((it_b->second).connec_idx[i_part]);
        }
        free((it_b->second).connec_idx);
      }

      if ((it_b->second).connec != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec[i_part] != NULL) free((it_b->second).connec[i_part]);
        }
        free((it_b->second).connec);
      }

      if ((it_b->second).elt_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).elt_global_num[i_part] != NULL) free((it_b->second).elt_global_num[i_part]);
        }
        free((it_b->second).elt_global_num);
      }

      if ((it_b->second).std_connec != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).std_connec[i_part] != NULL) free((it_b->second).std_connec[i_part]);
        }
        free((it_b->second).std_connec);
      }

      if ((it_b->second).std_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).std_global_num[i_part] != NULL) free((it_b->second).std_global_num[i_part]);
        }
        free((it_b->second).std_global_num);
      }

      if ((it_b->second).ho_std_connec != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).ho_std_connec[i_part] != NULL) free((it_b->second).ho_std_connec[i_part]);
        }
        free((it_b->second).ho_std_connec);
      }

      if ((it_b->second).ho_std_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).ho_std_global_num[i_part] != NULL) free((it_b->second).ho_std_global_num[i_part]);
        }
        free((it_b->second).ho_std_global_num);
      }

      it_b = svr_cwp.coupling[s].block.erase(it_b);
    }
  }

  // Field
  if (!svr_cwp.coupling[s].field.empty()) {
    std::map<std::string, t_field>::iterator it_f = svr_cwp.coupling[s].field.begin();
    while (it_f != svr_cwp.coupling[s].field.end()) {
      if ((it_f->second).data != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_f->second).data[i_part] != NULL) free((it_f->second).data[i_part]);
          (it_f->second).data[i_part] = NULL;
        }
        free((it_f->second).data);
      }

      if ((it_f->second).size != NULL) free((it_f->second).size);
      it_f = svr_cwp.coupling[s].field.erase(it_f);
    }
  }

  // Global Data
  if (!svr_cwp.coupling[s].global_data.empty()) {
    std::map<std::string, t_global_data>::iterator it_gd = svr_cwp.coupling[s].global_data.begin();
    while (it_gd != svr_cwp.coupling[s].global_data.end()) {
      if ((it_gd->second).send_data != NULL) free((it_gd->second).send_data);
      it_gd = svr_cwp.coupling[s].global_data.erase(it_gd);
    }
  }

  // Part Data
  if (!svr_cwp.coupling[s].part_data.empty()) {
    std::map<std::string, t_part_data>::iterator it_pd = svr_cwp.coupling[s].part_data.begin();
    while (it_pd != svr_cwp.coupling[s].part_data.end()) {

      if ((it_pd->second).n_elt != NULL) {
        free((it_pd->second).n_elt);
      }

      if ((it_pd->second).gnum_elt != NULL) {
        for (int i_part = 0; i_part < (it_pd->second).n_part; i_part++) {
          if ((it_pd->second).gnum_elt[i_part] != NULL) {
            free((it_pd->second).gnum_elt[i_part]);
          }
        }
        free((it_pd->second).gnum_elt);
      }

      it_pd = svr_cwp.coupling[s].part_data.erase(it_pd);
    }
  }

  // Properties
  if (!svr_cwp.coupling[s].property.empty()) {
    std::map<std::string, t_property>::iterator it_p = svr_cwp.coupling[s].property.begin();
    while (it_p != svr_cwp.coupling[s].property.end()) {
      if ((it_p->second).property_name  != NULL) free((it_p->second).property_name);
      if ((it_p->second).property_value != NULL) free((it_p->second).property_value);
      it_p = svr_cwp.coupling[s].property.erase(it_p);
    }
  }

  svr_cwp.coupling.erase(s);

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Properties_dump
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PROPERTIES_DUMP);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  char *properties = NULL;
  int size = CWP_Properties_str_dump(&properties); // .length + 1

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PROPERTIES_DUMP);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, (void*) &size, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, (void*) properties, size);

  // free
  free(properties);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Visu_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_VISU_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read frequency
  int freq;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &freq, sizeof(int));

  // read format
  CWP_Visu_format_t format;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &format, sizeof(CWP_Visu_format_t));

  // read format option
  char *format_option = (char *) malloc(sizeof(char));
  read_name(&format_option, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_VISU_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Visu_set(local_code_name,
               cpl_id,
               freq,
               format,
               format_option);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_VISU_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(format_option);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_State_update
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STATE_UPDATE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read state
  CWP_State_t state;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &state, sizeof(CWP_State_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STATE_UPDATE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_State_update(local_code_name,
                   state);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STATE_UPDATE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Time_step_beg
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_TIME_STEP_BEG);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read current time
  double current_time = - 1.0;
  CWP_transfer_readdata(svr->connected_socket, svr->max_msg_size, &current_time, sizeof(double));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_TIME_STEP_BEG);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Time_step_beg(local_code_name,
                    current_time);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_TIME_STEP_BEG);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Time_step_end
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_TIME_STEP_END);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_TIME_STEP_END);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Time_step_end(local_code_name);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_TIME_STEP_END);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Output_file_set
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi.intra_comms[0]);

  PDM_UNUSED(svr);
  printf("CWP: CWP_Output_file_set not implemented on server side\n");

  // // read code name
  // svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  // char *output_filename = (char *) malloc(sizeof(char));
  // read_name(&output_filename, svr);

  // // create FILE *
  // svr_cwp.output_file = NULL;

  // svr_cwp.output_file = fopen(output_filename, "a+");

  // // launch
  // CWP_Output_file_set(svr_cwp.output_file);

  // // free
  // free(output_filename);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_User_structure_set
(
  p_server                 svr
)
{
  PDM_UNUSED(svr);
  printf("CWP: CWP_User_structure_set not implemented in client/server mode\n");
}

void
CWP_server_User_structure_get
(
  p_server                 svr
)
{
  PDM_UNUSED(svr);
  printf("CWP: CWP_User_structure_get not implemented in client/server mode\n");
}

void
CWP_server_State_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STATE_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *code_name = (char *) malloc(sizeof(char));
  read_name(&code_name, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STATE_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int state = CWP_State_get(code_name);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STATE_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send state
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &state, sizeof(int));

  // free
  free(code_name);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Codes_nb_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CODES_NB_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_codes = CWP_Codes_nb_get();

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CODES_NB_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_codes, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Codes_list_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CODES_LIST_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // get number of codes
  int nb_codes = CWP_Codes_nb_get();

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CODES_LIST_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_codes, sizeof(int));

  // launch
  svr_cwp.code_names = (char **) CWP_Codes_list_get();
  for (int i = 0; i < nb_codes; i++) {
    write_name((char *) svr_cwp.code_names[i], svr);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Loc_codes_nb_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_LOC_CODES_NB_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_local_codes = CWP_Loc_codes_nb_get();

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_LOC_CODES_NB_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_local_codes, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Loc_codes_list_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_LOC_CODES_LIST_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // get number of codes
  int nb_local_codes = CWP_Loc_codes_nb_get();

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_LOC_CODES_LIST_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of codes
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_local_codes, sizeof(int));

  // launch
  svr_cwp.loc_code_names = (char **) CWP_Loc_codes_list_get();
  for (int i = 0; i < nb_local_codes; i++) {
    write_name((char *) svr_cwp.loc_code_names[i], svr);
  }

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Computed_tgts_bcast_enable
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_COMPUTED_TGTS_BCAST_ENABLE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_COMPUTED_TGTS_BCAST_ENABLE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Computed_tgts_bcast_enable(local_code_name,
                                 cpl_id,
                                 field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_COMPUTED_TGTS_BCAST_ENABLE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Involved_srcs_bcast_enable
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INVOLVED_SRCS_BCAST_ENABLE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INVOLVED_SRCS_BCAST_ENABLE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Involved_srcs_bcast_enable(local_code_name,
                                 cpl_id,
                                 field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INVOLVED_SRCS_BCAST_ENABLE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_uncomputed_tgts_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_tgts = CWP_N_uncomputed_tgts_get(local_code_name,
                                          cpl_id,
                                          field_id,
                                          i_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Uncomputed_tgts_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_UNCOMPUTED_TGTS_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_UNCOMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_tgts = CWP_N_uncomputed_tgts_get(local_code_name,
                                          cpl_id,
                                          field_id,
                                          i_part);

  const int *tgts = NULL;
  tgts = CWP_Uncomputed_tgts_get(local_code_name,
                                 cpl_id,
                                 field_id,
                                 i_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_UNCOMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) tgts, sizeof(int) * nb_tgts);

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_computed_tgts_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_COMPUTED_TGTS_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_COMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_tgts = CWP_N_computed_tgts_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_COMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Computed_tgts_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_COMPUTED_TGTS_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_COMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_tgts = CWP_N_computed_tgts_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  const int *tgts = NULL;
  tgts = CWP_Computed_tgts_get(local_code_name,
                               cpl_id,
                               field_id,
                               i_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_COMPUTED_TGTS_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_tgts, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) tgts, sizeof(int) * nb_tgts);

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_N_involved_srcs_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_INVOLVED_SRCS_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_INVOLVED_SRCS_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_srcs = CWP_N_involved_srcs_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_N_INVOLVED_SRCS_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of sources
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_srcs, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Involved_srcs_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INVOLVED_SRCS_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INVOLVED_SRCS_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int nb_srcs = CWP_N_involved_srcs_get(local_code_name,
                                        cpl_id,
                                        field_id,
                                        i_part);

  const int *srcs = NULL;
  srcs = CWP_Involved_srcs_get(local_code_name,
                               cpl_id,
                               field_id,
                               i_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_INVOLVED_SRCS_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of targets
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &nb_srcs, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) srcs, sizeof(int) * nb_srcs);

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Spatial_interp_weights_compute
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // receive coordinates if deformable mesh
  std::string s1(cpl_id);
  if (svr_cwp.coupling[s1].mesh_dynamic == CWP_DYNAMIC_MESH_DEFORMABLE) {
    for (int j_part = 0; j_part < svr_cwp.coupling[s1].n_part; j_part++) {
      CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s1].vtx_coord[j_part], sizeof(double) * 3 * svr_cwp.coupling[s1].n_vtx[j_part]);
      if (svr_cwp.coupling[s1].usr_coord != NULL) {
        CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s1].usr_coord[j_part], sizeof(double) * 3 * svr_cwp.coupling[s1].n_user_vtx[j_part]);
      }
    }
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Spatial_interp_weights_compute(local_code_name,
                                     cpl_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Spatial_interp_property_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read property name
  char *property_name = (char *) malloc(sizeof(char));
  read_name(&property_name, svr);

  // read property type
   CWP_Type_t property_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &property_type, sizeof(int));

  // read property value
  char *property_value = (char *) malloc(sizeof(char));
  read_name(&property_value, svr);

  // add name to free or not
  std::string s1(cpl_id);
  std::string s2(property_name);

  if (svr_cwp.coupling[s1].property.find(s2) == svr_cwp.coupling[s1].property.end()) {
    // not found thus create
    t_property property = t_property();
    svr_cwp.coupling[s1].property.insert(std::make_pair(s2, property));

    svr_cwp.coupling[s1].property[s2].property_name  = property_name;
    svr_cwp.coupling[s1].property[s2].property_type  = property_type;
    svr_cwp.coupling[s1].property[s2].property_value = property_value;

  } else {

    free(svr_cwp.coupling[s1].property[s2].property_value);
    svr_cwp.coupling[s1].property[s2].property_value = property_value;

    free(property_name);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Spatial_interp_property_set(local_code_name,
                                  cpl_id,
                                  svr_cwp.coupling[s1].property[s2].property_name,
                                  svr_cwp.coupling[s1].property[s2].property_type,
                                  svr_cwp.coupling[s1].property[s2].property_value);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_User_tgt_pts_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_USER_TGT_PTS_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_pts
  int n_pts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_pts, sizeof(int));

  // read coord
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].usr_coord == NULL) {
    svr_cwp.coupling[s].usr_coord      = (double **) malloc(sizeof(double *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].usr_global_num = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].usr_coord[j_part]      = NULL;
      svr_cwp.coupling[s].usr_global_num[j_part] = NULL;
    }
  }

  if (svr_cwp.coupling[s].mesh_dynamic == CWP_DYNAMIC_MESH_DEFORMABLE) {

    if (svr_cwp.coupling[s].n_user_vtx ==  NULL) {
      svr_cwp.coupling[s].n_user_vtx = (int *) malloc(sizeof(int) * svr_cwp.coupling[s].n_part);
    }

    svr_cwp.coupling[s].n_user_vtx[i_part] = n_pts;
  }

  svr_cwp.coupling[s].usr_coord[i_part] = (double *) malloc(sizeof(double) * 3 * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].usr_coord[i_part], sizeof(double) * 3 * n_pts);

  // read global_num
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_USER_TGT_PTS_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_User_tgt_pts_set(local_code_name,
                         cpl_id,
                         i_part,
                         n_pts,
                         svr_cwp.coupling[s].usr_coord[i_part],
                         NULL);
  }
  else {
    svr_cwp.coupling[s].usr_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].usr_global_num, sizeof(CWP_g_num_t) * n_pts);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_USER_TGT_PTS_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_User_tgt_pts_set(local_code_name,
                         cpl_id,
                         i_part,
                         n_pts,
                         svr_cwp.coupling[s].usr_coord[i_part],
                         svr_cwp.coupling[s].usr_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_USER_TGT_PTS_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_finalize
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FINALIZE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FINALIZE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Mesh_interf_finalize(local_code_name,
                           cpl_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FINALIZE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_vtx_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_VTX_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_pts
  int n_pts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_pts, sizeof(int));

  // read coord
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].vtx_coord == NULL) {
    svr_cwp.coupling[s].vtx_coord      = (double **) malloc(sizeof(double *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].vtx_global_num = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].vtx_coord[j_part]      = NULL;
      svr_cwp.coupling[s].vtx_global_num[j_part] = NULL;
    }
  }

  if (svr_cwp.coupling[s].mesh_dynamic == CWP_DYNAMIC_MESH_DEFORMABLE) {
    if (svr_cwp.coupling[s].n_vtx == NULL) {
      svr_cwp.coupling[s].n_vtx = (int *) malloc(sizeof(int) * svr_cwp.coupling[s].n_part);
    }

    svr_cwp.coupling[s].n_vtx[i_part] = n_pts;
  }

  svr_cwp.coupling[s].vtx_coord[i_part] = (double *) malloc(sizeof(double) * 3 * n_pts);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].vtx_coord[i_part], sizeof(double) * 3 * n_pts);

  // read global_num
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {
    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_VTX_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_vtx_set(local_code_name,
                            cpl_id,
                            i_part,
                            n_pts,
                            svr_cwp.coupling[s].vtx_coord[i_part],
                            NULL);
  }
  else {
    svr_cwp.coupling[s].vtx_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_pts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].vtx_global_num[i_part], sizeof(CWP_g_num_t) * n_pts);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_VTX_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_vtx_set(local_code_name,
                            cpl_id,
                            i_part,
                            n_pts,
                            svr_cwp.coupling[s].vtx_coord[i_part],
                            svr_cwp.coupling[s].vtx_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_VTX_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_add
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read block_type
  CWP_Block_t block_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_type, sizeof(CWP_Block_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int block_id = CWP_Mesh_interf_block_add(local_code_name,
                                           cpl_id,
                                           block_type);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send block identifier
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &block_id, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_std_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));;

  // read n_elts
  int n_elts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_elts, sizeof(int));

  // read n_vtx_elt
  CWP_Block_t block_type = CWP_std_block_type_get(local_code_name, cpl_id, block_id);
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &block_type, sizeof(int));
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  int n_vtx_elt;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_vtx_elt, sizeof(int));

  // read connectivity
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].block.find( block_id ) != svr_cwp.coupling[s].block.end()) {
    t_block block = t_block();
    svr_cwp.coupling[s].block.insert(std::make_pair(block_id, block));
  }

  if (svr_cwp.coupling[s].block[block_id].std_connec == NULL) {
    svr_cwp.coupling[s].block[block_id].std_connec      = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].std_global_num  = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t*) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].block[block_id].std_connec[j_part]     = NULL;
      svr_cwp.coupling[s].block[block_id].std_global_num[j_part] = NULL;
    }
  }

  svr_cwp.coupling[s].block[block_id].std_connec[i_part] = (int *) malloc(sizeof(int) * n_elts * n_vtx_elt);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].std_connec[i_part], sizeof(int) * n_elts * n_vtx_elt);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_block_std_set(local_code_name,
                                  cpl_id,
                                  i_part,
                                  block_id,
                                  n_elts,
                                  svr_cwp.coupling[s].block[block_id].std_connec[i_part],
                                  NULL);
  }
  else {
    svr_cwp.coupling[s].block[block_id].std_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].std_global_num[i_part], sizeof(CWP_g_num_t) * n_elts);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_block_std_set(local_code_name,
                                  cpl_id,
                                  i_part,
                                  block_id,
                                  n_elts,
                                  svr_cwp.coupling[s].block[block_id].std_connec[i_part],
                                  svr_cwp.coupling[s].block[block_id].std_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_std_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));;

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  int               n_elts     = -1;
  int              *connec     = NULL;
  CWP_g_num_t      *global_num = NULL;
  CWP_Mesh_interf_block_std_get(local_code_name,
                                cpl_id,
                                i_part,
                                block_id,
                                &n_elts,
                                &connec,
                                &global_num);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send n_elts
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_elts, sizeof(int));

  // send connectivity
  CWP_Block_t block_type = CWP_std_block_type_get(local_code_name, cpl_id, block_id);
  int n_vtx_elt = n_vtx_get(block_type);

  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_vtx_elt, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec, sizeof(int) * (n_vtx_elt * n_elts));

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_std_block_type_get
(
 p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STD_BLOCK_TYPE_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));;

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STD_BLOCK_TYPE_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Block_t block_type = CWP_std_block_type_get(local_code_name, cpl_id, block_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_STD_BLOCK_TYPE_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send block_type
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &block_type, sizeof(int));

    // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_f_poly_block_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // read n_elts
  int n_elts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_elts, sizeof(int));

  // read connectivity index
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].block.find( block_id ) != svr_cwp.coupling[s].block.end()) {
    t_block block = t_block();
    svr_cwp.coupling[s].block.insert(std::make_pair(block_id, block));
  }

  if (svr_cwp.coupling[s].block[block_id].connec_idx == NULL) {
    svr_cwp.coupling[s].block[block_id].connec_idx      = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].connec          = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].elt_global_num  = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t*) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].block[block_id].connec_idx[j_part]     = NULL;
      svr_cwp.coupling[s].block[block_id].connec[j_part]         = NULL;
      svr_cwp.coupling[s].block[block_id].elt_global_num[j_part] = NULL;
    }
  }

  svr_cwp.coupling[s].block[block_id].connec_idx[i_part] = (int *) malloc(sizeof(int) * (n_elts+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].connec_idx[i_part], sizeof(int) * (n_elts+1));

  // read connectivity
  svr_cwp.coupling[s].block[block_id].connec[i_part] = (int *) malloc(sizeof(int) * svr_cwp.coupling[s].block[block_id].connec_idx[i_part][n_elts]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].connec[i_part], sizeof(int) * svr_cwp.coupling[s].block[block_id].connec_idx[i_part][n_elts]);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_f_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     svr_cwp.coupling[s].block[block_id].connec_idx[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec[i_part],
                                     NULL);
  }
  else {
    svr_cwp.coupling[s].block[block_id].elt_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].elt_global_num[i_part], sizeof(CWP_g_num_t) * n_elts);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_f_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     svr_cwp.coupling[s].block[block_id].connec_idx[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec[i_part],
                                     svr_cwp.coupling[s].block[block_id].elt_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_f_poly_block_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int          n_elts     = -1;
  int         *connec_idx = NULL;
  int         *connec     = NULL;
  CWP_g_num_t *global_num = NULL;
  CWP_Mesh_interf_f_poly_block_get(local_code_name,
                                   cpl_id,
                                   i_part,
                                   block_id,
                                   &n_elts,
                                   &connec_idx,
                                   &connec,
                                   &global_num);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send n_elts
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_elts, sizeof(int));

  // send connectivity index
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_idx, sizeof(int) * (n_elts+1));

  // send connectivity
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec, sizeof(int) * connec_idx[n_elts]);

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_c_poly_block_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // read n_elts
  int n_elts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_elts, sizeof(int));

  // read n_faces
  int n_faces;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_faces, sizeof(int));

  // read connectivity faces index
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].block.find( block_id ) != svr_cwp.coupling[s].block.end()) {
    t_block block = t_block();
    svr_cwp.coupling[s].block.insert(std::make_pair(block_id, block));
  }

  if (svr_cwp.coupling[s].block[block_id].connec_faces_idx == NULL) {
    svr_cwp.coupling[s].block[block_id].connec_faces_idx = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].connec_faces     = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].connec_cells_idx = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].connec_cells     = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].cell_global_num  = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t*) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].block[block_id].connec_faces_idx[j_part] = NULL;
      svr_cwp.coupling[s].block[block_id].connec_faces[j_part]     = NULL;
      svr_cwp.coupling[s].block[block_id].connec_cells_idx[j_part] = NULL;
      svr_cwp.coupling[s].block[block_id].connec_cells[j_part]     = NULL;
      svr_cwp.coupling[s].block[block_id].cell_global_num[j_part]  = NULL;
    }
  }

  svr_cwp.coupling[s].block[block_id].connec_faces_idx[i_part] = (int *) malloc(sizeof(int) * (n_faces+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].connec_faces_idx[i_part], sizeof(int) * (n_faces+1));

  // read connectivity faces
  svr_cwp.coupling[s].block[block_id].connec_faces[i_part] = (int *) malloc(sizeof(int) * svr_cwp.coupling[s].block[block_id].connec_faces_idx[i_part][n_faces]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].connec_faces[i_part], sizeof(int) * svr_cwp.coupling[s].block[block_id].connec_faces_idx[i_part][n_faces]);

  // read connectivity index
  svr_cwp.coupling[s].block[block_id].connec_cells_idx[i_part] = (int *) malloc(sizeof(int) * (n_elts+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].connec_cells_idx[i_part], sizeof(int) * (n_elts+1));

  // read connectivity
  svr_cwp.coupling[s].block[block_id].connec_cells[i_part] = (int *) malloc(sizeof(int) * svr_cwp.coupling[s].block[block_id].connec_cells_idx[i_part][n_elts]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].connec_cells[i_part], sizeof(int) * svr_cwp.coupling[s].block[block_id].connec_cells_idx[i_part][n_elts]);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_c_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     n_faces,
                                     svr_cwp.coupling[s].block[block_id].connec_faces_idx[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec_faces[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec_cells_idx[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec_cells[i_part],
                                     NULL);
  }
  else {
    svr_cwp.coupling[s].block[block_id].cell_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].cell_global_num[i_part], sizeof(CWP_g_num_t) * n_elts);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_c_poly_block_set(local_code_name,
                                     cpl_id,
                                     i_part,
                                     block_id,
                                     n_elts,
                                     n_faces,
                                     svr_cwp.coupling[s].block[block_id].connec_faces_idx[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec_faces[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec_cells_idx[i_part],
                                     svr_cwp.coupling[s].block[block_id].connec_cells[i_part],
                                     svr_cwp.coupling[s].block[block_id].cell_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_c_poly_block_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int          n_elts           = -1;
  int          n_faces          = -1;
  int         *connec_faces_idx = NULL;
  int         *connec_faces     = NULL;
  int         *connec_cells_idx = NULL;
  int         *connec_cells     = NULL;
  CWP_g_num_t *global_num       = NULL;
  CWP_Mesh_interf_c_poly_block_get(local_code_name,
                                   cpl_id,
                                   i_part,
                                   block_id,
                                   &n_elts,
                                   &n_faces,
                                   &connec_faces_idx,
                                   &connec_faces,
                                   &connec_cells_idx,
                                   &connec_cells,
                                   &global_num);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send n_elts
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_elts, sizeof(int));

  // send n_faces
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_faces, sizeof(int));

  // send connectivity index
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_faces_idx, sizeof(int) * (n_faces+1));

  // send connectivity
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_faces, sizeof(int) * connec_faces_idx[n_faces]);

  // send connectivity index
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_cells_idx, sizeof(int) * (n_elts+1));

  // send connectivity
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec_cells, sizeof(int) * connec_cells_idx[n_elts]);

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }

  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_del
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_DEL);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_DEL);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Mesh_interf_del(local_code_name,
                      cpl_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_DEL);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].ijk_grid != NULL) {
    free(svr_cwp.coupling[s].ijk_grid);
    svr_cwp.coupling[s].ijk_grid = NULL;
  }

  // WARNING : might call this then time update but since coupling still exists
  // will try to write in a NULL pointer
  if (svr_cwp.coupling[s].mesh_dynamic == CWP_DYNAMIC_MESH_DEFORMABLE) {
    if (svr_cwp.coupling[s].n_user_vtx != NULL) {
      free(svr_cwp.coupling[s].n_user_vtx);
      svr_cwp.coupling[s].n_user_vtx = NULL;
    }
    if (svr_cwp.coupling[s].n_vtx != NULL) {
      free(svr_cwp.coupling[s].n_vtx);
      svr_cwp.coupling[s].n_vtx = NULL;
    }
  }

  // n_part
  if (svr_cwp.coupling[s].vtx_coord != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].vtx_coord[i_part] != NULL) {
        free(svr_cwp.coupling[s].vtx_coord[i_part]);
        svr_cwp.coupling[s].vtx_coord[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].vtx_coord);
    svr_cwp.coupling[s].vtx_coord = NULL;
  }

  if (svr_cwp.coupling[s].vtx_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].vtx_global_num[i_part] != NULL) {
        free(svr_cwp.coupling[s].vtx_global_num[i_part]);
        svr_cwp.coupling[s].vtx_global_num[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].vtx_global_num);
    svr_cwp.coupling[s].vtx_global_num = NULL;
  }

  if (svr_cwp.coupling[s].usr_coord != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].usr_coord[i_part] != NULL) {
        free(svr_cwp.coupling[s].usr_coord[i_part]);
        svr_cwp.coupling[s].usr_coord[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].usr_coord);
    svr_cwp.coupling[s].usr_coord = NULL;
  }

  if (svr_cwp.coupling[s].usr_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].usr_global_num[i_part] != NULL) free(svr_cwp.coupling[s].usr_global_num[i_part]);
    }
    free(svr_cwp.coupling[s].usr_global_num);
  }

  if (svr_cwp.coupling[s].face_edge_idx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_edge_idx[i_part] != NULL) {
        free(svr_cwp.coupling[s].face_edge_idx[i_part]);
        svr_cwp.coupling[s].face_edge_idx[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].face_edge_idx);
    svr_cwp.coupling[s].face_edge_idx = NULL;
  }

  if (svr_cwp.coupling[s].face_edge != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_edge[i_part] != NULL) {
        free(svr_cwp.coupling[s].face_edge[i_part]);
        svr_cwp.coupling[s].face_edge[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].face_edge);
    svr_cwp.coupling[s].face_edge = NULL;
  }

  if (svr_cwp.coupling[s].edge_vtx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].edge_vtx[i_part] != NULL) {
        free(svr_cwp.coupling[s].edge_vtx[i_part]);
        svr_cwp.coupling[s].edge_vtx[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].edge_vtx);
    svr_cwp.coupling[s].edge_vtx = NULL;
  }

  if (svr_cwp.coupling[s].face_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_global_num[i_part] != NULL) {
        free(svr_cwp.coupling[s].face_global_num[i_part]);
        svr_cwp.coupling[s].face_global_num[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].face_global_num);
    svr_cwp.coupling[s].face_global_num = NULL;
  }

  if (svr_cwp.coupling[s].cell_face_idx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].cell_face_idx[i_part] != NULL) {
        free(svr_cwp.coupling[s].cell_face_idx[i_part]);
        svr_cwp.coupling[s].cell_face_idx[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].cell_face_idx);
    svr_cwp.coupling[s].cell_face_idx = NULL;
  }

  if (svr_cwp.coupling[s].cell_face != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].cell_face[i_part] != NULL) {
        free(svr_cwp.coupling[s].cell_face[i_part]);
        svr_cwp.coupling[s].cell_face[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].cell_face);
    svr_cwp.coupling[s].cell_face = NULL;
  }

  if (svr_cwp.coupling[s].face_vtx_idx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_vtx_idx[i_part] != NULL) free(svr_cwp.coupling[s].face_vtx_idx[i_part]);
    }
    free(svr_cwp.coupling[s].face_vtx_idx);
    svr_cwp.coupling[s].face_vtx_idx = NULL;
  }

  if (svr_cwp.coupling[s].face_vtx != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].face_vtx[i_part] != NULL) {
        free(svr_cwp.coupling[s].face_vtx[i_part]);
        svr_cwp.coupling[s].face_vtx[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].face_vtx);
    svr_cwp.coupling[s].face_vtx = NULL;
  }

  if (svr_cwp.coupling[s].cell_global_num != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
      if (svr_cwp.coupling[s].cell_global_num[i_part] != NULL) {
        free(svr_cwp.coupling[s].cell_global_num[i_part]);
        svr_cwp.coupling[s].cell_global_num[i_part] = NULL;
      }
    }
    free(svr_cwp.coupling[s].cell_global_num);
    svr_cwp.coupling[s].cell_global_num = NULL;
  }

  // Block
  if (!svr_cwp.coupling[s].block.empty()) {
    std::map<int, t_block>::iterator it_b = svr_cwp.coupling[s].block.begin();
    while (it_b != svr_cwp.coupling[s].block.end()) {

      // n_part
      if ((it_b->second).connec_faces_idx != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_faces_idx[i_part] != NULL) {
            free((it_b->second).connec_faces_idx[i_part]);
            (it_b->second).connec_faces_idx[i_part] = NULL;
          }
        }
        free((it_b->second).connec_faces_idx);
        (it_b->second).connec_faces_idx = NULL;
      }

      if ((it_b->second).connec_faces != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_faces[i_part] != NULL) {
            free((it_b->second).connec_faces[i_part]);
            (it_b->second).connec_faces[i_part] = NULL;
          }
        }
        free((it_b->second).connec_faces);
        (it_b->second).connec_faces = NULL;
      }

      if ((it_b->second).connec_cells_idx != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_cells_idx[i_part] != NULL) {
            free((it_b->second).connec_cells_idx[i_part]);
            (it_b->second).connec_cells_idx[i_part] = NULL;
          }
        }
        free((it_b->second).connec_cells_idx);
        (it_b->second).connec_cells_idx = NULL;
      }

      if ((it_b->second).connec_cells != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_cells[i_part] != NULL) free((it_b->second).connec_cells[i_part]);
        }
        free((it_b->second).connec_cells);
      }

      if ((it_b->second).cell_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).cell_global_num[i_part] != NULL) {
            free((it_b->second).cell_global_num[i_part]);
            (it_b->second).cell_global_num[i_part] = NULL;
          }
        }
        free((it_b->second).cell_global_num);
        (it_b->second).cell_global_num = NULL;
      }

      if ((it_b->second).connec_idx != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec_idx[i_part] != NULL) {
            free((it_b->second).connec_idx[i_part]);
            (it_b->second).connec_idx[i_part] = NULL;
          }
        }
        free((it_b->second).connec_idx);
        (it_b->second).connec_idx = NULL;
      }

      if ((it_b->second).connec != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).connec[i_part] != NULL) {
            free((it_b->second).connec[i_part]);
            (it_b->second).connec[i_part] = NULL;
          }
        }
        free((it_b->second).connec);
        (it_b->second).connec = NULL;
      }

      if ((it_b->second).elt_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).elt_global_num[i_part] != NULL) {
            free((it_b->second).elt_global_num[i_part]);
            (it_b->second).elt_global_num[i_part] = NULL;
          }
        }
        free((it_b->second).elt_global_num);
        (it_b->second).elt_global_num = NULL;
      }

      if ((it_b->second).std_connec != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).std_connec[i_part] != NULL) {
            free((it_b->second).std_connec[i_part]);
            (it_b->second).std_connec[i_part] = NULL;
          }
        }
        free((it_b->second).std_connec);
        (it_b->second).std_connec = NULL;
      }

      if ((it_b->second).std_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).std_global_num[i_part] != NULL) {
            free((it_b->second).std_global_num[i_part]);
            (it_b->second).std_global_num[i_part] =  NULL;
          }
        }
        free((it_b->second).std_global_num);
        (it_b->second).std_global_num = NULL;
      }

      if ((it_b->second).ho_std_connec != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).ho_std_connec[i_part] != NULL) {
            free((it_b->second).ho_std_connec[i_part]);
            (it_b->second).ho_std_connec[i_part] = NULL;
          }
        }
        free((it_b->second).ho_std_connec);
        (it_b->second).ho_std_connec = NULL;
      }

      if ((it_b->second).ho_std_global_num != NULL) {
        for (int i_part = 0; i_part < svr_cwp.coupling[s].n_part; i_part++) {
          if ((it_b->second).ho_std_global_num[i_part] != NULL) {
            free((it_b->second).ho_std_global_num[i_part]);
            (it_b->second).ho_std_global_num[i_part] = NULL;
          }
        }
        free((it_b->second).ho_std_global_num);
        (it_b->second).ho_std_global_num = NULL;
      }

      it_b = svr_cwp.coupling[s].block.erase(it_b);
    }
  }

  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_from_cellface_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_cells
  int n_cells;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_cells, sizeof(int));

  // read connectivity cells index
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].cell_face_idx == NULL) {
    svr_cwp.coupling[s].cell_face_idx   = (int **) malloc(sizeof(int *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].cell_face       = (int **) malloc(sizeof(int *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].face_vtx_idx    = (int **) malloc(sizeof(int *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].face_vtx        = (int **) malloc(sizeof(int *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].cell_global_num = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].cell_face_idx[j_part]   = NULL;
      svr_cwp.coupling[s].cell_face[j_part]       = NULL;
      svr_cwp.coupling[s].face_vtx_idx[j_part]    = NULL;
      svr_cwp.coupling[s].face_vtx[j_part]        = NULL;
      svr_cwp.coupling[s].cell_global_num[j_part] = NULL;
    }
  }

  // read connectivity cells index
  svr_cwp.coupling[s].cell_face_idx[i_part] = (int *) malloc(sizeof(int) * (n_cells+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].cell_face_idx[i_part], sizeof(int) * (n_cells+1));

  // read connectivity cells
  svr_cwp.coupling[s].cell_face[i_part] = (int *) malloc(sizeof(int) * svr_cwp.coupling[s].cell_face_idx[i_part][n_cells]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].cell_face[i_part], sizeof(int) * svr_cwp.coupling[s].cell_face_idx[i_part][n_cells]);

  // read n_faces
  int n_faces;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_faces, sizeof(int));

  // read connectivity faces index
  svr_cwp.coupling[s].face_vtx_idx[i_part] = (int *) malloc(sizeof(int) * (n_faces+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].face_vtx_idx[i_part], sizeof(int) * (n_faces+1));

  // read connectivity faces
  svr_cwp.coupling[s].face_vtx[i_part] = (int *) malloc(sizeof(int) * svr_cwp.coupling[s].face_vtx_idx[i_part][n_faces]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].face_vtx[i_part], sizeof(int) * svr_cwp.coupling[s].face_vtx_idx[i_part][n_faces]);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_from_cellface_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_cells,
                                      svr_cwp.coupling[s].cell_face_idx[i_part],
                                      svr_cwp.coupling[s].cell_face[i_part],
                                      n_faces,
                                      svr_cwp.coupling[s].face_vtx_idx[i_part],
                                      svr_cwp.coupling[s].face_vtx[i_part],
                                      NULL);
  }
  else {
    svr_cwp.coupling[s].cell_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_cells);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].cell_global_num[i_part], sizeof(CWP_g_num_t) * n_cells);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_from_cellface_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_cells,
                                      svr_cwp.coupling[s].cell_face_idx[i_part],
                                      svr_cwp.coupling[s].cell_face[i_part],
                                      n_faces,
                                      svr_cwp.coupling[s].face_vtx_idx[i_part],
                                      svr_cwp.coupling[s].face_vtx[i_part],
                                      svr_cwp.coupling[s].cell_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_from_faceedge_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read n_cells
  int n_faces;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_faces, sizeof(int));

  // read connectivity faces index
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].face_edge_idx == NULL) {
    svr_cwp.coupling[s].face_edge_idx   = (int **) malloc(sizeof(int *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].face_edge       = (int **) malloc(sizeof(int *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].edge_vtx        = (int **) malloc(sizeof(int *) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].face_global_num = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].face_edge_idx[j_part]   = NULL;
      svr_cwp.coupling[s].face_edge[j_part]       = NULL;
      svr_cwp.coupling[s].edge_vtx[j_part]        = NULL;
      svr_cwp.coupling[s].face_global_num[j_part] = NULL;
    }
  }

  svr_cwp.coupling[s].face_edge_idx[i_part] = (int *) malloc(sizeof(int) * (n_faces+1));
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*)  svr_cwp.coupling[s].face_edge_idx[i_part], sizeof(int) * (n_faces+1));

  // read connectivity faces
  svr_cwp.coupling[s].face_edge[i_part] = (int *) malloc(sizeof(int) *  svr_cwp.coupling[s].face_edge_idx[i_part][n_faces]);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*)  svr_cwp.coupling[s].face_edge[i_part], sizeof(int) *  svr_cwp.coupling[s].face_edge_idx[i_part][n_faces]);

  // read n_edges
  int n_edges;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_edges, sizeof(int));

  // read connectivity edges
  svr_cwp.coupling[s].edge_vtx[i_part] = (int *) malloc(sizeof(int) * 2 * n_edges);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*)  svr_cwp.coupling[s].edge_vtx[i_part], sizeof(int) * 2 * n_edges);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_from_faceedge_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_faces,
                                      svr_cwp.coupling[s].face_edge_idx[i_part],
                                      svr_cwp.coupling[s].face_edge[i_part],
                                      n_edges,
                                      svr_cwp.coupling[s].edge_vtx[i_part],
                                      NULL);
  }
  else {
    svr_cwp.coupling[s].face_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_faces);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*)  svr_cwp.coupling[s].face_global_num[i_part], sizeof(CWP_g_num_t) * n_faces);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_from_faceedge_set(local_code_name,
                                      cpl_id,
                                      i_part,
                                      n_faces,
                                      svr_cwp.coupling[s].face_edge_idx[i_part],
                                      svr_cwp.coupling[s].face_edge[i_part],
                                      n_edges,
                                      svr_cwp.coupling[s].edge_vtx[i_part],
                                      svr_cwp.coupling[s].face_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_create
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_CREATE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read data_type
  CWP_Type_t data_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &data_type, sizeof(CWP_Type_t));

  // read storage
  CWP_Field_storage_t storage;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &storage, sizeof(CWP_Field_storage_t));

  // read n_component
  int n_component;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_component, sizeof(int));

  // read target_location
  CWP_Dof_location_t target_location;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &target_location, sizeof(CWP_Dof_location_t));

  // read exch_type
  CWP_Field_exch_t exch_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &exch_type, sizeof(CWP_Field_exch_t));

  // read visu_status
  CWP_Status_t visu_status;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &visu_status, sizeof(CWP_Status_t));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_CREATE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Field_create(local_code_name,
                   cpl_id,
                   field_id,
                   data_type,
                   storage,
                   n_component,
                   target_location,
                   exch_type,
                   visu_status);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_CREATE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // create occurence in map
  std::string s1(cpl_id);
  std::string s2(field_id);
  t_field field = t_field();
  svr_cwp.coupling[s1].field.insert(std::make_pair(s2, field));
  svr_cwp.coupling[s1].field[s2].data = (double **) malloc(sizeof(double*) * svr_cwp.coupling[s1].n_part);
  svr_cwp.coupling[s1].field[s2].size = (int *) malloc(sizeof(int) * svr_cwp.coupling[s1].n_part);
  for (int i_part = 0; i_part < svr_cwp.coupling[s1].n_part; i_part++) {
    svr_cwp.coupling[s1].field[s2].data[i_part] = NULL;
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_data_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DATA_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read map_type
  std::string s1(cpl_id);
  std::string s2(field_id);
  CWP_Field_map_t map_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &map_type, sizeof(CWP_Field_map_t));
  svr_cwp.coupling[s1].field[s2].map_type = map_type;

  // read data array
  int size;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &size, sizeof(int));
  svr_cwp.coupling[s1].field[s2].data[i_part] = (double *) realloc(svr_cwp.coupling[s1].field[s2].data[i_part],
                                                                   sizeof(double) * size);
  svr_cwp.coupling[s1].field[s2].size[i_part] = size;
  if (map_type == CWP_FIELD_MAP_SOURCE) {
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s1].field[s2].data[i_part], sizeof(double) * size);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DATA_SET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Field_data_set(local_code_name,
                     cpl_id,
                     field_id,
                     i_part,
                     map_type,
                     svr_cwp.coupling[s1].field[s2].data[i_part]);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DATA_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_dof_location_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DOF_LOCATION_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DOF_LOCATION_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int dof_location = CWP_Field_dof_location_get(local_code_name,
                                                       cpl_id,
                                                       field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DOF_LOCATION_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send degree of freedom location
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &dof_location, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_storage_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_STORAGE_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_STORAGE_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int storage_type = CWP_Field_storage_get(local_code_name,
                                           cpl_id,
                                           field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_STORAGE_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number storage type
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &storage_type, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_n_components_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_N_COMPONENTS_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_N_COMPONENTS_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int n_components = CWP_Field_n_components_get(local_code_name,
                                                cpl_id,
                                                field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_N_COMPONENTS_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of components
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_components, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_n_dof_get
(
  p_server                 svr
)
{
    // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_N_DOF_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_N_DOF_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  int n_dof = CWP_Field_n_dof_get(local_code_name,
                                  cpl_id,
                                  field_id,
                                  i_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_N_DOF_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send number of degrees of freedom
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_dof, sizeof(int));

  // free
  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_del
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DEL);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read field identifier
  char *field_id = (char *) malloc(sizeof(char));
  read_name(&field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DEL);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Field_del(local_code_name,
                cpl_id,
                field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_DEL);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  std::string s1(cpl_id);
  std::string s2(field_id);

  if (svr_cwp.coupling[s1].field[s2].size != NULL) {
    free(svr_cwp.coupling[s1].field[s2].size);
    svr_cwp.coupling[s1].field[s2].size = NULL;
  }

  if (svr_cwp.coupling[s1].field[s2].data != NULL) {
    for (int i_part = 0; i_part < svr_cwp.coupling[s1].n_part; i_part++) {
      if (svr_cwp.coupling[s1].field[s2].data[i_part] != NULL) free(svr_cwp.coupling[s1].field[s2].data[i_part]);
      svr_cwp.coupling[s1].field[s2].data[i_part] = NULL;
    }
    free(svr_cwp.coupling[s1].field[s2].data);
    svr_cwp.coupling[s1].field[s2].data = NULL;
  }

  svr_cwp.coupling[s1].field.erase(s2);

  free(local_code_name);
  free(cpl_id);
  free(field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_issend
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_ISSEND);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read source field identifier
  char *src_field_id = (char *) malloc(sizeof(char));
  read_name(&src_field_id, svr);

  // receive field data from client to update data pointer
  std::string s1(cpl_id);
  std::string s2(src_field_id);
  for (int j_part = 0; j_part < svr_cwp.coupling[s1].n_part; j_part++) {
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) (svr_cwp.coupling[s1].field[s2]).data[j_part], sizeof(double) * (svr_cwp.coupling[s1].field[s2]).size[j_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_ISSEND);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Field_issend(local_code_name,
                   cpl_id,
                   src_field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_ISSEND);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(src_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_irecv
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_IRECV);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read target field identifier
  char *tgt_field_id = (char *) malloc(sizeof(char));
  read_name(&tgt_field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_IRECV);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Field_irecv(local_code_name,
                  cpl_id,
                  tgt_field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_IRECV);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(tgt_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_wait_issend
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_WAIT_ISSEND);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read source field identifier
  char *src_field_id = (char *) malloc(sizeof(char));
  read_name(&src_field_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_WAIT_ISSEND);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Field_wait_issend(local_code_name,
                        cpl_id,
                        src_field_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_WAIT_ISSEND);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(src_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Field_wait_irecv
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_WAIT_IRECV);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read target field identifier
  char *tgt_field_id = (char *) malloc(sizeof(char));
  read_name(&tgt_field_id, svr);

  // read map_type
  CWP_Field_map_t map_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &map_type, sizeof(CWP_Field_map_t));

  // read n_part
  int n_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_part, sizeof(int));

  // read n_components
  int n_components;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_components, sizeof(int));

  // read n_entities
  int *n_entities = (int *) malloc(sizeof(int) * n_part);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) n_entities, sizeof(int) * n_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_WAIT_IRECV);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // launch
  CWP_Field_wait_irecv(local_code_name,
                       cpl_id,
                       tgt_field_id);

  double **data = (double **) malloc(sizeof(double *) * n_part);
  int *n_computed_tgts = (int *) malloc(sizeof(int) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    if (n_entities[i_part] != -1) {

      n_computed_tgts[i_part] = CWP_N_computed_tgts_get(local_code_name,
                                                        cpl_id,
                                                        tgt_field_id,
                                                        i_part);

    }
    data[i_part] = NULL;
  }

  // launch CWP_Field_data_get to retreive field and send to Client
  for (int i_part = 0; i_part < n_part; i_part++) {
    if (n_entities[i_part] != -1) {
      CWP_Field_data_get(local_code_name,
                         cpl_id,
                         tgt_field_id,
                         i_part,
                         map_type,
                         &data[i_part]);
    }
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_FIELD_WAIT_IRECV);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send data
  svr->state = CWP_SVRSTATE_SENDPGETDATA;
  // send sizes
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) n_computed_tgts, sizeof(int) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    if (n_entities[i_part] != -1) {
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) data[i_part], sizeof(double) * n_computed_tgts[i_part] * n_components);
    }
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(tgt_field_id);
  if (n_computed_tgts != NULL) free(n_computed_tgts);
  if (n_entities != NULL) free(n_entities);
  if (data != NULL) {
    free(data);
  }

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Interp_function_unset
(
  p_server                 svr
)
{
  // wait all ranks have receive msg
  MPI_Barrier(svr_mpi.intra_comms[0]);

  // read local code name
  svr->state = CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read source field identifier
  char *src_field_id = (char *) malloc(sizeof(char));
  read_name(&src_field_id, svr);

  // launch
  CWP_Field_interp_function_unset(local_code_name,
                                  cpl_id,
                                  src_field_id);

  // free
  free(local_code_name);
  free(cpl_id);
  free(src_field_id);

  svr->state = CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Interp_function_set
(
  p_server                 svr
)
{
  PDM_UNUSED(svr);
  // TO DO: create some standard interpolation function on server side for user to choose
}

void
CWP_server_Mesh_interf_block_ho_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // read n_elts
  int n_elts;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_elts, sizeof(int));

  // read order
  int order;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &order, sizeof(int));

  // read n_vtx_elt
  CWP_Block_t block_type = CWP_std_block_type_get(local_code_name, cpl_id, block_id);
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &block_type, sizeof(int));
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  int n_vtx_elt;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_vtx_elt, sizeof(int));

  // read connectivity
  std::string s(cpl_id);

  if (svr_cwp.coupling[s].block.find( block_id ) != svr_cwp.coupling[s].block.end()) {
    t_block block = t_block();
    svr_cwp.coupling[s].block.insert(std::make_pair(block_id, block));
  }

  if (svr_cwp.coupling[s].block[block_id].ho_std_connec == NULL) {
    svr_cwp.coupling[s].block[block_id].ho_std_connec      = (int **) malloc(sizeof(int*) * svr_cwp.coupling[s].n_part);
    svr_cwp.coupling[s].block[block_id].ho_std_global_num  = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t*) * svr_cwp.coupling[s].n_part);

    for (int j_part = 0; j_part < svr_cwp.coupling[s].n_part; j_part++) {
      svr_cwp.coupling[s].block[block_id].ho_std_connec[j_part]     = NULL;
      svr_cwp.coupling[s].block[block_id].ho_std_global_num[j_part] = NULL;
    }
  }


  svr_cwp.coupling[s].block[block_id].ho_std_connec[i_part] = (int *) malloc(sizeof(int) * n_elts * n_vtx_elt);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].ho_std_connec[i_part], sizeof(int) * n_elts * n_vtx_elt);

  // read global number
  int NULL_flag;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  // launch
  if (NULL_flag) {

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_block_ho_set(local_code_name,
                                  cpl_id,
                                  i_part,
                                  block_id,
                                  n_elts,
                                  order,
                                  svr_cwp.coupling[s].block[block_id].ho_std_connec[i_part],
                                  NULL);
  }
  else {
    svr_cwp.coupling[s].block[block_id].ho_std_global_num[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elts);
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].block[block_id].ho_std_global_num[i_part], sizeof(CWP_g_num_t) * n_elts);

    // send status msg
    MPI_Barrier(svr_mpi.intra_comms[0]);
    if (svr->flags & CWP_FLAG_VERBOSE) {
      t_message message;
      NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_SET);
      message.flag = CWP_SVR_LCH_BEGIN;
      CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
    }

    CWP_Mesh_interf_block_ho_set(local_code_name,
                                  cpl_id,
                                  i_part,
                                  block_id,
                                  n_elts,
                                  order,
                                  svr_cwp.coupling[s].block[block_id].ho_std_connec[i_part],
                                  svr_cwp.coupling[s].block[block_id].ho_std_global_num[i_part]);
  }

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_block_ho_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read i_part
  int i_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &i_part, sizeof(int));

  // read block_id
  int block_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_id, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  int               n_elts     = -1;
  int               order      = -1;
  int              *connec     = NULL;
  CWP_g_num_t      *global_num = NULL;
  CWP_Mesh_interf_block_ho_get(local_code_name,
                               cpl_id,
                               i_part,
                               block_id,
                               &n_elts,
                               &order,
                               &connec,
                               &global_num);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send n_elts
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_elts, sizeof(int));

  // send order
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &order, sizeof(int));

  // send connectivity
  CWP_Block_t block_type = CWP_std_block_type_get(local_code_name, cpl_id, block_id);
  int n_vtx_elt = ho_n_vtx_get(block_type, order);

  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &n_vtx_elt, sizeof(int));
  CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, connec, sizeof(int) * (n_vtx_elt * n_elts));

  // send global number
  int NULL_flag = 0;
  if (global_num == NULL) {
    NULL_flag = 1;
  }
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,(void*) &NULL_flag, sizeof(int));

  if (!NULL_flag) {
    CWP_transfer_writedata(svr->connected_socket, svr->max_msg_size, global_num, sizeof(CWP_g_num_t) * n_elts);
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Mesh_interf_ho_ordering_from_IJK_set
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_HO_ORDERING_FROM_IJK_SET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read block_type
  CWP_Block_t block_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &block_type, sizeof(CWP_Block_t));

  // read order
  int order;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &order, sizeof(int));

  // read n_nodes
  int n_nodes;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_nodes, sizeof(int));

  // read connectivity
  std::string s(cpl_id);
  int size;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &size, sizeof(int));
  svr_cwp.coupling[s].ijk_grid = (int *) malloc(sizeof(int) * size);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s].ijk_grid, sizeof(int) * size);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_HO_ORDERING_FROM_IJK_SET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Mesh_interf_ho_ordering_from_IJK_set(local_code_name,
                                           cpl_id,
                                           block_type,
                                           order,
                                           n_nodes,
                                           svr_cwp.coupling[s].ijk_grid);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_MESH_INTERF_HO_ORDERING_FROM_IJK_SET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Cpl_spatial_interp_algo_get
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_SPATIAL_INTERP_ALGO_GET);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_SPATIAL_INTERP_ALGO_GET);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Spatial_interp_t spatial_interp_algo = CWP_Cpl_spatial_interp_algo_get(local_code_name,
                                                                             cpl_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_CPL_SPATIAL_INTERP_ALGO_GET);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send block_type
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) &spatial_interp_algo, sizeof(CWP_Spatial_interp_t));

  // free
  free(local_code_name);
  free(cpl_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Global_data_issend
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_ISSEND);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

   // read global data identifier
  char *global_data_id = (char *) malloc(sizeof(char));
  read_name(&global_data_id, svr);

  // read s_send_entity
  size_t s_send_entity;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &s_send_entity, sizeof(size_t));

  // read send_stride
  int send_stride;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &send_stride, sizeof(int));

  // read n_send_entity
  int n_send_entity;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_send_entity, sizeof(int));

  // create occurence if global_data
  std::string s1(cpl_id);
  t_coupling coupling = svr_cwp.coupling[s1];
  std::string s2(global_data_id);
  t_global_data global_data = t_global_data();
  coupling.global_data.insert(std::make_pair(s2, global_data));

  // allocate array
  svr_cwp.coupling[s1].global_data[s2].send_data = malloc(s_send_entity * send_stride * n_send_entity);

  // read connectivity
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) svr_cwp.coupling[s1].global_data[s2].send_data, s_send_entity * send_stride * n_send_entity);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_ISSEND);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Global_data_issend(local_code_name,
                         cpl_id,
                         global_data_id,
                         s_send_entity,
                         send_stride,
                         n_send_entity,
                         svr_cwp.coupling[s1].global_data[s2].send_data);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_ISSEND);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(global_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Global_data_irecv
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_IRECV);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

   // read global data identifier
  char *global_data_id = (char *) malloc(sizeof(char));
  read_name(&global_data_id, svr);

  // read s_recv_entity
  size_t s_recv_entity;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &s_recv_entity, sizeof(size_t));

  // read recv_stride
  int recv_stride;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &recv_stride, sizeof(int));

  // read n_recv_entity
  int n_recv_entity;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_recv_entity, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_IRECV);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // create occurence if global_data
  // WARNING: created here and in send because functions should be run on different servers
  std::string s1(cpl_id);
  t_coupling coupling = svr_cwp.coupling[s1];
  std::string s2(global_data_id);
  t_global_data global_data = t_global_data();
  coupling.global_data.insert(std::make_pair(s2, global_data));

  // allocate array
  svr_cwp.coupling[s1].global_data[s2].s_recv_entity = s_recv_entity;
  svr_cwp.coupling[s1].global_data[s2].recv_stride = recv_stride;
  svr_cwp.coupling[s1].global_data[s2].n_recv_entity = n_recv_entity;
  svr_cwp.coupling[s1].global_data[s2].recv_data = malloc(s_recv_entity * recv_stride * n_recv_entity);

  CWP_Global_data_irecv(local_code_name,
                         cpl_id,
                         global_data_id,
                         s_recv_entity,
                         recv_stride,
                         n_recv_entity,
                         svr_cwp.coupling[s1].global_data[s2].recv_data);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_IRECV);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(global_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Global_data_wait_issend
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_WAIT_ISSEND);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

   // read global data identifier
  char *global_data_id = (char *) malloc(sizeof(char));
  read_name(&global_data_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_WAIT_ISSEND);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Global_data_wait_issend(local_code_name,
                              cpl_id,
                              global_data_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_WAIT_ISSEND);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  std::string s1(cpl_id);
  std::string s2(global_data_id);

  if (svr_cwp.coupling[s1].global_data[s2].send_data != NULL) free(svr_cwp.coupling[s1].global_data[s2].send_data);
  svr_cwp.coupling[s1].global_data[s2].send_data = NULL;

  svr_cwp.coupling[s1].global_data.erase(s2);

  // free
  free(local_code_name);
  free(cpl_id);
  free(global_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Global_data_wait_irecv
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_WAIT_IRECV);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

   // read global data identifier
  char *global_data_id = (char *) malloc(sizeof(char));
  read_name(&global_data_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_WAIT_IRECV);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Global_data_wait_irecv(local_code_name,
                             cpl_id,
                             global_data_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_GLOBAL_DATA_WAIT_IRECV);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send receive data
  // WARNING send here because in recv is not yet allocated
  svr->state=CWP_SVRSTATE_SENDPGETDATA;
  std::string s1(cpl_id);
  t_coupling coupling = svr_cwp.coupling[s1];
  std::string s2(global_data_id);
  CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, (void*) svr_cwp.coupling[s1].global_data[s2].recv_data,
                         svr_cwp.coupling[s1].global_data[s2].s_recv_entity * svr_cwp.coupling[s1].global_data[s2].recv_stride * svr_cwp.coupling[s1].global_data[s2].n_recv_entity);

  // free
  free(local_code_name);
  free(cpl_id);
  free(global_data_id);
  if (svr_cwp.coupling[s1].global_data[s2].recv_data != NULL) free(svr_cwp.coupling[s1].global_data[s2].recv_data);
  svr_cwp.coupling[s1].global_data[s2].recv_data = NULL;
  svr_cwp.coupling[s1].global_data.erase(s2);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Part_data_create
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_CREATE);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read partitionned data identifier
  char *part_data_id = (char *) malloc(sizeof(char));
  read_name(&part_data_id, svr);

  // read exch_type
  CWP_PartData_exch_t exch_type;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &exch_type, sizeof(int));

  // read n_part
  int n_part;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_part, sizeof(int));

  // create occurence if part_data
  std::string s1(cpl_id);
  std::string s2(part_data_id);

  t_coupling &coupling = svr_cwp.coupling[s1];

  t_part_data part_data = t_part_data();

  part_data.exch_type = exch_type;
  part_data.n_part    = n_part;

  // read n_elt
  int *n_elt = (int *) malloc(sizeof(int) * n_part);
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) n_elt, sizeof(int) * n_part);

  part_data.n_elt = n_elt;

  // read gnum
  CWP_g_num_t **gnum_elt = (CWP_g_num_t **) malloc(sizeof(CWP_g_num_t *) * n_part);
  for (int i_part = 0; i_part < n_part; i_part++) {
    gnum_elt[i_part] = (CWP_g_num_t *) malloc(sizeof(CWP_g_num_t) * n_elt[i_part]);
  }
  for (int i_part = 0; i_part < n_part; i_part++) {
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) gnum_elt[i_part], sizeof(CWP_g_num_t) * n_elt[i_part]);
  }

  part_data.gnum_elt = gnum_elt;

  coupling.part_data.insert(std::make_pair(s2, part_data));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_CREATE);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Part_data_create(local_code_name,
                       cpl_id,
                       part_data_id,
                       exch_type,
                       gnum_elt,
                       n_elt,
                       n_part);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_CREATE);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(part_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Part_data_del
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_DEL);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read partitionned data identifier
  char *part_data_id = (char *) malloc(sizeof(char));
  read_name(&part_data_id, svr);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_DEL);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Part_data_del(local_code_name,
                    cpl_id,
                    part_data_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_DEL);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free struct
  std::string s1(cpl_id);
  std::string s2(part_data_id);

  t_part_data &part_data = svr_cwp.coupling[s1].part_data[s2];

  if (part_data.gnum_elt != NULL) {
    for (int i = 0; i < part_data.n_part; i++) {
      if (part_data.gnum_elt[i] != NULL) {
        free(part_data.gnum_elt[i]);
      }
    }
    free(part_data.gnum_elt);
  }

  if (part_data.n_elt != NULL) {
    free(part_data.n_elt);
  }

  svr_cwp.coupling[s1].part_data.erase(s2);

  // free
  free(local_code_name);
  free(cpl_id);
  free(part_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Part_data_issend
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_ISSEND);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read partitionned data identifier
  char *part_data_id = (char *) malloc(sizeof(char));
  read_name(&part_data_id, svr);

  // read exch_id
  int exch_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &exch_id, sizeof(int));

  // read s_data
  size_t s_data;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &s_data, sizeof(size_t));

  // read n_components
  int n_components;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_components, sizeof(int));

  // read part1_to_part2_data
  std::string s1(cpl_id);
  std::string s2(part_data_id);

  t_part_data &part_data = svr_cwp.coupling[s1].part_data[s2];

  assert(part_data.exch_type == CWP_PARTDATA_SEND);

  void **send_data = (void **) malloc(sizeof(void *) * part_data.n_part);
  for (int i_part = 0; i_part < part_data.n_part; i_part++) {
    send_data[i_part] = malloc(s_data * n_components * part_data.n_elt[i_part]);
  }
  for (int i_part = 0; i_part < part_data.n_part; i_part++) {
    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) send_data[i_part], s_data * n_components * part_data.n_elt[i_part]);
  }

  part_data.data.insert(std::make_pair(exch_id, send_data)); // necessary?

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_ISSEND);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Part_data_issend(local_code_name,
                       cpl_id,
                       part_data_id,
                       exch_id,
                       s_data,
                       n_components,
                       send_data);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_ISSEND);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(part_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Part_data_irecv
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_IRECV);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read partitionned data identifier
  char *part_data_id = (char *) malloc(sizeof(char));
  read_name(&part_data_id, svr);

  // read exch_id
  int exch_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &exch_id, sizeof(int));

  // read s_data
  size_t s_data;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &s_data, sizeof(size_t));

  // read n_components
  int n_components;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &n_components, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_IRECV);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // allocate array
  std::string s1(cpl_id);
  std::string s2(part_data_id);

  t_part_data &part_data = svr_cwp.coupling[s1].part_data[s2];
  assert(part_data.exch_type == CWP_PARTDATA_RECV);

  void **recv_data = (void **) malloc(sizeof(void *) * part_data.n_part);
  for (int i_part = 0; i_part < part_data.n_part; i_part++) {
    recv_data[i_part] = malloc(s_data * n_components * part_data.n_elt[i_part]);
  }

  part_data.s_unit.insert(std::make_pair(exch_id, (int) s_data * n_components));
  part_data.data  .insert(std::make_pair(exch_id, recv_data));

  CWP_Part_data_irecv(local_code_name,
                      cpl_id,
                      part_data_id,
                      exch_id,
                      s_data,
                      n_components,
                      recv_data);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_IRECV);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(part_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Part_data_wait_issend
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_WAIT_ISSEND);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read partitionned data identifier
  char *part_data_id = (char *) malloc(sizeof(char));
  read_name(&part_data_id, svr);

  // read exch_id
  int exch_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &exch_id, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_WAIT_ISSEND);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  std::string s1(cpl_id);
  std::string s2(part_data_id);

  t_part_data &part_data = svr_cwp.coupling[s1].part_data[s2];
  assert(part_data.exch_type == CWP_PARTDATA_SEND);


  CWP_Part_data_wait_issend(local_code_name,
                            cpl_id,
                            part_data_id,
                            exch_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_WAIT_ISSEND);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // free
  free(local_code_name);
  free(cpl_id);
  free(part_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

void
CWP_server_Part_data_wait_irecv
(
  p_server                 svr
)
{
  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_WAIT_IRECV);
    message.flag = CWP_SVR_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // read local code name
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  char *local_code_name = (char *) malloc(sizeof(char));
  read_name(&local_code_name, svr);

  // read coupling identifier
  char *cpl_id = (char *) malloc(sizeof(char));
  read_name(&cpl_id, svr);

  // read partitionned data identifier
  char *part_data_id = (char *) malloc(sizeof(char));
  read_name(&part_data_id, svr);

  // read exch_id
  int exch_id;
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,(void*) &exch_id, sizeof(int));

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_WAIT_IRECV);
    message.flag = CWP_SVR_LCH_BEGIN;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  CWP_Part_data_wait_irecv(local_code_name,
                           cpl_id,
                           part_data_id,
                           exch_id);

  // send status msg
  MPI_Barrier(svr_mpi.intra_comms[0]);
  if (svr->flags & CWP_FLAG_VERBOSE) {
    t_message message;
    NEWMESSAGE(message, CWP_MSG_CWP_PART_DATA_WAIT_IRECV);
    message.flag = CWP_SVR_LCH_END;
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, &message, sizeof(t_message));
  }

  // send recv data
  std::string s1(cpl_id);
  std::string s2(part_data_id);

  t_part_data &part_data = svr_cwp.coupling[s1].part_data[s2];
  assert(part_data.exch_type == CWP_PARTDATA_RECV);

  void **data = part_data.data  [exch_id];
  int  s_unit = part_data.s_unit[exch_id];

  for (int i_part = 0; i_part < part_data.n_part; i_part++) {
    CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size, data[i_part],
                           s_unit * part_data.n_elt[i_part]);
  }

  part_data.s_unit.erase(exch_id);
  part_data.data  .erase(exch_id);

  // free
  free(local_code_name);
  free(cpl_id);
  free(part_data_id);

  svr->state=CWP_SVRSTATE_LISTENINGMSG;
}

/*============================================================================
 * Server function definitions
 *============================================================================*/

/* Create a server */

int
CWP_server_create
(
 MPI_Comm global_comm,
 uint16_t server_port,
 int      flags,
 p_server svr
)
{
  struct sockaddr_in *server_addr = (struct sockaddr_in *) malloc(sizeof(struct sockaddr_in));
  socklen_t max_msg_size_len;
  int vrai=1;

  memset(svr,0,sizeof(t_server));
  svr->port             = server_port;
  svr->flags            = flags;
  svr->server_endianess = CWP_transfer_endian_machine();

  svr_mpi.global_comm  = global_comm;

  // retreive hostname
  if(gethostname(svr->host_name,sizeof(svr->host_name)) != 0) {
    PDM_error(__FILE__, __LINE__, 0, "Could not get host name\n");
    return -1;
  }

  // tmp
  struct addrinfo *svr_info;

  char port_str[16] = {};
  sprintf(port_str, "%d", svr->port);

  const int status = getaddrinfo(svr->host_name, port_str, NULL, &svr_info); // hint
  char *dst = (char *) malloc(sizeof(char) * INET_ADDRSTRLEN);
  inet_ntop(AF_INET, svr_info->ai_addr->sa_data, dst, INET_ADDRSTRLEN);
  free(dst);
  freeaddrinfo(svr_info);

  CWP_UNUSED(status);

  // verbose
  if (svr_debug) {
    printf("CWP:Creating Server on %s port %i...\n", svr->host_name, svr->port);
  }

  // create socket (IPv4, binary data flux, unique protocol for AF_INET + SOCK_STREAM)
  svr->listen_socket = socket(AF_INET, SOCK_STREAM, 0);

  if (svr->listen_socket == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Could not create socket to listen\n");
    return -1;
  }

  // get maximum message size
  max_msg_size_len  = sizeof(int);
  svr->max_msg_size = 0;
  getsockopt(svr->listen_socket, SOL_SOCKET, SO_RCVBUF, (void *) &svr->max_msg_size, &max_msg_size_len);
  svr->max_msg_size = CWP_MSG_MAXMSGSIZE;

  // verbose
  if (svr_debug) {
    printf("CWP:Max message size:%i\n", svr->max_msg_size);
  }

  // ensure restarted process can reconnect to same socket
  if (setsockopt(svr->listen_socket, SOL_SOCKET, SO_REUSEADDR, &vrai, sizeof(int)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Setsockopt failed\n");
    return -1;
  }

  // fill data structure for ip adress + port
  server_addr->sin_family      = AF_INET;
  server_addr->sin_port        = htons(svr->port);
  server_addr->sin_addr.s_addr = INADDR_ANY;
  bzero(&(server_addr->sin_zero),8);

  // bind socket to ip adress + port
  if (bind(svr->listen_socket, (struct sockaddr *) server_addr, sizeof(struct sockaddr)) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Unable to bind socket\n");
    return -1;
  }

  // wait for client to connect
  if (listen(svr->listen_socket, SOMAXCONN) == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Unable to listen on socket\n");
    return -1;
  }

  svr->state = CWP_SVRSTATE_WAITCONN;

  // verbose
  if (svr_debug) {
    printf("CWP:Server created on %s port %i\n", svr->host_name, svr->port);
  }

  // free
  free(server_addr);

  return 0;
}

/* Kill a server */

int
CWP_server_kill
(
 p_server svr
)
{
  // verbose
  if (svr_debug) {
    printf("CWP:Server shutting down\n");
  }

  // shutdown
#ifdef WINDOWS
  shutdown(svr->listen_socket,SD_BOTH);
  shutdown(svr->connected_socket,SD_BOTH);
#else
  shutdown(svr->listen_socket,SHUT_RDWR);
  shutdown(svr->connected_socket,SHUT_RDWR);
#endif

  memset(svr,0,sizeof(t_server));

  return 0;
}

/* Message handler */

int
CWP_server_msg_handler
(
 p_server svr,
 p_message msg
)
{
  switch(msg->message_type) {

  case CWP_MSG_DIE:
    svr->state=CWP_SVRSTATE_TERMINATING;

    // verbose
    if (svr_debug) {
      printf("CWP: server recieved termination signal\n");
    }

    break;

  case CWP_MSG_CWP_INIT:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Init signal\n");
    }

    // launch
    CWP_server_Init(svr);

    break;

  case CWP_MSG_CWP_FINALIZE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Finalize signal\n");
    }

    // launch
    CWP_server_Finalize(svr);

    break;

  case CWP_MSG_CWP_PARAM_LOCK:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_lock signal\n");
    }

    // launch
    CWP_server_Param_lock(svr);

    break;

  case CWP_MSG_CWP_PARAM_UNLOCK:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_unlock signal\n");
    }

    // launch
    CWP_server_Param_unlock(svr);

    break;

  case CWP_MSG_CWP_PARAM_ADD:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_add signal\n");
    }

    // launch
    CWP_server_Param_add(svr);

    break;

  case CWP_MSG_CWP_PARAM_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_get signal\n");
    }

    // launch
    CWP_server_Param_get(svr);

    break;

  case CWP_MSG_CWP_PARAM_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_set signal\n");
    }

    // launch
    CWP_server_Param_set(svr);

    break;

  case CWP_MSG_CWP_PARAM_DEL:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_del signal\n");
    }

    // launch
    CWP_server_Param_del(svr);

    break;

  case CWP_MSG_CWP_PARAM_N_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_n_get signal\n");
    }

    // launch
    CWP_server_Param_n_get(svr);

    break;

  case CWP_MSG_CWP_PARAM_LIST_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_list_get signal\n");
    }

    // launch
    CWP_server_Param_list_get(svr);

    break;

  case CWP_MSG_CWP_PARAM_IS:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_is signal\n");
    }

    // launch
    CWP_server_Param_is(svr);

    break;

  case CWP_MSG_CWP_PARAM_REDUCE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Param_reduce signal\n");
    }

    // launch
    CWP_server_Param_reduce(svr);

    break;

  case CWP_MSG_CWP_CPL_CREATE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Cpl_create signal\n");
    }

    // launch
    CWP_server_Cpl_create(svr);

    break;

  case CWP_MSG_CWP_CPL_BARRIER:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Cpl_barrier signal\n");
    }

    // launch
    CWP_server_Cpl_barrier(svr);

    break;

  case CWP_MSG_CWP_CPL_DEL:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Cpl_del signal\n");
    }

    // launch
    CWP_server_Cpl_del(svr);

    break;

  case CWP_MSG_CWP_PROPERTIES_DUMP:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Properties_dump signal\n");
    }

    // launch
    CWP_server_Properties_dump(svr);

    break;

  case CWP_MSG_CWP_VISU_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Visu_set signal\n");
    }

    // launch
    CWP_server_Visu_set(svr);

    break;

  case CWP_MSG_CWP_STATE_UPDATE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_State_update signal\n");
    }

    // launch
    CWP_server_State_update(svr);

    break;

  case CWP_MSG_CWP_TIME_STEP_BEG:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Time_step_beg signal\n");
    }

    // launch
    CWP_server_Time_step_beg(svr);

    break;

  case CWP_MSG_CWP_TIME_STEP_END:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Time_step_end signal\n");
    }

    // launch
    CWP_server_Time_step_end(svr);

    break;

  case CWP_MSG_CWP_STATE_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_State_get signal\n");
    }

    // launch
    CWP_server_State_get(svr);

    break;

  case CWP_MSG_CWP_CODES_NB_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Codes_nb_get signal\n");
    }

    // launch
    CWP_server_Codes_nb_get(svr);

    break;

  case CWP_MSG_CWP_CODES_LIST_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Codes_list_get signal\n");
    }

    // launch
    CWP_server_Codes_list_get(svr);

    break;

  case CWP_MSG_CWP_LOC_CODES_NB_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Loc_codes_nb_get signal\n");
    }

    // launch
    CWP_server_Loc_codes_nb_get(svr);

    break;

  case CWP_MSG_CWP_LOC_CODES_LIST_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Loc_codes_list_get signal\n");
    }

    // launch
    CWP_server_Loc_codes_list_get(svr);

    break;

  case CWP_MSG_CWP_COMPUTED_TGTS_BCAST_ENABLE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Computed_tgts_bcast_enable signal\n");
    }

    // launch
    CWP_server_Computed_tgts_bcast_enable(svr);

    break;

  case CWP_MSG_CWP_INVOLVED_SRCS_BCAST_ENABLE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Involved_srcs_bcast_enable signal\n");
    }

    // launch
    CWP_server_Involved_srcs_bcast_enable(svr);

    break;

  case CWP_MSG_CWP_N_UNCOMPUTED_TGTS_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_N_uncomputed_tgts_get signal\n");
    }

    // launch
    CWP_server_N_uncomputed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_UNCOMPUTED_TGTS_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Uncomputed_tgts_get signal\n");
    }

    // launch
    CWP_server_Uncomputed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_N_COMPUTED_TGTS_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_N_computed_tgts_get signal\n");
    }

    // launch
    CWP_server_N_computed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_COMPUTED_TGTS_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Computed_tgts_get signal\n");
    }

    // launch
    CWP_server_Computed_tgts_get(svr);

    break;

  case CWP_MSG_CWP_N_INVOLVED_SRCS_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_N_involved_srcs_get\n");
    }

    // launch
    CWP_server_N_involved_srcs_get(svr);

    break;

  case CWP_MSG_CWP_INVOLVED_SRCS_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Involved_srcs_get signal\n");
    }

    // launch
    CWP_server_Involved_srcs_get(svr);

    break;

  case CWP_MSG_CWP_SPATIAL_INTERP_WEIGHTS_COMPUTE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Spatial_interp_weights_compute\n");
    }

    // launch
    CWP_server_Spatial_interp_weights_compute(svr);

    break;

  case CWP_MSG_CWP_SPATIAL_INTERP_PROPERTY_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Spatial_interp_property_set signal\n");
    }

    // launch
    CWP_server_Spatial_interp_property_set(svr);

    break;

  case CWP_MSG_CWP_USER_TGT_PTS_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_User_tgt_pts_set signal\n");
    }

    // launch
    CWP_server_User_tgt_pts_set(svr);

    break;

  case CWP_MSG_CWP_STD_BLOCK_TYPE_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_std_block_type_get signal\n");
    }

    // launch
    CWP_server_std_block_type_get(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_FINALIZE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_finalize signal\n");
    }

    // launch
    CWP_server_Mesh_interf_finalize(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_VTX_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_vtx_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_vtx_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_ADD:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_block_add signal\n");
    }

    // launch
    CWP_server_Mesh_interf_block_add(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_block_std_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_block_std_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_STD_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_block_std_get signal\n");
    }

    // launch
    CWP_server_Mesh_interf_block_std_get(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_f_poly_block_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_f_poly_block_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_F_POLY_BLOCK_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_f_poly_block_get signal\n");
    }

    // launch
    CWP_server_Mesh_interf_f_poly_block_get(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_c_poly_block_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_c_poly_block_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_C_POLY_BLOCK_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_c_poly_block_get signal\n");
    }

    // launch
    CWP_server_Mesh_interf_c_poly_block_get(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_DEL:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_del signal\n");
    }

    // launch
    CWP_server_Mesh_interf_del(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_FROM_CELLFACE_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_from_cellface_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_from_cellface_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_FROM_FACEEDGE_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_from_faceedge_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_from_faceedge_set(svr);

    break;

  case CWP_MSG_CWP_FIELD_CREATE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_create signal\n");
    }

    // launch
    CWP_server_Field_create(svr);

    break;

  case CWP_MSG_CWP_FIELD_DATA_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_data_set signal\n");
    }

    // launch
    CWP_server_Field_data_set(svr);

    break;

  case CWP_MSG_CWP_FIELD_DOF_LOCATION_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_dof_location_get signal\n");
    }

    // launch
    CWP_server_Field_dof_location_get(svr);

    break;

  case CWP_MSG_CWP_FIELD_STORAGE_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_storage_get signal\n");
    }

    // launch
    CWP_server_Field_storage_get(svr);

    break;

  case CWP_MSG_CWP_FIELD_N_DOF_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_n_dof_get signal\n");
    }

    // launch
    CWP_server_Field_n_dof_get(svr);

    break;

  case CWP_MSG_CWP_FIELD_N_COMPONENTS_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_n_components_get signal\n");
    }

    // launch
    CWP_server_Field_n_components_get(svr);

    break;

  case CWP_MSG_CWP_FIELD_DEL:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_del signal\n");
    }

    // launch
    CWP_server_Field_del(svr);

    break;

  case CWP_MSG_CWP_FIELD_ISSEND:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_CWP_Field_issend signal\n");
    }

    // launch
    CWP_server_Field_issend(svr);

    break;

  case CWP_MSG_CWP_FIELD_IRECV:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_irecv signal\n");
    }

    // launch
    CWP_server_Field_irecv(svr);

    break;

  case CWP_MSG_CWP_FIELD_WAIT_ISSEND:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_wait_issend signal\n");
    }

    // launch
    CWP_server_Field_wait_issend(svr);

    break;

  case CWP_MSG_CWP_FIELD_WAIT_IRECV:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Field_wait_irecv signal\n");
    }

    // launch
    CWP_server_Field_wait_irecv(svr);

    break;

  case CWP_MSG_CWP_INTERP_FUNCTION_UNSET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Interp_function_unset signal\n");
    }

    // launch
    CWP_server_Interp_function_unset(svr);

    break;

  case CWP_MSG_CWP_INTERP_FUNCTION_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Interp_function_set signal\n");
    }

    // launch
    CWP_server_Interp_function_set(svr);

    break;

  case CWP_MSG_CWP_OUTPUT_FILE_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Output_file_set signal\n");
    }

    // launch
    CWP_server_Output_file_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_block_ho_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_block_ho_set(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_BLOCK_HO_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_block_ho_get signal\n");
    }

    // launch
    CWP_server_Mesh_interf_block_ho_get(svr);

    break;

  case CWP_MSG_CWP_MESH_INTERF_HO_ORDERING_FROM_IJK_SET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Mesh_interf_ho_ordering_from_IJK_set signal\n");
    }

    // launch
    CWP_server_Mesh_interf_ho_ordering_from_IJK_set(svr);

    break;

  case CWP_MSG_CWP_CPL_SPATIAL_INTERP_ALGO_GET:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Cpl_spatial_interp_algo_get signal\n");
    }

    // launch
    CWP_server_Cpl_spatial_interp_algo_get(svr);

    break;

  case CWP_MSG_CWP_GLOBAL_DATA_ISSEND:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Global_data_issend signal\n");
    }

    // launch
    CWP_server_Global_data_issend(svr);

    break;

  case CWP_MSG_CWP_GLOBAL_DATA_IRECV:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Global_data_irecv signal\n");
    }

    // launch
    CWP_server_Global_data_irecv(svr);

    break;

  case CWP_MSG_CWP_GLOBAL_DATA_WAIT_ISSEND:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Global_data_wait_issend signal\n");
    }

    // launch
    CWP_server_Global_data_wait_issend(svr);

    break;

  case CWP_MSG_CWP_GLOBAL_DATA_WAIT_IRECV:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Global_data_wait_irecv signal\n");
    }

    // launch
    CWP_server_Global_data_wait_irecv(svr);

    break;

  case CWP_MSG_CWP_PART_DATA_CREATE:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Part_data_create signal\n");
    }

    // launch
    CWP_server_Part_data_create(svr);

    break;

  case CWP_MSG_CWP_PART_DATA_DEL:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Part_data_del signal\n");
    }

    // launch
    CWP_server_Part_data_del(svr);

    break;

  case CWP_MSG_CWP_PART_DATA_ISSEND:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Part_data_issend signal\n");
    }

    // launch
    CWP_server_Part_data_issend(svr);

    break;

  case CWP_MSG_CWP_PART_DATA_IRECV:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Part_data_irecv signal\n");
    }

    // launch
    CWP_server_Part_data_irecv(svr);

    break;

  case CWP_MSG_CWP_PART_DATA_WAIT_ISSEND:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Part_data_wait_issend signal\n");
    }

    // launch
    CWP_server_Part_data_wait_issend(svr);

    break;

  case CWP_MSG_CWP_PART_DATA_WAIT_IRECV:

    // verbose
    if (svr_debug) {
      printf("CWP: server received CWP_Part_data_wait_irecv signal\n");
    }

    // launch
    CWP_server_Part_data_wait_irecv(svr);

    break;

  default:
    PDM_error(__FILE__, __LINE__, 0, "Received unknown message type %i, terminating server\n", msg->message_type);
    return -1;
  }

  return 0;
}

/* Run a server */

int
CWP_server_run
(
 p_server svr
)
{
  struct sockaddr_in *client_addr = (struct sockaddr_in *) malloc(sizeof(struct sockaddr_in));
  socklen_t client_addr_len = sizeof(struct sockaddr_in);
  int il_sv_endian;
  int il_cl_endian;

  // accept client connexion
  svr->connected_socket = accept(svr->listen_socket, (struct sockaddr *) client_addr, &client_addr_len);

  // verbose
  if (svr_debug) {
    printf("CWP:got a connection from %s:%d\n",
    inet_ntoa(client_addr->sin_addr),ntohs(client_addr->sin_port));
  }

  // endianess
  svr->state=CWP_SVRSTATE_RECVPPUTDATA;
  if(CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,
           &il_cl_endian,sizeof(int))!=0) {
    PDM_error(__FILE__, __LINE__, 0, "Client endian read failed\n");
    svr->state=CWP_SVRSTATE_LISTENINGMSG;
    return -1;
  }
  il_cl_endian =  ntohl(il_cl_endian);

  svr->client_endianess = il_cl_endian;

  il_sv_endian = htonl(svr->server_endianess);

  if(CWP_transfer_writedata(svr->connected_socket,svr->max_msg_size,
           &il_sv_endian,sizeof(int))!=0) {
    PDM_error(__FILE__, __LINE__, 0, "Server endian send failed\n");
    svr->state=CWP_SVRSTATE_LISTENINGMSG;
    return -1;
  }

  // receive verbose state
  CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size, &svr->flags, sizeof(int));

  svr->state=CWP_SVRSTATE_LISTENINGMSG;

  if (svr_debug) {
    printf("Server : client endian %i server endian %i\n",svr->client_endianess,svr->server_endianess);
  }

  t_message msg;

  while (svr->state != CWP_SVRSTATE_TERMINATING) {

    CWP_transfer_readdata(svr->connected_socket,svr->max_msg_size,&msg, sizeof(t_message));

    if (CWP_server_msg_handler(svr,&msg) != 0) {
      PDM_error(__FILE__, __LINE__, 0, "Server message handling failed\n");
      return -1;
    }

  }

  // free
  free(client_addr);

  return 0;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
