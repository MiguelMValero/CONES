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
#include <unistd.h>

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <errno.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "transfer.h"

#include <pdm_error.h>
#include <pdm_mpi.h>
#include <pdm_logging.h>

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define SOCKET_TIMEOUT_CNT      60 //total waiting time 60*50ms=3s
#define SOCKET_TIMEOUT_INTERVAL 50 //50ms

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Private function interfaces
 *============================================================================*/

static void milli_sleep(unsigned int milli_seconds) {
  usleep(milli_seconds*1000);
}

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/* machine endianess */

int
CWP_transfer_endian_machine()
{
  int i = 0x12345678;
  if ( (int32_t) htonl(i) == i ) {return CWP_BIG;} else { return CWP_LITTLE;}
}

/* read */

int
CWP_transfer_readdata
(
 int   socket,
 int   batch_size,
 void* data_ptr,
 int   data_size
)
{
  int data_to_xfer=data_size;
  int xfer_size, read_bytes,data_read,timeout;
  unsigned char* buf=(unsigned char*) data_ptr;

  xfer_size = data_to_xfer;
  read_bytes=0;

  data_read=0;
  while( data_to_xfer >0 ) {
    xfer_size=data_to_xfer;
    if (xfer_size>batch_size) {
      xfer_size=batch_size;
    }
    /* poll the socket until data is found if timeout */
    read_bytes=0;
    timeout=SOCKET_TIMEOUT_CNT;
    while(read_bytes<=0 && timeout>0) {
      read_bytes=recv(socket,buf,xfer_size,MSG_PEEK);

      if(read_bytes<=0) {
        milli_sleep(SOCKET_TIMEOUT_INTERVAL);
      }
      timeout--;
    }
    if(timeout<=0 && read_bytes<=0) {
      PDM_error(__FILE__, __LINE__, 0, "Transfer Error, read data read %i/%i bytes and socket returned %i\n",
        data_read,data_size,read_bytes);
      return -1;
    }
    read_bytes=recv(socket,buf,xfer_size,0);
    xfer_size=read_bytes;
    buf+=xfer_size;
    data_read+=xfer_size;
    data_to_xfer-=xfer_size;
  }

  if (read_bytes!=xfer_size) {
    PDM_error(__FILE__, __LINE__, 0, "Transfer Error, read data was waiting for %i bytes and recieved %i bytes\n",
      data_size,data_read);
    return -1;
  }

  return 0;
}

/* write */

int
CWP_transfer_writedata
(
 int   socket,
 int   batch_size,
 void* data_ptr,
 int   data_size
)
{
  int data_to_xfer=data_size;
  int xfer_size, sent_bytes;
  unsigned char* buf=(unsigned char*) data_ptr;
  while( data_to_xfer >0 ) {
    xfer_size=data_to_xfer;
    if (xfer_size>batch_size) {
      xfer_size=batch_size;
    }
    sent_bytes=send(socket,buf,xfer_size,0);
    if(sent_bytes<=0 && xfer_size>0) {
      PDM_error(__FILE__, __LINE__, 0, "Transfer Error, write tried to send %i bytes, socket send returned %i\n",
        xfer_size,sent_bytes);
    return -1;
    }
    xfer_size=sent_bytes;
    buf+=xfer_size;
    data_to_xfer-=xfer_size;
  }

  return 0;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
