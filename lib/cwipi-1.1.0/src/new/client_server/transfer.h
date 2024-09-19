#ifndef __TRANSFER_H__
#define __TRANSFER_H__
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

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#define CWP_BIG    100
#define CWP_LITTLE 0

/*=============================================================================
 * Public function interfaces
 *============================================================================*/

/* machine endianess */

int
CWP_transfer_endian_machine
(
 void
);

/* read */

int
CWP_transfer_readdata
(
 int socket,
 int batch_size,
 void* dataptr,
 int data_size
);

/* write */

int
CWP_transfer_writedata
(
 int socket,
 int batch_size,
 void* dataptr,
 int data_size
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __TRANSFER_H__ */
