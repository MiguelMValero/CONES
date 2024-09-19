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

#include <bftc_mem.h>
#include "quickSort.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* Fonctions privees */

static int partitionner(int tableau[], int p, int r, int* index) {
  int pivot = tableau[p], i = p-1, j = r+1;
  int temp;
  while(1) {
    do
      j--;
    while(tableau[j] > pivot);
    do
      i++;
    while(tableau[i] < pivot);
    if(i<j) {
      temp = tableau[i];
      tableau[i] = tableau[j];
      tableau[j] = temp;
      if (index != NULL) {
        temp = index[i];
        index[i] = index[j];
        index[j] = temp;
      }
    }
    else
      return j;
  }
  return j;
}

/* Fonctions publiques */

void quickSort(int tableau[], int p, int r, int* index) {
  int q;

  if(p<r) {
    q = partitionner(tableau, p, r, index);
    quickSort(tableau, p, q, index);
    quickSort(tableau, q+1, r, index);
  }
  return;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */

