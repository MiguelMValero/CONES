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

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm_poly_surf_gen.h"
#include "pdm_config.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"


/*============================================================================
 * Local macro definitions
 *============================================================================*/


#define ABS(a)     ((a) <  0  ? -(a) : (a))
#define MIN(a,b)   ((a) > (b) ?  (b) : (a))

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Private function definitions
 *============================================================================*/


static double random01(void)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / ABS(rsigna - rsignb);
  double resultat = sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}


/*=============================================================================
 * Public function definitions
 *============================================================================*/


/**
 *
 * \brief Generate a distributed polygonal surface mesh
 *
 * \param [in]   pdm_comm         MPI communicator
 * \param [in]   xmin             Minimum x-coordinate
 * \param [in]   xmax             Maximum x-coordinate
 * \param [in]   ymin             Minimum y-coordinate
 * \param [in]   ymax             Maximum y-coordinate
 * \param [in]   have_random      Enable/disable randomization
 * \param [in]   init_random      Random seed
 * \param [in]   nx               Number of vertices in the x-direction
 * \param [in]   ny               Number of vertices in the y-direction
 * \param [out]  ng_face          Global number of faces
 * \param [out]  ng_vtx           Global number of vertices
 * \param [out]  ng_edge          Global number of edges
 * \param [out]  dn_vtx           Local number of vertices
 * \param [out]  dvtx_coord       Coordinates of local vertices (size = 3 * \ref dn_vtx)
 * \param [out]  dn_face          Local number of faces
 * \param [out]  dface_vtx_idx    Index of face-vertex connectivity (size = \ref dn_face + 1)
 * \param [out]  dface_vtx        Distributed face-vertex connectivity (size = \ref dface_vtx_idx[\ref dn_face])
 * \param [out]  dn_edge          Local number of edges
 * \param [out]  dedge_vtx        Distributed edge-vertex connectivity (size = 2 * \ref dn_edge)
 * \param [out]  dedge_face       Distributed edge-face connectivity (size = 2 * \ref dn_edge)
 * \param [out]  n_edge_group     Number of edge groups
 * \param [out]  dedge_group_idx  Index of dedge_group (size = \ref n_edge_group + 1)
 * \param [out]  dedge_group      Distributed lists of edges in each group (size = \ref dedge_group_idx[\ref n_edge_group])
 *
 */

void PDM_poly_surf_gen
(
PDM_MPI_Comm     pdm_comm,
double           xmin,
double           xmax,
double           ymin,
double           ymax,
int              have_random,
int              init_random,
PDM_g_num_t      nx,
PDM_g_num_t      ny,
PDM_g_num_t     *ng_face,
PDM_g_num_t     *ng_vtx,
PDM_g_num_t     *ng_edge,
int             *dn_vtx,
double         **dvtx_coord,
int             *dn_face,
int            **dface_vtx_idx,
PDM_g_num_t    **dface_vtx,
PDM_g_num_t    **dface_edge,
int             *dn_edge,
PDM_g_num_t    **dedge_vtx,
PDM_g_num_t    **dedge_face,
int             *n_edge_group,
int            **dedge_group_idx,
PDM_g_num_t    **dedge_group
)
{
  int verbose = 0;

  int n_rank;
  PDM_MPI_Comm_size(pdm_comm, &n_rank);

  int local_rank;
  PDM_MPI_Comm_rank(pdm_comm, &local_rank);

  const double coef_rand = 0.3;
  srand(init_random);

  /* Comptage des entites du maillage */
  /* -------------------------------- */

  /* nx et ny doivent etre pair et superieur a 4 */

  PDM_g_num_t nx1 = nx;
  PDM_g_num_t ny1 = ny;

  if (nx1 < 6)
    nx1 = 6;

  if (ny1 < 6)
    ny1 = 6;

  if (nx1 % 2 != 0)
    nx1 += 1;

  if (ny1 % 2 != 0)
    ny1 += 1;

  /* cotes d un octogone */

  double cote1;
  if (nx1 == 4)
    cote1 = (xmax - xmin) / (sqrt(2) + 1);
  else
    cote1 = (xmax - xmin) / ((nx1 - 2)/2 + ((nx1 - 2)/2 - 1) * sqrt(2) + sqrt(2));
  double cote2;
  if (ny1 == 4)
    cote2 = (ymax - ymin) / (sqrt(2) + 1);
  else
    cote2 = (ymax - ymin) / ((ny1 - 2)/2 + ((ny1 - 2)/2 - 1) * sqrt(2) + sqrt(2));

  /* Comptage des Sommets */

  PDM_g_num_t dn_vtxTotal = 0;

  PDM_g_num_t cpt = 0;
  PDM_g_num_t cpt_max;

  if (ny1 > 4)
    cpt_max = ny1 + (ny1-4)/2;
  else
    cpt_max = ny1;

  while(cpt < cpt_max) {
    if (cpt % 3 == 1 || cpt % 3 == 2) {
      dn_vtxTotal +=  nx1/2;
    }
    else if ((cpt % 3) == 0) {
      if ((cpt == (cpt_max-1)) || (cpt == 0)) {
        dn_vtxTotal++;
      }

      dn_vtxTotal++;

      for (PDM_g_num_t ix = 2; ix < nx1-1; ix++) {
        dn_vtxTotal++;
      }
      if ((cpt == (cpt_max-1)) || (cpt == 0)) {
        dn_vtxTotal++;
      }
    }
    cpt++;
  }

  /* Nombre de limites */

  *n_edge_group = 4;
  *dedge_group_idx = (int *) malloc(sizeof(int) * (*n_edge_group + 1));
  (*dedge_group_idx)[0] = 0;

  PDM_g_num_t *dn_edgeLim13Rank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank+1));
  dn_edgeLim13Rank[0] = 0;
  for (int i = 0; i < n_rank; i++)
    dn_edgeLim13Rank[i+1] = (PDM_g_num_t) (nx1 - 1) / n_rank;

  PDM_g_num_t _reste = (nx1 - 1) % n_rank;
  int reste =  (int) _reste;

  for (int i = 0; i < n_rank; i++) {
    if (i < reste) {
      dn_edgeLim13Rank[i+1] += 1;
    }
  }


  PDM_g_num_t dn_edgeLim13 = dn_edgeLim13Rank[local_rank + 1];

  for (int i = 0; i < n_rank; i++) {
    dn_edgeLim13Rank[i+1] = dn_edgeLim13Rank[i+1] + dn_edgeLim13Rank[i];
  }

  PDM_g_num_t *dn_edgeLim24Rank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank+1));
  dn_edgeLim24Rank[0] = 0;
  for (int i = 0; i < n_rank; i++)
    dn_edgeLim24Rank[i+1] = (PDM_g_num_t) (ny1 - 1) / n_rank;

  _reste = (ny1 - 1) % n_rank;
  reste =  (int) (_reste);

  for (int i = 0; i < n_rank; i++) {
    if (i < reste) {
      dn_edgeLim24Rank[i+1] += 1;
    }
  }

  PDM_g_num_t dn_edgeLim24 = dn_edgeLim24Rank[local_rank + 1];

  for (int i = 0; i < n_rank; i++) {
    dn_edgeLim24Rank[i+1] = dn_edgeLim24Rank[i+1] + dn_edgeLim24Rank[i];
  }

  *dedge_group = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (2*dn_edgeLim13 + 2*dn_edgeLim24));

  (*dedge_group_idx)[1] = (int) dn_edgeLim13;
  (*dedge_group_idx)[2] = (int) dn_edgeLim24;
  (*dedge_group_idx)[3] = (int) dn_edgeLim13;
  (*dedge_group_idx)[4] = (int) dn_edgeLim24;

  for (int i= 0; i < 4; i++)
    (*dedge_group_idx)[i+1] = (*dedge_group_idx)[i] + (*dedge_group_idx)[i+1];


  /* Construction des sommets */
  /* ------------------------ */

  PDM_g_num_t *dn_vtx_rank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank+1));

  dn_vtx_rank[0] = 0;
  for (int i = 0; i < n_rank; i++)
    dn_vtx_rank[i+1] = (PDM_g_num_t) dn_vtxTotal / n_rank;

  _reste = dn_vtxTotal % n_rank;
  reste =  (int) (_reste);

  for (int i = 0; i < n_rank; i++) {
    if (i < reste) {
      dn_vtx_rank[i+1] += 1;
    }
  }

  *dn_vtx = (int) dn_vtx_rank[local_rank + 1];

  for (int i = 0; i < n_rank; i++) {
    dn_vtx_rank[i+1] = dn_vtx_rank[i+1] + dn_vtx_rank[i];
  }

  *dvtx_coord = (double*) malloc (sizeof(double) * 3 * (*dn_vtx));

  /* Construction des elements et des aretes */
  /* --------------------------------------- */

  /* Comptage des elements et des aretes */

  PDM_g_num_t  dn_face_total = 0;

  /* Triangles */

  /* -- Premiere ligne */
  dn_face_total += nx1/2;

  /* -- Autres triangles (un a gauche un a droite */
  PDM_g_num_t nbLi = (ny1-4)/2;
  dn_face_total += 2*nbLi;

  /* -- Derniere ligne */
  dn_face_total += nx1/2;

  /* Quadrangles */
  PDM_g_num_t nx_quad = (nx1-4)/2;
  PDM_g_num_t ny_quad = (ny1-4)/2;
  dn_face_total += ny_quad * nx_quad;

  /* Polygones */
  PDM_g_num_t nx_poly = (nx1-2)/2;
  PDM_g_num_t ny_poly = (ny1-2)/2;
  dn_face_total += nx_poly * ny_poly;

  PDM_g_num_t dn_edgeTotal = (6 * nx_poly + 1) * ny_poly + nx_poly; /* Aretes touchant un polygone */
  dn_edgeTotal += nx1 + ny1; /* Aretes des cotes ne touchant pa sun polygone */

  PDM_g_num_t n_poly = nx_poly * ny_poly;
  PDM_g_num_t n_quad = nx_quad * ny_quad;
  PDM_g_num_t n_tri  = dn_face_total - n_quad - n_poly;

  /* Allocation des elts */

  PDM_g_num_t *dn_face_rank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank+1));
  PDM_g_num_t *dn_edge_rank = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (n_rank+1));

  dn_face_rank[0] = 0;
  for (int i = 0; i < n_rank; i++)
    dn_face_rank[i+1] = dn_face_total / n_rank ;

  _reste = dn_face_total % n_rank;
  reste =  (int) _reste ;
  for (int i = 0; i < n_rank; i++) {
    if (i < reste) {
      dn_face_rank[i+1] += 1;
    }
  }

  dn_edge_rank[0] = 0;
  for (int i = 0; i < n_rank; i++)
    dn_edge_rank[i+1] = dn_edgeTotal / n_rank ;

  _reste = dn_edgeTotal % n_rank;
  reste =  (int) _reste;
  for (int i = 0; i < n_rank; i++) {
    if (i < reste) {
      dn_edge_rank[i+1] += 1;
    }
  }

  *dn_face = (int) dn_face_rank[local_rank+1];
  *dn_edge = (int) dn_edge_rank[local_rank+1];

  for (int i = 0; i < n_rank; i++) {
    dn_face_rank[i+1] = dn_face_rank[i+1] + dn_face_rank[i];
    dn_edge_rank[i+1] = dn_edge_rank[i+1] + dn_edge_rank[i];
  }

  *dface_vtx_idx = (int *) malloc(sizeof(int) * ((*dn_face)+1));
  *dface_vtx     = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 8 * (*dn_face));
  *dface_edge    = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 8 * (*dn_face));
  *dedge_vtx     = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 2 * (*dn_edge));
  *dedge_face    = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * 2 * (*dn_edge));
  for (int i=0; i<(*dn_face)+1; i++)
    (*dface_vtx_idx)[i]=-1;
  for (int i=0; i<(8 * (*dn_face)); i++)
    (*dface_vtx)[i]=-1;
  for (int i=0; i<(8 * (*dn_face)); i++)
    (*dface_edge)[i]=-1;
  for (int i=0; i<(2 * (*dn_edge)); i++)
    (*dedge_vtx)[i]=-1;
  for (int i=0; i<(2 * (*dn_edge)); i++)
    (*dedge_face)[i]=-1;

  /* Calcul des entites */
  /* ------------------ */

  /* Calcul des coordonnees */

  PDM_g_num_t dn_vtx_tmp = 0;

  cpt = 0;
  if (ny1 > 4)
    cpt_max = ny1 + (ny1-4)/2;
  else
    cpt_max = ny1;

  double ycourant = ymin;
  const double eps = 1e-5;

  while (cpt < cpt_max) {
    if (cpt % 3 == 1 || cpt % 3 == 2) {
      for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
        if ((dn_vtx_rank[local_rank] <= dn_vtx_tmp) && (dn_vtx_rank[local_rank+1] > dn_vtx_tmp)) {
          PDM_g_num_t _local_idx = dn_vtx_tmp - dn_vtx_rank[local_rank];
          int local_idx = (int) _local_idx;
          (*dvtx_coord)[3*local_idx]   = xmin + ix * (1+sqrt(2)) * cote1;
          (*dvtx_coord)[3*local_idx+1] = ycourant;
          (*dvtx_coord)[3*local_idx+2] = 0.;
        }
        dn_vtx_tmp++;
      }
    }
    else if ((cpt % 3) == 0) {
      if ((cpt == (cpt_max-1)) || (cpt == 0)) {
        if ((dn_vtx_rank[local_rank] <= dn_vtx_tmp) && (dn_vtx_rank[local_rank+1] > dn_vtx_tmp)) {
          PDM_g_num_t _local_idx = dn_vtx_tmp - dn_vtx_rank[local_rank];
          int local_idx = (int) _local_idx;
          (*dvtx_coord)[3*local_idx]   = xmin;
          (*dvtx_coord)[3*local_idx+1] = ycourant;
          (*dvtx_coord)[3*local_idx+2] = 0.;
        }
        dn_vtx_tmp++;
      }

      if ((dn_vtx_rank[local_rank] <= dn_vtx_tmp) && (dn_vtx_rank[local_rank+1] > dn_vtx_tmp)) {
        PDM_g_num_t _local_idx = dn_vtx_tmp - dn_vtx_rank[local_rank];
        int local_idx = (int) _local_idx;
        (*dvtx_coord)[3*local_idx]   = xmin + cote1/sqrt(2);
        (*dvtx_coord)[3*local_idx+1] = ycourant;
        (*dvtx_coord)[3*local_idx+2] = 0.;
      }
      dn_vtx_tmp++;

      double xbase = xmin + cote1/sqrt(2);

      PDM_g_num_t nx11 = (nx1-2)/2;
      for (PDM_g_num_t ix = 0; ix < nx11; ix++) {
        if ((dn_vtx_rank[local_rank] <= dn_vtx_tmp) && (dn_vtx_rank[local_rank+1] > dn_vtx_tmp)) {
          PDM_g_num_t _local_idx = dn_vtx_tmp - dn_vtx_rank[local_rank];
          int local_idx = (int) _local_idx;
          (*dvtx_coord)[3*local_idx] = xbase + (ix+1) * cote1 + ix * cote1*sqrt(2);
          (*dvtx_coord)[3*local_idx+1] = ycourant;
          (*dvtx_coord)[3*local_idx+2] = 0.;
        }
        dn_vtx_tmp++;
        if (ix < (nx11 - 1)) {
          if ((dn_vtx_rank[local_rank] <= dn_vtx_tmp) && (dn_vtx_rank[local_rank+1] > dn_vtx_tmp)) {
            PDM_g_num_t _local_idx = dn_vtx_tmp - dn_vtx_rank[local_rank];
            int local_idx = (int) _local_idx;
            (*dvtx_coord)[3*local_idx] = xbase + (ix+1) * cote1 + (ix+1) * cote1*sqrt(2);
            (*dvtx_coord)[3*local_idx+1] = ycourant;
            (*dvtx_coord)[3*local_idx+2] = 0.;
          }
          dn_vtx_tmp++;
        }
      }

      if ((cpt == (cpt_max-1)) || (cpt == 0)) {
        if ((dn_vtx_rank[local_rank] <= dn_vtx_tmp) && (dn_vtx_rank[local_rank+1] > dn_vtx_tmp)) {
          PDM_g_num_t _local_idx = dn_vtx_tmp - dn_vtx_rank[local_rank];
          int local_idx = (int) _local_idx;
          (*dvtx_coord)[3*local_idx]   = xmax;
          (*dvtx_coord)[3*local_idx+1] = ycourant;
          (*dvtx_coord)[3*local_idx+2] = 0.;
        }
        dn_vtx_tmp++;
      }
    }
    cpt++;
    if ((cpt % 3 == 1) || (cpt % 3 == 0))
      ycourant += cote2/sqrt(2);
    else
      ycourant += cote2;
  }

  *ng_vtx = dn_vtx_tmp;

  /* Perturbation des coordonnees + Creation d'une courbure en Z */

  for (PDM_g_num_t ix = 0; ix <(*dn_vtx) ; ix++) {
    if (ABS(xmin-(*dvtx_coord)[3*ix]) > eps &&
        ABS(xmax-(*dvtx_coord)[3*ix]) > eps &&
        ABS(ymax-(*dvtx_coord)[3*ix+1]) > eps &&
        ABS(ymin-(*dvtx_coord)[3*ix+1]) > eps) {
      if (have_random != 0) {
        (*dvtx_coord)[3*ix]   += random01() * coef_rand * cote1;
        (*dvtx_coord)[3*ix+1] += random01() * coef_rand * cote2;
      }
    }
    //dvtx_coord[3*ix+2]=2*sin(3*dvtx_coord[3*ix+1])*sin(3*dvtx_coord[3*ix]);
    //dvtx_coord[3*ix+2]=20*sin(dvtx_coord[3*ix+1]/5.)*sin(dvtx_coord[3*ix]/5.);
    //(*dvtx_coord)[3*ix+2] = random01() * 1e-9 * PDM_MIN (cote1, cote2); // Perturbation 0 machine plan
    //printf("%12.5e\n", (*dvtx_coord)[3*ix+2]);

  }

  /* Construction simultanée des connectivites des elements et des aretes internes */
  /* ----------------------------------------------------------------------------- */

  int dn_face_tmp = 0;
  int dn_edge_tmp = 0;
  int dn_face_abs = 0;
  int dn_edge_abs = 0;

  PDM_g_num_t n1;
  PDM_g_num_t n2;
  PDM_g_num_t n3;
  PDM_g_num_t n4;
  PDM_g_num_t n5;
  PDM_g_num_t n6;
  PDM_g_num_t n7;
  PDM_g_num_t n8;

  /* Triangle */

  /* -- Premiere ligne */

  n1 = 1;
  n2 = n1+nx1;
  (*dface_vtx_idx)[0] = 0;

  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    PDM_g_num_t ix1 = 2 * ix;

    if ((dn_face_rank[local_rank] <= dn_face_abs) && (dn_face_rank[local_rank+1] > dn_face_abs)) {
      int ideb = (*dface_vtx_idx)[dn_face_tmp] ;
      (*dface_vtx)[ideb]   = n1+ix1;
      (*dface_vtx)[ideb+1] = n1+ix1+1;
      (*dface_vtx)[ideb+2] = n2+ix;

      if (ix == 0)
        (*dface_edge)[ideb]   = 1;
      else
        (*dface_edge)[ideb]   = n_tri + 4 + 6*(ix-1) + 2;

      if (ix == 0)
        (*dface_edge)[ideb+1]   = dn_edge_abs + 2;
      else
        (*dface_edge)[ideb+1]   = dn_edge_abs + 1;

      if (ix == (nx1/2 - 1))
        (*dface_edge)[ideb+2]   = dn_edge_abs + 2;
      else if (ix == ((nx1/2 - 1) - 1))
        (*dface_edge)[ideb+2]   = n_tri + 4 + 6*ix + 7;
      else
        (*dface_edge)[ideb+2]   = n_tri + 4 + 6*ix + 6;

      dn_face_tmp += 1;
      (*dface_vtx_idx)[dn_face_tmp] = ideb + 3;
    }
    dn_face_abs += 1;

    if (ix == 0) {
      if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
        (*dedge_vtx)[2*dn_edge_tmp]     = n2+ix;
        (*dedge_vtx)[2*dn_edge_tmp + 1] = n1+ix1;
        (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
        (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
        dn_edge_tmp += 1;
      }
      dn_edge_abs += 1;
    }

    if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
      (*dedge_vtx)[2*dn_edge_tmp]     = n1+ix1;
      (*dedge_vtx)[2*dn_edge_tmp + 1] = n1+ix1+1;
      (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
      (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
      dn_edge_tmp += 1;
    }
    dn_edge_abs += 1;

    if (ix == (nx1/2 - 1)) {
      if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
        (*dedge_vtx)[2*dn_edge_tmp]     = n1+ix1+1;
        (*dedge_vtx)[2*dn_edge_tmp + 1] = n2+ix;
        (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
        (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
        dn_edge_tmp += 1;
      }
      dn_edge_abs += 1;
    }

  }

  /* -- Autres triangles (un a gauche un a droite */

  n1 = 1 + nx1 + nx1/2;
  n2 = n1 + nx1/2;
  n3 = n2 + nx1-2;

  for (PDM_g_num_t itri = 0; itri < nbLi; itri++) {
    n4 = n1 + nx1/2 - 1;
    n5 = n2 + nx1-2 - 1;
    n6 = n3 + nx1/2 - 1;

    if ((dn_face_rank[local_rank] <= dn_face_abs) && (dn_face_rank[local_rank+1] > dn_face_abs)) {
      int ideb = (*dface_vtx_idx)[dn_face_tmp];
      (*dface_vtx)[ideb]   = n1;
      (*dface_vtx)[ideb+1] = n2;
      (*dface_vtx)[ideb+2] = n3;
      (*dface_edge)[ideb]     = dn_edge_abs + 1;
      (*dface_edge)[ideb+1]   = n_tri + 4 + (6 * nx_poly + 1) * itri + 4;
      if (itri == nbLi - 1)
        (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (itri+1) + 7;
      else
        (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (itri+1) + 6;
      dn_face_tmp += 1;
      (*dface_vtx_idx)[dn_face_tmp] = ideb + 3;

    }
    dn_face_abs += 1;

    if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
      (*dedge_vtx)[2*dn_edge_tmp]     = n3;
      (*dedge_vtx)[2*dn_edge_tmp + 1] = n1;
      (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
      (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
      dn_edge_tmp += 1;
    }
    dn_edge_abs += 1;

    if ((dn_face_rank[local_rank] <= dn_face_abs) && (dn_face_rank[local_rank+1] > dn_face_abs)) {
      int ideb = (*dface_vtx_idx)[dn_face_tmp];
      (*dface_vtx)[ideb]   = n4;
      (*dface_vtx)[ideb+1] = n6;
      (*dface_vtx)[ideb+2] = n5;
      (*dface_edge)[ideb]     = dn_edge_abs + 1;
      if (itri == nbLi - 1){
        (*dface_edge)[ideb+1]   = n_tri + 4 + (6 * nx_poly + 1) * (itri + 1) + (7 * nx_poly + 1) - 6;
      }
      else{
        (*dface_edge)[ideb+1]   = n_tri + 4 + (6 * nx_poly + 1) * (itri+2) - 5;
      }
      (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (itri + 1) - 3;
      dn_face_tmp += 1;
      (*dface_vtx_idx)[dn_face_tmp] = ideb + 3;
    }
    dn_face_abs += 1;

    if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
      (*dedge_vtx)[2*dn_edge_tmp]     = n4;
      (*dedge_vtx)[2*dn_edge_tmp + 1] = n6;
      (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
      (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
      dn_edge_tmp += 1;
    }
    dn_edge_abs += 1;

    n1 = n3 + nx1/2;
    n2 = n1 + nx1/2;
    n3 = n2 + nx1-2;
  }

  /* -- Derniere ligne */

  n2 = n1 + nx1/2;
  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    PDM_g_num_t ix1 = 2 * ix;
    if ((dn_face_rank[local_rank] <= dn_face_abs) && (dn_face_rank[local_rank+1] > dn_face_abs)) {
      int ideb = (*dface_vtx_idx)[dn_face_tmp] ;
      (*dface_vtx)[ideb]   = n1 + ix;
      (*dface_vtx)[ideb+1] = n2 + ix1 + 1;
      (*dface_vtx)[ideb+2] = n2 + ix1;

      if (ix == 0)
        (*dface_edge)[ideb]   = dn_edge_abs + 1;
      else if (ix == (nx1/2 - 1))
        (*dface_edge)[ideb]   = n_tri + 4 + (6 * nx_poly + 1) * (ny_poly - 1) + 7*(ix-1) + 4;
      else
        (*dface_edge)[ideb]   = n_tri + 4 + (6 * nx_poly + 1) * (ny_poly - 1) + 7*(ix-1) + 3;

      if (ix == 0)
        (*dface_edge)[ideb+1]   = dn_edge_abs + 2;
      else
        (*dface_edge)[ideb+1]   = dn_edge_abs + 1;

      if (ix == (nx1/2 - 1))
        (*dface_edge)[ideb+2]   = dn_edge_abs + 2;
      else if (ix == ((nx1/2 - 1) - 1))
        (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (ny_poly - 1) + 7*(nx_poly - 1) + 6;
      else
        (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (ny_poly - 1) + 7*ix + 5;
      dn_face_tmp += 1;
      (*dface_vtx_idx)[dn_face_tmp] = ideb + 3;
    }
    dn_face_abs += 1;

    if (ix == 0) {
      if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
        (*dedge_vtx)[2*dn_edge_tmp]     = n2 + ix1;
        (*dedge_vtx)[2*dn_edge_tmp + 1] = n1 + ix;
        (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
        (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
        dn_edge_tmp += 1;
      }
      dn_edge_abs += 1;
    }

    if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
      (*dedge_vtx)[2*dn_edge_tmp]     = n2 + ix1;
      (*dedge_vtx)[2*dn_edge_tmp + 1] = n2 + ix1 + 1;
      (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
      (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
      dn_edge_tmp += 1;
    }
    dn_edge_abs += 1;

    if (ix == (nx1/2 - 1)) {
      if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
        (*dedge_vtx)[2*dn_edge_tmp]     = n1 + ix;
        (*dedge_vtx)[2*dn_edge_tmp + 1] = n2 + ix1 + 1;
        (*dedge_face)[2*dn_edge_tmp]       = dn_face_abs;
        (*dedge_face)[2*dn_edge_tmp + 1]   = 0;
        dn_edge_tmp += 1;
      }
      dn_edge_abs += 1;
    }
  }

  /* Quadrangle */

  for (PDM_g_num_t iy = 0; iy < ny_quad; iy++) {
    for (PDM_g_num_t ix = 0; ix < nx_quad; ix++) {
      n1 = iy*(2*nx1-2) + nx1 + nx1/2 + 1 + ix + 1;
      n2 = iy*(2*nx1-2) + 2*nx1 + 1 + 2*ix + 1;
      n3 = n2 + 1;
      n4 = iy*(2*nx1-2) + 3*nx1 - 2 + 1 + ix + 1;
      if ((dn_face_rank[local_rank] <= dn_face_abs) && (dn_face_rank[local_rank+1] > dn_face_abs)) {
        int ideb = (*dface_vtx_idx)[dn_face_tmp];
        (*dface_vtx)[ideb]   = n1;
        (*dface_vtx)[ideb+1] = n3;
        (*dface_vtx)[ideb+2] = n4;
        (*dface_vtx)[ideb+3] = n2;
        (*dface_edge)[ideb]     = n_tri + 4 + (6 * nx_poly + 1) * iy     + 6*ix     + 3;

        if (ix == nx_quad - 1)
          (*dface_edge)[ideb+1]   = n_tri + 4 + (6 * nx_poly + 1) * iy     + 6*(ix+1) + 5;
        else
          (*dface_edge)[ideb+1]   = n_tri + 4 + (6 * nx_poly + 1) * iy     + 6*(ix+1) + 4;

        if (iy == ny_quad - 1)
          if (ix == nx_quad - 1)
            (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (iy+1) + 7*(ix+1) + 8;
          else
            (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (iy+1) + 7*(ix+1) + 7;
        else
          if (ix == nx_quad - 1)
            (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (iy+1) + 6*(ix+1) + 7;
          else
            (*dface_edge)[ideb+2]   = n_tri + 4 + (6 * nx_poly + 1) * (iy+1) + 6*(ix+1) + 6;

        if (iy == ny_quad - 1)
          (*dface_edge)[ideb+3]   = n_tri + 4 + (6 * nx_poly + 1) * (iy+1) + 7*ix     + 2;
        else
          (*dface_edge)[ideb+3]   = n_tri + 4 + (6 * nx_poly + 1) * (iy+1) + 6*ix     + 2;
        dn_face_tmp += 1;
        (*dface_vtx_idx)[dn_face_tmp] = ideb + 4;
      }
      dn_face_abs += 1;
    }
  }

  /* Polygones */

  PDM_g_num_t delta = 0;
  PDM_g_num_t ipoLy = 0;

  for (PDM_g_num_t iy = 0; iy < ny_poly; iy++) {
    if (iy == 1)
      delta += 2*nx1;
    else if (iy != 0)
      delta += 2*nx1-2;
    for (PDM_g_num_t ix = 0; ix < nx_poly; ix++) {

      ipoLy += 1;
      if (iy == 0)
        n1 = delta + 1 + 2*ix +1;
      else
        n1 = delta + 1 + 2*ix ;
      n2 = n1 + 1;
      n3 = iy*(2*nx1-2) + 1 + nx1 + ix;
      n4 = n3 + 1;
      n5 = iy*(2*nx1-2) + 1 + nx1 + nx1/2 + ix;
      n6 = n5 + 1;
      n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix;
      if (iy == (ny_poly - 1))
        n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix + 1;
      n8 = n7 + 1;

      PDM_g_num_t connecPoly[8];
      connecPoly[0]      = n1;
      connecPoly[1]      = n2;
      connecPoly[2]      = n4;
      connecPoly[3]      = n6;
      connecPoly[4]      = n8;
      connecPoly[5]      = n7;
      connecPoly[6]      = n5;
      connecPoly[7]      = n3;

      if ((dn_face_rank[local_rank] <= dn_face_abs) && (dn_face_rank[local_rank+1] > dn_face_abs)) {
        int ideb = (*dface_vtx_idx)[dn_face_tmp] ;
        for (int k = 0; k < 8; k++)
          (*dface_vtx)[ideb+k] = connecPoly[k];
        int id = 0;
        (*dface_edge)[ideb]     = dn_edge_abs + (++id);
        (*dface_edge)[ideb+1]   = dn_edge_abs + (++id);
        if (ix == (nx_poly - 1))
          (*dface_edge)[ideb+2]   = dn_edge_abs + (++id);
        else {
          (*dface_edge)[ideb+2]   = dn_edge_abs+(id+1) + 8;
          if (ix == (nx_poly-2))
            (*dface_edge)[ideb+2]= (*dface_edge)[ideb+2]+1;
          if (iy == (ny_poly-1))
            (*dface_edge)[ideb+2]= (*dface_edge)[ideb+2]+2;
        }
        (*dface_edge)[ideb+3]   = dn_edge_abs + (++id);
        if (iy == (ny_poly - 1))
          (*dface_edge)[ideb+4]   = dn_edge_abs + (++id);
        else {
          (*dface_edge)[ideb+4]   = dn_edge_abs+(id+1) + (6*(nx_poly-1)+1) +3;
          if (ix == (nx_poly-1))
            (*dface_edge)[ideb+4]= (*dface_edge)[ideb+4]-1;
          if (iy == (ny_poly-2))
            (*dface_edge)[ideb+4]= (*dface_edge)[ideb+4]+ix;
        }
        (*dface_edge)[ideb+5]   = dn_edge_abs + (++id);
        (*dface_edge)[ideb+6]   = dn_edge_abs + (++id);
        (*dface_edge)[ideb+7]   = dn_edge_abs + (++id);
        dn_face_tmp += 1;
        (*dface_vtx_idx)[dn_face_tmp] = ideb + 8;
      }
      dn_face_abs += 1;

      /* Definition de toutes les aretes internes */

      for (int k = 0; k < 8; k++) {
        if (!((k == 2) && (ix != (nx_poly - 1))) &&
            !((k == 4) && (iy != (ny_poly - 1)))) {


          if ((dn_edge_rank[local_rank] <= dn_edge_abs) && (dn_edge_rank[local_rank+1] > dn_edge_abs)) {
            (*dedge_vtx)[2*dn_edge_tmp]     = connecPoly[k];
            (*dedge_vtx)[2*dn_edge_tmp + 1] = connecPoly[(k + 1)%8];
            (*dedge_face)[2*dn_edge_tmp]       = n_tri + n_quad + ipoLy;
            if (k == 0) {
              if (iy == 0)
                (*dedge_face)[2*dn_edge_tmp + 1] = 0;
              else
                (*dedge_face)[2*dn_edge_tmp + 1] = n_tri + n_quad + ipoLy - nx_poly;
            }
            else if (k == 1) {
              if (iy == 0)
                (*dedge_face)[2*dn_edge_tmp + 1] = ix + 2;
              else {
                if (ix == (nx_poly - 1))
                  (*dedge_face)[2*dn_edge_tmp + 1] = nx1/2 + 2*(iy-1) + 2;
                else
                  (*dedge_face)[2*dn_edge_tmp + 1] = n_tri + (iy - 1)*nx_quad + ix + 1;
              }
            }
            else if (k == 2) {
              (*dedge_face)[2*dn_edge_tmp + 1] = 0;
            }
            else if (k == 3) {
              if (iy == (ny_poly - 1))
                (*dedge_face)[2*dn_edge_tmp + 1] = nx1/2 + 2*nbLi + ix +2;
              else {
                if (ix == (nx_poly - 1))
                  (*dedge_face)[2*dn_edge_tmp + 1] = nx1/2 + 2*(iy+1);
                else
                  (*dedge_face)[2*dn_edge_tmp + 1] = n_tri + iy*nx_quad + ix + 1;
              }
            }
            else if (k == 4) {
              (*dedge_face)[2*dn_edge_tmp + 1] = 0;
            }
            else if (k == 5) {
              if (iy == (ny_poly - 1))
                (*dedge_face)[2*dn_edge_tmp + 1] = nx1/2 + 2*nbLi + ix+1;
              else {
                if (ix == 0)
                  (*dedge_face)[2*dn_edge_tmp + 1] = nx1/2 + 2*(iy+1) - 1;
                else
                  (*dedge_face)[2*dn_edge_tmp + 1] = n_tri + iy*nx_quad + ix + 1 - 1;
              }
            }
            else if (k == 6) {
              if (ix == 0)
                (*dedge_face)[2*dn_edge_tmp + 1] = 0;
              else
                (*dedge_face)[2*dn_edge_tmp + 1] = n_tri + n_quad + ipoLy - 1;
            }
            else if (k == 7) {
              if (iy == 0)
                (*dedge_face)[2*dn_edge_tmp + 1] = ix + 1;
              else {
                if (ix == 0)
                  (*dedge_face)[2*dn_edge_tmp + 1] = nx1/2 + 2*(iy-1) + 1;
                else
                  (*dedge_face)[2*dn_edge_tmp + 1] = n_tri + (iy-1)*nx_quad + ix + 1 - 1;
              }
            }
            dn_edge_tmp += 1;
          }
          dn_edge_abs += 1;
        }
      }
    }
  }

  *ng_edge = dn_edge_abs;
  *ng_face = dn_face_abs;

  if (verbose && local_rank == 0) printf("gn_vtx = "PDM_FMT_G_NUM"\ngn_elt = "PDM_FMT_G_NUM"\n", *ng_vtx, *ng_face);

  /* Definition des limites */
  /* ---------------------- */

  int idxDedge_group = 0;

  /* - lim 1 - (triangles ligne bas + octogones bas) */

  PDM_g_num_t dn_edgeGroupAbs = 0;
  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    if ((dn_edgeLim13Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim13Rank[local_rank+1] > dn_edgeGroupAbs)) {
      (*dedge_group)[idxDedge_group++]  = 2 + ix;
    }
    ++dn_edgeGroupAbs;
  }

  for (PDM_g_num_t ix = 0; ix < nx_poly; ix++) {
    if ((dn_edgeLim13Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim13Rank[local_rank+1] > dn_edgeGroupAbs)) {
      (*dedge_group)[idxDedge_group++]  = nx1 + ny1 + (6 * ix) + 1;
    }
    ++dn_edgeGroupAbs;
  }

  /* - lim 2 - */

  dn_edgeGroupAbs = 0;
  for (PDM_g_num_t iy = 0; iy < ny1/2; iy++) {
    if ((dn_edgeLim24Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim24Rank[local_rank+1] > dn_edgeGroupAbs)) {
      if (iy == (ny1/2 - 1))
        (*dedge_group)[idxDedge_group++]  = nx1 + ny1;
      else
        (*dedge_group)[idxDedge_group++]  = nx1/2 + 2*iy + 2;
    }
    ++dn_edgeGroupAbs;
  }

  for (PDM_g_num_t iy = 0; iy < ny_poly; iy++) {
    if ((dn_edgeLim24Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim24Rank[local_rank+1] > dn_edgeGroupAbs)) {
      if (iy == ny_poly - 1)
        (*dedge_group)[idxDedge_group++]  =  nx1 + ny1 + ((6 * nx_poly) + 1) * (ny_poly - 1) + 7 * (nx_poly - 1) + 3;
      else
        (*dedge_group)[idxDedge_group++]  =  nx1 + ny1 + ((6 * nx_poly) + 1) * (iy+1) - 4;
    }
    ++dn_edgeGroupAbs;
  }

  /* - lim 3 - */

  dn_edgeGroupAbs = 0;
  for (PDM_g_num_t ix = 0; ix < nx1/2; ix++) {
    if ((dn_edgeLim13Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim13Rank[local_rank+1] > dn_edgeGroupAbs)) {
      (*dedge_group)[idxDedge_group++]  = nx1 + ny1 - nx1/2 + ix;
    }
    ++dn_edgeGroupAbs;
  }

  for (PDM_g_num_t ix = 0; ix < nx_poly; ix++) {
    if ((dn_edgeLim13Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim13Rank[local_rank+1] > dn_edgeGroupAbs)) {
      if (ix == (nx_poly - 1))
        (*dedge_group)[idxDedge_group++]  = nx1 + ny1 + ((6 * nx_poly) + 1) * (ny_poly - 1) + 7 * ix + 5;
      else
        (*dedge_group)[idxDedge_group++]  = nx1 + ny1 + ((6 * nx_poly) + 1) * (ny_poly - 1) + 7 * ix + 4;
    }
    ++dn_edgeGroupAbs;
  }

  /* - lim 4 - */

  dn_edgeGroupAbs = 0;
  for (PDM_g_num_t iy = 0; iy < ny1/2; iy++) {
    if ((dn_edgeLim24Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim24Rank[local_rank+1] > dn_edgeGroupAbs)) {
      if (iy == 0)
        (*dedge_group)[idxDedge_group++]  = 1;
      else
        (*dedge_group)[idxDedge_group++]  = 2 + nx1/2 + (2*(iy-1)) + 1;
    }
    ++dn_edgeGroupAbs;
  }

  for (PDM_g_num_t iy = 0; iy < ny_poly; iy++) {
    if ((dn_edgeLim24Rank[local_rank] <= dn_edgeGroupAbs) && (dn_edgeLim24Rank[local_rank+1] > dn_edgeGroupAbs)) {
      if (iy == ny_poly - 1)
        (*dedge_group)[idxDedge_group++]  =  nx1 + ny1 + ((6 * nx_poly) + 1) * iy + 6;
      else
        (*dedge_group)[idxDedge_group++]  =  nx1 + ny1 + ((6 * nx_poly) + 1) * iy + 5;
    }
    ++dn_edgeGroupAbs;
  }

  if (verbose){
    PDM_printf ("- dface_vtx_idx : \n");
    for (int i=0; i<(*dn_face)+1; i++)
      PDM_printf ("%d->%d  ", i+1, (*dface_vtx_idx)[i]);
    PDM_printf ("\n");
    PDM_printf ("- dvtx_coord : \n");
    for (int i = 0; i < (*dn_vtx); i++) {
      PDM_printf ("%d-> ", i+1);
      PDM_printf (" %f", (*dvtx_coord)[3*i]);
      PDM_printf (" %f", (*dvtx_coord)[3*i+1]);
      PDM_printf (" %f", (*dvtx_coord)[3*i+2]);
      PDM_printf ("\n");
    }
    PDM_printf ("- dface_vtx : \n");
    for (int i = 0; i <(*dn_face); i++) {
      PDM_printf ("%d-> ", i+1);
      for (int j = (*dface_vtx_idx)[i]; j < (*dface_vtx_idx)[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, (*dface_vtx)[j]);
      PDM_printf ("\n");
    }
    PDM_printf ("- dface_edge : \n");
    for (int i = 0; i < (*dn_face); i++) {
      PDM_printf ("%d-> ", i+1);
      for (int j = (*dface_vtx_idx)[i]; j < (*dface_vtx_idx)[i+1]; j++)
        PDM_printf (" "PDM_FMT_G_NUM, (*dface_edge)[j]);
      PDM_printf ("\n");
    }

    PDM_printf ("- dedge_vtx : \n");
    for (int i = 0; i < (*dn_edge); i++) {
      PDM_printf ("%d-> ", i+1);
      PDM_printf (" "PDM_FMT_G_NUM, (*dedge_vtx)[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, (*dedge_vtx)[2*i+1]);
      PDM_printf ("\n");
    }

    PDM_printf ("- dedge_face : \n");
    for (int i = 0; i < (*dn_edge); i++) {
      PDM_printf ("%d-> ", i);
      PDM_printf (" "PDM_FMT_G_NUM, (*dedge_face)[2*i]);
      PDM_printf (" "PDM_FMT_G_NUM, (*dedge_face)[2*i+1]);
      PDM_printf ("\n");
    }
  }


  free(dn_edgeLim13Rank);
  free(dn_edgeLim24Rank);
  free(dn_vtx_rank);
  free(dn_face_rank);
  free(dn_edge_rank);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
