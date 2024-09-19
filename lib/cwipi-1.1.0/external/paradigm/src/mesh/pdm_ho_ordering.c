/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_error.h"
#include "pdm_mesh_nodal.h"
#include "pdm_array.h"
#include "pdm_logging.h"

#include "pdm_ho_ordering.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Fortran function header
 *============================================================================*/

/*============================================================================
 * Local macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

typedef struct {

  int  order;
  int *user_to_ijk;
  int *ijk_to_user;

} _ordering_t;

typedef struct {

  char           *name;
  PDM_hash_tab_t *elt_ordering[PDM_MESH_NODAL_N_ELEMENT_TYPES];

} PDM_ho_ordering_t;

/*============================================================================
 * Global variable
 *============================================================================*/

/**
 * Storage of high-order node orderings
 */
static int key_max_order  = 4;
static int s_ho_orderings = 0;
static int n_ho_orderings = 0;
static PDM_ho_ordering_t **ho_orderings = NULL;

#define n_default_orderings 2

static const char *default_orderings_names[n_default_orderings] = {
  "PDM_HO_ORDERING_VTK",
  "PDM_HO_ORDERING_CGNS"
};

static const int _bar_o1_def_to_ijk[n_default_orderings][2] = {
  {0, 1}, // VTK LINE 3
  {0, 1}  // CGNS BAR_2
};

static const int _bar_o2_def_to_ijk[n_default_orderings][3] = {
  {0, 2, 1}, // VTK LAGRANGE_CURVE 68
  {0, 2, 1}  // CGNS BAR_3
};

static const int _bar_o3_def_to_ijk[n_default_orderings][4] = {
  {0, 3, 1, 2}, // VTK LAGRANGE_CURVE 68
  {0, 3, 1, 2}  // CGNS BAR_4
};


static const int _tria_o1_def_to_ijk[n_default_orderings][6] = {
  {0, 0,  1, 0,  0, 1}, // VTK TRIANGLE 5
  {0, 0,  1, 0,  0, 1}  // CGNS TRI_3
};

static const int _tria_o2_def_to_ijk[n_default_orderings][12] = {
  {0, 0,  2, 0,  0, 2,  1, 0,  1, 1,  0, 1}, // VTK LAGRANGE_TRIANGLE 69
  {0, 0,  2, 0,  0, 2,  1, 0,  1, 1,  0, 1}  // CGNS TRI_6
};

static const int _tria_o3_def_to_ijk[n_default_orderings][20] = {
  {0, 0,  3, 0,  0, 3,  1, 0,  2, 0,  2, 1,  1, 2,  0, 2,  0, 1,  1, 1}, // VTK LAGRANGE_TRIANGLE 69
  {0, 0,  3, 0,  0, 3,  1, 0,  2, 0,  2, 1,  1, 2,  0, 2,  0, 1,  1, 1}  // CGNS TRI_10
};


static const int _quad_o1_def_to_ijk[n_default_orderings][8] = {
  {0, 0,  1, 0,  1, 1,  0, 1}, // VTK QUAD 9
  {0, 0,  1, 0,  1, 1,  0, 1}  // CGNS QUAD_4
};

static const int _quad_o2_def_to_ijk[n_default_orderings][18] = {
  {0, 0,  2, 0,  2, 2,  0, 2,  1, 0,  2, 1,  1, 2,  0, 1,  1, 1}, // VTK LAGRANGE_QUADRILATERAL 70
  {0, 0,  2, 0,  2, 2,  0, 2,  1, 0,  2, 1,  1, 2,  0, 1,  1, 1}  // CGNS QUAD_9
};

static const int _quad_o3_def_to_ijk[n_default_orderings][32] = {
  {0, 0,  3, 0,  3, 3,  0, 3,  1, 0,  2, 0,  3, 1,  3, 2,  1, 3,  2, 3,  0, 1,  0, 2,  1, 1,  2, 1,  1, 2,  2, 2}, // VTK LAGRANGE_QUADRILATERAL 70
  {0, 0,  3, 0,  3, 3,  0, 3,  1, 0,  2, 0,  3, 1,  3, 2,  2, 3,  1, 3,  0, 2,  0, 1,  1, 1,  2, 1,  2, 2,  1, 2}  // CGNS QUAD_16
};


static const int _tetra_o1_def_to_ijk[n_default_orderings][12] = {
  {0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1}, // VTK TETRA 10
  {0, 0, 0,  1, 0, 0,  0, 1, 0,  0, 0, 1}  // CGNS TETRA_4
};

static const int _tetra_o2_def_to_ijk[n_default_orderings][30] = {
  {0, 0, 0,  2, 0, 0,  0, 2, 0,  0, 0, 2,  1, 0, 0,  1, 1, 0,  0, 1, 0,  0, 0, 1,  1, 0, 1,  0, 1, 1}, // VTK LAGRANGE_TETRAHEDRON 71
  {0, 0, 0,  2, 0, 0,  0, 2, 0,  0, 0, 2,  1, 0, 0,  1, 1, 0,  0, 1, 0,  0, 0, 1,  1, 0, 1,  0, 1, 1}  // CGNS TETRA_10
};

static const int _tetra_o3_def_to_ijk[n_default_orderings][60] = {
  {0, 0, 0,  3, 0, 0,  0, 3, 0,  0, 0, 3,
   1, 0, 0,  2, 0, 0,
   2, 1, 0,  1, 2, 0,
   0, 2, 0,  0, 1, 0,
   0, 0, 1,  0, 0, 2,
   2, 0, 1,  1, 0, 2,
   0, 2, 1,  0, 1, 2,
   1, 0, 1,
   1, 1, 1,
   0, 1, 1,
   1, 1, 0}, // VTK LAGRANGE_TETRAHEDRON 71
  {0, 0, 0,  3, 0, 0,  0, 3, 0,  0, 0, 3,  1, 0, 0,  2, 0, 0,  2, 1, 0,  1, 2, 0,  0, 2, 0,  0, 1, 0,  0, 0, 1,  0, 0, 2,  2, 0, 1,  1, 0, 2,  0, 2, 1,  0, 1, 2,  1, 1, 0,  1, 0, 1,  1, 1, 1,  0, 1, 1} // CGNS TETRA_20
};


static const int _pyramid_o1_def_to_ijk[n_default_orderings][15] = {
  {0, 0, 0,  1, 0, 0,  1, 1, 0,  0,  1,  0,
   0, 0, 1}, // VTK PYRAMID 14

  {0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,  0, 0, 1}  // CGNS PYRA_5
};

static const int _pyramid_o2_def_to_ijk[n_default_orderings][42] = {
  {0, 0, 0,  2, 0, 0,  2, 2, 0,  0, 2, 0,  0, 0, 2,
   1, 0, 0,  2, 1, 0,  1, 2, 0,  0, 1, 0,
   0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1,
   1, 1, 0}, // VTK LAGRANGE_PYRAMID 74

  {0, 0, 0,  2, 0, 0,  2, 2, 0,  0, 2, 0,  0, 0, 2,
   1, 0, 0,  2, 1, 0,  1, 2, 0,  0, 1, 0,
   0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1,
   1, 1, 0}  // CGNS PYRA_14
};

static const int _pyramid_o3_def_to_ijk[n_default_orderings][90] = {
  {0, 0, 0,  3, 0, 0,  3, 3, 0,  0, 3, 0,  0, 0, 3,
   1, 0, 0,  2, 0, 0,
   3, 1, 0,  3, 2, 0,
   2, 3, 0,  1, 3, 0,
   0, 2, 0,  0, 1, 0,
   0, 0, 1,  0, 0, 2,
   2, 0, 1,  1, 0, 2,
   2, 2, 1,  1, 1, 2,
   0, 2, 1,  0, 1, 2,
   1, 1, 0,  2, 1, 0,
   2, 2, 0,  1, 2, 0,
   1, 0, 1,  2, 1, 1,  1, 2, 1,  0, 1, 1,
   1, 1, 1}, // VTK LAGRANGE_PYRAMID 74

  {0, 0, 0,  3, 0, 0,  3, 3, 0,  0, 3, 0,  0, 0, 3,
   1, 0, 0,  2, 0, 0,
   3, 1, 0,  3, 2, 0,
   2, 3, 0,  1, 3, 0,
   0, 2, 0,  0, 1, 0,
   0, 0, 1,  0, 0, 2,
   2, 0, 1,  1, 0, 2,
   2, 2, 1,  1, 1, 2,
   0, 2, 1,  0, 1, 2,
   1, 1, 0,  2, 1, 0,
   2, 2, 0,  1, 2, 0,
   1, 0, 1,  2, 1, 1,  1, 2, 1,  0, 1, 1,
   1, 1, 1}  // CGNS PYRA_30
};


static const int _prism_o1_def_to_ijk[n_default_orderings][18] = {
  {0, 0, 0,  1, 0, 0,  0, 1, 0,
   0, 0, 1,  1, 0, 1,  0, 1, 1}, // VTK WEDGE 13

  {0, 0, 0,  1, 0, 0,  0, 1, 0,
   0, 0, 1,  1, 0, 1,  0, 1, 1}  // CGNS PENTA_6
};

static const int _prism_o2_def_to_ijk[n_default_orderings][54] = {
  {0, 0, 0,  2, 0, 0,  0, 2, 0,
   0, 0, 2,  2, 0, 2,  0, 2, 2,
   1, 0, 0,  1, 1, 0,  0, 1, 0,
   1, 0, 2,  1, 1, 2,  0, 1, 2,
   0, 0, 1,  2, 0, 1,  0, 2, 1,
   1, 0, 1,  1, 1, 1,  0, 1, 1}, // VTK LAGRANGE_WEDGE 73

  {0, 0, 0,  2, 0, 0,  0, 2, 0,
   0, 0, 2,  2, 0, 2,  0, 2, 2,
   1, 0, 0,  1, 1, 0,  0, 1, 0,
   0, 0, 1,  2, 0, 1,  0, 2, 1,
   1, 0, 2,  1, 1, 2,  0, 1, 2,
   1, 0, 1,  1, 1, 1,  0, 1, 1}  // CGNS PENTA_18
};

static const int _prism_o3_def_to_ijk[n_default_orderings][120] = {
  {0, 0, 0,  3, 0, 0,  0, 3, 0,
   0, 0, 3,  3, 0, 3,  0, 3, 3,

   1, 0, 0,  2, 0, 0,
   2, 1, 0,  1, 2, 0,
   0, 2, 0,  0, 1, 0,

   1, 0, 3,  2, 0, 3,
   2, 1, 3,  1, 2, 3,
   0, 2, 3,  0, 1, 3,

   0, 0, 1,  0, 0, 2,
   3, 0, 1,  3, 0, 2,
   0, 3, 1,  0, 3, 2,

   1, 1, 0,
   1, 1, 3,
   1, 0, 1,  2, 0, 1,  1, 0, 2,  2, 0, 2,
   2, 1, 1,  1, 2, 1,  2, 1, 2,  1, 2, 2,
   0, 1, 1,  0, 2, 1,  0, 1, 2,  0, 2, 2,
   1, 1, 1,
   1, 1, 2}, // VTK LAGRANGE_WEDGE 73

  {0, 0, 0,  3, 0, 0,  0, 3, 0,
   0, 0, 3,  3, 0, 3,  0, 3, 3,
   1, 0, 0,  2, 0, 0,
   2, 1, 0,  1, 2, 0,
   0, 2, 0,  0, 1, 0,
   0, 0, 1,  0, 0, 2,
   3, 0, 1,  3, 0, 2,
   0, 3, 1,  0, 3, 2,
   1, 0, 3,  2, 0, 3,
   2, 1, 3,  1, 2, 3,
   0, 2, 3,  0, 1, 3,
   1, 1, 0,
   1, 0, 1,  2, 0, 1,
   2, 0, 2,  1, 0, 2,
   2, 1, 1,  1, 2, 1,
   1, 2, 2,  2, 1, 2,
   0, 2, 1,  0, 1, 1,
   0, 1, 2,  0, 2, 2,
   1, 1, 3,
   1, 1, 1,
   1, 1, 2}  // CGNS PENTA_40
};


static const int _hexa_o1_def_to_ijk[n_default_orderings][24] = {
  {0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,  0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1}, // VTK HEXAHEDRON 12
  {0, 0, 0,  1, 0, 0,  1, 1, 0,  0, 1, 0,  0, 0, 1,  1, 0, 1,  1, 1, 1,  0, 1, 1}  // CGNS HEXA_8
};

static const int _hexa_o2_def_to_ijk[n_default_orderings][81] = {
  {0, 0, 0,  2, 0, 0,  2, 2, 0,  0, 2, 0,
   0, 0, 2,  2, 0, 2,  2, 2, 2,  0, 2, 2,
   1, 0, 0,  2, 1, 0,  1, 2, 0,  0, 1, 0,
   1, 0, 2,  2, 1, 2,  1, 2, 2,  0, 1, 2,
   0, 0, 1,  2, 0, 1,  0, 2, 1,  2, 2, 1,
   0, 1, 1,  2, 1, 1,
   1, 0, 1,  1, 2, 1,
   1, 1, 0,  1, 1, 2,
   1, 1, 1}, // VTK LAGRANGE_HEXAHEDRON 72

  {0, 0, 0,  2, 0, 0,  2, 2, 0,  0, 2, 0,  0, 0, 2,  2, 0, 2,  2, 2, 2,  0, 2, 2,
   1, 0, 0,  2, 1, 0,  1, 2, 0,  0, 1, 0,
   0, 0, 1,  2, 0, 1,  2, 2, 1,  0, 2, 1,
   1, 0, 2,  2, 1, 2,  1, 2, 2,  0, 1, 2,
   1, 1, 0,
   1, 0, 1,  2, 1, 1,  1, 2, 1,  0, 1, 1,
   1, 1, 2,
   1, 1, 1}  // CGNS HEXA_27
};

static const int _hexa_o3_def_to_ijk[n_default_orderings][192] = {
  {0, 0, 0,  3, 0, 0,  3, 3, 0,  0, 3, 0,  0, 0, 3,  3, 0, 3,  3, 3, 3,  0, 3, 3,
   1, 0, 0,  2, 0, 0,
   3, 1, 0,  3, 2, 0,
   1, 3, 0,  2, 3, 0,
   0, 1, 0,  0, 2, 0,
   1, 0, 3,  2, 0, 3,
   3, 1, 3,  3, 2, 3,
   1, 3, 3,  2, 3, 3,
   0, 1, 3,  0, 2, 3,

   0, 0, 1,  0, 0, 2,
   3, 0, 1,  3, 0, 2,
   0, 3, 1,  0, 3, 2,
   3, 3, 1,  3, 3, 2,

   0, 1, 1,  0, 2, 1,
   0, 1, 2,  0, 2, 2,
   3, 1, 1,  3, 2, 1,
   3, 1, 2,  3, 2, 2,

   1, 0, 1,  2, 0, 1,
   1, 0, 2,  2, 0, 2,
   1, 3, 1,  2, 3, 1,
   1, 3, 2,  2, 3, 2,

   1, 1, 0,  2, 1, 0,
   1, 2, 0,  2, 2, 0,

   1, 1, 3,  2, 1, 3,
   1, 2, 3,  2, 2, 3,

   1, 1, 1,  2, 1, 1,
   1, 2, 1,  2, 2, 1,
   1, 1, 2,  2, 1, 2,
   1, 2, 2,  2, 2, 2
   }, // VTK LAGRANGE_HEXAHEDRON 72

  {0, 0, 0,  3, 0, 0,  3, 3, 0,  0, 3, 0,  0, 0, 3,  3, 0, 3,  3, 3, 3,  0, 3, 3,
   1, 0, 0,  2, 0, 0,
   3, 1, 0,  3, 2, 0,
   2, 3, 0,  1, 3, 0,
   0, 2, 0,  0, 1, 0,
   0, 0, 1,  0, 0, 2,
   3, 0, 1,  3, 0, 2,
   3, 3, 1,  3, 3, 2,
   0, 3, 1,  0, 3, 2,
   1, 0, 3,  2, 0, 3,
   3, 1, 3,  3, 2, 3,
   2, 3, 3,  1, 3, 3,
   0, 2, 3,  0, 1, 3,
   1, 1, 0,  2, 1, 0,
   2, 2, 0,  1, 2, 0,
   1, 0, 1,  2, 0, 1,
   2, 0, 2,  1, 0, 2,
   3, 1, 1,  3, 2, 1,
   3, 2, 2,  3, 1, 2,
   2, 3, 1,  1, 3, 1,
   1, 3, 2,  2, 3, 2,
   0, 2, 1,  0, 1, 1,
   0, 1, 2,  0, 2, 2,
   1, 1, 3,  2, 1, 3,
   2, 2, 3,  1, 2, 3,
   1, 1, 1,  2, 1, 1,
   2, 2, 1,  1, 2, 1,
   1, 1, 2,  2, 1, 2,
   2, 2, 2,  1, 2, 2}  // CGNS HEXA_64
};

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*
static void _vtk_to_ijk_bar (const int order, int *user_to_ijk)
{
  user_to_ijk[0] = 0;
  user_to_ijk[1] = order;
  for (int i = 1; i < order; i++) {
    user_to_ijk[i+1] = 1;
  }
}
*/

static int _check_ijk_to_user (const int n_nodes, const int *ijk_to_user)
{
  for (int i = 0; i < n_nodes; i++) {
    if (ijk_to_user[i] < 0) {
      return 0;
    }
  }

  return 1;
}

static int *_compute_ijk_to_user_bar(const int order, const int *user_to_ijk)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_BARHO, order);

  int *ijk_to_user = PDM_array_const_int (n_nodes, -1);

  for (int i_user = 0; i_user < n_nodes; i_user++) {
    int i_ijk = user_to_ijk[i_user];

    ijk_to_user[i_ijk] = i_user;
  }

  return ijk_to_user;
}

static int *_compute_ijk_to_user_tria(const int order, const int *user_to_ijk)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TRIAHO, order);

  int *ijk_to_user = PDM_array_const_int (2 * n_nodes, -1);

  for (int i_user = 0; i_user < n_nodes; i_user++) {
    int i = user_to_ijk[2*i_user    ];
    int j = user_to_ijk[2*i_user + 1];
    int i_ijk = i + j*(order+1) - j*(j-1)/2;

    ijk_to_user[i_ijk] = i_user;
  }

  return ijk_to_user;
}

static int *_compute_ijk_to_user_quad(const int order, const int *user_to_ijk)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_QUADHO, order);

  int *ijk_to_user = PDM_array_const_int (2 * n_nodes, -1);

  for (int i_user = 0; i_user < n_nodes; i_user++) {
    int i = user_to_ijk[2*i_user    ];
    int j = user_to_ijk[2*i_user + 1];
    int i_ijk = i + j*(order+1);

    ijk_to_user[i_ijk] = i_user;
  }

  return ijk_to_user;
}


static int *_compute_ijk_to_user_tetra(const int order, const int *user_to_ijk)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TETRAHO, order);

  int *ijk_to_user = PDM_array_const_int (3 * n_nodes, -1);

  for (int i_user = 0; i_user < n_nodes; i_user++) {
    int i = user_to_ijk[3*i_user    ];
    int j = user_to_ijk[3*i_user + 1];
    int k = user_to_ijk[3*i_user + 2];
    int i_ijk = i +
      j*(order + 1 - k) - j*(j-1)/2 +
      (k*(k*(k - 3*order - 6) + 3*order*(order + 4) + 11)) / 6;

    ijk_to_user[i_ijk] = i_user;
  }

  return ijk_to_user;
}


static int *_compute_ijk_to_user_pyramid(const int order, const int *user_to_ijk)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_PYRAMIDHO, order);

  int *ijk_to_user = PDM_array_const_int (3 * n_nodes, -1);

  for (int i_user = 0; i_user < n_nodes; i_user++) {
    int i = user_to_ijk[3*i_user    ];
    int j = user_to_ijk[3*i_user + 1];
    int k = user_to_ijk[3*i_user + 2];
    int i_ijk = i +
      j*(order+1-k) +
      (k*(k*(2*k - 6*order - 9) + 6*order*(order + 3) + 13)) / 6;

    ijk_to_user[i_ijk] = i_user;
  }

  return ijk_to_user;
}


static int *_compute_ijk_to_user_prism(const int order, const int *user_to_ijk)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_PRISMHO, order);

  int *ijk_to_user = PDM_array_const_int (3 * n_nodes, -1);

  int n_nodes_iso_k = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_TRIAHO, order);

  for (int i_user = 0; i_user < n_nodes; i_user++) {
    int i = user_to_ijk[3*i_user    ];
    int j = user_to_ijk[3*i_user + 1];
    int k = user_to_ijk[3*i_user + 2];
    int i_ijk = i + j*(order+1) - j*(j-1)/2 + k*n_nodes_iso_k;

    ijk_to_user[i_ijk] = i_user;
  }

  return ijk_to_user;
}

static int *_compute_ijk_to_user_hexa(const int order, const int *user_to_ijk)
{
  int n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (PDM_MESH_NODAL_HEXAHO, order);

  int *ijk_to_user = PDM_array_const_int (3 * n_nodes, -1);

  for (int i_user = 0; i_user < n_nodes; i_user++) {
    int i = user_to_ijk[3*i_user    ];
    int j = user_to_ijk[3*i_user + 1];
    int k = user_to_ijk[3*i_user + 2];
    int i_ijk = i + (order+1)*(j + (order+1)*k);

    ijk_to_user[i_ijk] = i_user;
  }

  return ijk_to_user;
}


static int
_ho_ordering_create
(
 const char *name
 )
{
  int id = n_ho_orderings;
  ho_orderings[id] = malloc (sizeof(PDM_ho_ordering_t));

  PDM_ho_ordering_t *hoo = ho_orderings[id];
  hoo->name = malloc (sizeof(char) * (strlen(name) + 1));
  strcpy(hoo->name, name);

  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT;
       type < PDM_MESH_NODAL_N_ELEMENT_TYPES;
       type++) {
    hoo->elt_ordering[type] = NULL;
  }

  n_ho_orderings++;

  return id;
}


static void
_ordering_free
(
 _ordering_t *ord
 )
{
  if (ord == NULL) return;

  if (ord->user_to_ijk != NULL) free (ord->user_to_ijk);
  if (ord->ijk_to_user != NULL) free (ord->ijk_to_user);

  free (ord);
}

static void
_ho_ordering_free
(
 PDM_ho_ordering_t *hoo
 )
{
  if (hoo == NULL) return;

  for (PDM_Mesh_nodal_elt_t type = PDM_MESH_NODAL_POINT;
       type < PDM_MESH_NODAL_N_ELEMENT_TYPES;
       type++) {

    if (hoo->elt_ordering[type] != NULL) {

      int n_key = PDM_hash_tab_n_used_keys_get (hoo->elt_ordering[type]);
      PDM_g_num_t *keys = PDM_hash_tab_used_keys_get(hoo->elt_ordering[type]);
      for (int k = 0; k < n_key; k++) {
        int key = (int) keys[k];

        const int n_data = PDM_hash_tab_n_data_get (hoo->elt_ordering[type],
                                                    (void *) &key);

        _ordering_t **data = (_ordering_t **) PDM_hash_tab_data_get (hoo->elt_ordering[type],
                                                                     (void *) &key);
        for (int i = 0; i < n_data; i++) {
          _ordering_free (data[i]);
        }
      }

      PDM_hash_tab_free (hoo->elt_ordering[type]);
    }

  }

  free (hoo->name);
  free (hoo);
}




/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Initialize the structure that stores HO orderings
 */

void
PDM_ho_ordering_init
(
void
 )
{
  //printf(">> PDM_ho_ordering_init\n");
  if (ho_orderings == NULL) {

    s_ho_orderings = 2*n_default_orderings;
    n_ho_orderings = 0;
    ho_orderings = malloc (sizeof(PDM_ho_ordering_t *) * s_ho_orderings);


    int order, n_nodes;
    PDM_Mesh_nodal_elt_t t_elt;

    for (int idef = 0; idef < n_default_orderings; idef++) {

      t_elt = PDM_MESH_NODAL_BARHO;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _bar_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _bar_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _bar_o3_def_to_ijk[idef]);

      t_elt = PDM_MESH_NODAL_BARHO_BEZIER;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _bar_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _bar_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _bar_o3_def_to_ijk[idef]);



      t_elt = PDM_MESH_NODAL_TRIAHO;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tria_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tria_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tria_o3_def_to_ijk[idef]);

      t_elt = PDM_MESH_NODAL_TRIAHO_BEZIER;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tria_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tria_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tria_o3_def_to_ijk[idef]);



      t_elt = PDM_MESH_NODAL_QUADHO;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _quad_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _quad_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _quad_o3_def_to_ijk[idef]);



      t_elt = PDM_MESH_NODAL_TETRAHO;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tetra_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tetra_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _tetra_o3_def_to_ijk[idef]);


      t_elt = PDM_MESH_NODAL_PYRAMIDHO;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _pyramid_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _pyramid_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _pyramid_o3_def_to_ijk[idef]);


      t_elt = PDM_MESH_NODAL_PRISMHO;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _prism_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _prism_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _prism_o3_def_to_ijk[idef]);


      t_elt = PDM_MESH_NODAL_HEXAHO;

      order = 1;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _hexa_o1_def_to_ijk[idef]);

      order = 2;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _hexa_o2_def_to_ijk[idef]);

      order = 3;
      n_nodes = PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order);
      PDM_ho_ordering_user_to_ijk_add (default_orderings_names[idef],
                                       t_elt,
                                       order,
                                       n_nodes,
                                       _hexa_o3_def_to_ijk[idef]);
    }

    atexit (PDM_ho_ordering_free);
  }
}



/**
 * \brief Free the structure that stores HO orderings
 *
 * This function is automatically called upon exit of
 * the program which called \ref PDM_ho_ordering_init.
 *
 */

void
PDM_ho_ordering_free
(
void
 )
{
  //printf(">> PDM_ho_ordering_free\n");
  if (ho_orderings == NULL) return;

  for (int i = 0; i < n_ho_orderings; i++) {
    _ho_ordering_free (ho_orderings[i]);
  }

  free (ho_orderings);
}


/**
 * \brief Add a user-defined HO ordering from the locations
 * in the reference uvw-grid of a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 * \param[in] n_nodes      Number of nodes in the high-order element
 * \param[in] user_to_ijk  IJK-coordinates of the nodes in the high-order element
 *
 * \return                 Id of the HO ordering
 */

int
PDM_ho_ordering_user_to_ijk_add
(
 const char                 *name,
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                   n_nodes,
 const int                  *user_to_ijk
 )
{
  assert (n_nodes == PDM_Mesh_nodal_n_vtx_elt_get (t_elt, order));
  int id = PDM_ho_ordering_id_get (name);

  if (id < 0) {
    /* Create new user ordering */
    if (n_ho_orderings >= s_ho_orderings) {
      s_ho_orderings *= 2;
      ho_orderings = realloc (ho_orderings, sizeof(PDM_ho_ordering_t) * s_ho_orderings);
    }

    id = _ho_ordering_create(name);
  }
  PDM_ho_ordering_t *hoo = ho_orderings[id];

  if (hoo->elt_ordering[t_elt] == NULL) {
    hoo->elt_ordering[t_elt] = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT, (void *) &key_max_order);
  }

  PDM_hash_tab_t *elt_ordering = hoo->elt_ordering[t_elt];


  int key = order % key_max_order;

  const int n_data = PDM_hash_tab_n_data_get (elt_ordering, (void *) &key);

  _ordering_t **data = (_ordering_t **) PDM_hash_tab_data_get (elt_ordering, (void *) &key);
  for (int i = 0; i < n_data; i++) {
    if (data[i]->order == order) {
      PDM_error(__FILE__, __LINE__, 0, "Ordering '%s' already specified for elt %d at order %d\n",
                name, (int) t_elt, order);
    }
  }

  /* Compute ijk -> user */
  int *ijk_to_user = NULL;
  int elt_dim = 0;
  switch(t_elt) {
  case PDM_MESH_NODAL_BARHO:
  case PDM_MESH_NODAL_BARHO_BEZIER:
    ijk_to_user = _compute_ijk_to_user_bar(order, user_to_ijk);
    elt_dim = 1;
    break;
  case PDM_MESH_NODAL_TRIAHO:
  case PDM_MESH_NODAL_TRIAHO_BEZIER:
    ijk_to_user = _compute_ijk_to_user_tria(order, user_to_ijk);
    elt_dim = 2;
    break;
  case PDM_MESH_NODAL_QUADHO:
    ijk_to_user = _compute_ijk_to_user_quad(order, user_to_ijk);
    elt_dim = 2;
    break;
  case PDM_MESH_NODAL_TETRAHO:
    ijk_to_user = _compute_ijk_to_user_tetra(order, user_to_ijk);
    elt_dim = 3;
    break;
  case PDM_MESH_NODAL_PYRAMIDHO:
    ijk_to_user = _compute_ijk_to_user_pyramid(order, user_to_ijk);
    elt_dim = 3;
    break;
  case PDM_MESH_NODAL_PRISMHO:
    ijk_to_user = _compute_ijk_to_user_prism(order, user_to_ijk);
    elt_dim = 3;
    break;
  case PDM_MESH_NODAL_HEXAHO:
    ijk_to_user = _compute_ijk_to_user_hexa(order, user_to_ijk);
    elt_dim = 3;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt %d\n", (int) t_elt);
  }

  if (!_check_ijk_to_user (n_nodes, ijk_to_user)) {
    for (int i = 0; i < n_nodes; i++) {
      printf("ijk %3d -> user %3d\n", i, ijk_to_user[i]);
    }
    PDM_error(__FILE__, __LINE__, 0, "Invalid user_to_ijk for ordering '%s', t_elt %d, order %d\n",
              name, (int) t_elt, order);
  }

  int *_user_to_ijk = malloc (sizeof(int) * n_nodes * elt_dim);
  memcpy (_user_to_ijk, user_to_ijk, sizeof(int) * n_nodes * elt_dim);

  _ordering_t *_ord = malloc (sizeof(_ordering_t));
  _ord->order = order;
  _ord->user_to_ijk = _user_to_ijk;
  _ord->ijk_to_user = ijk_to_user;

  PDM_hash_tab_data_add (elt_ordering, (void *) &key, (void *) _ord);

  return id;
}



/**
 * \brief Get the node locations in the reference uvw-grid
 * for a user-defined HO ordering of a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 *
 * \return                 IJK-coordinates of the nodes in the high-order element
 */

int *
PDM_ho_ordering_user_to_ijk_get
(
 const char                  *name,
 const PDM_Mesh_nodal_elt_t   t_elt,
 const int                    order
 )
{
  int id = PDM_ho_ordering_id_get (name);

  if (id < 0) {
    PDM_error(__FILE__, __LINE__, 0, "Ordering '%s' not defined\n", name);
  }

  PDM_ho_ordering_t *hoo = ho_orderings[id];
  PDM_hash_tab_t *elt_ordering = hoo->elt_ordering[t_elt];

  if (elt_ordering == NULL) {
    PDM_error(__FILE__, __LINE__, 0,
              "Ordering '%s' not defined for elt type %d\n", name, (int) t_elt);
  }


  int key = order % key_max_order;

  const int n_data = PDM_hash_tab_n_data_get (elt_ordering, (void *) &key);

  _ordering_t **data = (_ordering_t **) PDM_hash_tab_data_get (elt_ordering, (void *) &key);
  for (int i = 0; i < n_data; i++) {

    if (data[i]->order == order) {
      return data[i]->user_to_ijk;
    }
  }

  PDM_error(__FILE__, __LINE__, 0,
            "Ordering '%s' not defined for elt type %d at order %d\n", name, (int) t_elt, order);


  return NULL;
}



/**
 * \brief Get the map from ijk HO ordering to a user-defined HO ordering
 * for a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 *
 * \return                 User-defined ordering of the nodes in the high-order element
 */

int *
PDM_ho_ordering_ijk_to_user_get
(
 const char                  *name,
 const PDM_Mesh_nodal_elt_t   t_elt,
 const int                    order
 )
{
  int id = PDM_ho_ordering_id_get (name);

  if (id < 0) {
    PDM_error(__FILE__, __LINE__, 0, "Ordering '%s' not defined\n", name);
  }

  PDM_ho_ordering_t *hoo = ho_orderings[id];
  PDM_hash_tab_t *elt_ordering = hoo->elt_ordering[t_elt];

  if (elt_ordering == NULL) {
    PDM_error(__FILE__, __LINE__, 0,
              "Ordering '%s' not defined for elt type %d\n", name, (int) t_elt);
  }


  int key = order % key_max_order;

  const int n_data = PDM_hash_tab_n_data_get (elt_ordering, (void *) &key);

  _ordering_t **data = (_ordering_t **) PDM_hash_tab_data_get (elt_ordering, (void *) &key);
  for (int i = 0; i < n_data; i++) {

    if (data[i]->order == order) {
      return data[i]->ijk_to_user;
    }
  }

  PDM_error(__FILE__, __LINE__, 0,
            "Ordering '%s' not defined for elt type %d at order %d\n", name, (int) t_elt, order);

  return NULL;
}



/**
 * \brief Compute the map from ijk HO ordering to a user-defined HO ordering
 * for a given element type of a given order
 *
 * \param[in] name         Name of the HO ordering
 * \param[in] t_elt        Element type
 * \param[in] order        Element order
 *
 * \return                 User-defined ordering of the nodes in the high-order element
 */

int *
PDM_ho_ordering_compute_ijk_to_user
(
 const PDM_Mesh_nodal_elt_t  t_elt,
 const int                   order,
 const int                  *user_to_ijk
 )
{
  int *ijk_to_user = NULL;
  switch(t_elt) {
  case PDM_MESH_NODAL_BARHO:
    ijk_to_user = _compute_ijk_to_user_bar(order, user_to_ijk);
    break;
  case PDM_MESH_NODAL_TRIAHO:
    ijk_to_user = _compute_ijk_to_user_tria(order, user_to_ijk);
    break;
  case PDM_MESH_NODAL_QUADHO:
    ijk_to_user = _compute_ijk_to_user_quad(order, user_to_ijk);
    break;
  case PDM_MESH_NODAL_TETRAHO:
    ijk_to_user = _compute_ijk_to_user_tetra(order, user_to_ijk);
    break;
  case PDM_MESH_NODAL_PYRAMIDHO:
    ijk_to_user = _compute_ijk_to_user_pyramid(order, user_to_ijk);
    break;
  case PDM_MESH_NODAL_PRISMHO:
    ijk_to_user = _compute_ijk_to_user_prism(order, user_to_ijk);
    break;
  case PDM_MESH_NODAL_HEXAHO:
    ijk_to_user = _compute_ijk_to_user_hexa(order, user_to_ijk);
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid t_elt %d\n", (int) t_elt);
  }

  return ijk_to_user;
}


/**
 * \brief Return the ID of an HO ordering
 *
 * \param[in] name         Name of the HO ordering
 *
 * \return   ID of the HO ordering (-1 if not found)
 */

int
PDM_ho_ordering_id_get
(
 const char *name
 )
{
  if (ho_orderings == NULL) {
    PDM_ho_ordering_init ();
  }

  int id = -1;

  for (int i = 0; i < n_ho_orderings; i++) {
    if (strcmp(ho_orderings[i]->name, name) == 0) {
      id = i;
      break;
    }
  }

  return id;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
