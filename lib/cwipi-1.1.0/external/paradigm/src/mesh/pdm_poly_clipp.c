/*
 * This file is an implementation of Greiner & Hormann algorithm with Foster
 * Overleft extension to remove degenerate cases
 *
 */

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_poly_clipp.h"
#include "pdm_poly_clipp_priv.h"
#include "pdm_binary_search.h"
#include "pdm_plane.h"
#include "pdm_polygon.h"
#include "pdm_priv.h"
#include "pdm_edges_intersect.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef	__cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function definition
 *============================================================================*/


/**
 *
 * \brief Replace default C modulo function (%) to taking into account
 *
 * \param [in]   val    Value
 * \param [in]   mod    Mod
 *
 * \return modulo
 */

static int
_modulo
(
 int val,
 int mod
)
{
  if (val >= 0) {
    return val % mod;
  }
  else {
    return val + mod * ((mod - val - 1)/mod);
  }
}


/**
 *
 * \brief Create a new vertex linked list
 *
 * \param [in]    coords   Vertex coordinates
 * \param [in]    charLgthVtxA      Characteristic length vertex
 * \param [in]    type         Type of point
 * \param [in]    gN       Global number of vertex
 * \param [in]    gNEdge   Global number of edge
 *
 * \return new \ref _vertex_poly_t object
 *
 */
static _vertex_poly_t *
_poly_clipp_new
(
const double *coords,
const PDM_g_num_t gN,
const PDM_g_num_t gNEdge
)
{
  _vertex_poly_t *vtxp = malloc (sizeof(_vertex_poly_t));

  vtxp->coords       = coords;
  vtxp->u            = 0.;
  vtxp->gN           = gN;
  vtxp->gNEdge       = gNEdge;
  vtxp->next         = vtxp;
  vtxp->previous     = vtxp;
  vtxp->used         = false;
  vtxp->isect        = false;
  vtxp->tag          = false;
  vtxp->neighbor     = NULL;
  vtxp->first        = vtxp;
  vtxp->cp           = NULL;
  return vtxp;
}


/**
 *
 * \brief Add a new vertex in linked list
 *
 * \param [in]    coords       Vertex coordinates
 * \param [in]    charLgthVtxA Characteristic length vertex
 * \param [in]    type         Type of point
 * \param [in]    gN           Global number of vertex
 * \param [in]    gNEdge       Global number of edge
 * \param [in]    loc          Location of new vertex
 * \param [in]    linked_vtxp  linked element
 *
 *
 * \return new added \ref _vertex_poly_t object
 *
 */

static _vertex_poly_t *
_poly_clipp_add
(
const double      *coords,
const PDM_g_num_t   gN,
const PDM_g_num_t   gNEdge,
const _poly_clipp_loc_t  loc,
_vertex_poly_t    *linked_vtxp
)
{
  _vertex_poly_t *vtxp = malloc (sizeof(_vertex_poly_t));

  vtxp->coords       = coords;
  vtxp->u            = 0.;
  vtxp->gN           = gN;
  vtxp->gNEdge       = gNEdge;
  vtxp->used         = false;
  vtxp->first        = linked_vtxp->first;

  if (loc == POLY_CLIPP_LOC_NEXT) {
    vtxp->previous       = linked_vtxp;
    vtxp->next           = linked_vtxp->next;
    linked_vtxp->next        = vtxp;
    vtxp->next->previous = vtxp;
  }
  else { //loc == POLY_CLIPP_LOC_PREVIOUS
    vtxp->next            = linked_vtxp;
    vtxp->previous        = linked_vtxp->previous;
    linked_vtxp->previous = vtxp;
    vtxp->previous->next  = vtxp;
  }

  vtxp->isect    = false;
  vtxp->tag      = false;
  vtxp->neighbor = NULL;
  return vtxp;
}


/**
 *
 * \brief Add a new intersection in linked list
 *
 * \param [in]    u            Parameterization in current edge
 * \param [in]    i            New point index
 * \param [in]    type         Type of point
 * \param [in]    loc          Location of new vertex
 * \param [in]    linked_vtxp  Linked element
 *
 *
 * \return new added \ref _vertex_poly_t object
 *
 */

static _vertex_poly_t *
_poly_clipp_intersect_add
(
const double             u,
const PDM_g_num_t         gN,
const _poly_clipp_loc_t  loc,
_vertex_poly_t          *linked_vtxp
)
{
  _vertex_poly_t *vtxp = malloc (sizeof(_vertex_poly_t));

  vtxp->coords       = NULL;
  vtxp->u            = u;
  vtxp->gN           = gN;
  vtxp->gNEdge       = -1;
  vtxp->used         = false;
  vtxp->first        = linked_vtxp->first;

  if (loc == POLY_CLIPP_LOC_NEXT) {
    vtxp->previous        = linked_vtxp;
    vtxp->next            = linked_vtxp->next;
    linked_vtxp->next     = vtxp;
    vtxp->next->previous  = vtxp;
  }
  else {//loc == POLY_CLIPP_LOC_PREVIOUS
    vtxp->next            = linked_vtxp;
    vtxp->previous        = linked_vtxp->previous;
    linked_vtxp->previous = vtxp;
    vtxp->previous->next  = vtxp;
  }

  vtxp->isect    = true;
  vtxp->tag      = true;
  vtxp->neighbor = NULL;
  return vtxp;
}


/**
 *
 * \brief Link neighbors
 *
 * \param [in] vtxpA  A intersection
 * \param [in] vtxpB  B intersection
 *
 */

static void
_poly_clipp_link_neighbors
(
 _vertex_poly_t * vtxpA,
 _vertex_poly_t * vtxpB
)
{
	if (vtxpB == NULL) abort();
  vtxpA->neighbor = vtxpB;
	if (vtxpA == NULL) abort();
  vtxpB->neighbor = vtxpA;
  vtxpA->isect    = true;
  vtxpB->isect    = true;
  vtxpA->tag      = true;
  vtxpB->tag      = true;
}


/**
 *
 * \brief Remove a vertex in linked list
 *
 * \param [in]  vtxp   Vertex to remove
 *
 * \return NULL
 *
 */

static _vertex_poly_t *
_poly_clipp_remove
(
_vertex_poly_t * vtxp
)
{
 if (vtxp != NULL) {
   vtxp->previous->next = vtxp->next;
   vtxp->next->previous = vtxp->previous;
 }
 free (vtxp);
 return NULL;
}


/**
 *
 * \brief Unset intersection vertex (Keep old neighbor in memory)
 *
 * \param [inout] current Current intersection
 *
 */

static void
_poly_clipp_unset_intersect
(
 _vertex_poly_t *current
)
{
  current->isect           = false;
  current->neighbor->isect = false;
}


/**
 *
 * \brief Copy a vertex poly list
 *
 * \return NULL
 *
 */


static _vertex_poly_t *
_poly_clipp_copy
(
_vertex_poly_t *vtxp
)
{
  _vertex_poly_t *vtxp_cp = NULL;
  _vertex_poly_t *cp_prev = NULL;
  _vertex_poly_t *cp_next = NULL;

  if (vtxp != NULL) {
    _vertex_poly_t *vtxp_curr = vtxp;

    do {

      _vertex_poly_t *cp_current = malloc(sizeof(_vertex_poly_t));
      if (vtxp_cp == NULL) {
        vtxp_cp = cp_current;
      }
      vtxp_curr->cp = cp_current;
      cp_current->coords = vtxp_curr->coords;
      cp_current->u = vtxp_curr->u;
      cp_current->first = vtxp_curr->first;
      cp_current->gN = vtxp_curr->gN;
      cp_current->gNEdge = vtxp_curr->gNEdge;
      if (cp_next != NULL) {
        cp_current->next = cp_next;
      }
      else {
        cp_current->next = cp_current;
      }
      if (cp_prev != NULL) {
        cp_current->previous = cp_prev;
      }
      else {
        cp_current->previous = cp_current;
      }
      cp_current->isect = vtxp_curr->isect;
      cp_current->tag = vtxp_curr->tag;
      cp_current->neighbor = vtxp_curr->neighbor;
      cp_current->used = vtxp_curr->used;

      if (cp_prev != NULL) {
        cp_prev->next = cp_current;
      }
      if (cp_next != NULL) {
        cp_next->previous = cp_current;
      }

      vtxp_curr = vtxp_curr->next;

      cp_prev = cp_current;
      cp_next = cp_current->next;

    } while (vtxp_curr != vtxp);



  }
  return vtxp_cp;
}


/**
 *
 * \brief Free a linked list
 *
 * \return NULL
 *
 */


//FIXME: Voir l'appel ou faire l'appel a _poly_clipp_free

static _vertex_poly_t *
_poly_clipp_free
(
_vertex_poly_t *vtxp
)
{

 while (vtxp->next != vtxp) {
   _poly_clipp_remove (vtxp->next);
 }
 return _poly_clipp_remove (vtxp);
}


/**
 *
 * \brief Perform vertex location (in, out) for no 'on' vertex
 *                (ray tracing algorithm)
 *
 * \param [inout] vtxpA        Subject polygon
 * \param [inout] vtxpB        Constraint polygon
 * \param [in]    face_vtxCooA  A vertex coordinates
 * \param [in]    face_vtxCooB  B vertex coordinates
 * \param [in]    n_vtxA        Number of A vertices
 * \param [in]    n_vtxB        Number of B vertices
 * \param [in]    nA           A normal
 * \param [in]    nB           B normal
 *
 */

static void
_location
(
_vertex_poly_t *vtxpA,
_vertex_poly_t *vtxpB,
double *face_vtxCooA,
double *face_vtxCooB,
int n_vtxA,
int n_vtxB,
double nA[3],
double nB[3]
)
{
  double *boundsA = PDM_polygon_bounds_get (n_vtxA, face_vtxCooA);
  double *boundsB = PDM_polygon_bounds_get (n_vtxB, face_vtxCooB);

  _vertex_poly_t *vtx_currB = vtxpB;
  for (int i = 0; i < n_vtxB; i++) {
    const double *_coo = face_vtxCooB + 3*i;
    PDM_polygon_status_t stat = PDM_polygon_point_in_new (_coo,
                                                      n_vtxA,
                                                      face_vtxCooA,
                                                      boundsA,
                                                      nA);
    if (stat == PDM_POLYGON_INSIDE) {
      vtx_currB->tag = true;
    }
    else if (stat == PDM_POLYGON_OUTSIDE) {
      vtx_currB->tag = false;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error PDM_poly_clipp : degenerated polygon\n");
      abort();
    }
    vtx_currB = vtx_currB->next;
  }

  _vertex_poly_t *vtx_currA = vtxpA;
  for (int i = 0; i < n_vtxA; i++) {
    const double *_coo = face_vtxCooA + 3*i;
    PDM_polygon_status_t stat = PDM_polygon_point_in_new (_coo,
                                                      n_vtxB,
                                                      face_vtxCooB,
                                                      boundsB,
                                                      nB);
    if (stat == PDM_POLYGON_INSIDE) {
      vtx_currA->tag = true;
    }
    else if (stat == PDM_POLYGON_OUTSIDE) {
      vtx_currA->tag = false;
    }
    else {
      PDM_error(__FILE__, __LINE__, 0, "Error PDM_poly_clipp : degenerated polygon\n");
      abort();
    }
    vtx_currA = vtx_currA->next;
  }

  free (boundsA);
  free (boundsB);

}


/**
 *
 * \brief Tag intersection vertices
 *
 * \param [inout] vtxpA        Subject polygon
 * \param [inout] vtxpB        Constraint polygon
 *
 */

static void
_tag_current
(
_vertex_poly_t *current
)
{

  _vertex_poly_t *neighbor = current->neighbor;
  _vertex_poly_t *cp = current->cp;
  _vertex_poly_t *neighbor_cp = current->neighbor->cp;

  assert(cp != NULL); // Original data have to be duplicate
  assert(neighbor_cp != NULL); // Original data have to be duplicate

  bool ent = true;
  bool ext = false;

  bool in = true;
  bool out = false;

  bool prev_on = cp->previous->isect;
  bool next_on = cp->next->isect;

  bool prev_in = !cp->previous->isect && cp->previous->tag;
  bool next_in = !cp->next->isect && cp->next->tag;

  bool prev_out = !cp->previous->isect && !cp->previous->tag;
  bool next_out = !cp->next->isect && !cp->next->tag;

  bool neighbor_prev_on = neighbor_cp->previous->isect;
  bool neighbor_next_on = neighbor_cp->next->isect;

  bool neighbor_prev_in = !neighbor_cp->previous->isect && neighbor_cp->previous->tag;
  bool neighbor_next_in = !neighbor_cp->next->isect && neighbor_cp->next->tag;

  bool neighbor_prev_out = !neighbor_cp->previous->isect && !neighbor_cp->previous->tag;
  bool neighbor_next_out = !neighbor_cp->next->isect && !neighbor_cp->next->tag;

  /*
   * on/on
   */

  if (prev_on && next_on) {

    if (neighbor_prev_on && neighbor_next_on) {
      _poly_clipp_unset_intersect (current);
      current->tag = in;
      neighbor->tag = in;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      _poly_clipp_unset_intersect (current);
      current->tag = in;
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      _poly_clipp_unset_intersect (current);
      current->tag = in;
      neighbor->tag = out;
    }

    else if (neighbor_prev_on && neighbor_next_out) {
      current->tag = ent;
      neighbor->tag = ext;
    }

    else if (neighbor_prev_on && neighbor_next_in) {
      current->tag = ext;
      neighbor->tag = ent;
    }

    else if (neighbor_prev_out && neighbor_next_on) {
      current->tag = ext;
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_on) {
      current->tag = ent;
      neighbor->tag = ext;
    }

    else if (neighbor_prev_out && neighbor_next_in) {
      current->tag = ext;
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_out) {
      current->tag = ent;
      neighbor->tag = ext;
    }

    else {

      abort();
    }

  }

  /*
   * on/out
   */

  else if (prev_on && next_out) {
    current->tag = ext;

    if (neighbor_prev_on && neighbor_next_on) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      neighbor->tag = out;
    }

    else if (neighbor_prev_on && neighbor_next_out) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_on && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_out && neighbor_next_on) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_out && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_out) {
      neighbor->tag = ext;
    }
  }

  /*
   * on/in
   */

  else if (prev_on && next_in) {
    current->tag = ent;

    if (neighbor_prev_on && neighbor_next_on) {
      neighbor->tag = out;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      neighbor->tag = out;
    }

    else if (neighbor_prev_on && neighbor_next_out) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_on && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_out && neighbor_next_on) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_out && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_out) {
      neighbor->tag = ext;
    }
  }

  /*
   * out/on
   */

  else if (prev_out && next_on) {
    current->tag = ent;

    if (neighbor_prev_on && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      neighbor->tag = out;
    }

    else if (neighbor_prev_on && neighbor_next_out) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_on && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_out && neighbor_next_on) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_out && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_out) {
      neighbor->tag = ext;
    }
  }


  /*
   * in/on
   */

  else if (prev_in && next_on) {
    current->tag = ext;

    if (neighbor_prev_on && neighbor_next_on) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      neighbor->tag = out;
      current->tag = in;
    }

    else if (neighbor_prev_on && neighbor_next_out) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_on && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_out && neighbor_next_on) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_out && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_out) {
      neighbor->tag = ext;
    }
  }

  /*
   * out/out
   */

  else if (prev_out && next_out) {

    if (neighbor_prev_on && neighbor_next_on) {
      _poly_clipp_unset_intersect (current);
      current->tag = out;
      neighbor->tag = in;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      _poly_clipp_unset_intersect (current);
      current->tag = out;
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      _poly_clipp_unset_intersect (current);
      current->tag = out;
      neighbor->tag = out;
    }
    else {
      current->tag = ext;
    }

  }

  /*
   * in/in
   */

  else if (prev_in && next_in) {

    if (neighbor_prev_on && neighbor_next_on) {
      _poly_clipp_unset_intersect (current);
      current->tag = in;
      neighbor->tag = in;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      _poly_clipp_unset_intersect (current);
      current->tag = in;
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      _poly_clipp_unset_intersect (current);
      current->tag = in;
      neighbor->tag = out;
    }

  }

  /*
   * out/in
   */

  else if (prev_out && next_in) {
    current->tag = ent;

    if (neighbor_prev_on && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      neighbor->tag = out;
    }

    else if (neighbor_prev_on && neighbor_next_out) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_on && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_out && neighbor_next_on) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_out && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_out) {
      neighbor->tag = ext;
    }
  }

  /*
   * in/out
   */

  else if (prev_in && next_out) {
    current->tag = ext;

    if (neighbor_prev_on && neighbor_next_on) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_in) {
      neighbor->tag = in;
    }

    else if (neighbor_prev_out && neighbor_next_out) {
      neighbor->tag = out;
    }

    else if (neighbor_prev_on && neighbor_next_out) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_on && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_out && neighbor_next_on) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_on) {
      neighbor->tag = ext;
    }

    else if (neighbor_prev_out && neighbor_next_in) {
      neighbor->tag = ent;
    }

    else if (neighbor_prev_in && neighbor_next_out) {
      neighbor->tag = ext;
    }
  }
}


/**
 *
 * \brief Tag intersection vertices
 *
 * \param [inout] vtxpA   Subject polygon
 * \param [inout] vtxpB   Constraint polygon
 * \param [in] nClippVtxA   Subject polygon
 * \param [in] nClippVtxB   Constraint polygon
 *
 */


//FIXME :Fonction _tag : Voir s'il faut faire le travail pour vtxpB

static void
_tag
(
_vertex_poly_t *vtxpA,
_vertex_poly_t *vtxpB,
const int       nClippVtxA,
const int       nClippVtxB
)
{
  PDM_UNUSED(nClippVtxA);
  PDM_UNUSED(nClippVtxB);

  _vertex_poly_t *cp_vtxpA = _poly_clipp_copy(vtxpA);
  _vertex_poly_t *cp_vtxpB = _poly_clipp_copy(vtxpB);

  _vertex_poly_t *vtx_currA = vtxpA;

  _vertex_poly_t * test =  vtxpA;
  do {
    test = test->next;
  } while (test != vtxpA);

  test = cp_vtxpA;
  do {
    test = test->next;
  } while (test != cp_vtxpA);

  test = vtxpB;
  do {
    test = test->next;
  } while (test != vtxpB);

  test = cp_vtxpB;
  do {
    test = test->next;
  } while (test != cp_vtxpB);

  do {
    _vertex_poly_t *nextA = vtx_currA->next;
    if (vtx_currA->isect) {
      _tag_current (vtx_currA);
    }

    vtx_currA = nextA;

  } while (vtx_currA != vtxpA);

  _poly_clipp_free (cp_vtxpA);
  _poly_clipp_free (cp_vtxpB);

}

/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 *
 * \brief Perform polygon clipping
 *
 * \param [in]    ei                   Edges intersection management
 * \param [in]    gNumA                Polygon A global number
 * \param [in]    n_vtxA                Number of polygon A vertices
 * \param [in]    faceToEdgeA          Polygon A face to edge connectivity
 * \param [in]    faceToVtxA           Polygon A face to vertex connectivity
 * \param [in]    face_vtxCooA          Polygon A vertex coordinates
 * \param [in]    face_vtxEpsA          Polygon A vertex characteristic length
 * \param [in]    gNumB                Polygon A global number
 * \param [in]    n_vtxB                Number of polygon B vertices
 * \param [in]    faceToEdgeB          Polygon B face to edge connectivity
 * \param [in]    faceToVtxB           Polygon B face to vertex connectivity
 * \param [in]    face_vtxCooB          Polygon B vertex coordinates
 * \param [in]    face_vtxEpsB          Polygon B vertex characteristic length
 * \param [in]    performed_t          Type of performed polygon :  PDM_POLY_CLIPP_CLIP,  PDM_POLY_CLIPP_REVERSE  // Perform clipped polygons , Perform opposit polygons
 * \param [out]   nPolyClippA           Number of clipped polygon
 * \param [out]   polyClippIdxA         Connectivity index for each polygon
 *                                     size = nPolyClipp + 1
 * \param [out]   polyClippConnecA     Connectivity of each clipped polygon
 *                                     size = polyClippIdx[nPolyClipp] (Vtx from A > 0 and vtx from B < 0)
 * \param [out]   polyClippCoordsA     Vertices coordinates of clipping polygon
 * \param [out]   nPolyClippB          Number of clipped polygon
 * \param [out]   polyClippIdxB         Connectivity index for each polygon
 *                                     size = nPolyClipp + 1
 * \param [out]   polyClippConnecB     Connectivity of each clipped polygon
 *                                     size = polyClippIdx[nPolyClipp] (Vtx from B > 0 and vtx from A < 0)
 * \param [out]   polyClippCoordsB     Vertices coordinates of clipping polygon
 *
 */

void
PDM_poly_clipp
(
PDM_edges_intersect_t  *ei,
PDM_g_num_t             gnum_boxA,
PDM_g_num_t             gnum_boxB,
const int               n_vtxA,
PDM_g_num_t            *faceToEdgeA,
PDM_g_num_t            *faceToVtxA,
double                 *face_vtxCooA,
const int               n_vtxB,
PDM_g_num_t            *faceToEdgeB,
PDM_g_num_t            *faceToVtxB,
double                 *face_vtxCooB,
PDM_poly_clipp_t        performed_t,
int                    *nPolyClippA,
int                   **polyClippIdxA,
PDM_g_num_t           **polyClippConnecA,
double                **polyClippCoordsA,
int                    *nPolyClippB,
int                   **polyClippIdxB,
PDM_g_num_t           **polyClippConnecB,
double                **polyClippCoordsB
)
{
  PDM_g_num_t *_faceToEdgeA = faceToEdgeA;
  PDM_g_num_t *_faceToVtxA  = faceToVtxA;
  double     *_face_vtxCooA = face_vtxCooA;

  PDM_g_num_t *_faceToEdgeB = faceToEdgeB;
  PDM_g_num_t *_faceToVtxB  = faceToVtxB;
  double     *_face_vtxCooB = face_vtxCooB;

  /*
   * Compute Normal
   *
   */

  double nA[3];
  double baryA[3];
  PDM_plane_normal (n_vtxA, face_vtxCooA, nA);
  PDM_plane_barycenter (n_vtxA, face_vtxCooA, baryA);

  double nB[3];
  PDM_plane_normal (n_vtxB, face_vtxCooB, nB);

  double dot = PDM_DOT_PRODUCT (nA, nB);

  bool revert = false;

  if (dot < 0) {
    revert = true;
  }

  /*
   * Reorient if necessary
   *
   */

  if (revert) {

    _faceToEdgeB = malloc (sizeof(PDM_g_num_t) * n_vtxB);
    _faceToVtxB  = malloc (sizeof(PDM_g_num_t) * n_vtxB);
    _face_vtxCooB = malloc (sizeof(double) * 3 * n_vtxB);

    int j = n_vtxB - 1;
    for (int i = 0; i < n_vtxB; i++) {
      _faceToEdgeB[i] = -faceToEdgeB[_modulo((j-1),n_vtxB)];
      _faceToVtxB[i] = faceToVtxB[j];
      for (int k = 0; k < 3; k++) {
        _face_vtxCooB[3*i+k] = face_vtxCooB[3*j+k];
      }
      j += -1;
    }
  }

  _face_vtxCooA = malloc (sizeof(double) * 3 * n_vtxA);
  for (int i = 0; i < n_vtxA; i++) {
    PDM_plane_projection (face_vtxCooA + 3 * i, baryA, nA, _face_vtxCooA + 3 * i);
  }
  PDM_plane_normal (n_vtxA, _face_vtxCooA, nA);

  if (revert) {
    for (int i = 0; i < n_vtxB; i++) {
      PDM_plane_projection (_face_vtxCooB + 3 * i, baryA, nA, _face_vtxCooB + 3 * i);
    }
  }
  else {
    _face_vtxCooB = malloc (sizeof(double) * 3 * n_vtxB);
    for (int i = 0; i < n_vtxB; i++) {
      PDM_plane_projection (face_vtxCooB + 3 * i, baryA, nA, _face_vtxCooB + 3 * i);
    }
  }



  /*
   * Create double linked list vertex structures
   *
   */

  _vertex_poly_t *vtxA = _poly_clipp_new (_face_vtxCooA,
                                          *_faceToVtxA,
                                          PDM_ABS(*_faceToEdgeA));

  _vertex_poly_t *vtxB = _poly_clipp_new (_face_vtxCooB,
                                          *_faceToVtxB,
                                          PDM_ABS(*_faceToEdgeB));

  _vertex_poly_t **vtxA_origin = malloc (sizeof(_vertex_poly_t *) *  n_vtxA);
  _vertex_poly_t **vtxB_origin = malloc (sizeof(_vertex_poly_t *) *  n_vtxB);

  vtxA_origin[0] = vtxA;
  vtxB_origin[0] = vtxB;

  int nClippVtxA = n_vtxA;

  for (int i = 1; i < n_vtxA; i++) {
    const double *_coo   = _face_vtxCooA + 3*i;
    PDM_g_num_t   *_gN     = _faceToVtxA +i;
    PDM_g_num_t   *_gNEdge = _faceToEdgeA +i;

    vtxA_origin[i] = _poly_clipp_add (_coo,
                                      *_gN,
                                      PDM_ABS(*_gNEdge),
                                      POLY_CLIPP_LOC_PREVIOUS,
                                      vtxA);
  }

  int nClippVtxB = n_vtxB;
  for (int i = 1; i < n_vtxB; i++) {
    const double *_coo = _face_vtxCooB + 3*i;
    PDM_g_num_t   *_gN     = _faceToVtxB +i;
    PDM_g_num_t   *_gNEdge = _faceToEdgeB +i;

    vtxB_origin[i] = _poly_clipp_add (_coo,
                                      *_gN,
                                      PDM_ABS(*_gNEdge),
                                      POLY_CLIPP_LOC_PREVIOUS,
                                      vtxB);
  }

  /*
   * Perform vertex location (in, out) for no 'on' vertex
   *    - First step  : In or out (ray tracing algorithm)
   *
   */

  //_location (vtxA, vtxB, _face_vtxCooA, _face_vtxCooB, n_vtxA, n_vtxB, nA, nB);
  _location (vtxA, vtxB, _face_vtxCooA, _face_vtxCooB, n_vtxA, n_vtxB, nA, nA);

  /*
   *   - Add intersection into linked list
   *   - Move to Intersection tag for Vertices located on Polygon
   *
   */

  for (int i1 = 0; i1 < n_vtxA; i1++) {
    _vertex_poly_t *vtx_currA = vtxA_origin[i1];
    _vertex_poly_t *nextA = vtxA_origin[(i1+1) % n_vtxA];

    for (int j1 = 0; j1 < n_vtxB; j1++) {
      _vertex_poly_t *vtx_currB = vtxB_origin[j1];
      _vertex_poly_t *nextB = vtxB_origin[(j1+1) % n_vtxB];

      /*
       * Look for intersection : (compute if not already stored)
       */

      int nData = 0;
      PDM_edges_intersect_res_t **_eir =
        PDM_edges_intersect_get (ei,
                                 PDM_EDGES_GET_FROM_AB,
                                 vtx_currA->gNEdge,
                                 vtx_currB->gNEdge,
                                 &nData);

      if (_eir != NULL) {

        assert (nData == 1);

        PDM_edges_intersect_res_t *eir = *_eir;

        /*
         * Get new points and modified vertices coming from intersections
         */

        PDM_line_intersect_t         tIntersect;

        PDM_g_num_t                   nGEdgeA;
        PDM_g_num_t                   originEdgeA;
        PDM_g_num_t                   endEdgeA;
        int                          nNewPointsA;
        PDM_edges_intersect_point_t *oNewPointsA;
        PDM_g_num_t                  *linkA;
        PDM_g_num_t                  *gNumVtxA;
        double                      *coordsA;
        double                      *uA;

        /*
         * Get intersection properties
         */

        PDM_edges_intersect_res_data_get (eir,
                                          PDM_EDGES_INTERSECT_MESHA,
                                          &nGEdgeA,
                                          &originEdgeA,
                                          &endEdgeA,
                                          &tIntersect,
                                          &nNewPointsA,
                                          &oNewPointsA,
                                          &linkA,
                                          &gNumVtxA,
                                          &coordsA,
                                          &uA);

        PDM_g_num_t                   nGEdgeB;
        PDM_g_num_t                   originEdgeB;
        PDM_g_num_t                   endEdgeB;
        int                          nNewPointsB;
        PDM_edges_intersect_point_t *oNewPointsB;
        PDM_g_num_t                  *linkB;
        PDM_g_num_t                  *gNumVtxB;
        double                      *coordsB;
        double                      *uB;

        PDM_edges_intersect_res_data_get (eir,
                                          PDM_EDGES_INTERSECT_MESHB,
                                          &nGEdgeB,
                                          &originEdgeB,
                                          &endEdgeB,
                                          &tIntersect,
                                          &nNewPointsB,
                                          &oNewPointsB,
                                          &linkB,
                                          &gNumVtxB,
                                          &coordsB,
                                          &uB);

        free (_eir);

        /*
         * Add new intersections vertex or switch vertex to intersection
         * if vertex is on subject polygon
         */

        for (int i = 0; i < nNewPointsA; i++) {

          if ((oNewPointsA[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA)
              || (oNewPointsA[i] == PDM_EDGES_INTERSECT_POINT_NEW)) {

            _vertex_poly_t *_currA = vtx_currA;
            _vertex_poly_t *_next_currA = _currA->next;

            double _u = uA[i];
            if (vtx_currA->gN != originEdgeA) {
              _u = 1. - _u;
            }

            while (_next_currA != nextA) {
              if ((_u < _next_currA->u) ||
                  ((gNumVtxA[i] != 0) && (gNumVtxA[i] == _next_currA->gN))) {
                break;
              }
              _currA = _next_currA;
              _next_currA = _currA->next;
            };

            if ((oNewPointsA[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) &&
                ((gNumVtxA[i] != 0) && (gNumVtxA[i] == _next_currA->gN))) {
              break;
            }

            _vertex_poly_t *interVtxA = _poly_clipp_intersect_add (_u,
                                                                   gNumVtxA[i],
                                                                   POLY_CLIPP_LOC_NEXT,
                                                                   _currA);
            interVtxA->coords = &(coordsA[3*i]);

            *nPolyClippA += 1;

            if (oNewPointsA[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {

              if (vtx_currB->gN == linkA[i]) {
                _poly_clipp_link_neighbors (interVtxA, vtx_currB);
              }

              else {
                assert (nextB->gN == linkA[i]);
                _poly_clipp_link_neighbors (interVtxA, nextB);
              }
            }

            else if (oNewPointsA[i] == PDM_EDGES_INTERSECT_POINT_NEW) {

              _vertex_poly_t *_currB = vtx_currB;
              _vertex_poly_t *_next_currB = _currB->next;

              double _v = uB[i];
              if (vtx_currB->gN != originEdgeB) {
                _v = 1. - _v;
              }
              while (_next_currB != nextB) {
                if (_v < _next_currB->u) {
                  break;
                }
                _currB = _next_currB;
                _next_currB = _currB->next;
              };

              assert (nNewPointsB == 1);

              _vertex_poly_t *interVtxB = _poly_clipp_intersect_add (_v,
                                                                     gNumVtxB[0],
                                                                     POLY_CLIPP_LOC_NEXT,
                                                                     _currB);
              interVtxB->coords = &(coordsB[0]);
              *nPolyClippB += 1;

              _poly_clipp_link_neighbors (interVtxA, interVtxB);

            }
          }

          else if (oNewPointsA[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {

            _vertex_poly_t *_currA = vtx_currA;
            if (_currA->gN != gNumVtxA[i]) {
              assert (nextA->gN == gNumVtxA[i]);
              _currA = nextA;
            }

            _vertex_poly_t *_currB = vtx_currB;
            if (_currB->gN != linkA[i]) {
              assert (nextB->gN == linkA[i]);
              _currB = nextB;
            }

            _poly_clipp_link_neighbors (_currA, _currB);

          }

        }


        for (int i = 0; i < nNewPointsB; i++) {

          if (oNewPointsB[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {

            _vertex_poly_t *_currB = vtx_currB;
            _vertex_poly_t *_next_currB = _currB->next;

            double _u = uB[i];
            if (vtx_currB->gN != originEdgeB) {
              _u = 1. - _u;
            }

            while (_next_currB != nextB) {
              if ((_u < _next_currB->u) ||
                  ((gNumVtxB[i] != 0) && (gNumVtxB[i] == _next_currB->gN))) {
                break;
              }
              _currB = _next_currB;
              _next_currB = _currB->next;
            };

            if ((gNumVtxB[i] != 0) && (gNumVtxB[i] == _next_currB->gN)) {
              break;
            }

            _vertex_poly_t *interVtxB = _poly_clipp_intersect_add (_u,
                                                                   gNumVtxB[i],
                                                                   POLY_CLIPP_LOC_NEXT,
                                                                   _currB);

            interVtxB->coords = &(coordsB[3*i]);
            *nPolyClippB += 1;

            if (vtx_currA->gN == linkB[i]) {
              _poly_clipp_link_neighbors (interVtxB, vtx_currA);
            }

            else {
              assert (nextA->gN == linkB[i]);
              _poly_clipp_link_neighbors (interVtxB, nextA);
            }
          }
        }
      }
    }
  }
  free (vtxA_origin);
  free (vtxB_origin);

  /*
   * Tag Intersection points
   *
   */

  if ((gnum_boxA == 19189330) && (gnum_boxB == 3206132) && 0) {
    printf("Ensemble des points avant tag (gn, isect, neighbour, tag) :\n");
    _vertex_poly_t *vtx_currA = vtxA;
    printf ("Polygon A : ");
    do {
      PDM_g_num_t neighbour = -1;
      if (vtx_currA->neighbor != NULL)
        neighbour = vtx_currA->neighbor->gN;
      if ( neighbour != -1)
        printf (" "PDM_FMT_G_NUM"-%d/"PDM_FMT_G_NUM"/%d", vtx_currA->gN, vtx_currA->isect, neighbour, vtx_currA->tag);
      else
        printf (" "PDM_FMT_G_NUM"-%d/NC/%d", vtx_currA->gN, vtx_currA->isect, vtx_currA->tag);
      vtx_currA = vtx_currA->next;
    } while (vtx_currA != vtxA);
    printf("\n");

    _vertex_poly_t *vtx_currB = vtxB;
    printf ("Polygon B : ");
    do {
      PDM_g_num_t neighbour = -1;
      if (vtx_currB->neighbor != NULL)
        neighbour = vtx_currB->neighbor->gN;
      if ( neighbour != -1)
        printf (" "PDM_FMT_G_NUM"-%d/"PDM_FMT_G_NUM"/%d", vtx_currB->gN, vtx_currB->isect, neighbour, vtx_currB->tag);
      else
        printf (" "PDM_FMT_G_NUM"-%d/NC/%d", vtx_currB->gN, vtx_currB->isect, vtx_currB->tag);
      vtx_currB = vtx_currB->next;
    } while (vtx_currB != vtxB);
    printf("\n");
  }

  _tag (vtxA, vtxB, nClippVtxA, nClippVtxB);

  if ((gnum_boxA == 19189330) && (gnum_boxB == 3206132) && 0) {
    //  if (0 == 1) {
    printf("Ensemble des points apres tag (gn, isect, neighbour, tag) :\n");
    _vertex_poly_t *vtx_currA = vtxA;
    printf ("Polygon A : ");
    do {
      PDM_g_num_t neighbour = -1;
      if (vtx_currA->neighbor != NULL)
        neighbour = vtx_currA->neighbor->gN;
      if ( neighbour != -1)
        printf (" "PDM_FMT_G_NUM"-%d/"PDM_FMT_G_NUM"/%d", vtx_currA->gN, vtx_currA->isect, neighbour, vtx_currA->tag);
      else
        printf (" "PDM_FMT_G_NUM"-%d/NC/%d", vtx_currA->gN, vtx_currA->isect, vtx_currA->tag);
      vtx_currA = vtx_currA->next;
    } while (vtx_currA != vtxA);
    printf("\n");

    _vertex_poly_t *vtx_currB = vtxB;
    printf ("Polygon B : ");
    do {
      PDM_g_num_t neighbour = -1;
      if (vtx_currB->neighbor != NULL)
        neighbour = vtx_currB->neighbor->gN;
      if ( neighbour != -1)
        printf (" "PDM_FMT_G_NUM"-%d/"PDM_FMT_G_NUM"/%d", vtx_currB->gN, vtx_currB->isect, neighbour, vtx_currB->tag);
      else
        printf (" "PDM_FMT_G_NUM"-%d/NC/%d", vtx_currB->gN, vtx_currB->isect, vtx_currB->tag);
      vtx_currB = vtx_currB->next;
    } while (vtx_currB != vtxB);
    printf("\n");
  }

  /*
   * Clipping
   */

  int nPolyPredicA = 2;
  int sPolyConnecA = n_vtxA + n_vtxB;

  *nPolyClippA = 0;
  *polyClippConnecA = NULL;
  *polyClippIdxA = NULL;
  int idx1 = 0;
  int idx2 = 0;

  int nPolyPredicB = 2;
  int sPolyConnecB = n_vtxA + n_vtxB;

  *nPolyClippB = 0;
  *polyClippConnecB = NULL;
  *polyClippIdxB = NULL;

  int sPolyCoordA = 3 * sPolyConnecA;
  int sPolyCoordB = 3 * sPolyConnecB;

  if (performed_t == PDM_POLY_CLIPP_CLIP) {

    *polyClippIdxA = malloc (sizeof(int) * (nPolyPredicA + 1));
    (*polyClippIdxA)[0] = 0;

    *polyClippIdxB = *polyClippIdxA;

    *polyClippConnecA = malloc (sizeof(PDM_g_num_t) * sPolyConnecA);
    *polyClippConnecB = malloc (sizeof(PDM_g_num_t) * sPolyConnecB);

    *polyClippCoordsA = malloc (sizeof(double) * sPolyCoordA);
    *polyClippCoordsB = *polyClippCoordsA;

    sPolyCoordB = sPolyCoordA;

    /* Look for First intersection or point A in B*/

    _vertex_poly_t *first = vtxA;

    while (!first->tag && !first->isect && first->next != vtxA) {
      first = first->next;
    }

    if (first->isect || first->tag) {

      int s_clipped_vtx = n_vtxA + n_vtxB;

      _vertex_poly_t **clipped_vtx = malloc (sizeof(_vertex_poly_t*) * s_clipped_vtx);
      int *clipped_multi = malloc (sizeof(int) * s_clipped_vtx);
      int *origin_vtx = malloc (sizeof(int) * s_clipped_vtx); //  1  : A, -1   : B
                                                              // 10  : A, -10  : B pour les points intersections
                                                              // 100 : A, -100 : B pour les points double (sans intersection)

      int *link_multi = malloc (sizeof(int) * s_clipped_vtx);
      _vertex_poly_t *curr = first;

      do {

        /* Clipping */

        int n_vtxClipp = 0;

        _vertex_poly_t *first_poly = curr;

        int onPolyA = 1;
        int onPolyB = -onPolyA;

        int n_intersection = 0;
        int n_multiplicate = 0;

        if (!curr->used && ((curr->isect) || ((!curr->isect) && (curr->tag)))) {

          int ipass = 0;
          do {

            if (n_vtxClipp >= s_clipped_vtx) {
              while (n_vtxClipp >= s_clipped_vtx) {
                s_clipped_vtx *= 2;
              }
              clipped_vtx =
                realloc (clipped_vtx, sizeof(_vertex_poly_t*) * s_clipped_vtx);
              clipped_multi =
                realloc (clipped_multi, sizeof(int) * s_clipped_vtx);
              origin_vtx =
                realloc (origin_vtx, sizeof(int) * s_clipped_vtx);
              link_multi =
                realloc (link_multi, sizeof(int) * s_clipped_vtx);
            }

            origin_vtx[n_vtxClipp] = onPolyA;
            if (curr->isect) {
              origin_vtx[n_vtxClipp] *= 10;
              n_intersection++;
            }
            else if (curr->neighbor != NULL) {
              origin_vtx[n_vtxClipp] *= 100;
              clipped_multi[n_multiplicate++] = n_vtxClipp;
            }

            link_multi[n_vtxClipp] = -1;

            clipped_vtx[n_vtxClipp++] = curr;
            curr->used = true;

            if ((!curr->isect && curr->tag) ||
                (curr->isect && curr->tag)) {
              curr = curr->next;
            }

            else if (curr->isect && !curr->tag) {
              curr = curr->neighbor->next;
              onPolyA = -onPolyA;
              onPolyB = -onPolyB;
            }

            else {
              printf("Poly clipp erreur : point exterieur dans un sous-poly\n");
              printf("polyA : "PDM_FMT_G_NUM" %d\n", gnum_boxA, n_vtxA);
              for (int i = 0; i < n_vtxA; i++) {
                printf("%20.16e %20.16e %20.16e\n", _face_vtxCooA[3*i], _face_vtxCooA[3*i+1], _face_vtxCooA[3*i+2]);
              }

              printf("polyB : "PDM_FMT_G_NUM" %d\n", gnum_boxB, n_vtxB);
              for (int i = 0; i < n_vtxB; i++) {
                printf("%20.16e %20.16e %20.16e\n", _face_vtxCooB[3*i], _face_vtxCooB[3*i+1], _face_vtxCooB[3*i+2]);
              }

              fflush(stdout);
              abort();
            }

            ipass++;

          } while (curr != first_poly &&
                   ((curr->neighbor != NULL) ? (curr->neighbor != first_poly) : true) &&
                   ipass < 100);

          if (ipass >= 100) {
            printf ("internal poly_clipp error : Erreur dans le parcours des polygones, boucle infinie\n");
            abort();
          }

          int new_n_multiplicate = 0;
          for (int i = 0; i < n_multiplicate; i++) {
            if (link_multi[clipped_multi[i]] == -1) {
              for (int j = i+1; j < n_multiplicate; j++) {
                if (clipped_vtx[clipped_multi[i]]->neighbor  == clipped_vtx[clipped_multi[j]]) {
                  link_multi[clipped_multi[i]] = clipped_multi[j];
                  link_multi[clipped_multi[j]] = clipped_multi[i];
                  new_n_multiplicate++;
                  break;
                }
              }
              if (link_multi[clipped_multi[i]] == -1) {
                origin_vtx[clipped_multi[i]] =
                  origin_vtx[clipped_multi[i]]/PDM_ABS(origin_vtx[clipped_multi[i]]);
              }
            }
          }

          n_multiplicate = new_n_multiplicate;

          // Verification des degenerescence + des points double pour sous-decoupage

          if (n_multiplicate == 0) {

            if (*nPolyClippA >= nPolyPredicA) {
              while (*nPolyClippA >= nPolyPredicA) {
                nPolyPredicA *= 2;
              }
              *polyClippIdxA = realloc (*polyClippIdxA, sizeof(int) * (nPolyPredicA + 1));
              *polyClippIdxB = *polyClippIdxA;
            }

            (*polyClippIdxA)[*nPolyClippA + 1] =
              (*polyClippIdxA)[*nPolyClippA] + n_vtxClipp;

            if ((*polyClippIdxA)[*nPolyClippA+1] >= sPolyConnecA) {
              while ((*polyClippIdxA)[*nPolyClippA+1] >= sPolyConnecA) {
                sPolyConnecA *= 2;
              }
              *polyClippConnecA = realloc (*polyClippConnecA, sizeof(PDM_g_num_t) *sPolyConnecA);
              sPolyConnecB = sPolyConnecA;
              *polyClippConnecB = realloc (*polyClippConnecB, sizeof(PDM_g_num_t) *sPolyConnecB);
            }

            if ((3 * ((*polyClippIdxA)[*nPolyClippA+1])) >= sPolyCoordA) {
              while ((3 * ((*polyClippIdxA)[*nPolyClippA+1])) >= sPolyCoordA) {
                sPolyCoordA *= 2;
              }
              *polyClippCoordsA = realloc (*polyClippCoordsA, sizeof(double) *sPolyCoordA);
              sPolyCoordB = sPolyCoordA;
              *polyClippCoordsB = *polyClippCoordsA;

           }

            idx1 = (*polyClippIdxA)[*nPolyClippA];
            idx2 = 3 * ((*polyClippIdxA)[*nPolyClippA]);
            int isDegenerated = 0;

            for (int i = 0; i < n_vtxClipp; i++) {
              int iprev = _modulo((i-1),n_vtxClipp);
              int inext = _modulo((i+1),n_vtxClipp);

              if ((clipped_vtx[iprev] == clipped_vtx[inext])
                  || ((clipped_vtx[iprev]->neighbor != NULL)
                      && (clipped_vtx[inext]->neighbor != NULL)
                      && ((clipped_vtx[iprev]->neighbor == clipped_vtx[inext])
                          || (clipped_vtx[inext]->neighbor == clipped_vtx[iprev])))) {
                isDegenerated = 1;
                clipped_vtx[iprev]->used = true;
                clipped_vtx[i]->used = true;
                clipped_vtx[inext]->used = true;
                break;
              }

              if (origin_vtx[i] > 0) {
                if (clipped_vtx[i]->neighbor != NULL) {
                  (*polyClippConnecA)[idx1  ] = clipped_vtx[i]->gN;
                  (*polyClippConnecB)[idx1++] = clipped_vtx[i]->neighbor->gN;
                }
                else {
                  (*polyClippConnecA)[idx1  ] = clipped_vtx[i]->gN;
                  (*polyClippConnecB)[idx1++] = -clipped_vtx[i]->gN;
                }
              }

              else {
                if (clipped_vtx[i]->neighbor != NULL) {
                  (*polyClippConnecB)[idx1  ] = clipped_vtx[i]->gN;
                  (*polyClippConnecA)[idx1++] = clipped_vtx[i]->neighbor->gN;
                }
                else {
                  (*polyClippConnecB)[idx1  ] = clipped_vtx[i]->gN;
                  (*polyClippConnecA)[idx1++] = -clipped_vtx[i]->gN;
                }
              }

              if (clipped_vtx[i]->coords != NULL) {
                for (int j = 0; j < 3; j++) {
                  (*polyClippCoordsA)[idx2++] = clipped_vtx[i]->coords[j];
                }
              }
              else if ((clipped_vtx[i]->neighbor != NULL) ?
                       (clipped_vtx[i]->neighbor->coords != NULL)  : false) {
                for (int j = 0; j < 3; j++) {
                  (*polyClippCoordsA)[idx2++] = clipped_vtx[i]->neighbor->coords[j];
                }
              }
              else {
                printf ("PDM_poly_clipp error : Pas de coordonnees definies pour le point\n");
                abort();
              }
              clipped_vtx[i]->used = true;
            }

            if (!isDegenerated) {
              *nPolyClippA += 1;
              *nPolyClippB = *nPolyClippA;
            }

          }

          else {


            int *used_vtx = PDM_array_zeros_int(n_vtxClipp);

            int i = 0;
            int inext = -1;

            do {
              inext = -1;
              int ipass2 = 0;
              int inext1;
              int ibeg = i;

              idx1 = (*polyClippIdxA)[*nPolyClippA];
              idx2 = 3 *(*polyClippIdxA)[*nPolyClippA];

              do {

                if (PDM_ABS(origin_vtx[i]) == 100) {
                  inext1  = _modulo(link_multi[i]+1, n_vtxClipp);
                  if (ipass2 == 0) {
                    ipass2++;
                    inext = _modulo((i+1),n_vtxClipp);
                  }
                }
                else {
                  inext1 = _modulo((i+1),n_vtxClipp);
                }

                // Copy tab

                if (idx1 >= sPolyConnecA) {
                  sPolyConnecA *= 2;

                  *polyClippConnecA = realloc (*polyClippConnecA, sizeof(PDM_g_num_t) *sPolyConnecA);
                  sPolyConnecB = sPolyConnecA;
                  *polyClippConnecB = realloc (*polyClippConnecB, sizeof(PDM_g_num_t) *sPolyConnecB);
                }

                if (idx2+3 >= sPolyCoordA) {
                  sPolyCoordA *= 2;

                  *polyClippCoordsA = realloc (*polyClippCoordsA, sizeof(double) *sPolyCoordA);
                  sPolyCoordB = sPolyCoordA;
                  *polyClippCoordsB = *polyClippCoordsA;
                }

                if (origin_vtx[i] > 0) {
                  if (clipped_vtx[i]->neighbor != NULL) {
                    (*polyClippConnecA)[idx1  ] = clipped_vtx[i]->gN;
                    (*polyClippConnecB)[idx1++] = clipped_vtx[i]->neighbor->gN;
                  }
                  else {
                    (*polyClippConnecA)[idx1  ] = clipped_vtx[i]->gN;
                    (*polyClippConnecB)[idx1++] = -clipped_vtx[i]->gN;
                  }
                }

                else {
                  if (clipped_vtx[i]->neighbor != NULL) {
                    (*polyClippConnecB)[idx1  ] = clipped_vtx[i]->gN;
                    (*polyClippConnecA)[idx1++] = clipped_vtx[i]->neighbor->gN;
                  }
                  else {
                    (*polyClippConnecB)[idx1  ] = clipped_vtx[i]->gN;
                    (*polyClippConnecA)[idx1++] = -clipped_vtx[i]->gN;
                  }
                }

                for (int j = 0; j < 3; j++) {
                  (*polyClippCoordsA)[idx2++] = clipped_vtx[i]->coords[j];
                }

                used_vtx[i]++;
                if (PDM_ABS(origin_vtx[i]) == 100) {
                  used_vtx[link_multi[i]]++;
                }
                i = inext1;

              }  while ((PDM_ABS(origin_vtx[ibeg]) == 100) ?
                        ((inext1 != ibeg) || (inext1 !=link_multi[ibeg])) : (inext1 != ibeg));

              *nPolyClippA += 1;
              (*polyClippIdxA)[*nPolyClippA] = idx1;
              i = inext;

            } while ((inext == -1) ? false : used_vtx[inext] == 0);

            free (used_vtx);
          }

          // Verification

        }

        /* Look for next polygon sortie si on a fait le tour */

        if ((curr == first) || (curr->neighbor == first)) {
          break;
        }

        if (curr == first_poly) {
          curr = curr->next;
        }
        else if (curr->neighbor == first_poly) {
          curr = curr->neighbor->next;
        }
        else {
          printf ("Erreur dans le parcours du polygone");
          abort();
        }

      } while ((curr != first) && (curr->neighbor != first)); // On blinde le teste de sortie

      free (clipped_vtx);
      free (clipped_multi);
      free (origin_vtx);
      free (link_multi);

    }

    else {

      /*
       * Check if the constraint polygon is inside the subject polygon
       */

      bool inside = true;

      _vertex_poly_t *curr = vtxA;

      idx1 = 0;
      idx2 = 0;

      do {
        if (!curr->tag) {
          inside = false;
          break;
        }

        (*polyClippConnecA)[idx1  ] = curr->gN;
        (*polyClippConnecB)[idx1++] = -curr->gN;

        for (int j = 0; j < 3; j++) {
          (*polyClippCoordsA)[idx2++] = curr->coords[j];
        }

        curr = curr->next;

      } while (curr != vtxA);


      if (inside) {
        *nPolyClippA = 1;
        *nPolyClippB = 1;
        (*polyClippIdxA)[0] = 0;
        (*polyClippIdxA)[1] = n_vtxA;
      }

      /*
       * Check if the subject polygon is inside the constraint polygon
       */

      else {

        curr = vtxB;

        idx1 = 0;
        idx2 = 0;

        inside = true;

        do {
          if (!curr->tag) {
            inside = false;
            break;
          }

          (*polyClippConnecA)[idx1  ] = -curr->gN;
          (*polyClippConnecB)[idx1++] = curr->gN;

          for (int j = 0; j < 3; j++) {
            (*polyClippCoordsA)[idx2++] = curr->coords[j];
          }

          curr = curr->next;
        } while (curr != vtxB);


        if (inside) {
          *nPolyClippA = 1;
          *nPolyClippB = 1;
          (*polyClippIdxA)[1] = n_vtxB;
        }

      }

    }

    /*
     * Update size
     */

    *polyClippIdxA =
            realloc (*polyClippIdxA, (sizeof(int) * (*nPolyClippA + 1)));

    *polyClippConnecA =
            realloc (*polyClippConnecA, (sizeof(PDM_g_num_t) * (*polyClippIdxA)[*nPolyClippA]));

    *polyClippCoordsA =
            realloc (*polyClippCoordsA, (sizeof(double) * 3 * (*polyClippIdxA)[*nPolyClippA]));

    *polyClippIdxB = *polyClippIdxA;
    *polyClippConnecB =
            realloc (*polyClippConnecB, (sizeof(PDM_g_num_t) * (*polyClippIdxB)[*nPolyClippB]));

    *polyClippCoordsB = *polyClippCoordsA;

  }

  /*
   * Reverse mod
   */

  else if (performed_t == PDM_POLY_CLIPP_REVERSE) {

    *polyClippIdxA = malloc (sizeof(int) * (nPolyPredicA + 1));
    (*polyClippIdxA)[0] = 0;

    *polyClippIdxB = malloc (sizeof(int) * (nPolyPredicB + 1));
    (*polyClippIdxB)[0] = 0;

    *polyClippCoordsA = malloc (sizeof(double) * sPolyCoordA);
    *polyClippCoordsB = malloc (sizeof(double) * sPolyCoordB);

    /*
     * A Reverse clipping
     */

    /* Look for first out vertex*/

    _vertex_poly_t *first_out = vtxA;
    while (first_out->isect || first_out->tag) {
      first_out = first_out->next;
    }

    /* Fist polygon */

    idx1 = 0;

    if (!first_out->isect && !first_out->tag) {

      int s_clipped_vtx = n_vtxA + n_vtxB;

      _vertex_poly_t **clipped_vtx = malloc (sizeof(_vertex_poly_t*) * s_clipped_vtx);

      _vertex_poly_t *curr = first_out;

      do {

        /* Clip */

        int n_vtxClipp = 0;

        _vertex_poly_t *first_poly_out = curr;

        if (!curr->used) {

          curr->used = true;

          bool direction = true;

          do {

            if (curr->neighbor != NULL) {

              if (curr->isect) {

                direction = !direction;
                curr = curr->neighbor;

              }

              else {

                if ((direction && curr->neighbor->tag) ||
                    (!direction && !curr->neighbor->tag)) {
                  direction = !direction;
                  curr = curr->neighbor;
                }

              }

            }

            do {

              if (n_vtxClipp >= s_clipped_vtx) {
                while (n_vtxClipp >= s_clipped_vtx) {
                  s_clipped_vtx *= 2;
                }
                clipped_vtx =
                realloc (clipped_vtx, sizeof(_vertex_poly_t*) * s_clipped_vtx);
              }

              clipped_vtx[n_vtxClipp++] = curr;
              curr->used = true;

              curr = direction ? curr->next : curr->previous;

            } while ((curr != first_poly_out) && (curr->neighbor == NULL));

          } while (curr != first_poly_out);

          if (n_vtxClipp > 2) {

            if (*nPolyClippA >= nPolyPredicA) {
              while (*nPolyClippA >= nPolyPredicA) {
                nPolyPredicA *= 2;
              }
              *polyClippIdxA = realloc (*polyClippIdxA, sizeof(int) * (nPolyPredicA + 1));
            }

            *nPolyClippA += 1;
            (*polyClippIdxA)[*nPolyClippA] =
              (*polyClippIdxA)[*nPolyClippA - 1] + n_vtxClipp;

            if ((*polyClippIdxA)[*nPolyClippA] >= sPolyConnecA) {
              while ((*polyClippIdxA)[*nPolyClippA] >= sPolyConnecA) {
                sPolyConnecA *= 2;
              }
              *polyClippConnecA = realloc (*polyClippConnecA, sizeof(PDM_g_num_t) *sPolyConnecA);
            }

            if ((3 * ((*polyClippIdxA)[*nPolyClippA])) >= sPolyCoordA) {
              while ((3 * ((*polyClippIdxA)[*nPolyClippA])) >= sPolyCoordA) {
                sPolyCoordA *= 2;
              }
              *polyClippCoordsA = realloc (*polyClippCoordsA, sizeof(double) *sPolyCoordA);
            }

            for (int i = 0; i < n_vtxClipp; i++) {
              if (clipped_vtx[i]->first == vtxA) {
                (*polyClippConnecA)[idx1++] = clipped_vtx[i]->gN;
              }
              else {
                (*polyClippConnecA)[idx1++] = clipped_vtx[i]->neighbor->gN;
              }
              for (int j = 0; j < 3; j++) {
                (*polyClippCoordsA)[idx2++] = clipped_vtx[i]->coords[j];
              }
            }
          }
        }

        do {
          curr = curr->next;
        } while (curr->isect && (curr != first_out));

      } while (curr != first_out);

      free (clipped_vtx);

    }

    /*
     * B Reverse clipping
     */

    /* Look for first out vertex*/

    first_out = vtxB;
    while (first_out->isect || first_out->tag) {
      first_out = first_out->next;
    }

    idx1 = 0;
    idx2 = 0;

    /* Fist polygon */

    if (!first_out->isect && !first_out->tag) {

      int s_clipped_vtx = n_vtxA + n_vtxB;

      _vertex_poly_t **clipped_vtx = malloc (sizeof(_vertex_poly_t*) * s_clipped_vtx);

      _vertex_poly_t *curr = first_out;

      do {

        /* Clip */

        int n_vtxClipp = 0;

        _vertex_poly_t *first_poly_out = curr;

        if (!curr->used) {

          curr->used = true;

          bool direction = true;

          do {

            if (curr->neighbor != NULL) {

              if (curr->isect) {
                direction = !direction;
                curr = curr->neighbor;
              }

              else {
                if ((direction && curr->neighbor->tag) ||
                    (!direction && !curr->neighbor->tag)) {
                  direction = !direction;
                  curr = curr->neighbor;
                }
              }
            }

            do {

              if (n_vtxClipp >= s_clipped_vtx) {
                while (n_vtxClipp >= s_clipped_vtx) {
                  s_clipped_vtx *= 2;
                }
                clipped_vtx =
                realloc (clipped_vtx, sizeof(_vertex_poly_t*) * s_clipped_vtx);
              }

              clipped_vtx[n_vtxClipp++] = curr;
              curr->used = true;

              curr = direction ? curr->next : curr->previous;

            } while ((curr != first_poly_out) && (curr->neighbor == NULL));

          } while (curr != first_poly_out);

          if (n_vtxClipp > 2) {

            if (*nPolyClippB >= nPolyPredicB) {
              while (*nPolyClippB >= nPolyPredicB) {
                nPolyPredicB *= 2;
              }
              *polyClippIdxB = realloc (*polyClippIdxB, sizeof(int) * (nPolyPredicB + 1));
            }

            *nPolyClippB += 1;
            (*polyClippIdxB)[*nPolyClippB] =
              (*polyClippIdxB)[*nPolyClippB - 1] + n_vtxClipp;

            if ((*polyClippIdxB)[*nPolyClippB] >= sPolyConnecB) {
              while ((*polyClippIdxB)[*nPolyClippB] >= sPolyConnecB) {
                sPolyConnecB *= 2;
              }
              *polyClippConnecB = realloc (*polyClippConnecB, sizeof(PDM_g_num_t) *sPolyConnecB);
            }

            if ((3 * ((*polyClippIdxB)[*nPolyClippB])) >= sPolyCoordB) {
              while ((3 * ((*polyClippIdxB)[*nPolyClippB])) >= sPolyCoordB) {
                sPolyCoordB *= 2;
              }
              *polyClippCoordsB = realloc (*polyClippCoordsB, sizeof(double) *sPolyCoordB);
            }

            for (int i = 0; i < n_vtxClipp; i++) {
              if (clipped_vtx[i]->first == vtxB) {
                (*polyClippConnecB)[idx1++] = clipped_vtx[i]->gN;
              }
              else {
                (*polyClippConnecB)[idx1++] = clipped_vtx[i]->neighbor->gN;
              }
              for (int j = 0; j < 3; j++) {
                (*polyClippCoordsB)[idx2++] = clipped_vtx[i]->coords[j];
              }
            }
          }
        }

        do {
          curr = curr->next;
        } while (curr->isect && (curr != first_out));

      } while (curr != first_out);

      free (clipped_vtx);
    }

    /*
     * Update size
     */

    *polyClippIdxA =
            realloc (*polyClippIdxA, (sizeof(int) * (*nPolyClippA + 1)));
    *polyClippConnecA =
            realloc (*polyClippConnecA, (sizeof(PDM_g_num_t) * (*polyClippIdxA)[*nPolyClippA]));
    *polyClippCoordsA =
            realloc (*polyClippCoordsA, (sizeof(double)* 3 * (*polyClippIdxA)[*nPolyClippA]));

    *polyClippIdxB =
            realloc (*polyClippIdxB, (sizeof(int) * (*nPolyClippB + 1)));
    *polyClippConnecB =
            realloc (*polyClippConnecB, (sizeof(PDM_g_num_t) * (*polyClippIdxB)[*nPolyClippB]));

    *polyClippCoordsB =
            realloc (*polyClippCoordsB, (sizeof(double)* 3 * (*polyClippIdxB)[*nPolyClippB]));
  }

  if (1 == 0) {
    printf ("Sous-polygones de A : %d\n", *nPolyClippA);
    for (int i1 = 0; i1 < *nPolyClippA; i1++) {
      printf ("%d :",i1);
      for (int j1 = (*polyClippIdxA)[i1]; j1 < (*polyClippIdxA)[i1+1]; j1++) {
        printf (" "PDM_FMT_G_NUM, (*polyClippConnecA)[j1]);
        printf ("=(%.2f",(*polyClippCoordsA)[3*j1]);
        printf (",%.2f)", (*polyClippCoordsA)[3*j1+1]);
        printf ("%d",j1);
      }
      printf("\n");
    }

    printf ("Sous-polygones de B : %d\n", *nPolyClippB);
    for (int i1 = 0; i1 < *nPolyClippB; i1++) {
      printf ("%d :",i1);
      for (int j1 = (*polyClippIdxB)[i1]; j1 < (*polyClippIdxB)[i1+1]; j1++) {
        printf (" "PDM_FMT_G_NUM,   (*polyClippConnecB)[j1]);
        printf ("=(%.2f", (*polyClippCoordsB)[3*j1]);
        printf (",%.2f)", (*polyClippCoordsB)[3*j1+1]);
        printf ("%d", j1);
      }
      printf("\n");
    }
  }

  /*
   * Free local memory
   */

  _poly_clipp_free (vtxA);
  _poly_clipp_free (vtxB);;

  if (_faceToEdgeA != faceToEdgeA) {
    free (_faceToEdgeA);
  }

  if (_faceToVtxA  != faceToVtxA) {
    free (_faceToVtxA);
  }

  if (_face_vtxCooA != face_vtxCooA) {
    free (_face_vtxCooA);
  }

  if (_faceToEdgeB != faceToEdgeB) {
    free (_faceToEdgeB);
  }

  if (_faceToVtxB != faceToVtxB) {
    free (_faceToVtxB);
  }

  if (_face_vtxCooB != face_vtxCooB) {
    free (_face_vtxCooB);
  }

}

void
PDM_vertex_poly_dump(
_vertex_poly_t  *v,
int  verbosity)
{
	printf (" v:%p:"PDM_FMT_G_NUM"  v->first:%p v->next:%p:"PDM_FMT_G_NUM" v->previous:%p:"PDM_FMT_G_NUM"\n", (void *)v, v->gN, (void *)v->first, (void *)v->next, v->next->gN, (void *)v->previous, v->previous->gN);
	_vertex_poly_t *vtx_curr = v;
    do {
      PDM_g_num_t neighbour = -1;
      if (vtx_curr->neighbor != NULL)
        neighbour = vtx_curr->neighbor->gN;
      if ( neighbour != -1)
	  {
        printf (" "PDM_FMT_G_NUM"/"PDM_FMT_G_NUM, vtx_curr->gN, neighbour);
		if (verbosity>1)
			printf ("/%d/%d", vtx_curr->isect, vtx_curr->tag);
	  }
      else
	  {
        printf (" "PDM_FMT_G_NUM"/NC", vtx_curr->gN);
		if (verbosity>1)
			printf ("/%d/%d", vtx_curr->isect, vtx_curr->tag);
	  }
      vtx_curr = vtx_curr->next;
    } while (vtx_curr != v->first);
    printf("\n");
}
#ifdef	__cplusplus
}
#endif
