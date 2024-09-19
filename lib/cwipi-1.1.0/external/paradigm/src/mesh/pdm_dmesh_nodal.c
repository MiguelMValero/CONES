
/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_mpi.h"
#include "pdm_dmesh_nodal.h"
#include "pdm_dmesh_nodal_priv.h"
#include "pdm_dmesh_nodal_elmts.h"
#include "pdm_part_to_block.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"
#include "pdm_quick_sort.h"
#include "pdm_geom_elem.h"
#include "pdm_dmesh_nodal_elements_utils.h"
#include "pdm_para_graph_dual.h"
#include "pdm_block_to_part.h"
#include "pdm_dconnectivity_transform.h"
#include "pdm_partitioning_algorithm.h"
#include "pdm_partitioning_nodal_algorithm.h"
#include "pdm_vtk.h"
#include "pdm_ho_ordering.h"
#include "pdm_array.h"
#include "pdm_distrib.h"

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

/*----------------------------------------------------------------------------
 * Maximum number of sections depending of section type
 *----------------------------------------------------------------------------*/

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Free vtx structure
 *
 * \param[inout]  vtx    Vertices
 *
 * \return        NULL
 *
 */

static
PDM_DMesh_nodal_vtx_t *
_vtx_free
(
 PDM_DMesh_nodal_vtx_t *vtx
)
{
  if (vtx != NULL) {
    if (vtx->distrib != NULL) {
      free (vtx->distrib);
      vtx->distrib = NULL;
    }

    if(vtx->owner == PDM_OWNERSHIP_KEEP) {
      if(vtx->_coords != NULL) {
        free(vtx->_coords);
        vtx->_coords = NULL;
      }

      if(vtx->dvtx_tag != NULL) {
        free(vtx->dvtx_tag);
        vtx->dvtx_tag = NULL;
      }

      if(vtx->dvtx_parent_g_num != NULL) {
        free(vtx->dvtx_parent_g_num);
        vtx->dvtx_parent_g_num = NULL;
      }
    }
    free (vtx);
  }
  return NULL;
}

static
PDM_dmesh_nodal_elmts_t*
_get_from_geometry_kind
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  PDM_dmesh_nodal_elmts_t* dmne = NULL;
  if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC){
    assert(dmesh_nodal->mesh_dimension == 3);
    dmne = dmesh_nodal->volumic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC){
    assert(dmesh_nodal->mesh_dimension >= 2);
    dmne = dmesh_nodal->surfacic;
  } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE){
    assert(dmesh_nodal->mesh_dimension >= 1);
    dmne = dmesh_nodal->ridge;
  } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER){
    dmne = dmesh_nodal->corner;
  } else {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom_kind in _get_from_geometry_kind \n");
  }
  return dmne;
}

/**
 *
 * \brief Initialize a mesh
 *
 * \param [inout]  mesh        Mesh
 * \param [in]     n_part      Number of partitions
 */

static void
_mesh_init
(
      PDM_dmesh_nodal_t  *dmesh_nodal,
const PDM_MPI_Comm        comm,
      int                 mesh_dimension,
      PDM_g_num_t         n_vtx,
      PDM_g_num_t         n_cell,
      PDM_g_num_t         n_face,
      PDM_g_num_t         n_edge
)
{
  int n_rank;
  int i_rank;

  PDM_MPI_Comm_size (comm, &n_rank);
  PDM_MPI_Comm_rank (comm, &i_rank);

  dmesh_nodal->comm                     = comm;
  dmesh_nodal->n_rank                   = n_rank;
  dmesh_nodal->i_rank                   = i_rank;

  dmesh_nodal->mesh_dimension           = mesh_dimension;
  dmesh_nodal->n_cell_abs               = n_cell;
  dmesh_nodal->n_face_abs               = n_face;
  dmesh_nodal->n_edge_abs               = n_edge;
  dmesh_nodal->n_vtx_abs                = n_vtx;

  dmesh_nodal->vtx                      = malloc(sizeof(PDM_DMesh_nodal_vtx_t ));
  dmesh_nodal->vtx->_coords             = NULL;
  dmesh_nodal->vtx->distrib             = NULL;
  dmesh_nodal->vtx->n_vtx               = 0;
  dmesh_nodal->vtx->dvtx_tag            = NULL;
  dmesh_nodal->vtx->dvtx_parent_g_num   = NULL;

  dmesh_nodal->volumic  = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 3, dmesh_nodal->n_cell_abs);
  dmesh_nodal->surfacic = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 2, dmesh_nodal->n_face_abs);
  dmesh_nodal->ridge    = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 1, dmesh_nodal->n_edge_abs);
  dmesh_nodal->corner   = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 0, dmesh_nodal->n_vtx_abs );

}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 * \brief Create a Mesh nodal structure
 *
 * \param [in]   n_part   Number of partition on the current process
 * \param [in]   comm     MPI communicator
 *
 * \return       New mesh nodal handle
 *
 */

PDM_dmesh_nodal_t*
PDM_DMesh_nodal_create
(
const PDM_MPI_Comm comm,
      int          mesh_dimension,
      PDM_g_num_t  n_vtx,
      PDM_g_num_t  n_cell,
      PDM_g_num_t  n_face,
      PDM_g_num_t  n_edge
)
{
  PDM_dmesh_nodal_t *mesh = (PDM_dmesh_nodal_t *) malloc (sizeof(PDM_dmesh_nodal_t));

  _mesh_init (mesh, comm, mesh_dimension, n_vtx, n_cell, n_face, n_edge);

  mesh->is_computed_g_extents = PDM_FALSE;

  return mesh;
}


/**
 * \brief Free a nodal mesh structure
 *
 * \param [in]  idx      Nodal mesh handle
 * \param [in]  partial  if partial is equal to 0, all data are removed.
 *                       Otherwise, results are kept.
 *
 * \return      NULL
 *
 */

void
PDM_DMesh_nodal_free
(
 PDM_dmesh_nodal_t *dmesh_nodal
)
{

  if (dmesh_nodal != NULL) {

    _vtx_free(dmesh_nodal->vtx);

    if(dmesh_nodal->volumic != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->volumic);
    }
    if(dmesh_nodal->surfacic != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->surfacic);
    }
    if(dmesh_nodal->ridge != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->ridge);
    }
    if(dmesh_nodal->corner != NULL){
      PDM_DMesh_nodal_elmts_free(dmesh_nodal->corner);
    }

    free(dmesh_nodal);
  }
}


/**
 * \brief Define partition vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 * \param [in]  n_vtx     Number of vertices
 * \param [in]  coords    Interlaced coordinates (size = 3 * \ref n_vtx)
 *
 */

void
PDM_DMesh_nodal_coord_set
(
       PDM_dmesh_nodal_t *dmesh_nodal,
 const int                n_vtx,
       PDM_real_t        *coords,
       PDM_ownership_t    owner
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  if (vtx->_coords != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Vertices are already defined\n");
  }

  /* Mapping memoire */

  vtx->n_vtx   = n_vtx;
  vtx->_coords = coords;
  vtx->owner   = owner;

  vtx->distrib = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (dmesh_nodal->n_rank + 1));

  PDM_g_num_t _n_vtx = n_vtx;
  PDM_MPI_Allgather((void *) &_n_vtx,
                    1,
                    PDM__PDM_MPI_G_NUM,
                    (void *) (&vtx->distrib[1]),
                    1,
                    PDM__PDM_MPI_G_NUM,
                    dmesh_nodal->comm);

  vtx->distrib[0] = 0;
  for (int i = 1; i < dmesh_nodal->n_rank + 1; i++) {
    vtx->distrib[i] +=  vtx->distrib[i-1];
  }
}


void
PDM_DMesh_nodal_vtx_tag_set
(
 PDM_dmesh_nodal_t *dmesh_nodal,
 int               *dvtx_tag
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  if (vtx->dvtx_tag != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "dvtx_tag are already defined\n");
  }

  vtx->dvtx_tag = dvtx_tag;
}

void
PDM_DMesh_nodal_vtx_parent_gnum_set
(
 PDM_dmesh_nodal_t *dmesh_nodal,
 PDM_g_num_t       *dvtx_parent_g_num
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  if (vtx->dvtx_parent_g_num != NULL) {
    PDM_error(__FILE__, __LINE__, 0, "dvtx_tag are already defined\n");
  }

  vtx->dvtx_parent_g_num = dvtx_parent_g_num;
}



void
PDM_DMesh_nodal_section_g_dims_get
(
  PDM_dmesh_nodal_t *dmesh_nodal,
  PDM_g_num_t       *n_cell_abs,
  PDM_g_num_t       *n_face_abs,
  PDM_g_num_t       *n_edge_abs,
  PDM_g_num_t       *n_vtx_abs
)
{
  *n_cell_abs = dmesh_nodal->n_cell_abs;
  *n_face_abs = dmesh_nodal->n_face_abs;
  *n_edge_abs = dmesh_nodal->n_edge_abs;
  *n_vtx_abs  = dmesh_nodal->n_vtx_abs;
}

/**
 * \brief  Return number of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Number of vertices
 *
 */

int
PDM_DMesh_nodal_n_vtx_get
(
  PDM_dmesh_nodal_t *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->n_vtx;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

double *
PDM_DMesh_nodal_vtx_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->_coords;
}

/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

int*
PDM_DMesh_nodal_vtx_tag_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->dvtx_tag;
}


/**
 * \brief  Return coordinates of vertices
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Coordinates of vertices
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_vtx_parent_gnum_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{

  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->dvtx_parent_g_num;
}

/**
 * \brief  Return number of sections
 *
 * \param [in]  hdl            Nodal mesh handle
 *
 * \return  Number of sections
 *
 */

int
PDM_DMesh_nodal_n_section_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  if(dmne){
    return dmne->n_section;
  } else {
    return 0;
  }
}

/**
 * \brief  Return sections identifier
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 *
 * \return  Blocks identifier
 *
 */
int *
PDM_DMesh_nodal_sections_id_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return dmne->sections_id;
}

/**
 * \brief  Return type of element of section
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  Type of section
 *
 */
PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_elt_type_get
(
        PDM_dmesh_nodal_t   *dmesh_nodal,
        PDM_geometry_kind_t  geom_kind,
  const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  if (id_section <= PDM_BLOCK_ID_BLOCK_POLY2D) { // std
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    return dmne->sections_std[_id_section]->t_elt;
  }
  assert(0); // only useful for std elements
  return (PDM_Mesh_nodal_elt_t) -1;
}

PDM_Mesh_nodal_elt_t
PDM_DMesh_nodal_section_type_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);

  PDM_Mesh_nodal_elt_t t_elt = PDM_MESH_NODAL_POLY_3D;

  if (id_section < PDM_BLOCK_ID_BLOCK_POLY2D) {

    t_elt = PDM_MESH_NODAL_POLY_3D;
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    PDM_DMesh_nodal_section_std_t *section = dmne->sections_std[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad section identifier\n");
    }

    t_elt = section->t_elt;
  }

  else if (id_section < PDM_BLOCK_ID_BLOCK_POLY3D) {

    t_elt = PDM_MESH_NODAL_POLY_2D;

  }

  else {

    t_elt = PDM_MESH_NODAL_POLY_3D;

  }

  return t_elt;
}

/**
 * \brief  Return distri of section
 *
 * \param [in] dmesh_nodal
 * \param [in] id_section   Block identifier
 *
 * \return  distri
 *
 */
PDM_g_num_t*
PDM_DMesh_nodal_section_distri_std_get
(
        PDM_dmesh_nodal_t   *dmesh_nodal,
        PDM_geometry_kind_t  geom_kind,
  const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  if (id_section <= PDM_BLOCK_ID_BLOCK_POLY2D) { // std
    int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;
    return dmne->sections_std[_id_section]->distrib;
  }
  assert(0); // only useful for std elements
  return NULL;
}

int
PDM_DMesh_nodal_section_add
(
      PDM_dmesh_nodal_t    *dmesh_nodal,
      PDM_geometry_kind_t   geom_kind,
const PDM_Mesh_nodal_elt_t  t_elt
)
{
  if( _get_from_geometry_kind(dmesh_nodal, geom_kind) == NULL) {
    if(geom_kind == PDM_GEOMETRY_KIND_VOLUMIC) {
      dmesh_nodal->volumic = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 3, dmesh_nodal->n_cell_abs);
    } else if( geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
      dmesh_nodal->surfacic = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 2, dmesh_nodal->n_face_abs);
    } else if( geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
      dmesh_nodal->ridge = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 1, dmesh_nodal->n_edge_abs);
    } else if( geom_kind == PDM_GEOMETRY_KIND_CORNER) {
      dmesh_nodal->corner = PDM_DMesh_nodal_elmts_create(dmesh_nodal->comm, 0, dmesh_nodal->n_vtx_abs);
    }
  }

  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  return PDM_DMesh_nodal_elmts_section_add(dmne, t_elt);
}

void
PDM_DMesh_nodal_update_ownership
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_ownership_t      owner
)
{
  dmesh_nodal->vtx->owner = owner;
  if (dmesh_nodal->volumic != NULL)
    PDM_DMesh_nodal_elmts_update_ownership(dmesh_nodal->volumic, owner);
  if (dmesh_nodal->surfacic != NULL)
    PDM_DMesh_nodal_elmts_update_ownership(dmesh_nodal->surfacic, owner);
  if (dmesh_nodal->ridge != NULL)
    PDM_DMesh_nodal_elmts_update_ownership(dmesh_nodal->ridge, owner);
  if (dmesh_nodal->corner != NULL)
    PDM_DMesh_nodal_elmts_update_ownership(dmesh_nodal->corner, owner);
}


void
PDM_DMesh_nodal_section_std_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const int                  n_elt,
      PDM_g_num_t         *connec,
      PDM_ownership_t      owner
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_section_std_set(dmne, id_section, n_elt, connec, owner);
}

void
PDM_DMesh_nodal_section_group_elmt_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  n_group_elmt,
      int                 *dgroup_elmt_idx,
      PDM_g_num_t         *dgroup_elmt,
      PDM_ownership_t      owner
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_group_set(dmne, n_group_elmt, dgroup_elmt_idx, dgroup_elmt, owner);
}

void
PDM_DMesh_nodal_section_group_elmt_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind,
 int                 *n_group_elmt,
 int                 **dgroup_elmt_idx,
 PDM_g_num_t         **dgroup_elmt
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  *n_group_elmt    = dmne->n_group_elmt;
  *dgroup_elmt_idx = dmne->dgroup_elmt_idx;
  *dgroup_elmt     = dmne->dgroup_elmt;
}

/**
 * \brief Return standard section description
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return  connect           Connectivity
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_section_std_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  return PDM_DMesh_nodal_elmts_section_std_get(dmne, id_section);
}

PDM_g_num_t *
PDM_DMesh_nodal_section_std_ho_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
      int                 *order,
const char               **ho_ordering
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);

  PDM_g_num_t *delt_vtx = PDM_DMesh_nodal_elmts_section_std_ho_get(dmne, id_section, order, ho_ordering);

  return delt_vtx;
}
/**
 * \brief Return standard section description
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return  connect           Connectivity
 *
 */

PDM_g_num_t *
PDM_DMesh_nodal_elmts_section_std_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
)
{
  if (dmn_elts == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  return section->_connec;
}


PDM_g_num_t *
PDM_DMesh_nodal_elmts_section_std_ho_get
(
      PDM_dmesh_nodal_elmts_t  *dmn_elts,
const int                       id_section,
      int                      *order,
const char                    **ho_ordering
)
{
  if (dmn_elts == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  int _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

  PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

  *order       = dmn_elts->sections_std[_id_section]->order;
  *ho_ordering = dmn_elts->sections_std[_id_section]->ho_ordering;

  if (section == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
  }

  return section->_connec;
}

/**
 * \brief Get number of section elements
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return      Number of elements
 *
 */

int
PDM_DMesh_nodal_section_n_elt_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return PDM_DMesh_nodal_elmts_section_n_elt_get(dmne, id_section);
}


/**
 * \brief Get number of section elements
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 *
 * \return      Number of elements
 *
 */

int
PDM_DMesh_nodal_elmts_section_n_elt_get
(
      PDM_dmesh_nodal_elmts_t *dmn_elts,
const int                      id_section
)
{
  int _id_section;

  if (id_section >= PDM_BLOCK_ID_BLOCK_POLY3D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY3D;

    PDM_DMesh_nodal_section_poly3d_t *section = dmn_elts->sections_poly3d[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad standard section identifier\n");
    }

    return section->n_elt;
  }

  else if (id_section >= PDM_BLOCK_ID_BLOCK_POLY2D) {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_POLY2D;

    PDM_DMesh_nodal_section_poly2d_t *section = dmn_elts->sections_poly2d[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polygon section identifier\n");
    }

    return section->n_elt;
  }

  else {

    _id_section = id_section - PDM_BLOCK_ID_BLOCK_STD;

    PDM_DMesh_nodal_section_std_t *section = dmn_elts->sections_std[_id_section];

    if (section == NULL) {
      PDM_error (__FILE__, __LINE__, 0, "Bad polyhedron section identifier\n");
    }

    return section->n_elt;
  }

}


/**
 * \brief Define a polygon section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of elements
 * \param [in]  connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [in]  connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_section_poly2d_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const PDM_l_num_t          n_elt,
      PDM_l_num_t         *connec_idx,
      PDM_g_num_t         *connec,
      PDM_ownership_t      owner
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_section_poly2d_set(dmne, id_section, n_elt, connec_idx, connec, owner);
}


/**
 * \brief Return a polygon section description
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out] connect_idx    Connectivity index (size = \ref n_elt + 1)
 * \param [out] connect        Connectivity (size = \ref connect_idx[\ref n_elt])
 *
 */

void
PDM_DMesh_nodal_section_poly2d_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
      PDM_l_num_t        **connec_idx,
      PDM_g_num_t        **connec
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  assert(dmne != NULL);
  PDM_DMesh_nodal_elmts_section_poly2d_get(dmne, id_section, connec_idx, connec);
}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [in]  n_elt          Number of polyhedra
 * \param [in]  n_face         Number of faces used to describe polyhedra
 * \param [in]  facvtx_idx     Index of face vertex connectivity
 * \param [in]  facvtx         Face vertex connectivity
 * \param [in]  cellfac_idx    Index of cell face connectivity
 * \param [in]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_set
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section,
const PDM_l_num_t          n_elt,
const PDM_l_num_t          n_face,
      PDM_l_num_t         *facvtx_idx,
      PDM_g_num_t         *facvtx,
      PDM_l_num_t         *cellfac_idx,
      PDM_g_num_t         *cellfac,
      PDM_ownership_t      owner
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  PDM_DMesh_nodal_elmts_section_poly3d_set(dmne, id_section, n_elt, n_face, facvtx_idx, facvtx, cellfac_idx, cellfac, owner);
}


/**
 * \brief Define a polyhedra section
 *
 * \param [in]  hdl            Distributed nodal mesh handle
 * \param [in]  id_section       Block identifier
 * \param [out]  n_face         Number of faces used to describe polyhedra
 * \param [out]  facvtx_idx     Index of face vertex connectivity
 * \param [out]  facvtx         Face vertex connectivity
 * \param [out]  cellfac_idx    Index of cell face connectivity
 * \param [out]  cellfac        Cell face connectivity
 *
 */

void
PDM_DMesh_nodal_section_poly3d_get
(
      PDM_dmesh_nodal_t    *dmesh_nodal,
      PDM_geometry_kind_t   geom_kind,
const int                   id_section,
      PDM_l_num_t          *n_face,
      PDM_l_num_t         **facvtx_idx,
      PDM_g_num_t         **facvtx,
      PDM_l_num_t         **cellfac_idx,
      PDM_g_num_t         **cellfac
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  PDM_DMesh_nodal_elmts_section_poly3d_get(dmne, id_section, n_face, facvtx_idx, facvtx, cellfac_idx, cellfac);
}


/**
 * \brief  Return total number of elements of a distributed mesh
 *
 * \param [in]  hdl       Distributed nodal mesh handle
 *
 * \return  Return number elements of a partition
 *
 */

PDM_g_num_t
PDM_dmesh_nodal_total_n_elmt_get
(
 PDM_dmesh_nodal_t   *dmesh_nodal,
 PDM_geometry_kind_t  geom_kind
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return PDM_DMesh_nodal_elmts_total_n_elmt_get(dmne);
}

/**
 * \brief  Return vtx distribution of a distributed mesh
 *
 * \param [in]  dmesh_nodal
 *
 * \return  Return vtx distribution
 *
 */
PDM_g_num_t*
PDM_dmesh_nodal_vtx_distrib_get
(
  PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  PDM_dmesh_nodal_t* mesh = (PDM_dmesh_nodal_t *) dmesh_nodal;
  return mesh->vtx->distrib;
}

/**
 * \brief  Return vtx distribution of a distributed mesh
 *
 * \param [in]  dmesh_nodal
 *
 * \return  Return vtx distribution
 *
 */
PDM_g_num_t*
PDM_dmesh_nodal_vtx_distrib_copy_get
(
  PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  PDM_dmesh_nodal_t* mesh = (PDM_dmesh_nodal_t *) dmesh_nodal;

  PDM_g_num_t* distrib = (PDM_g_num_t *) malloc( (dmesh_nodal->n_rank + 1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < dmesh_nodal->n_rank + 1; ++i) {
    distrib[i] = mesh->vtx->distrib[i];
  }

  return distrib;
}



PDM_g_num_t
PDM_dmesh_nodal_total_n_vtx_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->distrib[dmesh_nodal->n_rank];
}

/**
 *
 * \brief Setup global distribution of all elements register in current structure
 *
 * \param [inout]  mesh
 *
 * \return         Null
 *
 */
void
PDM_dmesh_nodal_generate_distribution
(
 PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  if(dmesh_nodal->volumic != NULL) {
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->volumic );
  }
  if(dmesh_nodal->surfacic != NULL){
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->surfacic);
  }
  if(dmesh_nodal->ridge != NULL){
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->ridge);
  }
  if(dmesh_nodal->corner != NULL){
     PDM_dmesh_nodal_elmts_generate_distribution(dmesh_nodal->corner);
  }
}



/**
*
* \brief Concatenates all element sections blocks
*
* \param [in]   dmesh_nodal
* \param [out]  section_idx        index of element section
* \param [out]  cat_delt_vtx_idx   index of element
* \param [out]  cat_delt_vtx       element vtx
*
 * \return     Number sections
*/
int PDM_concat_elt_sections
(
  PDM_dmesh_nodal_t  *dmesh_nodal,
  int               **section_idx,
  int               **cat_delt_vtx_idx,
  PDM_g_num_t       **cat_delt_vtx
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(section_idx);
  PDM_UNUSED(cat_delt_vtx_idx);
  PDM_UNUSED(cat_delt_vtx);
  abort(); // TODO Remove
  // 0. sizes
  int n_section = 0;
  // int n_section = PDM_DMesh_nodal_n_section_get(dmesh_nodal);
  // int dn_elt_vtx = 0;
  // *section_idx = (int*) malloc((n_section+1) * sizeof(int));
  // int* _section_idx = *section_idx;
  // int* n_vtx_by_elt_by_section = (int*) malloc(n_section * sizeof(int));
  // _section_idx[0] = 0;
  // int* n_elt_vtx_by_section = (int*) malloc(n_section * sizeof(int));
  // for (int i=0; i<n_section; ++i) {
  //   int n_elt_by_section = PDM_DMesh_nodal_section_n_elt_get(dmesh_nodal,i);
  //   _section_idx[i+1] = _section_idx[i] + n_elt_by_section;
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   n_vtx_by_elt_by_section[i] = PDM_Mesh_nodal_n_vertices_element(type,1); // 1: elements of order 1
  //   n_elt_vtx_by_section[i] = n_elt_by_section*n_vtx_by_elt_by_section[i];
  //   dn_elt_vtx += n_elt_vtx_by_section[i];
  // }

  // // 1. cat_delt_vtx_idx
  // int dn_elt = _section_idx[n_section];
  // *cat_delt_vtx_idx = (int*) malloc((dn_elt+1)* sizeof(int));
  // int* _cat_delt_vtx_idx = *cat_delt_vtx_idx;
  // _cat_delt_vtx_idx[0] = 0;
  // int pos_idx = 1;
  // for (int i=0; i<n_section; ++i) {
  //   int n_elt_by_section = _section_idx[i+1] - _section_idx[i];
  //   for (int j=0; j<n_elt_by_section; ++j) {
  //     _cat_delt_vtx_idx[pos_idx+j] = n_vtx_by_elt_by_section[i];
  //   }
  //   pos_idx += n_elt_by_section;
  // }
  // for (int i=1; i<dn_elt+1; ++i) {
  //   _cat_delt_vtx_idx[i] += _cat_delt_vtx_idx[i-1];
  // }

  // // 2. cat_delt_vtx
  // *cat_delt_vtx = (PDM_g_num_t *) malloc(dn_elt_vtx * sizeof(PDM_g_num_t));
  // PDM_g_num_t* _cat_delt_vtx = *cat_delt_vtx;
  // int pos = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_g_num_t* delt_vtx = PDM_DMesh_nodal_section_std_get(dmesh_nodal,i);
  //   for (int j=0; j<n_elt_vtx_by_section[i]; ++j) {
  //     _cat_delt_vtx[pos+j] = delt_vtx[j];
  //   }
  //   pos += n_elt_vtx_by_section[i];
  // }

  // // 3. free
  // free(n_elt_vtx_by_section);
  // free(n_vtx_by_elt_by_section);

  return n_section;
}


/**
*
* \brief Concatenates 3D element sections blocks
*
* \param [in]   dmesh_nodal
* \param [out]  section_idx        index of element section
* \param [out]  cat_delt_vtx_idx   index of element
* \param [out]  cat_delt_vtx       element vtx
*
 * \return     Number sections
*/
int PDM_concat_cell_sections
(
  PDM_dmesh_nodal_t  *dmesh_nodal,
  int               **cell_section_idx,
  int               **cat_dcell_vtx_idx,
  PDM_g_num_t       **cat_dcell_vtx
)
{
  PDM_UNUSED(dmesh_nodal);
  PDM_UNUSED(cell_section_idx);
  PDM_UNUSED(cat_dcell_vtx_idx);
  PDM_UNUSED(cat_dcell_vtx);
  int n_cell_section = -1;
  abort();
  // 0. sizes
  // int n_section = PDM_DMesh_nodal_n_section_get(dmesh_nodal);
  // int n_cell_section = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   if (PDM_Mesh_nodal_is_3D_element(type)) {
  //     ++n_cell_section;
  //   }
  // }

  // int dn_cell_vtx = 0;
  // *cell_section_idx = (int*) malloc((n_cell_section+1) * sizeof(int));
  // int* _cell_section_idx = *cell_section_idx;
  // int* n_vtx_by_cell_by_section = (int*) malloc(n_cell_section * sizeof(int));
  // _cell_section_idx[0] = 0;
  // int* n_cell_vtx_by_section = (int*) malloc(n_section * sizeof(int));
  // int cnt = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   if (PDM_Mesh_nodal_is_3D_element(type)) {
  //     int n_cell_by_section = PDM_DMesh_nodal_section_n_elt_get(dmesh_nodal,i);
  //     _cell_section_idx[cnt+1] = _cell_section_idx[cnt] + n_cell_by_section;
  //     n_vtx_by_cell_by_section[cnt] = PDM_Mesh_nodal_n_vertices_element(type,1); // 1: elements of order 1
  //     n_cell_vtx_by_section[i] = n_cell_by_section*n_vtx_by_cell_by_section[cnt];
  //     dn_cell_vtx += n_cell_vtx_by_section[i];
  //     ++cnt;
  //   }
  // }

  // // 1. cat_delt_vtx_idx
  // int dn_cell = _cell_section_idx[n_cell_section];
  // *cat_dcell_vtx_idx = (int*) malloc((dn_cell+1)* sizeof(int));
  // int* _cat_dcell_vtx_idx = *cat_dcell_vtx_idx;
  // _cat_dcell_vtx_idx[0] = 0;
  // int pos_idx = 1;
  // for (int i=0; i<n_cell_section; ++i) {
  //   int n_cell_by_section = _cell_section_idx[i+1] - _cell_section_idx[i];
  //   for (int j=0; j<n_cell_by_section; ++j) {
  //     _cat_dcell_vtx_idx[pos_idx+j] = n_vtx_by_cell_by_section[i];
  //   }
  //   pos_idx += n_cell_by_section;
  // }
  // for (int i=1; i<dn_cell+1; ++i) {
  //   _cat_dcell_vtx_idx[i] += _cat_dcell_vtx_idx[i-1];
  // }

  // // 2. cat_delt_vtx
  // *cat_dcell_vtx = (PDM_g_num_t *) malloc(dn_cell_vtx * sizeof(PDM_g_num_t));
  // PDM_g_num_t* _cat_dcell_vtx = *cat_dcell_vtx;
  // int pos = 0;
  // for (int i=0; i<n_section; ++i) {
  //   PDM_Mesh_nodal_elt_t type = PDM_DMesh_nodal_section_elt_type_get(dmesh_nodal,i);
  //   if (PDM_Mesh_nodal_is_3D_element(type)) {
  //     PDM_g_num_t* dcell_vtx = PDM_DMesh_nodal_section_std_get(dmesh_nodal,i);
  //     for (int j=0; j<n_cell_vtx_by_section[i]; ++j) {
  //       _cat_dcell_vtx[pos+j] = dcell_vtx[j];
  //     }
  //     pos += n_cell_vtx_by_section[i];
  //   }
  // }

  // // 3. free
  // free(n_cell_vtx_by_section);
  // free(n_vtx_by_cell_by_section);

  return n_cell_section;
}


/**
 * \brief  Compute elt->elt connectivity
 *
 * \param [out]  dual_graph_idx
 * \param [out]  dual_graph
 * \param [in]   dim Distributed nodal mesh handle
 *
 */
void
PDM_dmesh_nodal_dual_graph
(
  PDM_g_num_t   *vtx_dist,
  PDM_g_num_t   *elt_dist,
  int           *delt_vtx_idx,
  PDM_g_num_t   *delt_vtx,
  PDM_g_num_t  **delt_elt_idx,
  PDM_g_num_t  **delt_elt,
  PDM_MPI_Comm   comm
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // 0. transpose
  int* dvtx_elt_idx;
  PDM_g_num_t* dvtx_elt;

  PDM_dconnectivity_transpose(comm,
                              elt_dist, vtx_dist,
                              delt_vtx_idx,delt_vtx,
                              0, // not signed
                              &dvtx_elt_idx,&dvtx_elt);

  // 1. dual
  PDM_deduce_combine_connectivity_dual(comm,
                                       elt_dist, vtx_dist,
                                       delt_vtx_idx,delt_vtx,
                                       dvtx_elt_idx,dvtx_elt,
                                       0, // not signed
                                       delt_elt_idx,
                                       delt_elt);

}

/**
 * \brief  Return vertices distribution
 *
 * \param [in]  hdl  Distributed nodal mesh handle
 *
 * \return  A array of size \ref n_ranks + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_vtx_get
(
PDM_dmesh_nodal_t  *dmesh_nodal
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }

  PDM_DMesh_nodal_vtx_t *vtx = dmesh_nodal->vtx;

  return vtx->distrib;
}

/**
 * \brief  Return section distribution
 *
 * \param [in]  hdl        Distributed nodal mesh handle
 * \param [in]  id_section   Block identifier
 *
 * \return  A array of size \ref n_ranks + 1
 *
 */

const PDM_g_num_t *
PDM_DMesh_nodal_distrib_section_get
(
      PDM_dmesh_nodal_t   *dmesh_nodal,
      PDM_geometry_kind_t  geom_kind,
const int                  id_section
)
{
  if (dmesh_nodal == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad mesh nodal identifier\n");
  }
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmesh_nodal, geom_kind);
  return PDM_DMesh_nodal_elmts_distrib_section_get(dmne, id_section);
}

void
PDM_dmesh_nodal_transfer_to_new_dmesh_nodal
(
 PDM_dmesh_nodal_t   *dmn_in,
 PDM_dmesh_nodal_t   *dmn_out,
 PDM_geometry_kind_t  geom_kind,
 PDM_g_num_t         *dparent_vtx_distrib,
 PDM_g_num_t         *blk_parent_to_new_vtx_gnum
)
{
  PDM_dmesh_nodal_elmts_t* dmne_in = _get_from_geometry_kind(dmn_in, geom_kind);

  /* For all section we update connectivity */
  int n_section = dmne_in->n_section;
  for(int i = 0; i < n_section; ++i) {

    int id_section = dmne_in->sections_id[i];
    int                   n_elt     = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne_in, id_section);
    PDM_g_num_t          *delmt_vtx = PDM_DMesh_nodal_elmts_section_std_get  (dmne_in, id_section);
    PDM_Mesh_nodal_elt_t  t_elt     = PDM_DMesh_nodal_elmts_section_type_get (dmne_in, id_section);
    int n_vtx_per_elmt              = PDM_Mesh_nodal_n_vtx_elt_get           (t_elt    , 1);

    int n_connect = n_vtx_per_elmt * n_elt;

    PDM_block_to_part_t *btp_update_elmts_vtx = PDM_block_to_part_create (dparent_vtx_distrib,
                                                   (const PDM_g_num_t **) &delmt_vtx,
                                                                          &n_connect,
                                                                          1,
                                                                          dmn_in->comm);

    int stride_one = 1;
    PDM_g_num_t **tmp_delmt_vtx_new = NULL;
    PDM_block_to_part_exch (btp_update_elmts_vtx,
                             sizeof(PDM_g_num_t),
                             PDM_STRIDE_CST_INTERLACED,
                             &stride_one,
                    (void *) blk_parent_to_new_vtx_gnum,
                             NULL,
                  (void ***) &tmp_delmt_vtx_new);
    PDM_g_num_t *delmt_vtx_new = tmp_delmt_vtx_new[0];
    free(tmp_delmt_vtx_new);

    int id_section_post = PDM_DMesh_nodal_section_add(dmn_out,
                                                      geom_kind,
                                                      t_elt);

    PDM_DMesh_nodal_section_std_set(dmn_out,
                                    geom_kind,
                                    id_section_post,
                                    n_elt,
                                    delmt_vtx_new,
                                    PDM_OWNERSHIP_KEEP);

    PDM_block_to_part_free(btp_update_elmts_vtx);
  }
}



void
PDM_dmesh_nodal_transfer_to_new_dmesh_nodal_gen
(
 PDM_dmesh_nodal_t   *dmn_in,
 PDM_dmesh_nodal_t   *dmn_out,
 PDM_geometry_kind_t  geom_kind,
 PDM_g_num_t         *dparent_vtx_distrib,
 int                 *blk_parent_to_new_vtx_gnum_idx,
 PDM_g_num_t         *blk_parent_to_new_vtx_gnum
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn_in->comm, &i_rank);
  PDM_MPI_Comm_size(dmn_in->comm, &n_rank);
  PDM_dmesh_nodal_elmts_t* dmne_in = _get_from_geometry_kind(dmn_in, geom_kind);

  /* For all section we update connectivity */
  int n_section = dmne_in->n_section;
  for(int i_section = 0; i_section < n_section; ++i_section) {

    int id_section = dmne_in->sections_id[i_section];
    int                   n_elt     = PDM_DMesh_nodal_elmts_section_n_elt_get(dmne_in, id_section);
    PDM_g_num_t          *delmt_vtx = PDM_DMesh_nodal_elmts_section_std_get  (dmne_in, id_section);
    PDM_Mesh_nodal_elt_t  t_elt     = PDM_DMesh_nodal_elmts_section_type_get (dmne_in, id_section);
    int n_vtx_per_elmt              = PDM_Mesh_nodal_n_vtx_elt_get           (t_elt    , 1);

    int n_connect = n_vtx_per_elmt * n_elt;

    PDM_block_to_part_t *btp_update_elmts_vtx = PDM_block_to_part_create (dparent_vtx_distrib,
                                                   (const PDM_g_num_t **) &delmt_vtx,
                                                                          &n_connect,
                                                                          1,
                                                                          dmn_in->comm);

    int dn_vtx_old = dparent_vtx_distrib[i_rank+1] - dparent_vtx_distrib[i_rank];
    int         *dvtx_old_to_n = (int * ) malloc( dn_vtx_old * sizeof(int));
    for(int i = 0; i < dn_vtx_old; ++i) {
      dvtx_old_to_n[i] = blk_parent_to_new_vtx_gnum_idx[i+1] - blk_parent_to_new_vtx_gnum_idx[i];
    }

    // int stride_one = 1;
    int         **tmp_delmt_vtx_new_n = NULL;
    PDM_g_num_t **tmp_delmt_vtx_new = NULL;
    PDM_block_to_part_exch (btp_update_elmts_vtx,
                             sizeof(PDM_g_num_t),
                             PDM_STRIDE_VAR_INTERLACED,
                             dvtx_old_to_n,
                    (void *) blk_parent_to_new_vtx_gnum,
                             &tmp_delmt_vtx_new_n,
                  (void ***) &tmp_delmt_vtx_new);
    PDM_g_num_t *delmt_vtx_new = tmp_delmt_vtx_new[0];
    free(tmp_delmt_vtx_new);
    free(tmp_delmt_vtx_new_n[0]);
    free(tmp_delmt_vtx_new_n);
    free(dvtx_old_to_n);

    int id_section_post = PDM_DMesh_nodal_section_add(dmn_out,
                                                      geom_kind,
                                                      t_elt);

    PDM_DMesh_nodal_section_std_set(dmn_out,
                                    geom_kind,
                                    id_section_post,
                                    n_elt,
                                    delmt_vtx_new,
                                    PDM_OWNERSHIP_KEEP);

    PDM_block_to_part_free(btp_update_elmts_vtx);
  }
}




void
PDM_dmesh_nodal_dump_vtk
(
       PDM_dmesh_nodal_t   *dmn,
       PDM_geometry_kind_t  geom_kind,
 const char                *filename_patter
)
{
  /* TODO: add groups as scalar, elt-based fields */

  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(dmn->comm, &i_rank);
  PDM_MPI_Comm_size(dmn->comm, &n_rank);

  int* sections_id = PDM_DMesh_nodal_sections_id_get(dmn, geom_kind);
  int n_section    = PDM_DMesh_nodal_n_section_get(dmn, geom_kind);

  const char *field_name = "group";
  int n_field = 0;
  _pdm_dmesh_nodal_elts_t *dmne = NULL;
  int *delt_group_idx = NULL;
  int *delt_group     = NULL;
  double **field = NULL;
  if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
    dmne = dmn->ridge;
  } else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC && dmn->mesh_dimension == 3) {
    dmne = dmn->surfacic;
  } else if (geom_kind == PDM_GEOMETRY_KIND_CORNER) {
    dmne = dmn->corner;
  }

  PDM_g_num_t *distrib_elt = NULL;
  if (dmne != NULL) {
    distrib_elt = PDM_compute_uniform_entity_distribution(dmn->comm,
                                                          dmne->n_g_elmts);

    PDM_dgroup_entity_transpose(dmne->n_group_elmt,
                                dmne->dgroup_elmt_idx,
                                dmne->dgroup_elmt,
                (PDM_g_num_t *) distrib_elt,
                                &delt_group_idx,
                                &delt_group,
                                dmn->comm);
  }

  PDM_g_num_t shift = 0;
  for(int i_section = 0; i_section < n_section; ++i_section) {
    int id_section = sections_id[i_section];
    int order;
    const char *ho_ordering  = NULL;
    PDM_g_num_t *dconnec     = NULL;
    int         *dconnec_idx = NULL;
    const PDM_g_num_t    *delmt_distribution = PDM_DMesh_nodal_distrib_section_get(dmn, geom_kind, id_section);
    int                   n_elt              = PDM_DMesh_nodal_section_n_elt_get  (dmn, geom_kind, id_section);
    PDM_Mesh_nodal_elt_t  t_elt              = PDM_DMesh_nodal_section_type_get   (dmn, geom_kind, id_section);
    if (t_elt == PDM_MESH_NODAL_POLY_2D) {
      PDM_DMesh_nodal_section_poly2d_get(dmn,
                                         geom_kind,
                                         id_section,
                                         &dconnec_idx,
                                         &dconnec);
    }
    else {
      assert(t_elt != PDM_MESH_NODAL_POLY_3D);
      dconnec     = PDM_DMesh_nodal_section_std_ho_get (dmn, geom_kind, id_section, &order, &ho_ordering);
      int strid   = PDM_Mesh_nodal_n_vtx_elt_get(t_elt, order);
      dconnec_idx = PDM_array_new_idx_from_const_stride_int(strid, n_elt);
    }

    PDM_g_num_t *delmt_ln_to_gn = (PDM_g_num_t * ) malloc( (n_elt  ) * sizeof(PDM_g_num_t));
    for(int i = 0; i < n_elt; ++i) {
      delmt_ln_to_gn[i] = delmt_distribution[i_rank] + i + 1;
    }

    PDM_g_num_t *pvtx_ln_to_gn;
    int         *pcell_vtx_idx;
    int         *pcell_vtx;
    int          pn_vtx;
    PDM_part_dconnectivity_to_pconnectivity_sort_single_part(dmn->comm,
                                                             delmt_distribution,
                                                             dconnec_idx,
                                                             dconnec,
                                                             n_elt,
                                    (const PDM_g_num_t *)    delmt_ln_to_gn,
                                                            &pn_vtx,
                                                            &pvtx_ln_to_gn,
                                                            &pcell_vtx_idx,
                                                            &pcell_vtx);

    /*
     * Coordinates
     */
    PDM_g_num_t *vtx_distrib = PDM_dmesh_nodal_vtx_distrib_get(dmn);
    double      *dvtx_coord  = PDM_DMesh_nodal_vtx_get(dmn);
    // int          dn_vtx   = PDM_DMesh_nodal_n_vtx_get(dln->dmesh_nodal_in);
    // assert(dn_vtx == (vtx_distrib[i_rank+1]-vtx_distrib[i_rank]));
    double** tmp_pvtx_coord = NULL;
    PDM_part_dcoordinates_to_pcoordinates(dmn->comm,
                                          1,
                                          vtx_distrib,
                                          dvtx_coord,
                                          &pn_vtx,
                   (const PDM_g_num_t **) &pvtx_ln_to_gn,
                                          &tmp_pvtx_coord);

    double* pvtx_coord_out = tmp_pvtx_coord[0];


    /*
     * Groups
     */
    if (dmne != NULL) {
      for(int i = 0; i < n_elt; ++i) {
        delmt_ln_to_gn[i] += shift;
      }

      int **tmp_elt_group_idx = NULL;
      int **tmp_elt_group     = NULL;
      PDM_part_dentity_group_to_pentity_group(dmn->comm,
                                              1,
                                              distrib_elt,
                                              delt_group_idx,
                                              delt_group,
                                              &n_elt,
                      (const PDM_g_num_t **)  &delmt_ln_to_gn,
                                              &tmp_elt_group_idx,
                                              &tmp_elt_group);
      int *pelt_group_idx = tmp_elt_group_idx[0];
      int *pelt_group     = tmp_elt_group    [0];
      free (tmp_elt_group_idx);
      free (tmp_elt_group);

      n_field = 1;
      field = malloc (sizeof(double *) * n_field);
      field[0] = malloc (sizeof(double) * n_elt);
      for (int i = 0; i < n_elt; i++) {
        assert (pelt_group_idx[i+1] == pelt_group_idx[i] + 1);
        field[0][i] = (double) pelt_group[i];
      }
      free (pelt_group);
      free (pelt_group_idx);
    }

    /*
     *  Dump
     */
    char filename[999];
    sprintf(filename, "%s_section_%2.2d_%2.2d.vtk", filename_patter, i_section, i_rank);
    if (t_elt == PDM_MESH_NODAL_POLY_2D) {
      PDM_vtk_write_polydata(filename,
                             pn_vtx,
                             pvtx_coord_out,
                             pvtx_ln_to_gn,
                             n_elt,
                             pcell_vtx_idx,
                             pcell_vtx,
                             delmt_ln_to_gn,
                             NULL);
    }
    else {
      PDM_vtk_write_std_elements_ho(filename,
                                    order,
                                    pn_vtx,
                                    pvtx_coord_out,
                                    pvtx_ln_to_gn,
                                    t_elt,
                                    n_elt,
                                    pcell_vtx,
                                    delmt_ln_to_gn,
                                    n_field,
                  (const char   **) &field_name,
                  (const double **) field);
    }
    free(tmp_pvtx_coord);
    free(pvtx_ln_to_gn);
    free(pcell_vtx_idx);
    free(pcell_vtx);

    if (t_elt != PDM_MESH_NODAL_POLY_2D) {
      free(dconnec_idx);
    }
    free(delmt_ln_to_gn);

    free(pvtx_coord_out);

    shift += delmt_distribution[n_rank];

    if (dmne != NULL) {
      free (field[0]);
      free (field);
    }
  }

  if (dmne != NULL) {
    free (delt_group_idx);
    free (delt_group);
    free (distrib_elt);
  }
}



void
PDM_dmesh_nodal_reorder
(
 PDM_dmesh_nodal_t *dmesh_nodal,
 const char        *ordering_name
 )
{
  PDM_geometry_kind_t start = PDM_GEOMETRY_KIND_MAX;
  if (dmesh_nodal->mesh_dimension == 2) {
    start = PDM_GEOMETRY_KIND_SURFACIC;
  } else if (dmesh_nodal->mesh_dimension == 3) {
    start = PDM_GEOMETRY_KIND_VOLUMIC;
  }

  PDM_dmesh_nodal_elmts_t *dmne = NULL;

  for (PDM_geometry_kind_t geom_kind = start; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {

    if (geom_kind == PDM_GEOMETRY_KIND_RIDGE) {
      dmne = dmesh_nodal->ridge;
    } else if (geom_kind == PDM_GEOMETRY_KIND_SURFACIC) {
      dmne = dmesh_nodal->surfacic;
    } else if (geom_kind == PDM_GEOMETRY_KIND_VOLUMIC) {
      dmne = dmesh_nodal->volumic;
    }

    int *sections_id = PDM_DMesh_nodal_sections_id_get(dmesh_nodal, geom_kind);
    int n_section    = PDM_DMesh_nodal_n_section_get  (dmesh_nodal, geom_kind);

    for (int isection = 0; isection < n_section; isection++) {

      PDM_DMesh_nodal_elmts_section_std_ho_reorder(dmne,
                                                   sections_id[isection],
                                                   ordering_name);
    }
  }

}



PDM_part_mesh_nodal_elmts_t*
PDM_dmesh_nodal_to_part_mesh_nodal_elmts
(
 PDM_dmesh_nodal_t            *dmn,
 PDM_geometry_kind_t           geom_kind,
 int                           n_part,
 int                          *pn_vtx,
 PDM_g_num_t                 **vtx_ln_to_gn,
 int                          *pn_elmt,
 PDM_g_num_t                 **elmt_ln_to_gn,
 PDM_g_num_t                 **pparent_entitity_ln_to_gn
)
{
  PDM_dmesh_nodal_elmts_t* dmne = _get_from_geometry_kind(dmn, geom_kind);
  return PDM_dmesh_nodal_elmts_to_part_mesh_nodal_elmts(dmne, n_part, pn_vtx, vtx_ln_to_gn, pn_elmt, NULL, elmt_ln_to_gn, pparent_entitity_ln_to_gn);
}

const double *
PDM_dmesh_nodal_global_extents_get
(
 PDM_dmesh_nodal_t         *dmn
)
{
  if (dmn->is_computed_g_extents == PDM_FALSE) {
    //TO DO Bastien: if high-order, compute extents using Bézier control points

    double l_min[3] = { HUGE_VAL};
    double l_max[3] = {-HUGE_VAL};

    for (int i = 0; i < dmn->vtx->n_vtx; i++) {
      for (int j = 0; j < 3; j++) {
        double x = dmn->vtx->_coords[3*i + j];
        l_min[j] = PDM_MIN(l_min[j], x);
        l_max[j] = PDM_MAX(l_max[j], x);
      }
    }

    PDM_MPI_Allreduce(l_min, dmn->g_extents,   3, PDM_MPI_DOUBLE, PDM_MPI_MIN, dmn->comm);
    PDM_MPI_Allreduce(l_max, dmn->g_extents+3, 3, PDM_MPI_DOUBLE, PDM_MPI_MAX, dmn->comm);

    dmn->is_computed_g_extents = PDM_TRUE;
  }

  return dmn->g_extents;
}


int
PDM_dmesh_nodal_have_ho
(
 PDM_dmesh_nodal_t         *dmesh_nodal
)
{
  int have_ho_volumic  = 0;
  int have_ho_surfacic = 0;
  int have_ho_ridge    = 0;
  int have_ho_corner   = 0;
  if(dmesh_nodal->volumic  != NULL) {
    have_ho_volumic = PDM_dmesh_nodal_elmts_have_ho(dmesh_nodal->volumic);
  }
  if(dmesh_nodal->surfacic != NULL) {
    have_ho_surfacic = PDM_dmesh_nodal_elmts_have_ho(dmesh_nodal->surfacic);
  }
  if(dmesh_nodal->ridge    != NULL) {
    have_ho_ridge = PDM_dmesh_nodal_elmts_have_ho(dmesh_nodal->ridge);
  }
  if(dmesh_nodal->corner   != NULL) {
    have_ho_corner = PDM_dmesh_nodal_elmts_have_ho(dmesh_nodal->corner);
  }

  return have_ho_volumic || have_ho_surfacic || have_ho_ridge || have_ho_corner;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
