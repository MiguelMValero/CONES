/*
 * \file
 */

#ifndef __PDM_DMESH_PARTITIONING_H__
#define __PDM_DMESH_PARTITIONING_H__

#include <stdio.h>
#include "pdm.h"
#include "pdm_mpi.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

typedef struct _dmesh_partitioning_t PDM_dmesh_partitioning_t;

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/

typedef enum {
  PDM_PART_NONE               = 0x0000000, /*=* empty flag, all values set to false */
  PDM_PART_NULL               = 0xFFFFFFF, /*=* unsignificant flags value           */
  PDM_PART_CELL_VTX           = 0x0000001, /*=*  0                                  */
  PDM_PART_FACE_CELL          = 0x0000002, /*=*  1                                  */
  PDM_PART_CELL_FACE          = 0x0000004, /*=*  2                                  */
  PDM_PART_EDGE_VTX           = 0x0000008, /*=*  3                                  */
  PDM_PART_FACE_VTX           = 0x0000010, /*=*  4                                  */
  PDM_PART_CELL_WEIGHT        = 0x0000020, /*=*  5                                  */
  PDM_PART_EDGE_WEIGHT        = 0x0000040, /*=*  6                                  */
  PDM_PART_CELL_TAG           = 0x0000080, /*=*  7                                  */
  PDM_PART_FACE_TAG           = 0x0000100, /*=*  8                                  */
  PDM_PART_VTX_TAG            = 0x0000200, /*=*  9                                  */
  PDM_PART_FACE_GROUP         = 0x0000400, /*=*  10                                 */
  PDM_PART_CELL_GROUP         = 0x0000800, /*=*  11                                 */
  PDM_PART_VTX_GROUP          = 0x0001000, /*=*  12                                 */
  PDM_PART_EDGE_GROUP         = 0x0002000, /*=*  13                                 */
  PDM_PART_CELL_PART          = 0x0004000, /*=*  14                                 */
  PDM_PART_VTX_PART           = 0x0008000, /*=*  15                                 */
  PDM_PART_GRAPH_COMM_FACE    = 0x0010000, /*=*  16                                 */
  PDM_PART_GRAPH_COMM_VTX     = 0x0020000, /*=*  17                                 */
  PDM_PART_GRAPH_COMM_EDGE    = 0x0040000, /*=*  18                                 */
  PDM_PART_VTX_COORD          = 0x0080000, /*=*  19                                 */
  PDM_PART_CELL_LN_TO_GN      = 0x0100000, /*=*  20                                 */
  PDM_PART_FACE_LN_TO_GN      = 0x0200000, /*=*  21                                 */
  PDM_PART_VTX_LN_TO_GN       = 0x0400000, /*=*  22                                 */
  PDM_PART_FACEGROUP_LN_TO_GN = 0x0800000, /*=*  23                                 */
  PDM_PART_CELLGROUP_LN_TO_GN = 0x1000000, /*=*  24                                 */
  PDM_PART_VTXGROUP_LN_TO_GN  = 0x2000000, /*=*  25                                 */
  PDM_PART_UNUSED6            = 0x4000000, /*=*  26                                 */
  PDM_PART_OWNDATA            = 0x8000000, /*=*  27                                 */
} PDM_partitioning_option_t;

/**
 * \enum PDM_part_split_t
 * \brief Split method
 *
 */

typedef enum {
  PDM_PARTITIONING_WITH_PARMETIS = 1,
  PDM_PARTITIONING_WITH_PTSCOTCH = 2,
  PDM_PARTITIONING_WITH_HILBERT  = 3
} PDM_partitioning_method_t;

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/**
 *
 * \brief Create a new \ref PDM_dmesh_partitioning_t object
 *
 * \param [in]   comm          MPI Communicator
 * \param [in]   split_method  Graph splitting method
 *
 * \return       Pointer to a new \ref PDM_dmesh_partitioning_t object
 */
PDM_dmesh_partitioning_t *
PDM_dmesh_partitioning_create
(
 const PDM_MPI_Comm              comm,
 const PDM_partitioning_method_t split_method
);

/**
 *
 * \brief Compute partitioning
 *
 * \param [in] dmp            Pointer to \ref PDM_dmesh_partitioning_t object
 * \param [in] input_flags    ?
 * \param [in] queries_flags  ?
 *
 */
void
PDM_dmesh_partitioning_compute
(
 PDM_dmesh_partitioning_t *dmp,
 const int                 input_flags,
 const int                 queries_flags
);

/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in] dmp            Pointer to \ref PDM_dmesh_partitioning_t object
 * \param [in] dmesh_id       Dmesh identifier
 *
 */
void
PDM_dmesh_partitioning_set_from_dmesh
(
 PDM_dmesh_partitioning_t *dmp,
 const int                 dmesh_id
);

/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
// void
// PDM_dmesh_partitioning_set_from_dmesh_nodal_id
// (
//  const int dmesh_partitioning_id,
//  const int dmesh_nodal_id
// );

/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
// void
// PDM_dmesh_partitioning_get_part_mesh_nodal_id
// (
//  const int dmesh_partitioning_id,
// );
// return mesh_nodal_id


/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
void
PDM_dmesh_partitioning_part_get
(
 PDM_dmesh_partitioning_t  *dmp,
 const int                  part_id,
 const int                  input_field_key,
      void                **field
);


/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
// void
// PDM_dmesh_partitioning_part_get
// (
//  const int   dmesh_partitioning_id,
//  const int   part_id,
//  const int   input_field_key,
//       void **field
//       void **field_idx
// );

/**
 *
 * \brief Setup structure with a dmesh
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
void
PDM_dmesh_partitioning_get
(
 PDM_dmesh_partitioning_t   *dmp,
 const int                   input_field_key,
       void               ***field
);
/**
 *
 * \brief Free
 *
 * \param [in]   id           Identifier
 *
 */
void
PDM_dmesh_partitioning_free
(
 PDM_dmesh_partitioning_t *dmp
);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif  /* __PDM_DMESH_PARTITIONING_H__ */
