#ifndef PDM_PART_COARSE_MESH_PRIV_H
#define	PDM_PART_COARSE_MESH_PRIV_H

/*
 * File:   pdm_part_coarse_mesh_priv.h
 * Author: jmagnene
 *
 * Created on July 12, 2016, 10:54 AM
 */

#include "pdm_part_priv.h"
#include "pdm_timer.h"
#include "pdm_fortran_to_c_string.h"
#include "pdm_mpi.h"

#ifdef	__cplusplus
extern "C" {
#endif

/*============================================================================
 * Type definitions
 *============================================================================*/

/**
 * \struct _coarse_part_t
 * \brief  Coarse partition object
 *
 * _part_t defines a coarse partition
 *
 */

typedef struct  {

  _part_t      *part;        //Coarse mesh

  int           n_coarse_cell_wanted;    /*!< Number Cell wanted for agglomeration     */

  // int *cell_weight;    /*!< Integer weight for graoh partitionning  */
  // int *face_weight;    /*!< Number Cell wanted for agglomeration     */

  int *coarse_cell_cell_idx;    //Array of indexes of the connected partitions (size : n_coarse_cell + 1)

  int *coarse_cell_cell;       //Partitioning array (size : coarse_cell_cell_idx[n_coarse_cell])

  int *coarse_face_group_to_fine_face_group; //Coarse face group - fine face group connectivity (size = face_group_idx[n_face_group])

  int *coarse_face_to_fine_face; //Coarse face - fine face connectivity (size = nCoarseFace)

  int *coarse_vtx_to_fine_vtx;   //Coarse vertex - fine vertex connectivity (size = nCoarseVtx)

  void *specific_data;       /*!< Specific data      */

} _coarse_part_t;

/**
 * \struct PDM_part_renum_fct_t
 *
 * \brief  Function pointer used to define a renumbering specific to agglomeration array
 *
 */

typedef void (*PDM_coarse_mesh_renum_fct_t) (_coarse_part_t  *part);



/**
 * \struct _coarse_part_t
 * \brief  Coarse partition object
 *
 * _coarse_mesh_t defines a coarse mesh
 *
 */

struct _coarse_mesh_t {

  /* Partitions */

  int n_part;        /*!< Number of partitions to define
                                      on this process */

  int method;       /*!< Partitioning method */
  int n_total_part;       /*!< Total number of partitions */
  int n_face_group;   /*!< Number of boundaries */

  int have_cell_tag;
  int have_face_tag;
  int have_vtx_tag;
  int have_cell_weight;
  int have_face_weight;
  int have_face_group;

  void *specific_data;

  PDM_coarse_mesh_renum_fct_t specific_func;       /*!< Specific function  */

  /* Reordering */
  int        renum_face_method;               /*!< Renumbering face method       */
  int        renum_cell_method;               /*!< Renumbering cell method       */
  int        n_property_cell;                   /*!< Size of cells properties      */
  int        n_property_face;                   /*!< Size of faces properties      */
  const int* renum_properties_cell;           /*!< Renumbering cells properties  */
  const int* renum_properties_face;           /*!< Renumbering faces properties  */

  //TIMER

  PDM_timer_t *timer;             /*!< Timer */

  double times_elapsed [18];          /*!< Elapsed times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu[18];             /*!< CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_u[18];           /*!< User CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */

  double times_cpu_s[18];          /*!< Systeme CPU times :
                                      - Total,
                                      - build dualgraph,
                                      - split graph
                                      - build meshes partition */


  /* Communicator */

  PDM_MPI_Comm  comm;   /*!< Communicator */

  _part_t        **part_ini;               /*!< Partition: fine mesh                            */

  _coarse_part_t **part_res;               //Coarse mesh



} ;


/**
 * \struct PDM_part_renum_fct_t
 *
 * \brief  Function pointer used to define a coarse mesh method
 *
 */
typedef void (*PDM_coarse_mesh_fct_t) (struct _coarse_mesh_t  *cm,
                                       const int       i_part,
                                       int             *n_coarse_cell_computed,
                                       int             *cell_cell_idx,
                                       int             *cell_cell,
                                       int             *cell_part);

/**
 * \struct _coarse_mesh_method_t
 * \brief coarse mesh method
 *
 */

typedef struct _renum_method_t {

  char                  *name;  /*!< Name of method          */
  PDM_coarse_mesh_fct_t   fct;  /*!< Agglometration function */

} _coarse_mesh_method_t;


/*============================================================================
 * Private function definitions
 *============================================================================*/


/**
 *
 * \brief Return an initialized coarse part object
 *
 */

static inline _coarse_part_t *
_coarse_part_create
(
void
 )
{
  _coarse_part_t *cp = (_coarse_part_t *) malloc(sizeof(_coarse_part_t));
  cp->part = _part_create();

  cp->coarse_cell_cell = NULL;

  cp->coarse_cell_cell_idx = NULL;

  cp->coarse_face_group_to_fine_face_group = NULL;

  cp->coarse_face_to_fine_face = NULL;

  cp->coarse_vtx_to_fine_vtx = NULL;

  cp->specific_data = NULL;

  return cp;

}


/*============================================================================
 * Public function prototypes
 *============================================================================*/



/**
 *
 * \brief Add a new coarse mesh method
 *
 * \param [in]      name          Mesh entity to renumber
 * \param [in]      fct           Function
 *
 */

int
PDM_coarse_mesh_method_add
(
 const char                 *name,     /*!< Name          */
 PDM_coarse_mesh_fct_t       fct       /*!< Function      */
 );


/**
 *
 * \brief Get index of a coarse mesh method from it's name
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

int
PDM_coarse_mesh_method_idx_get
(
const char *name
 );


/**
 *
 * \brief Get name of a coarse mesh method from it's index
 *
 * \param [in]  name   Name of the method
 *
 * \return Index (-1 if not found)
 */

void
PDM_coarse_mesh_method_name_get_cf
(
 const int  idx,
 char      *name,
 int       *l_name
 );

char *
PDM_coarse_mesh_method_name_get
(
const int id
 );


/**
 *
 * \brief Get the number of coarse mesh method
 *
 * \return Number of methods
 *
 */

int
PDM_coarse_mesh_method_n_get
(
void
);


/**
 *
 * \brief Purge coarse mesh methods catalog
 *
 */

void
PDM_coarse_mesh_method_purge
(
void
 );


/**
 *
 * \brief Load local coarse mesh methods
 *
 */

void
PDM_coarse_mesh_method_load_local
(
void
 );


/**
 *
 * \brief Return coarse mesh object from its identifier
 *
 * \param [in]   cmId        Coarse mesh identifier
 *
 */

// _coarse_mesh_t *
// PDM_part_coarse_mesh_get_from_id
// (
//  int  cmId
//  );


#ifdef	__cplusplus
}
#endif

#endif	/* PDM_PART_COARSE_MESH_PRIV_H */
