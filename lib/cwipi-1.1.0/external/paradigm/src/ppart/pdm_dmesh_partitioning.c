
/*============================================================================
 * Parallel partitioning
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  System headers
 *----------------------------------------------------------------------------*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_config.h"
#include "pdm_priv.h"
#include "pdm_timer.h"
#include "pdm_mpi.h"
#include "pdm_mpi_ext_dependencies.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_sort.h"
#include "pdm_binary_search.h"

#include "pdm_dmesh_partitioning.h"
#include "pdm_para_graph_dual.h"

/*----------------------------------------------------------------------------
 *  Optional headers
 *----------------------------------------------------------------------------*/

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

/**
 * \struct _dmesh_partitioning_t
 * \brief  Define a global numberring
 *
 */
typedef struct {


  int set_flags;
  int input_flags;
  int query_flags;

  PDM_g_num_t** p_face_ln_to_gn;

  // int* ownership = NULL; // 64 et on met le code dedans


} _dmesh_partitioning_t;

/*============================================================================
 * Global variable
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*=============================================================================
 * Public function definitions
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
)
{
  PDM_UNUSED(comm);
  PDM_UNUSED(split_method);

  _dmesh_partitioning_t *_dmp = (_dmesh_partitioning_t *) malloc(sizeof(_dmesh_partitioning_t));

  return (PDM_dmesh_partitioning_t *) _dmp;
}


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
)
{
  PDM_UNUSED(queries_flags);
  PDM_UNUSED(input_flags);
  PDM_UNUSED(dmp);

  // printf("PDM_dmesh_partitioning_compute::dmesh_partitioning_id :: %d \n", dmesh_partitioning_id);
  printf("PDM_dmesh_partitioning_compute::input_flags   :: %d \n", input_flags);
  printf("PDM_dmesh_partitioning_compute::queries_flags :: %d \n", queries_flags);

  // _dmesh_partitioning_t* _dmp = _get_from_id(dmesh_partitioning_id);

  /* Parse input flags - to know wich data is available to build a partiotionning */

  if( PDM_HASFLAG(input_flags, PDM_PART_FACE_CELL) ) {

    /* Sanity check */
    // assert(_dmp->dface_cell != NULL);

    // Rebuild CELL_FACE for graph computation
    // dual
    // Split

  } else if ( PDM_HASFLAG(input_flags, PDM_PART_CELL_FACE) ) {

    /* Sanity check */
    // assert(_dmp->dcell_face != NULL);

    // dual
    // Split


  } else if ( PDM_HASFLAG(input_flags, PDM_PART_CELL_VTX) ) {

    /* Sanity check */
    // assert(_dmp->dcell_vtx != NULL);
    // DG Case
    // dmesh_nodal_id
    // Il faut faire un appel a dmesh_nodal pour calculer le
    // face_cell (ou le cell_face) pour pouvoir splitter

  } else {

    PDM_error(__FILE__, __LINE__, 0,
              "PDM_dmesh_partitioning_compute error : None of input data can split graph. ");
  }

  /* Compute and split */


  /* Post-treatment */
  // D'abord ln_to_gn > connectivité > renum > graph_comm

  // renum_entiies ()
  //

  // Multigrille ou uniform_refinement


  /*
   * Le probleme est qu'il faut le faire dans un ordre précis si il y a une dependances -_-
   *   - Par exemple : On doit reconstruire cell_face > face_vtx
   *     Il y a une dependances entre les demandes comment trié ?
   */
  // if()()
  // do_..
  // do_..
  // do_..

  // while( calc < demande ){
  // if compute alreay done
  // demande +=1
  // sinon
  // Je cherche dans le truc qui correspond a ma demande
  // }



}


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
)
{
  PDM_UNUSED(dmp);
  PDM_UNUSED(dmesh_id);
  // _dmesh_partitioning_t* _dmp = _get_from_id(dmesh_partitioning_id);
  // ---> Depuis l'interface du dmesh
  // Recopie
  // PDM_dmesh_dims_get(dmesh_id, &)
}

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
 *
 * \param [in]   dmpartitioning_id        ppart identifier
 *
 */
// void
// PDM_dmesh_partitioning_set
// (
//  const int     dmesh_partitioning_id,
//  const int     input_field_key,
//        void ***field
// )
// {
//   // int owner      = PDM_HASFLAG  (dmm->set_flags, input_field_key);
//   // assert(owner != 0);

// }

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
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
)
{
  PDM_UNUSED(dmp);
  PDM_UNUSED(field);
  int field_key  = input_field_key;
  //int owner      = PDM_HASFLAG  (field_key, PDM_PART_OWNDATA);
  int field_name = PDM_UNSETFLAG(field_key, PDM_PART_OWNDATA); // Il reste que le nom

  // Search the correct field who want
  switch (field_name) {
    case PDM_PART_FACE_CELL:
    {
      printf("PDM_dmesh_partitioning_get::PDM_PART_FACE_CELL\n");

      break;
    }
    case PDM_PART_CELL_FACE:
    {
      printf("PDM_dmesh_partitioning_get::PDM_PART_CELL_FACE\n");

      break;
    }
    default:
    {
      printf("PDM_dmesh_partitioning_get::NOT FOUND\n");
      break;
    }
  }
}

/**
 *
 * \brief Return _dmesh_partitioning_t object from it identifier
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
)
{
  PDM_UNUSED(dmp);
  PDM_UNUSED(part_id);
  PDM_UNUSED(field);
  int field_key  = input_field_key;
  //  int owner      = PDM_HASFLAG  (field_key, PDM_PART_OWNDATA);
  int field_name = PDM_UNSETFLAG(field_key, PDM_PART_OWNDATA); // Il reste que le nom

  // Search the correct field who want
  switch (field_name) {
    case PDM_PART_FACE_CELL:
    {
      printf("PDM_dmesh_partitioning_part_get::PDM_PART_FACE_CELL\n");

      break;
    }
    case PDM_PART_CELL_FACE:
    {
      printf("PDM_dmesh_partitioning_part_get::PDM_PART_CELL_FACE\n");

      break;
    }
    default:
    {
      printf("PDM_dmesh_partitioning_part_get::NOT FOUND\n");
      break;
    }
  }
}

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
)
{
  _dmesh_partitioning_t* _dmp = (_dmesh_partitioning_t *) dmp;

  free (_dmp);
  _dmp = NULL;
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
