/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_binary_search.h"
#include "pdm_io.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"
#include "pdm_vtk.h"
#include "pdm_logging.h"
#include "pdm_gnum.h"
#include "pdm_array.h"
#include "pdm_distrib.h"
#include "pdm_multipart.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_block_to_block.h"

#include "pdm_reader_stl.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type
 *============================================================================*/

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*=============================================================================
 * Private function prototypes
 *============================================================================*/

static void
_read_distributed_stl
(
 PDM_MPI_Comm   comm,
 const char    *filename,
 int          **dface_vtx_idx,
 PDM_g_num_t  **dface_vtx,
 double       **dvtx_coord,
 PDM_g_num_t  **distrib_vtx,
 PDM_g_num_t  **distrib_face
)
{
  int debug = 0;

  PDM_UNUSED(comm);
  PDM_UNUSED(filename);
  PDM_UNUSED(dface_vtx_idx);
  PDM_UNUSED(dface_vtx);
  PDM_UNUSED(dvtx_coord);
  PDM_UNUSED(distrib_vtx);
  PDM_UNUSED(distrib_face);

  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // assert(n_rank == 1);
  int n_vtx  = 0;
  int n_face = 0;
  double      *face_normal    = NULL;
  double      *vtx_coord      = NULL;
  int         *face_vtx_idx   = NULL;
  int         *face_vtx_n     = NULL;
  // PDM_g_num_t *face_vtx       = NULL;
  int *face_vtx       = NULL;
  double      *face_vtx_coord = NULL;

  if (i_rank == 0) {
    FILE *f = fopen(filename, "r");

    if (f == NULL) {
      PDM_error(__FILE__, __LINE__, 0, "Could not read file %s\n", filename);
    }

    char header[999];
    fscanf(f, "%[^\n]", header);

    if (debug) {
      printf("header = %s \n", header);
    }

    int end = 0;
    while(end == 0) {

      int stat = fscanf(f, "%s", header);

      // printf("header = %s \n", header);

      if(strcmp(header, "vertex") == 0) {
        n_vtx += 1;
      }
      if(strcmp(header, "facet") == 0) {
        n_face += 1;
      }

      // if(strcmp(header, "OpenSCAD_Model") == 0) {
      //   end = 1;
      // }
      if(stat == EOF) {
        end = 1;
      }
    }

    if (debug) {
      printf("n_face = %i \n", n_face);
      printf("n_vtx  = %i \n", n_vtx );
    }


    face_normal    = (double * ) malloc(3 * n_face    * sizeof(double));
    vtx_coord      = (double * ) malloc(3 * n_vtx     * sizeof(double));
    face_vtx_idx   = (int    * ) malloc( (n_face + 1) * sizeof(int   ));
    face_vtx_n     = (int    * ) malloc( (n_face + 1) * sizeof(int   ));
    face_vtx       = (int    * ) malloc(3 * n_face    * sizeof(int   ));
    face_vtx_coord = (double * ) malloc(9 * n_face    * sizeof(double));


    fseek(f, 0, SEEK_SET);
    fscanf(f, "%[^\n]", header);

    int i_vtx  = 0;
    int i_face = 0;
    end = 0;
    face_vtx_idx[0] = 0;
    while(end == 0) {

      int stat = fscanf(f, "%s", header);

      // printf("header = %s \n", header);

      if(strcmp(header, "vertex") == 0) {
        fscanf(f, "%lf %lf %lf \n", &vtx_coord[3*i_vtx  ],
                                          &vtx_coord[3*i_vtx+1],
                                          &vtx_coord[3*i_vtx+2]);
        // printf("%f %f %f \n", vtx_coord[3*i_vtx  ],
        //                                   vtx_coord[3*i_vtx+1],
        //                                   vtx_coord[3*i_vtx+2]);

        face_vtx[face_vtx_idx[i_face]] = i_vtx+1;//(PDM_g_num_t) i_vtx+1;

        face_vtx_coord[3*face_vtx_idx[i_face]  ] = vtx_coord[3*i_vtx  ];
        face_vtx_coord[3*face_vtx_idx[i_face]+1] = vtx_coord[3*i_vtx+1];
        face_vtx_coord[3*face_vtx_idx[i_face]+2] = vtx_coord[3*i_vtx+2];

        face_vtx_idx[i_face]++;
        i_vtx += 1;

        // if(i_vtx == 12) {
        //   return;
        // }
      }
      if(strcmp(header, "facet") == 0) {
        face_vtx_idx[i_face+1] = face_vtx_idx[i_face];
        i_face += 1;
      }

      // if(strcmp(header, "OpenSCAD_Model") == 0) {
      // if(strcmp(header, "endsolid") == 0) {
      //   end = 1;
      // }
      if(stat == EOF) {
        end = 1;
      }

    }


    if (debug) {
      const char outfilename[999] = "debug_stl.vtk";
      PDM_vtk_write_point_cloud(outfilename,
                                n_vtx,
                                vtx_coord,
                                NULL,
                                NULL);

      const char outfilenamep[999] = "debug_stl_face_vtx.vtk";
      PDM_vtk_write_polydata(outfilenamep,
                             n_vtx,
                             vtx_coord,
                             NULL,
                             n_face,
                             face_vtx_idx,
                             face_vtx,
                             NULL,
                             NULL);
    }

    for(int i = 0; i < n_face; ++i) {
      face_vtx_n[i] = face_vtx_idx[i+1] - face_vtx_idx[i];
    }

    fclose(f);
  }

  // Re-repart surface among all procs
  PDM_g_num_t* init_distrib_vtx  = PDM_compute_entity_distribution(comm, n_vtx );
  PDM_g_num_t* init_distrib_face = PDM_compute_entity_distribution(comm, n_face);

  PDM_g_num_t _n_g_face = n_face;
  PDM_MPI_Bcast(&_n_g_face, 1, PDM__PDM_MPI_G_NUM, 0, comm);

  // PDM_g_num_t* _distrib_vtx  = PDM_compute_uniform_entity_distribution(comm, n_vtx );
  PDM_g_num_t* _distrib_face = PDM_compute_uniform_entity_distribution(comm, _n_g_face);

  PDM_block_to_block_t *btb_face = PDM_block_to_block_create(init_distrib_face,
                                                             _distrib_face,
                                                             comm);

  int dn_face_end = _distrib_face[i_rank+1] - _distrib_face[i_rank];
  int *tmp_dface_vtx_n = (int * ) malloc( dn_face_end * sizeof(int));
  // PDM_g_num_t *tmp_dface_vtx   = NULL;

  // PDM_block_to_block_exch(btb_face,
  //                         sizeof(int),//sizeof(PDM_g_num_t),
  //                         PDM_STRIDE_VAR_INTERLACED,
  //                         -1,
  //                         face_vtx_n,
  //               (void *)  face_vtx,
  //                         tmp_dface_vtx_n,
  //               (void **) &tmp_dface_vtx);

  // PDM_log_trace_array_int   (face_vtx_n    , n_face, "face_vtx_n : ");
  // PDM_log_trace_array_double(face_vtx_coord, 3 * 3 * n_face, "face_vtx_coord : ");

  double *dface_vtx_coord = NULL;
  PDM_block_to_block_exch(btb_face,
                          3 * sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          -1,
                          face_vtx_n,
                (void *)  face_vtx_coord,
                          tmp_dface_vtx_n,
                (void **) &dface_vtx_coord);


  PDM_block_to_block_free(btb_face);

  // *distrib_vtx  = _distrib_vtx;
  *distrib_face = _distrib_face;

  int *tmp_dface_vtx_idx = (int * ) malloc( (dn_face_end + 1) * sizeof(int));
  tmp_dface_vtx_idx[0] = 0;
  for(int i = 0; i < dn_face_end; ++i) {
    tmp_dface_vtx_idx[i+1] = tmp_dface_vtx_idx[i] + tmp_dface_vtx_n[i];
  }

  // PDM_log_trace_array_long(tmp_dface_vtx_idx, dn_face_end+1, "tmp_dface_vtx_idx ::");

  /*
   * Create gnum
   */
  PDM_gen_gnum_t* gen_gnum = PDM_gnum_create(3, 1, PDM_TRUE, 1.e-6, comm, PDM_OWNERSHIP_USER);

  double *char_length = malloc( tmp_dface_vtx_idx[dn_face_end] * sizeof(double));

  for (int i = 0; i < tmp_dface_vtx_idx[dn_face_end]; ++i) {
    char_length[i] = HUGE_VAL;//1.e-6;
  }

  double tol = 1e-6;
  const double eps_base = 1e-12;
  for(int i = 0; i < dn_face_end; ++i) {
    // double tmp_char_lenght = 1e30;     /* Big value */
    int n_vtx_on_face = tmp_dface_vtx_idx[i+1] - tmp_dface_vtx_idx[i];
    double *fvc = dface_vtx_coord + tmp_dface_vtx_idx[i]*3;
    double *chl = char_length     + tmp_dface_vtx_idx[i];
    for (int j = 0; j < n_vtx_on_face; j++) {
      int ivtx1 = j;
      int ivtx2 = (j+1) % n_vtx_on_face;

      double length2 = 0.;
      for (int k = 0; k < 3; k++) {
        double delta = fvc[3*ivtx1 + k] - fvc[3*ivtx2 + k];
        length2 += delta * delta;
      }

      chl[ivtx1] = PDM_MIN(chl[ivtx1], length2);
      chl[ivtx2] = PDM_MIN(chl[ivtx2], length2);
    }
    // for(int idx_vtx = tmp_dface_vtx_idx[i]; idx_vtx < tmp_dface_vtx_idx[i+1]; ++idx_vtx) {
    //   int idx_vtx1 = dface_vtx[idx_vtx] -1;
    //   double x1 = dface_vtx_coord[3*idx_vtx1  ];
    //   double y1 = dface_vtx_coord[3*idx_vtx1+1];
    //   double z1 = dface_vtx_coord[3*idx_vtx1+2];
    //   int idx_vtx2 = dface_vtx[(idx_vtx +1) % n_vtx_on_face] -1;
    //   double x2 = dface_vtx_coord[3*idx_vtx2  ];
    //   double y2 = dface_vtx_coord[3*idx_vtx2+1];
    //   double z2 = dface_vtx_coord[3*idx_vtx2+2];
    // }
  }

  for (int i = 0; i < tmp_dface_vtx_idx[dn_face_end]; ++i) {
    char_length[i] = PDM_MAX(eps_base, tol*sqrt(char_length[i]));
  }

  PDM_gnum_set_from_coords(gen_gnum,
                           0,
                           tmp_dface_vtx_idx[dn_face_end],
                           dface_vtx_coord,
                           char_length);

  PDM_gnum_compute(gen_gnum);
  free(char_length);

  PDM_g_num_t* _dface_vtx = PDM_gnum_get(gen_gnum, 0);
  PDM_gnum_free(gen_gnum);

  // PDM_log_trace_array_long(_dface_vtx, tmp_dface_vtx_idx[dn_face_end], "_dface_vtx");

  /*
   * Unified vtx
   */
  PDM_part_to_block_t* ptb_vtx = PDM_part_to_block_create(PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                          PDM_PART_TO_BLOCK_POST_CLEANUP,
                                                          1.,
                                                          &_dface_vtx,
                                                          NULL,
                                                          &tmp_dface_vtx_idx[dn_face_end],
                                                          1,
                                                          comm);

  PDM_g_num_t* _distrib_vtx = PDM_part_to_block_distrib_index_get(ptb_vtx);

  if (debug) {
    PDM_log_trace_array_long(init_distrib_vtx, n_rank+1, "init_distrib_vtx : ");
    PDM_log_trace_array_long(_distrib_vtx    , n_rank+1, "_distrib_vtx : ");
  }

  // int    *dvtx_coord_n  = NULL;
  double *_dvtx_coord   = NULL;
  PDM_part_to_block_exch(ptb_vtx,
                         3 * sizeof(double),
                         PDM_STRIDE_CST_INTERLACED,
                         1,
                         NULL,
                (void**) &dface_vtx_coord,
                         NULL,
                (void**) &_dvtx_coord);

  free(dface_vtx_coord);

  if (debug) {
    char outfilename[999];
    sprintf(outfilename, "debug_stl_dvtx_coord_%i.vtk", i_rank);
    int dn_vtx = _distrib_vtx[i_rank+1] - _distrib_vtx[i_rank];
    PDM_vtk_write_point_cloud(outfilename,
                              dn_vtx,
                              _dvtx_coord,
                              NULL,
                              NULL);
  }

  PDM_g_num_t *_tmp_distrib_vtx = malloc( (n_rank+1) * sizeof(PDM_g_num_t));
  for(int i = 0; i < n_rank+1; ++i) {
    _tmp_distrib_vtx[i] = _distrib_vtx[i];
  }

  PDM_part_to_block_free(ptb_vtx);

  *dvtx_coord   = _dvtx_coord;
  *dface_vtx    = _dface_vtx;

  *distrib_vtx  = _tmp_distrib_vtx;
  *distrib_face = _distrib_face;


  // free(_dvtx_coord);
  // free(_dface_vtx);

  free(tmp_dface_vtx_idx);
  // free(tmp_dface_vtx);

  free(tmp_dface_vtx_n  );
  // free(tmp_dvtx_coord  );

  free(init_distrib_vtx );
  free(init_distrib_face);

  if(i_rank == 0) {
    free(face_normal  );
    free(vtx_coord    );
    free(face_vtx_idx );
    free(face_vtx     );
    free(face_vtx_n   );
    free(face_vtx_coord);
  }
}


static
PDM_dmesh_nodal_t *
_create_dmesh_nodal
(
 PDM_MPI_Comm  comm,
 PDM_g_num_t  *distrib_vtx,
 PDM_g_num_t  *distrib_face,
 double       *dvtx_coord,
 PDM_g_num_t  *dface_vtx
)
{
  int i_rank, n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  /*
   *  Create dmesh nodal
   */
  PDM_g_num_t gn_vtx  = distrib_vtx [n_rank];
  PDM_g_num_t gn_face = distrib_face[n_rank];
  int dn_vtx  = distrib_vtx [i_rank+1] - distrib_vtx [i_rank];
  int dn_face = distrib_face[i_rank+1] - distrib_face[i_rank];

  PDM_dmesh_nodal_t *dmn = PDM_DMesh_nodal_create(comm,
                                                  2,
                                                  gn_vtx,
                                                  0,  // gn_cell
                                                  gn_face,
                                                  0); // gn_edge

  PDM_DMesh_nodal_coord_set(dmn,
                            dn_vtx,
                            dvtx_coord,
                            PDM_OWNERSHIP_KEEP);

  // dmn->surfacic->n_g_elmts = gn_face;
  int id_section = PDM_DMesh_nodal_section_add(dmn,
                                               PDM_GEOMETRY_KIND_SURFACIC,
                                               PDM_MESH_NODAL_TRIA3);

  PDM_DMesh_nodal_section_std_set(dmn,
                                  PDM_GEOMETRY_KIND_SURFACIC,
                                  id_section,
                                  dn_face,
                                  dface_vtx,
                                  PDM_OWNERSHIP_KEEP);


  PDM_dmesh_nodal_generate_distribution(dmn);

  return dmn;
}


/*=============================================================================
 * Public function prototypes
 *============================================================================*/


/**
 *
 * \brief Create a dmesh nodal from a file in ASCII STL mesh format
 *
 * \param[in]  comm                MPI communicator
 * \param[in]  filename            Filename
 *
 * \return Pointer to PDM_dmesh_nodal object
 *
 */

PDM_dmesh_nodal_t *
PDM_reader_stl_dmesh_nodal
(
 PDM_MPI_Comm   comm,
 const char    *filename
 )
{
  /* Read distributed mesh */
  int          *dface_vtx_idx = NULL;
  PDM_g_num_t  *dface_vtx     = NULL;
  double       *dvtx_coord    = NULL;
  PDM_g_num_t  *distrib_vtx   = NULL;
  PDM_g_num_t  *distrib_face  = NULL;

  _read_distributed_stl(comm,
                        filename,
                        &dface_vtx_idx,
                        &dface_vtx,
                        &dvtx_coord,
                        &distrib_vtx,
                        &distrib_face);

  /* Split mesh */
  /* !!! We assume we only have triangles */
  if (dface_vtx_idx != NULL){
    free(dface_vtx_idx);
  }

  PDM_dmesh_nodal_t *dmn = _create_dmesh_nodal(comm,
                                               distrib_vtx,
                                               distrib_face,
                                               dvtx_coord,
                                               dface_vtx);
  free(distrib_vtx);
  free(distrib_face);

  return dmn;
}
