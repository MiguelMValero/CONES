/*
  This file is part of the CWIPI library.

  Copyright (C) 2021-2023  ONERA

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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "cwipi.h"
#include "cwp.h"
#include "cwp_priv.h"

#include "pdm_dmesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_elmts.h"
#include "pdm_poly_surf_gen.h"
#include "pdm_part.h"
#include "pdm_mpi_node_first_rank.h"
#include "pdm_error.h"
#include "pdm_timer.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"

#include "pdm_multipart.h"
#include "pdm_gnum.h"

#include "pdm_array.h"
#include "pdm_writer.h"
#include "pdm_printf.h"
#include "pdm_logging.h"
#include "pdm_error.h"

#include "pdm_part_extension.h"
#include "pdm_vtk.h"

#include "pdm_generate_mesh.h"
#include "pdm_reader_gamma.h"


#define ABS(a) ((a) <  0  ? -(a) : (a))

/*----------------------------------------------------------------------
 *
 * Display usage
 *
 * parameters:
 *   exit code           <-- Exit code
 *---------------------------------------------------------------------*/

static void
_usage(int exit_code) {
  printf("\n"
         "  Usage: \n\n"
         "  -n           <> Number of vertices in band length.\n\n"
         "  -no_random      Disable mesh randomization\n\n"
         "  -n_proc_data <> Number of processes where there are data \n\n"
         "  -h              this message.\n\n");

  exit(exit_code);
}


/*----------------------------------------------------------------------
 *
 * Read args from the command line
 *
 * parameters:
 *   nVertex             <-- Number of vertices in bandwidth
 *   randLevel           <-- Random level
 *---------------------------------------------------------------------*/

static void
_read_args
(
  int                    argc,
  char                 **argv,
  int                   *n_vtx_seg1,
  int                   *n_vtx_seg2,
  double                *length,
  double                *separation_x,
  double                *separation_y,
  double                *separation_z,
  int                   *deform,
  double                *tolerance,
  int                   *randomize,
  int                   *nProcData,
  PDM_split_dual_t      *part_method,
  char                 **output_filename,
  int                   *verbose,
  int                   *use_gnum,
  int                   *n_part,
  int                   *n_block,
  char                 **filename
)
{
  int i = 1;

  // Parse and check command line
  while (i < argc) {
    if (strcmp(argv[i], "-h") == 0) {
      _usage(EXIT_SUCCESS);
    }
    else if (strcmp(argv[i], "-n") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n1") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg1 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n2") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_vtx_seg2 = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-length") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *length = atof(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-sep") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepx") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_x = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepy") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_y = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-sepz") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *separation_z = atof(argv[i]);
    }

    else if (strcmp(argv[i], "-output") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *output_filename = argv[i];
      }
    }
    else if (strcmp(argv[i], "-def") == 0) {
      *deform = 1;
    }
    else if (strcmp(argv[i], "-no_random") == 0) {
      *randomize = 0;
    }
    else if (strcmp(argv[i], "-n_proc_data") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *nProcData = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-t") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else
        *tolerance = atof(argv[i]);
    }
    else if (strcmp(argv[i], "-pt-scotch") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
    }
    else if (strcmp(argv[i], "-parmetis") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_PARMETIS;
    }
    else if (strcmp(argv[i], "-hilbert") == 0) {
      *part_method = PDM_SPLIT_DUAL_WITH_HILBERT;
    }
    else if (strcmp(argv[i], "-v") == 0) {
      *verbose = 1;
    }
    else if (strcmp(argv[i], "-no_gnum") == 0) {
      *use_gnum = 0;
    }
    else if (strcmp(argv[i], "-n_part") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_part = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-n_block") == 0) {
      i++;
      if (i >= argc) {
        _usage(EXIT_FAILURE);
      }
      else {
        *n_block = atoi(argv[i]);
      }
    }
    else if (strcmp(argv[i], "-f") == 0) {
      i++;
      if (i >= argc)
        _usage(EXIT_FAILURE);
      else {
        *filename = argv[i];
      }
    }
    else
      _usage(EXIT_FAILURE);
    i++;
  }
}

static void
_read_mesh
(
 const int                  current_rank_has_mesh,
 const PDM_MPI_Comm         comm,
 const int                  n_part,
 const PDM_split_dual_t     part_method,
 const char                *filename,
       int                 *n_block,
       int               ***n_elt,
       int              ****elt_vtx,
       PDM_g_num_t      ****elt_g_num,
       CWP_Block_t        **elt_type,
       int                **n_vtx,
       double            ***vtx_coord,
       PDM_g_num_t       ***vtx_g_num
 )
{
  assert(current_rank_has_mesh); // TODO : bcast n_block + block_type?

  *n_vtx     = malloc(sizeof(int          ) * n_part);
  *vtx_coord = malloc(sizeof(double      *) * n_part);
  *vtx_g_num = malloc(sizeof(PDM_g_num_t *) * n_part);

  PDM_dmesh_nodal_t *dmn = PDM_reader_gamma_dmesh_nodal(comm,
                                                        filename,
                                                        0,
                                                        0);

  int n_zone = 1;
  PDM_multipart_t *mpart = PDM_multipart_create(n_zone,
                                                &n_part,
                                                PDM_FALSE,
                                                part_method,
                                                PDM_PART_SIZE_HOMOGENEOUS,
                                                NULL,
                                                comm,
                                                PDM_OWNERSHIP_KEEP);

  PDM_multipart_set_reordering_options(mpart,
                                       -1,
                                       "PDM_PART_RENUM_CELL_NONE",
                                       NULL,
                                       "PDM_PART_RENUM_FACE_NONE");

  PDM_multipart_dmesh_nodal_set(mpart, 0, dmn);
  PDM_multipart_compute(mpart);


  PDM_part_mesh_nodal_t *pmn = NULL;
  PDM_multipart_get_part_mesh_nodal(mpart,
                                    0,
                                    &pmn,
                                    PDM_OWNERSHIP_USER);

  PDM_part_mesh_nodal_elmts_t *pmne = PDM_part_mesh_nodal_part_mesh_nodal_elmts_get(pmn,
                                                                                    PDM_GEOMETRY_KIND_VOLUMIC);

  const int *sections_id = PDM_part_mesh_nodal_elmts_sections_id_get(pmne);

  int _n_block = PDM_part_mesh_nodal_elmts_n_section_get(pmne);
  *n_block = _n_block;



  *n_elt          = malloc(sizeof(int          *) * _n_block);
  *elt_vtx        = malloc(sizeof(int         **) * _n_block);
  *elt_g_num      = malloc(sizeof(PDM_g_num_t **) * _n_block);
  *elt_type       = malloc(sizeof(CWP_Block_t  *) * _n_block);


  for (int iblock = 0; iblock < _n_block; iblock++) {
    (*n_elt)    [iblock] = PDM_array_zeros_int(n_part);
    (*elt_vtx)  [iblock] = malloc(sizeof(int         *) * n_part);
    (*elt_g_num)[iblock] = malloc(sizeof(PDM_g_num_t *) * n_part);

    (*elt_type)[iblock] = (CWP_Block_t) PDM_part_mesh_nodal_elmts_section_type_get(pmne,
                                                                                   sections_id[iblock]);
  }


  for (int ipart = 0; ipart < n_part; ipart++) {

    (*n_vtx)[ipart]         = PDM_part_mesh_nodal_n_vtx_get    (pmn, ipart);
    double      *_vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, ipart);
    PDM_g_num_t *_vtx_g_num = PDM_part_mesh_nodal_vtx_g_num_get(pmn, ipart);

    (*vtx_coord)[ipart] = malloc(sizeof(double) * (*n_vtx)[ipart] * 3);
    memcpy((*vtx_coord)[ipart],
           _vtx_coord,
           sizeof(double) * (*n_vtx)[ipart] * 3);

    (*vtx_g_num)[ipart] = malloc(sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);
    memcpy((*vtx_g_num)[ipart],
           _vtx_g_num,
           sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);


    for (int iblock = 0; iblock < _n_block; iblock++) {
      int         *connec              = NULL;
      PDM_g_num_t *numabs              = NULL;
      int         *parent_num          = NULL;
      PDM_g_num_t *parent_entity_g_num = NULL;
      PDM_part_mesh_nodal_elmts_section_std_get(pmne,
                                                sections_id[iblock],
                                                ipart,
                                                &connec,
                                                &numabs,
                                                &parent_num,
                                                &parent_entity_g_num,
                                                PDM_OWNERSHIP_USER);

      (*n_elt)[iblock][ipart] = PDM_part_mesh_nodal_elmts_section_n_elt_get(pmne,
                                                                            sections_id[iblock],
                                                                            ipart);

      if ((*elt_type)[iblock] == CWP_BLOCK_CELL_PRISM6) {
        // only for mesh from vtk ...
        for (int i = 0; i < (*n_elt)[iblock][ipart]; i++) {

          int tmp = connec[6*i+2];
          connec[6*i+2] = connec[6*i+1];
          connec[6*i+1] = tmp;

          tmp = connec[6*i+5];
          connec[6*i+5] = connec[6*i+4];
          connec[6*i+4] = tmp;
        }
      }


      (*elt_vtx)  [iblock][ipart] = connec;
      (*elt_g_num)[iblock][ipart] = numabs;

      if (parent_num != NULL) {
        free(parent_num);
      }
      if (parent_entity_g_num != NULL) {
        free(parent_entity_g_num);
      }
    }

  }

  PDM_part_mesh_nodal_free(pmn);
  PDM_multipart_free(mpart);
  PDM_DMesh_nodal_free(dmn);
}



static void
_gen_mesh
(
 const int                   active_rank,
 const PDM_MPI_Comm          comm,
 const int                   n_part,
 const int                   n_block,
 const PDM_split_dual_t      part_method,
 const PDM_g_num_t           n,
 const double                xmin,
 const double                ymin,
 const double                zmin,
 const double                length,
       int                ***n_elt,
       int               ****elt_vtx,
       PDM_g_num_t       ****elt_g_num,
       int                 **n_vtx,
       double             ***vtx_coord,
       PDM_g_num_t        ***vtx_g_num
 )
{
  PDM_Mesh_nodal_elt_t elt_type = PDM_MESH_NODAL_HEXA8;
  int elt_vtx_n = PDM_Mesh_nodal_n_vtx_elt_get(elt_type, 1);

  *n_vtx     = malloc(sizeof(int          ) * n_part);
  *vtx_coord = malloc(sizeof(double      *) * n_part);
  *vtx_g_num = malloc(sizeof(PDM_g_num_t *) * n_part);


  *n_elt     = malloc(sizeof(int          *) * n_block);
  *elt_vtx   = malloc(sizeof(int         **) * n_block);
  *elt_g_num = malloc(sizeof(PDM_g_num_t **) * n_block);

  if (active_rank) {
    PDM_part_mesh_nodal_t *pmn = PDM_generate_mesh_parallelepiped(comm,
                                                                  elt_type,
                                                                  1,
                                                                  NULL,
                                                                  xmin,
                                                                  ymin,
                                                                  zmin,
                                                                  length,
                                                                  length,
                                                                  length,
                                                                  n,
                                                                  n,
                                                                  n,
                                                                  n_part,
                                                                  part_method);


    for (int iblock = 0; iblock < n_block; iblock++) {
      (*n_elt)    [iblock] = PDM_array_zeros_int(n_part);
      (*elt_vtx)  [iblock] = malloc(sizeof(int         *) * n_part);
      (*elt_g_num)[iblock] = malloc(sizeof(PDM_g_num_t *) * n_part);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {

      int part_n_elt = PDM_part_mesh_nodal_section_n_elt_get(pmn, 0, ipart);

      int         *connec              = NULL;
      PDM_g_num_t *numabs              = NULL;
      int         *parent_num          = NULL;
      PDM_g_num_t *parent_entity_g_num = NULL;
      PDM_part_mesh_nodal_section_std_get(pmn,
                                          0,
                                          ipart,
                                          &connec,
                                          &numabs,
                                          &parent_num,
                                          &parent_entity_g_num,
                                          PDM_OWNERSHIP_KEEP);

      int block_n_elt = part_n_elt / n_block;

      for (int iblock = 0; iblock < n_block-1; iblock++) {
        (*n_elt)[iblock][ipart] = block_n_elt;
      }

      (*n_elt)[n_block-1][ipart] = part_n_elt - (n_block-1)*block_n_elt;


      int idx = 0;
      for (int iblock = 0; iblock < n_block; iblock++) {
        (*elt_vtx)[iblock][ipart] = malloc(sizeof(int) * elt_vtx_n * (*n_elt)[iblock][ipart]);

        memcpy((*elt_vtx)[iblock][ipart],
               connec + elt_vtx_n * idx,
               sizeof(int) * elt_vtx_n * (*n_elt)[iblock][ipart]);

        (*elt_g_num)[iblock][ipart] = malloc(sizeof(PDM_g_num_t) * (*n_elt)[iblock][ipart]);
        memcpy((*elt_g_num)[iblock][ipart],
               numabs + idx,
               sizeof(PDM_g_num_t) * (*n_elt)[iblock][ipart]);

        idx += (*n_elt)[iblock][ipart];
      }

      (*n_vtx)[ipart]         = PDM_part_mesh_nodal_n_vtx_get(pmn, ipart);
      double      *_vtx_coord = PDM_part_mesh_nodal_vtx_coord_get(pmn, ipart);
      PDM_g_num_t *_vtx_g_num = PDM_part_mesh_nodal_vtx_g_num_get(pmn, ipart);

      (*vtx_coord)[ipart] = malloc(sizeof(double) * (*n_vtx)[ipart] * 3);
      memcpy((*vtx_coord)[ipart],
             _vtx_coord,
             sizeof(double) * (*n_vtx)[ipart] * 3);

      (*vtx_g_num)[ipart] = malloc(sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);
      memcpy((*vtx_g_num)[ipart],
             _vtx_g_num,
             sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);

    }

    PDM_part_mesh_nodal_free(pmn);
  }
  else {

    for (int iblock = 0; iblock < n_block; iblock++) {
      (*n_elt)    [iblock] = PDM_array_zeros_int(n_part);
      (*elt_vtx)  [iblock] = malloc(sizeof(int         *) * n_part);
      (*elt_g_num)[iblock] = malloc(sizeof(PDM_g_num_t *) * n_part);
    }

    for (int ipart = 0; ipart < n_part; ipart++) {
      for (int iblock = 0; iblock < n_block; iblock++) {
        (*elt_vtx)  [iblock][ipart] = malloc(sizeof(int)         * elt_vtx_n * (*n_elt)[iblock][ipart]);
        (*elt_g_num)[iblock][ipart] = malloc(sizeof(PDM_g_num_t) * elt_vtx_n);
      }

      (*n_vtx)    [ipart] = 0;
      (*vtx_coord)[ipart] = malloc(sizeof(double     ) * (*n_vtx)[ipart] * 3);
      (*vtx_g_num)[ipart] = malloc(sizeof(PDM_g_num_t) * (*n_vtx)[ipart]);
    }

  }
}



/*----------------------------------------------------------------------
 *
 * Main : surface coupling test : P1P0_P0P1
 *
 *---------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
  // Read args from command line
  int n_part                      = 1;
  int n_block                     = 2;
  int n_vtx_seg1                  = 4;
  int n_vtx_seg2                  = 4;
  int randomize                   = 1;
  int n_proc_data                 = -1;

#ifdef PDM_HAVE_PARMETIS
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PARMETIS;
#else
#ifdef PDM_HAVE_PTSCOTCH
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_PTSCOTCH;
#else
  PDM_split_dual_t part_method    = PDM_SPLIT_DUAL_WITH_HILBERT;
#endif
#endif
  int verbose                     = 0;

  double length                   = 20.;
  int    deform                   = 0;

  double      separation_x        = 2.;
  double      separation_y        = 0.;
  double      separation_z        = 0.;

  double      tolerance           = 1e-2;

  char* output_filename           = NULL;
  // int filedump                    = 0;

  int         use_gnum            = 1;
  char       *mesh_filename       = NULL;

  _read_args(argc,
             argv,
             &n_vtx_seg1,
             &n_vtx_seg2,
             &length,
             &separation_x,
             &separation_y,
             &separation_z,
             &deform,
             &tolerance,
             &randomize,
             &n_proc_data,
             &part_method,
             &output_filename,
             &verbose,
             &use_gnum,
             &n_part,
             &n_block,
             &mesh_filename);

  // if (output_filename !=NULL) {
  //   filedump = 1;
  // }

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int i_rank;
  int comm_world_size;

  MPI_Comm_rank(MPI_COMM_WORLD, &i_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_world_size);

  assert (comm_world_size > 1);

  if (n_proc_data == 1) {
    n_proc_data = 2;
  }




  // Initialize CWIPI
  int n_code = 1;
  int code_id;
  const char **code_name = malloc(sizeof(char *) * n_code);
  const char **coupled_code_name = malloc(sizeof(char *) * n_code);
  CWP_Status_t is_active_rank = CWP_STATUS_ON;

  int n_vtx_seg;
  if (i_rank < comm_world_size / 2) {
    code_id = 1;
    code_name[0] = "code1";
    coupled_code_name[0] = "code2";
    n_vtx_seg = n_vtx_seg1;
  }
  else {
    code_id = 2;
    code_name[0] = "code2";
    coupled_code_name[0] = "code1";
    n_vtx_seg = n_vtx_seg2;
  }



  MPI_Comm *intra_comm = malloc(sizeof(MPI_Comm) * n_code);

  CWP_Init(MPI_COMM_WORLD,
           n_code,
           (const char **) code_name,
           is_active_rank,
           intra_comm);


  if (verbose && i_rank == 0) {
    printf("CWIPI Init OK\n");
  }


  // Create coupling
  const char *coupling_name = "c_new_api_multiblock";

  CWP_Cpl_create(code_name[0],
                 coupling_name,
                 coupled_code_name[0],
                 CWP_INTERFACE_VOLUME,
                 CWP_COMM_PAR_WITH_PART,
                 CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_OCTREE,
                 // CWP_SPATIAL_INTERP_FROM_LOCATION_MESH_LOCATION_LOCATE_ALL_TGT,
                 n_part,
                 CWP_DYNAMIC_MESH_STATIC,
                 CWP_TIME_EXCH_USER_CONTROLLED);

  CWP_Visu_set(code_name[0], coupling_name, 1, CWP_VISU_FORMAT_ENSIGHT, "text");

  if (verbose && i_rank == 0) {
    printf("Create coupling OK\n");
  }

  // Define mesh
  PDM_MPI_Comm mesh_comm = PDM_MPI_mpi_2_pdm_mpi_comm((void *) intra_comm);

  // int _n_proc_data = n_proc_data;
  // if (n_proc_data > 0) {
  //   if (code_id == 1) {
  //     _n_proc_data /= 2;
  //   }
  //   else {
  //     _n_proc_data -= n_proc_data / 2;
  //   }
  // }
  int current_rank_has_mesh = 1;//_set_rank_has_mesh(intra_comm[0], _n_proc_data, &mesh_comm);

  // int true_n_proc_data;
  // MPI_Reduce(&current_rank_has_mesh, &true_n_proc_data, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  // if (i_rank == 0) {
  //   printf("nb procs with mesh data = %d\n", true_n_proc_data);
  // }


  double xmin = -0.5 * length;
  double ymin = -0.5 * length;
  double zmin = -0.5 * length;
  //  int init_random = (int) time(NULL);
  int init_random = 5;

  PDM_MPI_Comm code_mesh_comm;
  PDM_MPI_Comm_split(mesh_comm, code_id, i_rank, &code_mesh_comm);

  if (code_id == 2) {
    init_random++;
    xmin += separation_x;
    ymin += separation_y;
    zmin += separation_z;
  }

  int          **n_elt          = NULL;
  int         ***elt_vtx        = NULL;
  PDM_g_num_t ***elt_g_num      = NULL;
  CWP_Block_t   *elt_type       = NULL;
  int           *n_vtx          = NULL;
  double       **vtx_coord      = NULL;
  PDM_g_num_t  **vtx_g_num      = NULL;
  if (mesh_filename == NULL) {
    _gen_mesh(current_rank_has_mesh,
              code_mesh_comm,
              n_part,
              n_block,
              part_method,
              n_vtx_seg,
              xmin,
              ymin,
              zmin,
              length,
              &n_elt,
              &elt_vtx,
              &elt_g_num,
              &n_vtx,
              &vtx_coord,
              &vtx_g_num);

    elt_type = malloc(sizeof(CWP_Block_t) * n_block);
    for (int i = 0; i < n_block; i++) {
      elt_type[i] = CWP_BLOCK_CELL_HEXA8;
    }
  }
  else {
    _read_mesh(current_rank_has_mesh,
               code_mesh_comm,
               n_part,
               part_method,
               (const char *) mesh_filename,
               &n_block,
               &n_elt,
               &elt_vtx,
               &elt_g_num,
               &elt_type,
               &n_vtx,
               &vtx_coord,
               &vtx_g_num);

    if (code_id == 2) {
      for (int ipart = 0; ipart < n_part; ipart++) {
        for (int i = 0; i < n_vtx[ipart]; i++) {
          vtx_coord[ipart][3*i  ] += separation_x;
          vtx_coord[ipart][3*i+1] += separation_y;
          vtx_coord[ipart][3*i+2] += separation_z;
        }
      }
    }
  }

  if (0) {
    for (int iblock = 0; iblock < n_block; iblock++) {
      for (int ipart = 0; ipart < n_part; ipart++) {
        char filename[999];
        sprintf(filename, "%s_block%d_part%d_rank%d.vtk", code_name[0], iblock, ipart, i_rank);
        PDM_vtk_write_std_elements(filename,
                                   n_vtx    [ipart],
                                   vtx_coord[ipart],
                                   vtx_g_num[ipart],
            (PDM_Mesh_nodal_elt_t) elt_type [iblock],
                                   n_elt    [iblock][ipart],
                                   elt_vtx  [iblock][ipart],
                                   elt_g_num[iblock][ipart],
                                   0,
                                   NULL,
                                   NULL);
      }
    }
  }


  // Set interface mesh
  PDM_g_num_t *_vtx_g_num  = NULL;
  PDM_g_num_t *_elt_g_num = NULL;

  int *block_id = malloc(sizeof(int) * n_block);
  for (int iblock = 0; iblock < n_block; iblock++) {
    block_id[iblock] = CWP_Mesh_interf_block_add(code_name[0],
                                                 coupling_name,
                                                 elt_type[iblock]);
  }

  for (int ipart = 0; ipart < n_part; ipart++) {

    if (use_gnum) {
      _vtx_g_num = vtx_g_num[ipart];
    }

    CWP_Mesh_interf_vtx_set(code_name[0],
                            coupling_name,
                            ipart,
                            n_vtx[ipart],
                            vtx_coord[ipart],
                            _vtx_g_num);

    for (int iblock = 0; iblock < n_block; iblock++) {

      if (use_gnum) {
        _elt_g_num = elt_g_num[iblock][ipart];
      }

      CWP_Mesh_interf_block_std_set(code_name[0],
                                    coupling_name,
                                    ipart,
                                    block_id[iblock],
                                    n_elt   [iblock][ipart],
                                    elt_vtx [iblock][ipart],
                                    _elt_g_num);
    }
  }

  CWP_Mesh_interf_finalize(code_name[0],
                           coupling_name);

  if (verbose && i_rank == 0) {
    printf("Set mesh OK\n");
  }


  // Create and set fields
  double **send_val1 = NULL;
  double **recv_val1 = NULL;
  double **send_val2 = NULL;
  double **recv_val2 = NULL;

  const char *field_name1 = "coo";
  const char *field_name2 = "gnum";

  // TODO : cell-based field in multiblock??

  int *pn_cell = malloc(sizeof(int) * n_part);
  for (int ipart = 0; ipart < n_part; ipart++) {
    pn_cell[ipart] = 0;
    for (int iblock = 0; iblock < n_block; iblock++) {
      pn_cell[ipart] += n_elt[iblock][ipart];
    }
  }

  if (code_id == 1) {
    send_val1 = malloc(sizeof(double *) * n_part);
    send_val2 = malloc(sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      send_val1[ipart] = malloc(sizeof(double) * n_vtx[ipart] * 3);
      memcpy(send_val1[ipart], vtx_coord[ipart], sizeof(double) * n_vtx[ipart] * 3);

      send_val2[ipart] = malloc(sizeof(double) * pn_cell[ipart]);
      int idx = 0;
      for (int iblock = 0; iblock < n_block; iblock++) {
        for (int ielt = 0; ielt < n_elt[iblock][ipart]; ielt++) {
          send_val2[ipart][idx++] = (double) elt_g_num[iblock][ipart][ielt];
          // send_val2[ipart][idx++] = elt_type[iblock];
        }
      }
    }
  } else {
    recv_val1 = malloc(sizeof(double *) * n_part);
    recv_val2 = malloc(sizeof(double *) * n_part);
    for (int ipart = 0; ipart < n_part; ipart++) {
      recv_val1[ipart] = malloc(sizeof(double) * n_vtx[ipart] * 3);
      recv_val2[ipart] = malloc(sizeof(double) * pn_cell[ipart]);//n_vtx[ipart]);
    }
  }

  CWP_Status_t visu_status = CWP_STATUS_ON;
  MPI_Barrier(MPI_COMM_WORLD);

  if (code_id == 1) {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     3,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,
                     CWP_FIELD_EXCH_SEND,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    for (int ipart = 0; ipart < n_part; ipart++) {
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name1,
                         ipart,
                         CWP_FIELD_MAP_SOURCE,
                         send_val1[ipart]);

      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         ipart,
                         CWP_FIELD_MAP_SOURCE,
                         send_val2[ipart]);
    }
  } else {
    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name1,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     3,
                     CWP_DOF_LOCATION_NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Field_create(code_name[0],
                     coupling_name,
                     field_name2,
                     CWP_DOUBLE,
                     CWP_FIELD_STORAGE_INTERLACED,
                     1,
                     CWP_DOF_LOCATION_CELL_CENTER,//NODE,
                     CWP_FIELD_EXCH_RECV,
                     visu_status);

    CWP_Time_step_beg(code_name[0],
                      0.0);

    for (int ipart = 0; ipart < n_part; ipart++) {
      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name1,
                         ipart,
                         CWP_FIELD_MAP_TARGET,
                         recv_val1[ipart]);

      CWP_Field_data_set(code_name[0],
                         coupling_name,
                         field_name2,
                         ipart,
                         CWP_FIELD_MAP_TARGET,
                         recv_val2[ipart]);
    }
  }

  if (verbose && i_rank == 0) {
    printf("Fields OK\n");
  }

  // Perform geometric algorithm
  CWP_Spatial_interp_property_set(code_name[0], coupling_name, "tolerance", CWP_DOUBLE, "0.01");
  CWP_Spatial_interp_weights_compute(code_name[0], coupling_name);

  //  Exchange interpolated fields
  if (code_id == 1) {
    CWP_Field_issend(code_name[0], coupling_name, field_name1);
    CWP_Field_issend(code_name[0], coupling_name, field_name2);
  }
  else {
    CWP_Field_irecv(code_name[0], coupling_name, field_name1);
    CWP_Field_irecv(code_name[0], coupling_name, field_name2);
  }

  if (code_id == 1) {
    CWP_Field_wait_issend(code_name[0], coupling_name, field_name1);
    CWP_Field_wait_issend(code_name[0], coupling_name, field_name2);
  }
  else {
    CWP_Field_wait_irecv(code_name[0], coupling_name, field_name1);
    CWP_Field_wait_irecv(code_name[0], coupling_name, field_name2);
  }

  //  Check
  double max_err = 0.;
  int    n_err   = 0;
  if (code_id == 2) {

    for (int ipart = 0; ipart < n_part; ipart++) {
      int        n_located = CWP_N_computed_tgts_get(code_name[0], coupling_name, field_name1, ipart);
      const int *located   = CWP_Computed_tgts_get  (code_name[0], coupling_name, field_name1, ipart);

      if (0) {
        log_trace("recv_val / coord = \n");
        for (int i = 0 ; i < n_located; i++) {
          int ivtx = located[i] - 1;
          log_trace("%d ("PDM_FMT_G_NUM"): %f %f %f / %f %f %f\n",
                    located[i],
                    vtx_g_num[ipart][located[i]-1],
                    recv_val1[ipart][3*i], recv_val1[ipart][3*i+1], recv_val1[ipart][3*i+2],
                    vtx_coord[ipart][3*ivtx], vtx_coord[ipart][3*ivtx+1], vtx_coord[ipart][3*ivtx+2]);
        }
      }

      for (int i = 0 ; i < n_located ; i++) {
        int wrong = 0;

        for (int j = 0; j < 3; j++) {
          double err = ABS (recv_val1[ipart][3*i + j] - vtx_coord[ipart][3 * (located[i] -1) + j]);
          if (err > 1.e-4) {
            wrong = 1;
            printf("[%d] !! vtx "PDM_FMT_G_NUM" %d err = %g (coord#%d = %f, recv = %f)\n",
                   i_rank, vtx_g_num[ipart][(located[i] - 1)],
                    located[i], err, j,
                    vtx_coord[ipart][3*(located[i]-1) + j], recv_val1[ipart][3*i + j]);
          }
          if (err > max_err) {
            max_err = err;
          }
        }

        if (wrong) {
          n_err++;
        }
      }
    }

  }

  double global_max_err = 0.;
  MPI_Reduce(&max_err, &global_max_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("Max error = %g\n", global_max_err);
  }

  int global_n_err = 0.;
  MPI_Reduce(&n_err, &global_n_err, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
  if (i_rank == 0) {
    printf("N error = %d\n", global_n_err);
  }

  CWP_Time_step_end(code_name[0]);

  //  Delete interface mesh
  CWP_Mesh_interf_del(code_name[0], coupling_name);

  //  Delete coupling
  CWP_Cpl_del(code_name[0], coupling_name);

  // Free memory
  free(code_name);
  free(coupled_code_name);
  free(intra_comm);

  if (current_rank_has_mesh) {
    for (int iblock = 0; iblock < n_block; iblock++) {
      for (int ipart = 0; ipart < n_part; ipart++) {
        free(elt_vtx  [iblock][ipart]);
        free(elt_g_num[iblock][ipart]);
      }
      free(n_elt    [iblock]);
      free(elt_vtx  [iblock]);
      free(elt_g_num[iblock]);
    }
    for (int ipart = 0; ipart < n_part; ipart++) {
      free(vtx_g_num[ipart]);
      free(vtx_coord[ipart]);
    }
  }
  free(n_elt    );
  free(elt_vtx  );
  free(elt_g_num);
  free(elt_type );
  free(n_vtx    );
  free(vtx_coord);
  free(vtx_g_num);
  free(pn_cell  );

  free(block_id);


  if (code_id == 1) {
    for (int ipart = 0; ipart < n_part; ipart++) {
      free(send_val1[ipart]);
      free(send_val2[ipart]);
    }
    free(send_val1);
    free(send_val2);
  }
  else {
    for (int ipart = 0; ipart < n_part; ipart++) {
      free(recv_val1[ipart]);
      free(recv_val2[ipart]);
    }
    free(recv_val1);
    free(recv_val2);
  }

  //  Finalize CWIPI
  CWP_Finalize();

  // Finalize MPI
  MPI_Finalize();

  return EXIT_SUCCESS;
}
