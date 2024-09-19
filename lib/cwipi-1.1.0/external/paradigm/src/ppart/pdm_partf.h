#include "pdmf.h"

  integer, parameter :: PDM_part_SPLIT_PARMETIS = 1
  integer, parameter :: PDM_part_SPLIT_PTSCOTCH = 2
  integer, parameter :: PDM_part_SPLIT_HILBERT = 3

  integer, parameter :: PDM_PART_RENUM_FACE_RANDOM        = 1
  integer, parameter :: PDM_PART_RENUM_FACE_NONE          = 2
  integer, parameter :: PDM_PART_RENUM_FACE_LEXICOGRAPHIC = 3

  integer, parameter :: PDM_part_RENUM_CELL_HILBERT = 1
  integer, parameter :: PDM_part_RENUM_CELL_RANDOM = 2
  integer, parameter :: PDM_part_RENUM_CELL_NONE = 3
  integer, parameter :: PDM_part_RENUM_CELL_CUTHILL = 4

interface

 !================================================================================
 !
 ! \brief Build a initial partitioning
 !
 !  Build a initial partitioning from :
 !      - Cell block distribution with implicit global numbering
 !         (the first cell is the first cell of the first process and
 !          the latest cell is the latest cell of the latest process)
 !      - Face block distribution with implicit global numbering
 !      - Vertex block distribution with implicit global numbering
 !  To repart an existing partition use \ref PDM_part_repart function
 !
 ! \param [out]  ppart_id        ppart identifier
 ! \param [in]   pt_comm        MPI Comminicator
 ! \param [in]   method         Choice between (1 for ParMETIS or 2 for PT-Scotch)
 ! \param [in]   n_part          Number of partition to build on this process
 ! \param [in]   dn_cell         Number of distributed cells
 ! \param [in]   dn_face         Number of distributed faces
 ! \param [in]   dn_vtx          Number of distributed vertices
 ! \param [in]   n_face_group     Number of face groups
 ! \param [in]   dcell_faceIdx   Distributed cell face connectivity index or NULL
 !                              (size : dn_cell + 1, numbering : 0 to n-1)
 ! \param [in]   dcell_face      Distributed cell face connectivity or NULL
 !                              (size : dface_vtx_idx[dn_cell], numbering : 1 to n)
 ! \param [in]   dcell_tag       Cell tag (size : n_cell) or NULL
 ! \param [in]   dcell_weight    Cell weight (size : n_cell) or NULL
 ! \param [in]   dcell_part      Distributed cell partitioning
 !                              (size = dn_cell) or NULL (No partitioning if != NULL)
 ! \param [in]   dface_cell      Distributed face cell connectivity or NULL
 !                              (size : 2 * dn_face, numbering : 1 to n)
 ! \param [in]   dface_vtx_idx    Distributed face to vertex connectivity index
 !                              (size : dn_face + 1, numbering : 0 to n-1)
 ! \param [in]   dface_vtx       Distributed face to vertex connectivity
 !                              (size : dface_vtx_idx[dn_face], numbering : 1 to n)
 ! \param [in]   dface_tag       Distributed face tag (size : dn_face)
 !                              or NULL
 ! \param [in]   dvtx_coord      Distributed vertex coordinates
 !                              (size : 3!dn_vtx)
 ! \param [in]   dvtx_tag        Distributed vertex tag (size : dn_vtx) or NULL
 ! \param [in]   dface_group_idx  Index of distributed faces list of each group
 !                              (size = n_face_group + 1) or NULL
 ! \param [in]   dface_group     distributed faces list of each group
 !                              (size = dface_group[dface_group_idx[n_face_group]],
 !                              numbering : 1 to n)
 !                              or NULL
 !================================================================================

    subroutine pdm_part_create(ppart_id, &
                          pt_comm, &
                          method,  &
                          n_part, &
                          dn_cell, &
                          dn_face, &
                          dn_vtx,&
                          n_face_group, &
                          have_dcell_face,&
                          dcell_faceIdx,&
                          dcell_face,&
                          have_dcell_tag,&
                          dcell_tag, &
                          have_dcell_weight,&
                          dcell_weight,&
                          have_dcell_part,&
                          dcell_part,&
                          have_dface_cell,&
                          dface_cell,&
                          dface_vtx_idx,&
                          dface_vtx,&
                          have_dface_tag,&
                          dface_tag,&
                          dvtx_coord,&
                          have_dvtx_tag,&
                          dvtx_tag,&
                          dface_group_idx,&
                          dface_group)

    implicit none

    integer                     ::  ppart_id
    integer                     ::  pt_comm
    integer                     ::  method
    integer                     ::  n_part
    integer                     ::  dn_cell
    integer                     ::  dn_face
    integer                     ::  dn_vtx
    integer                     ::  n_face_group
    integer                     ::  have_dcell_face
    integer                     ::  dcell_faceIdx(*)
    integer (kind = PDM_g_num_s) :: dcell_face(*)
    integer                     ::  have_dcell_tag
    integer                     ::  dcell_tag(*)
    integer                     ::  have_dcell_weight
    integer                     ::  dcell_weight(*)
    integer                     ::  have_dcell_part
    integer                     ::  dcell_part(*)
    integer                     ::  have_dface_cell
    integer (kind = PDM_g_num_s) :: dface_cell(*)
    integer                     ::  dface_vtx_idx(*)
    integer (kind = PDM_g_num_s) :: dface_vtx(*)
    integer                     ::  have_dface_tag
    integer                     ::  dface_tag(*)
    double precision            ::  dvtx_coord(*)
    integer                     ::  have_dvtx_tag
    integer                     ::  dvtx_tag(*)
    integer                     ::  dface_group_idx(*)
    integer (kind = PDM_g_num_s) :: dface_group(*)

  end subroutine


 !==============================================================================
 !
 ! \brief Return a mesh partition dimensions
 !
 ! \param [in]   ppart_id            ppart identifier
 ! \param [in]   i_part              Current partition
 ! \param [out]  n_cell              Number of cells
 ! \param [out]  n_face              Number of faces
 ! \param [out]  n_face_part_bound     Number of partitioning boundary faces
 ! \param [out]  n_vtx               Number of vertices
 ! \param [out]  n_proc              Number of processus
 ! \param [out]  n_total_part             Total number of partitions
 ! \param [out]  scell_face          Size of cell-face connectivity
 ! \param [out]  sface_vtx           Size of face-vertex connectivity
 ! \param [out]  sFacePartBound     Size of face_part_bound array
 ! \param [out]  sface_group         Size of face_group array
 !
 !==============================================================================


   subroutine pdm_part_part_dim_get (ppart_id, &
                                  i_part, &
                                  n_cell, &
                                  n_face, &
                                  n_face_part_bound, &
                                  n_vtx, &
                                  n_proc, &
                                  n_total_part, &
                                  scell_face, &
                                  sface_vtx, &
                                  sface_group)
     implicit none

     integer :: ppart_id
     integer :: i_part
     integer :: n_cell
     integer :: n_face
     integer :: n_face_part_bound
     integer :: n_vtx
     integer :: n_proc
     integer :: n_total_part
     integer :: scell_face
     integer :: sface_vtx
     integer :: sface_group

   end subroutine pdm_part_part_dim_get


 !==============================================================================
 !
 ! \brief Return a mesh partition
 !
 ! \param [in]   ppart_id            ppart identifier
 ! \param [in]   i_part              Current partition
 ! \param [out]  cell_tag            Cell tag (size = n_cell)
 ! \param [out]  cell_face_idx        Cell to face connectivity index
 !                                  (size = n_cell + 1, numbering : 0 to n-1)
 ! \param [out]  cell_face           Cell to face connectivity
 !                                 (size = cell_face_idx[n_cell] = lcell_face
 !                                  numbering : 1 to n)
 ! \param [out]  cell_ln_to_gn         Cell local numbering to global numbering
 !                                  (size = n_cell, numbering : 1 to n)
 ! \param [out]  face_tag            Face tag (size = n_face)
 ! \param [out]  face_cell           Face to cell connectivity
 !                                  (size = 2 * n_face, numbering : 1 to n)
 ! \param [out]  face_vtx_idx         Face to Vertex connectivity index
 !                                  (size = n_face + 1, numbering : 0 to n-1)
 ! \param [out]  face_vtx            Face to Vertex connectivity
 !                                  (size = faceVertexIdx[n_face], numbering : 1 to n)
 ! \param [out]  face_ln_to_gn         Face local numbering to global numbering
 !                                  (size = n_face, numbering : 1 to n)
 ! \param [out]  face_part_bound_proc_idx Partitioning boundary faces block
 !                                    distribution from processus
 !                                    (size = n_proc + 1)
 ! \param [out]  face_part_bound_part_idx Partitioning boundary faces block
 !                                    distribution from partition
 !                                    (size = n_total_part + 1)
 ! \param [out]  face_part_bound      Partitioning boundary faces
 !                                  (size = 4 * n_face_part_bound)
 !                                    sorted by processus,
 !                                    sorted by partition in each processus, and
 !                                    sorted by absolute face number in each partition
 !                                  For each face :
 !                                     - Face local number
 !                                       (numbering : 1 to n)
 !                                     - Connected process
 !                                       (numbering : 0 to n-1)
 !                                     - Connected Partition
 !                                       on the connected process
 !                                       (numbering :1 to n)
 !                                     - Connected face local number
 !                                       in the connected partition
 !                                       (numbering :1 to n)
 ! \param [out]  vtx_tag             Vertex tag (size = nVertex)
 ! \param [out]  vtx                Vertex coordinates (size = 3 * nVertex)
 ! \param [out]  vtx_ln_to_gn          Vertex local numbering to global numbering
 !                                  (size = n_vtx, numbering : 1 to n)
 ! \param [out]  face_group_idx       Face group index
 !                                  (size = n_face_group + 1, numbering : 1 to n-1)
 ! \param [out]  face_group          faces for each group
 !                                  (size = face_group_idx[n_face_group] = lFaceGroup,
 !                                   numbering : 1 to n)
 ! \param [out]  face_group_ln_to_gn    Faces global numbering for each group
 !                                  (size = face_group_idx[n_face_group] = lFaceGroup,
 !                                  numbering : 1 to n)
 !==============================================================================

 subroutine pdm_part_part_val_get (ppart_id, &
                                i_part, &
                                cell_tag, &
                                cell_face_idx, &
                                cell_face, &
                                cell_ln_to_gn, &
                                face_tag, &
                                face_cell, &
                                face_vtx_idx, &
                                face_vtx, &
                                face_ln_to_gn, &
                                face_part_bound_proc_idx, &
                                face_part_bound_part_idx, &
                                face_part_bound, &
                                vtx_tag, &
                                vtx, &
                                vtx_ln_to_gn, &
                                face_group_idx, &
                                face_group, &
                                face_group_ln_to_gn)

   implicit none

   integer                       :: ppart_id
   integer                       :: i_part
   integer                       :: cell_tag(*)
   integer                       :: cell_face_idx(*)
   integer                       :: cell_face(*)
   integer (kind = PDM_g_num_s)  :: cell_ln_to_gn(*)
   integer                       :: face_tag(*)
   integer                       :: face_cell(*)
   integer                       :: face_vtx_idx(*)
   integer                       :: face_vtx(*)
   integer (kind = PDM_g_num_s)  :: face_ln_to_gn(*)
   integer                       :: face_part_bound_proc_idx(*)
   integer                       :: face_part_bound_part_idx(*)
   integer                       :: face_part_bound(*)
   integer                       :: vtx_tag(*)
   double precision              :: vtx(*)
   integer (kind = PDM_g_num_s)  :: vtx_ln_to_gn(*)
   integer                       :: face_group_idx(*)
   integer                       :: face_group(*)
   integer (kind = PDM_g_num_s)  :: face_group_ln_to_gn(*)

 end subroutine pdm_part_part_val_get

 !==============================================================================
 !
 ! \brief Free ppart
 !
 ! \param [in]   ppart_id        ppart identifier
 !
 !==============================================================================

 subroutine pdm_part_free (ppart_id)

   implicit none

   integer                       :: ppart_id

 end subroutine pdm_part_free

 !==============================================================================
 !
 ! \brief Return times
 !
 ! \param [in]   ppart_id     ppart identifier
 ! \param [out]  elapsed     elapsed times (size = 4)
 ! \param [out]  cpu         cpu times (size = 4)
 ! \param [out]  cpu_user    user cpu times (size = 4)
 ! \param [out]  cpu_sys     system cpu times (size = 4)
 !
 !==============================================================================

 subroutine pdm_part_time_get (ppart_id, &
                            elapsed, &
                            cpu, &
                            cpu_user, &
                            cpu_sys)
   implicit none

   integer           :: ppart_id
   double precision  :: elapsed
   double precision  :: cpu
   double precision  :: cpu_user
   double precision  :: cpu_sys

 end subroutine pdm_part_time_get

 !==============================================================================
 !
 ! \brief Return statistic
 !
 ! \param [in]   ppart_id                        ppart identifier
 ! \param [out]  cells_average                  average of cells number
 ! \param [out]  cells_median                   median of cells number
 ! \param [out]  cells_std_deviation            standard deviation of cells number
 ! \param [out]  cells_min                      minimum of cells nummber
 ! \param [out]  cells_max                      maximum of cells nummber
 ! \param [out]  bound_part_faces_average       average of partitioning boundary faces
 ! \param [out]  bound_part_faces_median        median of partitioning boundary faces
 ! \param [out]  bound_part_faces_std_deviation standard deviation of partitioning boundary faces
 ! \param [out]  bound_part_faces_min           minimum of partitioning boundary faces
 ! \param [out]  bound_part_faces_max           maximum of partitioning boundary faces
 !
 !==============================================================================

 subroutine pdm_part_stat_get (ppart_id,  &
                            cells_average, &
                            cells_median, &
                            cells_std_deviation, &
                            cells_min,  &
                            cells_max, &
                            bound_part_faces_average, &
                            bound_part_faces_median,  &
                            bound_part_faces_std_deviation, &
                            bound_part_faces_min, &
                            bound_part_faces_max, &
                            bound_part_faces_sum)

   implicit none

   integer      :: ppart_id
   integer      :: cells_average
   integer      :: cells_median
   double precision   :: cells_std_deviation
   integer      :: cells_min
   integer      :: cells_max
   integer      :: bound_part_faces_average
   integer      :: bound_part_faces_median
   double precision   :: bound_part_faces_std_deviation
   integer      :: bound_part_faces_min
   integer      :: bound_part_faces_max
   integer      :: bound_part_faces_sum

 end subroutine pdm_part_stat_get


end interface
