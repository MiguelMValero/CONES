#ifndef __PDM_CONFIG_H__
#define __PDM_CONFIG_H__

! #undef PDM_HAVE_FORTRAN_MPI_MODULE
! #undef PDM_HAVE_PARMETIS
#define PDM_HAVE_OPENMP
! #undef PDM_HAVE_CUDA
! #undef PDM_HAVE_ANISO_AGGLO
#define PDM_HAVE_PTSCOTCH
#define PDM_HAVE_GETRUSAGE
#define PDM_HAVE_GETTIMEOFDAY
#define PDM_LONG_G_NUM
! #undef PDM_HAVE_LAPACK
! #undef PDM_HAVE_BLAS
! #undef PDM_HAVE_MKL
! #undef PDM_USE_MULTIPART

#define PDM_LONG_G_NUM_BOOL ""
#define PDM_VERSION_MAJOR "2"
#define PDM_VERSION_MINOR "4"
#define PDM_VERSION_PATCH "1"
#define PDM_VERSION       "2.4.1"
#define PDM_MESH_DIR      "/home/miguel/cwipi-1.1.0_1/cwipi-1.1.0/build/external/paradigm/test/meshes/"
#define PDM_MESH_DIR_F    "/home/miguel/cwipi-1.1.0_1/cwipi-1.1.0/build/external/paradigm"

#endif /*__PDM_CONFIG_H__*/
