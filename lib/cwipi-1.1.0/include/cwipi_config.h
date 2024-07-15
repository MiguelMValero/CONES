#ifndef __CWP_CONFIG_H__
#define __CWP_CONFIG_H__

/*------------------------------------------------------------------------------
 * Version (CWIPI_* variable are deprecated)
 *-----------------------------------------------------------------------------*/

#define CWP_VERSION       "1.1.0"
#define CWIPI_VERSION       "1.1.0"
#define CWP_DEF_VERSION       "1.1.0"

#define CWP_VERSION_MAJOR 1
#define CWP_VERSION_MINOR 1
#define CWP_VERSION_PATCH 0

#define CWP_MAJOR_VERSION 1
#define CWP_MINOR_VERSION 1
#define CWP_RELEASE_VERSION 0

/* #undef HAVE_MPI */

/* #undef HAVE_MPI_IO */

/* #undef HAVE_MPI_ONE_SIDED */

/* #undef HAVE_SPACE_BASIS */

/* #undef CWP_HAVE_FORTRAN_MPI_MODULE */

/* #undef CWP_HAVE_BLAS */

/* #undef CWP_HAVE_LAPACK */

/* #undef CWP_HAVE_MKL */

#define CWP_HAVE_NOT_FORTRAN_IN_C 1

#define CWP_MESH_DIR "/home/miguel/cwipi-1.1.0_1/cwipi-1.1.0/build/tests/meshes/"


#endif /* __CWP_CONFIG_H__ */
