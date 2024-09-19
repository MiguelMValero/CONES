/*
 * \file
 */

#ifndef __PDM_MPI_H__
#define __PDM_MPI_H__

/*============================================================================
 * Bibliotheque de messagerie
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include "pdm.h"
/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

typedef int PDM_MPI_Request;
typedef int PDM_MPI_Win;
typedef int PDM_MPI_Comm;
typedef int PDM_MPI_Datatype;
typedef int PDM_MPI_File;

typedef long long PDM_MPI_Offset;

typedef long PDM_MPI_Aint;
typedef int PDM_MPI_Fint;

#define PDM_MPI_UNDEFINED MPI_UNDEFINED
#define PDM_MPI_IN_PLACE  MPI_IN_PLACE


enum  {

  PDM_MPI_SUCCESS,
  PDM_MPI_ERR_BUFFER,
  PDM_MPI_ERR_COUNT,
  PDM_MPI_ERR_TYPE,
  PDM_MPI_ERR_TAG,
  PDM_MPI_ERR_COMM,
  PDM_MPI_ERR_RANK,
  PDM_MPI_ERR_ROOT,
  PDM_MPI_ERR_TRUNCATE,
  PDM_MPI_ERR_GROUP,
  PDM_MPI_ERR_OP,
  PDM_MPI_ERR_REQUEST,
  PDM_MPI_ERR_TOPOLOGY,
  PDM_MPI_ERR_DIMS,
  PDM_MPI_ERR_ARG,
  PDM_MPI_ERR_UNKNOWN,
  PDM_MPI_ERR_OTHER,
  PDM_MPI_ERR_INTERN,
  PDM_MPI_ERR_IN_STATUS,
  PDM_MPI_ERR_PENDING,
  PDM_MPI_MAX_ERROR_STRING,
  PDM_MPI_ERR_ACCESS,
  PDM_MPI_ERR_AMODE,
  PDM_MPI_ERR_BAD_FILE,
  PDM_MPI_ERR_CONVERSION,
  PDM_MPI_ERR_DUP_DATAREP,
  PDM_MPI_ERR_FILE_EXISTS,
  PDM_MPI_ERR_FILE_IN_USE,
  PDM_MPI_ERR_FILE,
  PDM_MPI_ERR_INFO_KEY,
  PDM_MPI_ERR_INFO_NOKEY,
  PDM_MPI_ERR_INFO_VALUE,
  PDM_MPI_ERR_IO,
  PDM_MPI_ERR_NO_MEM,
  PDM_MPI_ERR_NOT_SAME,
  PDM_MPI_ERR_NO_SPACE,
  PDM_MPI_ERR_NO_SUCH_FILE,
  PDM_MPI_ERR_QUOTA,
  PDM_MPI_ERR_READ_ONLY,
  PDM_MPI_ERR_UNSUPPORTED_DATAREP,
  PDM_MPI_ERR_UNSUPPORTED_OPERATION,
  PDM_MPI_ERR_WIN,
  PDM_MPI_ERR_LASTCODE,
  PDM_MPI_ERR_ASSERT,
  PDM_MPI_ERR_BASE,
  PDM_MPI_ERR_DISP,
  PDM_MPI_ERR_KEYVAL,
  PDM_MPI_ERR_LOCKTYPE,
  PDM_MPI_ERR_RMA_CONFLICT,
  PDM_MPI_ERR_RMA_SYNC,
  PDM_MPI_ERR_SIZE

};

enum {

  PDM_MPI_MODE_CREATE,
  PDM_MPI_MODE_RDONLY,
  PDM_MPI_MODE_WRONLY,
  PDM_MPI_MODE_RDWR ,
  PDM_MPI_MODE_DELETE_ON_CLOSE,
  PDM_MPI_MODE_UNIQUE_OPEN,
  PDM_MPI_MODE_EXCL,
  PDM_MPI_MODE_APPEND,
  PDM_MPI_MODE_SEQUENTIAL,
  PDM_MPI_DISPLACEMENT_CURRENT,
  PDM_MPI_SEEK_SET,
  PDM_MPI_SEEK_CUR,
  PDM_MPI_SEEK_END,
  PDM_MPI_MODE_WRONLY_APPEND,
  PDM_MPI_MODE_WRONLY_CREATE

};

typedef enum {

  PDM_MPI_MAX,
  PDM_MPI_MIN,
  PDM_MPI_SUM,
  PDM_MPI_MINLOC,
  PDM_MPI_MAXLOC,
  PDM_MPI_OP_NULL

} PDM_MPI_Op;

enum {
  PDM_MPI_BYTE                   = -1,
  PDM_MPI_PACKED                 = -2,
  PDM_MPI_CHAR                   = -3,
  PDM_MPI_SHORT                  = -4,
  PDM_MPI_INT                    = -5,
  PDM_MPI_LONG                   = -6,
  PDM_MPI_FLOAT                  = -7,
  PDM_MPI_DOUBLE                 = -8,
  PDM_MPI_LONG_DOUBLE            = -9,
  PDM_MPI_UNSIGNED_CHAR          = -10,
  PDM_MPI_UNSIGNED_SHORT         = -11,
  PDM_MPI_UNSIGNED_LONG          = -12,
  PDM_MPI_UNSIGNED               = -13,
  PDM_MPI_FLOAT_INT              = -14,
  PDM_MPI_DOUBLE_INT             = -15,
  PDM_MPI_LONG_DOUBLE_INT        = -16,
  PDM_MPI_LONG_INT               = -17,
  PDM_MPI_SHORT_INT              = -18,
  PDM_MPI_2INT                   = -19,
  PDM_MPI_CHARACTER              = -20,
  PDM_MPI_INTEGER                = -21,
  PDM_MPI_REAL                   = -22,
  PDM_MPI_DOUBLE_PRECISION       = -23,
  PDM_MPI_DATATYPE_NULL          = -24,
  PDM_MPI_INT8_T                 = -25,
  PDM_MPI_INT16_T                = -26,
  PDM_MPI_INT32_T                = -27,
  PDM_MPI_INT64_T                = -28,
  PDM_MPI_UINT8_T                = -29,
  PDM_MPI_UINT16_T               = -30,
  PDM_MPI_UINT32_T               = -31,
  PDM_MPI_UINT64_T               = -32,
  PDM_MPI_UNSIGNED_LONG_LONG     = -33,
};

enum {
  PDM_MPI_COMM_NULL  = -1,
  PDM_MPI_COMM_WORLD = -2
};

enum {
  PDM_MPI_FILE_NULL = -1
};

enum {
  PDM_MPI_REQUEST_NULL  = -123456789
};

enum {
  PDM_MPI_WIN_NULL  = -1
};

enum {
  PDM_MPI_SPLIT_SHARED  = 1,
  PDM_MPI_SPLIT_NUMA    = 2
};

typedef struct _pdm_mpi_win_shared_t PDM_mpi_win_shared_t;

/*============================================================================
 * Prototype des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Init(int *argc, char ***argv);

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Finalize(void);

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_mpi_2_pdm_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm(void *);

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_free_mpi_comm(void *mpi_comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_open (wrapping de la fonction MPI_File_open)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_open(PDM_MPI_Comm comm, char *filename, int amode,
                  PDM_MPI_File *fh);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_close (wrapping de la fonction MPI_File_close)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_close(PDM_MPI_File *fh);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_seek (wrapping de la fonction MPI_File_seek)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_seek(PDM_MPI_File, PDM_MPI_Offset, int);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_size (wrapping de la fonction MPI_File_get_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_size(PDM_MPI_File, PDM_MPI_Offset *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_position (wrapping de la fonction MPI_File_get_position)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_position(PDM_MPI_File, PDM_MPI_Offset *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_set_view (wrapping de la fonction MPI_File_set_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_set_view(PDM_MPI_File, PDM_MPI_Offset, PDM_MPI_Datatype,
                      PDM_MPI_Datatype, const char *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_view (wrapping de la fonction MPI_File_get_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_view(PDM_MPI_File, PDM_MPI_Offset *,
                 PDM_MPI_Datatype *, PDM_MPI_Datatype *, char *);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at (wrapping de la fonction MPI_File_read_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                     int count, PDM_MPI_Datatype datatype, int *n_octet_lus);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at_all (wrapping de la fonction MPI_File_read_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                         int count, PDM_MPI_Datatype datatype, int *n_octet_lus);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at (wrapping de la fonction MPI_File_write_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                      int count, PDM_MPI_Datatype datatype, int *n_octet_ecrits);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at_all (wrapping de la fonction MPI_File_write_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_ecrits);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read (wrapping de la fonction MPI_File_read)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read(PDM_MPI_File fh, void *buf, int count,
                  PDM_MPI_Datatype datatype, int *n_octet_lus);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_all (wrapping de la fonction MPI_File_read_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_all(PDM_MPI_File fh, void *buf, int count,
                      PDM_MPI_Datatype datatype, int *n_octet_lus);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write (wrapping de la fonction MPI_File_write)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write(PDM_MPI_File fh, void *buf, int count,
                   PDM_MPI_Datatype datatype, int *n_octet_ecrits);

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_all (wrapping de la fonction MPI_File_write_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_all(PDM_MPI_File fh, void *buf, int count,
                       PDM_MPI_Datatype datatype, int *n_octet_ecrits);

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Igather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm, PDM_MPI_Request *reqauest);

/*----------------------------------------------------------------------------
 * PDM_MPI_Gatherv (wrapping de la fonction MPI_Gatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int *recvcounts, int *displs,
                PDM_MPI_Datatype recvtype, int root, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Recv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
             int tag, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Irecv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
              int tag, PDM_MPI_Comm comm, PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Send (wrapping de la fonction MPI_Send)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Send(void *buf, int count, PDM_MPI_Datatype datatype, int dest,
             int tag, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Issend (wrapping de la fonction MPI_Issend)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Issend(const void *buf, int count, PDM_MPI_Datatype datatype, int dest, int tag,
               PDM_MPI_Comm comm, PDM_MPI_Request *request);


/*----------------------------------------------------------------------------
 * PDM_MPI_Wait (wrapping de la fonction MPI_Wait)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Wait(PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Test (wrapping de la fonction MPI_Test)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Test(PDM_MPI_Request *request, int *flag);

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_hindexed (wrapping de la fonction MPI_Type_hindexed)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_create_hindexed(int count, const int array_of_blocklengths[],
                      const PDM_MPI_Aint array_of_displacements[],
                      PDM_MPI_Datatype oldtype, PDM_MPI_Datatype *newtype);

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_create_contiguous (wrapping de la fonction MPI_Type_create_contiguous)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_create_contiguous(int               count,
                                   PDM_MPI_Datatype  oldtype,
                                   PDM_MPI_Datatype *newtype);
/*----------------------------------------------------------------------------
 * PDM_MPI_Type_commit (wrapping de la fonction MPI_Type_commit)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_commit(PDM_MPI_Datatype *datatype);

/*----------------------------------------------------------------------------
 * MPI_Type_size (wrapping de la fonction MPI_Type_commit)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_size(PDM_MPI_Datatype datatype, int *size);

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_free (wrapping de la fonction MPI_Type_free)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_free(PDM_MPI_Datatype *datatype);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_f2c (wrapping de la fonction MPI_comm_f2c)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_Comm_f2c(PDM_MPI_Fint comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_c2f (wrapping de la fonction MPI_comm_c2f)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Fint PDM_MPI_Comm_c2f(PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Scatter (wrapping de la fonction MPI_Scatter)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scatter(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
                int root, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Barrier (wrapping de la fonction MPI_Barrier)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Barrier(PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Wtime (wrapping de la fonction MPI_Barrier)
 *
 *----------------------------------------------------------------------------*/

double PDM_MPI_Wtime(void);

/*----------------------------------------------------------------------------
 * PDM_MPI_Bcast (wrapping de la fonction MPI_Bcast)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Bcast(void *buffer, int count, PDM_MPI_Datatype datatype,
              int root, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_IBcast (wrapping de la fonction MPI_IBcast)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ibcast(void *buffer, int count, PDM_MPI_Datatype datatype,
                   int root, PDM_MPI_Comm comm, PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                  void *recvbuf, int recvcount,
                  PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgatherv (wrapping de la fonction MPI_Allgatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts,
                   int *displs, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce (wrapping de la fonction MPI_Reduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		   PDM_MPI_Datatype datatype, PDM_MPI_Op op,
		   int root, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce_scatter (wrapping de la fonction MPI_Reduce_scatter)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *counts,
                           PDM_MPI_Datatype datatype, PDM_MPI_Op op,
                           PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Allreduce (wrapping de la fonction MPI_Allreduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm);


/*----------------------------------------------------------------------------
 * PDM_MPI_Scan (wrapping de la fonction MPI_Scan)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm);

int PDM_MPI_Exscan(const void *sendbuf, void *recvbuf, int count,
                   PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm);

int PDM_MPI_Iscan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm,
             PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoall (wrapping de la fonction MPI_Alltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                 void *recvbuf, int recvcount,
                 PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Alltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                      void *recvbuf, int recvcount,
                      PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm, PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoallv (wrapping de la fonction MPI_Alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                  PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                  int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

int PDM_MPI_Alltoallv_l(void *sendbuf, int *sendcounts, size_t *sdispls,
                  PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                  size_t *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoallv (wrapping de la fonction MPI_Ialltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                  PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                  int *rdispls, PDM_MPI_Datatype recvtype,
                  PDM_MPI_Comm comm, PDM_MPI_Request *request);


/*----------------------------------------------------------------------------
 * PDM_MPI_Get_ialltoallv (Implemtation of alltoall like with window )
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Get_ialltoallv(PDM_MPI_Win       win_send,
                           PDM_MPI_Win       win_recv,
                           void             *sendbuf,
                           int              *sendcounts,
                           int              *sdispls,
                           PDM_MPI_Datatype  sendtype,
                           void             *recvbuf,
                           int              *recvcounts,
                           int              *rdispls,
                           PDM_MPI_Datatype  recvtype,
                           PDM_MPI_Comm      comm);


/*----------------------------------------------------------------------------
 * PDM_MPI_Win_allocate (wrapping de la fonction MPI_Win_allocate)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_allocate(PDM_MPI_Aint  size,
                         int           disp_unit,
                         PDM_MPI_Comm  comm,
                         void         *baseptr,
                         PDM_MPI_Win  *win);


/*----------------------------------------------------------------------------
 * PDM_MPI_Win_free (wrapping de la fonction MPI_Win_free)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_free(PDM_MPI_Win  *win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_fence (wrapping de la fonction MPI_Win_fence)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Win_fence(int assert, PDM_MPI_Win win);

/*----------------------------------------------------------------------------
 * PDM_MPI_Error_string (wrapping de la fonction MPI_Error_string)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Error_string(int errorcode, char *string, int *resultlen);


/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_rank(PDM_MPI_Comm comm, int *rank);


/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_size(PDM_MPI_Comm comm, int *size);


/*----------------------------------------------------------------------------
 * PDM_MPI_get_max_error_string
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_get_max_error_string(void);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_free
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_free(PDM_MPI_Comm *comm);


/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_split(PDM_MPI_Comm comm, int color, int key, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_dup
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_dup(PDM_MPI_Comm comm, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_split_type_numa(PDM_MPI_Comm comm, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_split_type(PDM_MPI_Comm comm, int split_type, PDM_MPI_Comm *newcomm);

/*----------------------------------------------------------------------------
 * PDM_mpi_Win_allocate_shared_get
 *
 *----------------------------------------------------------------------------*/
PDM_mpi_win_shared_t*
PDM_mpi_win_shared_create(PDM_MPI_Aint          size,
                          int                   disp_unit,
                          PDM_MPI_Comm          comm);

void* PDM_mpi_win_shared_get(PDM_mpi_win_shared_t *wins);

void PDM_mpi_win_shared_free(PDM_mpi_win_shared_t *wins);

PDM_MPI_Comm PDM_MPI_get_group_of_master(PDM_MPI_Comm comm, PDM_MPI_Comm sub_comm);

int PDM_mpi_win_shared_lock_all(int assert, PDM_mpi_win_shared_t* win);
int PDM_mpi_win_shared_unlock_all(PDM_mpi_win_shared_t* win);
int PDM_mpi_win_shared_sync(PDM_mpi_win_shared_t* win);

/*----------------------------------------------------------------------------
 * PDM_MPI_rand_tag_get
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_get_attr_tag_ub(PDM_MPI_Comm comm, void *attribute_val, int *flag);
int PDM_MPI_Rand_tag            (PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Dist_graph_create_adjacent
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Dist_graph_create_adjacent(PDM_MPI_Comm  comm_old,
                                             int     indegree,
                                       const int     sources[],
                                             int     outdegree,
                                       const int     destinations[],
                                       int           reorder,
                                       PDM_MPI_Comm *comm_dist_graph);


/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_allgather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                               void *recvbuf, int recvcount,
                               PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_allgatherv (wrapping de la fonction MPI_Neighbor_allgatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_allgatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                   void *recvbuf, int *recvcounts,
                   int *displs, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);


/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_alltoall (wrapping de la fonction MPI_Neighbor_alltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                              void *recvbuf, int recvcount,
                              PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Ialltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ineighbor_alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                               void *recvbuf, int recvcount,
                               PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm, PDM_MPI_Request *request);

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_alltoallv (wrapping de la fonction MPI_Neighbor_alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                               PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                               int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm);

/*----------------------------------------------------------------------------
 * PDM_MPI_Ineighbor_alltoallv (wrapping de la fonction MPI_Ineighbor_alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ineighbor_alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                                PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                                int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm,
                                PDM_MPI_Request *request);

void
PDM_MPI_setup_hybrid_dist_comm_graph
(
  PDM_MPI_Comm   comm,
  PDM_MPI_Comm  *comm_shared_out,
  PDM_MPI_Comm  *comm_dist_graph_out,
  int           *n_degree,
  int          **neighbor
);

int
PDM_MPI_Dist_graph_neighbors_count
(
  PDM_MPI_Comm  comm,
  int          *n_degree_in,
  int          *n_degree_out,
  int          *is_weighted
);

int
PDM_MPI_Dist_graph_neighbors
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *sources,
  int            n_degree_out,
  int           *destinations
);

void
PDM_MPI_setup_dist_graph_from_neighbor_in
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *neighbor_in,
  PDM_MPI_Comm  *comm_dist_graph_out
);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_MPI_H__ */
