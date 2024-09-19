/*============================================================================
 * Encapsulation de MPI
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <assert.h>
#include <sys/time.h>
#ifdef __linux__
#include <sys/syscall.h> //Non portable mettre un ifdef
#endif

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_priv.h"
#include "pdm_mpi_priv.h"
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/


/*============================================================================
 * Definition des variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Indirection sur le code d'erreur
 *----------------------------------------------------------------------------*/

static const int mpi_err[] = {

  MPI_SUCCESS,
  MPI_ERR_BUFFER,
  MPI_ERR_COUNT,
  MPI_ERR_TYPE,
  MPI_ERR_TAG,
  MPI_ERR_COMM,
  MPI_ERR_RANK,
  MPI_ERR_ROOT,
  MPI_ERR_TRUNCATE,
  MPI_ERR_GROUP,
  MPI_ERR_OP,
  MPI_ERR_REQUEST,
  MPI_ERR_TOPOLOGY,
  MPI_ERR_DIMS,
  MPI_ERR_ARG,
  MPI_ERR_UNKNOWN,
  MPI_ERR_OTHER,
  MPI_ERR_INTERN,
  MPI_ERR_IN_STATUS,
  MPI_ERR_PENDING,
  MPI_MAX_ERROR_STRING,



  MPI_ERR_ACCESS,
  MPI_ERR_AMODE,
  MPI_ERR_BAD_FILE,
  MPI_ERR_CONVERSION,
  MPI_ERR_DUP_DATAREP,
  MPI_ERR_FILE_EXISTS,
  MPI_ERR_FILE_IN_USE,
  MPI_ERR_FILE,
  MPI_ERR_INFO_KEY,
  MPI_ERR_INFO_NOKEY,
  MPI_ERR_INFO_VALUE,
  MPI_ERR_IO,
  MPI_ERR_NO_MEM,
  MPI_ERR_NOT_SAME,
  MPI_ERR_NO_SPACE,
  MPI_ERR_NO_SUCH_FILE,
  MPI_ERR_QUOTA,
  MPI_ERR_READ_ONLY,
  MPI_ERR_UNSUPPORTED_DATAREP,
  MPI_ERR_UNSUPPORTED_OPERATION,
  MPI_ERR_WIN,
  MPI_ERR_LASTCODE,
  MPI_ERR_ASSERT,
  MPI_ERR_BASE,
  MPI_ERR_DISP,
  MPI_ERR_KEYVAL,
  MPI_ERR_LOCKTYPE,
  MPI_ERR_RMA_CONFLICT,
  MPI_ERR_RMA_SYNC,
  MPI_ERR_SIZE

};

/*----------------------------------------------------------------------------
 * Indirection sur le mode du fichier
 *----------------------------------------------------------------------------*/

static const int mpi_file_mode[] = {

MPI_MODE_CREATE,
MPI_MODE_RDONLY,
MPI_MODE_WRONLY,
MPI_MODE_RDWR,
MPI_MODE_DELETE_ON_CLOSE,
MPI_MODE_UNIQUE_OPEN,
MPI_MODE_EXCL,
MPI_MODE_APPEND,
MPI_MODE_SEQUENTIAL,
MPI_DISPLACEMENT_CURRENT,
MPI_SEEK_SET,
MPI_SEEK_CUR,
MPI_SEEK_END,
MPI_MODE_WRONLY | MPI_MODE_APPEND,
MPI_MODE_WRONLY | MPI_MODE_CREATE

};

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Datatype -> MPI_Datatype
 *----------------------------------------------------------------------------*/

static const MPI_Datatype mpi_datatype_cste[] = {

  MPI_BYTE,
  MPI_PACKED,
  MPI_CHAR,
  MPI_SHORT,
  MPI_INT,
  MPI_LONG,
  MPI_FLOAT,
  MPI_DOUBLE,
  MPI_LONG_DOUBLE,
  MPI_UNSIGNED_CHAR,
  MPI_UNSIGNED_SHORT,
  MPI_UNSIGNED_LONG,
  MPI_UNSIGNED,
  MPI_FLOAT_INT,
  MPI_DOUBLE_INT,
  MPI_LONG_DOUBLE_INT,
  MPI_LONG_INT,
  MPI_SHORT_INT,
  MPI_2INT,
  MPI_CHARACTER,
  MPI_INTEGER,
  MPI_REAL,
  MPI_DOUBLE_PRECISION,
  MPI_DATATYPE_NULL,
  MPI_INT8_T,
  MPI_INT16_T,
  MPI_INT32_T,
  MPI_INT64_T,
  MPI_UINT8_T,
  MPI_UINT16_T,
  MPI_UINT32_T,
  MPI_UINT64_T,
  MPI_UNSIGNED_LONG_LONG
};

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Op -> MPI_Op
 *----------------------------------------------------------------------------*/

static const MPI_Op mpi_op[] = {

  MPI_MAX,
  MPI_MIN,
  MPI_SUM,
  MPI_MINLOC,
  MPI_MAXLOC,
  MPI_OP_NULL

};


/*----------------------------------------------------------------------------
 * Indirection constantes PDM_MPI_File ->constantes MPI_File
 *----------------------------------------------------------------------------*/

static const MPI_File mpi_file_cste[] = {

  MPI_FILE_NULL

};

/*----------------------------------------------------------------------------
 * Indirection constantes PDM_MPI_Comm ->constantes MPI_Comm
 *----------------------------------------------------------------------------*/

static const MPI_Comm mpi_comm_cste[] = {

  MPI_COMM_NULL,
  MPI_COMM_WORLD

};

/*----------------------------------------------------------------------------
 * Indirection constantes PDM_MPI_Request ->constantes MPI_Request
 *----------------------------------------------------------------------------*/

static const MPI_Request mpi_request_cste[] = {

  MPI_REQUEST_NULL,

};

/*----------------------------------------------------------------------------
 * Indirection constantes PDM_MPI_Request ->constantes MPI_Request
 *----------------------------------------------------------------------------*/

static const MPI_Win mpi_win_cste[] = {

  MPI_WIN_NULL,

};


/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_File -> MPI_File
 * stockage dans un tableau
 *----------------------------------------------------------------------------*/

static MPI_File **mpi_file   = NULL; /* Tableau de stockage */
static int       l_mpi_file = 0;     /* Taille du tableau */
static int       n_mpi_file = 0;     /* Nombre de fichiers stockes */

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Comm -> MPI_Comm
 * stockage dans un tableau
 *----------------------------------------------------------------------------*/

static MPI_Comm **mpi_comm   = NULL; /* Tableau de stockage */
static int       l_mpi_comm = 0;     /* Taille du tableau */
static int       n_mpi_comm = 0;     /* Nombre de communicateurs stockes */

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Request -> MPI_Request
 * stockage dans un tableau
 *----------------------------------------------------------------------------*/

static MPI_Request **mpi_request = NULL; /* Tableau de stockage */
static int       l_mpi_request = 0;   /* Taille du tableau */
static int       n_mpi_request = 0;   /* Nombre de communicateurs stockes */

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Win -> MPI_Win
 * stockage dans un tableau
 *----------------------------------------------------------------------------*/

static MPI_Win **mpi_win = NULL; /* Tableau de stockage */
static int       l_mpi_win = 0;   /* Taille du tableau */
static int       n_mpi_win = 0;   /* Nombre de communicateurs stockes */

/*----------------------------------------------------------------------------
 * Indirection PDM_MPI_Datatype -> MPI_Datatype
 * stockage dans un tableau des types utilisateurs
 *----------------------------------------------------------------------------*/

static MPI_Datatype **mpi_datatype   = NULL; /* Tableau de stockage */
static int           l_mpi_datatype = 0;     /* Taille du tableau */
static int           n_mpi_datatype = 0;     /* Nombre de communicateurs stockes */

/*============================================================================
 * Defintion des fonctions pprivees
 *============================================================================*/

/*----------------------------------------------------------------------------
 * mpi_err -> pdm_mpi_err
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

static int _mpi_2_pdm_mpi_err(int code_mpi)
{
  int code = PDM_MPI_ERR_OTHER;
  switch(code_mpi) {
  case MPI_SUCCESS:
    code = PDM_MPI_SUCCESS;
  break;
  case MPI_ERR_BUFFER:
    code = PDM_MPI_ERR_BUFFER;
      break;
  case MPI_ERR_COUNT:
    code = PDM_MPI_ERR_COUNT;
      break;
  case MPI_ERR_TYPE:
    code = PDM_MPI_ERR_TYPE;
    break;
  case MPI_ERR_TAG:
    code = PDM_MPI_ERR_TAG;
    break;
  case MPI_ERR_COMM:
    code = PDM_MPI_ERR_COMM;
    break;
  case MPI_ERR_RANK:
    code = PDM_MPI_ERR_RANK;
    break;
  case MPI_ERR_ROOT:
    code = PDM_MPI_ERR_ROOT;
    break;
  case MPI_ERR_TRUNCATE:
    code = PDM_MPI_ERR_TRUNCATE;
    break;
  case MPI_ERR_GROUP:
    code = PDM_MPI_ERR_GROUP;
    break;
  case MPI_ERR_OP:
    code = PDM_MPI_ERR_OP;
    break;
  case MPI_ERR_REQUEST:
    code = PDM_MPI_ERR_REQUEST;
    break;
  case MPI_ERR_TOPOLOGY:
    code = PDM_MPI_ERR_TOPOLOGY;
    break;
  case MPI_ERR_DIMS:
    code = PDM_MPI_ERR_DIMS;
    break;
  case MPI_ERR_ARG:
    code = PDM_MPI_ERR_ARG;
    break;
  case MPI_ERR_UNKNOWN:
    code = PDM_MPI_ERR_UNKNOWN;
    break;
  case MPI_ERR_OTHER:
    code = PDM_MPI_ERR_OTHER;
    break;
  case MPI_ERR_INTERN:
    code = PDM_MPI_ERR_INTERN;
    break;
  case MPI_ERR_IN_STATUS:
    code = PDM_MPI_ERR_IN_STATUS;
    break;
  case MPI_ERR_PENDING:
    code = PDM_MPI_ERR_PENDING;
    break;



  case MPI_ERR_ACCESS:
    code = PDM_MPI_ERR_ACCESS;
    break;
  case MPI_ERR_AMODE:
    code = PDM_MPI_ERR_AMODE;
    break;
  case MPI_ERR_BAD_FILE:
    code = PDM_MPI_ERR_BAD_FILE;
    break;
  case MPI_ERR_CONVERSION:
    code = PDM_MPI_ERR_CONVERSION;
    break;
  case MPI_ERR_DUP_DATAREP:
    code = PDM_MPI_ERR_DUP_DATAREP;
    break;
  case MPI_ERR_FILE_EXISTS:
    code = PDM_MPI_ERR_FILE_EXISTS;
    break;
  case MPI_ERR_FILE_IN_USE:
    code = PDM_MPI_ERR_FILE_IN_USE;
    break;
  case MPI_ERR_FILE:
    code = PDM_MPI_ERR_FILE;
    break;
  case MPI_ERR_INFO_KEY:
    code = PDM_MPI_ERR_INFO_KEY;
    break;
  case MPI_ERR_INFO_NOKEY:
    code = PDM_MPI_ERR_INFO_NOKEY;
    break;
  case MPI_ERR_INFO_VALUE:
    code = PDM_MPI_ERR_INFO_VALUE;
    break;
  case MPI_ERR_IO:
    code = PDM_MPI_ERR_IO;
    break;
  case MPI_ERR_NO_MEM:
    code = PDM_MPI_ERR_NO_MEM;
    break;
  case MPI_ERR_NOT_SAME:
    code = PDM_MPI_ERR_NOT_SAME;
    break;
  case MPI_ERR_NO_SPACE:
    code = PDM_MPI_ERR_NO_SPACE;
    break;
  case MPI_ERR_NO_SUCH_FILE:
    code = PDM_MPI_ERR_NO_SUCH_FILE;
    break;
  case MPI_ERR_QUOTA:
    code = PDM_MPI_ERR_QUOTA;
    break;
  case MPI_ERR_READ_ONLY:
    code = PDM_MPI_ERR_READ_ONLY;
    break;
  case MPI_ERR_UNSUPPORTED_DATAREP:
    code = PDM_MPI_ERR_UNSUPPORTED_DATAREP;
    break;
  case MPI_ERR_UNSUPPORTED_OPERATION:
    code = PDM_MPI_ERR_UNSUPPORTED_OPERATION;
    break;
  case MPI_ERR_WIN:
    code = PDM_MPI_ERR_WIN;
    break;
  case MPI_ERR_LASTCODE:
    code = PDM_MPI_ERR_LASTCODE;
    break;



  case MPI_ERR_ASSERT:
    code = PDM_MPI_ERR_ASSERT;
    break;
  case MPI_ERR_BASE:
    code = PDM_MPI_ERR_BASE;
    break;
  case MPI_ERR_DISP:
    code = PDM_MPI_ERR_DISP;
    break;
  case MPI_ERR_KEYVAL:
    code = PDM_MPI_ERR_KEYVAL;
    break;
  case MPI_ERR_LOCKTYPE:
    code = PDM_MPI_ERR_LOCKTYPE;
    break;
  case MPI_ERR_RMA_CONFLICT:
    code = PDM_MPI_ERR_RMA_CONFLICT;
    break;
  case MPI_ERR_RMA_SYNC:
    code = PDM_MPI_ERR_RMA_SYNC;
    break;
  case MPI_ERR_SIZE:
    code = PDM_MPI_ERR_SIZE;

  }
  return code;
}

/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

static MPI_Comm _pdm_mpi_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm)
{

  /* Traitement des communicateurs predefinis */

  if (pdm_mpi_comm < 0)
    return mpi_comm_cste[-pdm_mpi_comm - 1];

  /* Traitement des communicateurs utilisateurs */

  else {
    if (pdm_mpi_comm < l_mpi_comm)
      return *(mpi_comm[pdm_mpi_comm]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_comm :"
            " pdm_mpi_comm '%d' non valide\n", pdm_mpi_comm);
      abort();
      return MPI_COMM_NULL;
    }
  }
}

/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_comm
 *
 * MPI_Comm -> PDM_MPI_Comm
 *----------------------------------------------------------------------------*/

static PDM_MPI_Comm _mpi_2_pdm_mpi_comm(MPI_Comm comm)
{

  /* Traitement des communicateurs predefinis */

  if (comm == MPI_COMM_NULL)
    return PDM_MPI_COMM_NULL;

  else if (comm == MPI_COMM_WORLD)
    return PDM_MPI_COMM_WORLD;

  /* Traitement des communicateurs utilisateurs */

  else {

    /* Recherche du communicateur MSG correspondant au communicateur MPI */


    if (mpi_comm != NULL) {
      for (int i = 0; i < l_mpi_comm; i++)
        if (mpi_comm[i] != NULL)
          if (*(mpi_comm[i]) == comm)
            return (PDM_MPI_Comm) i;
    }

    /* Si non trouve cree un nouveau communicateur MSG */

    if (mpi_comm == NULL) {
      l_mpi_comm = 4;
      mpi_comm = (MPI_Comm **) malloc(sizeof(MPI_Comm *) * l_mpi_comm);
      for (int i = 0; i < l_mpi_comm; i++)
        mpi_comm[i] = NULL;
    }

    if (l_mpi_comm <= n_mpi_comm) {
      int  p_l_mpi_comm = l_mpi_comm;
      l_mpi_comm = 2 * l_mpi_comm;
      mpi_comm = (MPI_Comm **) realloc((void*) mpi_comm,
                                       l_mpi_comm *
                                       sizeof(MPI_Comm *));
      for (int i = p_l_mpi_comm; i < l_mpi_comm; i++)
        mpi_comm[i] = NULL;
    }

    /* Recherche de la premiere place libre pour stocker le fichier */

    int i = 0;
    while (mpi_comm[i] != NULL)
      i++;

    mpi_comm[i] = (MPI_Comm *) malloc(sizeof(MPI_Comm));
    *(mpi_comm[i]) = comm;
    n_mpi_comm += 1;

    return (PDM_MPI_Comm) i;
  }
}


/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_request
 *
 * PDM_MPI_Request -> MPI_Request
 *----------------------------------------------------------------------------*/

static MPI_Request _pdm_mpi_2_mpi_request(PDM_MPI_Request pdm_mpi_request)
{

  /* Traitement des communicateurs predefinis */

  if (pdm_mpi_request < 0)
    return mpi_request_cste[-pdm_mpi_request - 1];

  /* Traitement des communicateurs utilisateurs */

  else {
    if (pdm_mpi_request < l_mpi_request)
      return *(mpi_request[pdm_mpi_request]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_request :"
            " pdm_mpi_request '%d' non valide\n", pdm_mpi_request);
      abort();
      return MPI_REQUEST_NULL;
    }
  }
}


/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_request
 *
 * PDM_MPI_win -> MPI_win
 *----------------------------------------------------------------------------*/

static MPI_Win _pdm_mpi_2_mpi_win(PDM_MPI_Win pdm_mpi_win)
{

  /* Traitement des communicateurs predefinis */

  if (pdm_mpi_win < 0)
    return mpi_win_cste[-pdm_mpi_win - 1];

  /* Traitement des communicateurs utilisateurs */

  else {
    if (pdm_mpi_win < l_mpi_win)
      return *(mpi_win[pdm_mpi_win]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_win :"
            " pdm_mpi_win '%d' non valide\n", pdm_mpi_win);
      abort();
      return MPI_WIN_NULL;
    }
  }
}


/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_request
 *
 * MPI_Request -> PDM_MPI_Request
 *----------------------------------------------------------------------------*/

// static PDM_MPI_Request _mpi_2_pdm_mpi_request(MPI_Request request)
// {

//   /* Traitement des communicateurs predefinis */

//   if (request == MPI_REQUEST_NULL) {
//     return PDM_MPI_REQUEST_NULL;
//   }

//   /* Traitement des communicateurs utilisateurs */

//   else {

//     /* Recherche du communicateur MSG correspondant au communicateur MPI */

//     if (mpi_request != NULL) {
//       for (int i = 0; i < l_mpi_request; i++)
//         if (mpi_request[i] != NULL)
//           if (*(mpi_request[i]) == request) {
//             return (PDM_MPI_Request) i;
//           }
//     }

//     /* Si non trouve cree un nouveau communicateur MSG */

//     if (mpi_request == NULL) {
//       l_mpi_request = 4;
//       mpi_request = (MPI_Request **) malloc(sizeof(MPI_Request *) * l_mpi_request);
//       for (int i = 0; i < l_mpi_request; i++)
//         mpi_request[i] = NULL;
//     }

//     if (l_mpi_request <= n_mpi_request) {
//       int  p_l_mpi_request = l_mpi_request;
//       l_mpi_request = 2 * l_mpi_request;
//       mpi_request = (MPI_Request **) realloc((void*) mpi_request,
//                                              l_mpi_request *
//                                              sizeof(MPI_Request *));
//       for (int i = p_l_mpi_request; i < l_mpi_request; i++)
//         mpi_request[i] = NULL;
//     }

//     /* Recherche de la premiere place libre pour stocker le fichier */

//     int i = 0;
//     while (mpi_request[i] != NULL)
//       i++;

//     mpi_request[i] = (MPI_Request *) malloc(sizeof(MPI_Request));
//     *(mpi_request[i]) = request;
//     n_mpi_request += 1;

//     return (PDM_MPI_Request) i;
//   }
// }



/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_request
 *
 * MPI_Request -> PDM_MPI_Request
 *----------------------------------------------------------------------------*/

static PDM_MPI_Request _mpi_2_pdm_mpi_request_add(MPI_Request request)
{

  /* Traitement des communicateurs predefinis */

  if (request == MPI_REQUEST_NULL) {
    return PDM_MPI_REQUEST_NULL;
  }

  /* Traitement des communicateurs utilisateurs */

  else {


    /* On stocke le request */

    if (mpi_request == NULL) {
      l_mpi_request = 4;
      mpi_request = (MPI_Request **) malloc(sizeof(MPI_Request *) * l_mpi_request);
      for (int i = 0; i < l_mpi_request; i++)
        mpi_request[i] = NULL;
    }

    if (l_mpi_request <= n_mpi_request) {
      int  p_l_mpi_request = l_mpi_request;
      l_mpi_request = 2 * l_mpi_request;
      mpi_request = (MPI_Request **) realloc((void*) mpi_request,
                                             l_mpi_request *
                                             sizeof(MPI_Request *));
      for (int i = p_l_mpi_request; i < l_mpi_request; i++)
        mpi_request[i] = NULL;
    }

    /* Recherche de la premiere place libre pour stocker le fichier */

    int i = 0;
    while (mpi_request[i] != NULL)
      i++;

    mpi_request[i] = (MPI_Request *) malloc(sizeof(MPI_Request));
    *(mpi_request[i]) = request;
    n_mpi_request += 1;

    return (PDM_MPI_Request) i;
  }
}


/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_win
 *
 * MPI_Win -> PDM_MPI_win
 *----------------------------------------------------------------------------*/

static PDM_MPI_Win _mpi_2_pdm_mpi_win_add(MPI_Win win)
{

  /* Traitement des communicateurs predefinis */

  if (win == MPI_WIN_NULL) {
    return PDM_MPI_WIN_NULL;
  }

  /* Traitement des communicateurs utilisateurs */

  else {


    /* On stocke le win */

    if (mpi_win == NULL) {
      l_mpi_win = 4;
      mpi_win = (MPI_Win **) malloc(sizeof(MPI_Win *) * l_mpi_win);
      for (int i = 0; i < l_mpi_win; i++)
        mpi_win[i] = NULL;
    }

    if (l_mpi_win <= n_mpi_win) {
      int  p_l_mpi_win = l_mpi_win;
      l_mpi_win = 2 * l_mpi_win;
      mpi_win = (MPI_Win **) realloc((void*) mpi_win,
                                             l_mpi_win *
                                             sizeof(MPI_Win *));
      for (int i = p_l_mpi_win; i < l_mpi_win; i++)
        mpi_win[i] = NULL;
    }

    /* Recherche de la premiere place libre pour stocker le fichier */

    int i = 0;
    while (mpi_win[i] != NULL)
      i++;

    mpi_win[i] = (MPI_Win *) malloc(sizeof(MPI_Win));
    *(mpi_win[i]) = win;
    n_mpi_win += 1;

    return (PDM_MPI_Win) i;
  }
}

/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_datatype
 *
 * PDM_MPI_Datatype -> MPI_Datatype
 *----------------------------------------------------------------------------*/

static MPI_Datatype _pdm_mpi_2_mpi_datatype(PDM_MPI_Datatype pdm_mpi_datatype)
{

  /* Traitement des MPI_Datatype connus  */

  if (pdm_mpi_datatype < 0)
    return mpi_datatype_cste[-pdm_mpi_datatype - 1];

  /* Traitement des MPI_Datatype utilisateurs  */

  else {
    if (pdm_mpi_datatype < l_mpi_datatype)
      return *(mpi_datatype[pdm_mpi_datatype]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_datatype :"
            " pdm_mpi_datatype '%d' non valide\n", pdm_mpi_datatype);
      abort();
      return MPI_DATATYPE_NULL;
    }
  }
}


/*----------------------------------------------------------------------------
 * _mpi_2_pdm_mpi_data
 *
 * MPI_Datatype -> PDM_MPI_Datatype
 *----------------------------------------------------------------------------*/

static PDM_MPI_Datatype _mpi_2_pdm_mpi_datatype(MPI_Datatype datatype)
{

  /* Traitement des communicateurs predefinis */

  if (datatype == MPI_BYTE)
    return PDM_MPI_BYTE;
  else if (datatype == MPI_PACKED)
    return  PDM_MPI_PACKED;
  else if (datatype == MPI_CHAR)
    return  PDM_MPI_CHAR;
  else if (datatype == MPI_SHORT)
    return  PDM_MPI_SHORT;
  else if (datatype == MPI_INT)
    return  PDM_MPI_INT;
  else if (datatype == MPI_LONG)
    return  PDM_MPI_LONG;
  else if (datatype == MPI_FLOAT)
    return  PDM_MPI_FLOAT;
  else if (datatype == MPI_DOUBLE)
    return  PDM_MPI_DOUBLE;
  else if (datatype == MPI_LONG_DOUBLE)
    return  PDM_MPI_LONG_DOUBLE;
  else if (datatype == MPI_UNSIGNED_CHAR)
    return  PDM_MPI_UNSIGNED_CHAR;
  else if (datatype == MPI_UNSIGNED_SHORT)
    return  PDM_MPI_UNSIGNED_SHORT;
  else if (datatype == MPI_UNSIGNED_LONG)
    return  PDM_MPI_UNSIGNED_LONG;
  else if (datatype == MPI_UNSIGNED_LONG_LONG)
    return  PDM_MPI_UNSIGNED_LONG_LONG;
  else if (datatype == MPI_UNSIGNED)
    return PDM_MPI_UNSIGNED;
  else if (datatype == MPI_FLOAT_INT)
    return  PDM_MPI_FLOAT_INT;
  else if (datatype == MPI_DOUBLE_INT)
    return  PDM_MPI_DOUBLE_INT;
  else if (datatype == MPI_LONG_DOUBLE_INT)
    return PDM_MPI_LONG_DOUBLE_INT;
  else if (datatype == MPI_LONG_INT)
    return PDM_MPI_LONG_INT;
  else if (datatype == MPI_SHORT_INT)
    return PDM_MPI_SHORT_INT;
  else if (datatype == MPI_2INT)
    return PDM_MPI_2INT;
  else if (datatype == MPI_CHARACTER)
    return PDM_MPI_CHARACTER;
  else if (datatype == MPI_INTEGER)
    return PDM_MPI_INTEGER;
  else if (datatype == MPI_REAL)
    return PDM_MPI_REAL;
  else if (datatype == MPI_DOUBLE_PRECISION)
    return PDM_MPI_DOUBLE_PRECISION;
  else if (datatype == MPI_DATATYPE_NULL)
    return PDM_MPI_DATATYPE_NULL;
  else if (datatype == MPI_INT8_T)
    return PDM_MPI_INT8_T;
  else if (datatype == MPI_INT16_T)
    return PDM_MPI_INT16_T;
  else if (datatype == MPI_INT32_T)
    return PDM_MPI_INT32_T;
  else if (datatype == MPI_INT64_T)
    return PDM_MPI_INT64_T;
  else if (datatype == MPI_UINT8_T)
    return PDM_MPI_UINT8_T;
  else if (datatype == MPI_UINT16_T)
    return PDM_MPI_UINT16_T;
  else if (datatype == MPI_UINT32_T)
    return PDM_MPI_UINT32_T;
  else if (datatype == MPI_UINT64_T)
    return PDM_MPI_UINT64_T;

  /* Traitement des communicateurs utilisateurs */

  else {

    /* Recherche du datatype MSG correspondant au datatype MPI */

    if (mpi_datatype != NULL) {
      for (int i = 0; i < l_mpi_datatype; i++)
        if (mpi_datatype[i] != NULL)
          if (*(mpi_datatype[i]) == datatype)
            return (PDM_MPI_Datatype) i;
    }

    /* Si non trouve cree un nouveau datatype MSG */

    if (mpi_datatype == NULL) {
      l_mpi_datatype = 4;
      mpi_datatype = (MPI_Datatype **)
        malloc(sizeof(MPI_Datatype *) * l_mpi_datatype);
      for (int i = 0; i < l_mpi_datatype; i++)
        mpi_datatype[i] = NULL;
    }

    if (l_mpi_datatype <= n_mpi_datatype) {
      int  p_l_mpi_datatype = l_mpi_datatype;
      l_mpi_datatype = 2 * l_mpi_datatype;
      mpi_datatype = (MPI_Datatype **) realloc((void*) mpi_datatype,
                                       l_mpi_datatype *
                                       sizeof(MPI_Datatype *));
      for (int i = p_l_mpi_datatype; i < l_mpi_datatype; i++)
        mpi_datatype[i] = NULL;
    }

    /* Recherche de la premiere place libre pour stocker le fichier */

    int i = 0;
    while (mpi_datatype[i] != NULL)
      i++;

    mpi_datatype[i] = (MPI_Datatype *) malloc(sizeof(MPI_Datatype));
    *(mpi_datatype[i]) = datatype;
    n_mpi_datatype += 1;

    return (PDM_MPI_Datatype) i;
  }
}



/*----------------------------------------------------------------------------
 * PDM_MPI_File_Create
 *
 * PDM_MPI_File -> MPI_File
 *----------------------------------------------------------------------------*/

static PDM_MPI_File _pdm_mpi_file_create(void)
{

  /* Si non trouve, on cree un nouveau fichier MSG */

  if (mpi_file == NULL) {
    l_mpi_file = 4;
      mpi_file = (MPI_File **) malloc(sizeof(MPI_File *) * l_mpi_file);
      for (int i = 0; i < l_mpi_file; i++)
        mpi_file[i] = NULL;
  }

  if (l_mpi_file <= n_mpi_file) {
    int  p_l_mpi_file = l_mpi_file;
    l_mpi_file = 2 * l_mpi_file;
    mpi_file = (MPI_File **) realloc((void*) mpi_file,
                                     l_mpi_file *
                                     sizeof(MPI_File *));
    for (int i = p_l_mpi_file; i < l_mpi_file; i++)
      mpi_file[i] = NULL;
  }

  /* Recherche de la premiere place libre pour stocker le fichier */

  int i = 0;
  while (mpi_file[i] != NULL)
    i++;

  mpi_file[i] = (MPI_File *) malloc(sizeof(MPI_File));
  n_mpi_file += 1;
  return i;
}

/*----------------------------------------------------------------------------
 * _pdm_mpi_2_mpi_file
 *
 * PDM_MPI_File -> MPI_File
 *----------------------------------------------------------------------------*/

static MPI_File _pdm_mpi_2_mpi_file(PDM_MPI_File pdm_mpi_file)
{

  /* Traitement des MPI_File connus  */

  if (pdm_mpi_file < 0)
    return mpi_file_cste[-pdm_mpi_file - 1];

  /* Traitement des MPI_File utilisateurs  */

  else {
    if (pdm_mpi_file < l_mpi_file)
      return *(mpi_file[pdm_mpi_file]);
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_file :"
              " pdm_mpi_file '%d' non valide\n", pdm_mpi_file);
      abort();
      return MPI_FILE_NULL;
    }
  }
}

/*============================================================================
 * Defintion des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Init(int *argc, char ***argv)
{
  return MPI_Init(argc, argv);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Init
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

int PDM_MPI_Finalize (void)
{
  if (mpi_file != NULL) {
    for (int i = 0; i < l_mpi_file; i++) {
      if (mpi_file[i] != NULL) {
        MPI_File_close(mpi_file[i]);
        mpi_file[i] = NULL;
      }
    }
    free(mpi_file);
    l_mpi_file = 0;
    n_mpi_file = 0;
  }
  if (mpi_comm != NULL) {
    for (int i = 0; i < l_mpi_comm; i++) {
      if (mpi_comm[i] != NULL) {
        MPI_Comm_free(mpi_comm[i]);
        free (mpi_comm[i]);
        mpi_comm[i] = NULL;
      }
    }
    free(mpi_comm);
    l_mpi_comm = 0;
    n_mpi_comm = 0;
  }

  if (mpi_request != NULL) {
    for (int i = 0; i < l_mpi_request; i++) {
      if (mpi_request[i] != NULL) {
        MPI_Request_free(mpi_request[i]);
        mpi_request[i] = NULL;
      }
    }
    free(mpi_request);
    l_mpi_request = 0;
    n_mpi_request = 0;
  }

  if (mpi_win != NULL) {
    for (int i = 0; i < l_mpi_win; i++) {
      if (mpi_win[i] != NULL) {
        MPI_Win_free(mpi_win[i]);
        mpi_win[i] = NULL;
      }
    }
    free(mpi_win);
    l_mpi_win = 0;
    n_mpi_win = 0;
  }

  if (mpi_datatype != NULL) {
    for (int i = 0; i < l_mpi_datatype; i++) {
      if (mpi_datatype[i] != NULL) {
        MPI_Type_free(mpi_datatype[i]);
        mpi_datatype[i] = NULL;
      }
    }
    free(mpi_datatype);
    l_mpi_datatype = 0;
    n_mpi_datatype = 0;
  }
  return MPI_Finalize();
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_2_mpi_comm(PDM_MPI_Comm pdm_mpi_comm)
{

  /* Traitement des communicateurs predefinis */

  if (pdm_mpi_comm < 0)
    return (void *) &mpi_comm_cste[-pdm_mpi_comm - 1];

  /* Traitement des communicateurs utilisateurs */

  else {
    if (pdm_mpi_comm < l_mpi_comm)
      return (void *) mpi_comm[pdm_mpi_comm];
    else {
      PDM_error(__FILE__, __LINE__, 0,"_pdm_mpi_2_mpi_comm :"
            " pdm_mpi_comm '%d' non valide\n", pdm_mpi_comm);
      abort();
      return NULL;
    }
  }
}

/*----------------------------------------------------------------------------
 * pdm_mpi_2_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

void *PDM_MPI_free_mpi_comm(void *pt_mpi_comm)
{

  MPI_Comm *comm = (MPI_Comm *) pt_mpi_comm;
  MPI_Comm_free (comm);
  return NULL;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_mpi_2_pdm_mpi_comm
 *
 * PDM_MPI_Comm -> MPI_Comm
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_mpi_2_pdm_mpi_comm(void *pt_mpi_comm)
{

  MPI_Comm _mpi_comm = *((MPI_Comm *) pt_mpi_comm);
  return _mpi_2_pdm_mpi_comm(_mpi_comm);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_File_open (wrapping de la fonction MPI_File_open)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_open(PDM_MPI_Comm comm, char *filename, int amode, PDM_MPI_File *fh)
{

  *fh = _pdm_mpi_file_create();

  char *hints = getenv("PDM_IO_HINTS");

  MPI_Info hints_mpi = MPI_INFO_NULL;

  if (hints != NULL) {

    MPI_Info_create (&hints_mpi);

    char *cp_hints = malloc (sizeof(char *) * (strlen(hints) + 1));
    char *name = malloc (sizeof(char *) * (strlen(hints) + 1));
    char *value = malloc (sizeof(char *) * (strlen(hints) + 1));
    strcpy (cp_hints, hints);

    char *pch;
    char *str2 = cp_hints;

    do {
      pch = strtok (str2,"=");
      str2 = NULL;
      if (pch != NULL) {
        strcpy(name, pch);
        pch = strtok (str2, ":");
        if (pch == NULL) {
          PDM_printf ("Error PDM_MPI_File_open : No value for hint \"%s\"."
                  " Check \"PDM_IO_HINTS\" environment variable\n", name);
          exit(1);
        }
        else {
          strcpy(value, pch);
          MPI_Info_set (hints_mpi, name, value);
          PDM_printf ("MPI/IO hint \"%s\" = \"%s\"\n", name, value);
        }
      }
    } while (pch != NULL);

    free (cp_hints);
    free (name);
    free (value);

  }

  int code = MPI_File_open(_pdm_mpi_2_mpi_comm(comm),
                           filename,
                           mpi_file_mode[amode],
                           hints_mpi,
                           mpi_file[*fh]);

  if (hints != NULL) {

    MPI_Info_free(&hints_mpi);

  }

  if (code != MPI_SUCCESS) {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_close (wrapping de la fonction MPI_File_close)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_close(PDM_MPI_File *fh)
{
  int code =  MPI_File_close(mpi_file[*fh]);

  free(mpi_file[*fh]);

  mpi_file[*fh] = NULL;
  n_mpi_file -= 1;

  if (code != MPI_SUCCESS) {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_seek (wrapping de la fonction MPI_File_seek)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_seek(PDM_MPI_File fh, PDM_MPI_Offset offset, int whence)
{
  int code = MPI_File_seek(_pdm_mpi_2_mpi_file(fh),
                           (MPI_Offset) offset,
                           whence);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_size (wrapping de la fonction MPI_File_get_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_size(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  MPI_Offset _tmp_offset;
  int code = MPI_File_get_size(_pdm_mpi_2_mpi_file(fh),
                           (MPI_Offset*) &_tmp_offset);
  *offset = _tmp_offset;
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_position (wrapping de la fonction MPI_File_get_position)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_position(PDM_MPI_File fh, PDM_MPI_Offset *offset)
{
  MPI_Offset _tmp_offset;
  int code = MPI_File_get_position(_pdm_mpi_2_mpi_file(fh),
                                    &_tmp_offset);
  *offset = _tmp_offset;
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_set_view (wrapping de la fonction MPI_File_set_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_set_view(PDM_MPI_File fh, PDM_MPI_Offset disp, PDM_MPI_Datatype etype,
	              PDM_MPI_Datatype filetype, const char *datarep)
{
  int code = MPI_File_set_view(_pdm_mpi_2_mpi_file(fh),
                               (MPI_Offset) disp,
                               _pdm_mpi_2_mpi_datatype(etype),
                               _pdm_mpi_2_mpi_datatype(filetype),
                               datarep,
                               MPI_INFO_NULL);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_get_view (wrapping de la fonction MPI_File_get_view)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_get_view(PDM_MPI_File fh, PDM_MPI_Offset *disp,
                      PDM_MPI_Datatype *etype, PDM_MPI_Datatype *filetype, char *datarep)
{

  MPI_Datatype mpi_etype;
  MPI_Datatype mpi_filetype;
  MPI_Offset _disp = (MPI_Offset) *disp;

  int code = MPI_File_get_view(_pdm_mpi_2_mpi_file(fh),
                               &_disp,
                               &mpi_etype,
                               &mpi_filetype,
                               datarep);

  *etype    = _mpi_2_pdm_mpi_datatype(mpi_etype);
  *filetype = _mpi_2_pdm_mpi_datatype(mpi_filetype);
  *disp     = (PDM_MPI_Offset) _disp;

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at (wrapping de la fonction MPI_File_read_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                     int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code = MPI_File_read_at(_pdm_mpi_2_mpi_file(fh),
                              (MPI_Offset) offset,
                              buf,
                              count,
                              _pdm_mpi_2_mpi_datatype(datatype),
                              &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_at_all (wrapping de la fonction MPI_File_read_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code = MPI_File_read_at_all(_pdm_mpi_2_mpi_file(fh),
                                  (MPI_Offset) offset,
                                  buf,
                                  count,
                                  _pdm_mpi_2_mpi_datatype(datatype),
                                  &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at (wrapping de la fonction MPI_File_write_at)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                      int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  MPI_Offset _offset = (MPI_Offset) offset;
  int code = MPI_File_write_at(_pdm_mpi_2_mpi_file(fh),
                               _offset,
                               buf,
                               count,
                               _pdm_mpi_2_mpi_datatype(datatype),
                               &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_at_all (wrapping de la fonction MPI_File_write_at_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_at_all(PDM_MPI_File fh, PDM_MPI_Offset offset, void *buf,
                          int count, PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  MPI_Offset _offset = (MPI_Offset) offset;
  int code = MPI_File_write_at_all(_pdm_mpi_2_mpi_file(fh),
                                   (MPI_Offset) _offset,
                                   buf,
                                   count,
                                   _pdm_mpi_2_mpi_datatype(datatype),
                                   &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read (wrapping de la fonction MPI_File_read)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read(PDM_MPI_File fh, void *buf, int count,
                  PDM_MPI_Datatype datatype, int *n_octet_lus)
{

  MPI_Status status;

  int code =  MPI_File_read(_pdm_mpi_2_mpi_file(fh),
                            buf,
                            count,
                            _pdm_mpi_2_mpi_datatype(datatype),
                            &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_read_all (wrapping de la fonction MPI_File_read_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_read_all(PDM_MPI_File fh, void *buf, int count,
                      PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  int code = MPI_File_read_all(_pdm_mpi_2_mpi_file(fh),
                                buf,
                                count,
                                _pdm_mpi_2_mpi_datatype(datatype),
                                &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write (wrapping de la fonction MPI_File_write)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write(PDM_MPI_File fh, void *buf, int count,
                   PDM_MPI_Datatype datatype, int *n_octet_lus)
{
  MPI_Status status;

  int code =  MPI_File_write(_pdm_mpi_2_mpi_file(fh),
                             buf,
                             count,
                             _pdm_mpi_2_mpi_datatype(datatype),
                             &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_File_write_all (wrapping de la fonction MPI_File_write_all)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_File_write_all(PDM_MPI_File fh, void *buf, int count,
                       PDM_MPI_Datatype datatype, int *n_octet_lus)

{

  MPI_Status status;

  int code =  MPI_File_write_all(_pdm_mpi_2_mpi_file(fh),
                                 buf,
                                 count,
                                 _pdm_mpi_2_mpi_datatype(datatype),
                                 &status);

  if (code == MPI_SUCCESS)
    MPI_Get_count(&status, MPI_BYTE, n_octet_lus);
  else {
    char buffer[MPI_MAX_ERROR_STRING];
    int  buffer_len;

    MPI_Error_string(code, buffer, &buffer_len);

    PDM_error(__FILE__, __LINE__, 0, "%s\n", buffer);

    abort();
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gather (wrapping de la fonction MPI_Gather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm)
{
  int code = MPI_Gather(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                        recvbuf, recvcount, _pdm_mpi_2_mpi_datatype(recvtype),
                        root, _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Igather (wrapping de la fonction MPI_Igather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Igather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
               void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
               int root, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code = MPI_Igather(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                        recvbuf, recvcount, _pdm_mpi_2_mpi_datatype(recvtype),
                        root, _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Gatherv (wrapping de la fonction MPI_Gatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Gatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int *recvcounts, int *displs,
                PDM_MPI_Datatype recvtype, int root, PDM_MPI_Comm comm)
{
  int code = MPI_Gatherv(sendbuf,
                         sendcount,
                         _pdm_mpi_2_mpi_datatype(sendtype),
                         recvbuf,
                         recvcounts,
                         displs,
                         _pdm_mpi_2_mpi_datatype(recvtype),
                         root,
                         _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Recv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
             int tag, PDM_MPI_Comm comm)
{
  int code =  MPI_Recv(buf, count, _pdm_mpi_2_mpi_datatype(datatype), source,
                       tag, _pdm_mpi_2_mpi_comm(comm), MPI_STATUS_IGNORE);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Recv (wrapping de la fonction MPI_Recv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Irecv(void *buf, int count, PDM_MPI_Datatype datatype, int source,
              int tag, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code =  MPI_Irecv(buf, count, _pdm_mpi_2_mpi_datatype(datatype), source,
                       tag, _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  assert(code == 0);
  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  assert(_mpi_request != MPI_REQUEST_NULL);

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Send (wrapping de la fonction MPI_Send)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Send(void *buf, int count, PDM_MPI_Datatype datatype, int dest,
             int tag, PDM_MPI_Comm comm)
{
  int code = MPI_Send(buf, count, _pdm_mpi_2_mpi_datatype(datatype), dest,
                      tag, _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Issend (wrapping de la fonction MPI_Issend)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Issend(const void *buf, int count, PDM_MPI_Datatype datatype, int dest, int tag,
               PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  int code = MPI_Issend(buf, count, _pdm_mpi_2_mpi_datatype(datatype), dest,
                        tag, _pdm_mpi_2_mpi_comm(comm), &_mpi_request);

  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  assert(code == 0);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Wait (wrapping de la fonction MPI_Wait)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Wait(PDM_MPI_Request *request)

{
  MPI_Request _request = _pdm_mpi_2_mpi_request(*request);
  int code = MPI_Wait(&_request, MPI_STATUS_IGNORE);
  assert(code == 0);

  free(mpi_request[*request]);
  mpi_request[*request] = NULL;
  n_mpi_request += -1;
  *request = PDM_MPI_REQUEST_NULL;

  if (n_mpi_request == 0) {
    free(mpi_request);
    mpi_request = NULL;

    l_mpi_request = 0;
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Test (wrapping de la fonction MPI_Test)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Test(PDM_MPI_Request *request, int *flag)
{
  MPI_Request _request = _pdm_mpi_2_mpi_request(*request);

  // Test was already done
  if(_request == MPI_REQUEST_NULL) {
    *flag = 1;
    return _mpi_2_pdm_mpi_err(MPI_SUCCESS);
  }
  int code = MPI_Test(&_request, flag, MPI_STATUS_IGNORE);

  if(*flag == 0) {
    return code; // Message was not ready
  }

  free(mpi_request[*request]);
  mpi_request[*request] = NULL;
  n_mpi_request += -1;
  *request = PDM_MPI_REQUEST_NULL;

  if (n_mpi_request == 0) {
    free(mpi_request);
    mpi_request = NULL;

    l_mpi_request = 0;
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_hindexed (wrapping de la fonction MPI_Type_hindexed)
 *
 *----------------------------------------------------------------------------*/


int PDM_MPI_Type_create_hindexed (int count,
                              const int array_of_blocklengths[],
                              const PDM_MPI_Aint array_of_displacements[],
                              PDM_MPI_Datatype oldtype,
                              PDM_MPI_Datatype *newtype)
{
  MPI_Datatype mpi_newtype;
  MPI_Aint *_array_of_displacements = malloc (sizeof(MPI_Aint) * count);

  for (int i = 0; i < count; i++) {
    _array_of_displacements[i] = array_of_displacements[i];
  }

  int code = MPI_Type_create_hindexed(count,
                               array_of_blocklengths,
                               _array_of_displacements,
                               _pdm_mpi_2_mpi_datatype(oldtype),
                               &mpi_newtype);

  *newtype = _mpi_2_pdm_mpi_datatype(mpi_newtype);
  free (_array_of_displacements);
  return _mpi_2_pdm_mpi_err(code);

}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_hindexed (wrapping de la fonction MPI_Type_hindexed)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_create_contiguous(int               count,
                                   PDM_MPI_Datatype  old_datatype,
                                   PDM_MPI_Datatype *newtype)
{
  MPI_Datatype mpi_newtype;
  int code = MPI_Type_contiguous(count,
                                 _pdm_mpi_2_mpi_datatype(old_datatype),
                                 &mpi_newtype);
  assert(code == 0);

  *newtype = _mpi_2_pdm_mpi_datatype(mpi_newtype);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_commit (wrapping de la fonction MPI_Type_commit)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_commit(PDM_MPI_Datatype *datatype)
{
  MPI_Datatype mpi_type = _pdm_mpi_2_mpi_datatype(*datatype);
  int code =  MPI_Type_commit(&mpi_type);
  *datatype = _mpi_2_pdm_mpi_datatype(mpi_type);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * MPI_Type_size (wrapping de la fonction MPI_Type_commit)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Type_size(PDM_MPI_Datatype datatype, int *size)
{
  return MPI_Type_size(_pdm_mpi_2_mpi_datatype(datatype), size);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Type_free (wrapping de la fonction MPI_Type_free)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Type_free(PDM_MPI_Datatype *datatype)
{
  MPI_Datatype mpi_type = _pdm_mpi_2_mpi_datatype(*datatype);
  int code = MPI_Type_free(&mpi_type);
  free(mpi_datatype[*datatype]);
  mpi_datatype[*datatype] = NULL;
  *datatype = PDM_MPI_DATATYPE_NULL;
  n_mpi_datatype += -1;
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_f2c (wrapping de la fonction MPI_comm_f2c)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Comm PDM_MPI_Comm_f2c(PDM_MPI_Fint comm)
{

  /* Conversion Fortran vers C */

  MPI_Comm _mpi_comm = MPI_Comm_f2c(comm);
  PDM_MPI_Comm c_comm = _mpi_2_pdm_mpi_comm(_mpi_comm);
  return c_comm;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_c2f (wrapping de la fonction MPI_comm_c2f)
 *
 *----------------------------------------------------------------------------*/

PDM_MPI_Fint PDM_MPI_Comm_c2f(PDM_MPI_Comm comm)
{

  /* Conversion Fortran vers C */

  MPI_Comm _mpi_comm = _pdm_mpi_2_mpi_comm(comm);
  PDM_MPI_Fint f_comm = (PDM_MPI_Fint) MPI_Comm_c2f(_mpi_comm);
  return f_comm;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Scatter (wrapping de la fonction MPI_Scatter)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scatter(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                void *recvbuf, int recvcount, PDM_MPI_Datatype recvtype,
                int root, PDM_MPI_Comm comm)
{
  int code = MPI_Scatter(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                         recvbuf, recvcount, _pdm_mpi_2_mpi_datatype(recvtype),
                         root, _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Barrier (wrapping de la fonction MPI_Barrier)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Barrier(PDM_MPI_Comm comm)
{
  int code =  MPI_Barrier(_pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Wtime (wrapping de la fonction MPI_Wtime)
 *
 *----------------------------------------------------------------------------*/

double PDM_MPI_Wtime(void)
{

  return MPI_Wtime();
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Bcast (wrapping de la fonction MPI_Bcast)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Bcast(void *buffer, int count, PDM_MPI_Datatype datatype,
                  int root, PDM_MPI_Comm comm)
{
  int code = MPI_Bcast(buffer,
                       count,
                       _pdm_mpi_2_mpi_datatype(datatype),
                       root,
                       _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_IBcast (wrapping de la fonction MPI_IBcast)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ibcast(void *buffer, int count, PDM_MPI_Datatype datatype,
                   int root, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;

  int code = MPI_Ibcast(buffer,
                        count,
                        _pdm_mpi_2_mpi_datatype(datatype),
                        root,
                        _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                      void *recvbuf, int recvcount,
                      PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code =  MPI_Allgather(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                            recvbuf, recvcount,
                            _pdm_mpi_2_mpi_datatype(recvtype),
                            _pdm_mpi_2_mpi_comm(comm));
  assert(code == MPI_SUCCESS);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allgatherv (wrapping de la fonction MPI_Allgatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allgatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                       void *recvbuf, int *recvcounts,
                       int *displs, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Allgatherv(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                            recvbuf, recvcounts, displs,
                            _pdm_mpi_2_mpi_datatype(recvtype), _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce (wrapping de la fonction MPI_Reduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Reduce(void *sendbuf, void *recvbuf, int count,
		   PDM_MPI_Datatype datatype, PDM_MPI_Op op,
		   int root, PDM_MPI_Comm comm)
{
  int code = MPI_Reduce(sendbuf, recvbuf, count,
                           _pdm_mpi_2_mpi_datatype(datatype),
                           mpi_op[op], root,
                           _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Reduce_scatter (wrapping de la fonction MPI_Reduce_scatter)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *counts,
                           PDM_MPI_Datatype datatype, PDM_MPI_Op op,
                           PDM_MPI_Comm comm)
{
  int code = MPI_Reduce_scatter(sendbuf, recvbuf, counts,
                                _pdm_mpi_2_mpi_datatype(datatype),
                                mpi_op[op],
                                _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Allreduce (wrapping de la fonction MPI_Allreduce)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                  PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm)
{
  int code = MPI_Allreduce(sendbuf, recvbuf, count,
                           _pdm_mpi_2_mpi_datatype(datatype),
                           mpi_op[op],
                           _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Scan (wrapping de la fonction MPI_Scan)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Scan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm)
{
  int code = MPI_Scan(sendbuf, recvbuf, count,
                      _pdm_mpi_2_mpi_datatype(datatype),
                      mpi_op[op],
                      _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Exscan (wrapping de la fonction MPI_Exscan)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Exscan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm)
{
  int code = MPI_Exscan(sendbuf, recvbuf, count,
                        _pdm_mpi_2_mpi_datatype(datatype),
                        mpi_op[op],
                        _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}



int PDM_MPI_Iscan(const void *sendbuf, void *recvbuf, int count,
             PDM_MPI_Datatype datatype, PDM_MPI_Op op, PDM_MPI_Comm comm,
             PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;

  int code = MPI_Iscan(sendbuf, recvbuf, count,
                      _pdm_mpi_2_mpi_datatype(datatype),
                      mpi_op[op],
                      _pdm_mpi_2_mpi_comm(comm), &_mpi_request);

  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoall (wrapping de la fonction MPI_Alltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                 void *recvbuf, int recvcount,
                 PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Alltoall(sendbuf, sendcount,
                          _pdm_mpi_2_mpi_datatype(sendtype),
                          recvbuf, recvcount,
                          _pdm_mpi_2_mpi_datatype(recvtype),
                          _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Ialltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                 void *recvbuf, int recvcount,
                 PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;

  int code = MPI_Ialltoall(sendbuf, sendcount,
                          _pdm_mpi_2_mpi_datatype(sendtype),
                          recvbuf, recvcount,
                          _pdm_mpi_2_mpi_datatype(recvtype),
                          _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Alltoallv (wrapping de la fonction MPI_Alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                      PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                      int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Alltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           _pdm_mpi_2_mpi_datatype(sendtype),
                           recvbuf,
                           recvcounts,
                           rdispls,
                           _pdm_mpi_2_mpi_datatype(recvtype),
                           _pdm_mpi_2_mpi_comm(comm));

  return _mpi_2_pdm_mpi_err(code);
}


int PDM_MPI_Alltoallv_l(void *sendbuf, int *sendcounts, size_t *sdispls,
                      PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                      size_t *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_SUCCESS;

  int size;
  MPI_Comm_size(_pdm_mpi_2_mpi_comm(comm), &size);

  INT_MAX;
  int coeff = 4;
  int large = 0;
  //  for (int i = 0; i < size; i++) {
  if ((sdispls[size-1] > (size_t) (INT_MAX/coeff)) || (rdispls[size-1] > (size_t) (INT_MAX/coeff))) {
    large = 1;
  }
  //}

  int s_large = 0;
  MPI_Allreduce (&large, &s_large, 1, MPI_INT,  MPI_SUM,  _pdm_mpi_2_mpi_comm(comm));
  large = s_large;


  if (!large) {

    int *_sdispls = malloc(sizeof(int) * size);
    int *_rdispls = malloc(sizeof(int) * size);

    for (int i = 0; i < size; i++) {
      _sdispls[i] = (int) sdispls[i];
      _rdispls[i] = (int) rdispls[i];
    }

    MPI_Alltoallv(sendbuf,
                  sendcounts,
                  _sdispls,
                  _pdm_mpi_2_mpi_datatype(sendtype),
                  recvbuf,
                  recvcounts,
                  _rdispls,
                  _pdm_mpi_2_mpi_datatype(recvtype),
                  _pdm_mpi_2_mpi_comm(comm));

    free (_sdispls);
    free (_rdispls);
  }

  else {

    MPI_Request *request_r = malloc(sizeof(MPI_Request) * size);
    MPI_Request *request_s = malloc(sizeof(MPI_Request) * size);

    int size_sendType;
    MPI_Type_size(_pdm_mpi_2_mpi_datatype(sendtype), &size_sendType);

    int size_recvType;
    MPI_Type_size(_pdm_mpi_2_mpi_datatype(recvtype), &size_recvType);

    for (int i = 0; i < size; i++) {
      if (recvcounts[i] != 0) {
        void *buf = (void *) ((unsigned char*) recvbuf + rdispls[i] * size_recvType);
        code = MPI_Irecv(buf, recvcounts[i], _pdm_mpi_2_mpi_datatype(recvtype), i,
                         0, _pdm_mpi_2_mpi_comm(comm), request_r + i);
        if (code != MPI_SUCCESS) {
          break;
        }
      }
      if (sendcounts[i] != 0) {
        void *buf = (void *) ((unsigned char*) sendbuf + sdispls[i] * size_sendType);
        code = MPI_Issend(buf, sendcounts[i], _pdm_mpi_2_mpi_datatype(sendtype), i,
                          0, _pdm_mpi_2_mpi_comm(comm), request_s + i);
        if (code != MPI_SUCCESS) {
          break;
        }
      }
    }

    if (code != MPI_SUCCESS) {
      return _mpi_2_pdm_mpi_err(code);
    }

    for (int i = 0; i < size; i++) {
      if (recvcounts[i] != 0) {
        code = MPI_Wait(request_r + i, MPI_STATUS_IGNORE);
      }
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    if (code != MPI_SUCCESS) {
      return _mpi_2_pdm_mpi_err(code);
    }

    for (int i = 0; i < size; i++) {
      if (sendcounts[i] != 0) {
        code = MPI_Wait(request_s + i, MPI_STATUS_IGNORE);
      }
      if (code != MPI_SUCCESS) {
        break;
      }
    }

    free (request_r);
    free (request_s);
  }

  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoallv (wrapping de la fonction MPI_Ialltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ialltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                       PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                       int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm,
                       PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  // double t1 = MPI_Wtime();
  int code = MPI_Ialltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           _pdm_mpi_2_mpi_datatype(sendtype),
                           recvbuf,
                           recvcounts,
                           rdispls,
                           _pdm_mpi_2_mpi_datatype(recvtype),
                           _pdm_mpi_2_mpi_comm(comm), &_mpi_request);

  // double dt = MPI_Wtime() - t1;

  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);

  return _mpi_2_pdm_mpi_err(code);
}

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
                           PDM_MPI_Comm      comm)
{
  int code = 0;

  /*
   * Exchange in target view the correct displacement for MPI_Get
   */
  int n_rank, i_rank;
  MPI_Comm_size(_pdm_mpi_2_mpi_comm(comm), &n_rank);
  MPI_Comm_rank(_pdm_mpi_2_mpi_comm(comm), &i_rank);

  int *target_disp = (int *) malloc(n_rank * sizeof(int));

  MPI_Alltoall(sdispls    , 1, MPI_INT,
               target_disp, 1, MPI_INT, _pdm_mpi_2_mpi_comm(comm));


  // double t1 = MPI_Wtime();
  for(int i = 0; i < n_rank; ++i) {

    int   origin_data_size = -1;
    MPI_Type_size(_pdm_mpi_2_mpi_datatype(recvtype), &origin_data_size);

    int            origin_displ = rdispls[i]; // + recvcounts[i] *
    unsigned char *origin_addr  = (unsigned char *) recvbuf + origin_displ * origin_data_size;
    int            origin_count = recvcounts[i];

    if(origin_count > 0 && i != i_rank) {
      MPI_Get(origin_addr,
              origin_count,
              _pdm_mpi_2_mpi_datatype(recvtype),
              i,
              (MPI_Aint) target_disp[i],
              origin_count,
              _pdm_mpi_2_mpi_datatype(sendtype),
              _pdm_mpi_2_mpi_win(win_send));
    }

  }

  // double dt = MPI_Wtime() - t1;

  // t1 = MPI_Wtime();
  int   origin_data_size = -1;
  MPI_Type_size(_pdm_mpi_2_mpi_datatype(recvtype), &origin_data_size);

  int            origin_displ = rdispls[i_rank]; // + recvcounts[i_rank] *
  unsigned char *origin_addr  = (unsigned char *) recvbuf + origin_displ * origin_data_size;
  int            origin_count = recvcounts[i_rank];
  MPI_Get(origin_addr,
          origin_count,
          _pdm_mpi_2_mpi_datatype(recvtype),
          i_rank,
          (MPI_Aint) target_disp[i_rank],
          origin_count,
          _pdm_mpi_2_mpi_datatype(sendtype),
          _pdm_mpi_2_mpi_win(win_send));

  // dt = MPI_Wtime() - t1;

  PDM_UNUSED(win_recv  );
  PDM_UNUSED(sendcounts);
  PDM_UNUSED(sendbuf   ); // Implicit in win_send

  free(target_disp);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_allocate (wrapping de la fonction MPI_Win_allocate)
 *
 *----------------------------------------------------------------------------*/
int PDM_MPI_Win_allocate(PDM_MPI_Aint  size,
                         int           disp_unit,
                         PDM_MPI_Comm  comm,
                         void         *baseptr,
                         PDM_MPI_Win  *win)
{
  MPI_Win _mpi_win;
  int code = MPI_Win_allocate(size,
                              disp_unit,
                              MPI_INFO_NULL,
                              _pdm_mpi_2_mpi_comm(comm),
                              baseptr,
                              &_mpi_win);
  *win = _mpi_2_pdm_mpi_win_add(_mpi_win);

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_free (wrapping de la fonction MPI_Win_free)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Win_free(PDM_MPI_Win *win)

{
  MPI_Win _win = _pdm_mpi_2_mpi_win(*win);
  int code = MPI_Win_free(&_win);

  free(mpi_win[*win]);
  mpi_win[*win] = NULL;
  n_mpi_win += -1;
  *win = PDM_MPI_WIN_NULL;

  if (n_mpi_win == 0) {
    free(mpi_win);
    mpi_win = NULL;

    l_mpi_win = 0;
  }

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Win_fence (wrapping de la fonction MPI_Win_fence)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Win_fence(int assert, PDM_MPI_Win win)

{
  MPI_Win _win = _pdm_mpi_2_mpi_win(win);
  int code = MPI_Win_fence(assert, _win);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Error_string (wrapping de la fonction MPI_Error_string)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Error_string(int errorcode, char *string, int *resultlen)
{
   int code = MPI_Error_string(mpi_err[errorcode], string, resultlen);
   return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_rank (wrapping de la fonction MPI_Comm_rank)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_rank(PDM_MPI_Comm comm, int *rank)
{
  int code = MPI_Comm_rank(_pdm_mpi_2_mpi_comm(comm), rank);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_size (wrapping de la fonction MPI_Comm_size)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_size(PDM_MPI_Comm comm, int *size)
{
  int code = MPI_Comm_size(_pdm_mpi_2_mpi_comm(comm), size);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_get_max_error_string
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_get_max_error_string(void)
{
  return MPI_MAX_ERROR_STRING;
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_free
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_free(PDM_MPI_Comm *comm)
{
 int code = 0;
  if ((*comm != PDM_MPI_COMM_NULL) || (*comm != PDM_MPI_COMM_WORLD)) {

    MPI_Comm mpi_comm_loc = _pdm_mpi_2_mpi_comm(*comm);
    code = MPI_Comm_free(&mpi_comm_loc);

    free(mpi_comm[*comm]);
    mpi_comm[*comm] = NULL;
    n_mpi_comm += -1;
    return _mpi_2_pdm_mpi_err(code);
  }
  *comm = PDM_MPI_COMM_NULL;
  return code;
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_split(PDM_MPI_Comm comm, int color, int key, PDM_MPI_Comm *newcomm)
{
  MPI_Comm _newcomm;
  int code = MPI_Comm_split(_pdm_mpi_2_mpi_comm(comm), color, key, &_newcomm);
  *newcomm = _mpi_2_pdm_mpi_comm(_newcomm);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_dup
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Comm_dup(PDM_MPI_Comm comm, PDM_MPI_Comm *newcomm)
{
  MPI_Comm _newcomm;
  int code = MPI_Comm_dup(_pdm_mpi_2_mpi_comm(comm), &_newcomm);
  *newcomm = _mpi_2_pdm_mpi_comm(_newcomm);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split_type_numa // Non portable mettre un ifdef
 *
 *----------------------------------------------------------------------------*/
int
PDM_MPI_Comm_split_type_numa
(
 PDM_MPI_Comm comm,
 PDM_MPI_Comm *comm_numa
)
{
  PDM_MPI_Comm comm_node;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_SHARED, &comm_node);

  int i_rank_node;
  PDM_MPI_Comm_rank(comm_node, &i_rank_node);

  int i_cpu;
  int i_numa;
#ifdef __linux__
  syscall(SYS_getcpu, &i_cpu, &i_numa, NULL);
#else
  printf("PDM_MPI_Comm_split_type_numa : appel a SYS_getcpu commente car non portable : a reintroduire après tests dans CMake\n");
  abort();
#endif


  /* Sur le shared on split par numa */
  int code = PDM_MPI_Comm_split(comm_node, i_numa, i_rank_node, comm_numa);

  // *comm_numa = comm_node;
  // PDM_MPI_Comm_free(&comm_node);
  // int code = 0;
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Comm_split_type
 *
 *-------------------s---------------------------------------------------------*/
int PDM_MPI_Comm_split_type(PDM_MPI_Comm comm, int split_type, PDM_MPI_Comm *newcomm)
{
  int i_rank;
  MPI_Comm_rank(_pdm_mpi_2_mpi_comm(comm), &i_rank);

  // PDM_MPI_Comm _newcomm;
  int code = 0;
  if(split_type == PDM_MPI_SPLIT_SHARED) {
    MPI_Comm comm_shared;
    code = MPI_Comm_split_type(_pdm_mpi_2_mpi_comm(comm), MPI_COMM_TYPE_SHARED, i_rank /* Key */,
                               MPI_INFO_NULL, &comm_shared);
    *newcomm = _mpi_2_pdm_mpi_comm(comm_shared);
  } else if(split_type == PDM_MPI_SPLIT_NUMA) {
    PDM_MPI_Comm_split_type_numa(comm, newcomm);
  } else {
    PDM_error(__FILE__, __LINE__, 0,"PDM_MPI_Comm_split_type :"
            " split_type '%d' non valide\n", split_type);
    abort();
  }
  // *newcomm = _mpi_2_pdm_mpi_comm(_newcomm);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_allocate_shared_get
 *
 *----------------------------------------------------------------------------*/
PDM_mpi_win_shared_t*
PDM_mpi_win_shared_create(PDM_MPI_Aint size,
                          int          disp_unit,
                          PDM_MPI_Comm comm)
{
  PDM_mpi_win_shared_t* wins = (PDM_mpi_win_shared_t*) malloc(sizeof(PDM_mpi_win_shared_t));

  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  MPI_Info info;
  MPI_Info_create( &info );
  // MPI_Info_set(info, "no_locks", "true");
  // MPI_Info_set( info, "alloc_shared_noncontig", "true" );

  wins->win = MPI_WIN_NULL;
  wins->ptr = NULL;
  int res = 0;
  if(i_rank == 0) {
    res = MPI_Win_allocate_shared(size * disp_unit, disp_unit, info, _pdm_mpi_2_mpi_comm(comm), &wins->ptr , &wins->win);
  } else {
    res = MPI_Win_allocate_shared(0, disp_unit, info , _pdm_mpi_2_mpi_comm(comm), &wins->ptr , &wins->win );
    MPI_Aint size_0;
    int disp_0;
    MPI_Win_shared_query(wins->win, 0, &size_0, &disp_0, &wins->ptr);
  }
  MPI_Info_free(&info);
  assert(res == PDM_MPI_SUCCESS);
  return wins;
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_get
 *
 *----------------------------------------------------------------------------*/
void* PDM_mpi_win_shared_get(PDM_mpi_win_shared_t *wins){
  return wins->ptr;
}


/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_free
 *
 *----------------------------------------------------------------------------*/
void PDM_mpi_win_shared_free(PDM_mpi_win_shared_t *wins){
  MPI_Win_free(&wins->win);
  wins->ptr = NULL;
  free(wins);
}


/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_lock_all
 *
 *----------------------------------------------------------------------------*/
int PDM_mpi_win_shared_lock_all(int assert, PDM_mpi_win_shared_t* win)
{
  int code = MPI_Win_lock_all(assert, win->win);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_unlock_all
 *
 *----------------------------------------------------------------------------*/
int PDM_mpi_win_shared_unlock_all(PDM_mpi_win_shared_t* win)
{
  int code = MPI_Win_unlock_all(win->win);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_mpi_win_shared_sync
 *
 *----------------------------------------------------------------------------*/
int PDM_mpi_win_shared_sync(PDM_mpi_win_shared_t* win)
{
  int code = MPI_Win_sync(win->win);
  return _mpi_2_pdm_mpi_err(code);
}

// ------------------------------------------------------------------
PDM_MPI_Comm PDM_MPI_get_group_of_master(PDM_MPI_Comm comm, PDM_MPI_Comm sub_comm)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int i_rank_sub;
  PDM_MPI_Comm_rank(sub_comm, &i_rank_sub);

  PDM_MPI_Comm master_of_sub_comm;
  int res = 0;
  if(i_rank_sub == 0){
    res = PDM_MPI_Comm_split(comm, 0, i_rank, &master_of_sub_comm);
  } else {
    res = PDM_MPI_Comm_split(comm, PDM_MPI_UNDEFINED, i_rank, &master_of_sub_comm);
    assert(master_of_sub_comm == PDM_MPI_COMM_NULL);
  }
  assert(res == PDM_MPI_SUCCESS);
  return master_of_sub_comm;
}

// ------------------------------------------------------------------
int PDM_MPI_Comm_get_attr_tag_ub(PDM_MPI_Comm comm, void *attribute_val, int *flag)
{
  int code = MPI_Comm_get_attr(_pdm_mpi_2_mpi_comm(comm), MPI_TAG_UB, attribute_val, flag);
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_rand_tag_get
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Rand_tag (PDM_MPI_Comm comm)
{
  struct timeval t;
  gettimeofday(&t, NULL);

  long ltag = t.tv_usec + 1000000 * t.tv_sec;

  MPI_Bcast (&ltag, 1, MPI_LONG, 0, _pdm_mpi_2_mpi_comm(comm));

  void  *max_tag_tmp;
  int flag;

  // Mandatory to call with PDM_MPI_COMM_WORLD becuase only this one keep attributes (openMPI implemntation for exemple)
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &max_tag_tmp, &flag);
  long max_tag = (long) (*((int *) max_tag_tmp));

  // printf("max_tag = %li | ltag = %li \n", max_tag, ltag);

  return (int) (ltag % max_tag);
}



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
                                       PDM_MPI_Comm *newcomm)
{
  MPI_Comm _newcomm;
  const int *weight_in  = MPI_UNWEIGHTED;
  const int *weight_out = MPI_UNWEIGHTED;
  int code = MPI_Dist_graph_create_adjacent(_pdm_mpi_2_mpi_comm(comm_old),
                                            indegree,
                                            sources,
                                            weight_in,
                                            outdegree,
                                            destinations,
                                            weight_out,
                                            MPI_INFO_NULL,
                                            reorder,
                                            &_newcomm);

  *newcomm = _mpi_2_pdm_mpi_comm(_newcomm);
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Allgather (wrapping de la fonction MPI_Allgather)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_allgather(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                               void *recvbuf, int recvcount,
                               PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code =  MPI_Neighbor_allgather(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                                     recvbuf, recvcount,
                                     _pdm_mpi_2_mpi_datatype(recvtype),
                                     _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_allgatherv (wrapping de la fonction MPI_Neighbor_allgatherv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_allgatherv(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                                void *recvbuf, int *recvcounts,
                                int *displs, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Neighbor_allgatherv(sendbuf, sendcount, _pdm_mpi_2_mpi_datatype(sendtype),
                                     recvbuf, recvcounts, displs,
                                     _pdm_mpi_2_mpi_datatype(recvtype), _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}


/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_alltoall (wrapping de la fonction MPI_Neighbor_alltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                              void *recvbuf, int recvcount,
                              PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Neighbor_alltoall(sendbuf, sendcount,
                                   _pdm_mpi_2_mpi_datatype(sendtype),
                                   recvbuf, recvcount,
                                   _pdm_mpi_2_mpi_datatype(recvtype),
                                   _pdm_mpi_2_mpi_comm(comm));
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Ialltoall (wrapping de la fonction MPI_Ialltoall)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ineighbor_alltoall(void *sendbuf, int sendcount, PDM_MPI_Datatype sendtype,
                               void *recvbuf, int recvcount,
                               PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm, PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;

  int code = MPI_Ineighbor_alltoall(sendbuf, sendcount,
                                    _pdm_mpi_2_mpi_datatype(sendtype),
                                    recvbuf, recvcount,
                                    _pdm_mpi_2_mpi_datatype(recvtype),
                                    _pdm_mpi_2_mpi_comm(comm), &_mpi_request);
  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);
  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Neighbor_alltoallv (wrapping de la fonction MPI_Neighbor_alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Neighbor_alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                               PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                               int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm)
{
  int code = MPI_Neighbor_alltoallv(sendbuf,
                           sendcounts,
                           sdispls,
                           _pdm_mpi_2_mpi_datatype(sendtype),
                           recvbuf,
                           recvcounts,
                           rdispls,
                           _pdm_mpi_2_mpi_datatype(recvtype),
                           _pdm_mpi_2_mpi_comm(comm));

  return _mpi_2_pdm_mpi_err(code);
}

/*----------------------------------------------------------------------------
 * PDM_MPI_Ineighbor_alltoallv (wrapping de la fonction MPI_Ineighbor_alltoallv)
 *
 *----------------------------------------------------------------------------*/

int PDM_MPI_Ineighbor_alltoallv(void *sendbuf, int *sendcounts, int *sdispls,
                                PDM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts,
                                int *rdispls, PDM_MPI_Datatype recvtype, PDM_MPI_Comm comm,
                                PDM_MPI_Request *request)
{
  MPI_Request _mpi_request = MPI_REQUEST_NULL;
  // double t1 = MPI_Wtime();
  int code = MPI_Ineighbor_alltoallv(sendbuf,
                                     sendcounts,
                                     sdispls,
                                     _pdm_mpi_2_mpi_datatype(sendtype),
                                     recvbuf,
                                     recvcounts,
                                     rdispls,
                                     _pdm_mpi_2_mpi_datatype(recvtype),
                                     _pdm_mpi_2_mpi_comm(comm), &_mpi_request);

  *request = _mpi_2_pdm_mpi_request_add(_mpi_request);

  return _mpi_2_pdm_mpi_err(code);
}

void
PDM_MPI_setup_hybrid_dist_comm_graph
(
  PDM_MPI_Comm   comm,
  PDM_MPI_Comm  *comm_shared_out,
  PDM_MPI_Comm  *comm_dist_graph_out,
  int           *n_degree,
  int          **neighbor
)
{
  int i_rank;
  int n_rank;

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  // Shared
  PDM_MPI_Comm comm_shared;
  PDM_MPI_Comm_split_type(comm, PDM_MPI_SPLIT_NUMA, &comm_shared);

  int n_rank_in_shm, i_rank_in_shm;
  PDM_MPI_Comm_rank (comm_shared, &i_rank_in_shm);
  PDM_MPI_Comm_size (comm_shared, &n_rank_in_shm);

  PDM_MPI_Comm comm_master_of_shm = PDM_MPI_get_group_of_master(comm, comm_shared);

  int i_rank_master_of_shm = -1;
  int n_rank_master_of_shm;
  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Comm_rank(comm_master_of_shm, &i_rank_master_of_shm);
    PDM_MPI_Comm_size(comm_master_of_shm, &n_rank_master_of_shm);
  }
  PDM_MPI_Bcast(&n_rank_master_of_shm, 1, PDM_MPI_INT, 0, comm_shared);
  PDM_MPI_Bcast(&i_rank_master_of_shm, 1, PDM_MPI_INT, 0, comm_shared);

  PDM_mpi_win_shared_t* wnuma_by_numa_n = PDM_mpi_win_shared_create(n_rank_master_of_shm, sizeof(int), comm_shared);
  int *numa_by_numa_n  = PDM_mpi_win_shared_get(wnuma_by_numa_n);
  PDM_mpi_win_shared_lock_all (0, wnuma_by_numa_n);


  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    PDM_MPI_Allgather(&n_rank_in_shm, 1, PDM_MPI_INT,
                      numa_by_numa_n, 1, PDM_MPI_INT, comm_master_of_shm);
  }
  PDM_mpi_win_shared_sync(wnuma_by_numa_n);
  PDM_MPI_Barrier(comm_shared);

  int n_tot_numa = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    n_tot_numa += numa_by_numa_n[i];
  }
  /*
   * Create idx  and  gid of each numa
   */
  PDM_mpi_win_shared_t* wnuma_core_gid    = PDM_mpi_win_shared_create(n_tot_numa               , sizeof(int), comm_shared);
  PDM_mpi_win_shared_t* wnuma_by_numa_idx = PDM_mpi_win_shared_create(n_rank_master_of_shm+1, sizeof(int), comm_shared);
  int *numa_core_gid    = PDM_mpi_win_shared_get(wnuma_core_gid);
  int *numa_by_numa_idx = PDM_mpi_win_shared_get(wnuma_by_numa_idx);
  PDM_mpi_win_shared_lock_all (0, wnuma_core_gid);
  PDM_mpi_win_shared_lock_all (0, wnuma_by_numa_idx);


  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    numa_by_numa_idx[0] = 0;
    for(int i = 0; i < n_rank_master_of_shm; ++i) {
      numa_by_numa_idx[i+1] = numa_by_numa_idx[i] + numa_by_numa_n[i];
    }
  }
  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_by_numa_idx);

  numa_core_gid[numa_by_numa_idx[i_rank_master_of_shm]+i_rank_in_shm] = i_rank;

  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_core_gid);

  /*
   *  Exchange of the global numbering of rank for each NUMA
   */
  if(comm_master_of_shm != PDM_MPI_COMM_NULL) {
    int *lnuma_core_gid = malloc(n_rank_in_shm * sizeof(int));
    for(int i = 0; i < n_rank_in_shm; ++i) {
      lnuma_core_gid[i] = numa_core_gid[numa_by_numa_idx[i_rank_master_of_shm]+i];
    }
    PDM_MPI_Allgatherv(lnuma_core_gid, n_rank_in_shm, PDM_MPI_INT,
                       numa_core_gid , numa_by_numa_n, numa_by_numa_idx, PDM_MPI_INT, comm_master_of_shm);
    free(lnuma_core_gid);
  }
  PDM_MPI_Barrier(comm_shared);
  PDM_mpi_win_shared_sync(wnuma_core_gid);

  /*
   * Computation of degree_in
   */
  int *send_n   = malloc(  n_rank    * sizeof(int));
  int *recv_n   = malloc(  n_rank    * sizeof(int));
  int *send_idx = malloc( (n_rank+1) * sizeof(int));
  int *recv_idx = malloc( (n_rank+1) * sizeof(int));

  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
    recv_n[i] = 0;
  }

  int n_degrees_in = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int lid_rank = (j - numa_by_numa_idx[i]) % n_rank_in_shm; // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        n_degrees_in++;
      }
    }
  }

  int* neighbor_in = malloc( (n_degrees_in ) * sizeof(int));
  n_degrees_in = 0;
  for(int i = 0; i < n_rank_master_of_shm; ++i) {
    for(int j = numa_by_numa_idx[i]; j < numa_by_numa_idx[i+1]; ++j) {
      int gid_rank = numa_core_gid[j];
      int lid_rank = (j - numa_by_numa_idx[i]) % n_rank_in_shm;  // Donc numero de numa dans le group
      if(lid_rank == i_rank_in_shm){
        neighbor_in[n_degrees_in++] = gid_rank;
      }
    }
  }


  for(int i = 0; i < n_degrees_in; ++i) {
    send_n[neighbor_in[i]]++;
  }

  send_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    send_n[i] = 0;
  }

  int *send_cur_i_rank = malloc(send_idx[n_rank] * sizeof(int));

  for(int i = 0; i < n_degrees_in; ++i) {
    int idx_write = send_idx[neighbor_in[i]] + send_n[neighbor_in[i]]++;
    send_cur_i_rank[idx_write] = i_rank;
  }


  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }
  int *recv_opp_i_rank = malloc(recv_idx[n_rank] * sizeof(int));

  PDM_MPI_Alltoallv(send_cur_i_rank, send_n, send_idx, PDM_MPI_INT,
                    recv_opp_i_rank, recv_n, recv_idx, PDM_MPI_INT, comm);


  int n_degrees_out = recv_idx[n_rank];
  int *neighbor_out = recv_opp_i_rank; // Already sort normaly

  free(send_n);
  free(recv_n);
  free(send_idx);
  free(recv_idx);
  free(send_cur_i_rank);


  PDM_MPI_Comm comm_dist_graph;
  PDM_MPI_Dist_graph_create_adjacent(comm,
                                     n_degrees_in,
                                     neighbor_in,
                                     n_degrees_out,
                                     neighbor_out,
                                     0,
                                     &comm_dist_graph);

  PDM_mpi_win_shared_unlock_all(wnuma_by_numa_n);
  PDM_mpi_win_shared_unlock_all(wnuma_core_gid);
  PDM_mpi_win_shared_unlock_all(wnuma_by_numa_idx);
  PDM_mpi_win_shared_free(wnuma_by_numa_n);
  PDM_mpi_win_shared_free(wnuma_core_gid);
  PDM_mpi_win_shared_free(wnuma_by_numa_idx);

  free(recv_opp_i_rank);

  *comm_shared_out     = comm_shared;
  *comm_dist_graph_out = comm_dist_graph;

  *n_degree = n_degrees_in;
  *neighbor = neighbor_in;
}




int
PDM_MPI_Dist_graph_neighbors_count
(
  PDM_MPI_Comm  comm,
  int          *n_degree_in,
  int          *n_degree_out,
  int          *is_weighted
)
{
  int code = MPI_Dist_graph_neighbors_count(_pdm_mpi_2_mpi_comm(comm),
                                            n_degree_in,
                                            n_degree_out,
                                            is_weighted);
  return _mpi_2_pdm_mpi_err(code);
}



int
PDM_MPI_Dist_graph_neighbors
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *sources,
  int            n_degree_out,
  int           *destinations
)
{

  int *weight_in  = NULL;
  int *weight_out = NULL;
  int code = MPI_Dist_graph_neighbors(_pdm_mpi_2_mpi_comm(comm),
                                      n_degree_in,
                                      sources,
                                      weight_in,
                                      n_degree_out,
                                      destinations,
                                      weight_out);
  return _mpi_2_pdm_mpi_err(code);
}




void
PDM_MPI_setup_dist_graph_from_neighbor_in
(
  PDM_MPI_Comm   comm,
  int            n_degree_in,
  int           *neighbor_in,
  PDM_MPI_Comm  *comm_dist_graph_out
)
{
  int i_rank;
  int n_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  int *send_n   = malloc(  n_rank    * sizeof(int));
  int *recv_n   = malloc(  n_rank    * sizeof(int));
  int *send_idx = malloc( (n_rank+1) * sizeof(int));
  int *recv_idx = malloc( (n_rank+1) * sizeof(int));

  for(int i = 0; i < n_rank; ++i) {
    send_n[i] = 0;
    recv_n[i] = 0;
  }

  PDM_MPI_Comm_rank(comm, &i_rank);
  PDM_MPI_Comm_size(comm, &n_rank);

  for(int i = 0; i < n_degree_in; ++i) {
    send_n[neighbor_in[i]]++;
  }

  send_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    send_idx[i+1] = send_idx[i] + send_n[i];
    send_n[i] = 0;
  }

  int *send_cur_i_rank = malloc(send_idx[n_rank] * sizeof(int));

  for(int i = 0; i < n_degree_in; ++i) {
    int idx_write = send_idx[neighbor_in[i]] + send_n[neighbor_in[i]]++;
    send_cur_i_rank[idx_write] = i_rank;
  }


  PDM_MPI_Alltoall(send_n, 1, PDM_MPI_INT,
                   recv_n, 1, PDM_MPI_INT, comm);

  recv_idx[0] = 0;
  for(int i = 0; i < n_rank; ++i) {
    recv_idx[i+1] = recv_idx[i] + recv_n[i];
  }
  int *recv_opp_i_rank = malloc(recv_idx[n_rank] * sizeof(int));

  PDM_MPI_Alltoallv(send_cur_i_rank, send_n, send_idx, PDM_MPI_INT,
                    recv_opp_i_rank, recv_n, recv_idx, PDM_MPI_INT, comm);


  int n_degrees_out = recv_idx[n_rank];
  int *neighbor_out = recv_opp_i_rank; // Already sort normaly

  free(send_n);
  free(recv_n);
  free(send_idx);
  free(recv_idx);
  free(send_cur_i_rank);

  PDM_MPI_Dist_graph_create_adjacent(comm,
                                     n_degree_in,
                                     neighbor_in,
                                     n_degrees_out,
                                     neighbor_out,
                                     0,
                                     comm_dist_graph_out);
}




#ifdef __cplusplus
}
#endif /* __cplusplus */
