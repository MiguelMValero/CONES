/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_line.h"
#include "pdm_plane.h"
#include "pdm_edges_intersect.h"
#include "pdm_edges_intersect_priv.h"
#include "pdm_hash_tab.h"
#include "pdm_part_to_block.h"
#include "pdm_block_to_part.h"
#include "pdm_mpi.h"
#include "pdm_morton.h"
#include "pdm_array.h"
#include "pdm_printf.h"
#include "pdm_error.h"

/*=============================================================================
 * Macro definitions
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Type
 *============================================================================*/


/*=============================================================================
 * Static global variables
 *============================================================================*/


static const   char* _typeInter[] = {"PDM_LINE_INTERSECT_UNDEF",  /*!< No intersection */
                                     "PDM_LINE_INTERSECT_NO",       /*!< No intersection */
                                     "PDM_LINE_INTERSECT_YES",       /*!< Intersection */
                                     "PDM_LINE_INTERSECT_ON_LINE"};  /*!< On line  */

static const   char* _typeInter2[] = {"PDM_EDGES_INTERSECT_POINT_NEW",
                                      "PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB",
                                      "PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA",
                                      "PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB"};

/*=============================================================================
 * Static function definitions
 *============================================================================*/


/**
 *
 * \brief Replace default C modulo function (%) to taking into account
 *
 * \param [in]   val    Value
 * \param [in]   mod    Mod
 *
 * \return modulo
 */

static int
_modulo
(
 int val,
 int mod
)
{
  if (val >= 0) {
    return val % mod;
  }
  else {
    return val + mod * ((mod - val - 1)/mod);
  }
}


/**
 *
 * \brief Create a new \ref _edges_intersect_res_t object
 *
 * \param [in]   nGEdgeA    Global number mesh A edge
 * \param [in]   nGEdgeB    Global number mesh B edge
 * \param [in]   nNewPointA Number of new point for mesh A
 * \param [in]   nNewPointB Number of new point for mesh B
 *
 * \return      A new \ref _edges_intersect_res_t
 */

static _edges_intersect_res_t *
_edges_intersect_res_create
(
const PDM_g_num_t   nGEdgeA,
const PDM_g_num_t   nGEdgeB,
const int          nNewPointsA,
const int          nNewPointsB
)
{
  _edges_intersect_res_t *newInter = malloc (sizeof(_edges_intersect_res_t));

  newInter->nGEdgeA     = nGEdgeA;
  newInter->nGEdgeB     = nGEdgeB;
  newInter->originEdgeA = -1;
  newInter->originEdgeB = -1;
  newInter->endEdgeA = -1;
  newInter->endEdgeB = -1;
  newInter->tIntersect  = PDM_LINE_INTERSECT_UNDEF;

  newInter->nNewPointsA = nNewPointsA;
  newInter->uA          = malloc (sizeof(double) * nNewPointsA);
  newInter->coordsA     = malloc (sizeof(double) * 3 * nNewPointsA);
  newInter->linkA       = malloc (sizeof(PDM_g_num_t) * nNewPointsA);
  newInter->gNumA       = malloc (sizeof(PDM_g_num_t) * nNewPointsA);
  newInter->oNewPointsA = malloc (sizeof(PDM_edges_intersect_point_t) * nNewPointsA);

  newInter->nNewPointsB = nNewPointsB;
  newInter->uB          = malloc (sizeof(double) * nNewPointsB);
  newInter->coordsB     = malloc (sizeof(double) * 3 * nNewPointsB);
  newInter->linkB       = malloc (sizeof(PDM_g_num_t) * nNewPointsB);
  newInter->gNumB       = malloc (sizeof(PDM_g_num_t) * nNewPointsB);
  newInter->oNewPointsB = malloc (sizeof(PDM_edges_intersect_point_t) * nNewPointsB);

  return newInter;
}


/**
 *
 * \brief Free a \ref _edges_intersect_res_t object
 *
 * \param [in]   eir  Edges intersection results object
 *
 * \return NULL
 */

static _edges_intersect_res_t *
_edges_intersect_res_free (_edges_intersect_res_t *eir)
{
  if (eir != NULL)
  {
	if (eir->linkA != NULL) {
		free (eir->gNumA);
		free (eir->linkA);
		free (eir->uA);
		free (eir->coordsA);
		free (eir->oNewPointsA);
		eir->gNumA      = NULL;
		eir->linkA      = NULL;
		eir->uA           = NULL;
		eir->oNewPointsA  = NULL;
		eir->coordsA      = NULL;
	}
	eir->nNewPointsB      = 0;
	if (eir->linkB != NULL) {
		free (eir->gNumB);
		free (eir->linkB);
		free (eir->uB);
		free (eir->coordsB);
		free (eir->oNewPointsB);
		eir->gNumB      = NULL;
		eir->linkB      = NULL;
		eir->uB           = NULL;
		eir->oNewPointsB  = NULL;
		eir->coordsB      = NULL;
	}
	free (eir);
  }
  return NULL;
}

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Create a new \ref PDM_edges_intersect_t object
 *
 * \param [in]   maxGNEdgeA    Max global number of edges in mesh A
 * \param [in]   maxGNEdgeN    Max global number of edges in mesh B
 * \param [in]   vtxCarLengthTol Absolute tolerance for characteristic length
 * \param [in]   sMSGComm        size of mpicomm
 *
 * \return      A new \ref PDM_edges_intersect_t
 */

PDM_edges_intersect_t *
PDM_edges_intersect_create
(
const PDM_g_num_t maxGNEdgeA,
const PDM_g_num_t maxGNEdgeB,
const double     vtxCarLengthTol,
const PDM_MPI_Comm   comm
)
{
  int sMSGComm;
  PDM_MPI_Comm_size (comm, &sMSGComm);

  _edges_intersect_t *ei = malloc( sizeof(_edges_intersect_t));

  PDM_g_num_t _keyMax = (maxGNEdgeA + maxGNEdgeB) / sMSGComm + 1;

  int keyMax = (int) _keyMax;

  ei->ht = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                (void *) &keyMax);

  PDM_g_num_t _keyMaxA = maxGNEdgeA / sMSGComm +1;
  int keyMaxA = (int) _keyMaxA;
  ei->htA = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                (void *) &keyMaxA);

  PDM_g_num_t _keyMaxB = maxGNEdgeB / sMSGComm + 1;
  int keyMaxB = (int) _keyMaxB;
  ei->htB = PDM_hash_tab_create (PDM_HASH_TAB_KEY_INT,
                                (void *) &keyMaxB);

  ei->maxGNEdgeA = maxGNEdgeA;
  ei->maxGNEdgeB = maxGNEdgeB;
  ei->vtxCarLengthTol = vtxCarLengthTol;
  ei->comm = comm;
  ei->sMSGComm = sMSGComm;

  return (PDM_edges_intersect_t *) ei;
}


/**
 *
 * \brief Get result of the intersection
 *
 * \param [in]   ei               Current edges intersection pointer
 * \param [in]   get_t            Type of key to return data
 * \param [in]   nGEdgeA          Global number of meshA edge
 * \param [in]   nGEdgeB          Global number of meshB edge
 * \param [out]  n_intersect      Number of intersections
 *
 * \return    Result of the intersection
 */

PDM_edges_intersect_res_t **
PDM_edges_intersect_get
(
PDM_edges_intersect_t  *ei,
PDM_edges_get_t         get_t,
const PDM_g_num_t       nGEdgeA,
const PDM_g_num_t       nGEdgeB,
int                    *n_intersect
)
{

  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;
  *n_intersect = 0;

  if (get_t == PDM_EDGES_GET_FROM_AB) {

    PDM_hash_tab_t *ht = _ei->ht;
    PDM_g_num_t _key = (nGEdgeA + nGEdgeB) / _ei->sMSGComm;
    int key = (int) _key;

    const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    _edges_intersect_res_t ** datas =
              (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht,
                                                                 (void *) &key);
    /**********************************************
     * Look for requested intersection            *
     **********************************************/

    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if ((data->nGEdgeA == nGEdgeA) && (data->nGEdgeB == nGEdgeB)) {
        *n_intersect = 1;
        _edges_intersect_res_t **_datas_edge =
            malloc(sizeof(_edges_intersect_res_t *));
        _datas_edge[0] = data;
        return (PDM_edges_intersect_res_t **) _datas_edge;
      }
    }

    return NULL;

  }

  else if (get_t == PDM_EDGES_GET_FROM_A) {

    PDM_hash_tab_t *ht = _ei->htA;
    PDM_g_num_t _key = nGEdgeA / _ei->sMSGComm;
    int key = (int) _key;

    const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    _edges_intersect_res_t ** datas =
              (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht,
                                                                 (void *) &key);
    /**********************************************
     * Look for requested intersection            *
     **********************************************/

    int _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data == NULL) {
        PDM_error(__FILE__, __LINE__, 0, "les donnees de la table de hashage ne sont pas renseignees\n");
        abort();
      }
      if (data->nGEdgeA == nGEdgeA) {
        _nData += 1;
      }
    }

    _edges_intersect_res_t **_datas_edge =
            malloc(sizeof(_edges_intersect_res_t *) * _nData);

    _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeA == nGEdgeA) {
        _datas_edge[_nData++] = data;
      }
    }
    *n_intersect = _nData;
    return (PDM_edges_intersect_res_t **) _datas_edge;
  }

  else if (get_t == PDM_EDGES_GET_FROM_B) {

    PDM_hash_tab_t *ht = _ei->htB;
    PDM_g_num_t _key = nGEdgeB / _ei->sMSGComm;
    int key = (int) _key;

    const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    _edges_intersect_res_t ** datas =
              (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht,
                                                                 (void *) &key);
    /**********************************************
     * Look for requested intersection            *
     **********************************************/

    int _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeB == nGEdgeB) {
        _nData += 1;
      }
    }

    _edges_intersect_res_t **_datas_edge =
            malloc(sizeof(_edges_intersect_res_t *) * _nData);

    _nData = 0;
    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *data = datas[i];
      if (data->nGEdgeB == nGEdgeB) {
        _datas_edge[_nData++] = data;
      }
    }
    *n_intersect = _nData;
    return (PDM_edges_intersect_res_t **) _datas_edge;
  }

  return NULL;
}


/**
 *
 * \brief Perform an intersection between a meshA edge and a meshB edge
 *
 * \param [in]   ei               Current edges intersection pointer
 * \param [in]   nGEdgeA          Global number of meshA edge
 * \param [in]   nGVtxA           Global number of edgeA vertices
 * \param [in]   charLgthVtxA[2]  Characteristic length of edgeA vertices
 * \param [in]   coordsVtxA[6]    Coordinates of edgeA vertices
 * \param [in]   nGEdgeB          Global number of meshB edge
 * \param [in]   nGVtxB           Global number of edgeB vertices
 * \param [in]   charLgthVtxB[2]  Characteristic length of edgeB vertices
 * \param [in]   coordsVtxB[6]    Coordinates of edgeB vertices
 *
 * \return    Result of the intersection
 */

PDM_edges_intersect_res_t *
PDM_edges_intersect_add
(
PDM_edges_intersect_t       *ei,
const PDM_g_num_t             nGEdgeA,
const PDM_g_num_t             nGVtxA[2],
const double                 charLgthVtxA[2],
const double                 coordsVtxA[6],
const PDM_g_num_t             nGEdgeB,
const PDM_g_num_t             nGVtxB[2],
const double                 charLgthVtxB[2],
const double                 coordsVtxB[6]
)
{


  int vb = 0;//1;
  if (vb)   {
    PDM_printf ("==== PDM_edges_intersect_add ==== \n");
	  PDM_printf ("--- nGEdgeA:"PDM_FMT_G_NUM", nGVtxA:"PDM_FMT_G_NUM"-"PDM_FMT_G_NUM" \n charLgthVtxA:%12.5e-%12.5e, coordsVtxA:%12.5e-%12.5e-%12.5e-%12.5e-%12.5e-%12.5e\n",
                nGEdgeA, nGVtxA[0], nGVtxA[1],
                charLgthVtxA[0], charLgthVtxA[1], coordsVtxA[0],
                coordsVtxA[1], coordsVtxA[2],
                coordsVtxA[3], coordsVtxA[4], coordsVtxA[5]);
	  PDM_printf ("--- nGEdgeB:"PDM_FMT_G_NUM", nGVtxB:"PDM_FMT_G_NUM"-"PDM_FMT_G_NUM" \n charLgthVtxB:%12.5e-%12.5e, coordsVtxB:%12.5e-%12.5e-%12.5e-%12.5e-%12.5e-%12.5e \n",
                nGEdgeB, nGVtxB[0], nGVtxB[1],
                charLgthVtxB[0], charLgthVtxB[1],
                coordsVtxB[0], coordsVtxB[1], coordsVtxB[2],
                coordsVtxB[3], coordsVtxB[4], coordsVtxB[5]);
  }
  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;

  PDM_hash_tab_t *ht = _ei->ht;
  PDM_hash_tab_t *htA = _ei->htA;
  PDM_hash_tab_t *htB = _ei->htB;

  const double minMin = 1e-12;

  double _charLgthVtxA[2] = {PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxA[0], minMin),
                             PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxA[1], minMin)};

  double _charLgthVtxB[2] = {PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxB[0], minMin),
                             PDM_MIN (_ei->vtxCarLengthTol * charLgthVtxB[1], minMin)};

  PDM_g_num_t _key  = (nGEdgeA + nGEdgeB) / _ei->sMSGComm;
  PDM_g_num_t _keyA = nGEdgeA / _ei->sMSGComm;
  PDM_g_num_t _keyB = nGEdgeB / _ei->sMSGComm;

  int key = (int) _key;
  int keyA = (int) _keyA;
  int keyB = (int) _keyB;

  const int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
  _edges_intersect_res_t ** datas = (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht,(void *) &key);


  /**********************************************
   * Check if intersection is already preformed *
   **********************************************/

  for (int i = 0; i < nData; i++) {
    _edges_intersect_res_t *data = datas[i];
    if ((nGEdgeA == data->nGEdgeA) && (nGEdgeB == data->nGEdgeB)){
      if (vb) {
        PDM_printf ("intersection is already preformed\n");
        PDM_printf ("==== PDM_edges_intersect_add ==== terminated ====\n");
      }
      return (PDM_edges_intersect_res_t *) data;
    }
  }

  /******************************
   * Perform a new intersection *
   ******************************/

  /*
   * Line-line intersection
   */

  double u1;
  double v1;

#if 1
  PDM_line_intersect_t tIntersect =
    PDM_line_intersection_mean_square (coordsVtxA, &(coordsVtxA[3]),
                                       coordsVtxB, &(coordsVtxB[3]),
                                       &u1, &v1);
#else
  // !!! Works only in the xy-plane
  PDM_line_intersect_t tIntersect =
    PDM_line_intersection_2drobust (coordsVtxA, &(coordsVtxA[3]),
                                    coordsVtxB, &(coordsVtxB[3]),
                                    &u1, &v1);
#endif
  /*  PDM_LINE_INTERSECT_UNDEF   = -1,  !< No intersection */
  /*  PDM_LINE_INTERSECT_NO      = 0,  !< No intersection */
  /*  PDM_LINE_INTERSECT_YES     = 1,  !< Intersection */
  /*  PDM_LINE_INTERSECT_ON_LINE = 2,  !< On line  */

  bool isInitialOnLine = false;
  if (tIntersect == PDM_LINE_INTERSECT_ON_LINE) {
    isInitialOnLine = true;
  }

  /*
   * Initialization
   */

  _edges_intersect_res_t *newInter = NULL;

  double vA[3];
  double vB[3];

  for (int i = 0; i < 3; i++) {
    vA[i] = coordsVtxA[3 + i] - coordsVtxA[i];
    vB[i] = coordsVtxB[3 + i] - coordsVtxB[i];
  }

  double vA_norm2 = PDM_DOT_PRODUCT (vA, vA);
  double vB_norm2 = PDM_DOT_PRODUCT (vB, vB);

  double vA_norm = PDM_MODULE (vA);
  double vB_norm = PDM_MODULE (vB);

  bool isSetted = false;

  /*
   * Resolution of inconsistencies
   */

  if ((tIntersect == PDM_LINE_INTERSECT_NO) ||
      (tIntersect == PDM_LINE_INTERSECT_YES)) {

    double A1B1[3] = {coordsVtxB[0] - coordsVtxA[0],
                      coordsVtxB[1] - coordsVtxA[1],
                      coordsVtxB[2] - coordsVtxA[2]};
    double mA1B1 = PDM_MODULE (A1B1);

    double A1B2[3] = {coordsVtxB[3+0] - coordsVtxA[0],
                      coordsVtxB[3+1] - coordsVtxA[1],
                      coordsVtxB[3+2] - coordsVtxA[2]};
    double mA1B2 = PDM_MODULE (A1B2);

    double B1A2[3] = {coordsVtxA[3+0] - coordsVtxB[0],
                      coordsVtxA[3+1] - coordsVtxB[1],
                      coordsVtxA[3+2] - coordsVtxB[2]};
    double mB1A2 = PDM_MODULE (B1A2);

    double B2A2[3] = {coordsVtxA[3+0] - coordsVtxB[3+0],
                      coordsVtxA[3+1] - coordsVtxB[3+1],
                      coordsVtxA[3+2] - coordsVtxB[3+2]};
    double mB2A2 = PDM_MODULE (B2A2);

    int isA1B1 = (mA1B1 < _charLgthVtxA[0]) &&
      (mA1B1 < _charLgthVtxB[0]);

    int isA1B2 = (mA1B2 < _charLgthVtxA[0]) &&
      (mA1B2 < _charLgthVtxB[1]);

    int isA2B1 = (mB1A2 < _charLgthVtxA[1]) &&
      (mB1A2 < _charLgthVtxB[0]);

    int isA2B2 = (mB2A2 < _charLgthVtxA[1]) &&
      (mB2A2 < _charLgthVtxB[1]);


    /* double dVtxA[2] = {PDM_ABS (u1 * vA_norm), */
    /*                    PDM_ABS ((1. - u1) * vA_norm)}; */

    /* double dVtxB[2] = {PDM_ABS (v1 * vB_norm), */
    /*                    PDM_ABS ((1. - v1) * vB_norm)}; */

    /*
     * Check if A vertex and B vertex are the same
     */

    /* int isA1B1 = (dVtxA[0] < _charLgthVtxA[0]) && */
    /*              (dVtxB[0] < _charLgthVtxB[0]); */

    /* int isA1B2 = (dVtxA[0] < _charLgthVtxA[0]) && */
    /*              (dVtxB[1] < _charLgthVtxB[1]); */

    /* int isA2B1 = (dVtxA[1] < _charLgthVtxA[1]) && */
    /*              (dVtxB[0] < _charLgthVtxB[0]); */

    /* int isA2B2 = (dVtxA[1] < _charLgthVtxA[1]) && */
    /*              (dVtxB[1] < _charLgthVtxB[1]); */


    /*
     * Check if A vertex is on B Edge and vice versa
     */

    double closestA1EdgeB[3];
    double tA1EdgeB;
    double d2A1EdgeB = PDM_line_distance (coordsVtxA,
                                          coordsVtxB,
                                          coordsVtxB + 3,
                                          &tA1EdgeB,
                                          closestA1EdgeB);

    double closestA2EdgeB[3];
    double tA2EdgeB;
    double d2A2EdgeB = PDM_line_distance (coordsVtxA + 3,
                                          coordsVtxB,
                                          coordsVtxB + 3,
                                          &tA2EdgeB,
                                          closestA2EdgeB);

    double closestB1EdgeA[3];
    double tB1EdgeA;
    double d2B1EdgeA = PDM_line_distance (coordsVtxB,
                                          coordsVtxA,
                                          coordsVtxA + 3,
                                          &tB1EdgeA,
                                          closestB1EdgeA);

    double closestB2EdgeA[3];
    double tB2EdgeA;
    double d2B2EdgeA = PDM_line_distance (coordsVtxB + 3,
                                          coordsVtxA,
                                          coordsVtxA + 3,
                                          &tB2EdgeA,
                                          closestB2EdgeA);

    double _deltaEdgA = _charLgthVtxA[1] - _charLgthVtxA[0];
    double _deltaEdgB = _charLgthVtxB[1] - _charLgthVtxB[0];

    double A1closestEdgeB[3] = {closestA1EdgeB[0] - coordsVtxA[0],
                                closestA1EdgeB[1] - coordsVtxA[1],
                                closestA1EdgeB[2] - coordsVtxA[2]};
    double dA1closestEdgeB = PDM_MODULE (A1closestEdgeB);

    double A2closestEdgeB[3] = {closestA2EdgeB[0] - coordsVtxA[3+0],
                                closestA2EdgeB[1] - coordsVtxA[3+1],
                                closestA2EdgeB[2] - coordsVtxA[3+2]};
    double dA2closestEdgeB = PDM_MODULE (A2closestEdgeB);


    double B1closestEdgeA[3] = {closestB1EdgeA[0] - coordsVtxB[0],
                                closestB1EdgeA[1] - coordsVtxB[1],
                                closestB1EdgeA[2] - coordsVtxB[2]};
    double dB1closestEdgeA = PDM_MODULE (B1closestEdgeA);

    double B2closestEdgeA[3] = {closestB2EdgeA[0] - coordsVtxB[3+0],
                                closestB2EdgeA[1] - coordsVtxB[3+1],
                                closestB2EdgeA[2] - coordsVtxB[3+2]};
    double dB2closestEdgeA = PDM_MODULE (B2closestEdgeA);

    int isA1OnEdgeB = (d2A1EdgeB < ((_charLgthVtxB[0] + tA1EdgeB * _deltaEdgB) *
                                    (_charLgthVtxB[0] + tA1EdgeB * _deltaEdgB))) &&
                      (dA1closestEdgeB < _charLgthVtxA[0]) && !isA1B1 && !isA1B2;

    int isA2OnEdgeB = (d2A2EdgeB < ((_charLgthVtxB[0] + tA2EdgeB * _deltaEdgB) *
                                    (_charLgthVtxB[0] + tA2EdgeB * _deltaEdgB))) &&
                      (dA2closestEdgeB < _charLgthVtxA[1]) && !isA2B1 && !isA2B2;

    int isB1OnEdgeA = (d2B1EdgeA < ((_charLgthVtxA[0] + tB1EdgeA * _deltaEdgA) *
                                    (_charLgthVtxA[0] + tB1EdgeA * _deltaEdgA))) &&
                      (dB1closestEdgeA < _charLgthVtxB[0]) && !isA1B1 && !isA2B1;

    int isB2OnEdgeA = (d2B2EdgeA < ((_charLgthVtxA[0] + tB2EdgeA * _deltaEdgA) *
                                    (_charLgthVtxA[0] + tB2EdgeA * _deltaEdgA))) &&
                      (dB2closestEdgeA < _charLgthVtxB[1]) && !isA1B2 && !isA2B2;

    if (vb) {
      PDM_printf ("isA1B1:%d isA1B2:%d isA2B1:%d isA2B2:%d isA1OnEdgeB:%d isA2OnEdgeB:%d "
                  "isB1OnEdgeA:%d isB2OnEdgeA:%d\n",
                  isA1B1, isA1B2, isA2B1, isA2B2,
                  isA1OnEdgeB, isA2OnEdgeB, isB1OnEdgeA, isB2OnEdgeA);
    }

    /*
     * Check if A Edge and B are the same line
     */

    if ((isA1B1 && isA2B2) || (isA2B1 && isA1B2)) {
      tIntersect = PDM_LINE_INTERSECT_ON_LINE;
    }

    if ((isA1OnEdgeB && isA2OnEdgeB) || (isB1OnEdgeA && isB2OnEdgeA)) {
      tIntersect = PDM_LINE_INTERSECT_ON_LINE;
    }

    if ((isA1B1 && isA2OnEdgeB) || (isA1B1 && isB2OnEdgeA) ||
        (isA2B1 && isA1OnEdgeB) || (isA2B1 && isB2OnEdgeA) ||
        (isA1B2 && isA2OnEdgeB) || (isA1B2 && isB1OnEdgeA) ||
        (isA2B2 && isA1OnEdgeB) || (isA2B2 && isB1OnEdgeA)) {
      tIntersect = PDM_LINE_INTERSECT_ON_LINE;
    }

    if (tIntersect != PDM_LINE_INTERSECT_ON_LINE) {

      /*
       * Define intersection for same vertex inconsistencies
       */

      if (isA1B1 || isA1B2 || isA2B1 || isA2B2) {

        int nNewPointsA = 1;
        int nNewPointsB = 1;

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

        int iA = 0;
        int iB = 0;

        if (isA1B1) {
          iA = 0;
          iB = 0;
        }
        else if (isA1B2) {
          iA = 0;
          iB = 1;
        }
        else if (isA2B1) {
          iA = 1;
          iB = 0;
        }
        else if (isA2B2) {
          iA = 1;
          iB = 1;
        }

        newInter->uA[0] = (double) iA;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        newInter->uB[0] = (double) iB;
        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }

        newInter->linkA[0] = nGVtxB[iB];
        newInter->linkB[0] = nGVtxA[iA];
        newInter->gNumA[0] = nGVtxA[iA];
        newInter->gNumB[0] = nGVtxB[iB];

        for (int i = 0; i < 3; i++) {
          double _comp = (coordsVtxB[3*iB+i] + coordsVtxA[3*iA+i]) / 2;
          newInter->coordsA[i] = _comp;
          newInter->coordsB[i] = _comp;
        }
        isSetted = true;

      }

      /*
       * Define intersection for vertex on edge inconsistencies
       */

      else if (isA1OnEdgeB) {
        int nNewPointsA = 0;
        int nNewPointsB = 1;

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

        double B1ClosestA1[3] = {tA1EdgeB * vB[0],
                                 tA1EdgeB * vB[1],
                                 tA1EdgeB * vB[2]};

         newInter->uB[0] = PDM_DOT_PRODUCT (B1ClosestA1, vB) / vB_norm2;
        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsB[i] = closestA1EdgeB[i];
        }

        newInter->linkB[0] = nGVtxA[0];
        newInter->gNumB[0] = 0;

        isSetted = true;

      }

      else if (isA2OnEdgeB) {
        int nNewPointsA = 0;
        int nNewPointsB = 1;

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

        double B1ClosestA2[3] = {tA2EdgeB * vB[0],
                                 tA2EdgeB * vB[1],
                                 tA2EdgeB * vB[2]};

        newInter->uB[0] = PDM_DOT_PRODUCT (B1ClosestA2, vB) / vB_norm2;
        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsB[i] = closestA2EdgeB[i];
        }

        newInter->linkB[0] = nGVtxA[1];
        newInter->gNumB[0] = 0;

        isSetted = true;

      }

      else if (isB1OnEdgeA) {

        int nNewPointsA = 1;
        int nNewPointsB = 0;

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;
        double A1ClosestB1[3] = {tB1EdgeA *vA[0],
                                 tB1EdgeA *vA[1],
                                 tB1EdgeA *vA[2]};

        newInter->uA[0] = PDM_DOT_PRODUCT (A1ClosestB1, vA) / vA_norm2;

        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = closestB1EdgeA[i];
        }

        newInter->linkA[0] = nGVtxB[0];
        newInter->gNumA[0] = 0;

        isSetted = true;

      }

      else if (isB2OnEdgeA) {

        int nNewPointsA = 1;
        int nNewPointsB = 0;

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

        double A1ClosestB2[3] = {tB2EdgeA * vA[0],
                                 tB2EdgeA * vA[1],
                                 tB2EdgeA * vA[2]};

        newInter->uA[0] = PDM_DOT_PRODUCT (A1ClosestB2, vA) / vA_norm2;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = closestB2EdgeA[i];
        }

        newInter->linkA[0] = nGVtxB[1];
        newInter->gNumA[0] = 0;

        isSetted = true;
      }
    }
  }

  if (tIntersect == PDM_LINE_INTERSECT_ON_LINE) {

    if (isInitialOnLine) {

      double A1B1[3] = {coordsVtxB[0] - coordsVtxA[0],
                        coordsVtxB[1] - coordsVtxA[1],
                        coordsVtxB[2] - coordsVtxA[2]};
      double A1B2[3] = {coordsVtxB[3] - coordsVtxA[0],
                        coordsVtxB[4] - coordsVtxA[1],
                        coordsVtxB[5] - coordsVtxA[2]};
      double B1A1[3] = {coordsVtxA[0] - coordsVtxB[0],
                        coordsVtxA[1] - coordsVtxB[1],
                        coordsVtxA[2] - coordsVtxB[2]};
      double B1A2[3] = {coordsVtxA[3] - coordsVtxB[0],
                        coordsVtxA[4] - coordsVtxB[1],
                        coordsVtxA[5] - coordsVtxB[2]};

      double cp[3];
      PDM_CROSS_PRODUCT(cp, A1B1, vA);
      int isInBallB1 = (PDM_MODULE (cp) / vA_norm) < _charLgthVtxB[0];
      PDM_CROSS_PRODUCT(cp, A1B2, vA);
      int isInBallB2 = (PDM_MODULE (cp) / vA_norm) < _charLgthVtxB[1];
      PDM_CROSS_PRODUCT(cp, B1A1, vB);
      int isInBallA1 = (PDM_MODULE (cp) / vB_norm) < _charLgthVtxA[0];
      PDM_CROSS_PRODUCT(cp, B1A2, vB);
      int isInBallA2 = (PDM_MODULE (cp) / vB_norm) < _charLgthVtxA[1];

      int isSameLine = (isInBallB1 && isInBallB2 && isInBallA1 && isInBallA2);

      if (!isSameLine) {
        tIntersect = PDM_LINE_INTERSECT_NO;
        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                0,
                                                0);
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        isSetted = true;

      }
    }

    if (tIntersect == PDM_LINE_INTERSECT_ON_LINE) {

      double A1B1[3] = {coordsVtxB[0] - coordsVtxA[0],
                        coordsVtxB[1] - coordsVtxA[1],
                        coordsVtxB[2] - coordsVtxA[2]};
      double mA1B1 = PDM_MODULE (A1B1);

      double A1B2[3] = {coordsVtxB[3+0] - coordsVtxA[0],
                        coordsVtxB[3+1] - coordsVtxA[1],
                        coordsVtxB[3+2] - coordsVtxA[2]};
      double mA1B2 = PDM_MODULE (A1B2);

      double B1A2[3] = {coordsVtxA[3+0] - coordsVtxB[0],
                        coordsVtxA[3+1] - coordsVtxB[1],
                        coordsVtxA[3+2] - coordsVtxB[2]};
      double mB1A2 = PDM_MODULE (B1A2);

      double B2A2[3] = {coordsVtxA[3+0] - coordsVtxB[3+0],
                        coordsVtxA[3+1] - coordsVtxB[3+1],
                        coordsVtxA[3+2] - coordsVtxB[3+2]};
      double mB2A2 = PDM_MODULE (B2A2);

      int isA1B1 = (mA1B1 < _charLgthVtxA[0]) &&
        (mA1B1 < _charLgthVtxB[0]);

      int isA1B2 = (mA1B2 < _charLgthVtxA[0]) &&
        (mA1B2 < _charLgthVtxB[1]);

      int isA2B1 = (mB1A2 < _charLgthVtxA[1]) &&
        (mB1A2 < _charLgthVtxB[0]);

      int isA2B2 = (mB2A2 < _charLgthVtxA[1]) &&
        (mB2A2 < _charLgthVtxB[1]);

      int nNewPointsA = 0;
      int nNewPointsB = 0;

      double closestB2EdgeA[3];
      double tB2EdgeA;

      PDM_line_distance (coordsVtxB + 3,
                         coordsVtxA,
                         coordsVtxA + 3,
                         &tB2EdgeA,
                         closestB2EdgeA);

      double A1ClosestB2[3] = {tB2EdgeA * vA[0],
                               tB2EdgeA * vA[1],
                               tB2EdgeA * vA[2]};

      double closestA2EdgeB[3];
      double tA2EdgeB;
      PDM_line_distance (coordsVtxA + 3,
                         coordsVtxB,
                         coordsVtxB + 3,
                         &tA2EdgeB,
                         closestA2EdgeB);

      double B1ClosestA2[3] = {tA2EdgeB * vB[0],
                               tA2EdgeB * vB[1],
                               tA2EdgeB * vB[2]};

      double closestB1EdgeA[3];
      double tB1EdgeA;
      PDM_line_distance (coordsVtxB,
                         coordsVtxA,
                         coordsVtxA + 3,
                         &tB1EdgeA,
                         closestB1EdgeA);

      double A1ClosestB1[3] = {tB1EdgeA * vA[0],
                               tB1EdgeA * vA[1],
                               tB1EdgeA * vA[2]};

      double closestA1EdgeB[3];
      double tA1EdgeB;
      PDM_line_distance (coordsVtxA,
                         coordsVtxB,
                         coordsVtxB + 3,
                         &tA1EdgeB,
                         closestA1EdgeB);

      double B1ClosestA1[3] = {tA1EdgeB * vB[0],
                               tA1EdgeB * vB[1],
                               tA1EdgeB * vB[2]};

      if ((isA1B1 && isA2B2) || (isA2B1 && isA1B2)) {

        nNewPointsA = 2;
        nNewPointsB = 2;

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

        if (isA1B1 && isA2B2) {
          newInter->uA[0] = 0.;
          newInter->uA[1] = 1.;
          newInter->uB[0] = 0.;
          newInter->uB[1] = 1.;

          newInter->linkA[0] = nGVtxB[0];
          newInter->linkA[1] = nGVtxB[1];
          newInter->gNumA[0] = nGVtxA[0];
          newInter->gNumA[1] = nGVtxA[1];

          newInter->linkB[0] = nGVtxA[0];
          newInter->linkB[1] = nGVtxA[1];
          newInter->gNumB[0] = nGVtxB[0];
          newInter->gNumB[1] = nGVtxB[1];

          double _coords1[3];
          double _coords2[3];

          for (int i = 0; i < 3; i++) {
            _coords1[i] = (coordsVtxB[i] + coordsVtxA[i]) / 2;
            _coords2[i] = (coordsVtxB[3+i] + coordsVtxA[3+i]) / 2;
          }

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[i]     = _coords1[i];
            newInter->coordsB[i]     = _coords1[i];
            newInter->coordsA[3 + i] = _coords2[i];
            newInter->coordsB[3 + i] = _coords2[i];
          }
        }

        else { //  (isA2B1 && isA1B2))
          newInter->uA[0] = 0.;
          newInter->uA[1] = 1.;
          newInter->uB[0] = 0.;
          newInter->uB[1] = 1.;

          newInter->linkA[0] = nGVtxB[1];
          newInter->linkA[1] = nGVtxB[0];
          newInter->gNumA[0] = nGVtxA[0];
          newInter->gNumA[1] = nGVtxA[1];

          newInter->linkB[0] = nGVtxA[1];
          newInter->linkB[1] = nGVtxA[0];
          newInter->gNumB[0] = nGVtxB[0];
          newInter->gNumB[1] = nGVtxB[1];

          double _coords1[3];
          double _coords2[3];

          for (int i = 0; i < 3; i++) {
            _coords1[i] = (coordsVtxB[3+i] + coordsVtxA[  i]) / 2;
            _coords2[i] = (coordsVtxB[  i] + coordsVtxA[3+i]) / 2;
          }

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[i]     = _coords1[i];
            newInter->coordsB[i]     = _coords2[i];
            newInter->coordsA[3 + i] = _coords2[i];
            newInter->coordsB[3 + i] = _coords1[i];
          }
        }

        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
          newInter->uA[1] = 1 - newInter->uA[1];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
          newInter->uB[1] = 1 - newInter->uB[1];
        }

      }

      else if (isA1B1) {
        nNewPointsA = 1;
        nNewPointsB = 1;

        double u = PDM_DOT_PRODUCT (A1ClosestB2, vA) / vA_norm2;
        double v = PDM_DOT_PRODUCT (B1ClosestA2, vB) / vB_norm2;

        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->linkA[0] = nGVtxB[0];
        newInter->linkB[0] = nGVtxA[0];
        newInter->gNumA[0] = nGVtxA[0];
        newInter->gNumB[0] = nGVtxB[0];

        newInter->uA[0] = 0.;
        newInter->uB[0] = 0.;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[i] + coordsVtxA[i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i];
          newInter->coordsB[i] = _coords[i];
        }

        if (nNewPointsA == 2) {
          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[1];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[3 + i] =  closestB2EdgeA[i];
          }
        }
        if (nNewPointsB == 2) {
          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[1];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA2EdgeB[i];
          }
        }
      }

      else if (isA2B2) {

        nNewPointsA = 1;
        nNewPointsB = 1;

        double u = PDM_DOT_PRODUCT (A1ClosestB1, vA) / vA_norm2;
        double v = PDM_DOT_PRODUCT (B1ClosestA1, vB) / vB_norm2;

        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->linkA[0] = nGVtxB[1];
        newInter->linkB[0] = nGVtxA[1];
        newInter->gNumA[0] = nGVtxA[1];
        newInter->gNumB[0] = nGVtxB[1];

        newInter->uA[0] = 1.;
        newInter->uB[0] = 1.;
        if (newInter->originEdgeA == nGVtxA[1]) {
           newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[3+i] + coordsVtxA[3+i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i];
          newInter->coordsB[i] = _coords[i];
        }

        if (nNewPointsA == 2) {
          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[0];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[3 + i] =  closestB1EdgeA[i];
          }
        }

        if (nNewPointsB == 2) {
          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[0];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA1EdgeB[i];
          }
        }
      }

      else if (isA1B2) {

        nNewPointsA = 1;
        nNewPointsB = 1;

        double u = PDM_DOT_PRODUCT (A1ClosestB1, vA) / vA_norm2;
        double v = PDM_DOT_PRODUCT (B1ClosestA2, vB) / vB_norm2;

        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->linkA[0] = nGVtxB[1];
        newInter->linkB[0] = nGVtxA[0];
        newInter->gNumA[0] = nGVtxA[0];
        newInter->gNumB[0] = nGVtxB[1];

        newInter->uA[0] = 1.;
        newInter->uB[0] = 0.;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[3+i] + coordsVtxA[i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i];
          newInter->coordsB[i] = _coords[i];
        }

        if (nNewPointsA == 2) {

          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[0];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[3 + i] =  closestB1EdgeA[i];
          }
        }

        if (nNewPointsB == 2) {

          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[1];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA2EdgeB[i];
          }
        }

      }

      else if (isA2B1) {

        nNewPointsA = 1;
        nNewPointsB = 1;

        double u = PDM_DOT_PRODUCT (A1ClosestB2, vA) / vA_norm2;
        double v = PDM_DOT_PRODUCT (B1ClosestA1, vB) / vA_norm2;

        if ((u <= 1.) && (u >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v <= 1.) && (v >= 0.)) {
          nNewPointsB += 1;
        }

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);
        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        newInter->linkA[0] = nGVtxB[0];
        newInter->linkB[0] = nGVtxA[1];
        newInter->gNumA[0] = nGVtxA[1];
        newInter->gNumB[0] = nGVtxB[0];

        newInter->uA[0] = 1.;
        newInter->uB[0] = 0.;
        if (newInter->originEdgeA == nGVtxA[1]) {
          newInter->uA[0] = 1 - newInter->uA[0];
        }

        if (newInter->originEdgeB == nGVtxB[1]) {
          newInter->uB[0] = 1 - newInter->uB[0];
        }

        newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
        newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

        double _coords[3];

        for (int i = 0; i < 3; i++) {
          _coords[i] = (coordsVtxB[i] + coordsVtxA[3+i]) / 2;
        }

        for (int i = 0; i < 3; i++) {
          newInter->coordsA[i] = _coords[i];
          newInter->coordsB[i] = _coords[i];
        }

        if (nNewPointsA == 2) {
          newInter->oNewPointsA[1] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[1] = u;
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[1] = 1 - newInter->uA[1];
          }

          newInter->linkA[1] = nGVtxB[1];
          newInter->gNumA[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[3 + i] =  closestA2EdgeB[i];
          }
        }

        if (nNewPointsB == 2) {
          newInter->oNewPointsB[1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[1] = v;
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[1] = 1 - newInter->uB[1];
          }

          newInter->linkB[1] = nGVtxA[0];
          newInter->gNumB[1] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3 + i] =  closestA1EdgeB[i];
          }
        }

      }

      else {
        nNewPointsA = 0;
        nNewPointsB = 0;

        double u[2] = {PDM_DOT_PRODUCT (A1ClosestB1, vA) / vA_norm2,
                       PDM_DOT_PRODUCT (A1ClosestB2, vA) / vA_norm2};
        double v[2] = {PDM_DOT_PRODUCT (B1ClosestA1, vB) / vB_norm2,
                       PDM_DOT_PRODUCT (B1ClosestA2, vB) / vB_norm2};

        if ((u[0] <= 1) && (u[0] >= 0.)) {
          nNewPointsA += 1;
        }

        if ((u[1] <= 1.) && (u[1] >= 0.)) {
          nNewPointsA += 1;
        }

        if ((v[0] <= 1.) && (v[0] >= 0.)) {
          nNewPointsB += 1;
        }

        if ((v[1] <= 1.) && (v[1] >= 0.)) {
          nNewPointsB += 1;
        }

        newInter = _edges_intersect_res_create (nGEdgeA,
                                                nGEdgeB,
                                                nNewPointsA,
                                                nNewPointsB);

        newInter->tIntersect = tIntersect;
        newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
        newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
        newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
        newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

        nNewPointsA = 0;
        nNewPointsB = 0;

        if ((u[0] <= 1.) && (u[0] >= 0.)) {

          newInter->oNewPointsA[nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;
          newInter->uA[nNewPointsA] = u[0];
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[nNewPointsA] = 1 - newInter->uA[nNewPointsA];
          }

          newInter->linkA[nNewPointsA] = nGVtxB[0];
          newInter->gNumA[nNewPointsA] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[3*nNewPointsA + i] =  closestB1EdgeA[i];
          }
          nNewPointsA += 1;
        }

        if ((u[1] <= 1.) && (u[1] >= 0.)) {
          newInter->oNewPointsA[nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA;

          newInter->uA[nNewPointsA] = u[1];
          if (newInter->originEdgeA == nGVtxA[1]) {
            newInter->uA[nNewPointsA] = 1 - newInter->uA[nNewPointsA];
          }

          newInter->linkA[nNewPointsA] = nGVtxB[1];
          newInter->gNumA[nNewPointsA] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsA[3*nNewPointsA + i] =  closestB2EdgeA[i];
          }
          nNewPointsA += 1;
        }

        if ((v[0] <= 1.) && (v[0] >= 0.)) {
          newInter->oNewPointsB[nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[nNewPointsB] = v[0];
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[nNewPointsB] = 1 - newInter->uB[nNewPointsB];
          }

          newInter->linkB[nNewPointsB] = nGVtxA[0];
          newInter->gNumB[nNewPointsB] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3*nNewPointsB + i] =  closestA1EdgeB[i];
          }
          nNewPointsB += 1;
        }

        if ((v[1] <= 1.) && (v[1] >= 0.)) {
          newInter->oNewPointsB[nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB;

          newInter->uB[nNewPointsB] = v[1];
          if (newInter->originEdgeB == nGVtxB[1]) {
            newInter->uB[nNewPointsB] = 1 - newInter->uB[nNewPointsB];
          }

          newInter->linkB[nNewPointsB] = nGVtxA[1];
          newInter->gNumB[nNewPointsB] = 0;

          for (int i = 0; i < 3; i++) {
            newInter->coordsB[3*nNewPointsB + i] =  closestA2EdgeB[i];
          }
          nNewPointsB += 1;
        }
      }

      isSetted = true;
    }
  }

  /*
   * Storage intersection
   */

  if (!isSetted && (tIntersect == PDM_LINE_INTERSECT_YES)) {

    int nNewPointsA = 1;
    int nNewPointsB = 1;

    newInter = _edges_intersect_res_create (nGEdgeA,
                                            nGEdgeB,
                                            nNewPointsA,
                                            nNewPointsB);
    newInter->tIntersect = tIntersect;
    newInter->originEdgeA = PDM_MIN (nGVtxA[0], nGVtxA[1]);
    newInter->originEdgeB = PDM_MIN (nGVtxB[0], nGVtxB[1]);
    newInter->endEdgeA = PDM_MAX (nGVtxA[0], nGVtxA[1]);
    newInter->endEdgeB = PDM_MAX (nGVtxB[0], nGVtxB[1]);

    newInter->oNewPointsA[0] = PDM_EDGES_INTERSECT_POINT_NEW;

    newInter->oNewPointsB[0] = PDM_EDGES_INTERSECT_POINT_NEW;

    newInter->linkA[0] = -1;
    newInter->linkB[0] = -1;
    newInter->gNumA[0] = 0;
    newInter->gNumB[0] = 0;

    newInter->uA[0] = u1;
    newInter->uB[0] = v1;

    if (newInter->originEdgeA == nGVtxA[1]) {
      newInter->uA[0] = 1 - newInter->uA[0];
    }

    if (newInter->originEdgeB == nGVtxB[1]) {
      newInter->uB[0] = 1 - newInter->uB[0];
    }

    for (int i = 0; i < 3; i++) {
      newInter->coordsA[i] =  coordsVtxA[i] + u1 * vA[i];
    }

    for (int i = 0; i < 3; i++) {
      newInter->coordsB[i] =  coordsVtxB[i] + v1 * vB[i];
    }

    isSetted = true;
  }

  if( newInter != NULL) {
	  PDM_hash_tab_data_add (ht, (void *) &key, newInter);
	  PDM_hash_tab_data_add (htA, (void *) &keyA, newInter);
	  PDM_hash_tab_data_add (htB, (void *) &keyB, newInter);
  }
  if (vb) {
    if (newInter != NULL)
      PDM_edges_intersect_res_dump((PDM_edges_intersect_res_t *)newInter);
    PDM_printf ("==== PDM_edges_intersect_add ==== terminated ====\n");
  }

  return (PDM_edges_intersect_res_t *) newInter;
}


/**
 *
 * \brief Free \ref PDM_edges_intersect_t object
 *
 * \param [in]   ei               Current edges intersection pointer
 *
 * \return     NULL
 */

PDM_edges_intersect_t *
PDM_edges_intersect_free
(
PDM_edges_intersect_t *ei
)
{
  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;

  int *keyMax;

  keyMax = (int *) PDM_hash_tab_keyMax_get (_ei->ht);
  for (int i = 0; i < *keyMax; i++) {
    int nData = PDM_hash_tab_n_data_get (_ei->ht, &i);
    _edges_intersect_res_t **eir =
            (_edges_intersect_res_t **) PDM_hash_tab_data_get (_ei->ht, &i);

    for (int j = 0; j < nData; j++) {
      _edges_intersect_res_free (eir[j]);
    }
  }

  PDM_hash_tab_free (_ei->ht);
  PDM_hash_tab_free (_ei->htA);
  PDM_hash_tab_free (_ei->htB);

  free (_ei);

  return NULL;
}


/**
 *
 * \brief Get data intersection
 *
 * \param [in]  eir             Current edges intersection result pointer
 * \param [in]  mesh            Origin mesh \ref PDM_edges_intersect_MESHA or
 *                                          \ref PDM_edges_intersect_MESHB
 * \param [out] nGEdge          Global number of meshA edge
 * \param [out] originEdge      Global number of edge origin
 * \param [out] tIntersect      Intersection type
 * \param [out] nNewPoints      Number of intersection points
 * \param [out] oNewPoints      Origin of intersection points
 * \param [out] link            Linked vertex in linked mesh
 * \param [out] gNum            Global number in overlay mesh
 * \param [out] coords          Coordinates of intersection point
 * \param [out] u               Parameter of the intersections in edges
 *                              parametric coordinates
 *
 */

void
PDM_edges_intersect_res_data_get
(
PDM_edges_intersect_res_t   *eir,
PDM_edges_intersect_mesh_t   mesh,
PDM_g_num_t                  *nGEdge,
PDM_g_num_t                  *originEdge,
PDM_g_num_t                  *endEdge,
PDM_line_intersect_t        *tIntersect,
int                         *nNewPoints,
PDM_edges_intersect_point_t **oNewPoints,
PDM_g_num_t                  **link,
PDM_g_num_t                  **gNum,
double                      **coords,
double                      **u
)
{

  if (eir == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_edges_intersect_res_data_get : The Input PDM_edges_intersect_res_t is NULL.\n");
    abort();
}

  _edges_intersect_res_t *_eir = (_edges_intersect_res_t *) eir;

  if (mesh == PDM_EDGES_INTERSECT_MESHA) {
    *nGEdge          = _eir->nGEdgeA;
    *originEdge      = _eir->originEdgeA;
    *endEdge         = _eir->endEdgeA;
    *tIntersect      = _eir->tIntersect;
    *nNewPoints      = _eir->nNewPointsA;
    *oNewPoints      = _eir->oNewPointsA;
    *link            = _eir->linkA;
    *gNum            = _eir->gNumA;
    *coords          = _eir->coordsA;
    *u               = _eir->uA;
  }
  if (mesh == PDM_EDGES_INTERSECT_MESHB) {
    *nGEdge          = _eir->nGEdgeB;
    *originEdge      = _eir->originEdgeB;
    *endEdge         = _eir->endEdgeB;
    *tIntersect      = _eir->tIntersect;
    *nNewPoints      = _eir->nNewPointsB;
    *oNewPoints      = _eir->oNewPointsB;
    *link            = _eir->linkB;
    *gNum            = _eir->gNumB;
    *coords          = _eir->coordsB;
    *u               = _eir->uB;
  }
}


/**
 *
 * \brief Perform edges intersection from two polygons
 *
 * \param [in]    intersect            Edges intersection management
 * \param [in]    n_vtxA                Number of polygon A vertices
 * \param [in]    faceToEdgeA          Polygon A face to edge connectivity
 * \param [in]    faceToVtxA           Polygon A face to vertex connectivity
 * \param [in]    face_vtxCooA          Polygon A vertex coordinates
 * \param [in]    face_vtxEpsA          Polygon A vertex characteristic length
 * \param [in]    gNumB                Polygon B global number
 * \param [in]    n_vtxB                Number of polygon B vertices
 * \param [in]    faceToEdgeB          Polygon B face to edge connectivity
 * \param [in]    faceToVtxB           Polygon B face to vertex connectivity
 * \param [in]    face_vtxCooB          Polygon B vertex coordinates
 * \param [in]    face_vtxEpsB          Polygon B vertex characteristic length
 *
 */

void
PDM_edges_intersect_poly_add
(
PDM_edges_intersect_t  *ei,
const int               n_vtxA,
PDM_g_num_t            *faceToEdgeA,
PDM_g_num_t            *faceToVtxA,
double                 *face_vtxCooA,
double                 *face_vtxEpsA,
const int               n_vtxB,
PDM_g_num_t            *faceToEdgeB,
PDM_g_num_t            *faceToVtxB,
double                 *face_vtxCooB,
double                 *face_vtxEpsB
)
{
  int vb = 0;
  if (vb) {
    PDM_printf ("==== PDM_edges_intersect_poly_add ==== \n");
  }

  PDM_g_num_t *_faceToEdgeA = faceToEdgeA;
  PDM_g_num_t *_faceToVtxA  = faceToVtxA;
  double     *_face_vtxCooA = face_vtxCooA;
  double     *_face_vtxEpsA = face_vtxEpsA;

  PDM_g_num_t *_faceToEdgeB = faceToEdgeB;
  PDM_g_num_t *_faceToVtxB  = faceToVtxB;
  double     *_face_vtxCooB = face_vtxCooB;
  double     *_face_vtxEpsB = face_vtxEpsB;

  /*
   * Compute Normal
   */

  double nA[3];
  double baryA[3];
  PDM_plane_normal (n_vtxA, face_vtxCooA, nA);
  PDM_plane_barycenter (n_vtxA, face_vtxCooA, baryA);

  double nB[3];
  PDM_plane_normal (n_vtxB, face_vtxCooB, nB);

  double dot = PDM_DOT_PRODUCT (nA, nB);

  bool revert = false;

   if (dot < 0) {
     revert = true;
   }

  /*
   * Reorient if necessary
   */

  if (revert) {

    _faceToEdgeB = malloc (sizeof(PDM_g_num_t) * n_vtxB);
    _faceToVtxB  = malloc (sizeof(PDM_g_num_t) * n_vtxB);
    _face_vtxCooB = malloc (sizeof(double) * 3 * n_vtxB);
    _face_vtxEpsB = malloc (sizeof(double) * n_vtxB);

    int j = n_vtxB - 1;
    for (int i = 0; i < n_vtxB; i++) {
      _faceToEdgeB[i] = faceToEdgeB[j];
      _faceToEdgeB[i] = -faceToEdgeB[_modulo((j-1),n_vtxB)];
      _faceToVtxB[i] = faceToVtxB[j];
      for (int k = 0; k < 3; k++) {
        _face_vtxCooB[3*i+k] = face_vtxCooB[3*j+k];
      }
      _face_vtxEpsB[i] = face_vtxEpsB[j];
      j += -1;
    }

  }

  _face_vtxCooA = malloc (sizeof(double) * 3 * n_vtxA);
  for (int i = 0; i < n_vtxA; i++) {
    PDM_plane_projection (face_vtxCooA + 3 * i, baryA, nA, _face_vtxCooA + 3 * i);
  }
  //PDM_plane_normal (n_vtxA, _face_vtxCooA, nA);

  if (revert) {
    for (int i = 0; i < n_vtxB; i++) {
      PDM_plane_projection (_face_vtxCooB + 3 * i, baryA, nA, _face_vtxCooB + 3 * i);
    }
  }
  else {
    _face_vtxCooB = malloc (sizeof(double) * 3 * n_vtxB);
    for (int i = 0; i < n_vtxB; i++) {
      PDM_plane_projection (face_vtxCooB + 3 * i, baryA, nA, _face_vtxCooB + 3 * i);
    }
  }


  /*
   *   Compute Edges Intersection :
   *   - First step : Remove case with vertex located on 2 two different edges
   *   - Second step : Compute other intersection
   */

  //
  // FIXME: faire un test pour verifier que le resultat de 2 intersections d'une
  // arete de A de B donnent toujours des resultats differents
  // Faire l'inverse : Detection  des aretes a pb .
  // Echange MPI des aretes a pb broad cast apres un premier clipping
  // Mise a jour des clipping concernes
  //

  _edges_intersect_res_t **vtxAOnEdgeBEir = malloc(sizeof(_edges_intersect_res_t *) * n_vtxA);
  for (int i = 0; i < n_vtxA; i++) {
    vtxAOnEdgeBEir[i] = NULL;
  }

  _edges_intersect_res_t **vtxBOnEdgeAEir = malloc(sizeof(_edges_intersect_res_t *) * n_vtxB);
  for (int i = 0; i < n_vtxB; i++) {
    vtxBOnEdgeAEir[i] = NULL;
  }

  for (int i = 0; i < n_vtxA; i++) {

    int inext = (i + 1) % n_vtxA;

    PDM_g_num_t nGVtxA[2]   = {_faceToVtxA[i], _faceToVtxA[inext]};
    double charLgthVtxA[2] = {_face_vtxEpsA[i], _face_vtxEpsA[inext]};
    double coordsVtxA[6]   = {_face_vtxCooA[3*i],
                              _face_vtxCooA[3*i+1],
                              _face_vtxCooA[3*i+2],
                              _face_vtxCooA[3*inext],
                              _face_vtxCooA[3*inext+1],
                              _face_vtxCooA[3*inext+2]};

    for (int j = 0; j < n_vtxB; j++) {
      int jnext = (j + 1) % n_vtxB;
      PDM_g_num_t nGVtxB[2]   = {_faceToVtxB[j],
                                 _faceToVtxB[jnext]};
      double charLgthVtxB[2] = {_face_vtxEpsB[j], _face_vtxEpsB[jnext]};
      double coordsVtxB[6]   = {_face_vtxCooB[3*j],
                                _face_vtxCooB[3*j+1],
                                _face_vtxCooB[3*j+2],
                                _face_vtxCooB[3*jnext],
                                _face_vtxCooB[3*jnext+1],
                                _face_vtxCooB[3*jnext+2]};

      /*
       * Perform intersection
       */

      if (vb) {
        PDM_printf ("_faceToEdgeA[i]:%d, nGVtxA:%d-%d, charLgthVtxA:%12.5e-%12.5e, "
                    "coordsVtxA:%12.5e-%12.5e-%12.5e-%12.5e-%12.5e-%12.5e\n",
                    _faceToEdgeA[i], nGVtxA[0],
                    nGVtxA[1], charLgthVtxA[0], charLgthVtxA[1], coordsVtxA[0],
                    coordsVtxA[1], coordsVtxA[2], coordsVtxA[3], coordsVtxA[4], coordsVtxA[5]);
        PDM_printf ("_faceToEdgeB[j]:%d, nGVtxB:%d-%d, charLgthVtxB:%12.5e-%12.5e, "
                    "coordsVtxB:%12.5e-%12.5e-%12.5e-%12.5e-%12.5e-%12.5e \n",
                    _faceToEdgeB[j], nGVtxB[0],
                    nGVtxB[1], charLgthVtxB[0], charLgthVtxB[1], coordsVtxB[0],
                    coordsVtxB[1], coordsVtxB[2], coordsVtxB[3], coordsVtxB[4], coordsVtxB[5]);
      }

      PDM_edges_intersect_res_t *eir =
        PDM_edges_intersect_add (ei,
                                 PDM_ABS(_faceToEdgeA[i]),
                                 nGVtxA,
                                 charLgthVtxA,
                                 coordsVtxA,
                                 PDM_ABS(_faceToEdgeB[j]),
                                 nGVtxB,
                                 charLgthVtxB,
                                 coordsVtxB);


      _edges_intersect_res_t *_eir = (_edges_intersect_res_t *) eir;

      /*
       * Get intersection properties if intersection!!!
       */

      if (eir != NULL) {

        PDM_line_intersect_t         tIntersect;

        PDM_g_num_t                   nGEdgeA;
        PDM_g_num_t                   originEdgeA;
        PDM_g_num_t                   endEdgeA;
        int                          nNewPointsA;
        PDM_edges_intersect_point_t *oNewPointsA;
        PDM_g_num_t                  *linkA;
        PDM_g_num_t                  *gNumA;
        double                      *coordsA;
        double                      *uA;

        PDM_edges_intersect_res_data_get (eir,
                                          PDM_EDGES_INTERSECT_MESHA,
                                          &nGEdgeA,
                                          &originEdgeA,
                                          &endEdgeA,
                                          &tIntersect,
                                          &nNewPointsA,
                                          &oNewPointsA,
                                          &linkA,
                                          &gNumA,
                                          &coordsA,
                                          &uA);

        PDM_g_num_t                   nGEdgeB;
        PDM_g_num_t                   originEdgeB;
        PDM_g_num_t                   endEdgeB;
        int                          nNewPointsB;
        PDM_edges_intersect_point_t *oNewPointsB;
        PDM_g_num_t                  *linkB;
        PDM_g_num_t                  *gNumB;
        double                      *coordsB;
        double                      *uB;

        PDM_edges_intersect_res_data_get (eir,
                                          PDM_EDGES_INTERSECT_MESHB,
                                          &nGEdgeB,
                                          &originEdgeB,
                                          &endEdgeB,
                                          &tIntersect,
                                          &nNewPointsB,
                                          &oNewPointsB,
                                          &linkB,
                                          &gNumB,
                                          &coordsB,
                                          &uB);

        /*
         * Remove inconsistencies :
         * Check if a vertex is not on two different edges
         *   - case 1 : B vertex on two 2 different A edges
         *   - case 2 : A vertex on two 2 different B edges
         */

        for (int k = 0; k < nNewPointsA; k++) {

          if (oNewPointsA[k] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
            int ind = j;
            if (linkA[k] != nGVtxB[0]) {
              ind = jnext;
              assert(linkA[k] == nGVtxB[1]);
            }

            /*
             * If B vertex is already on an other A edge :
             * Update intersections to merge vertex case
             */

            if (vtxBOnEdgeAEir[ind] != NULL) {

              PDM_edges_intersect_res_t *preEir = (PDM_edges_intersect_res_t *) vtxBOnEdgeAEir[ind];
              _edges_intersect_res_t *_preEir = (_edges_intersect_res_t *) preEir;

              PDM_g_num_t                  preNGEdgeA;
              PDM_g_num_t                  preExtEdgeA[2];
              PDM_line_intersect_t         preTIntersectA;
              int                          preNNewPointsA;
              PDM_edges_intersect_point_t *preONewPointsA;
              PDM_g_num_t                 *preLinkA;
              PDM_g_num_t                 *preGNumA;
              double                      *preCoordsA;
              double                      *preUA;

              PDM_edges_intersect_res_data_get (preEir,
                                                PDM_EDGES_INTERSECT_MESHA,
                                                &preNGEdgeA,
                                                &(preExtEdgeA[0]),
                                                &(preExtEdgeA[1]),
                                                &preTIntersectA,
                                                &preNNewPointsA,
                                                &preONewPointsA,
                                                &preLinkA,
                                                &preGNumA,
                                                &preCoordsA,
                                                &preUA);

              PDM_g_num_t                  preNGEdgeB;
              PDM_g_num_t                  preExtEdgeB[2];
              PDM_line_intersect_t         preTIntersectB;
              int                          preNNewPointsB;
              PDM_edges_intersect_point_t *preONewPointsB;
              PDM_g_num_t                 *preLinkB;
              PDM_g_num_t                 *preGNumB;
              double                      *preCoordsB;
              double                      *preUB;

              PDM_edges_intersect_res_data_get (preEir,
                                                PDM_EDGES_INTERSECT_MESHB,
                                                &preNGEdgeB,
                                                &(preExtEdgeB[0]),
                                                &(preExtEdgeB[1]),
                                                &preTIntersectB,
                                                &preNNewPointsB,
                                                &preONewPointsB,
                                                &preLinkB,
                                                &preGNumB,
                                                &preCoordsB,
                                                &preUB);

              if (nGEdgeA != preNGEdgeA) {
                printf("Warning : B vtx already on Edge A new old : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" \n" , linkA[k], nGEdgeA, preNGEdgeA);
                printf("TODO : Revoir le cas et creer deux points d'intersection"
                       " au lieu de fusionner avec le sommet commun\n.");
                //FIXME : Au lieu de fusionner avec le sommet commun, il faut creer deux vrais intersections
                //sans prendre en compte le eps relatif
                /*
                 * Look for common vertex
                 */

                int common_isom = -1;

                for (int k1 = 0; k1 < 2; k1++) {
                  for (int k2 = 0; k2 < 2; k2++) {
                    if (nGVtxA[k1] == preExtEdgeA[k2]) {
                      common_isom = k1;
                      break;
                    }
                  }
                }

                if (common_isom == -1) {
                  PDM_error(__FILE__, __LINE__, 0, "Probleme simplication sommet proche de 2 aretes :\n"
                            "les 2 aretes n'ont pas de sommet commun\n");
                  abort();
                }

                double coords_new[3];
                for (int k1 = 0; k1 < 3; k1++) {
                  coords_new[k1] = _face_vtxCooB[3*ind+k1];
                }

                /*
                 * Update preEir and eir to merge vertex case
                 *   - From A : switch to PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB
                 *   - From B : Add a new modified A vertex
                 *   - Move vertices to half-distance.
                 */

                /* From A */

                for (int k1 = 0; k1 < _preEir->nNewPointsA; k1++) {
                  if ((_preEir->oNewPointsA[k1] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) &&
                      (_preEir->linkA[k1] == linkA[k])) {

                    _preEir->oNewPointsA[k1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

                    _preEir->gNumA[k1]     = nGVtxA[common_isom];

                    if (_preEir->originEdgeA == nGVtxA[common_isom]) {
                      _preEir->uA[k1] = 0;
                    }
                    else {
                      _preEir->uA[k1] = 1.;
                    }
                    // Le point d'intersection est le sommet B

                    for (int k2 = 0; k2 < 3; k2++) {
                      _preEir->coordsA[3*k1+k2] = coords_new[k2];
                    }
                  }
                }

                _eir->oNewPointsA[k] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

                _eir->gNumA[k]     = nGVtxA[common_isom];

                if (_eir->originEdgeA == nGVtxA[common_isom]) {
                  _eir->uA[k] = 0;
                }
                else {
                  _eir->uA[k] = 1.;
                }
                for (int k2 = 0; k2 < 3; k2++) {
                _eir->coordsA[3*k+k2] = coords_new[k2];
                }

                /* From B */

                _preEir->oNewPointsB =
                  realloc(_preEir->oNewPointsB,
                          (_preEir->nNewPointsB + 1) * sizeof(PDM_edges_intersect_point_t));
                oNewPointsB = _preEir->oNewPointsB;

                _preEir->linkB = realloc(_preEir->linkB, sizeof(PDM_g_num_t) * (_preEir->nNewPointsB + 1));
                linkB = _preEir->linkB;

                _preEir->gNumB = realloc(_preEir->gNumB, sizeof(PDM_g_num_t) * (_preEir->nNewPointsB + 1));
                gNumB = _preEir->gNumB;

                _preEir->uB  = realloc(_preEir->uB, sizeof(double) * (_preEir->nNewPointsB + 1));
                uB = _preEir->uB;

                _preEir->coordsB  =
                  realloc(_preEir->coordsB, sizeof(double) * 3 * (_preEir->nNewPointsB + 1));
                coordsB = _preEir->coordsB;

                _preEir->oNewPointsB[_preEir->nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                _preEir->linkB[_preEir->nNewPointsB] = nGVtxA[common_isom];
                _preEir->gNumB[_preEir->nNewPointsB] = linkA[k];

                if (_preEir->originEdgeB == linkA[k]) {
                  _preEir->uB[_preEir->nNewPointsB] = 0;
                }
                else {
                  _preEir->uB[_preEir->nNewPointsB] = 1.;
                }
                for (int k2 = 0; k2 < 3; k2++) {
                  _preEir->coordsB[3*_preEir->nNewPointsB+k2] = coords_new[k2];
                }

                _preEir->nNewPointsB += 1 ;

                _eir->oNewPointsB = realloc(_eir->oNewPointsB,
                                            (_eir->nNewPointsB + 1) * sizeof(PDM_edges_intersect_point_t));
                oNewPointsB = _eir->oNewPointsB;

                _eir->linkB = realloc(_eir->linkB, sizeof(PDM_g_num_t) * (_eir->nNewPointsB + 1));
                linkB = _eir->linkB;

                _eir->gNumB = realloc(_eir->gNumB, sizeof(PDM_g_num_t) * (_eir->nNewPointsB + 1));
                gNumB = _eir->gNumB;

                _eir->uB  = realloc(_eir->uB, sizeof(double) * (_eir->nNewPointsB + 1));
                uB = _eir->uB;

                _eir->coordsB  = realloc(_eir->coordsB, sizeof(double) * 3 * (_eir->nNewPointsB + 1));
                coordsB = _eir->coordsB;

                _eir->oNewPointsB[_eir->nNewPointsB] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                _eir->linkB[_eir->nNewPointsB] = nGVtxA[common_isom];
                _eir->gNumB[_eir->nNewPointsB] = linkA[k];

                if (_eir->originEdgeB == linkA[k]) {
                  _eir->uB[_eir->nNewPointsB] = 0;
                }
                else {
                  _eir->uB[_eir->nNewPointsB] = 1.;
                }
                for (int k2 = 0; k2 < 3; k2++) {
                  _eir->coordsB[3*_eir->nNewPointsB+k2] = coords_new[k2];
                }

                _eir->nNewPointsB += 1 ;

              }
            }

            /*
             * If point B is not on an other edge : tag it
             */

            else {
              vtxBOnEdgeAEir[ind] = _eir;
            }
          }
        }

        for (int k = 0; k < nNewPointsB; k++) {

          if (oNewPointsB[k] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {

            int ind = i;
            if (linkB[k] != nGVtxA[0]) {
              ind = inext;
              assert(linkB[k] == nGVtxA[1]);
            }

            /*
             * If A vertex is already on an other B edge :
             * Update intersections to merge vertex case
             */

            if (vtxAOnEdgeBEir[ind] != NULL) {

              PDM_edges_intersect_res_t *preEir = (PDM_edges_intersect_res_t *) vtxAOnEdgeBEir[ind];
              _edges_intersect_res_t *_preEir = (_edges_intersect_res_t *) preEir;

              PDM_g_num_t                  preNGEdgeA;
              PDM_g_num_t                  preExtEdgeA[2];
              PDM_line_intersect_t         preTIntersectA;
              int                          preNNewPointsA;
              PDM_edges_intersect_point_t *preONewPointsA;
              PDM_g_num_t                 *preLinkA;
              PDM_g_num_t                 *preGNumA;
              double                      *preCoordsA;
              double                      *preUA;

              PDM_edges_intersect_res_data_get (preEir,
                                                PDM_EDGES_INTERSECT_MESHA,
                                                &preNGEdgeA,
                                                &(preExtEdgeA[0]),
                                                &(preExtEdgeA[1]),
                                                &preTIntersectA,
                                                &preNNewPointsA,
                                                &preONewPointsA,
                                                &preLinkA,
                                                &preGNumA,
                                                &preCoordsA,
                                                &preUA);

              PDM_g_num_t                  preNGEdgeB;
              PDM_g_num_t                  preExtEdgeB[2];
              PDM_line_intersect_t         preTIntersectB;
              int                          preNNewPointsB;
              PDM_edges_intersect_point_t *preONewPointsB;
              PDM_g_num_t                 *preLinkB;
              PDM_g_num_t                 *preGNumB;
              double                      *preCoordsB;
              double                      *preUB;

              PDM_edges_intersect_res_data_get (preEir,
                                                PDM_EDGES_INTERSECT_MESHB,
                                                &preNGEdgeB,
                                                &(preExtEdgeB[0]),
                                                &(preExtEdgeB[1]),
                                                &preTIntersectB,
                                                &preNNewPointsB,
                                                &preONewPointsB,
                                                &preLinkB,
                                                &preGNumB,
                                                &preCoordsB,
                                                &preUB);

              if (nGEdgeB != preNGEdgeB) {
                printf("Warning : A vtx already on Edge B new old : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n" , linkB[k], nGEdgeB, preNGEdgeB);
                printf("TODO : Revoir le cas et creer deux points d'intersection"
                       " au lieu de fusionner avec le sommet commun\n.");

                /*
                 * Look for common vertex
                 */

                int common_isom = -1;

                for (int k1 = 0; k1 < 2; k1++) {
                  for (int k2 = 0; k2 < 2; k2++) {
                    if (nGVtxB[k1] == preExtEdgeB[k2]) {
                      common_isom = k1;
                      break;
                    }
                  }
                }

                if (common_isom == -1) {
                  PDM_error(__FILE__, __LINE__, 0, "Probleme simplication sommet proche de 2 aretes :\n"
                            "les 2 aretes n'ont pas de sommet commun\n");
                  abort();
                }

                double coords_new[3];
                for (int k1 = 0; k1 < 3; k1++) {
                  coords_new[k1] = _face_vtxCooA[3*ind+k1];
                }

                /*
                 * Update preEir and eir to merge vertex case
                 *   - From B : switch to PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB
                 *   - From A : Add a new modified A vertex
                 *   - Move vertices to half-distance.
                 */

                /* From B */

                for (int k1 = 0; k1 < _preEir->nNewPointsB; k1++) {
                  if ((_preEir->oNewPointsB[k1] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) &&
                    (_preEir->linkB[k1] == linkB[k])) {

                    _preEir->oNewPointsB[k1] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

                    _preEir->gNumB[k1]     = nGVtxB[common_isom];

                    if (_preEir->originEdgeB == nGVtxB[common_isom]) {
                      _preEir->uB[k1] = 0;
                    }
                    else {
                      _preEir->uB[k1] = 1.;
                    }

                    // Le point d'intersection est le sommet A

                    for (int k2 = 0; k2 < 3; k2++) {
                      _preEir->coordsB[3*k1+k2] = coords_new[k2];
                    }
                  }
                }

                _eir->oNewPointsB[k] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;

                _eir->gNumB[k]     = nGVtxB[common_isom];

                if (_eir->originEdgeB == nGVtxB[common_isom]) {
                  _eir->uB[k] = 0;
                }
                else {
                  _eir->uB[k] = 1.;
                }
                for (int k2 = 0; k2 < 3; k2++) {
                  _eir->coordsB[3*k+k2] = coords_new[k2];
                }

                /* From A */

                _preEir->oNewPointsA =
                  realloc(_preEir->oNewPointsA,
                          (_preEir->nNewPointsA + 1) * sizeof(PDM_edges_intersect_point_t));
                oNewPointsA = _preEir->oNewPointsA;

                _preEir->linkA = realloc(_preEir->linkA, sizeof(PDM_g_num_t) * (_preEir->nNewPointsA + 1));
                linkA = _preEir->linkA;

                _preEir->gNumA = realloc(_preEir->gNumA, sizeof(PDM_g_num_t) * (_preEir->nNewPointsA + 1));
                gNumA = _preEir->gNumA;

                _preEir->uA  = realloc(_preEir->uA, sizeof(double) * (_preEir->nNewPointsA + 1));
                uA = _preEir->uA;

                _preEir->coordsA  =
                  realloc(_preEir->coordsA, sizeof(double) * 3 * (_preEir->nNewPointsA + 1));
                coordsA = _preEir->coordsA;

                _preEir->oNewPointsA[_preEir->nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                _preEir->linkA[_preEir->nNewPointsA] = nGVtxB[common_isom];
                _preEir->gNumA[_preEir->nNewPointsA] = linkB[k];

                if (_preEir->originEdgeA == linkB[k]) {
                  _preEir->uA[_preEir->nNewPointsA] = 0;
                }
                else {
                  _preEir->uA[_preEir->nNewPointsA] = 1.;
                }
                for (int k2 = 0; k2 < 3; k2++) {
                _preEir->coordsA[3*_preEir->nNewPointsA+k2] = coords_new[k2];
                }

                _preEir->nNewPointsA += 1 ;

                _eir->oNewPointsA = realloc(_eir->oNewPointsA,
                                            (_eir->nNewPointsA + 1) * sizeof(PDM_edges_intersect_point_t));
                oNewPointsA = _eir->oNewPointsA;

                _eir->linkA = realloc(_eir->linkA, sizeof(PDM_g_num_t) * (_eir->nNewPointsA + 1));
                linkA = _eir->linkA;

                _eir->gNumA = realloc(_eir->gNumA, sizeof(PDM_g_num_t) * (_eir->nNewPointsA + 1));
                gNumA = _eir->gNumA;

                _eir->uA  = realloc(_eir->uA, sizeof(double) * (_eir->nNewPointsA + 1));
                uA = _eir->uA;

                _eir->coordsA  = realloc(_eir->coordsA, sizeof(double) * 3 * (_eir->nNewPointsA + 1));
                coordsA = _eir->coordsA;

                _eir->oNewPointsA[_eir->nNewPointsA] = PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB;
                _eir->linkA[_eir->nNewPointsA] = nGVtxB[common_isom];
                _eir->gNumA[_eir->nNewPointsA] = linkB[k];

                if (_eir->originEdgeA == linkB[k]) {
                  _eir->uA[_eir->nNewPointsA] = 0;
                }
                else {
                  _eir->uA[_eir->nNewPointsA] = 1.;
                }
                for (int k2 = 0; k2 < 3; k2++) {
                  _eir->coordsA[3*_eir->nNewPointsA+k2] = coords_new[k2];
                }

                _eir->nNewPointsA += 1 ;

              }
            }

            /*
             * If point B is not on an other edge : tag it
             */

            else {
              vtxAOnEdgeBEir[ind] = _eir;
            }
          }
        }
      }
    }
  }


  /*
   * Free local memory
   */


  free (vtxAOnEdgeBEir);

  free (vtxBOnEdgeAEir);

  if (_faceToEdgeA != faceToEdgeA) {
    free (_faceToEdgeA);
  }

  if (_faceToVtxA  != faceToVtxA) { /* EQU == */
    free (_faceToVtxA);
  }

  if (_face_vtxCooA != face_vtxCooA) { /* EQU == */
    free (_face_vtxCooA);
  }

  if (_face_vtxEpsA != face_vtxEpsA) { /* EQU == */
    free (_face_vtxEpsA);
  }

  if (_faceToEdgeB != faceToEdgeB) {
    free (_faceToEdgeB);
  }

  if (_faceToVtxB != faceToVtxB) { /* EQU == */
    free (_faceToVtxB);
  }

  if (_face_vtxCooB != face_vtxCooB) { /* EQU == */
    free (_face_vtxCooB);
  }

  if (_face_vtxEpsB != face_vtxEpsB) { /* EQU == */
    free (_face_vtxEpsB);
  }
  if (vb) {
    PDM_printf ("==== PDM_edges_intersect_poly_add ==== terminated ====\n");
  }

}


/**
 *
 * \brief Remove inconsistencies between processes
 *
 * \param [in]   ei           Current edges intersection pointer
 * \param [in]   nAbsVtxA     Absolute number of vertices in initial A mesh
 * \param [in]   nAbsVtxB     Absolute number of vertices in initial B mesh
 * \param [out]  nAbsNewVtxA  Absolute number of vertices in A mesh after intersections
 * \param [out]  nAbsNewVtxB  Absolute number of vertices in B mesh after intersections
 *
 */

void
PDM_edges_intersect_synchronize
(
PDM_edges_intersect_t       *ei,
PDM_g_num_t             nAbsVtxA,
PDM_g_num_t             nAbsVtxB,
PDM_g_num_t            *nAbsNewVtxA,
PDM_g_num_t            *nAbsNewVtxB
)
{
  int vb = 0;
  if (vb) {
    PDM_printf ("==== PDM_edges_intersect_synchronize ====\n");
  }

  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;

  int sMSGComm;
  PDM_MPI_Comm_size (_ei->comm, &sMSGComm);
  int lastRank = sMSGComm - 1;

  int i_rank;
  PDM_MPI_Comm_rank (_ei->comm, &i_rank);

  /*
   * Part to block hash table
   */

  PDM_hash_tab_t *ht = _ei->ht;

  int keyMax  = * ((int *) PDM_hash_tab_keyMax_get (ht));
  // int *nDataKey = malloc (sizeof(int) * keyMax);

  /* for (int key = 0; key < keyMax; key++) { */
  /*   nDataKey[key] = 0; */
  /* } */

  int n_procData = 0;
  for (int key = 0; key < keyMax; key++) {
    int nData = PDM_hash_tab_n_data_get (ht, (void *) &key);
    //    nDataKey[key] = nData;
    n_procData += nData;
  }

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 1\n"); */
  /*   fflush(stdout); */
  /* } */


  PDM_g_num_t *keys        = malloc (sizeof(PDM_g_num_t) * n_procData);

  int        *tIntersects = malloc (sizeof(int) * n_procData);

  PDM_g_num_t *gNumEdgeA   = malloc (sizeof(PDM_g_num_t) * n_procData);
  PDM_g_num_t *gNumEdgeB   = malloc (sizeof(PDM_g_num_t) * n_procData);
  int        *nNewPointsA = malloc (sizeof(int) * n_procData);
  PDM_edges_intersect_point_t *oNewPointsA =
    malloc (sizeof(PDM_edges_intersect_point_t) * 2 * n_procData);

  PDM_g_num_t *connectPointA = malloc (sizeof(PDM_g_num_t) * 2 * n_procData);
  PDM_g_num_t *gNumA = malloc (sizeof(PDM_g_num_t) * 2 * n_procData);
  double *uPointA = malloc (sizeof(double) * 2 * n_procData);
  double *coordsPointA = malloc (sizeof(double) * 6 * n_procData);

  int        *nNewPointsB = malloc (sizeof(int) * n_procData);
  PDM_edges_intersect_point_t *oNewPointsB =
    malloc (sizeof(PDM_edges_intersect_point_t) * 2 * n_procData);

  PDM_g_num_t *connectPointB = malloc (sizeof(PDM_g_num_t) * 2 * n_procData);
  PDM_g_num_t *gNumB = malloc (sizeof(PDM_g_num_t) * 2 * n_procData);
  double *uPointB = malloc (sizeof(double) * 2 * n_procData);
  double *coordsPointB = malloc (sizeof(double) * 6 * n_procData);

  n_procData = 0;
  int idxA = 0;
  int idxB = 0;
  for (int key = 0; key < keyMax; key++) {

    _edges_intersect_res_t ** datas =
            (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht,
                                                               (void *) &key);

    int nData = PDM_hash_tab_n_data_get (ht, &key);

    for (int i = 0; i < nData; i++) {
      _edges_intersect_res_t *_data = datas[i];

      keys[n_procData]        = _data->nGEdgeA + _data->nGEdgeB;
      tIntersects[n_procData] = _data->tIntersect;

      gNumEdgeA[n_procData  ] = _data->nGEdgeA;
      gNumEdgeB[n_procData  ] = _data->nGEdgeB;
      nNewPointsA[n_procData] = _data->nNewPointsA;
      for (int j = 0; j < nNewPointsA[n_procData]; j++) {
        oNewPointsA[idxA] = _data->oNewPointsA[j];

        PDM_g_num_t gNumInB = _data->linkA[j];
        for (int k = 0; k < 3; k++) {
          coordsPointA[3*idxA + k] = _data->coordsA[3*j + k];
        }

        connectPointA[idxA] = gNumInB;
        gNumA[idxA] = _data->gNumA[j];
        uPointA[idxA] = _data->uA[j];
        idxA += 1;
      }

      nNewPointsB[n_procData] = _data->nNewPointsB;
      for (int j = 0; j < nNewPointsB[n_procData]; j++) {
        oNewPointsB[idxB] = _data->oNewPointsB[j];

        PDM_g_num_t gNumInA = _data->linkB[j];
        for (int k = 0; k < 3; k++) {
          coordsPointB[3*idxB + k] = _data->coordsB[3*j + k];
        }

        connectPointB[idxB] = gNumInA;
        gNumB[idxB] = _data->gNumB[j];
        uPointB[idxB] = _data->uB[j];
        idxB += 1;
      }

      n_procData++;
    }
  }

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 2\n"); */
  /*   fflush(stdout); */
  /* } */
  if (vb) {
    PDM_printf ("\nTable de hachage des aretes avant synchronisation, keyMax : %d\n", keyMax);

    n_procData = 0;
    idxA = 0;
    idxB = 0;
    for (int key = 0; key < keyMax; key++) {

      _edges_intersect_res_t ** datas =
        (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht,
                                                           (void *) &key);

      int nData = PDM_hash_tab_n_data_get (ht, &key);
      PDM_printf ("\n --> key : %d    nData : %d\n", key, nData);

      for (int i = 0; i < nData; i++) {
        _edges_intersect_res_t *_data = datas[i];
        PDM_printf ("\n      - idata: %d\n", i);
        PDM_printf ("       nGEdgeA :%d  nGEdgeB :%d  tIntersect : %d\n",
                    _data->nGEdgeA, _data->nGEdgeB, _data->tIntersect);

        PDM_printf ("       nNewPointsA[%d] :%d\n", n_procData, nNewPointsA[n_procData]);
        for (int j = 0; j < nNewPointsA[n_procData]; j++) {
          PDM_printf ("        oNewPointsA[%d] :%d\n", idxA, oNewPointsA[idxA]);
          PDM_printf ("        cPointA[%d] : %12.5e  %12.5e %12.5e\n",
                      idxA ,coordsPointA[3*idxA], coordsPointA[3*idxA+1], coordsPointA[3*idxA+2]);
          PDM_printf ("        gnumA[%d] :"PDM_FMT_G_NUM"\n", idxA, _data->gNumA[j]);
          PDM_printf ("        linkA[%d] :"PDM_FMT_G_NUM"\n", idxA, connectPointA[idxA]);
          PDM_printf ("        uPointA[%d] :%12.5e\n", idxA, uPointA[idxA]);
          idxA += 1;
        }

        PDM_printf ("       nNewPointsB[%d] :%d\n", n_procData, nNewPointsB[n_procData]);
        for (int j = 0; j < nNewPointsB[n_procData]; j++) {
          PDM_printf ("        oNewPointsB[%d] :%d\n", idxB, oNewPointsB[idxB]);

          PDM_printf ("        cPointB[%d] : %12.5e  %12.5e %12.5e\n", idxB
                      , coordsPointB[3*idxB], coordsPointB[3*idxB+1], coordsPointB[3*idxB+2]);

          PDM_printf ("        gnumB[%d] :"PDM_FMT_G_NUM"\n", idxB, _data->gNumB[j]);
          PDM_printf ("        linkB[%d] :"PDM_FMT_G_NUM"\n", idxB, connectPointB[idxB]);
          PDM_printf ("        uPointB[%d] :%12.5e\n", idxB, uPointB[idxB]);
          idxB += 1;
        }

        n_procData++;
      }
    }
  }

  /*
   * TODO : - A single exchange to optimize communications !
   *        - Add a revert param to PDM_part_to_block_exch :
   *          To perform block to part communication
   *          (less expensive than PDM_block_to_part))
   */

  PDM_part_to_block_t *ptb = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                     PDM_PART_TO_BLOCK_POST_MERGE,
                                                     1.,
                                                     (PDM_g_num_t **) &keys,
                                                       NULL,
                                                     &n_procData,
                                                     1,
                                                     _ei->comm);

  PDM_g_num_t *block_gnum = (PDM_g_num_t *) PDM_part_to_block_block_gnum_get (ptb);
  int n_elt_block = PDM_part_to_block_n_elt_block_get (ptb);

  /*
   * A info
   */

  int *stride_one   = malloc (sizeof(PDM_g_num_t) * n_procData);
  for (int i = 0; i < n_procData; i++) {
    stride_one[i] = 1;
  }

  int *b_tIntersects = NULL;
  int *b_stride_one = NULL;

  PDM_part_to_block_exch (ptb,
                          sizeof(int),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &stride_one,
                          (void **) &tIntersects,
                          &b_stride_one,
                          (void **) &b_tIntersects);

  free (b_stride_one);

  PDM_g_num_t *b_gNumEdgeA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &stride_one,
                         (void **) &gNumEdgeA,
                         &b_stride_one,
                         (void **) &b_gNumEdgeA);


  free (b_stride_one);
  PDM_g_num_t *b_gNumEdgeB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &stride_one,
                         (void **) &gNumEdgeB,
                         &b_stride_one,
                         (void **) &b_gNumEdgeB);

  PDM_edges_intersect_point_t *b_oNewPointsA  = NULL;
  int                         *b_nNewPointsA;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsA,
                         (void **)&oNewPointsA,
                         &b_nNewPointsA,
                         (void **)&b_oNewPointsA);

  free (b_stride_one);

  free(b_nNewPointsA);

  PDM_g_num_t *b_connectPointA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsA,
                         (void **)&connectPointA,
                         &b_nNewPointsA,
                         (void **)&b_connectPointA);

  free(b_nNewPointsA);

  PDM_g_num_t *b_gNumA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsA,
                         (void **)&gNumA,
                         &b_nNewPointsA,
                         (void **)&b_gNumA);

  free(b_nNewPointsA);

  double *b_uPointA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsA,
                         (void **)&uPointA,
                         &b_nNewPointsA,
                         (void **)&b_uPointA);

  free(b_nNewPointsA);

  for (int k = 0; k < n_procData; k++) {
    nNewPointsA[k] = 3 * nNewPointsA[k];
  }

  double *b_coordsPointA = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsA,
                         (void **)&coordsPointA,
                         &b_nNewPointsA,
                         (void **)&b_coordsPointA);

  for (int k = 0; k < n_procData; k++) {
    nNewPointsA[k] = nNewPointsA[k] / 3;
  }

  free (b_nNewPointsA);
  PDM_part_to_block_exch (ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &stride_one,
                         (void **) &nNewPointsA,
                         &b_stride_one,
                         (void **) &b_nNewPointsA);

  // int sum1=0;
  // for (int k = 0; k < n_elt_block; k++) {
  //   sum1+=b_stride_one[k];
  // }
  free (b_stride_one);

  /*
   * B info
   */

  PDM_edges_intersect_point_t *b_oNewPointsB  = NULL;
  int                         *b_nNewPointsB;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsB,
                         (void **)&oNewPointsB,
                         &b_nNewPointsB,
                         (void **)&b_oNewPointsB);

  free(b_nNewPointsB);

  PDM_g_num_t *b_connectPointB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsB,
                         (void **)&connectPointB,
                         &b_nNewPointsB,
                         (void **)&b_connectPointB);

  free(b_nNewPointsB);

  PDM_g_num_t *b_gNumB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(PDM_g_num_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsB,
                         (void **)&gNumB,
                         &b_nNewPointsB,
                         (void **)&b_gNumB);

  free(b_nNewPointsB);

  double *b_uPointB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsB,
                         (void **)&uPointB,
                         &b_nNewPointsB,
                         (void **)&b_uPointB);

  free(b_nNewPointsB);


  for (int k = 0; k < n_procData; k++) {
    nNewPointsB[k] = 3 * nNewPointsB[k];
  }

  double *b_coordsPointB = NULL;
  PDM_part_to_block_exch (ptb,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &nNewPointsB,
                         (void **)&coordsPointB,
                         &b_nNewPointsB,
                         (void **)&b_coordsPointB);

  for (int k = 0; k < n_procData; k++) {
    nNewPointsB[k] = nNewPointsB[k] / 3;
  }

  free (b_nNewPointsB);
  PDM_part_to_block_exch (ptb,
                         sizeof(int),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &stride_one,
                         (void **) &nNewPointsB,
                         &b_stride_one,
                         (void **) &b_nNewPointsB);

  free (stride_one);
  free (tIntersects);

  free (gNumEdgeA);
  free (gNumEdgeB);
  free (oNewPointsA);
  free (connectPointA);
  free (gNumA);
  free (coordsPointA);
  free (uPointA);

  free (oNewPointsB);
  free (connectPointB);
  free (gNumB);
  free (coordsPointB);
  free (uPointB);

  free (nNewPointsA);
  free (nNewPointsB);


  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 3\n"); */
  /*   fflush(stdout); */
  /* } */


  /*
   * Remove inconsistencies
   *   - Fill true intersection properties
   */

  int *b_stride_one_idx = PDM_array_new_idx_from_sizes_int(b_stride_one, n_elt_block);
  int *b_nNewPointsA_idx = PDM_array_new_idx_from_sizes_int(b_nNewPointsA, b_stride_one_idx[n_elt_block]);
  int *b_nNewPointsB_idx = PDM_array_new_idx_from_sizes_int(b_nNewPointsB, b_stride_one_idx[n_elt_block]);
  free (b_stride_one);

  int *tag = PDM_array_zeros_int(b_stride_one_idx[n_elt_block]);

  int *b_stride_one_idx_true = malloc(sizeof(int) * (n_elt_block + 1));

  int *b_tIntersects_true = malloc(sizeof(int) * b_stride_one_idx[n_elt_block]);
  PDM_g_num_t *b_gNumEdgeA_true = malloc(sizeof(PDM_g_num_t) * b_stride_one_idx[n_elt_block]);
  PDM_g_num_t *b_gNumEdgeB_true = malloc(sizeof(PDM_g_num_t) * b_stride_one_idx[n_elt_block]);

  int *b_nNewPointsA_true = malloc(sizeof(int) * b_stride_one_idx[n_elt_block]);

  PDM_edges_intersect_point_t *b_oNewPointsA_true =
          malloc(sizeof(PDM_edges_intersect_point_t) * 2 * b_stride_one_idx[n_elt_block]);
  PDM_g_num_t *b_connectPointA_true =
          malloc(sizeof(PDM_g_num_t) * 2 * b_stride_one_idx[n_elt_block]);
  PDM_g_num_t *b_gNumA_true =
          malloc(sizeof(PDM_g_num_t) * 2 * b_stride_one_idx[n_elt_block]);
  double *b_uPointA_true =
          malloc(sizeof(double) * 2 * b_stride_one_idx[n_elt_block]);
  double *b_coordsPointA_true =
          malloc(sizeof(double) * 6 * b_stride_one_idx[n_elt_block]);

  int *b_nNewPointsB_true = malloc(sizeof(int) * b_stride_one_idx[n_elt_block]);
  PDM_edges_intersect_point_t *b_oNewPointsB_true =
          malloc(sizeof(PDM_edges_intersect_point_t) * 2 * b_stride_one_idx[n_elt_block]);
  PDM_g_num_t *b_connectPointB_true =
          malloc(sizeof(PDM_g_num_t) * 2 * b_stride_one_idx[n_elt_block]);
  PDM_g_num_t *b_gNumB_true =
          malloc(sizeof(PDM_g_num_t) * 2 * b_stride_one_idx[n_elt_block]);
  double *b_uPointB_true = malloc(sizeof(double) * 2 * b_stride_one_idx[n_elt_block]);
  double *b_coordsPointB_true = malloc(sizeof(double) * 6 * b_stride_one_idx[n_elt_block]);

  int idx_true   = 0;
  int idx_newPtA = 0;
  int idx_newPtB = 0;

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 4\n"); */
  /*   fflush(stdout); */
  /* } */

  for (int i = 0; i < n_elt_block; i++) {
    b_stride_one_idx_true[i] = idx_true;
    PDM_g_num_t gNum = block_gnum[i];

    if (vb) {
      printf("\n\n+++ Cle : "PDM_FMT_G_NUM" (begin) %d %d\n", gNum, b_stride_one_idx_true[i], b_stride_one_idx[i]);
    }

    for (int j = b_stride_one_idx[i]; j < b_stride_one_idx[i+1]; j++) {
      if (tag[j] == 0) {
        tag[j] = 1;

        for (int k = j+1; k < b_stride_one_idx[i+1]; k++) {
          if ((tag[k] == 0) && (b_gNumEdgeA[j] == b_gNumEdgeA[k]) && (b_gNumEdgeB[j] == b_gNumEdgeB[k])) {
            tag[j] += 1;
            if (tag[j] > 2) {
              PDM_error(__FILE__, __LINE__, 0, "intersection calculée + de 2 fois : "
                        "Erreur ou pas erreur le max est peut-être en fait 4\n");
              abort();
            }
            tag[k] = -1;
            if (vb) {
              PDM_printf ("\n+++ b_gNumEdgeA[j] : "PDM_FMT_G_NUM" - b_gNumEdgeB[j] : "
                          PDM_FMT_G_NUM" / b_gNumEdgeA[k] : "PDM_FMT_G_NUM" - b_gNumEdgeB[k] : "PDM_FMT_G_NUM" $$$$$\n",
                        b_gNumEdgeA[j], b_gNumEdgeB[j], b_gNumEdgeA[k], b_gNumEdgeB[k]);

              printf("+++b_tIntersects[j] : %s / b_tIntersects[k] : %s\n"
                     ,_typeInter[b_tIntersects[j]+1] ,_typeInter[b_tIntersects[k]+1]);
            }
            b_tIntersects_true[idx_true] = PDM_MAX (b_tIntersects[j], b_tIntersects[k]);
            b_gNumEdgeA_true[idx_true] = b_gNumEdgeA[j];
            b_gNumEdgeB_true[idx_true] = b_gNumEdgeB[j];

            /*
             * Merge new points for A
             */

            b_nNewPointsA_true[idx_true] = PDM_MAX (b_nNewPointsA[j], b_nNewPointsA[k]);
            if ((b_nNewPointsA[k] == 2) && (b_nNewPointsA[j] == 2)) {
              int link[2] = {-1, -1};
              for (int k1 = 0; k1 < b_nNewPointsA[j]; k1++) {
                for (int k2 = 0; k2 < b_nNewPointsA[k]; k2++) {
                  if (vb) {
                    PDM_printf ("\nk1:%d k2:%d b_connectPointA[b_nNewPointsA_idx[j]+k1]:"
                                "%d b_connectPointA[b_nNewPointsA_idx[k]+k2]:%d ",
                                k1, k2, b_connectPointA[b_nNewPointsA_idx[j]+k1],
                                b_connectPointA[b_nNewPointsA_idx[k]+k2]);
                  }
                  if (b_connectPointA[b_nNewPointsA_idx[j]+k1] ==
                      b_connectPointA[b_nNewPointsA_idx[k]+k2]) {
                    if (vb) {
                      PDM_printf ("link[%d] = %d", k1, k2);
                    }
                    link[k1] = k2;
                    break;
                  }
                }
              }

              if ((link[0] == -1) || (link[1] == -1)) {
                printf("Sortie en erreur : Incoherence sur les deux sommets de B"
                          " devant être ajoutés à A.\n");
                printf("b_tIntersects[j], b_tIntersects[k] : %d %d\n", b_tIntersects[j], b_tIntersects[k]);
                printf("b_connectPointA[b_nNewPointsA_idx[j]=["PDM_FMT_G_NUM" "PDM_FMT_G_NUM"]\n",
                       b_connectPointA[b_nNewPointsA_idx[j]], b_connectPointA[b_nNewPointsA_idx[j]+1]);
                printf("b_connectPointA[b_nNewPointsA_idx[k]=["PDM_FMT_G_NUM" "PDM_FMT_G_NUM"]\n",
                       b_connectPointA[b_nNewPointsA_idx[k]], b_connectPointA[b_nNewPointsA_idx[k]+1]);
                printf("b_coordsPointA[j1] %12.5e %12.5e %12.5e\n", b_coordsPointA[3*b_nNewPointsA_idx[j]],
                       b_coordsPointA[3*b_nNewPointsA_idx[j]+1], b_coordsPointA[3*b_nNewPointsA_idx[j]+2]);
                printf("b_coordsPointA[j2] %12.5e %12.5e %12.5e\n", b_coordsPointA[3*(b_nNewPointsA_idx[j]+1)],
                       b_coordsPointA[3*(b_nNewPointsA_idx[j]+1)+1], b_coordsPointA[3*(b_nNewPointsA_idx[j]+1)+2]);
                printf("b_coordsPointA[k1] %12.5e %12.5e %12.5e\n", b_coordsPointA[3*b_nNewPointsA_idx[k]],
                       b_coordsPointA[3*b_nNewPointsA_idx[k]+1], b_coordsPointA[3*b_nNewPointsA_idx[k]+2]);
                printf("b_coordsPointA[k2] %12.5e %12.5e %12.5e\n", b_coordsPointA[3*(b_nNewPointsA_idx[k]+1)],
                       b_coordsPointA[3*(b_nNewPointsA_idx[k]+1)+1], b_coordsPointA[3*(b_nNewPointsA_idx[k]+1)+2]);
                abort();
              }

              b_oNewPointsA_true[idx_newPtA] = PDM_MAX (b_oNewPointsA[b_nNewPointsA_idx[j]],
                                                        b_oNewPointsA[b_nNewPointsA_idx[k] + link[0]]);

              int i_true1;

              if (b_oNewPointsA_true[idx_newPtA] == b_oNewPointsA[b_nNewPointsA_idx[j]]) {
                i_true1 = b_nNewPointsA_idx[j];
              }
              else {
                i_true1 = b_nNewPointsA_idx[k] + link[0];
              }

              b_oNewPointsA_true[idx_newPtA + 1] = PDM_MAX (b_oNewPointsA[b_nNewPointsA_idx[j] + 1],
                                                            b_oNewPointsA[b_nNewPointsA_idx[k] + link[1]]);

              int i_true2;

              if (b_oNewPointsA_true[idx_newPtA + 1] == b_oNewPointsA[b_nNewPointsA_idx[j]]) {
                i_true2 = b_nNewPointsA_idx[j] + 1;
              }
              else {
                i_true2 = b_nNewPointsA_idx[k] + link[1];
              }
              b_connectPointA_true[idx_newPtA    ] = b_connectPointA[i_true1];
              b_connectPointA_true[idx_newPtA + 1] = b_connectPointA[i_true2];

              b_gNumA_true[idx_newPtA    ] = b_gNumA[i_true1];
              b_gNumA_true[idx_newPtA + 1] = b_gNumA[i_true2];

              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true1 + k1];
                b_coordsPointA_true[3*(idx_newPtA+1) + k1] = b_coordsPointA[3*i_true2 + k1];
              }

              b_uPointA_true[idx_newPtA] = b_uPointA[i_true1];
              b_uPointA_true[idx_newPtA+1] = b_uPointA[i_true2];


              idx_newPtA += 2;

            }

            else if ((b_nNewPointsA[k] == 1) && (b_nNewPointsA[j] == 1)) {
              b_oNewPointsA_true[idx_newPtA] = PDM_MAX (b_oNewPointsA[b_nNewPointsA_idx[j]],
                                                        b_oNewPointsA[b_nNewPointsA_idx[k]]);
              int i_true1;
              if (b_oNewPointsA_true[idx_newPtA] == b_oNewPointsA[b_nNewPointsA_idx[j]]) {
                i_true1 = b_nNewPointsA_idx[j];
              }
              else {
                i_true1 = b_nNewPointsA_idx[k];
              }

              b_connectPointA_true[idx_newPtA    ] = b_connectPointA[i_true1];
              b_gNumA_true[idx_newPtA    ] = b_gNumA[i_true1];
              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true1 + k1];
              }

              b_uPointA_true[idx_newPtA] = b_uPointA[i_true1];

              idx_newPtA += 1;
            }

            else if (b_nNewPointsA[k] != b_nNewPointsA[j]) {
              int idx_sup;
              int idx_inf;
              if (b_nNewPointsA[j] > b_nNewPointsA[k]) {
                idx_sup = j;
                idx_inf = k;
              }
              else {
                idx_sup = k;
                idx_inf = j;
              }

              if (b_nNewPointsA[idx_inf] == 0) {

                int i_true = b_nNewPointsA_idx[idx_sup];
                b_oNewPointsA_true[idx_newPtA] = b_oNewPointsA[i_true];

                b_connectPointA_true[idx_newPtA] = b_connectPointA[i_true];
                b_gNumA_true[idx_newPtA] = b_gNumA[i_true];
                for (int k1 = 0; k1 < 3; k1++) {
                  b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true + k1];
                }

                b_uPointA_true[idx_newPtA] = b_uPointA[i_true];
                idx_newPtA += 1;

              }

              else {

                if (vb) {
                  printf("b_nNewPointsA[idx_inf], b_nNewPointsA[idx_sup] : %d %d\n",
                         b_nNewPointsA[idx_inf], b_nNewPointsA[idx_sup]);
                }
                assert (b_nNewPointsA[idx_inf] == 1);
                assert (b_nNewPointsA[idx_sup] == 2);

                /*
                 * if "idx_inf" is "Vertex A on vertex B" :
                 *    - Look for the related point in "idx_sup"
                 *    - copy "idx_inf"
                 *    - copy the second "idx_sup" vertex
                 * Otherwise :
                 *    - copy "idx_sup" vertices
                 *
                 */

                int i_true[2] = {-1, -1};

                if (b_oNewPointsA[b_nNewPointsA_idx[idx_inf]]
                    == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {

                  i_true[0] = b_nNewPointsA_idx[idx_inf];
                  if (b_connectPointA[b_nNewPointsA_idx[idx_sup]] == gNum) {
                    i_true[1] = b_nNewPointsA_idx[idx_sup];
                  }
                  else if (b_connectPointA[b_nNewPointsA_idx[idx_sup] + 1] == gNum) {
                    i_true[1] = b_nNewPointsA_idx[idx_sup] + 1;
                  }
                  assert (i_true[1] != -1);

                }

                else {
                  i_true[0] = b_nNewPointsA_idx[idx_sup];
                  i_true[1] = b_nNewPointsA_idx[idx_sup] + 1;
                }

                for (int k2 = 0; k2 < 2; k2++) {
                  b_oNewPointsA_true[idx_newPtA] = b_oNewPointsA[i_true[k2]];

                  b_connectPointA_true[idx_newPtA] = b_connectPointA[i_true[k2]];
                  b_gNumA_true[idx_newPtA] = b_gNumA[i_true[k2]];
                  for (int k1 = 0; k1 < 3; k1++) {
                    b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3 * i_true[k2] + k1];
                  }

                  b_uPointA_true[idx_newPtA] = b_uPointA[i_true[k2]];
                  idx_newPtA += 1;
                }
              }
            }

            else {
              b_nNewPointsA_true[idx_true] = 0;
            }

            /*
             * Merge new points for B
             */

            b_nNewPointsB_true[idx_true] = PDM_MAX (b_nNewPointsB[j], b_nNewPointsB[k]);
            if (vb >= 1) {
              PDM_printf ("b_nNewPointsB[k]:%d  b_nNewPointsB[j]:%d \n", b_nNewPointsB[k], b_nNewPointsB[j]);
            }

            if ((b_nNewPointsB[k] == 2) && (b_nNewPointsB[j] == 2)) {
              int link[2] = {-1, -1};
              for (int k1 = 0; k1 < b_nNewPointsB[j]; k1++) {
                for (int k2 = 0; k2 < b_nNewPointsB[k]; k2++) {
                  if (b_connectPointB[b_nNewPointsB_idx[j]+k1] ==
                      b_connectPointB[b_nNewPointsB_idx[k]+k2]) {
                    link[k1] = k2;
                    break;
                  }
                }
              }

              if ((link[0] == -1) || (link[1] == -1)) {
                printf("link[0], link[1] : %d %d\n", link[0], link[1]);
                PDM_error(__FILE__, __LINE__, 0,
                          "Sortie en erreur : Incoherence sur les deux sommets de B "
                          "devant être ajoutés à A.\n");
                abort();
              }

              b_oNewPointsB_true[idx_newPtB] = PDM_MAX (b_oNewPointsB[b_nNewPointsB_idx[j]],
                                                        b_oNewPointsB[b_nNewPointsB_idx[k] + link[0]]);
              int i_true1 = b_nNewPointsB_idx[j];
              if (b_oNewPointsB_true[idx_newPtB] == b_oNewPointsB[b_nNewPointsB_idx[k] + link[0]]) {
                i_true1 = b_nNewPointsB_idx[k] + link[0];
              }

              b_oNewPointsB_true[idx_newPtB + 1] = PDM_MAX (b_oNewPointsB[b_nNewPointsB_idx[j] + 1],
                                                            b_oNewPointsB[b_nNewPointsB_idx[k] + link[1]]);
              int i_true2 = b_nNewPointsB_idx[j] + 1;
              if (b_oNewPointsB_true[idx_newPtB + 1] == b_oNewPointsB[b_nNewPointsB_idx[k] + link[1]]) {
                i_true2 = b_nNewPointsB_idx[k] + link[1];
              }

              b_connectPointB_true[idx_newPtB    ] = b_connectPointB[i_true1];
              b_connectPointB_true[idx_newPtB + 1] = b_connectPointB[i_true2];

              b_gNumB_true[idx_newPtB    ] = b_gNumB[i_true1];
              b_gNumB_true[idx_newPtB + 1] = b_gNumB[i_true2];

              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true1 + k1];
                b_coordsPointB_true[3*(idx_newPtB+1) + k1] = b_coordsPointB[3*i_true2 + k1];
              }

              b_uPointB_true[idx_newPtB ] = b_uPointB[i_true1];
              b_uPointB_true[idx_newPtB+1] = b_uPointB[i_true2];

              idx_newPtB += 2;

            }

            else if ((b_nNewPointsB[k] == 1) && (b_nNewPointsB[j] == 1)) {
              b_oNewPointsB_true[idx_newPtB] = PDM_MAX (b_oNewPointsB[b_nNewPointsB_idx[j]],
                                                        b_oNewPointsB[b_nNewPointsB_idx[k]]);
              int i_true1 = b_nNewPointsB_idx[j];
              if (b_oNewPointsB_true[idx_newPtB] == b_oNewPointsB[b_nNewPointsB_idx[k]]) {
                i_true1 = b_nNewPointsB_idx[k];
              }

              b_connectPointB_true[idx_newPtB    ] = b_connectPointB[i_true1];
              b_gNumB_true[idx_newPtB    ] = b_gNumB[i_true1];
              for (int k1 = 0; k1 < 3; k1++) {
                b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true1 + k1];
              }

              b_uPointB_true[idx_newPtB] = b_uPointB[i_true1];

              idx_newPtB += 1;
            }

            else if (b_nNewPointsB[k] != b_nNewPointsB[j]) {
              int idx_sup;
              int idx_inf;
              if (b_nNewPointsB[j] > b_nNewPointsB[k]) {
                idx_sup = j;
                idx_inf = k;
              }
              else {
                idx_sup = k;
                idx_inf = j;
              }

              if (b_nNewPointsB[idx_inf] == 0) {

                int i_true = b_nNewPointsB_idx[idx_sup];
                b_oNewPointsB_true[idx_newPtB] = b_oNewPointsB[i_true];

                b_connectPointB_true[idx_newPtB] = b_connectPointB[i_true];
                b_gNumB_true[idx_newPtB] = b_gNumB[i_true];
                for (int k1 = 0; k1 < 3; k1++) {
                  b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true + k1];
                }

                /* for (int k1 = 0; k1 < 2; k1++) { */
                /*   b_uPointB_true[2*idx_newPtB + k1] = b_uPointB[2*i_true + k1]; */
                /* } */
                b_uPointB_true[idx_newPtB] = b_uPointB[i_true];
                idx_newPtB += 1;

              }

              else {

                assert (b_nNewPointsB[idx_inf] == 1);
                assert (b_nNewPointsB[idx_sup] == 2);

                /*
                 * if "idx_inf" is "Vertex A on vertex B" :
                 *    - Look for the related point in "idx_sup"
                 *    - copy "idx_inf"
                 *    - copy the second "idx_sup" vertex
                 * Otherwise :
                 *    - copy "idx_sup" vertices
                 *
                 */

                int i_true[2] = {-1, -1};

                if (b_oNewPointsB[b_nNewPointsB_idx[idx_inf]]
                    == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {

                  i_true[0] = b_nNewPointsB_idx[idx_inf];
                  if (b_connectPointB[b_nNewPointsB_idx[idx_sup]] == gNum) {
                    i_true[1] = b_nNewPointsB_idx[idx_sup];
                  }
                  else if (b_connectPointB[b_nNewPointsB_idx[idx_sup] + 1] == gNum) {
                    i_true[1] = b_nNewPointsB_idx[idx_sup] + 1;
                  }
                  assert (i_true[1] != -1);

                }

                else {
                  i_true[0] = b_nNewPointsB_idx[idx_sup];
                  i_true[1] = b_nNewPointsB_idx[idx_sup] + 1;
                }

                for (int k2 = 0; k2 < 2; k2++) {
                  b_oNewPointsB_true[idx_newPtB] = b_oNewPointsB[i_true[k2]];

                  b_connectPointB_true[idx_newPtB] = b_connectPointB[i_true[k2]];
                  b_gNumB_true[idx_newPtB] = b_gNumB[i_true[k2]];
                  for (int k1 = 0; k1 < 3; k1++) {
                    b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3 * i_true[k2] + k1];
                  }

                  /* for (int k1 = 0; k1 < 2; k1++) { */
                  /*   b_uPointB_true[2*idx_newPtB + k1] = b_uPointB[2 * i_true[k2] + k1]; */
                  /* } */
                  b_uPointB_true[idx_newPtB] = b_uPointB[i_true[k2]];
                  idx_newPtB += 1;
                }
              }
            }

            else {
              b_nNewPointsB_true[idx_true] = 0;
            }

            idx_true += 1;

          }
        }

        /*
         * Copy intersection if only one computation
         */

        if (tag[j] == 1) {

          b_tIntersects_true[idx_true] = b_tIntersects[j];
          b_gNumEdgeA_true[idx_true] = b_gNumEdgeA[j];
          b_gNumEdgeB_true[idx_true] = b_gNumEdgeB[j];

          b_nNewPointsA_true[idx_true] = b_nNewPointsA[j];

          for (int k = 0; k < b_nNewPointsA[j]; k++) {
            int i_true = b_nNewPointsA_idx[j] + k ;
            b_oNewPointsA_true[idx_newPtA] = b_oNewPointsA[i_true];

            b_connectPointA_true[idx_newPtA] = b_connectPointA[i_true];
            b_gNumA_true[idx_newPtA] = b_gNumA[i_true];
            for (int k1 = 0; k1 < 3; k1++) {
              b_coordsPointA_true[3*idx_newPtA + k1] = b_coordsPointA[3*i_true + k1];
            }

            b_uPointA_true[idx_newPtA] = b_uPointA[i_true];

            idx_newPtA += 1;
          }

          b_nNewPointsB_true[idx_true] = b_nNewPointsB[j];

          for (int k = 0; k < b_nNewPointsB[j]; k++) {
            int i_true = b_nNewPointsB_idx[j] + k ;

            b_oNewPointsB_true[idx_newPtB] = b_oNewPointsB[i_true];
            b_connectPointB_true[idx_newPtB] = b_connectPointB[i_true];
            b_gNumB_true[idx_newPtB] = b_gNumB[i_true];
            for (int k1 = 0; k1 < 3; k1++) {
              b_coordsPointB_true[3*idx_newPtB + k1] = b_coordsPointB[3*i_true + k1];
            }

            b_uPointB_true[idx_newPtB] = b_uPointB[i_true];

            idx_newPtB += 1;
          }

          idx_true += 1;

        }
      }
    }
    if (vb) {
      printf("+++ Cle : "PDM_FMT_G_NUM" (fin)\n", gNum);
    }
  }


  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 5\n"); */
  /*   fflush(stdout); */
  /* } */

  b_stride_one_idx_true[n_elt_block] = idx_true;

   if (vb) {
     printf("\n\n------------\n Bilan\n");

     for (int i = 0; i < n_elt_block; i++) {
       PDM_g_num_t gNum = block_gnum[i];
       printf("\n\n+++ Cle : "PDM_FMT_G_NUM" (begin) plage init : %d -> %d, plage true : %d -> %d\n",
              gNum,
              b_stride_one_idx[i], b_stride_one_idx[i+1],
              b_stride_one_idx_true[i], b_stride_one_idx_true[i+1]);
       for (int j = b_stride_one_idx[i]; j < b_stride_one_idx[i+1]; j++) {
         printf ("   - Init gnumA gnumB : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",b_gNumEdgeA[j], b_gNumEdgeB[j]);
       }
       for (int j = b_stride_one_idx_true[i]; j < b_stride_one_idx_true[i+1]; j++) {
         printf ("   - Init b_gNumEdgeA_true, b_gNumEdgeB_true : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",b_gNumEdgeA_true[j], b_gNumEdgeB_true[j]);
       }
       printf("+++ Cle : "PDM_FMT_G_NUM" (end)\n", gNum);

     }
   }

  free (tag);
  free (b_nNewPointsB_idx);
  free (b_nNewPointsA_idx);

  free (b_connectPointA);
  free (b_connectPointB);
  free (b_gNumA);
  free (b_gNumB);
  free (b_coordsPointA);
  free (b_coordsPointB);
  free (b_gNumEdgeA);
  free (b_gNumEdgeB);
  free (b_nNewPointsA);
  free (b_nNewPointsB);
  free (b_oNewPointsA);
  free (b_oNewPointsB);
  free (b_stride_one_idx);
  free (b_tIntersects);

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 6\n"); */
  /*   fflush(stdout); */
  /* } */

  /*
   * Perform absolute number of new points
   *  - grace a b_oNewPointsA_true et b_oNewPointsA_true
   *  - verifier qu'on a le meme nombre de nouveaux points de chaque côté !
   *  - mpi_all_reduce ou equivalent
   *
   */

  PDM_g_num_t nNewPtsA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      nNewPtsA += 1;
    }
  }

  PDM_g_num_t beg_nNewPtsA = 0;
  PDM_MPI_Request request1;
  PDM_MPI_Iscan(&nNewPtsA, &beg_nNewPtsA, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request1);

  PDM_g_num_t nNewPtsB = 0;
  for (int i = 0; i < idx_newPtB; i++) {
    if (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      nNewPtsB += 1;
    }
  }

  PDM_g_num_t beg_nNewPtsB = 0;
  PDM_MPI_Request request2;
  PDM_MPI_Iscan(&nNewPtsB, &beg_nNewPtsB, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request2);

  int nNewPtsFromBForA = 0;
  int nPtsFromBForA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if ((b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) ||
        (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      nPtsFromBForA += 1;
      if (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
        nNewPtsFromBForA += 1;
      }
    }
  }
  PDM_UNUSED(nPtsFromBForA);
  PDM_UNUSED(nNewPtsFromBForA);

  PDM_MPI_Wait (&request1);

  PDM_g_num_t end_nNewPtsA = beg_nNewPtsA + nAbsVtxA + 1;
  beg_nNewPtsA += -nNewPtsA + 1 + nAbsVtxA;

  PDM_g_num_t beg_nPtsFromBForA = end_nNewPtsA;

  PDM_MPI_Bcast (&beg_nPtsFromBForA, 1, PDM__PDM_MPI_G_NUM, lastRank, _ei->comm);

  nNewPtsA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      b_gNumA_true[i] = beg_nNewPtsA + nNewPtsA;
      nNewPtsA += 1;
    }
  }

  int nNewPtsFromAForB = 0;
  int nPtsFromAForB = 0;
  for (int i = 0; i < idx_newPtB; i++) {
    if ((b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) ||
        (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      nPtsFromAForB += 1;
      if (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {
        nNewPtsFromAForB += 1;
      }
    }
  }
  PDM_UNUSED(nNewPtsFromAForB);

  PDM_MPI_Wait (&request2);
  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 7\n"); */
  /*   fflush(stdout); */
  /* } */

  //PDM_g_num_t *b_gNumPointsB_true = malloc(sizeof(PDM_g_num_t) * idx_newPtB);
  PDM_g_num_t end_nNewPtsB = beg_nNewPtsB + nAbsVtxB + 1;
  beg_nNewPtsB += -nNewPtsB + 1 + nAbsVtxB;

  PDM_g_num_t beg_nPtsFromAForB = end_nNewPtsB;

  PDM_MPI_Bcast (&beg_nPtsFromAForB, 1, PDM__PDM_MPI_G_NUM, lastRank, _ei->comm);

  nNewPtsB = 0;
  for (int i = 0; i < idx_newPtB; i++) {
    if (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_NEW) {
      b_gNumB_true[i] = beg_nNewPtsB + nNewPtsB;
      nNewPtsB += 1;
    }
  }

  /*
   * New points from B for A
   *    - Compute absolute number
   *    - Synchronize coordinates
   */

  double *b_cNewPointsA_true_pack = malloc (sizeof(double) * 3 * nPtsFromBForA);
  PDM_edges_intersect_point_t *b_oNewPointsA_true_pack =
    malloc (sizeof(PDM_edges_intersect_point_t) * nPtsFromBForA);
  PDM_g_num_t *b_lNewPointsA_true_pack = malloc (sizeof(PDM_g_num_t) * nPtsFromBForA);

  nPtsFromBForA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if ((b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) ||
        (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      for (int k = 0; k < 3; k++) {
        b_cNewPointsA_true_pack[3*nPtsFromBForA + k] = b_coordsPointA_true[3*i+k];
      }
      b_oNewPointsA_true_pack[nPtsFromBForA] = b_oNewPointsA_true[i];
      b_lNewPointsA_true_pack[nPtsFromBForA] = b_connectPointA_true[i];
      nPtsFromBForA += 1;
    }
  }

  fflush(stdout);

  PDM_part_to_block_t *ptbBForA = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                            PDM_PART_TO_BLOCK_POST_MERGE,
                                                            1.,
                                                            (PDM_g_num_t **) &b_lNewPointsA_true_pack, // Max des points selectionnes
                                                            NULL,
                                                            &nPtsFromBForA,
                                                            1,
                                                            _ei->comm);

  int n_BForA_gnum = PDM_part_to_block_n_elt_block_get (ptbBForA);

  int *b_stride_packA = PDM_array_const_int(nPtsFromBForA, 3);

  int *b_b_stride_packA = NULL;
  double *b_b_cNewPointsA_true_pack = NULL;
  PDM_part_to_block_exch (ptbBForA,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          1,
                          &b_stride_packA,
                          (void **)&b_cNewPointsA_true_pack,
                          &b_b_stride_packA,
                          (void **)&b_b_cNewPointsA_true_pack);

  free (b_cNewPointsA_true_pack);
  free (b_b_stride_packA);
  for (int i = 0; i < nPtsFromBForA; i++) {
    b_stride_packA[i] = 1;
  }

  PDM_edges_intersect_point_t *b_b_oNewPointsA_true_pack = NULL;
  PDM_part_to_block_exch (ptbBForA,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &b_stride_packA,
                         (void **)&b_oNewPointsA_true_pack,
                         &b_b_stride_packA,
                         (void **)&b_b_oNewPointsA_true_pack);

  free (b_oNewPointsA_true_pack);
  free (b_stride_packA);

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 8\n"); */
  /*   fflush(stdout); */
  /* } */

  /* Synchronize coordinates for A */

  int *b_b_idx_packA = malloc(sizeof(int) * (n_BForA_gnum +1));
  b_b_idx_packA[0] = 0;
  for (int i = 0; i < n_BForA_gnum; i++) {
    b_b_idx_packA[i+1] = b_b_idx_packA[i] + b_b_stride_packA[i];
    b_b_stride_packA[i] = 1;
  }

  PDM_g_num_t n_active_BForA_gnum = 0;
  for (int i = 0; i < n_BForA_gnum; i++) {
    int beg = b_b_idx_packA[i];
    int end = b_b_idx_packA[i+1];
    int cpt = 0;
    int state = 0;
    PDM_edges_intersect_point_t eip;

    for (int j = beg; j < end; j++) {
      if (state == 0) {
        state = 1;
        cpt   = 1;
        eip   = b_b_oNewPointsA_true_pack[j];
        if (eip == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
          n_active_BForA_gnum += 1;
        }
      }
      else {
        if (eip != b_b_oNewPointsA_true_pack[j]) {
          PDM_error(__FILE__, __LINE__, 0, "Error PDM_edges_intersection : Inconsistencies : "
                    "A B point is on a A vertex and on a A edge");
          abort();
        }
        cpt += 1;

        for (int k = 0; k < 3; k++) {
          b_b_cNewPointsA_true_pack[3 * beg + k] += b_b_cNewPointsA_true_pack[3 * j + k];
        }
      }
    }
    for (int k = 0; k < 3; k++) {
      b_b_cNewPointsA_true_pack[3 * beg + k] = b_b_cNewPointsA_true_pack[3 * beg + k] / cpt;
    }
  }

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 9\n"); */
  /*   fflush(stdout); */
  /* } */

  /*
   * Compress A coordinates
   */

  int idx = 0;
  for (int i = 0; i < n_BForA_gnum; i++) {
    int beg = b_b_idx_packA[i];
    int end = b_b_idx_packA[i+1];
    if (end > beg) {
      for (int k = 0; k < 3; k++) {
        b_b_cNewPointsA_true_pack[idx++] = b_b_cNewPointsA_true_pack[3*beg + k];
      }
    }
  }

  /*
   * Define and copy A Numabs
   */

  PDM_g_num_t *b_b_gNumVtxFromBForA = malloc (sizeof(PDM_g_num_t) * b_b_idx_packA[n_BForA_gnum]);

  PDM_g_num_t beg_n_BForA_gnum;

  PDM_MPI_Iscan (&n_active_BForA_gnum, &beg_n_BForA_gnum, 1,
                 PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request1);

  PDM_MPI_Wait (&request1);

  //beg_n_BForA_gnum += -n_active_BForA_gnum + beg_nPtsFromBForA + 1;
  beg_n_BForA_gnum += -n_active_BForA_gnum + beg_nPtsFromBForA;

  idx = 0;

  for (int i = 0; i < n_BForA_gnum; i++) {
    int beg = b_b_idx_packA[i];
    int end = b_b_idx_packA[i+1];
    if (end > beg) {
      if (b_b_oNewPointsA_true_pack[beg] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
        b_b_gNumVtxFromBForA[idx++] = -1;
      }
      else {
        b_b_gNumVtxFromBForA[idx++] = beg_n_BForA_gnum++;
      }
    }
  }

  free (b_b_oNewPointsA_true_pack);
  free (b_b_idx_packA);

  beg_n_BForA_gnum += -1;

  *nAbsNewVtxA = beg_n_BForA_gnum;

  PDM_MPI_Bcast(nAbsNewVtxA, 1, PDM__PDM_MPI_G_NUM,
                lastRank,
                _ei->comm);

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 10\n"); */
  /*   fflush(stdout); */
  /* } */

  /*
   * New points from A for B  :
   *   - Compute absolute number
   *   - Synchronize coordinates
   */

  double *b_cNewPointsB_true_pack = malloc (sizeof(double) * 3 * nPtsFromAForB);
  PDM_edges_intersect_point_t *b_oNewPointsB_true_pack =
          malloc (sizeof(PDM_edges_intersect_point_t) * nPtsFromAForB);
  PDM_g_num_t *b_lNewPointsB_true_pack = malloc (sizeof(PDM_g_num_t) * nPtsFromAForB);

  nPtsFromAForB = 0;
  for (int i = 0; i < idx_newPtB; i++) {
    if ((b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) ||
        (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      for (int k = 0; k < 3; k++) {
        b_cNewPointsB_true_pack[3*nPtsFromAForB + k] = b_coordsPointB_true[3*i+k];
      }
      b_oNewPointsB_true_pack[nPtsFromAForB] = b_oNewPointsB_true[i];
      b_lNewPointsB_true_pack[nPtsFromAForB] = b_connectPointB_true[i];
      nPtsFromAForB += 1;
    }
  }

  PDM_part_to_block_t *ptbAForB = PDM_part_to_block_create (PDM_PART_TO_BLOCK_DISTRIB_ALL_PROC,
                                                            PDM_PART_TO_BLOCK_POST_MERGE,
                                                            1.,
                                                            (PDM_g_num_t **) &b_lNewPointsB_true_pack,
                                                            NULL,
                                                            &nPtsFromAForB,
                                                            1,
                                                            _ei->comm);

  int n_AForB_gnum = PDM_part_to_block_n_elt_block_get (ptbAForB);

  int *b_stride_packB = PDM_array_const_int(nPtsFromAForB, 3);

  int *b_b_stride_packB = NULL;
  double *b_b_cNewPointsB_true_pack = NULL;
  PDM_part_to_block_exch (ptbAForB,
                         sizeof(double),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &b_stride_packB,
                         (void **)&b_cNewPointsB_true_pack,
                         &b_b_stride_packB,
                         (void **)&b_b_cNewPointsB_true_pack);

  free (b_cNewPointsB_true_pack);

  free (b_b_stride_packB);
  for (int i = 0; i < nPtsFromAForB; i++) {
    b_stride_packB[i] = 1;
  }

  PDM_edges_intersect_point_t *b_b_oNewPointsB_true_pack = NULL;
  PDM_part_to_block_exch (ptbAForB,
                         sizeof(PDM_edges_intersect_point_t),
                         PDM_STRIDE_VAR_INTERLACED,
                         1,
                         &b_stride_packB,
                         (void **)&b_oNewPointsB_true_pack,
                         &b_b_stride_packB,
                         (void **)&b_b_oNewPointsB_true_pack);
  free (b_stride_packB);
  free (b_oNewPointsB_true_pack);


  int *b_b_idx_packB = malloc(sizeof(int) * (n_AForB_gnum +1));
  b_b_idx_packB[0] = 0;
  for (int i = 0; i < n_AForB_gnum; i++) {
    b_b_idx_packB[i+1] = b_b_idx_packB[i] + b_b_stride_packB[i];
    b_b_stride_packB[i] = 1;
  }

  /*
   * Synchronize B coordinates
   */

  PDM_g_num_t n_active_AForB_gnum = 0;
  for (int i = 0; i < n_AForB_gnum; i++) {
    int beg = b_b_idx_packB[i];
    int end = b_b_idx_packB[i+1];
    int cpt = 0;
    int state = 0;
    PDM_edges_intersect_point_t eip;

    int cpte = 0;
    for (int j = beg; j < end; j++) {
      if (state == 0) {
        state = 1;
        cpt   = 1;
        eip   = b_b_oNewPointsB_true_pack[j];
        if (eip == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {
          n_active_AForB_gnum += 1;
        }
      }
      else {
        if (eip != b_b_oNewPointsB_true_pack[j]) {
          PDM_error(__FILE__, __LINE__, 0,
                    "Error PDM_edges_intersection : Inconsistencies : "
                    "A B point is on a A vertex and on a A edge");
          cpte++;
        }
        cpt += 1;

        for (int k = 0; k < 3; k++) {
          b_b_cNewPointsB_true_pack[3 * beg + k] += b_b_cNewPointsB_true_pack[3 * j + k];
        }
      }
    }
    if (cpte >= 1)
      abort();

    for (int k = 0; k < 3; k++) {
      b_b_cNewPointsB_true_pack[3 * beg + k] = b_b_cNewPointsB_true_pack[3 * beg + k] / cpt;
    }
  }

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 11\n"); */
  /*   fflush(stdout); */
  /* } */

  /*
   * Compress B coordinates
   */

  idx = 0;
  for (int i = 0; i < n_AForB_gnum; i++) {
    int beg = b_b_idx_packB[i];
    int end = b_b_idx_packB[i+1];
    if (end > beg) {
      for (int k = 0; k < 3; k++) {
        b_b_cNewPointsB_true_pack[idx++] = b_b_cNewPointsB_true_pack[3*beg + k];
      }
    }
  }

  /*
   * Define and copy B Numabs
   */

  PDM_g_num_t *b_b_gNumVtxFromAForB = malloc (sizeof(PDM_g_num_t) * b_b_idx_packB[n_AForB_gnum]);

  PDM_g_num_t beg_n_AForB_gnum;

  PDM_MPI_Iscan (&n_active_AForB_gnum, &beg_n_AForB_gnum, 1,
             PDM__PDM_MPI_G_NUM, PDM_MPI_SUM, _ei->comm, &request1);

  PDM_MPI_Wait (&request1);

  //beg_n_AForB_gnum += -n_active_AForB_gnum + beg_nPtsFromAForB + 1;
  beg_n_AForB_gnum += -n_active_AForB_gnum + beg_nPtsFromAForB;

  idx = 0;

  for (int i = 0; i < n_AForB_gnum; i++) {
    int beg = b_b_idx_packB[i];
    int end = b_b_idx_packB[i+1];
    if (end > beg) {
      if (b_b_oNewPointsB_true_pack[beg] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB) {
        b_b_gNumVtxFromAForB[idx++] = -1;
      }
      else {
        b_b_gNumVtxFromAForB[idx++] = beg_n_AForB_gnum++;
      }
    }
  }

  free (b_b_idx_packB);
  free (b_b_oNewPointsB_true_pack);

  beg_n_AForB_gnum += -1;

  if (lastRank == i_rank) {
    *nAbsNewVtxB = beg_n_AForB_gnum;
  }

  PDM_MPI_Bcast(nAbsNewVtxB, 1,
                PDM__PDM_MPI_G_NUM,
                lastRank,
                _ei->comm);

  /*
   * block to part the true intersection and absolute number
   *   - For A
   *   - For B
   */

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 12\n"); */
  /*   fflush(stdout); */
  /* } */

  PDM_g_num_t    *blockDistribIdxA = PDM_part_to_block_distrib_index_get (ptbBForA);
  PDM_g_num_t    *ptbBForA_block_gnum = PDM_part_to_block_block_gnum_get (ptbBForA);

  PDM_block_to_part_t *btpBForA = PDM_block_to_part_create (blockDistribIdxA,
                                                            (const PDM_g_num_t **) &b_lNewPointsA_true_pack,
                                                            &nPtsFromBForA,
                                                            1,
                                                            _ei->comm);

  int btpBForA_n_elt_block = blockDistribIdxA[i_rank+1] - blockDistribIdxA[i_rank];
  int *_b_b_stride_packA = PDM_array_zeros_int(btpBForA_n_elt_block);

  for (int i = 0; i < n_BForA_gnum; i++) {
    int idx_block = PDM_block_to_part_gnum_idx_get(btpBForA, ptbBForA_block_gnum[i]);
    _b_b_stride_packA[idx_block] = b_b_stride_packA[i];
  }

  free (b_b_stride_packA);

  int    **part_strideA = NULL;
  double **b_cNewPointsA_true_merge = NULL;
  PDM_g_num_t  **b_cNewPointsA_true_gnum = NULL;

  PDM_block_to_part_exch (btpBForA,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          _b_b_stride_packA,
                          (void *) b_b_gNumVtxFromBForA,
                          &part_strideA,
                          (void ***) &b_cNewPointsA_true_gnum);

  for (int i = 0; i < btpBForA_n_elt_block; i++) {
    _b_b_stride_packA[i] *= 3;
  }

  free (part_strideA[0]);
  free (part_strideA);

  PDM_block_to_part_exch (btpBForA,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          _b_b_stride_packA,
                           (void *) b_b_cNewPointsA_true_pack,
                          &part_strideA,
                          (void ***) &b_cNewPointsA_true_merge);

  free (part_strideA[0]);
  free (part_strideA);

  nPtsFromBForA = 0;
  for (int i = 0; i < idx_newPtA; i++) {
    if ((b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) ||
        (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      for (int k = 0; k < 3; k++) {
        b_coordsPointA_true[3*i+k] = (*b_cNewPointsA_true_merge)[3*nPtsFromBForA + k];
      }
      if (b_oNewPointsA_true[i] == PDM_EDGES_INTERSECT_POINT_VTXB_ON_EDGEA) {
        b_gNumA_true[i] = (*b_cNewPointsA_true_gnum)[nPtsFromBForA];
      }
      nPtsFromBForA += 1;
    }
  }

  PDM_g_num_t    *blockDistribIdxB = PDM_part_to_block_distrib_index_get (ptbAForB);
  PDM_g_num_t    *ptbAForB_block_gnum = PDM_part_to_block_block_gnum_get (ptbAForB);

  PDM_block_to_part_t *btpAForB = PDM_block_to_part_create (blockDistribIdxB,
                                                           (const PDM_g_num_t **) &b_lNewPointsB_true_pack,
                                                           &nPtsFromAForB,
                                                            1,
                                                            _ei->comm);

  free (b_lNewPointsB_true_pack);

  int btpAForB_n_elt_block = blockDistribIdxB[i_rank+1] - blockDistribIdxB[i_rank];
  int *_b_b_stride_packB = PDM_array_zeros_int(btpAForB_n_elt_block);

  for (int i = 0; i < n_AForB_gnum; i++) {
    int idx_block = PDM_block_to_part_gnum_idx_get(btpAForB, ptbAForB_block_gnum[i]);
    _b_b_stride_packB[idx_block] = b_b_stride_packB[i];
  }

  int    **part_strideB = NULL;
  double **b_cNewPointsB_true_merge = NULL;
  PDM_g_num_t  **b_cNewPointsB_true_gnum = NULL;

  PDM_block_to_part_exch (btpAForB,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          _b_b_stride_packB,
                          (void *) b_b_gNumVtxFromAForB,
                          &part_strideB,
                          (void ***) &b_cNewPointsB_true_gnum);

  for (int i = 0; i < btpAForB_n_elt_block; i++) {
    _b_b_stride_packB[i] *= 3;
  }

  free (part_strideB[0]);
  free (part_strideB);

  PDM_block_to_part_exch (btpAForB,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          _b_b_stride_packB,
                           (void *) b_b_cNewPointsB_true_pack,
                          &part_strideB,
                          (void ***) &b_cNewPointsB_true_merge);

  free (part_strideB[0]);
  free (part_strideB);
  nPtsFromAForB = 0;

  for (int i = 0; i < idx_newPtB; i++) {
    if ((b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) ||
        (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_VTXB)) {
      for (int k = 0; k < 3; k++) {
        b_coordsPointB_true[3*i+k] = (*b_cNewPointsB_true_merge)[3*nPtsFromAForB + k];
      }
      if (b_oNewPointsB_true[i] == PDM_EDGES_INTERSECT_POINT_VTXA_ON_EDGEB) {
        b_gNumB_true[i] = (*b_cNewPointsB_true_gnum)[nPtsFromAForB];
      }
      nPtsFromAForB += 1;
    }
  }

  free (b_cNewPointsB_true_gnum[0]);
  free (b_cNewPointsB_true_gnum);

  free (b_cNewPointsB_true_merge[0]);
  free (b_cNewPointsB_true_merge);

  free (b_cNewPointsA_true_gnum[0]);
  free (b_cNewPointsA_true_gnum);

  free (b_cNewPointsA_true_merge[0]);
  free (b_cNewPointsA_true_merge);

  PDM_block_to_part_free (btpBForA);
  PDM_part_to_block_free (ptbBForA);

  if (b_lNewPointsA_true_pack != NULL)
    free (b_lNewPointsA_true_pack);

  PDM_block_to_part_free (btpAForB);
  PDM_part_to_block_free (ptbAForB);

  free (b_uPointA);
  free (b_uPointB);

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 13\n"); */
  /*   fflush(stdout); */
  /* } */

  /*
   * Update _edges_intersect_structure (pdm_block_to_part)
   *  transfer :
   *
   * b_stride_one_idx_true;
   * b_tIntersects_true;
   * b_gNumEdgeA_true;
   * b_nNewPointsA_true;
   * b_oNewPointsA_true;
   * b_connectPointA_true;
   * b_uPointA_true;
   * b_coordsPointA_true;
   * b_gNumA_true;
   * b_nNewPointsB_true;
   * b_oNewPointsB_true;
   * b_connectPointB_true;
   * b_uPointB_true;
   * b_coordsPointB_true;
   * b_gNumB_true;
   *
   */

  PDM_g_num_t *blockDistribIdx = PDM_part_to_block_distrib_index_get (ptb);

  PDM_block_to_part_t *btp = PDM_block_to_part_create (blockDistribIdx,
                                                       (const PDM_g_num_t **) &keys,
                                                       &n_procData,
                                                       1,
                                                       _ei->comm);
  free (keys);

  if (vb) {
    printf("\n\n------------\n Bilan 2\n");

    for (int i = 0; i < n_elt_block; i++) {
      PDM_g_num_t gNum = block_gnum[i];

      printf("\n\n+++ Cle : "PDM_FMT_G_NUM" (begin)  plage true : %d -> %d\n",
             gNum,
             b_stride_one_idx_true[i], b_stride_one_idx_true[i+1]);

      for (int j = b_stride_one_idx_true[i]; j < b_stride_one_idx_true[i+1]; j++) {
        printf ("   - Init b_gNumEdgeA_true, b_gNumEdgeB_true : "PDM_FMT_G_NUM" "PDM_FMT_G_NUM"\n",b_gNumEdgeA_true[j], b_gNumEdgeB_true[j]);
      }
      printf("+++ Cle : "PDM_FMT_G_NUM" (end)\n", gNum);

    }
  }

  int n_elt_block_from_ptb = n_elt_block;
  n_elt_block = blockDistribIdx[i_rank+1] - blockDistribIdx[i_rank];

  int *b_stride_one_true = PDM_array_zeros_int(n_elt_block);


  for (int i = 0; i < n_elt_block_from_ptb; i++) {
    int idx_block = PDM_block_to_part_gnum_idx_get(btp, block_gnum[i]);
    b_stride_one_true [idx_block] = b_stride_one_idx_true[i+1] - b_stride_one_idx_true[i];
  }

  int **r_stride_one_true;
  int **r_tIntersects_true;

  if (vb) {
    printf("blockDistribIdx %d : ", n_elt_block);
    for (int i = 0; i < sMSGComm + 1; i++) {
      printf(" "PDM_FMT_G_NUM, blockDistribIdx[i]);
    }
    printf("\n");


    printf("b_tIntersects_true :");
    for (int i = 0; i < b_stride_one_idx_true[n_elt_block_from_ptb]; i++) {
      printf (" %d : ", b_tIntersects_true[i]);
    }
    printf("\n");
  }

  PDM_block_to_part_exch (btp,
                          sizeof(int),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stride_one_true,
                          (void *) b_tIntersects_true,
                          &r_stride_one_true,
                          (void ***) &r_tIntersects_true);
  free (*r_stride_one_true);
  free (r_stride_one_true);

  PDM_g_num_t **r_gNumEdgeA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stride_one_true,
                          (void *) b_gNumEdgeA_true,
                          &r_stride_one_true,
                          (void ***) &r_gNumEdgeA_true);
  free (*r_stride_one_true);
  free (r_stride_one_true);

  PDM_g_num_t **r_gNumEdgeB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stride_one_true,
                          (void *) b_gNumEdgeB_true,
                          &r_stride_one_true,
                          (void ***) &r_gNumEdgeB_true);
  free (*r_stride_one_true);
  free (r_stride_one_true);

  int **r_nNewPointsA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(int),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stride_one_true,
                          (void *) b_nNewPointsA_true,
                          &r_stride_one_true,
                          (void ***) &r_nNewPointsA_true);
  free (*r_stride_one_true);
  free (r_stride_one_true);

  int *b_stridePtsADep_true = PDM_array_zeros_int(n_elt_block);

  for (int i = 0; i < n_elt_block_from_ptb; i++) {
    int idx_block = PDM_block_to_part_gnum_idx_get(btp, block_gnum[i]);
    b_stridePtsADep_true[idx_block] = 0;
    for (int j = b_stride_one_idx_true[i]; j < b_stride_one_idx_true[i+1]; j++) {
      b_stridePtsADep_true[idx_block] += b_nNewPointsA_true[j];
    }
  }

  int **r_stridePtsADep_true = NULL;
  PDM_edges_intersect_point_t **r_oNewPointsA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_edges_intersect_point_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsADep_true,
                          (void *) b_oNewPointsA_true,
                          &r_stridePtsADep_true,
                          (void ***) &r_oNewPointsA_true);

  free (*r_stridePtsADep_true);
  free (r_stridePtsADep_true);
  PDM_g_num_t **r_connectPointA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsADep_true,
                          (void *) b_connectPointA_true,
                          &r_stridePtsADep_true,
                          (void ***) &r_connectPointA_true);

  free (*r_stridePtsADep_true);
  free (r_stridePtsADep_true);
  double **r_uPointA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsADep_true,
                          (void *) b_uPointA_true,
                          &r_stridePtsADep_true,
                          (void ***) &r_uPointA_true);
  free (b_uPointA_true);
  free (*r_stridePtsADep_true);
  free (r_stridePtsADep_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsADep_true[i] *= 3;
  }
  double **r_coordsPointA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsADep_true,
                          (void *) b_coordsPointA_true,
                          &r_stridePtsADep_true,
                          (void ***) &r_coordsPointA_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsADep_true[i] = b_stridePtsADep_true[i]/3;
  }

  free (*r_stridePtsADep_true);
  free (r_stridePtsADep_true);
  PDM_g_num_t **r_gNumA_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsADep_true,
                          (void *) b_gNumA_true,
                          &r_stridePtsADep_true,
                          (void ***) &r_gNumA_true);

  int **r_nNewPointsB_true = NULL;
  PDM_block_to_part_exch (btp,
                           sizeof(int),
                           PDM_STRIDE_VAR_INTERLACED,
                           b_stride_one_true,
                           (void *) b_nNewPointsB_true,
                           &r_stride_one_true,
                           (void ***) &r_nNewPointsB_true);

  free (b_stride_one_true);

  int *b_stridePtsBDep_true = PDM_array_zeros_int(n_elt_block);


  for (int i = 0; i < n_elt_block_from_ptb; i++) {
    int idx_block = PDM_block_to_part_gnum_idx_get(btp, block_gnum[i]);
    b_stridePtsBDep_true[idx_block] = 0;
    for (int j = b_stride_one_idx_true[i]; j < b_stride_one_idx_true[i+1]; j++) {
      b_stridePtsBDep_true[idx_block] += b_nNewPointsB_true[j];
    }
  }

  int **r_stridePtsBDep_true = NULL;
  PDM_edges_intersect_point_t **r_oNewPointsB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_edges_intersect_point_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsBDep_true,
                          (void *) b_oNewPointsB_true,
                          &r_stridePtsBDep_true,
                          (void ***) &r_oNewPointsB_true);

  free (*r_stridePtsBDep_true);
  free (r_stridePtsBDep_true);
  PDM_g_num_t **r_connectPointB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsBDep_true,
                          (void *) b_connectPointB_true,
                          &r_stridePtsBDep_true,
                          (void ***) &r_connectPointB_true);

  free (*r_stridePtsBDep_true);
  free (r_stridePtsBDep_true);
  double **r_uPointB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsBDep_true,
                          (void *) b_uPointB_true,
                          &r_stridePtsBDep_true,
                          (void ***) &r_uPointB_true);
  free (b_uPointB_true);
  free (*r_stridePtsBDep_true);
  free (r_stridePtsBDep_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsBDep_true[i] *= 3;
  }
  double **r_coordsPointB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(double),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsBDep_true,
                          (void *) b_coordsPointB_true,
                          &r_stridePtsBDep_true,
                          (void ***) &r_coordsPointB_true);
  for (int i = 0; i < n_elt_block; i++) {
    b_stridePtsBDep_true[i] = b_stridePtsBDep_true[i]/3;
  }

  free (*r_stridePtsBDep_true);
  free (r_stridePtsBDep_true);
  PDM_g_num_t **r_gNumB_true = NULL;
  PDM_block_to_part_exch (btp,
                          sizeof(PDM_g_num_t),
                          PDM_STRIDE_VAR_INTERLACED,
                          b_stridePtsBDep_true,
                          (void *) b_gNumB_true,
                          &r_stridePtsBDep_true,
                          (void ***) &r_gNumB_true);

  free (b_stridePtsBDep_true);

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 14\n"); */
  /*   fflush(stdout); */
  /* } */

  /*
   * Cleanup
   */

  PDM_part_to_block_free (ptb);
  PDM_block_to_part_free (btp);

  free (b_connectPointA_true);
  free (b_connectPointB_true);
  free (b_gNumA_true);
  free (b_gNumB_true);
  free (b_coordsPointA_true);
  free (b_coordsPointB_true);
  free (b_gNumEdgeA_true);
  free (b_gNumEdgeB_true);
  free (b_nNewPointsA_true);
  free (b_nNewPointsB_true);
  free (b_oNewPointsA_true);
  free (b_oNewPointsB_true);
  free (b_stride_one_idx_true);
  free (b_tIntersects_true);

  free (b_b_gNumVtxFromAForB);
  free (_b_b_stride_packB);
  free (b_b_stride_packB);
  free (b_b_cNewPointsB_true_pack);

  free (b_b_gNumVtxFromBForA);
  free (_b_b_stride_packA);
  free (b_b_cNewPointsA_true_pack);

  free (b_stridePtsADep_true);

  /*
   * Update edges intersection structure structures
   *
   *  keys : absolu ou local ?
   *  Comment relier au mieux la structure _edges_intersect_res_t* avec
   *  les tableau r_
   *
   */

  int idxData = 0;

  int *r_stride_one_idx_true     = malloc (sizeof(int) * (n_procData + 1));
  int *r_stridePtsADep_idx_true = malloc (sizeof(int) * (n_procData + 1));
  int *r_stridePtsBDep_idx_true = malloc (sizeof(int) * (n_procData + 1));

  r_stride_one_idx_true[0]     = 0;
  r_stridePtsADep_idx_true[0] = 0;
  r_stridePtsBDep_idx_true[0] = 0;

  for (int i = 0; i < n_procData; i++) {
    r_stride_one_idx_true[i+1] = r_stride_one_idx_true[i] + r_stride_one_true[0][i];
    r_stridePtsADep_idx_true[i+1] = r_stridePtsADep_idx_true[i] + r_stridePtsADep_true[0][i];
    r_stridePtsBDep_idx_true[i+1] = r_stridePtsBDep_idx_true[i] + r_stridePtsBDep_true[0][i];
  }

  for (int key = 0; key < keyMax; key++) {

    //  int keyloc

    _edges_intersect_res_t ** datas =
      (_edges_intersect_res_t **) PDM_hash_tab_data_get (ht,
                                                         (void *) &key);

    if (vb >= 3) printf("\n\n ++> key : %d\n", key);

    int nData = PDM_hash_tab_n_data_get (ht, &key);

    for (int i = 0; i < nData; i++) {

      _edges_intersect_res_t *_data = datas[i];

      int isFound = 0;

      int idxA_ = r_stridePtsADep_idx_true[idxData];
      int idxB_ = r_stridePtsBDep_idx_true[idxData];

      if (vb >= 3) PDM_printf ("\n      - idata: %d\n", i);
      if (vb >= 3) {
        PDM_printf ("       nGEdgeA :%d  nGEdgeB :%d  tIntersect : %d\n",
                    _data->nGEdgeA,
                    _data->nGEdgeB,
                    _data->tIntersect);}

      for (int k = r_stride_one_idx_true[idxData]; k < r_stride_one_idx_true[idxData+1]; k++) {

        if ((_data->nGEdgeA != r_gNumEdgeA_true[0][k]) ||
            (_data->nGEdgeB != r_gNumEdgeB_true[0][k])){
          idxA_ += r_nNewPointsA_true[0][k];
          idxB_ += r_nNewPointsB_true[0][k];
          continue;
        }

        isFound = 1;
        _data->tIntersect = (PDM_line_intersect_t) r_tIntersects_true[0][k];
        const int _nNewPointsA = r_nNewPointsA_true[0][k];

        if (_data->nNewPointsA != _nNewPointsA) {
          _data->nNewPointsA = _nNewPointsA;
          _data->uA          = realloc (_data->uA, sizeof(double) * _nNewPointsA);
          _data->coordsA     = realloc (_data->coordsA, sizeof(double) * 3 * _nNewPointsA);
          _data->linkA       = realloc (_data->linkA, sizeof(PDM_g_num_t) * _nNewPointsA);
          _data->gNumA       = realloc (_data->gNumA, sizeof(PDM_g_num_t) * _nNewPointsA);
          _data->oNewPointsA =
            realloc (_data->oNewPointsA, sizeof(PDM_edges_intersect_point_t) * _nNewPointsA);
        }

        _data->nNewPointsA = _nNewPointsA;
        if (vb >= 3) PDM_printf ("       nNewPointsA :%d\n",  _data->nNewPointsA);

        for (int k1 = 0; k1 < _nNewPointsA; k1++) {
          _data->oNewPointsA[k1] = r_oNewPointsA_true[0][idxA_+k1];
          _data->gNumA[k1]       = r_gNumA_true[0][idxA_+k1];
          _data->linkA[k1]       = r_connectPointA_true[0][idxA_+k1];
          _data->uA[k1]          = r_uPointA_true[0][idxA_+k1];
          for (int k2 = 0; k2 < 3; k2++) {
            _data->coordsA[3*k1+k2] = r_coordsPointA_true[0][3*(idxA_+k1)+k2];
          }
          if (vb >= 3) {
            PDM_printf ("        oNewPointsA :%d\n", _data->oNewPointsA[k1]);
            PDM_printf ("        cPointA : %12.5e  %12.5e %12.5e\n"
                        , _data->coordsA[3*k1]
                        , _data->coordsA[3*k1+1]
                        , _data->coordsA[3*k1+2]);
            PDM_printf ("        gNumA :"PDM_FMT_G_NUM"\n", _data->gNumA[k1]);
            PDM_printf ("        linkA :"PDM_FMT_G_NUM"\n", _data->linkA[k1]);
            PDM_printf ("        uPointA :%12.5e\n", _data->uA[k1]);
          }
        }
        const int _nNewPointsB = r_nNewPointsB_true[0][k];

        if (_data->nNewPointsB != _nNewPointsB) {
          _data->nNewPointsB = _nNewPointsB;
          _data->uB          = realloc (_data->uB, sizeof(double) * _nNewPointsB);
          _data->coordsB     = realloc (_data->coordsB, sizeof(double) * 3 * _nNewPointsB);
          _data->linkB       = realloc (_data->linkB, sizeof(PDM_g_num_t) * _nNewPointsB);
          _data->gNumB       = realloc (_data->gNumB, sizeof(PDM_g_num_t) * _nNewPointsB);
          _data->oNewPointsB =
            realloc (_data->oNewPointsB, sizeof(PDM_edges_intersect_point_t) * _nNewPointsB);
        }

        _data->nNewPointsB = _nNewPointsB;
        if (vb >= 3) PDM_printf ("       nNewPointsB :%d\n",  _data->nNewPointsB);

        for (int k1 = 0; k1 < _nNewPointsB; k1++) {
          _data->oNewPointsB[k1] = r_oNewPointsB_true[0][idxB_+k1];
          _data->gNumB[k1]       = r_gNumB_true[0][idxB_+k1];
          _data->linkB[k1]       = r_connectPointB_true[0][idxB_+k1];
          _data->uB[k1]          = r_uPointB_true[0][idxB_+k1];
          for (int k2 = 0; k2 < 3; k2++) {
            _data->coordsB[3*k1+k2] = r_coordsPointB_true[0][3*(idxB_+k1)+k2];
          }
          if (vb >= 3) {
            PDM_printf ("        oNewPointsB :%d\n", _data->oNewPointsB[k1]);
            PDM_printf ("        cPointB : %12.5e  %12.5e %12.5e\n"
                        , _data->coordsB[3*k1]
                        , _data->coordsB[3*k1+1]
                        , _data->coordsB[3*k1+2]);
            PDM_printf ("        gNumB :"PDM_FMT_G_NUM"\n", _data->gNumB[k1]);
            PDM_printf ("        linkB :"PDM_FMT_G_NUM"\n", _data->linkB[k1]);
            PDM_printf ("        uPointB :%12.5e\n", _data->uB[k1]);
          }
        }
        break;
      }

      if (!isFound) {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_edges_intersections :"
                         "No Data found to update current intersection\n");
        abort();
      }
      idxData += 1;
    }


  }

  /* PDM_MPI_Barrier (_ei->comm); */
  /* if (i_rank == 0) { */
  /*   printf("-- step 15\n"); */
  /*   fflush(stdout); */
  /* } */

  /*
   * Cleanup
   */

  free (r_stride_one_idx_true);
  free (r_stridePtsADep_idx_true);
  free (r_stridePtsBDep_idx_true);

  free (r_stride_one_true[0]);
  free (r_tIntersects_true[0]);
  free (r_gNumEdgeA_true[0]);
  free (r_gNumEdgeB_true[0]);

  free (r_stridePtsADep_true[0]);
  free (r_nNewPointsA_true[0]);
  free (r_oNewPointsA_true[0]);
  free (r_connectPointA_true[0]);
  free (r_gNumA_true[0]);
  free (r_gNumA_true);
  free (r_uPointA_true[0]);
  free (r_coordsPointA_true[0]);

  free (r_stridePtsBDep_true[0]);
  free (r_nNewPointsB_true[0]);
  free (r_oNewPointsB_true[0]);
  free (r_connectPointB_true[0]);
  free (r_gNumB_true[0]);
  free (r_gNumB_true);
  free (r_uPointB_true[0]);
  free (r_coordsPointB_true[0]);

  free (r_stride_one_true);
  free (r_tIntersects_true);
  free (r_gNumEdgeA_true);
  free (r_gNumEdgeB_true);

  free (r_stridePtsADep_true);
  free (r_nNewPointsA_true);
  free (r_oNewPointsA_true);
  free (r_connectPointA_true);
  free (r_uPointA_true);
  free (r_coordsPointA_true);

  free (r_stridePtsBDep_true);
  free (r_nNewPointsB_true);
  free (r_oNewPointsB_true);
  free (r_connectPointB_true);
  free (r_uPointB_true);
  free (r_coordsPointB_true);

  if (vb >= 1) PDM_printf ("==== PDM_edges_intersect_synchronize ==== terminated ====\n");

}

/**
 *
 * \brief dump
 *
 * \param [in]   ei           Current edges intersection pointer
 *
 */

void
PDM_edges_intersect_dump
(
PDM_edges_intersect_t       *ei
)
{

  _edges_intersect_t *_ei = (_edges_intersect_t *) ei;

  PDM_printf ("_ei->ht : %d\n", _ei->ht);
  PDM_printf ("_ei->ht->keyMax : %d\n", *((int *) PDM_hash_tab_keyMax_get(_ei->ht)));
  PDM_hash_tab_dump(_ei->ht);
  PDM_printf ("_ei->htA : %d\n", _ei->htA);
  PDM_printf ("_ei->htA->keyMax : %d\n", *((int *) PDM_hash_tab_keyMax_get(_ei->htA)));
  PDM_hash_tab_dump(_ei->htA);
  PDM_printf ("_ei->htB : %d\n", _ei->htB);
  PDM_printf ("_ei->htB->keyMax : %d\n", *((int *) PDM_hash_tab_keyMax_get(_ei->htB)));
  PDM_hash_tab_dump(_ei->htB);

}

void
PDM_edges_intersect_res_dump
(
PDM_edges_intersect_res_t       *eir
)
{


  _edges_intersect_res_t *_eir = (_edges_intersect_res_t *) eir;
	PDM_printf ("\n\n---\n--- Intersection - nGEdgeA : %d, nGEdgeB : %d begin\n",
              _eir->nGEdgeA, _eir->nGEdgeB );


  PDM_printf ("---  * Line intersection type : %s\n", _typeInter[_eir->tIntersect+1]);
	PDM_printf ("---  * Number of new points in mesh A : %d\n", _eir->nNewPointsA);

  for (int k1 = 0; k1 < _eir->nNewPointsA; k1++) {
    PDM_printf ("---     - Point             : %d\n", k1);
    PDM_printf ("---        Type             : %s\n", _typeInter2[_eir->oNewPointsA[k1]]);
    if ((_eir->oNewPointsA[k1] == 2) || (_eir->oNewPointsA[k1] == 3)) {
      PDM_printf ("---        merged vertex A  : "PDM_FMT_G_NUM"\n", _eir->gNumA[k1]);
    }
    if ((_eir->oNewPointsA[k1] == 2) || (_eir->oNewPointsA[k1] == 3)) {
      PDM_printf ("---        linked vertex B  : %d\n",  _eir->linkA[k1]);
    }
    //if ((_eir->oNewPointsA[k1] == 0) || (_eir->oNewPointsA[k1] == 2)) {
    PDM_printf ("---        Parametric coord from "PDM_FMT_G_NUM" : %12.5e\n",
                _eir->originEdgeA, _eir->uA[k1]);
      //}
    PDM_printf ("---        Coordinates      : %12.5e %12.5e %12.5e\n",
                _eir->coordsA[3*k1], _eir->coordsA[3*k1+1], _eir->coordsA[3*k1+2]);
  }

	PDM_printf ("---  * Number of new points in mesh B : %d\n", _eir->nNewPointsB);

  for (int k1 = 0; k1 < _eir->nNewPointsB; k1++) {
    PDM_printf ("---     - Point             : %d\n", k1);
    PDM_printf ("---        Type             : %s\n",  _typeInter2[_eir->oNewPointsB[k1]]);
    if ((_eir->oNewPointsB[k1] == 1) || (_eir->oNewPointsB[k1] == 3)) {
      PDM_printf ("---        merged vertex B  : "PDM_FMT_G_NUM"\n", _eir->gNumB[k1]);
    }
    if  ((_eir->oNewPointsB[k1] == 1) || (_eir->oNewPointsB[k1] == 3)) {
      PDM_printf ("---        linked vertex A  : %d\n",  _eir->linkB[k1]);
    }
    //if ((_eir->oNewPointsB[k1] == 0) || (_eir->oNewPointsB[k1] == 1)) {
    PDM_printf ("---        Parametric coord from "PDM_FMT_G_NUM" : %12.5e\n", _eir->originEdgeB, _eir->uB[k1]);
      //}
    PDM_printf ("---        Coordinates      : %12.5e %12.5e %12.5e\n",
                _eir->coordsB[3*k1], _eir->coordsB[3*k1+1], _eir->coordsB[3*k1+2]);
  }

	PDM_printf ("--- Intersection - nGEdgeA : %d, nGEdgeB : %d end\n\n", _eir->nGEdgeA, _eir->nGEdgeB );
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
