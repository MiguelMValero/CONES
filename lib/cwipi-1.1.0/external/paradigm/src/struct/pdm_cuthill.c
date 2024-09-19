/*============================================================================
 * Hilbert encoding for 2D or 3D coordinates.
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_cuthill.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define _MIN(a,b)   ((a) < (b) ?  (a) : (b))  /* Minimum of a et b */

#define _MAX(a,b)   ((a) > (b) ?  (a) : (b))  /* Maximum of a et b */

/*=============================================================================
 * Static global variables
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief Compute Bandwidth of a graph
 *
 * \param [in,out]  node_num            The number of nodes.
 * \param [in,out]  adj_row[ADJ_NUM]    Information about row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ
 * \param [in,out]  adj                 The adjacency structure. For each row, it contains the column indices of the nonzero entries.
 * \param [out]     ADJ_BANDWIDTH,      The bandwidth of the adjacency matrix.
 */

static int
_adj_bandwidth
(
int node_num,
int adj_row[],
int adj[]
)
{
  int band_hi;
  int band_lo;
  int i,j,col;
  int value;

  band_lo = 0;
  band_hi = 0;

  for ( i = 0; i < node_num; i++ )
  {
    for ( j = adj_row[i]; j <= adj_row[i+1]-1; j++ )
    {
      col = adj[j-1] - 1;
      band_lo = _MAX( band_lo, i - col );
      band_hi = _MAX( band_hi, col - i );
    }
  }

  value = band_lo + 1 + band_hi;

  return value;

}

/**
 *
 * \brief TODOUX
 *
 */

static void
_level_set
(
int  root,
int  adj_row[],
int  adj[],
int  mask[],
int *level_num,
int  level_row[],
int  level[]
)
{
  int i;
  int iccsze;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;

  mask[root-1] = 0;
  level[0] = root;
  *level_num = 0;
  lvlend = 0;
  iccsze = 1;
  /*
   *  LBEGIN is the pointer to the beginning of the current level, and
   *  LVLEND points to the end of this level.
   */
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = iccsze;
    *level_num = *level_num + 1;
    level_row[*level_num-1] = lbegin;
   /*
    *   Generate the next level by finding all the masked neighbors of nodes
    *   in the current level.
    */
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = level[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          iccsze = iccsze + 1;
          level[iccsze-1] = nbr;
          mask[nbr-1] = 0;
        }
      }
    }
    /*
     *   Compute the current level width (the number of nodes encountered.)
     *   If it is positive, generate the next level.
     */
    lvsize = iccsze - lvlend;

    if ( lvsize <= 0 )
    {
      break;
    }
  }

  level_row[*level_num] = lvlend + 1;
  /**  Reset MASK to 1 for the nodes in the level structure. **/
  for ( i = 0; i < iccsze; i++ )
  {
    mask[level[i]-1] = 1;
  }

  return;
}

/**
 *
 * \brief ROOT_FIND finds a pseudo-peripheral node.
 *
 */

static void
_root_find
(
int *root,
int adj_row[],
int adj[],
int mask[],
int *level_num,
int level_row[],
int level[]
)
{
  int iccsze;
  int j;
  int jstrt;
  int k;
  int kstop;
  int kstrt;
  int level_num2;
  int mindeg;
  int nabor;
  int ndeg;
  int node;

  /** Determine the level structure rooted at ROOT. **/

  _level_set ( *root, adj_row, adj, mask, level_num,
              level_row, level);

    /** Count the number of nodes in this level structure. **/
  iccsze = level_row[*level_num] - 1;

  /*  Extreme case:
   *    A complete graph has a level set of only a single level.
   *    Every node is equally good (or bad).
   */

  if ( *level_num == 1 )
  {
    return;
  }
  /*   Extreme case:
   *     A "line graph" 0--0--0--0--0 has every node in its only level.
   *     By chance, we've stumbled on the ideal root.
   */
  if ( *level_num == iccsze )
  {
    return;
  }
  /*
   *   Pick any node from the last level that has minimum degree
   *   as the starting point to generate a new level set.
   */
  for ( ; ; )
  {
    mindeg = iccsze;

    jstrt = level_row[*level_num-1];
    *root = level[jstrt-1];

    if ( jstrt < iccsze )
    {
      for ( j = jstrt; j <= iccsze; j++ )
      {
        node = level[j-1];
        ndeg = 0;
        kstrt = adj_row[node-1];
        kstop = adj_row[node] - 1;

        for ( k = kstrt; k <= kstop; k++ )
        {
          nabor = adj[k-1];
          if ( 0 < mask[nabor-1] )
          {
            ndeg = ndeg + 1;
          }
        }

        if ( ndeg < mindeg )
        {
          *root = node;
          mindeg = ndeg;
        }
      }
    }

   /** Generate the rooted level structure associated with this node. **/
   _level_set( *root, adj_row, adj, mask, &level_num2,
                level_row, level);

   /** If the number of levels did not increase, accept the new ROOT. **/

    if ( level_num2 <= *level_num )
    {
      break;
    }

    *level_num = level_num2;
    /*
     *   In the unlikely case that ROOT is one endpoint of a line graph,
     *   we can exit now.
     */
    if ( iccsze <= *level_num )
    {
      break;
    }
  }

  return;

}

/**
 *
 * \brief Reverse an integer array
   \param [in]     The number of entries in the array.
   \param [in,out] int A(N), the array to be reversed.
 *
 */

static void
_i4vec_reverse
(
int n,
int a[]
)
{
  int i;
  int j;

  for ( i = 0; i < n / 2; i++ )
  {
    j        = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = j;
  }

  return;
}

/**
 *
 * \brief TODOUX
 *
 */

static void
_degree
(
int root,
int adj_row[],
int adj[],
int mask[],
int deg[],
int *iccsze,
int ls[]
)
{
  int i;
  int ideg;
  int j;
  int jstop;
  int jstrt;
  int lbegin;
  int lvlend;
  int lvsize;
  int nbr;
  int node;

  /** The sign of ADJ_ROW(I) is used to indicate if node I has been considered. **/

  ls[0] = root;
  adj_row[root-1] = -adj_row[root-1];
  lvlend = 0;
  *iccsze = 1;
  /*
   *   LBEGIN is the pointer to the beginning of the current level, and
   *   LVLEND points to the end of this level.
   */
  for ( ; ; )
  {
    lbegin = lvlend + 1;
    lvlend = *iccsze;
    /*
     *   Find the degrees of nodes in the current level,
     *   and at the same time, generate the next level.
     */
    for ( i = lbegin; i <= lvlend; i++ )
    {
      node = ls[i-1];
      jstrt = -adj_row[node-1];
      jstop = abs ( adj_row[node] ) - 1;
      ideg = 0;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          ideg = ideg + 1;

          if ( 0 <= adj_row[nbr-1] )
          {
            adj_row[nbr-1] = -adj_row[nbr-1];
            *iccsze = *iccsze + 1;
            ls[*iccsze-1] = nbr;
          }
        }
      }
      deg[node-1] = ideg;
    }

    /** Compute the current level width. **/

    lvsize = *iccsze - lvlend;

    /**  If the current level width is nonzero, generate another level. **/

    if ( lvsize == 0 )
    {
      break;
    }
  }

  /** Reset ADJ_ROW to its correct sign and return. **/

  for ( i = 0; i < *iccsze; i++ )
  {
    node = ls[i] - 1;
    adj_row[node] = -adj_row[node];
  }

  return;
}

/**
 *
 * \brief Compute reverse Cut-Hill Mac-Kee ordering
 *
 * \param [in,out]  node_num            The number of nodes.
 * \param [in,out]  adj_num[NODE_NUM+1] The number of adjacency entries
 * \param [in,out]  adj_row[ADJ_NUM]    Information about row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ
 * \param [in,out]  adj                 The adjacency structure. For each row, it contains the column indices of the nonzero entries.
 * \param [out]     perm                The RCM ordering
 */

static void
_rcm
(
int root,
int adj_row[],
int adj[],
int mask[],
int perm[],
int *iccsze,
int node_num
)
{

  int fnbr;
  int i;
  int j;
  int jstop;
  int jstrt;
  int k;
  int l;
  int lbegin;
  int lnbr;
  int lperm;
  int lvlend;
  int nbr;
  int node;

  /** If node_num out of bounds, something is wrong. **/

  if ( node_num < 1 )
  {
    PDM_printf( "\n");
    PDM_printf( "RCM - Fatal error!\n");
    PDM_printf( "  Unacceptable input value of NODE_NUM = %d \n ", node_num);
    exit ( 1 );
  }

  /** If the root is out of bounds, something is wrong. **/

  if ( root < 1 || node_num < root )
  {
    PDM_printf("\n");
    PDM_printf("RCM - Fatal error!\n");
    PDM_printf("  Unacceptable input value of ROOT = %d \n ",root);
    PDM_printf("  Acceptable values are between 1 and %d inclusive \n",node_num);
    exit ( 1 );
  }

  /** Allocate memory for the degree array. **/
  int *deg = (int *) malloc( sizeof(int) * node_num); // deg = new int[node_num];

  /** Find the degrees of the nodes in the component specified by MASK and ROOT. **/

  _degree ( root, adj_row, adj, mask, deg, iccsze, perm);

  /** If the connected component size is less than 1, something is wrong. **/

  if ( *iccsze < 1 )
  {
    PDM_printf("\n");
    PDM_printf("RCM - Fatal error!\n");
    PDM_printf("  Connected component size ICCSZE returned from DEGREE as %d \n ", *iccsze);
    exit ( 1 );
  }

  /** Set the mask value for the root. **/

  mask[root-1] = 0;

  /** If the connected component is a singleton, there is no ordering necessary. **/

  if ( *iccsze == 1 )
  {
    free(deg);// delete [] deg;
    return;
  }
  /*   Carry out the reordering.
   *
   *   LBEGIN and LVLEND point to the beginning and
   *  the end of the current level respectively.
   */
  lvlend = 0;
  lnbr = 1;

  while ( lvlend < lnbr )
  {
    lbegin = lvlend + 1;
    lvlend = lnbr;

    for ( i = lbegin; i <= lvlend; i++ )
    {

      /** For each node in the current level...**/

      node = perm[i-1];
      jstrt = adj_row[node-1];
      jstop = adj_row[node] - 1;

      /* Find the unnumbered neighbors of NODE.
       *   FNBR and LNBR point to the first and last neighbors
       *   of the current node in PERM.
       */
      fnbr = lnbr + 1;

      for ( j = jstrt; j <= jstop; j++ )
      {
        nbr = adj[j-1];

        if ( mask[nbr-1] != 0 )
        {
          lnbr = lnbr + 1;
          mask[nbr-1] = 0;
          perm[lnbr-1] = nbr;
        }
      }

      /** If no neighbors, skip to next node in this level. **/

      if ( lnbr <= fnbr )
      {
        continue;
      }
     /*
      *   Sort the neighbors of NODE in increasing order by degree.
      *   Linear insertion is used.
      */
      k = fnbr;

      while ( k < lnbr )
      {
        l = k;
        k = k + 1;
        nbr = perm[k-1];

        while ( fnbr < l )
        {
          lperm = perm[l-1];

          if ( deg[lperm-1] <= deg[nbr-1] )
          {
            break;
          }

          perm[l] = lperm;
          l = l - 1;
        }
        perm[l] = nbr;
      }
    }
  }

  /*
  *   We now have the Cuthill-McKee ordering.
  *   Reverse it to get the Reverse Cuthill-McKee ordering.
  */
  _i4vec_reverse ( *iccsze, perm );

  /**  Free memory. **/
  free(deg);// delete [] deg;


  return;
}


/**
 *
 * \brief Compute reverse Cut-Hill Mac-Kee ordering
 *
 * \param [in,out]  node_num            The number of nodes.
 * \param [in,out]  adj_num[NODE_NUM+1] The number of adjacency entries
 * \param [in,out]  adj_row[ADJ_NUM]    Information about row I is stored in entries ADJ_ROW(I) through ADJ_ROW(I+1)-1 of ADJ
 * \param [in,out]  adj                 The adjacency structure. For each row, it contains the column indices of the nonzero entries.
 * \param [out]     perm                The RCM ordering
 */

void
PDM_genrcm
(
int node_num,
int adj_row[],
int adj[],
int perm[]
)
{
  int i;
  int iccsze;
  int level_num;
  int num;
  int root;

  int *level_row  = (int *) malloc(sizeof(int) * (node_num + 1)); //level_row = new int[node_num+1];
  int *mask       = PDM_array_const_int(node_num, 1); //mask = new int[node_num];

  num = 1;

  for ( i = 0; i < node_num; i++ )
  {
    /*
     * For each masked connected component...
     */
    // printf(" mask[%i] = %i \n", i, mask[i]);
    if ( mask[i] != 0 )
    {
      root = i + 1;
      /*
       *  Find a pseudo-peripheral node ROOT.  The level structure found by
       *   ROOT_FIND is stored starting at PERM(NUM).
       */
      _root_find(&root, adj_row, adj, mask, &level_num,
                level_row, perm+num-1);
      /*
       *   RCM orders the component using ROOT as the starting node.
       */
      _rcm(root, adj_row, adj, mask, perm+num-1, &iccsze, node_num );

      num = num + iccsze;

      /*
       *   We can stop once every node is in one of the connected components.
       */
      if ( node_num < num )
      {
        free(level_row); //delete [] level_row;
        free(mask);      //delete [] mask;
        return;
      }
    }
  }

  free(level_row); //delete [] level_row;
  free(mask);      //delete [] mask;

  return;
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

/**
 *
 * \brief Compute bandwidth of a mesh partition
 *
 * \param [in] n_elm          Number of elements to reorder
 * \param [in] dual_graph_idx Element to element connectivity indexes (size=n_elm+1)
 * \param [in] dual_graph     Element to element connectivity (size=dual_graph_idx[n_elm])
 */

int
PDM_cuthill_checkbandwidth
(
 int                n_elm,
 int               *dual_graph_idx,
 int               *dual_graph
)
{

  /** Do a copy since graph seems to be modified (?) **/
  int *dual_graph_idx_tmp = (int *) malloc((n_elm + 1) * sizeof(int));
  int *dual_graph_tmp     = (int *) malloc(dual_graph_idx[n_elm] * sizeof(int));

  /** Offset Graph and Arr **/
  for (int i = 0; i < n_elm; i++){
    dual_graph_idx_tmp[i] = dual_graph_idx[i]+1;
  }
  dual_graph_idx_tmp[n_elm] = dual_graph_idx[n_elm];
  for (int i = 0; i < dual_graph_idx[n_elm]; i++){
    dual_graph_tmp[i] = dual_graph[i]+1;
  }

  int dualBandWidth = _adj_bandwidth(dual_graph_idx_tmp[n_elm], dual_graph_idx_tmp, dual_graph_tmp);

  free(dual_graph_idx_tmp);
  free(dual_graph_tmp);
  return dualBandWidth;
}

/**
 *
 * \brief Compute reverse CutHill-MecKee reordering
 *
 * \param [in] n_elm          Number of elements to reorder
 * \param [in] dual_graph_idx Element to element connectivity indexes (size=n_elm+1)
 * \param [in] dual_graph     Element to element connectivity (size=dual_graph_idx[n_elm])
 *
 * \param [out] perm          Array of permutations
 */

void
PDM_cuthill_generate
(
 int                n_elm,
 int               *dual_graph_idx,
 int               *dual_graph,
 int               *perm
)
{

  /** Do a copy since graph seems to be modified (?) **/
  int *dual_graph_idx_tmp = (int *) malloc((n_elm + 1) * sizeof(int));
  int *dual_graph_tmp     = (int *) malloc(dual_graph_idx[n_elm] * sizeof(int));

  /** Offset Graph and Arr **/
  for (int i = 0; i < n_elm; i++){
    dual_graph_idx_tmp[i] = dual_graph_idx[i]+1;
  }
  dual_graph_idx_tmp[n_elm] = dual_graph_idx[n_elm];
  for (int i = 0; i < dual_graph_idx[n_elm]; i++){
    dual_graph_tmp[i] = dual_graph[i]+1;
  }

  /** Apply rcm to current Graph **/
  PDM_genrcm(n_elm, dual_graph_idx_tmp, dual_graph_tmp, perm);

  /** Offset Permutation and Graph arrays **/
  for (int i = 0; i < n_elm; i++)
  {
    perm[i] = perm[i]-1;
  }

  /** Verbose **/
  if (0 == 1) {
      PDM_printf("\n Contenu de perm : \n");
      for(int i = 0; i < n_elm; i++) {
        PDM_printf(" %d ", perm[i]);
    }
    PDM_printf("\n");
  }
  /** Free **/
  free(dual_graph_idx_tmp);
  free(dual_graph_tmp);
}

#ifdef  __cplusplus
}
#endif
