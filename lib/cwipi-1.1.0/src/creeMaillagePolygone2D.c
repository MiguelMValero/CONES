/*
  This file is part of the CWIPI library. 

  Copyright (C) 2011  ONERA

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
#include <mpi.h>

#include <stdio.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "creeMaillagePolygone2D.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
  Attention :
  -----------
  Les fichiers inclus par <metis.h> posent problème à minima sous Linux du à
  des redéclarations de fonctions de <stdlib.h>. On recopie dans ce cas le
  minimum nécessaire issu de <metis.h> (correspondant à METIS 4.0)
*/

#ifdef HAVE_METIS
typedef int idxtype;
void METIS_PartGraphRecursive(int *, idxtype *, idxtype *, idxtype *, idxtype *,
                              int *, int *, int *, int *, int *, idxtype *);
void METIS_PartGraphKway(int *, idxtype *, idxtype *, idxtype *, idxtype *,
                         int *, int *, int *, int *, int *, idxtype *);
#endif

#ifdef __cplusplus
}
#endif

#define ABS(a)     ((a) <  0  ? -(a) : (a))
#define MIN(a,b)   ((a) > (b) ?  (b) : (a))

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


static double random01(void)
{
  int sign;
  int rsigna = rand();
  int rsignb = rand();
  sign = (rsigna - rsignb) / ABS(rsigna - rsignb);
  double resultat =   sign*((double)rand())/((double)RAND_MAX);
  return resultat;
}

void creeMaillagePolygone2D(int order,
                            MPI_Comm localComm,
                            double xmin,
                            double xmax,
                            double ymin,
                            double ymax,
                            int initRandom,
                            int nx,
                            int ny,
                            int *nVertex,
                            double **meshCoords,
                            int *nElts,
                            int **eltsConnecPointer,
                            int **eltsConnec)
{

  int nRank;
  MPI_Comm_size(localComm, &nRank);

  int localRank;
  MPI_Comm_rank(localComm, &localRank);

  int *globalVertexNum = NULL;
  int *globalEltNum = NULL;

  if (localRank == 0) {
    const double coefRand = 0.3;
    //const double coefRand = .3;
    //const double coefRand = 0.45;
    srand(initRandom);

    /* Construction des coordonnes */
    /* --------------------------- */

    /* nx et ny doivent etre pair et superieur a 4 */

    int nx1 = nx;
    int ny1 = ny;

    if (nx1 < 4)
      nx1 = 4;

    if (ny1 < 4)
      ny1 = 4;

    if (nx1 % 2 != 0)
      nx1 += 1;

    if (ny1 % 2 != 0)
      ny1 += 1;

    /* cotes d un octogone */

    double cote1;
    if (nx1 == 4)
      cote1 = (xmax - xmin) / (sqrt(2) + 1);
    else
      cote1 = (xmax - xmin) / ((nx1 - 2)/2 + ((nx1 - 2)/2 - 1) * sqrt(2) + sqrt(2));
    double cote2;
    if (ny1 == 4)
      cote2 = (ymax - ymin) / (sqrt(2) + 1);
    else
      cote2 = (ymax - ymin) / ((ny1 - 2)/2 + ((ny1 - 2)/2 - 1) * sqrt(2) + sqrt(2));

    /* Sommets */

    *nVertex = 0;

    /* Allocation temporaire (surdimensionnee) */

    *meshCoords = (double*) malloc (sizeof(double) * 3 * (nx1 * ny1));

    int cpt = 0;
    int cptMax;
    if (ny1 > 4)
      cptMax = ny1 + (ny1-4)/2;
    else
      cptMax = ny1;

    double ycourant = ymin;
    const double eps = 1e-5;
    while(cpt < cptMax) {
      if (cpt % 3 == 1 || cpt % 3 == 2) {
        for (int ix = 0; ix < nx1/2; ix++) {
          (*meshCoords)[3*(*nVertex)]   = xmin + ix * (1+sqrt(2)) * cote1;
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }
      }
      else if ((cpt % 3) == 0) {
        if ((cpt == (cptMax-1)) || (cpt == 0)) {
          (*meshCoords)[3*(*nVertex)]   = xmin;
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }

        (*meshCoords)[3*(*nVertex)]   = xmin + cote1/sqrt(2);
        (*meshCoords)[3*(*nVertex)+1] = ycourant;
        (*meshCoords)[3*(*nVertex)+2] = 0.;
        (*nVertex)++;

        for (int ix = 2; ix < nx1-1; ix++) {
          if (ix % 2 == 0)
            (*meshCoords)[3*(*nVertex)] = (*meshCoords)[3*((*nVertex)-1)]+cote1;
          else
            (*meshCoords)[3*(*nVertex)] = (*meshCoords)[3*((*nVertex)-1)]+cote1*sqrt(2);
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }
        if ((cpt == (cptMax-1)) || (cpt == 0)) {
          (*meshCoords)[3*(*nVertex)]   = xmax;
          (*meshCoords)[3*(*nVertex)+1] = ycourant;
          (*meshCoords)[3*(*nVertex)+2] = 0.;
          (*nVertex)++;
        }
      }
      cpt++;
      if ((cpt % 3 == 1) || (cpt % 3 == 0))
        ycourant += cote2/sqrt(2);
      else
        ycourant += cote2;
    }

    *meshCoords = (double *) realloc(*meshCoords,3*(*nVertex) * sizeof(double)); 
    for (int ix = 0; ix <(*nVertex) ; ix++) {
      if (ABS(xmin-(*meshCoords)[3*ix]) > eps &&
        ABS(xmax-(*meshCoords)[3*ix]) > eps &&
          ABS(ymax-(*meshCoords)[3*ix+1]) > eps &&
          ABS(ymin-(*meshCoords)[3*ix+1]) > eps) {
        (*meshCoords)[3*ix] += random01() * coefRand * cote1;
        (*meshCoords)[3*ix+1] += random01() * coefRand * cote2;

      }
      //(*meshCoords)[3*ix+2]=2*sin(3*(*meshCoords)[3*ix+1])*sin(3*(*meshCoords)[3*ix]);
      //(*meshCoords)[3*ix+2]=20*sin((*meshCoords)[3*ix+1]/5.)*sin((*meshCoords)[3*ix]/5.);
      (*meshCoords)[3*ix+2]=0.;
    }


    /* Construction des blocs d'elements */
    /* --------------------------------- */

    *nElts = 0;

    /* Allocation avec surestimation de la taille */
    /* Optimisation : Peut-etre reajuste */

    *eltsConnecPointer = (int*) malloc(sizeof(int) * (nx1 * ny1 + 1));

    *eltsConnec = (int*) malloc(sizeof(int) * 8 * nx1 * ny1);

    (*eltsConnecPointer)[0] = 0;

    for(int itype = 0; itype < 3; itype++) {

      int itype1 = itype;

      if (order != 1)
        itype1 = 2 - itype;

      int n1;
      int n2;
      int n3;

      if (itype1 == 0) {
        /* Triangles */

        /* -- Premiere ligne */
        n1 = 1;
        n2 = n1+nx1;

        for (int ix = 0; ix < nx1/2; ix++) {
          int ix1 = 2 * ix;
          int ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n1+ix1;
          (*eltsConnec)[ideb+1] = n1+ix1+1;
          (*eltsConnec)[ideb+2] = n2+ix;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;
        }

        /* -- Autres triangles (un a gauche un a droite */
        int nbLi = (ny1-4)/2;
        n1   = 1 + nx1 + nx1/2;
        n2   = n1 + nx1/2;
        n3   = n2 + nx1-2;

        for (int itri = 0; itri < nbLi; itri++) {
          int n4 = n1 + nx1/2 - 1;
          int n5 = n2 + nx1-2 - 1;
          int n6 = n3 + nx1/2 - 1;
          int ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n1;
          (*eltsConnec)[ideb+1] = n2;
          (*eltsConnec)[ideb+2] = n3;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;

          ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n4;
          (*eltsConnec)[ideb+1] = n6;
          (*eltsConnec)[ideb+2] = n5;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;
          n1 = n3 + nx1/2;
          n2 = n1 + nx1/2;
          n3 = n2 + nx1-2;
        }

        /* -- Derniere ligne */
        n2 = n1+nx1/2;

        for (int ix = 0; ix < nx1/2; ix++) {
          int ix1 = 2 * ix;
          int ideb = (*eltsConnecPointer)[*nElts] ;
          (*eltsConnec)[ideb]   = n1+ix;
          (*eltsConnec)[ideb+1] = n2+ix1+1;
          (*eltsConnec)[ideb+2] = n2+ix1;
          *nElts += 1;
          (*eltsConnecPointer)[*nElts] = ideb + 3;
        }
      }

      else if (itype1 == 1) {

        /* Quadrangles */
        int nxQuad = (nx1-4)/2;
        int nyQuad = (ny1-4)/2;

        for (int iy = 0; iy < nyQuad; iy++) {
          for (int ix = 0; ix < nxQuad; ix++) {
            n1 = iy*(2*nx1-2) + nx1 + nx1/2 + 1 + ix + 1;
            n2 = iy*(2*nx1-2) + 2*nx1 + 1 + 2*ix + 1;
            n3 = n2 + 1;
            int n4 = iy*(2*nx1-2) + 3*nx1 - 2 + 1 + ix + 1 ;
            int ideb = (*eltsConnecPointer)[*nElts] ;
            (*eltsConnec)[ideb]   = n1;
            (*eltsConnec)[ideb+1] = n3;
            (*eltsConnec)[ideb+2] = n4;
            (*eltsConnec)[ideb+3] = n2;
            *nElts += 1;
            (*eltsConnecPointer)[*nElts] = ideb + 4;
          }
        }
      }

      else if (itype1 == 2) {

        /* Polygones */
        int nxPoly = (nx1-2)/2;
        int nyPoly = (ny1-2)/2;

        int delta = 0;

        for (int iy = 0; iy < nyPoly; iy++) {

          if (iy == 1)
            delta += 2*nx1;
          else if (iy !=0)
            delta += 2*nx1-2;
          for (int ix = 0; ix < nxPoly; ix++) {
            if (iy == 0)
              n1 = delta + 1 + 2*ix +1;
            else
              n1 = delta + 1 + 2*ix ;
            n2 = n1 + 1;
            n3 = iy*(2*nx1-2) + 1 + nx1 + ix;
            int n4 = n3 + 1;
            int n5 = iy*(2*nx1-2) + 1 + nx1 + nx1/2 + ix;
            int n6 = n5 + 1;
            int n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix;
            if (iy == (nyPoly - 1))
              n7 = iy*(2*nx1-2) + 1 + 2*nx1 + 2*ix + 1;
            int n8 = n7 + 1;
            int ideb = (*eltsConnecPointer)[*nElts] ;
            (*eltsConnec)[ideb]   = n1;
            (*eltsConnec)[ideb+1] = n2;
            (*eltsConnec)[ideb+2] = n4;
            (*eltsConnec)[ideb+3] = n6;
            (*eltsConnec)[ideb+4] = n8;
            (*eltsConnec)[ideb+5] = n7;
            (*eltsConnec)[ideb+6] = n5;
            (*eltsConnec)[ideb+7] = n3;
            *nElts += 1;
            (*eltsConnecPointer)[*nElts] = ideb + 8;
          }
        }
      }
    }
  }

  /* Decoupage du maillage par Metis */
  /* ------------------------------- */

  if (nRank > 1) {

#ifdef HAVE_METIS
    int *neighbourPointer = NULL;
    int *neighbour = NULL;
    int *downConnectivity = NULL;
    int *edges = NULL;
    int *edgeToFace = NULL;

    if (localRank == 0) {

      int nEdges = 0;
      int localNVertex = 0;
      int localNElts = 0;
      double *localCoords = NULL;

      /* Connectivite descendante */
      /* ------------------------ */

      /* Construction de la connectivite descendante et de la table des voisins
         A sortir dans une fonction propore*/

      /* - Construction des aretes + connectivite descendante */
      /* Numerotation de 1 a n */
      /* tri des aretes par le min des sommets */

      int *aretes = (int*) malloc(sizeof(int) * (*eltsConnecPointer)[*nElts]*2);

      downConnectivity =  (int) malloc(sizeof(int) * (*eltsConnecPointer)[*nElts]);

      int **triAre = (int**) malloc(sizeof(int*) * (*nVertex));

      int nDefaultAreVertex = 8;
      int *nAreVertex = (int*) malloc(sizeof(int) * (*nVertex));

      for (int i = 0; i < *nVertex; i++) {
        nAreVertex[i] = nDefaultAreVertex;
        triAre[i] = (int*) malloc(sizeof(int) * nAreVertex[i]);
        for (int j = 0; j < nAreVertex[i]; j++) {
          triAre[i][j] = -1;
        }
      }

      /* - Construction des aretes - Boucle sur les elements  */

      int iare = 0;
      for (int i = 0; i < *nElts; i++) {
        for (int j = (*eltsConnecPointer)[i]; j < (*eltsConnecPointer)[i+1]; j++) {
          aretes[2*iare]  = (*eltsConnec)[j];
          if (j != ((*eltsConnecPointer)[i+1] - 1))
            aretes[2*iare+1] = (*eltsConnec)[j+1];
          else
            aretes[2*iare+1] = (*eltsConnec)[(*eltsConnecPointer)[i]];
          downConnectivity[j] = iare+1;

          int minVertex = MIN(aretes[2*iare], aretes[2*iare+1]) - 1;
          int k = 1;
          while ((k < nAreVertex[minVertex]) && (triAre[minVertex][k-1] > -1))
            k++;
          if (k ==  nAreVertex[minVertex]){
            nAreVertex[minVertex] *= 2;
            triAre[minVertex] = (int*) realloc(triAre[minVertex], nAreVertex[minVertex] * sizeof(int));
            triAre[minVertex][k-1] = iare+1;
            for(int k1 = k ; k1 <nFacVertex[minVertex];k1++ )
              triFac[minVertex][k1] = -1;
          }
          else if (triAre[minVertex][k-1] == -1)
            triAre[minVertex][k-1] = iare+1;
          iare++;
        }
      }

      /* - Elimination des doublons - Boucle sur les sommets  */

      int *ptAretesCompactees = (int*) malloc(sizeof(int) * (*eltsConnecPointer)[*nElts]);

      for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++)
        ptAretesCompactees[i] = -1;

      int iareCompactee = 1;
      for (int i = 0; i < *nVertex; i++) {
        int j = 0;
        while(triAre[i][j] > -1) {
          int k = j+1;
          int iare1 = triAre[i][j] - 1;
          if (ptAretesCompactees[iare1] == -1) {
            ptAretesCompactees[iare1] = iareCompactee++;
            if (k < nAreVertex[i]) {
              while(triAre[i][k] > - 1) {
                int iare2 = triAre[i][k] - 1;
                if (ptAretesCompactees[iare2] == -1)
                  if (((aretes[2*iare1] == aretes[2*iare2]) &&
                       (aretes[2*iare1+1] == aretes[2*iare2+1])) ||
                      ((aretes[2*iare1] == aretes[2*iare2+1]) &&
                       (aretes[2*iare1+1] == aretes[2*iare2])))
                    ptAretesCompactees[iare2] = ptAretesCompactees[iare1];
                k += 1;
              }
            }
          }
          j +=1;
        }
        if (triAre[i] != NULL)
          free(triAre[i]);
      }

      nEdges = iareCompactee-1;
      if (triAre != NULL)
        free(triAre);

      if (nAreVertex != NULL)
        free(nAreVertex);


      /* - Renumerotation connectivite descendante - Boucle sur la connectivite des éléments */

      for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++)
        downConnectivity[i] = ptAretesCompactees[downConnectivity[i]-1];


      /* - Compression du tableau de description des aretes - Boucle sur la connectivitédes éléments */

      //edges = NULL;
      //BFT_MALLOC(edges, nEdges*2, int) ;

      //for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++) {
      //  edges[2*(ptAretesCompactees[i]-1)]   = aretes[2*i];
      //  edges[2*(ptAretesCompactees[i]-1)+1] = aretes[2*i+1];
      // }

      if (aretes != NULL)
        free(aretes);

      if (ptAretesCompactees != NULL)
        free(ptAretesCompactees);

      /* - Rangement des elements par couple suivant l'arete commune */

      edgeToFace = (int*) malloc(sizeof(int) * nEdges*2);  

      for (int i = 0; i < nEdges*2; i++)
        edgeToFace[i] = -1;

      for (int i = 0; i < *nElts; i++) {
        for (int j = (*eltsConnecPointer)[i]; j < (*eltsConnecPointer)[i+1]; j++) {
          if (edgeToFace[2*(downConnectivity[j]-1)] == -1 )
            edgeToFace[2*(downConnectivity[j]-1)] = i;
          else if (edgeToFace[2*(downConnectivity[j]-1)+1] == -1 )
            edgeToFace[2*(downConnectivity[j]-1)+1] = i;
          else {
            printf("Arete a plus de 2 facettes !!!!\n");
            exit(1);
          }
        }
      }

      if (downConnectivity != NULL)
        free(downConnectivity);

      /* - Creation de la table des voisins (Numerotation de 1 a n)
         Le voisin d'une face de bord est temporairement marque a -1 */

      int *tmpVoisins = (int*) malloc(sizeof(int) *  (*eltsConnecPointer)[*nElts]);

      for (int i = 0; i < (*eltsConnecPointer)[*nElts]; i++)
        tmpVoisins[i] = -2;

      for (int i = 0; i < nEdges; i++) {
        int elt1 = edgeToFace[2*i];
        int elt2 = edgeToFace[2*i+1];
        if (elt1 != -1) {
          int j = 0;
          while(tmpVoisins[(*eltsConnecPointer)[elt1] + j++] != -2);
          tmpVoisins[(*eltsConnecPointer)[elt1] + --j] = elt2;
          }
        if (elt2 != -1) {
          int j = 0;
          while(tmpVoisins[(*eltsConnecPointer)[elt2] + j++] != -2);
          tmpVoisins[(*eltsConnecPointer)[elt2] + --j] = elt1;
        }
      }

      free(edgeToFace);

      /* - Filtrage des faces de bords */

      neighbour = (int) malloc(sizeof(int) * (*eltsConnecPointer)[*nElts]);;

      neighbourPointer = (int) malloc(sizeof(int) *  (*nElts + 1));

      int nVoisin = 0;
      neighbourPointer[0] = nVoisin;
      for (int i = 0; i < *nElts; i++) {
        for (int j = (*eltsConnecPointer)[i]; j < (*eltsConnecPointer)[i+1]; j++) {
          if (tmpVoisins[j] != -1)
            neighbour[nVoisin++] = tmpVoisins[j];
        }
        neighbourPointer[i+1] = nVoisin;
      }
      
      neighbour = (int*) realloc((void *) neighbour, nVoisin * sizeof(int));

      if (tmpVoisins != NULL)
        free(tmpVoisins);

      /* Decoupage du maillage par Metis */
      /* ------------------------------- */

      int     wgtflag    = 0 ; /* Pas de pondération pour les faces ou cellules */
      int     numflag    = 0 ; /* Numérotation de 1 à n (type Fortran) */
      int     options[5] = {0, 3, 1, 1, 0} ; /* Par défaut si options[0] = 0 */
      int     edgecut    = 0 ; /* <-- nombre de faces sur la partition */

      int *numDomElt = (int*) malloc(sizeof(int) * (*nElts));

      assert(sizeof(idxtype) == sizeof(int));
      if (nRank < 8)

        METIS_PartGraphRecursive(nElts,
                                 (idxtype*)neighbourPointer,
                                 (idxtype*)neighbour,
                                 NULL,       /* vwgt   : poids des cellules */
                                 NULL,       /* adjwgt : poids des faces    */
                                 &wgtflag,
                                 &numflag,
                                 &nRank,
                                 options,
                                 &edgecut,
                                 numDomElt) ;

      else

        METIS_PartGraphKway(nElts,
                            (idxtype*)neighbourPointer,
                            (idxtype*)neighbour,
                            NULL,       /* vwgt   : poids des cellules */
                            NULL,       /* adjwgt : poids des faces    */
                            &wgtflag,
                            &numflag,
                            &nRank,
                            options,
                            &edgecut,
                            numDomElt) ;

      if (neighbour != NULL)
        free(neighbour);

      if (neighbourPointer != NULL)
        free(neighbourPointer);

      int  localEltsSize              = (*nElts)/nRank;   /* Estimation du nombre d'elements locaux
                                                             pour allocation memoire */
      int  localEltsConnecSize        = 6 * (*nElts)/nRank; /* Estimation de la taille
                                                             de la connectivite locale
                                                             pour allocation memoire */


      int *localEltsConnecPointer = NULL;

      localEltsConnecPointer = (int *) malloc(sizeof(int) * (localEltsSize+1));
      globalEltNum = (int *) malloc(sizeof(int) * localEltsSize);
      
      int *localEltsConnec = (int *) malloc(sizeof(int) * localEltsConnecSize);
      int *globalVertexNum = (int *) malloc(sizeof(int) * localEltsConnecSize);

      /* Pour chaque proc construction du maillage local envoi des donnees
         On finit par le proc 0 */

      for (int irank = nRank-1; irank >= 0; irank--) {
        int idom = irank;
        localNElts = 0;
        for (int ielt = 0; ielt < (*nElts); ielt++) {
          if (numDomElt[ielt] == idom) {
            if (localEltsSize <= localNElts) {
              localEltsSize *= 2;
              localEltsConnecPointer = (int *) realloc((void *) localEltsConnecPointer, 
                                                       sizeof(int) * (localEltsSize+1)); 
              globalEltNum = (int *) realloc((void *) globalEltNum, 
                                                       sizeof(int) * localEltsSize); 
            }
            globalEltNum[localNElts++] = ielt+1;
          }
        }

        int tmpSize = 0;
        localEltsConnecPointer[0] = 0;
        for (int ielt = 0; ielt < localNElts; ielt++) {
          for (int i = (*eltsConnecPointer)[globalEltNum[ielt]-1]; i < (*eltsConnecPointer)[globalEltNum[ielt]]; i++) {
            if (localEltsConnecSize <= tmpSize) {
              localEltsConnecSize *= 2;
              localEltsConnec = (int *) realloc((void *) localEltsConnec, 
                                                       sizeof(int) * localEltsConnecSize); 
              globalVertexNum = (int *) realloc((void *) globalVertexNum, 
                                                       sizeof(int) * localEltsConnecSize); 
            }
            globalVertexNum[tmpSize]   = (*eltsConnec)[i];
            localEltsConnec[tmpSize++] = (*eltsConnec)[i];
          }
          localEltsConnecPointer[ielt+1] = tmpSize;
        }

        /* Tri et elimination des doublons */

        localNVertex = 0;
        _quickSort(globalVertexNum, 0, tmpSize-1);
        int i = 0;
        int val = globalVertexNum[i];
        localNVertex  = 1;

        do {
          i++;
          while (( i < tmpSize && (val == globalVertexNum[i])))
            i++;
          if (i < tmpSize) {
            val = globalVertexNum[i];
            globalVertexNum[localNVertex++] = val;
          }
        } while (i < tmpSize);

        /* Renumerotation de la connectivite */

        int *tmpRenum = (int *) malloc (sizeof(int) * val);

        for (int i = 0; i < val; i++)
          tmpRenum[i] = -1;

        for (int i = 0; i < localNVertex; i++)
          tmpRenum[globalVertexNum[i]-1] = i;

        for (int i = 0; i < localEltsConnecPointer[localNElts]; i++)
          localEltsConnec[i] = tmpRenum[localEltsConnec[i]-1]+1;

        if (tmpRenum != NULL)
          free(tmpRenum);

        /* Coordonnees des sommets locaux */

        if (irank == nRank-1)
          localCoords  = (double *) malloc(sizeof(double) * 3 * localNVertex);
        else
          localCoords  = (double *) realloc(localCoords, sizeof(double) * 3 * localNVertex);

        for (int i = 0; i < localNVertex; i++) {
          localCoords[3*i]   = (*meshCoords)[3*(globalVertexNum[i]-1)];
          localCoords[3*i+1] = (*meshCoords)[3*(globalVertexNum[i]-1)+1];
          localCoords[3*i+2] = (*meshCoords)[3*(globalVertexNum[i]-1)+2];
        }

        /* Envoi des donnees maillages si different du proc 0 */

        if (irank != 0) {

          /* Envoi les infos concernant les sommets */
          int ierror = MPI_SUCCESS;
          ierror = MPI_Send(&localNVertex, 1, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS) {
            printf("Erreur MPI\n");
            exit(1);
          }

          ierror = MPI_Send(localCoords, 3*localNVertex, MPI_DOUBLE, irank, 0, localComm);
          if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
          }

          ierror = MPI_Send(globalVertexNum, localNVertex, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
          }

          /* Envoi les infos concernant les elements */

          ierror = MPI_Send(&localNElts, 1, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
          }

          ierror = MPI_Send(globalEltNum, localNElts, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
          }

          ierror = MPI_Send(localEltsConnecPointer, localNElts+1, MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
          }

          ierror = MPI_Send(localEltsConnec, localEltsConnecPointer[localNElts], MPI_INT, irank, 0, localComm);
          if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
          }
        }

      }

      /* Suppression des donnees globales du maillage */

      if (numDomElt != NULL)
        free(numDomElt);

      if ((*meshCoords) != NULL)
        free(*meshCoords);

      if ((*eltsConnecPointer) != NULL)
        free(*eltsConnecPointer);

      if ((*eltsConnec) != NULL)
        free(*eltsConnec);

      /* Contenu du maillage local sur proc 0 */


      localCoords = (double *) realloc(localCoords, sizeof(double) * 3 * localNVertex);
      localEltsConnec = (int *) realloc(localEltsConnec, sizeof(int) * localEltsConnecPointer[localNElts]);
      localEltsConnecPointer = (int *) realloc(localEltsConnecPointer, sizeof(int) * (localNElts + 1));
      globalVertexNum = (int *) realloc(globalVertexNum, sizeof(int) * localNVertex);
      globalEltNum = (int *) realloc(globalEltNum, sizeof(int) * localNElts);

      *nVertex = localNVertex;
      *nElts = localNElts;
      *meshCoords = localCoords;
      *eltsConnec = localEltsConnec;
      *eltsConnecPointer = localEltsConnecPointer;

    }

    /* Reception des donnees maillage si different du proc 0 */

    if (localRank != 0) {

      /* Reception des infos concernant les sommets */

      int ierror = MPI_SUCCESS;
      ierror = MPI_Recv(nVertex, 1, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
      }

      *meshCoords = (double *) malloc(sizeof(double) * 3*(*nVertex));

      ierror = MPI_Recv(*meshCoords, 3*(*nVertex), MPI_DOUBLE, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
      }

      globalVertexNum = (int *) malloc(sizeof(int) * (*nVertex));

      ierror = MPI_Recv(globalVertexNum, *nVertex, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
      }

      /* Reception des infos concernant les elements */

      ierror = MPI_Recv(nElts, 1, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
      }

      globalEltNum = (int *) malloc(sizeof(int) * (*nElts));

      ierror = MPI_Recv(globalEltNum, *nElts, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
      }

      *eltsConnecPointer = (int *) malloc(sizeof(int) * ((*nElts)+1));

      ierror = MPI_Recv(*eltsConnecPointer, (*nElts)+1, MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
      }

      *eltsConnec = (int *) malloc(sizeof(int) * (*nElts));

      ierror = MPI_Recv(*eltsConnec, (*eltsConnecPointer)[(*nElts)], MPI_INT, 0, MPI_ANY_TAG, localComm, &status);
      if (ierror != MPI_SUCCESS){
            printf("Erreur MPI\n");
            exit(1);
      }
    }
#else
    printf("To make parallel test : Build cwipi with metis \n");
    exit(1);
#endif
  }

  if (globalEltNum != NULL)
    free(globalEltNum);

  if (globalVertexNum != NULL)
    free(globalVertexNum);

}


void PROCF(creemaillagepolygone2d_f, CREEMAILLAGEPOLYGONE2D_F)(int *order,
							       MPI_Fint* localFComm,
                                                               /*MPI_Comm *localComm,*/
                                                               double *xmin,
                                                               double *xmax,
                                                               double *ymin,
                                                               double *ymax,
                                                               int *initRandom,
                                                               int *nx,
                                                               int *ny,
                                                               int *nVertex,
                                                               double *meshCoords_f,
                                                               int *nElts,
                                                               int *lEltsConnecPointer_f,
                                                               int *eltsConnecPointer_f,
                                                               int *eltsConnec_f
                                                               )
{
  MPI_Comm localComm = MPI_Comm_f2c(*localFComm);
  int nVertex_f = *nVertex;
  int nElts_f = *nElts;

  double *meshCoords = NULL;
  int    *eltsConnecPointer = NULL;
  int    *eltsConnec = NULL;

  creeMaillagePolygone2D(*order,
                         localComm,
                         *xmin,
                         *xmax,
                         *ymin,
                         *ymax,
                         *initRandom,
                         *nx,
                         *ny,
                         nVertex,
                         &meshCoords,
                         nElts,
                         &eltsConnecPointer,
                         &eltsConnec);

  if (nVertex_f < *nVertex) {
    printf("Augmenter le nombre de sommets Fortran a : %i \n", *nVertex);
    exit(1);
  }

  if (nElts_f < *nElts) {
    printf("Augmenter le nombre d'elements a : %i \n", *nElts);
    exit(1);
  }

  if (*lEltsConnecPointer_f < eltsConnecPointer[*nElts]) {
    printf("Augmenter la taille du tableau de connectivite a : %i \n", eltsConnecPointer[*nElts]);
    exit(1);
  }

  for(int i = 0; i < 3*(*nVertex); i++)
    meshCoords_f[i] = meshCoords[i];

  for(int i = 0; i < *nElts + 1; i++)
    eltsConnecPointer_f[i] = eltsConnecPointer[i];

  for(int i = 0; i < eltsConnecPointer[*nElts]; i++)
    eltsConnec_f[i] = eltsConnec[i];

  free(meshCoords);
  free(eltsConnecPointer);
  free(eltsConnec);
}



#ifdef __cplusplus
}
#endif /* __cplusplus */
