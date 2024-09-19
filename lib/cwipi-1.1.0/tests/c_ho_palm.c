/*
  This file is part of the CWIPI library. 

  Copyright (C) 2017  ONERA

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
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "grid_mesh.h"
#include <mpi.h>

#include "cwipi.h"







void _goto(FILE *f,char* word) 
{ char test[1000];
  int r;
  while(strcmp(test,word)!=0) {
      r = fscanf(f, "%s",test);
      if (r == EOF) 
        return EXIT_FAILURE;
      
 //     printf("%s\n",&test[0]);
  }
}

void _generate_gmsh_mesh(char* geofile,int localCommSize,int order){
  char *s = (char*) malloc(1000*sizeof(char));
  int len=0;
  for(len = 0; geofile[len] != '\0'; ++len);
  
  char* filename=(char*)malloc(sizeof(char)*(len-4));
  for (int i=0;i<len-4;i++) filename[i]=geofile[i];

  if(localCommSize>1) 
    sprintf(s, "gmsh %s -2 -format msh -order %i  -part %i -part_split -o %s.msh",geofile,order,localCommSize,filename);
  else
    sprintf(s, "gmsh %s -2 -format msh -order %i  -o %s.msh",geofile,order,filename);    
  system(s);
  free(s);
  free(filename);
}


typedef struct elType elType;
/* Declare the struct with integer members x, y */
struct elType {
   int    nNodes;
   char*    descri;
};









/*----------------------------------------------------------------------
 *                                                                     
 * Read mesh dimension                                             
 *                                                                     
 * parameters:
 *   f                   <-- Mesh file                 
 *   dimension           --> Dimension                   
 *   nb_Vertex             <-- number of vertices
 *   nElements           <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index  
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/
 static int _read_mesh_dim(FILE *f, 
                      int* nb_Vertex, 
                      int* nb_Elts,
                      int* nBlock,
                      int** nElBlock,
                      int** typeBlock) {

  elType* elementType = (elType*)malloc(sizeof(elType*)*100);
  
elementType[1] = (elType){2,"line"};
elementType[2] = (elType){3,"triangle"};
elementType[3] = (elType){4,"quadrangle"};
elementType[4] = (elType){4,"tetrahedron"};
elementType[5] = (elType){8,"hexahedron"};
elementType[6] = (elType){6,"prism"};
elementType[7] = (elType){5,"pyramid"};
elementType[8] = (elType){3,"second order line (2 nodes associated with the vertices and 1 with the edge)"};
elementType[9] = (elType){6,"second order triangle (3 nodes associated with the vertices and 3 with the edges)"};
elementType[10] = (elType){9,"second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)"};
elementType[11] = (elType){10,"second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)"};
elementType[12] = (elType){27,"second order hexahedron (8 nodes associated with the vertices, 12 with the edges 6 with the faces and 1 with the volume)"};
elementType[13] = (elType){18,"second order prism (6 nodes associated with the vertices], 9 with the edges and 3 with the quadrangular faces)"};
elementType[14] = (elType){14,"second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)"};
elementType[15] = (elType){1,"point"};
elementType[16] = (elType){8,"second order quadrangle (4 nodes associated with the vertices and 4 with the edges)"};
elementType[17] = (elType){20,"second order hexahedron (8 nodes associated with the vertices and 12 with the edges)"};
elementType[18] = (elType){15,"second order prism (6 nodes associated with the vertices and 9 with the edges)"};
elementType[19] = (elType){13,"second order pyramid (5 nodes associated with the vertices and 8 with the edges)"};
elementType[20] = (elType){9,"third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)"};
elementType[21] = (elType){10,"third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)"};
elementType[22] = (elType){12,"fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)"};
elementType[23] = (elType){15,"fourth order triangle (3 nodes associated with the vertices, 9 with the edges 3 with the face)"};
elementType[24] = (elType){15,"fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)"};
elementType[25] = (elType){21,"fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges 6 with the face)"};
elementType[26] = (elType){4,"third order edge (2 nodes associated with the vertices 2 internal to the edge)"};
elementType[27] = (elType){5,"fourth order edge (2 nodes associated with the vertices 3 internal to the edge)"};
elementType[28] = (elType){6,"fifth order edge (2 nodes associated with the vertices 4 internal to the edge)"};
elementType[29] = (elType){20,"third order tetrahedron (4 nodes associated with the vertices 12 with the edges 4 with the faces)"};
elementType[30] = (elType){35,"fourth order tetrahedron (4 nodes associated with the vertices 18 with the edges 12 with the faces 1 in the volume)"};
elementType[31] = (elType){56,"fifth order tetrahedron (4 nodes associated with the vertices 24 with the edges 24 with the faces 4 in the volume)"};
elementType[92] = (elType){64,"third order hexahedron (8 nodes associated with the vertices 24 with the edges 24 with the faces 8 in the volume)"};
elementType[93] = (elType){125,"fourth order hexahedron (8 nodes associated with the vertices 36 with the edges 54 with the faces 27 in the volume)"};

  int i, j, r;


  char test[1000];
  
  _goto(f,"$Nodes");
   // _goto(f,"$EndNodes");
  int nv,nEl;
  double poubd;
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nb_Vertex);
  
  printf("nb_Vertex %i\n",*nb_Vertex);
//while(1==1){}
  int dimension =3;

  for(int block = 1; block<(*nBlock)+1;block++) {
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&nv);

    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl);  
      nEl=nEl-1;
      r = fscanf(f, "%lf",&poubd);
      r = fscanf(f, "%lf",&poubd);
      r = fscanf(f, "%lf",&poubd);       

      if (r == EOF) 
        return EXIT_FAILURE;
    }
  }

  int IelType,nEl2,poub;

  _goto(f,"$Elements");
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nb_Elts);
  printf("nb_Elts %i\n",*nb_Elts);
  printf("nBlock %i\n",*nBlock);

  *nElBlock = (int*)malloc(sizeof(int)*(*nBlock));
  *typeBlock = (int*)malloc(sizeof(int)*(*nBlock));

  int block1 =0;
  for( int block=0;block<(*nBlock);block++) {
    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&IelType);  
    r = fscanf(f, "%i",&nv);
    //To use with Paraview
    if(IelType == 1 || IelType == 15) {/**nBlock=*nBlock-1;*//*block=block-1;*/*nb_Elts=*nb_Elts-nv;
       for(int s=0;s<nv*(1+elementType[IelType].nNodes);s++) {
         r = fscanf(f, "%i",&poub);
       }
    }
    else {
    (*nElBlock)[block1] = nv;
    (*typeBlock)[block1] = IelType;

    int size_el;
    size_el = elementType[IelType].nNodes;   

    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl2); 
      for(int jv=0;jv<size_el;jv++) { 
        r = fscanf(f, "%i",&poub);
      }

      if (r == EOF) 
        return EXIT_FAILURE;
    }
    block1++;
    }
   // free(connec);
  }  
  *nBlock = block1;
    printf("nb_Elts %i\n",*nb_Elts);
  return EXIT_SUCCESS;
}


int tabSearch2(int value,int** gnum,int gnum_size) {
 /* Déclarations */

 int POS;   /* position de la valeur */
 int I;     /* indice courant */
 int INF, MIL, SUP; /* limites du champ de recherche */
 
 /* Initialisation des limites du domaine de recherche */
 INF=0;
 SUP=gnum_size-1;
 /* Recherche de la position de la valeur */
 POS=INF;
 while (POS<SUP)
        {
          if(value == (*gnum)[POS]) break;
          POS++;
        }
  if(POS==SUP+1) POS=-1;
  return POS;
}



int tabSearch(int value,int** gnum,int gnum_size) {
 /* Déclarations */

 int POS;   /* position de la valeur */
 int I;     /* indice courant */
 int INF, MIL, SUP; /* limites du champ de recherche */
 
 /* Initialisation des limites du domaine de recherche */
 INF=0;
 SUP=gnum_size-1;
 /* Recherche de la position de la valeur */
 POS=-1;
 while ((INF<=SUP) && (POS==-1))
        {
         MIL=(SUP+INF)/2;
         //printf("value %i size %i MIL %i\n",value,MIL,gnum_size);
         if (value < (*gnum)[MIL])
               SUP=MIL-1;
         else if (value > (*gnum)[MIL])
               INF=MIL+1;
         else
               POS=MIL;
        }
 
  /* Edition du résultat */
 if (POS==-1) {
     printf("La valueeur recherchée ne se trouve pas "
            "dans le tableau. %i %i\n",value,gnum_size);
     return -1;
 }
 else {
     return POS;
 }
}


void tricroissant( int* a, int b) {
    int ind_min=0;
    int i=ind_min;
    int x=ind_min;
    int j=ind_min;
 
    for(i=ind_min;i<b;i++)
    {
        for(j=ind_min+1;j<b;j++)
        {
            if(a[i]<a[j])
            {
                x=a[i];
                a[i]=a[j];
                a[j]=x;
                j--;
                }
 
        }
 
        }
 
    x=a[ind_min];
    for(i=ind_min;i<b;i++)
    a[i]=a[i+1];
    a[b-1]=x;
 
}







/*----------------------------------------------------------------------
 *                                                                     
 * Read mesh dimension                                             
 *                                                                     
 * parameters:
 *   f                   <-- Mesh file                 
 *   dimension           --> Dimension                   
 *   nb_Vertex             <-- number of vertices
 *   nElements           <-- number of elements
 *   nConnecVertex       <-- size of connectivity
 *   coords              --> vertices coordinates
 *   connecPointer       --> connectivity index  
 *   connec              --> connectivity
 *---------------------------------------------------------------------*/
 static int _read_mesh(FILE *f, 
                      int* nb_Vertex, 
                      int* nb_Elts,
                      int* nBlock,
                      int** nElBlock,
                      int** typeBlock,
                      int** connec,
                      int** eltsConnecPointer,
                      double *coords,
                      int *gnum_coords) {

  elType* elementType = (elType*)malloc(sizeof(elType*)*100);
  
elementType[1] = (elType){2,"line"};
elementType[2] = (elType){3,"triangle"};
elementType[3] = (elType){4,"quadrangle"};
elementType[4] = (elType){4,"tetrahedron"};
elementType[5] = (elType){8,"hexahedron"};
elementType[6] = (elType){6,"prism"};
elementType[7] = (elType){5,"pyramid"};
elementType[8] = (elType){3,"second order line (2 nodes associated with the vertices and 1 with the edge)"};
elementType[9] = (elType){6,"second order triangle (3 nodes associated with the vertices and 3 with the edges)"};
elementType[10] = (elType){9,"second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)"};
elementType[11] = (elType){10,"second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)"};
elementType[12] = (elType){27,"second order hexahedron (8 nodes associated with the vertices, 12 with the edges 6 with the faces and 1 with the volume)"};
elementType[13] = (elType){18,"second order prism (6 nodes associated with the vertices], 9 with the edges and 3 with the quadrangular faces)"};
elementType[14] = (elType){14,"second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)"};
elementType[15] = (elType){1,"point"};
elementType[16] = (elType){8,"second order quadrangle (4 nodes associated with the vertices and 4 with the edges)"};
elementType[17] = (elType){20,"second order hexahedron (8 nodes associated with the vertices and 12 with the edges)"};
elementType[18] = (elType){15,"second order prism (6 nodes associated with the vertices and 9 with the edges)"};
elementType[19] = (elType){13,"second order pyramid (5 nodes associated with the vertices and 8 with the edges)"};
elementType[20] = (elType){9,"third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)"};
elementType[21] = (elType){10,"third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)"};
elementType[22] = (elType){12,"fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)"};
elementType[23] = (elType){15,"fourth order triangle (3 nodes associated with the vertices, 9 with the edges 3 with the face)"};
elementType[24] = (elType){15,"fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)"};
elementType[25] = (elType){21,"fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges 6 with the face)"};
elementType[26] = (elType){4,"third order edge (2 nodes associated with the vertices 2 internal to the edge)"};
elementType[27] = (elType){5,"fourth order edge (2 nodes associated with the vertices 3 internal to the edge)"};
elementType[28] = (elType){6,"fifth order edge (2 nodes associated with the vertices 4 internal to the edge)"};
elementType[29] = (elType){20,"third order tetrahedron (4 nodes associated with the vertices 12 with the edges 4 with the faces)"};
elementType[30] = (elType){35,"fourth order tetrahedron (4 nodes associated with the vertices 18 with the edges 12 with the faces 1 in the volume)"};
elementType[31] = (elType){56,"fifth order tetrahedron (4 nodes associated with the vertices 24 with the edges 24 with the faces 4 in the volume)"};
elementType[92] = (elType){64,"third order hexahedron (8 nodes associated with the vertices 24 with the edges 24 with the faces 8 in the volume)"};
elementType[93] = (elType){125,"fourth order hexahedron (8 nodes associated with the vertices 36 with the edges 54 with the faces 27 in the volume)"};

  int i, j, r;


  char test[1000];
  
  _goto(f,"$Nodes");
   // _goto(f,"$EndNodes");
  int nv,nEl;
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nb_Vertex);
  printf("nb_Vertex1 %i\n",*nb_Vertex);

  int dimension =3;

  int i_el=0;
  for(int block = 1; block<(*nBlock)+1;block++) {
    r = fscanf(f, "%i",&nv);  
  //  printf("nv %i\n",nv);
    r = fscanf(f, "%i",&nv);  
  //  printf("nv %i\n",nv);
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&nv);
 //   printf("nv %i\n",nv);  
 //   r = fscanf(f, "%i",&nv);  

    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",gnum_coords+i_el);  
      //gnum_coords[i_el]=nEl;
      r = fscanf(f, "%lf",coords+3*i_el);
      r = fscanf(f, "%lf",coords+3*i_el+1);
      r = fscanf(f, "%lf",coords+3*i_el+2);       
      printf("gnum_coords %i block %i nv %i x %f y %f z %f\n",gnum_coords[i_el],block,nv,coords[3*i_el],coords[3*i_el+1],coords[3*i_el+2]);
      if (r == EOF) 
        return EXIT_FAILURE;
      i_el++;
    }
  }
 // while(1==1){}
  int IelType,nEl2,poub;
  
  _goto(f,"$Elements");
  r = fscanf(f, "%i",nBlock);
  r = fscanf(f, "%i",nb_Elts);
  printf("nb_Elts %i\n",*nb_Elts);
    printf("nBlock %i\n",*nBlock);
    
  int** toto = *connec;
  
  int iop = *nb_Elts;
  
  int block1 = 0;

  for( int block=0;block<(*nBlock);block++) {
  

    r = fscanf(f, "%i",&nv);  
    r = fscanf(f, "%i",&nv);
    r = fscanf(f, "%i",&IelType);  
    r = fscanf(f, "%i",&nv);
    //To use with Paraview
    if(IelType == 1 || IelType == 15) {/**nBlock=*nBlock-1;*//*block=block-1;*/*nb_Elts=*nb_Elts-nv;
       for(int s=0;s<nv*(1+elementType[IelType].nNodes);s++) {
         r = fscanf(f, "%i",&poub);
       }
    }
    else {
    (*nElBlock)[block1] = nv;

    (*typeBlock)[block1] = IelType;

    int size_el;
    if(IelType!=2 && IelType!=3) printf("IelType %i nv %i\n",IelType,nv);
    size_el = elementType[IelType].nNodes;  
    

     
  //  printf("IelType %i nv %i size_el %i block %i block1 %i nBlock %i\n",IelType,nv,size_el,block,block1,*nBlock); 
   //   if(jv==1) while(1==1){}
  //  if(block1==0)while(1==1){}
  
    eltsConnecPointer[block1][0]=0;
  
    for (int i = 0; i < nv; i++) {
      r = fscanf(f, "%i",&nEl2); 
      for(int jv=0;jv<size_el;jv++) { 
        r = fscanf(f, "%i",connec[block1]+size_el*i+jv);
        
        connec[block1][size_el*i+jv] =  1+tabSearch2(connec[block1][size_el*i+jv],&gnum_coords,*nb_Vertex);
        //printf("connect %i \n",connec[block1][size_el*i+jv]);
     //   printf("%i nEl2 %i block %i nv %i connec %i %i %i size_el %i IelType %i\n",
      //  *nb_Elts,nEl2,block,nv,i,jv,connec[block1][size_el*i+jv],size_el,IelType);
      }
      
     eltsConnecPointer[block1][i+1]=eltsConnecPointer[block1][i]+size_el;
      if (r == EOF) 
        return EXIT_FAILURE;
    }
    block1++;
    }
   // free(connec);
  }  
  *nBlock = block1;
  printf("nb_Elts %i *nBlock %i\n",*nb_Elts,*nBlock);
  

  return EXIT_SUCCESS;
}



int sizeForType(int type) {
  if(type == 1)  return 2;
  if(type == 2)  return 3;
  if(type == 3)  return 4;
  if(type == 9)  return 6;
  if(type == 15) return 1;
 }


/*----------------------------------------------------------------------
 *                                                                     
 * Main : linear coupling test                                         
 *
 *---------------------------------------------------------------------*/
 
int main
(
 int    argc,    /* Nombre d'arguments dans la ligne de commandes */
 char  *argv[]   /* Tableau des arguments de la ligne de commandes */
)
{

  FILE *outputFile;

  MPI_Init(&argc, &argv);

  int rank;
  int commWorldSize;
  

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &commWorldSize);


  /* Initialization
   * -------------- */
  double *times_init = NULL;

  char *codeName;
  int codeId;
  char *codeCoupledName;

  if (rank < commWorldSize / 2) {
    codeName = "code1";
    codeId = 1;
    codeCoupledName = "code2";
  }
  else {
    codeName = "code2";
    codeId = 2;
    codeCoupledName = "code1";
  }

  char* fileName = (char *) malloc(sizeof(char) * 256);
  sprintf(fileName,"c_surf_cpl_palm_%4.4d.txt",rank);

  outputFile = fopen(fileName,"w");

  free(fileName);

  cwipi_set_output_listing(outputFile);

  MPI_Comm localComm;
  cwipi_init(MPI_COMM_WORLD,
             codeName ,
             &localComm);

  /* Output redirection
   * ------------------ */

  int currentRank;
  int localCommSize;

  MPI_Comm_rank(localComm, &currentRank);
  MPI_Comm_size(localComm, &localCommSize);



  char* geofile = CWP_MESH_DIR"sphere_palm.geo";
  int order = 2;

  if(currentRank==0) {
    _generate_gmsh_mesh(geofile,localCommSize,order);
  }
  MPI_Status status;
  MPI_Barrier(MPI_COMM_WORLD);
  
  int len=0;
  for(len = 0; geofile[len] != '\0'; ++len);
  
  int lenfile;
  
  if((1+currentRank)>=1000 && (1+currentRank)<10000) lenfile=len+5;
  else if((1+currentRank)>=100) lenfile=len+4;
  else if((1+currentRank)>=10)  lenfile=len+3;  
  else lenfile=len+2;

  if(localCommSize==1) {
   lenfile ==len;
  }

  char* root=(char*)malloc(sizeof(char)*(len-4));
  for (int i=0;i<len-4;i++) root[i]=geofile[i];   

  char* filename=(char*)malloc(sizeof(char)*(lenfile));    
  if(localCommSize!=1) sprintf(filename,"%s_%i.msh",root,1+currentRank);
  else sprintf(filename,"%s.msh",root);

  MPI_Barrier(MPI_COMM_WORLD);
  printf("filename %s %i\n",filename,currentRank);

  //while(1==1){}
  FILE* meshFile;
  meshFile = fopen(filename, "r");

  int nElts,nBlock,nb_Vertex;
 
 
  int* typeBlock= NULL;
  int* nElBlock= NULL;

  _read_mesh_dim(meshFile,
              &nb_Vertex, 
              &nElts,
              &nBlock,
              &(nElBlock),
              &(typeBlock));


  double* coords = (double*) malloc(3 * nb_Vertex * sizeof(double));
  int* gnum_coord = (int*) malloc(nb_Vertex * sizeof(int));
  int** eltsConnec = (int**)malloc(nBlock*sizeof(int*));
  int** eltsConnecPointer = (int**)malloc(nBlock*sizeof(int*));


  
  printf("currentRank %i nb_Vertex %i nElts %i nBlock %i nElBlock[0] %i typeBlock[0] %i\n",
  currentRank, nb_Vertex,nElts,nBlock,nElBlock[0],typeBlock[0]);
  
  for(int b=0;b<nBlock;b++) {
    //printf("Partition %i Block %i size = %i rank %i\n",i_part,b,nElBlock[b],rank);
    eltsConnec[b]=(int*)malloc(sizeForType(typeBlock[b])*nElBlock[b]*sizeof(int));
    eltsConnecPointer[b]=(int *) malloc(sizeof(int) * (nElBlock[b] + 1));
  }
  
  fclose(meshFile);
  meshFile = fopen(filename, "r");

  _read_mesh(meshFile,
            &nb_Vertex, 
            &nElts,
            &nBlock,
            &(nElBlock),
            &(typeBlock),
            eltsConnec,
            eltsConnecPointer,
            coords,
            gnum_coord);


  fprintf(outputFile, "  Surface coupling test : location in tria P2\n");
  fprintf(outputFile, "\n");

  fprintf(outputFile, "\nDump after initialization\n");
  fprintf(outputFile, "-------------------------\n");
  cwipi_dump_application_properties();

  if (rank == 0)
    printf("        Create coupling\n");
  
  cwipi_solver_type_t solver_type;
  
  solver_type = CWIPI_SOLVER_CELL_VERTEX;

  /* Coupling creation
   * ----------------- */

  const int postFreq = 1;
  
  cwipi_create_coupling("c_surf_cpl_palm",                                // Coupling id
                        CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING, // Coupling type
                        codeCoupledName,                           // Coupled application id
                        2,                                         // Geometric entities dimension
                        1e-1,                                       // Geometric tolerance
                        CWIPI_STATIC_MESH,                         // Mesh type
                        solver_type,                               // Solver type
                        postFreq,                                  // Postprocessing frequency
                        "EnSight Gold",                            // Postprocessing format
                        "text");                                   // Postprocessing option

  cwipi_ho_define_mesh("c_surf_cpl_palm",
                       nb_Vertex,
                       nElts,
                       order,
                       coords,
                       eltsConnecPointer[0],
                       eltsConnec[0]);


  const int n_node = 6;
  
  int *ijk = malloc(sizeof(int)*2*n_node);

  ijk[0] = 0;
  ijk[1] = 0;

  ijk[2] = 2;
  ijk[3] = 0;

  ijk[4] = 0;
  ijk[5] = 2;
  
  ijk[6] = 1;
  ijk[7] = 0;

  ijk[8] = 1;
  ijk[9] = 1;

  ijk[10] = 0;
  ijk[11] = 1;
  
  cwipi_ho_ordering_from_IJK_set ("c_surf_cpl_palm",
                                  CWIPI_FACE_TRIAHO,
                                  n_node,
                                  ijk);

  double *sendValues = NULL;
  double *recvValues = NULL;
  
  sendValues = (double *) malloc(sizeof(double) * nb_Vertex);
  recvValues = (double *) malloc(sizeof(double) * nb_Vertex);

  for (int i = 0; i < nb_Vertex; i++) {
    if (rank == 0)
      sendValues[i] = coords[3 * i];
    else
      sendValues[i] = coords[3 * i + 1];
  }

  int nNotLocatedPoints = 0;
  cwipi_exchange_status_t status_exch = cwipi_exchange("c_surf_cpl_palm",
                                                  "echange1",
                                                  1,
                                                  1,     // n_step
                                                  0.1,   // physical_time
                                                  "cooX",
                                                  sendValues,
                                                  "cooY",
                                                  recvValues,
                                                  &nNotLocatedPoints);




  /* Coupling deletion
   * ----------------- */

  if (rank == 0)
    printf("        Delete coupling\n");

  cwipi_delete_coupling("c_surf_cpl_palm");

  


  /* Finalize
   * -------- */
  
  cwipi_finalize();
  
  MPI_Finalize();

  fclose(outputFile);
  
  return EXIT_SUCCESS;

}





