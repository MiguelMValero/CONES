/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_mpi.h"
#include "pdm_writer.h"
#include "pdm_writer_priv.h"
#include "pdm_writer_ensight.h"
#include "pdm_binary_search.h"
#include "pdm_writer_ensight_case.h"
#include "pdm_io.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_array.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_mesh_nodal.h"

/*=============================================================================
 * Definitions des macro
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* Longueur max d'une ligne pour le format Ensight Gold
 * auxquels il faut ajouter le `\n' et le `\0'
 * pour l'affectation dans la chaîne réceptrice */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Particularisation de la structure CS
 *----------------------------------------------------------------------------*/

typedef struct {

  PDM_writer_ensight_case_t *ensight_case; /* Gestion du fichier case */
  PDM_io_file_t          *f_unit_geom;  /* Unite du fichier de géométrie */
  int                        n_time_step;  /* Nombre de pas de temps */
  int                        n_part_ecr;   /* Nombre de parts ensight
                                                écrites dans le fichier chr.geo */
  int                      append;  /* 1 == 1 seul fichier par viraible et 
                                          par geometrie sinon 1 pour chaque instant*/
} PDM_writer_ensight_t;


/*----------------------------------------------------------------------------
 * Particularisation de la structure PDM_writer_geom
 *----------------------------------------------------------------------------*/

typedef struct {

  int num_part; /* Numero de part associe a la geometrie */

} PDM_writer_geom_ensight_t;


/*----------------------------------------------------------------------------
 * Particularisation de la structure PDM_writer_var
 *----------------------------------------------------------------------------*/

typedef struct {

  // int f_unit; /* Unite du fichier associe à la variable */
  PDM_io_file_t *f_unit; /* Unite du fichier associe à la variable */

} PDM_writer_var_ensight_t;

/*=============================================================================
 * Variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Longueur max d'une chaine
 *----------------------------------------------------------------------------*/

//static const int _l_max_chaine_ens = 1024;

/*----------------------------------------------------------------------------
 * Longueur max d'une chaine
 *----------------------------------------------------------------------------*/

static const char  *_ensight_type_name[10] = {"point",
                                              "bar2",
                                              "tria3",
                                              "quad4",
                                              "nsided",
                                              "tetra4",
                                              "pyramid5",
                                              "penta6",
                                              "hexa8",
                                              "nfaced"};

/*=============================================================================
 * Fonctions privées
 *============================================================================*/

/*----------------------------------------------------------------------------
 *
 * Fonction max
 *
 * parameters:
 *   a           <-- Premiere valeur
 *   b           <-- Seconde valeur
 *
 * return:
 *   max(a, b)
 *
 *----------------------------------------------------------------------------*/

static inline PDM_g_num_t
_max
(
 PDM_g_num_t a,
 PDM_g_num_t b
)
{
  return ((a) > (b) ? (a) : (b));
}

static inline int
_max_int
(
 int a,
 int b
)
{
  return ((a) > (b) ? (a) : (b));
}

/*----------------------------------------------------------------------------
 * Write string to a text or C binary EnSight Gold file
 *
 * parameters:
 *   cs <-- Structure cs courante
 *   f  <-- File to write to
 *   s  <-- String to write
 *----------------------------------------------------------------------------*/

static void
_ecr_string
(
 PDM_writer_t     *cs,
 PDM_io_file_t *f_unit_geom,
 const char       *s
 )
{
  size_t  i;
  char  buf[82];

  if (cs->fmt_fic == PDM_WRITER_FMT_ASCII) {
    strncpy(buf, s, 80);
    buf[80] = '\0';
    buf[81] = '\n';
    PDM_io_fmt_data_set(f_unit_geom,
                          1,
                          PDM_IO_T_CHAR,
                          "%c");
    size_t s_buf =  strlen(buf);
    PDM_io_global_write(f_unit_geom,
                            (PDM_l_num_t) sizeof(char),
                            (PDM_l_num_t) s_buf,
                            buf);
  }

  else if (cs->fmt_fic == PDM_WRITER_FMT_BIN) {
    strncpy(buf, s, 80);
    buf[80] = '\0';
    for (i = strlen(buf); i < 80; i++)
      buf[i] = ' ';
    PDM_io_global_write(f_unit_geom, sizeof(char), 80, buf);
  }
}

/*----------------------------------------------------------------------------
 * Ecriture d'un entier
 *
 * parameters:
 *   cs <-- Structure cs courante
 *   f  <-- File to write to
 *   n  <-- Integer value to write
 *----------------------------------------------------------------------------*/

inline static void
_ecr_int
(
 PDM_writer_t     *cs,
 PDM_io_file_t *f_unit_geom,
 int32_t           n
 )
{
  if (cs->fmt_fic == PDM_WRITER_FMT_ASCII) {
    PDM_io_fmt_data_set(f_unit_geom, 10, PDM_IO_T_INT, "%10d");
  }
  PDM_io_global_write(f_unit_geom, sizeof(int32_t), 1, &n);
}

/*----------------------------------------------------------------------------
 * Ecriture d'un tableau de flottant avec tri suivant la numérotation
 * absolue
 *
 * parameters:
 *   cs     <-- Structure cs courante
 *   f      <-- File to write to
 *   values <-- Integer value to write
 *----------------------------------------------------------------------------*/

static void
_ecr_entrelace_float
(
 PDM_writer_t                 *cs,
 const PDM_writer_status_t     s_ecr_n_valeur,
 PDM_io_file_t             *f_unit_geom,
 const PDM_stride_t  t_comp,
 const PDM_l_num_t            *n_comp,
 const PDM_l_num_t             n_valeur,
 const PDM_g_num_t            *indirection,
 const float                  *valeurs
 )
{

  if (s_ecr_n_valeur == PDM_WRITER_ON) {

    PDM_g_num_t n_val_abs_loc = 0;
    PDM_g_num_t n_val_abs     = 0;

    for (int i = 0; i < n_valeur; i++) {
      n_val_abs_loc = _max(n_val_abs_loc, indirection[i]);
    }

    PDM_MPI_Allreduce(&n_val_abs_loc, &n_val_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, cs->pdm_mpi_comm);

    int32_t n_val_abs_32_t = (int32_t) n_val_abs;

    if (cs->fmt_fic == PDM_WRITER_FMT_ASCII) {
      char  buf[12];
      int n_val = sprintf(buf, "%10d", n_val_abs_32_t);
      PDM_io_fmt_data_set(f_unit_geom, 1, PDM_IO_T_CHAR, "%c");
      PDM_io_global_write(f_unit_geom, sizeof(char), n_val, buf);
    }

    else if (cs->fmt_fic == PDM_WRITER_FMT_BIN) {
      PDM_io_global_write(f_unit_geom, sizeof(int32_t), 1, &n_val_abs_32_t);
    }
  }

  PDM_io_fmt_data_set(f_unit_geom,
                          12,
                          PDM_IO_T_FLOAT,
                          "%12.5e");
  PDM_io_par_interlaced_write(f_unit_geom,
                              t_comp,
                              n_comp,
                              sizeof(float),
                              n_valeur,
                              indirection,
                              (void *) valeurs);
}


/*----------------------------------------------------------------------------
 * Ecriture d'un tableau de flottant avec tri suivant la numérotation
 * absolue
 *
 * parameters:
 *   cs     <-- Structure cs courante
 *   f      <-- File to write to
 *   values <-- Integer value to write
 *----------------------------------------------------------------------------*/

static void
_ecr_entrelace_int
(
 PDM_writer_t                 *cs,
 const PDM_writer_status_t     s_ecr_n_valeur,
 PDM_io_file_t             *f_unit_geom,
 const PDM_stride_t  t_comp,
 const PDM_l_num_t            *n_comp,
 const PDM_l_num_t             n_valeur,
 const PDM_g_num_t            *indirection,
 const int32_t                *valeurs
 )
{
  if (s_ecr_n_valeur == PDM_WRITER_ON) {

    PDM_g_num_t n_val_abs_loc = 0;
    PDM_g_num_t n_val_abs     = 0;

    for (int i = 0; i < n_valeur; i++) {
      n_val_abs_loc = _max(n_val_abs_loc, indirection[i]);
    }

    PDM_MPI_Allreduce(&n_val_abs_loc, &n_val_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, cs->pdm_mpi_comm);

    int32_t n_val_abs_32_t = (int32_t) n_val_abs;

    if (cs->fmt_fic == PDM_WRITER_FMT_ASCII) {
      char  buf[12];
      int n_val = sprintf(buf, "%10d", n_val_abs_32_t);
      PDM_io_fmt_data_set(f_unit_geom, 1, PDM_IO_T_CHAR, "%c");
      PDM_io_global_write(f_unit_geom, sizeof(char), n_val, buf);
    }

    else if (cs->fmt_fic == PDM_WRITER_FMT_BIN) {
      PDM_io_global_write(f_unit_geom, sizeof(int32_t), 1, &n_val_abs_32_t);
    }
  }

  PDM_io_fmt_data_set(f_unit_geom,
                          10,
                          PDM_IO_T_INT,
                          "%10d");

  PDM_io_par_interlaced_write(f_unit_geom,
                              t_comp,
                              n_comp,
                              sizeof(int32_t),
                              n_valeur,
                              indirection,
                              (void *) valeurs);
}

/*----------------------------------------------------------------------------
 * Ecriture de l'entete d'un fichier geom Ensight
 *
 * parameters:
 *   cs           <-- Pointeur sur la structure cs courante
 *   f_unit_geom  <-- Unité du fichier PDM_io
 *----------------------------------------------------------------------------*/

static void
_geom_entete_ecr
(
 PDM_writer_t     *cs,
 PDM_io_file_t *f_unit_geom,
 PDM_writer_ensight_t *PDM_writer_ensight
 )
{
  if (cs->fmt_fic == PDM_WRITER_FMT_BIN)
    _ecr_string(cs, f_unit_geom, "C Binary");

  if (PDM_writer_ensight->append) {
//  if ((PDM_writer_ensight->append) && (cs->topologie != PDM_WRITER_TOPO_CST)) {
    char buff_append[81];
    strncpy(buff_append, "BEGIN TIME STEP", 80);
    buff_append[80] = '\0';  
    _ecr_string(cs,
               f_unit_geom,
               buff_append);
  }

  /* 1st description line */
  {
    char buf[81] = "";
    if (cs->nom_sortie != NULL)
      strncpy(buf, cs->nom_sortie, 80);
    buf[80] = '\0';
    _ecr_string(cs, f_unit_geom, buf);
  }
  /* 2nd description line */
  _ecr_string(cs, f_unit_geom, "Nothing");
  _ecr_string(cs, f_unit_geom, "node id assign");
  _ecr_string(cs, f_unit_geom, "element id assign");
}



/*----------------------------------------------------------------------------
 * Calcule la numerotation aboslue des faces des polyèdres
 *
 * parameters :
 *   geom            <-- Geometrie associee
 *
 *----------------------------------------------------------------------------*/

static void
_calcul_numabs_face_poly3d
(
 PDM_writer_geom_t  *geom,
 const int           iblock,
 PDM_g_num_t       **numabs_face
)
{

  /* Calcul de la taille du bloc du processus courant */

  int n_procs = 0;
  PDM_MPI_Comm_size(geom->pdm_mpi_comm,
                &n_procs);

  int i_proc = 0;
  PDM_MPI_Comm_rank(geom->pdm_mpi_comm,
                &i_proc);

  PDM_g_num_t *d_elt_proc =
          (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_procs + 1));

  /* Calcul du nombre d'elements abs du bloc
     repartis sur l'ensemble des processus */

  PDM_g_num_t max_loc = 0;
  PDM_g_num_t max_abs = 0;
  PDM_l_num_t n_elt_proc = 0;

  const int n_part = PDM_part_mesh_nodal_n_part_get(geom->mesh_nodal);

  for (int i = 0; i < n_part; i++) {

    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, iblock, i);
    PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                      iblock,
                                                                      i,
                                                                      PDM_OWNERSHIP_KEEP);

    n_elt_proc += n_elt;
    for (int j = 0; j < n_elt; j++) {
      max_loc = _max((PDM_g_num_t) numabs_block[j], max_loc);
    }
  }
  PDM_UNUSED(n_elt_proc);

  PDM_MPI_Allreduce(&max_loc, &max_abs, 1,
                    PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, geom->pdm_mpi_comm);

  /* Tri a l'aide d'un tableau defini par blocs continus
     repartis sur l'ensemble des processus */

  PDM_g_num_t div_entiere = max_abs / n_procs;
  PDM_g_num_t div_reste   = max_abs % n_procs;

  d_elt_proc[0] = 1;
  for (int i = 0; i < n_procs; i++) {
    d_elt_proc[i+1] =  div_entiere;
    if (i < div_reste)
      d_elt_proc[i+1] += 1;
  }

  /* Calcul de la repartition des elements sur les processus */

  for (int j = 0; j < n_procs; j++) {
    d_elt_proc[j+1] += d_elt_proc[j] ;
  }

  /* Allocation des tableaux pour echanges MPI */

  int *send_buff_n   = (int *) malloc (sizeof(int) * n_procs);
  int *send_buff_idx = (int *) malloc (sizeof(int) * n_procs);

  int *recv_buff_n   = (int *) malloc (sizeof(int) * n_procs);
  int *recv_buff_idx = (int *) malloc (sizeof(int) * n_procs);

  /* Calcul du nombre total d'elements du bloc */

  PDM_l_num_t n_elt_loc_total = 0;

  for (int j = 0; j < n_part; j++) {

    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, iblock, j);

    n_elt_loc_total += n_elt;
  }

  /* Comptage du nombre d'elements a envoyer a chaque processus */

  for (int j = 0; j < n_procs; j++) {
    send_buff_n[j]   = 0;
    send_buff_idx[j] = 0;
    recv_buff_n[j]   = 0;
    recv_buff_idx[j] = 0;
  }

  for (int j = 0; j < n_part; j++) {

    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, iblock, j);
    PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                      iblock,
                                                                      j,
                                                                      PDM_OWNERSHIP_KEEP);

    for (int k = 0; k < n_elt; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long(numabs_block[k],
                                                        d_elt_proc,
                                                        n_procs+1);
      send_buff_n[i_elt_proc] += 1;
    }

  }


  send_buff_idx[0] = 0;
  for (int j = 1; j < n_procs; j++) {
    send_buff_idx[j] = send_buff_idx[j-1] + send_buff_n[j-1];
  }

  /* Determination du nombre d'elements recu de chaque processus */

  PDM_MPI_Alltoall(send_buff_n,
                   1,
                   PDM_MPI_INT,
                   recv_buff_n,
                   1,
                   PDM_MPI_INT,
                   geom->pdm_mpi_comm);

  recv_buff_idx[0] = 0;
  for(int j = 1; j < n_procs; j++) {
    recv_buff_idx[j] = recv_buff_idx[j-1] + recv_buff_n[j-1];
  }

  /* Transmission : des numeros absolus des éléments + nb de faces
     Comme les valeurs n'ont pas le même type on stocke dans un tableau de char */

  const int n_octet_exch = sizeof(int) + sizeof(PDM_g_num_t); /* Nb d'octet échangés */

  for (int j = 0; j < n_procs; j++) {
    send_buff_idx[j] = send_buff_idx[j] * n_octet_exch;
    recv_buff_idx[j] = recv_buff_idx[j] * n_octet_exch;
    recv_buff_n[j]   = recv_buff_n[j]   * n_octet_exch;
  }

  unsigned char *send_buff_data =
    (unsigned char *) malloc(sizeof(unsigned char) * n_elt_loc_total * n_octet_exch);
  unsigned char *recv_buff_data =
    (unsigned char *) malloc(sizeof(unsigned char) * (recv_buff_idx[n_procs - 1] +
                                                      recv_buff_n[n_procs - 1]) * n_octet_exch);

  PDM_array_reset_int(send_buff_n, n_procs, 0);


  unsigned char *current_data = send_buff_data;
  for (int j = 0; j < n_part; j++) {
      int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, iblock, j);
      PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                        iblock,
                                                                        j,
                                                                        PDM_OWNERSHIP_KEEP);

      PDM_l_num_t   n_face;
      PDM_l_num_t  *facvtx_idx;
      PDM_l_num_t  *facvtx;
      PDM_l_num_t  *cellfac_idx;
      PDM_l_num_t  *cellfac;
      PDM_g_num_t  *face_ln_to_gn;
      PDM_g_num_t  *cell_ln_to_gn;
      int          *parent_num;
      PDM_g_num_t  *parent_entity_g_num;
      PDM_part_mesh_nodal_section_poly3d_get(geom->mesh_nodal,
                                             iblock,
                                             j,
                                             &n_face,
                                             &face_ln_to_gn,
                                             &facvtx_idx,
                                             &facvtx,
                                             &cell_ln_to_gn,
                                             &cellfac_idx,
                                             &cellfac,
                                             &parent_num,
                                             &parent_entity_g_num,
                                             PDM_OWNERSHIP_KEEP);

    for (int k = 0; k < n_elt; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long (numabs_block[k],
                                                         d_elt_proc,
                                                         n_procs+1);

      PDM_g_num_t *current_data_long =
        (PDM_g_num_t *) (current_data + send_buff_idx[i_elt_proc] + send_buff_n[i_elt_proc]);
      *current_data_long = numabs_block[k];
      int *current_dataInt =
        (int *) (current_data + send_buff_idx[i_elt_proc] + send_buff_n[i_elt_proc] + sizeof(PDM_g_num_t));
      *current_dataInt = cellfac_idx[k+1] - cellfac_idx[k];

      send_buff_n[i_elt_proc] += n_octet_exch;

    }
  }

  PDM_MPI_Alltoallv((void *) send_buff_data,
                    send_buff_n,
                    send_buff_idx,
                    PDM_MPI_BYTE,
                    (void *) recv_buff_data,
                    recv_buff_n,
                    recv_buff_idx,
                    PDM_MPI_BYTE,
                    geom->pdm_mpi_comm);

  /* Tri des éléments locaux détermination */

  int n_elt_recv = (recv_buff_idx[n_procs-1] + recv_buff_n[n_procs-1]) / n_octet_exch;

  PDM_g_num_t *face_abs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * (d_elt_proc[i_proc+1] - d_elt_proc[i_proc] + 1));
  PDM_g_num_t n_face_proc = 0;

  face_abs[0] = 0;
  int current_octet = 0;
  for (int i = 0; i < n_elt_recv; i++) {
    PDM_g_num_t *recv_buff_data_long = (PDM_g_num_t *) (recv_buff_data + current_octet);
    int       *recv_buff_data_int = (int *) (recv_buff_data + current_octet + sizeof(PDM_g_num_t));
    PDM_g_num_t  gCel  = *recv_buff_data_long;
    PDM_g_num_t _lCel  = gCel - d_elt_proc[i_proc]; // local numbering
    int          lCel  = (int) _lCel; // local numbering
    face_abs[lCel+1] = (PDM_g_num_t) *recv_buff_data_int;
    n_face_proc     += face_abs[lCel+1];
    current_octet += n_octet_exch;
  }

  /* Echange de la somme des faces des polyèdres stockes sur chaque processus  */

  PDM_g_num_t *n_face_procs = (PDM_g_num_t *) malloc (sizeof(PDM_g_num_t) * (n_procs+1));

  PDM_MPI_Allgather((void *) &n_face_proc,
                1,
                PDM__PDM_MPI_G_NUM,
                (void *) (n_face_procs + 1),
                1,
                PDM__PDM_MPI_G_NUM,
                geom->pdm_mpi_comm);

  n_face_procs[0] = 1;
  for (int j = 1; j < n_procs + 1; j++) {
    n_face_procs[j] += n_face_procs[j-1];
  }

  /* Determination la numerotation absolue des faces
     independante du parallelisme */

  for (int i = 1; i < d_elt_proc[i_proc+1] - d_elt_proc[i_proc] + 1; i++) {
    face_abs[i] = face_abs[i] + face_abs[i-1];
  }

  for (int i = 0; i < d_elt_proc[i_proc+1] - d_elt_proc[i_proc] + 1; i++) {
    face_abs[i] += n_face_procs[i_proc];
  }

  free(n_face_procs);

  /* Retour à l'envoyeur de la numérotation absolue */

  /* On remplit le buffer de reception qui devient le buffer d'envoi
     Le buffer d'envoi devient lui le buffer de reception */

  current_octet = 0;
  for (int i = 0; i < n_elt_recv; i++) {
    PDM_g_num_t *recv_buff_data_long = (PDM_g_num_t *) (recv_buff_data + current_octet);
    int         *recv_buff_data_int  = (int *        ) (recv_buff_data + current_octet + sizeof(PDM_g_num_t));
    PDM_g_num_t  g_cell  = *recv_buff_data_long;
    PDM_g_num_t _l_cell  = g_cell - d_elt_proc[i_proc];
    int        l_cell    = (int) _l_cell; // local numbering
    *recv_buff_data_long = face_abs[l_cell];
    *recv_buff_data_int = -1; /* Pas d'info dans la deuxieme partie du buffer */
    current_octet += n_octet_exch;
  }

  free(face_abs);

  PDM_MPI_Alltoallv((void *) recv_buff_data,
                    recv_buff_n,
                    recv_buff_idx,
                    PDM_MPI_BYTE,
                    (void *) send_buff_data,
                    send_buff_n,
                    send_buff_idx,
                    PDM_MPI_BYTE,
                    geom->pdm_mpi_comm);

  /* On Stocke l'information recue */

  PDM_array_reset_int(send_buff_n, n_procs, 0);

  current_data = send_buff_data;
  for (int j = 0; j < n_part; j++) {

    int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, iblock, j);
    PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                      iblock,
                                                                      j,
                                                                      PDM_OWNERSHIP_KEEP);


    for (int k = 0; k < n_elt; k++) {
      const int i_elt_proc = PDM_binary_search_gap_long (numabs_block[k],
                                                         d_elt_proc,
                                                         n_procs+1);

      PDM_g_num_t *current_data_long =
        (PDM_g_num_t *) (current_data + send_buff_idx[i_elt_proc] + send_buff_n[i_elt_proc]);
      numabs_face[j][k] = (PDM_g_num_t) *current_data_long;
      send_buff_n[i_elt_proc] += n_octet_exch;
    }
  }


  /* Liberation memoire */

  free(send_buff_n);
  free(send_buff_idx);
  free(send_buff_data);
  free(recv_buff_n);
  free(recv_buff_idx);
  free(recv_buff_data);
  free(d_elt_proc);

}


/*----------------------------------------------------------------------------
 * Fermeture du fichier de géométrie
 *
 * parameters:
 *   cs           <-- Pointeur sur la structure cs courante
 *----------------------------------------------------------------------------*/

static void
_geom_close(PDM_writer_t *cs)
{
  int rank = 0;
  PDM_MPI_Comm_rank(cs->pdm_mpi_comm, &rank);

  PDM_writer_ensight_t *PDM_writer_ensight = (PDM_writer_ensight_t *) cs->sortie_fmt;

  if (PDM_writer_ensight->f_unit_geom != NULL) {
    PDM_io_close(PDM_writer_ensight->f_unit_geom);
    double t_cpu;
    double t_elapsed;
    PDM_io_timer_total_get(PDM_writer_ensight->f_unit_geom, &t_cpu, &t_elapsed);
    const char * nom_fichier = PDM_io_file_name_get(PDM_writer_ensight->f_unit_geom);
    if (1 == 0) {
      if (rank == 0) {
        PDM_printf("Temps elapsed d'ecriture du fichier '%s' : %12.5e s\n", nom_fichier, t_elapsed);
        PDM_printf("Temps cpu d'ecriture du fichier '%s' : %12.5e s\n", nom_fichier, t_cpu);
      }
    }
    PDM_io_free(PDM_writer_ensight->f_unit_geom);
  }
  PDM_writer_ensight->f_unit_geom = NULL;
}


/*----------------------------------------------------------------------------
 * Fermeture d'un fichier de variable
 *
 * parameters:
 *   cs           <-- Pointeur sur la structure cs courante
 *----------------------------------------------------------------------------*/

static void
_var_close(PDM_writer_var_ensight_t *var, const int rank)
{
  if (var->f_unit != NULL) {
    PDM_io_close(var->f_unit);
    double t_cpu;
    double t_elapsed;
    PDM_io_timer_total_get(var->f_unit, &t_cpu, &t_elapsed);
    const char * nom_fichier = PDM_io_file_name_get(var->f_unit);
    if (1 == 0) {
      if (rank == 0) {
        PDM_printf("Temps elapsed d'ecriture du fichier '%s' : %12.5e s\n", nom_fichier, t_elapsed);
        PDM_printf("Temps cpu d'ecriture du fichier '%s' : %12.5e s\n", nom_fichier, t_cpu);
      }
    }
    PDM_io_free(var->f_unit);
    var->f_unit = NULL;
  }
}


/*----------------------------------------------------------------------------
 * Fermeture des fichiers de variable
 *
 * parameters:
 *   cs           <-- Pointeur sur la structure cs courante
 *----------------------------------------------------------------------------*/

static void
_vars_close(PDM_writer_t *cs)
{
  int rank = 0;
  PDM_MPI_Comm_rank(cs->pdm_mpi_comm, &rank);

  if (cs->var_tab != NULL) {
    const int n_ind = cs->var_tab->n_var;

    for (int i = 0; i < n_ind; i++) {
      PDM_writer_var_t *var = cs->var_tab->var[i];
      if (var != NULL) {
        PDM_writer_var_ensight_t *_var_ensight = (PDM_writer_var_ensight_t *) var->var_fmt;
        _var_close(_var_ensight, rank);
      }
    }
  }
}


/*=============================================================================
 * Fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Particularise la creation d'un objet CS (Cedre Sortie)
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_create
(
PDM_writer_t *cs
)
{

  cs->sortie_fmt = malloc(sizeof(PDM_writer_ensight_t));

  PDM_writer_ensight_t *_PDM_writer_ensight = (PDM_writer_ensight_t *) cs->sortie_fmt;
  _PDM_writer_ensight->f_unit_geom = NULL;
  const int restart = (int) cs->st_reprise;
  _PDM_writer_ensight->append = 0;

  for (int i = 0; i < cs->n_options; i++) {
    if(!strcmp(cs->options[i].nom, "append")) {
      if(!strcmp(cs->options[i].val, "1")) {
        _PDM_writer_ensight->append = 1;
      }
    }
  }

  _PDM_writer_ensight->ensight_case = PDM_writer_ensight_case_cree(cs->nom_sortie,
                                                                   restart,
                                                                   cs->rep_sortie,
                                                                   cs->topologie,
                                                                   &(cs->cst_global_var_tab),
                                                                   _PDM_writer_ensight->append);
  _PDM_writer_ensight->n_time_step = 0;

}

/*----------------------------------------------------------------------------
 * Particularise la liberation d'un objet CS (Cedre Sortie)
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_free
(
PDM_writer_t *cs
)
{
  _vars_close(cs);
  _geom_close(cs);
  PDM_writer_ensight_t *PDM_writer_ensight = (PDM_writer_ensight_t *) cs->sortie_fmt;
  PDM_writer_ensight_case_lib(PDM_writer_ensight->ensight_case);
  free(cs->sortie_fmt);
}

/*----------------------------------------------------------------------------
 * Particularise le debut d'increment
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_step_beg
(
PDM_writer_t *cs
)
{

  PDM_writer_ensight_t *PDM_writer_ensight = (PDM_writer_ensight_t *) cs->sortie_fmt;
  PDM_writer_ensight_case_time_step_add(PDM_writer_ensight->ensight_case,
                                        cs->physical_time);
  PDM_writer_ensight->n_time_step += 1;
}

/*----------------------------------------------------------------------------
 * Particularise la fin d'increment
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_step_end
(
PDM_writer_t *cs
)
{

  int rank = 0;
  PDM_MPI_Comm_rank(cs->pdm_mpi_comm, &rank);
  PDM_writer_ensight_t *PDM_writer_ensight = (PDM_writer_ensight_t *) cs->sortie_fmt;
  _geom_close(cs);
  _vars_close(cs);
  PDM_writer_ensight_case_write(PDM_writer_ensight->ensight_case,
                        rank);
}


/*----------------------------------------------------------------------------
 * Particularise la creation de la geometrie
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_geom_create
(
PDM_writer_geom_t *geom
)
{
  geom->geom_fmt = malloc(sizeof(PDM_writer_geom_ensight_t));
  PDM_writer_geom_ensight_t *_geom_ensight = (PDM_writer_geom_ensight_t *) geom->geom_fmt;
  _geom_ensight->num_part = geom->_cs->geom_tab->n_geom;
}


/*----------------------------------------------------------------------------
 * Particularise la creation d'une variable
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *      var          <-- Nouvelle variable
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_var_create
(
 PDM_writer_var_t *var
)
{
  var->var_fmt = malloc(sizeof(PDM_writer_var_ensight_t));
  PDM_writer_var_ensight_t *_var_ensight = (PDM_writer_var_ensight_t *) var->var_fmt;
  _var_ensight->f_unit = NULL;

  PDM_writer_ensight_t *PDM_writer_ensight = (PDM_writer_ensight_t *) var->_cs->sortie_fmt;
  PDM_writer_ensight_case_var_cree(PDM_writer_ensight->ensight_case,
                                   var->nom_var,
                                   var->dim,
                                   var->st_dep_tps,
                                   var->loc);
}


/*----------------------------------------------------------------------------
 * Particularise l'ecriture de la geometrie
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_geom_write
(
 PDM_writer_geom_t    *geom
)
{
  PDM_writer_t* _cs = (PDM_writer_t*) geom->_cs;
  PDM_writer_ensight_t *PDM_writer_ensight = (PDM_writer_ensight_t *) _cs->sortie_fmt;
  PDM_io_file_t *f_unit_geom = PDM_writer_ensight->f_unit_geom;
  /* Premier passage : Ouverture du fichier + Ecriture entête */

  if (f_unit_geom == NULL) {

    const char* geom_file_name = 
    PDM_writer_ensight_case_geo_file_name_get(PDM_writer_ensight->ensight_case);

    PDM_io_file_t *unite = NULL;
    PDM_l_num_t              ierr;
    PDM_io_fmt_t              PDM_io_fmt;

    if (_cs->fmt_fic == PDM_WRITER_FMT_BIN) {
      PDM_io_fmt = PDM_IO_FMT_BIN;
    }

    else {
      PDM_io_fmt = PDM_IO_FMT_TXT;
    }

    PDM_io_mod_t io_mod = PDM_IO_MOD_WRITE;

    if (PDM_writer_ensight->append) {
      io_mod = PDM_IO_MOD_APPEND;      
    }

    PDM_io_open(geom_file_name,
                PDM_io_fmt,
                PDM_IO_SUFF_MAN,
                "",
                PDM_IO_BACKUP_OFF,
                _cs->acces,
                io_mod,
                PDM_IO_NATIVE,
                _cs->pdm_mpi_comm,
                _cs->prop_noeuds_actifs,
                &unite,
                &ierr);

    _geom_entete_ecr(_cs,
                     unite,
                     PDM_writer_ensight);


    PDM_writer_ensight->f_unit_geom = unite;
    f_unit_geom = PDM_writer_ensight->f_unit_geom;
    PDM_writer_ensight->n_part_ecr = 0;
  }

  /* Ecriture de la part associee a la structure geom courante */

  /* Ecriture de l'entete de la part */

  PDM_writer_ensight->n_part_ecr += 1;
  _ecr_string(_cs, f_unit_geom, "part");
  int32_t _n_part = PDM_writer_ensight->n_part_ecr;
  _ecr_int(_cs, f_unit_geom, _n_part);
  if (geom->nom_geom != NULL)
    _ecr_string(_cs, f_unit_geom, geom->nom_geom);
  else
    _ecr_string(_cs, f_unit_geom, "unnamed");

  /* Calcul du nombre total de sommets */

  int n_som_proc = 0;

  const int n_part = PDM_part_mesh_nodal_n_part_get(geom->mesh_nodal);

  for (int i_part = 0; i_part < n_part; i_part++) {
    const int n_vtx = PDM_part_mesh_nodal_n_vtx_get(geom->mesh_nodal, i_part);
    n_som_proc += n_vtx;
  }

  /* Concatenation des coordonnees et ecriture */

  _ecr_string(_cs, f_unit_geom, "coordinates");

  float *coord_tmp = (float *) malloc(n_som_proc * sizeof(float));
  PDM_g_num_t *numabs_tmp =
    (PDM_g_num_t *) malloc(n_som_proc * sizeof(PDM_g_num_t));
  PDM_writer_status_t s_ecr_n_val;
  for (int idim = 0; idim < 3; idim++) {
    if (idim == 0)
      s_ecr_n_val = PDM_WRITER_ON;
    else
      s_ecr_n_val = PDM_WRITER_OFF;
    n_som_proc = 0;

    for (int i_part = 0; i_part < n_part; i_part++) {

      const int n_vtx = PDM_part_mesh_nodal_n_vtx_get(geom->mesh_nodal, i_part);
      const double *vtx = PDM_part_mesh_nodal_vtx_coord_get(geom->mesh_nodal, i_part);
      const PDM_g_num_t *numabs =
                    PDM_part_mesh_nodal_vtx_g_num_get(geom->mesh_nodal, i_part);

      for (int i = 0; i < n_vtx; i++) {
        coord_tmp[n_som_proc+i] = (float) vtx[3*i+idim];
        numabs_tmp[n_som_proc+i] = (PDM_g_num_t) numabs[i];
      }
      n_som_proc += n_vtx;
    }

    PDM_l_num_t n_comp = 1;

    _ecr_entrelace_float(_cs,
                         s_ecr_n_val,
                         f_unit_geom,
                         PDM_STRIDE_CST_INTERLACED,
                         &n_comp,
                         n_som_proc,
                         numabs_tmp,
                         coord_tmp);
  }

  free(coord_tmp);
  free(numabs_tmp);

  /* Ecriture des blocs standard */
  const int n_blocks   = PDM_part_mesh_nodal_n_section_get  (geom->mesh_nodal);
  // const int *blocks_id = PDM_part_mesh_nodal_sections_id_get(geom->mesh_nodal);

  for (int ibloc = 0; ibloc < n_blocks; ibloc++) {

    PDM_writer_elt_geom_t t_elt =
            (PDM_writer_elt_geom_t) PDM_part_mesh_nodal_section_elt_type_get(geom->mesh_nodal, ibloc);

    /* Check geometry kind */
    PDM_geometry_kind_t geom_kind = PDM_Mesh_nodal_geom_kind_from_elt_type((PDM_Mesh_nodal_elt_t) t_elt);

    if (geom_kind != geom->geom_kind) {
      PDM_error(__FILE__, __LINE__, 0, "Multiple geometry kinds are not supported\n");
    }

    /* Type de bloc */

    _ecr_string(_cs,
                f_unit_geom,
                _ensight_type_name[t_elt]);


    if (   (t_elt != PDM_WRITER_POLY_2D)
        && (t_elt != PDM_WRITER_POLY_3D)) {

      PDM_g_num_t max_loc = 0;
      PDM_g_num_t max_abs = 0;

      PDM_l_num_t n_elt_proc = 0;

      for (int i = 0; i < n_part; i++) {
        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);
        PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                          ibloc,
                                                                          i,
                                                                          PDM_OWNERSHIP_KEEP);
        n_elt_proc += n_elt;
        for (int j = 0; j < n_elt; j++) {
          max_loc = _max((PDM_g_num_t) numabs_block[j], max_loc);
        }
      }

      PDM_MPI_Allreduce(&max_loc, &max_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _cs->pdm_mpi_comm);

//      int32_t n_g_elt = (int32_t) max_abs;

      int n_comp = 0;
      switch (t_elt) {

      case PDM_WRITER_POINT    :
        n_comp = 1;
        break;
      case PDM_WRITER_BAR2     :
        n_comp = 2;
        break;
      case PDM_WRITER_TRIA3    :
        n_comp = 3;
        break;
      case PDM_WRITER_QUAD4    :
        n_comp = 4;
        break;
      case PDM_WRITER_TETRA4   :
        n_comp = 4;
        break;
      case PDM_WRITER_PYRAMID5 :
        n_comp = 5;
        break;
      case PDM_WRITER_PRISM6   :
        n_comp = 6;
        break;
      case PDM_WRITER_HEXA8    :
        n_comp = 8;
        break;

      default :
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_ensight_geom_ecr : Type d'element inconnu\n");
        abort();

      }

      /* Copie de la connectivité en numérotation absolue */

      numabs_tmp = (PDM_g_num_t *) malloc(n_elt_proc * sizeof(PDM_g_num_t));
      int32_t *connec_tmp = (int32_t *) malloc(n_elt_proc * n_comp * sizeof(int32_t));

      n_elt_proc = 0;
      for (int i = 0; i < n_part; i++) {
        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);
        PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                          ibloc,
                                                                          i,
                                                                          PDM_OWNERSHIP_KEEP);

        PDM_l_num_t  *connec;
        PDM_g_num_t  *numabs;
        int          *parent_num;
        PDM_g_num_t  *parent_entity_g_num;
        PDM_part_mesh_nodal_section_std_get(geom->mesh_nodal,
                                            ibloc,
                                            i,
                                            &connec,
                                            &numabs,
                                            &parent_num,
                                            &parent_entity_g_num,
                                            PDM_OWNERSHIP_KEEP);

        const PDM_g_num_t *g_num_vtx = PDM_part_mesh_nodal_vtx_g_num_get(geom->mesh_nodal,
                                                                          i);
        for (int j = 0; j < n_elt; j++) {
          numabs_tmp[n_elt_proc] = numabs_block[j];
          for (int k = 0; k < n_comp; k++) {
            int isom = connec[j * n_comp + k] - 1;
            int32_t isom_g = (int32_t) g_num_vtx[isom];
            connec_tmp[n_elt_proc * n_comp + k] = isom_g;
          }
          n_elt_proc += 1;
        }
      }

      /* Ecriture */

      _ecr_entrelace_int(_cs,
                         PDM_WRITER_ON,
                         f_unit_geom,
                         PDM_STRIDE_CST_INTERLACED,
                         &n_comp,
                         n_elt_proc,
                         numabs_tmp,
                         connec_tmp);

      free(numabs_tmp);
      free(connec_tmp);

    }

    else if (t_elt == PDM_WRITER_POLY_2D) {

      PDM_g_num_t max_loc = 0;
      PDM_g_num_t max_abs = 0;

      PDM_l_num_t n_elt_proc = 0;
      int l_connec = 0;

      for (int i = 0; i < n_part; i++) {

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);
        PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                          ibloc,
                                                                          i,
                                                                          PDM_OWNERSHIP_KEEP);

        PDM_l_num_t  *connec_idx;
        PDM_l_num_t  *connec;

        PDM_part_mesh_nodal_section_poly2d_get(geom->mesh_nodal,
                                               ibloc,
                                               i,
                                               &connec_idx,
                                               &connec,
                                               PDM_OWNERSHIP_KEEP);

        if (n_elt > 0) {
          n_elt_proc += n_elt;
          l_connec += connec_idx[n_elt];
          for (int j = 0; j < n_elt; j++) {
            max_loc = _max((PDM_g_num_t) numabs_block[j], max_loc);
          }
        }

      }

      PDM_MPI_Allreduce(&max_loc, &max_abs, 1,
                        PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _cs->pdm_mpi_comm);

      int32_t n_g_elt = (int32_t) max_abs;

      _ecr_int(_cs,
               f_unit_geom,
               n_g_elt);

      /* Copie de la connectivité en numérotation absolue */

      numabs_tmp = (PDM_g_num_t *) malloc(n_elt_proc * sizeof(PDM_g_num_t));
      int32_t *connec_tmp = (int32_t *) malloc(l_connec * sizeof(int32_t));
      int32_t *n_comp_tmp = (int32_t *) malloc(n_elt_proc * sizeof(int32_t));

      n_elt_proc = 0;
      l_connec = 0;

      for (int i = 0; i < n_part; i++) {

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);
        PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                          ibloc,
                                                                          i,
                                                                          PDM_OWNERSHIP_KEEP);

        PDM_l_num_t  *connec_idx;
        PDM_l_num_t  *connec;

        PDM_part_mesh_nodal_section_poly2d_get(geom->mesh_nodal,
                                               ibloc,
                                               i,
                                               &connec_idx,
                                               &connec,
                                               PDM_OWNERSHIP_KEEP);

        const PDM_g_num_t *g_num_vtx =
              PDM_part_mesh_nodal_vtx_g_num_get(geom->mesh_nodal, i);

        for (int j = 0; j < n_elt; j++) {
          numabs_tmp[n_elt_proc] = (PDM_g_num_t) numabs_block[j];
          n_comp_tmp[n_elt_proc] = (int32_t) (connec_idx[j+1] - connec_idx[j]);
          for (int k = connec_idx[j]; k < connec_idx[j+1]; k++) {
            int isom = connec[k] - 1;
            int32_t isom_g = (int32_t) g_num_vtx[isom];
            connec_tmp[l_connec++] = isom_g;
          }
          n_elt_proc += 1;
        }
      }

      /* Ecriture du nombre de sommets */

      int n_comp_cste = 1;

      _ecr_entrelace_int(_cs,
                         PDM_WRITER_OFF,
                         f_unit_geom,
                         PDM_STRIDE_CST_INTERLACED,
                         &n_comp_cste,
                         n_elt_proc,
                         numabs_tmp,
                         n_comp_tmp);

      /* Ecriture de la connectivité */

      PDM_l_num_t *n_comp_tmp2;
      if (sizeof(PDM_l_num_t) == sizeof(int32_t)) {
        n_comp_tmp2 = (PDM_l_num_t *) n_comp_tmp;
      }
      else {
        PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_ensight_geom_ecr : sizeof(int32_t) != sizeof(PDM_l_num_t)\n");
        abort();
      }

      _ecr_entrelace_int(_cs,
                         PDM_WRITER_OFF,
                         f_unit_geom,
                         PDM_STRIDE_VAR_INTERLACED,
                         n_comp_tmp2,
                         n_elt_proc,
                         numabs_tmp,
                         connec_tmp);

      free(numabs_tmp);
      free(connec_tmp);
      free(n_comp_tmp);

    }

    else {
      /* Nombre total d'éléments du bloc */

      PDM_g_num_t max_loc = 0;
      PDM_g_num_t max_abs = 0;

      PDM_l_num_t n_elt_proc = 0;
      PDM_l_num_t n_face_proc = 0;
      int l_connec = 0;
      for (int i = 0; i < n_part; i++) {

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);
        PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                          ibloc,
                                                                          i,
                                                                          PDM_OWNERSHIP_KEEP);

        PDM_l_num_t   n_face;
        PDM_l_num_t  *facvtx_idx;
        PDM_l_num_t  *facvtx;
        PDM_l_num_t  *cellfac_idx;
        PDM_l_num_t  *cellfac;
        PDM_g_num_t  *face_ln_to_gn;
        PDM_g_num_t  *cell_ln_to_gn;
        int          *parent_num;
        PDM_g_num_t  *parent_entity_g_num;
        PDM_part_mesh_nodal_section_poly3d_get(geom->mesh_nodal,
                                               ibloc,
                                               i,
                                               &n_face,
                                               &face_ln_to_gn,
                                               &facvtx_idx,
                                               &facvtx,
                                               &cell_ln_to_gn,
                                               &cellfac_idx,
                                               &cellfac,
                                               &parent_num,
                                               &parent_entity_g_num,
                                               PDM_OWNERSHIP_KEEP);

        n_elt_proc += n_elt;
        for (int k = 0; k < n_elt; k++) {
          n_face_proc += cellfac_idx[k+1] - cellfac_idx[k];
          for (int j = cellfac_idx[k]; j < cellfac_idx[k+1]; j++) {
            int ifac = PDM_ABS(cellfac[j]) - 1;//cellfac[j] - 1;
            l_connec += facvtx_idx[ifac+1] - facvtx_idx[ifac];
          }
          max_loc = _max( (PDM_g_num_t) numabs_block[k], max_loc);
        }
      }

      PDM_MPI_Allreduce(&max_loc, &max_abs, 1, PDM__PDM_MPI_G_NUM, PDM_MPI_MAX, _cs->pdm_mpi_comm);

      int32_t n_g_elt = (int32_t) max_abs;

      _ecr_int(_cs,
               f_unit_geom,
               n_g_elt);

      /* Allocation du buffer int_32 au plus grand tableau à écrire (connectivité) */

      int32_t *buff_int32 = (int32_t *) malloc(sizeof(int32_t) * l_connec);
      PDM_g_num_t *numabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_face_proc);
      int n_comp_cste = 1;

      /* Ecriture du nombre de faces de chaque cellule */

      n_elt_proc = 0;
      for (int i = 0; i < n_part; i++) {

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);
        PDM_g_num_t *numabs_block = PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                                          ibloc,
                                                                          i,
                                                                          PDM_OWNERSHIP_KEEP);

        PDM_l_num_t   n_face;
        PDM_l_num_t  *facvtx_idx;
        PDM_l_num_t  *facvtx;
        PDM_l_num_t  *cellfac_idx;
        PDM_l_num_t  *cellfac;
        PDM_g_num_t  *face_ln_to_gn;
        PDM_g_num_t  *cell_ln_to_gn;
        int          *parent_num;
        PDM_g_num_t  *parent_entity_g_num;
        PDM_part_mesh_nodal_section_poly3d_get(geom->mesh_nodal,
                                               ibloc,
                                               i,
                                               &n_face,
                                               &face_ln_to_gn,
                                               &facvtx_idx,
                                               &facvtx,
                                               &cell_ln_to_gn,
                                               &cellfac_idx,
                                               &cellfac,
                                               &parent_num,
                                               &parent_entity_g_num,
                                               PDM_OWNERSHIP_KEEP);

        for (int k = 0; k < n_elt; k++) {
          buff_int32[n_elt_proc] =  cellfac_idx[k+1] - cellfac_idx[k];
          numabs[n_elt_proc] =  (PDM_g_num_t) numabs_block[k];
          n_elt_proc += 1;
        }
      }

      _ecr_entrelace_int(_cs,
                         PDM_WRITER_OFF,
                         f_unit_geom,
                         PDM_STRIDE_CST_INTERLACED,
                         &n_comp_cste,
                         n_elt_proc,
                         numabs,
                         buff_int32);

      /* Calcul d'une numérotation absolue pour l'ensemble des faces de tous les polyèdres */

      PDM_g_num_t **numabs_face =
        (PDM_g_num_t **) malloc(sizeof(PDM_g_num_t *) * n_part);
      for (int i = 0; i < n_part; i++) {

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);

        numabs_face[i] = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_elt);
      }

      _calcul_numabs_face_poly3d(geom,
                                 ibloc,
                                 numabs_face);

      /* Ecriture du nombre de sommets de chaque face de chaque cellule */

      n_face_proc = 0;
      for (int i = 0; i < n_part; i++) {

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);

        PDM_l_num_t   n_face;
        PDM_l_num_t  *facvtx_idx;
        PDM_l_num_t  *facvtx;
        PDM_l_num_t  *cellfac_idx;
        PDM_l_num_t  *cellfac;
        PDM_g_num_t  *face_ln_to_gn;
        PDM_g_num_t  *cell_ln_to_gn;
        int          *parent_num;
        PDM_g_num_t  *parent_entity_g_num;
        PDM_part_mesh_nodal_section_poly3d_get(geom->mesh_nodal,
                                               ibloc,
                                               i,
                                               &n_face,
                                               &face_ln_to_gn,
                                               &facvtx_idx,
                                               &facvtx,
                                               &cell_ln_to_gn,
                                               &cellfac_idx,
                                               &cellfac,
                                               &parent_num,
                                               &parent_entity_g_num,
                                               PDM_OWNERSHIP_KEEP);

        for (int k = 0; k < n_elt; k++) {
          for (int j = cellfac_idx[k]; j < cellfac_idx[k+1]; j++) {
            int ifac = PDM_ABS(cellfac[j]) - 1;//;cellfac[j] - 1;
            buff_int32[n_face_proc] =
              (int32_t) (facvtx_idx[ifac+1] -
                         facvtx_idx[ifac]);
            numabs[n_face_proc] = numabs_face[i][k] + j - cellfac_idx[k];
            n_face_proc += 1;
          }
        }
      }

      for (int i = 0; i < n_part; i++) {
        free(numabs_face[i]);
      }
      free(numabs_face);

      _ecr_entrelace_int(_cs,
                         PDM_WRITER_OFF,
                         f_unit_geom,
                         PDM_STRIDE_CST_INTERLACED,
                         &n_comp_cste,
                         n_face_proc,
                         numabs,
                         buff_int32);

      /* Copie de la connectivité des faces en numérotation absolue */

      PDM_l_num_t *n_comp_tmp = (PDM_l_num_t *) malloc(sizeof(PDM_l_num_t) * n_face_proc);
      for (int i = 0; i < n_face_proc; i++) {
        n_comp_tmp[i] = (PDM_l_num_t) buff_int32[i];
      }

      l_connec = 0;
      for (int i = 0; i < n_part; i++) {

        int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal, ibloc, i);

        PDM_l_num_t   n_face;
        PDM_l_num_t  *facvtx_idx;
        PDM_l_num_t  *facvtx;
        PDM_l_num_t  *cellfac_idx;
        PDM_l_num_t  *cellfac;
        PDM_g_num_t  *face_ln_to_gn;
        PDM_g_num_t  *cell_ln_to_gn;
        int          *parent_num;
        PDM_g_num_t  *parent_entity_g_num;
        PDM_part_mesh_nodal_section_poly3d_get(geom->mesh_nodal,
                                               ibloc,
                                               i,
                                               &n_face,
                                               &face_ln_to_gn,
                                               &facvtx_idx,
                                               &facvtx,
                                               &cell_ln_to_gn,
                                               &cellfac_idx,
                                               &cellfac,
                                               &parent_num,
                                               &parent_entity_g_num,
                                               PDM_OWNERSHIP_KEEP);

        const PDM_g_num_t *g_num_vtx =
              PDM_part_mesh_nodal_vtx_g_num_get(geom->mesh_nodal, i);

        for (int k = 0; k < n_elt; k++) {
          for (int j = cellfac_idx[k]; j < cellfac_idx[k+1]; j++) {
            int ifac = PDM_ABS(cellfac[j]) - 1;//cellfac[j] - 1;
            for (int j2 = facvtx_idx[ifac]; j2 < facvtx_idx[ifac+1]; j2++) {
              int isom = facvtx[j2] - 1;
              buff_int32[l_connec] =  (int32_t) g_num_vtx[isom];
              l_connec += 1;
            }
          }
        }
      }

      _ecr_entrelace_int(_cs,
                         PDM_WRITER_OFF,
                         f_unit_geom,
                         PDM_STRIDE_VAR_INTERLACED,
                         n_comp_tmp,
                         n_face_proc,
                         numabs,
                         buff_int32);

      /* Libération mémoire */

      free(numabs);
      free(buff_int32);
      free(n_comp_tmp);

    }

  }

//  if ((PDM_writer_ensight->append) && (_cs->topologie != PDM_WRITER_TOPO_CST)) {
  if (PDM_writer_ensight->append) {
    char buff_append[81];
    strncpy(buff_append, "END TIME STEP", 80);
    buff_append[80] = '\0';  
    _ecr_string(_cs,
                f_unit_geom,
                buff_append);
  }

}


/*----------------------------------------------------------------------------
 * Particularise la libération de la geometrie
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_geom_free
(
 PDM_writer_geom_t    *geom
)
{
  if (geom->geom_fmt != NULL)
    free(geom->geom_fmt);
}


/*----------------------------------------------------------------------------
 * Particularise l'ecriture de la geometrie
 *
 * parameters :
 *      cs           <-> objet Cedre Sortie
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_var_write
(
 PDM_writer_var_t        *var
)
{
  PDM_writer_t              *cs            = (PDM_writer_t *) var->_cs;
  PDM_writer_ensight_t      *PDM_writer_ensight    = (PDM_writer_ensight_t *) cs->sortie_fmt;
  PDM_writer_var_ensight_t  *_var_ensight  = (PDM_writer_var_ensight_t *) var->var_fmt;

  /* Ouverture du fichier au premier passage  */

  const char* file_name = PDM_writer_ensight_case_var_file_name_get(PDM_writer_ensight->ensight_case,
                                                            var->nom_var);

  PDM_io_file_t *unite = NULL;
  PDM_l_num_t ierr;
  PDM_io_fmt_t PDM_io_fmt;

  if (cs->fmt_fic == PDM_WRITER_FMT_BIN) {
    PDM_io_fmt = PDM_IO_FMT_BIN;
  }

  else {
    PDM_io_fmt = PDM_IO_FMT_TXT;
  }


  PDM_io_mod_t io_mod = PDM_IO_MOD_WRITE;

  if (PDM_writer_ensight->append) {
    io_mod = PDM_IO_MOD_APPEND;      
  }


  PDM_io_open(file_name,
              PDM_io_fmt,
              PDM_IO_SUFF_MAN,
              "",
              PDM_IO_BACKUP_OFF,
              cs->acces,
              io_mod,
              PDM_IO_NATIVE,
              cs->pdm_mpi_comm,
              cs->prop_noeuds_actifs,
              &unite,
              &ierr);

  _var_ensight->f_unit = unite;

  /* Ecriture de l'entête */

  if (PDM_writer_ensight->append) {
//  if (PDM_writer_ensight->append && (var->st_dep_tps == PDM_WRITER_ON)) {

    char buff_append[81];
    strncpy(buff_append, "BEGIN TIME STEP", 80);
    buff_append[80] = '\0';  
    _ecr_string(cs,
                unite,
                buff_append);
  }

  char buff_entete[81];

  if (var->st_dep_tps == PDM_WRITER_ON) {
    for (int i = 0; i < 81; i++)
      buff_entete[i] = ' ';
    snprintf(buff_entete, 80, "%s (time values: %d, %g)",
             var->nom_var, PDM_writer_ensight->n_time_step, cs->physical_time);
  }

  else
    strncpy(buff_entete, var->nom_var, 80);

  buff_entete[80] = '\0';

  _ecr_string(cs,
              unite,
              buff_entete);

  /* Boucle sur les géométries */

  const int n_ind = cs->geom_tab->n_geom;

  for (int i1 = 0; i1 < n_ind; i1++) {
    int igeom = i1;

    PDM_writer_geom_t *geom = cs->geom_tab->geom[igeom];

    const int n_part = PDM_part_mesh_nodal_n_part_get(geom->mesh_nodal);

    if ((geom != NULL) && (var->_val[igeom] != NULL)) {
      PDM_writer_geom_ensight_t *_geom_ensight = (PDM_writer_geom_ensight_t *) geom->geom_fmt;
      int                num_part      = _geom_ensight->num_part;

      /* Ecriture de l'entête de la part */

      _ecr_string(cs, unite, "part");
      _ecr_int(cs, unite, num_part);

      if (var->loc == PDM_WRITER_VAR_VERTICES) {

        /* Ecriture des valeurs du mot "coordinates" */

        _ecr_string(cs, unite, "coordinates");

        /* Ecriture des valeurs aux sommets */

       int n_som_proc = 0;
        for (int i = 0; i < n_part; i++) {
          const int n_vertices = PDM_part_mesh_nodal_n_vtx_get(geom->mesh_nodal, i);
          n_som_proc += n_vertices;
        }

        float *buff = (float *) malloc(sizeof(float) * n_som_proc);
        PDM_g_num_t *numabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_som_proc);

        n_som_proc = 0;
        for (int i = 0; i < n_part; i++) {
          const int n_vertices = PDM_part_mesh_nodal_n_vtx_get(geom->mesh_nodal, i);
          const PDM_g_num_t *gnum = PDM_part_mesh_nodal_vtx_g_num_get(geom->mesh_nodal, i);
          for (int j = 0; j < n_vertices; j++) {
            numabs[n_som_proc++] = (PDM_g_num_t) gnum[j];
          }
        }

        PDM_writer_status_t s_ecr_n_val = PDM_WRITER_OFF ;
        for (int k = 0; k < (int) var->dim; k++) {
          s_ecr_n_val = PDM_WRITER_OFF;
          n_som_proc = 0;
          int comp_a_ecrire;
          comp_a_ecrire = k;
          if (var->dim == 9)
            comp_a_ecrire = 3 * (k % 3) + k / 3;


          for (int i = 0; i < n_part; i++) {
            const int n_vertices = PDM_part_mesh_nodal_n_vtx_get(geom->mesh_nodal, i);
            for (int j = 0; j < n_vertices; j++) {
              buff[n_som_proc++] = (float) var->_val[igeom][i][j*var->dim + comp_a_ecrire];
            }
          }

          PDM_l_num_t un = 1;
          _ecr_entrelace_float(cs,
                               s_ecr_n_val,
                               unite,
                               PDM_STRIDE_CST_INTERLACED,
                               &un,
                               n_som_proc,
                               numabs,
                               buff);

        }

        free(buff);
        free(numabs);

      }

      else if (var->loc == PDM_WRITER_VAR_ELEMENTS) {

        /* Allocation du buffer */
        const int n_blocks   = PDM_part_mesh_nodal_n_section_get  (geom->mesh_nodal);
        // const int *blocks_id = PDM_part_mesh_nodal_sections_id_get(geom->mesh_nodal);

        int n_elt_max_bloc = 0;

        for (int iblock = 0; iblock < n_blocks; iblock++) {

          int n_elt_bloc = 0;
          for (int i = 0; i < n_part; i++) {

            int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal,
                                                              iblock, i);
            n_elt_bloc += n_elt;
          }

          n_elt_max_bloc = _max_int(n_elt_max_bloc, n_elt_bloc);
        }

        float       *buff = (float *) malloc(sizeof(float) * n_elt_max_bloc);
        PDM_g_num_t *numabs = (PDM_g_num_t *) malloc(sizeof(PDM_g_num_t) * n_elt_max_bloc);

        /* Boucle sur les blocs standard */

        int *ideb = PDM_array_zeros_int(n_part);

        for (int iblock = 0; iblock < n_blocks; iblock++) {

          PDM_writer_elt_geom_t t_elt = (PDM_writer_elt_geom_t) PDM_part_mesh_nodal_section_elt_type_get(geom->mesh_nodal, iblock);

          /* Ecriture du Type de bloc */

          _ecr_string(cs,
                      unite,
                      _ensight_type_name[t_elt]);

         /* Construction de l'indirection */

          int n_val_buff = 0;
          for (int i = 0; i < n_part; i++) {
            int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal,
                                                              iblock,
                                                              i);

           PDM_g_num_t  *numabs_block =
           PDM_part_mesh_nodal_section_g_num_get(geom->mesh_nodal,
                                                 iblock,
                                                 i,
                                                 PDM_OWNERSHIP_KEEP);
            for (int j = 0; j < n_elt; j++) {
              numabs[n_val_buff++] = (PDM_g_num_t) numabs_block[j];
            }
          }

          PDM_writer_status_t s_ecr_n_val = PDM_WRITER_OFF;

          for (int k = 0; k < (int) var->dim; k++) {
            n_val_buff = 0;
            int comp_a_ecrire;
            comp_a_ecrire = k;
            if (var->dim == 9)
              comp_a_ecrire = 3 * (k % 3) + k / 3;
            for (int i = 0; i < n_part; i++) {
              int n_elt = PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal,
                                                                iblock,
                                                                i);

              for (int j = 0; j < n_elt; j++) {
                buff[n_val_buff++] = (float) (var->_val[igeom][i][(ideb[i] + j)*var->dim + comp_a_ecrire]);
              }
            }

            PDM_l_num_t un = 1;

            _ecr_entrelace_float(cs,
                                 s_ecr_n_val,
                                 unite,
                                 PDM_STRIDE_CST_INTERLACED,
                                 &un,
                                 n_val_buff,
                                 numabs,
                                 buff);
          }
          for (int i = 0; i < n_part; i++) {

            ideb[i] += PDM_part_mesh_nodal_section_n_elt_get(geom->mesh_nodal,
                                                             iblock,
                                                             i);
          }
        }

        /* Libération */

        free(ideb);
        free(buff);
        free(numabs);

      }

      else if (var->loc == PDM_WRITER_VAR_PARTICLES) {

        /* Ecriture des valeurs aux particules */

        PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_ensight_var_ecr    : Ecriture des variables aux particules indisponible\n");
        abort();

      }
    }
  }

  if (PDM_writer_ensight->append) {
//  if (PDM_writer_ensight->append && (var->st_dep_tps == PDM_WRITER_ON)) {
    char buff_append[81];
    strncpy(buff_append, "END TIME STEP", 80);
    buff_append[80] = '\0';  
    _ecr_string(cs,
               unite,
               buff_append);
  }

  int rank = 0;
  PDM_MPI_Comm_rank(cs->pdm_mpi_comm,
                &rank);
  _var_close(_var_ensight, rank);
}


/*----------------------------------------------------------------------------
 * Particularise la libération d'une variable
 *
 * parameters :
 *      var          <-> Objet variable à libérer
 *   id_geom         <-- Identificateur de l'objet geometrique
 *
 *----------------------------------------------------------------------------*/

void
PDM_writer_ensight_var_free
(
 PDM_writer_var_t     *var
)
{
  int rank = 0;
  PDM_MPI_Comm_rank(var->_cs->pdm_mpi_comm,
                &rank);
  _var_close((PDM_writer_var_ensight_t *) var->var_fmt, rank);
  if (var->var_fmt != NULL)
    free(var->var_fmt);
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
