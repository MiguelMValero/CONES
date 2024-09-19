/*============================================================================
 * Sorties Cedre
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_writer.h"
#include "pdm_writer_priv.h"
#include "pdm_writer_ensight.h"
#include "pdm_binary_search.h"
#include "pdm.h"
#include "pdm_priv.h"
#include "pdm_fortran_to_c_string.h"
#include "pdm_remove_blank.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"
#include "pdm_part_mesh_nodal_priv.h"
#include "pdm_array.h"
#include "pdm_logging.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des macros locales
 *============================================================================*/

/*============================================================================
 * Definition des types locaux
 *============================================================================*/

/*----------------------------------------------------------------------------
 * NOMBRE DE BLOCS MAX
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_DEB_ID_BLOC_STD    = 0,
  PDM_WRITER_DEB_ID_BLOC_POLY2D = 1000000,
  PDM_WRITER_DEB_ID_BLOC_POLY3D = 2000000

} PDM_writer_deb_id_bloc_t;

/*============================================================================
 * Variables globales
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Stockage des formats
 *----------------------------------------------------------------------------*/

static PDM_writer_fmt_t **fmt_tab = NULL;

/*----------------------------------------------------------------------------
 * Nombre d'objets cs stockes dans fmt_tab
 *----------------------------------------------------------------------------*/

static const int n_intern_fmt = 1;
static       int s_fmt_tab    = 0;
static       int n_fmt_tab    = 0;

/*============================================================================
 * Definition des fonctions privees
 *============================================================================*/

// Wrap mesh writer for training

void
writer_wrapper
(
 const PDM_MPI_Comm     comm,
 const char            *folder,
 const char            *file,
 int                    n_part,
 int                   *n_vtx,
 double               **coords,
 PDM_g_num_t          **vtx_ln_to_gn,
 int                   *n_elt,
 int                  **elt_vtx_idx,
 int                  **elt_vtx,
 PDM_g_num_t          **elt_ln_to_gn,
 PDM_writer_elt_geom_t  cell_t,
 int                   *n_face,
 int                  **cell_face_idx,
 int                  **cell_face,
 const char            *format,
 int                    n_elt_field,
 const char           **elt_field_name,
 double              ***elt_field_values,
 int                    n_vtx_field,
 const char           **vtx_field_name,
 double              ***vtx_field_values
)
{
  int i_rank;
  PDM_MPI_Comm_rank(comm, &i_rank);

  int is_3d_nodal = (((int) cell_t) >= 0);
  int is_2d       = ((cell_face_idx == NULL) || (cell_face == NULL)) && (!is_3d_nodal);

  PDM_writer_t *wrt = PDM_writer_create(format,
                                        PDM_WRITER_FMT_BIN,
                                        PDM_WRITER_TOPO_VARIABLE,
                                        PDM_WRITER_OFF,
                                        folder,
                                        file,
                                        comm,
                                        PDM_IO_KIND_MPI_SIMPLE,
                                        1.,
                                        NULL);

  int id_geom = PDM_writer_geom_create(wrt,
                                       file,
                                       n_part);

  int id_var_part = PDM_writer_var_create(wrt,
                                          PDM_WRITER_OFF,
                                          PDM_WRITER_VAR_SCALAR,
                                          PDM_WRITER_VAR_ELEMENTS,
                                          "i_part");

  int id_var_elt_gnum = PDM_writer_var_create(wrt,
                                              PDM_WRITER_OFF,
                                              PDM_WRITER_VAR_SCALAR,
                                              PDM_WRITER_VAR_ELEMENTS,
                                              "elt_gnum");

  // elt based field
  int *id_var_elt_field = NULL;
  if (n_elt_field > 0) {
    id_var_elt_field = malloc(sizeof(int) * n_elt_field);

    for (int i = 0; i < n_elt_field; i++) {
      id_var_elt_field[i] = PDM_writer_var_create(wrt,
                                                  PDM_WRITER_OFF,
                                                  PDM_WRITER_VAR_SCALAR,
                                                  PDM_WRITER_VAR_ELEMENTS,
                                                  elt_field_name[i]);
    }
  }

  // node based field
  int *id_var_vtx_field = NULL;
  if (n_vtx_field > 0) {
    id_var_vtx_field = malloc(sizeof(int) * n_vtx_field);

    for (int i = 0; i < n_vtx_field; i++) {
      id_var_vtx_field[i] = PDM_writer_var_create(wrt,
                                                  PDM_WRITER_OFF,
                                                  PDM_WRITER_VAR_SCALAR,
                                                  PDM_WRITER_VAR_VERTICES,
                                                  vtx_field_name[i]);
    }
  }

  PDM_writer_step_beg(wrt, 0.);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_writer_geom_coord_set(wrt,
                              id_geom,
                              i_part,
                              n_vtx       [i_part],
                              coords      [i_part],
                              vtx_ln_to_gn[i_part],
                              PDM_OWNERSHIP_USER);

    if (is_2d) {
      PDM_writer_geom_faces_facesom_add(wrt,
                                        id_geom,
                                        i_part,
                                        n_elt       [i_part],
                                        elt_vtx_idx [i_part],
                                        NULL,
                                        elt_vtx     [i_part],
                                        elt_ln_to_gn[i_part]);
    } else {
      if (is_3d_nodal) {
        int id_bloc = PDM_writer_geom_bloc_add(wrt,
                                               id_geom,
                                               cell_t,
                                               PDM_OWNERSHIP_USER);
        PDM_writer_geom_bloc_std_set(wrt,
                                     id_geom,
                                     id_bloc,
                                     i_part,
                                     n_elt       [i_part],
                                     elt_vtx     [i_part],
                                     elt_ln_to_gn[i_part]);
      } else {
        PDM_writer_geom_cell3d_cellface_add(wrt,
                                            id_geom,
                                            i_part,
                                            n_elt         [i_part],
                                            n_face        [i_part],
                                            elt_vtx_idx   [i_part],
                                            NULL,
                                            elt_vtx       [i_part],
                                            cell_face_idx [i_part],
                                            NULL,
                                            cell_face     [i_part],
                                            elt_ln_to_gn  [i_part]);

      }
    }
  }

  PDM_writer_geom_write(wrt, id_geom);

  for (int i_part = 0; i_part < n_part; i_part++) {
    PDM_real_t *val_part = malloc(sizeof(PDM_real_t) * n_elt[i_part]);
    PDM_real_t *val_gnum = malloc(sizeof(PDM_real_t) * n_elt[i_part]);

    for (int i_face = 0; i_face < n_elt[i_part]; i_face++) {
      val_part[i_face] = i_rank*n_part + i_part;
      val_gnum[i_face] = elt_ln_to_gn[i_part][i_face];
    }

    PDM_writer_var_set(wrt,
                       id_var_part,
                       id_geom,
                       i_part,
                       val_part);
    free(val_part);

    PDM_writer_var_set(wrt,
                       id_var_elt_gnum,
                       id_geom,
                       i_part,
                       val_gnum);
    free(val_gnum);
  }


  PDM_writer_var_write(wrt, id_var_part);
  PDM_writer_var_write(wrt, id_var_elt_gnum);

  // elt based field
  if (n_elt_field > 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_real_t *val = malloc(sizeof(PDM_real_t) * n_elt[i_part]);

      for (int i_field = 0; i_field < n_elt_field; i_field++) {
        for (int i_elt = 0; i_elt < n_elt[i_part]; i_elt++) {
          val[i_elt] = elt_field_values[i_field][i_part][i_elt];
        }

        PDM_writer_var_set(wrt,
                           id_var_elt_field[i_field],
                           id_geom,
                           i_part,
                           val);
      }
      free(val);
    }

    for (int i_field = 0; i_field < n_elt_field; i_field++) {
      PDM_writer_var_write(wrt, id_var_elt_field[i_field]);
    }
  }

  // vtx based field
  if (n_vtx_field > 0) {
    for (int i_part = 0; i_part < n_part; i_part++) {
      PDM_real_t *val = malloc(sizeof(PDM_real_t) * n_vtx[i_part]);

      for (int i_field = 0; i_field < n_vtx_field; i_field++) {
        for (int i_vtx = 0; i_vtx < n_vtx[i_part]; i_vtx++) {
          val[i_vtx] = vtx_field_values[i_field][i_part][i_vtx];
        }

        PDM_writer_var_set(wrt,
                           id_var_vtx_field[i_field],
                           id_geom,
                           i_part,
                           val);
      }
      free(val);
    }

    for (int i_field = 0; i_field < n_vtx_field; i_field++) {
      PDM_writer_var_write(wrt, id_var_vtx_field[i_field]);
    }
  }

  PDM_writer_step_end(wrt);

  if (n_elt_field > 0) {
    free(id_var_elt_field);
  }

  if (n_vtx_field > 0) {
    free(id_var_vtx_field);
  }

  PDM_writer_free(wrt);

}

/**
 *
 * \brief Create a \ref _PDM_writer_geom_tab_t object
 *
 * \param [in]  size   Initial size of the array of geometries
 *
 * \return    Pointer to a new \ref _PDM_writer_geom_tab_t object
 *
 */

static _PDM_writer_geom_tab_t *
_pdm_writer_geom_tab_create
(
 const int size
 )
{
  _PDM_writer_geom_tab_t *geom_tab = (_PDM_writer_geom_tab_t *) malloc(sizeof(_PDM_writer_geom_tab_t));

  geom_tab->n_geom = 0;
  geom_tab->s_geom = size;

  geom_tab->geom = (PDM_writer_geom_t **) malloc(sizeof(PDM_writer_geom_t *) * geom_tab->s_geom);
  for (int i = 0; i < geom_tab->s_geom; i++) {
    geom_tab->geom[i] = NULL;
  }

  return geom_tab;
}


/**
 *
 * \brief Add a geometry
 *
 * \param [in] geom_tab   Pointer to \ref _PDM_writer_geom_tab_t object
 * \param [in] geom       Pointer to \ref PDM_writer_geom_t object
 *
 */

static int
_pdm_writer_geom_tab_add
(
 _PDM_writer_geom_tab_t *geom_tab,
 PDM_writer_geom_t      *geom
 )
{
  assert (geom_tab != NULL);

  if (geom_tab->n_geom >= geom_tab->s_geom) {
    geom_tab->s_geom = PDM_MAX(2*geom_tab->s_geom, geom_tab->n_geom+1);
    geom_tab->geom = (PDM_writer_geom_t **) realloc(geom_tab->geom, sizeof(PDM_writer_geom_t *) * geom_tab->s_geom);

    for (int i = geom_tab->n_geom+1; i < geom_tab->s_geom; i++) {
      geom_tab->geom[i] = NULL;
    }
  }

  int id_geom = geom_tab->n_geom;
  geom_tab->n_geom++;

  geom_tab->geom[id_geom] = geom;


  return id_geom;
}


/**
 *
 * \brief Create a \ref _PDM_writer_var_tab_t object
 *
 * \param [in]  size   Initial size of the array of variables
 *
 * \return    Pointer to a new \ref _PDM_writer_var_tab_t object
 *
 */

static _PDM_writer_var_tab_t *
_pdm_writer_var_tab_create
(
 const int size
 )
{
  _PDM_writer_var_tab_t *var_tab = (_PDM_writer_var_tab_t *) malloc(sizeof(_PDM_writer_var_tab_t));

  var_tab->n_var = 0;
  var_tab->s_var = size;

  var_tab->var = (PDM_writer_var_t **) malloc(sizeof(PDM_writer_var_t *) * var_tab->s_var);
  for (int i = 0; i < var_tab->s_var; i++) {
    var_tab->var[i] = NULL;
  }

  return var_tab;
}


/**
 *
 * \brief Add a variable
 *
 * \param [in] var_tab   Pointer to \ref _PDM_writer_var_tab_t object
 * \param [in] var       Pointer to \ref PDM_writer_var_t object
 *
 */

static int
_pdm_writer_var_tab_add
(
 _PDM_writer_var_tab_t *var_tab,
 PDM_writer_var_t      *var
 )
{
  assert (var_tab != NULL);

  if (var_tab->n_var >= var_tab->s_var) {
    var_tab->s_var = PDM_MAX(2*var_tab->s_var, var_tab->n_var+1);
    var_tab->var = (PDM_writer_var_t **) realloc(var_tab->var, sizeof(PDM_writer_var_t *) * var_tab->s_var);

    for (int i = var_tab->n_var+1; i < var_tab->s_var; i++) {
      var_tab->var[i] = NULL;
    }
  }

  int id_var = var_tab->n_var;
  var_tab->n_var++;

  var_tab->var[id_var] = var;


  return id_var;
}


/**
 *
 * \brief Create a \ref _PDM_writer_name_map_tab_t object
 *
 * \param [in]  size   Initial size of the array of name maps
 *
 * \return    Pointer to a new \ref _PDM_writer_name_map_tab_t object
 *
 */

static _PDM_writer_name_map_tab_t *
_pdm_writer_name_map_tab_create
(
 const int size
 )
{
  _PDM_writer_name_map_tab_t *name_map_tab = (_PDM_writer_name_map_tab_t *) malloc(sizeof(_PDM_writer_name_map_tab_t));

  name_map_tab->n_name_map = 0;
  name_map_tab->s_name_map = size;

  name_map_tab->name_map = (PDM_writer_name_map_t **) malloc(sizeof(PDM_writer_name_map_t *) * name_map_tab->s_name_map);
  for (int i = 0; i < name_map_tab->s_name_map; i++) {
    name_map_tab->name_map[i] = NULL;
  }

  return name_map_tab;
}


/**
 *
 * \brief Add a name map
 *
 * \param [in] name_map_tab   Pointer to \ref _PDM_writer_name_map_tab_t object
 * \param [in] name_map       Pointer to \ref PDM_writer_name_map_t object
 *
 */

static int
_pdm_writer_name_map_tab_add
(
 _PDM_writer_name_map_tab_t *name_map_tab,
 PDM_writer_name_map_t      *name_map
 )
{
  assert (name_map_tab != NULL);

  if (name_map_tab->n_name_map >= name_map_tab->s_name_map) {
    name_map_tab->s_name_map = PDM_MAX(2*name_map_tab->s_name_map, name_map_tab->n_name_map+1);
    name_map_tab->name_map = (PDM_writer_name_map_t **) realloc(name_map_tab->name_map, sizeof(PDM_writer_name_map_t *) * name_map_tab->s_name_map);

    for (int i = name_map_tab->n_name_map+1; i < name_map_tab->s_name_map; i++) {
      name_map_tab->name_map[i] = NULL;
    }
  }

  int id_name_map = name_map_tab->n_name_map;
  name_map_tab->n_name_map++;

  name_map_tab->name_map[id_name_map] = name_map;


  return id_name_map;
}


/**
 * \brief Initialize a \ref PDM_writer_geom_t object
 *
 * \param [in] geom            Pointer to \ref PDM_writer_geom_t object
 * \param [in] n_part          Number of partition
 * \param [in] comm            MPI communicator
 *
 */

static void
_geom_init
(
PDM_writer_geom_t  *geom,
const int           n_part,
const PDM_MPI_Comm  comm
)
{
  geom->nom_geom = NULL;

  int mesh_dimension = 3; // ?

  geom->_mesh_nodal = PDM_part_mesh_nodal_create(mesh_dimension, n_part, comm);
  geom->mesh_nodal = geom->_mesh_nodal; 

  geom->geom_fmt       = NULL;
  geom->_cs            = NULL;
  geom->pdm_mpi_comm   = comm;

  geom->s_section = 10;
  geom->section_owner = malloc(sizeof(PDM_ownership_t) * geom->s_section);

  geom->n_part = n_part;
  geom->_face_vtx_idx  = malloc(sizeof(int *) * n_part);
  geom->_cell_face_idx = malloc(sizeof(int *) * n_part);
  for (int i = 0; i < n_part; i++) {
    geom->_face_vtx_idx [i] = NULL;
    geom->_cell_face_idx[i] = NULL;
  }

  geom->geom_kind = PDM_GEOMETRY_KIND_MAX; // not yet defined
}


/**
 * \brief Initialize a \ref PDM_writer_var_t object
 *
 * \param [in] var            Pointer to \ref PDM_writer_var_t object
 *
 */

static void
_var_init
(
PDM_writer_var_t *var
)
{
  var->nom_var    = NULL;                    /* Nom de la geometrie */
  var->st_dep_tps = PDM_WRITER_OFF;          /* Variable en temps */
  var->dim        = PDM_WRITER_VAR_CST;     /* Dimension de la variable */
  var->loc        = PDM_WRITER_VAR_VERTICES;  /* Dimension de la variable */
  var->_val       = NULL;                    /* Valeurs de la variable */
  var->var_fmt    = NULL;                    /* Description propre au format fmt */
  var->_cs        = NULL;
}


/**
 *
 * \brief Parse la chaine options pour construire la structure CS correspondante
 *
 * \param [in]  options_str           Options_str : chaine provenant de cs_cree
 * \param [out] n_options             Nombre d'options trouvees
 * \param [out] options               Liste des options parsees
 *
 */

static void
_parse_options
(
 const char           *options_str,
 int                  *n_options,
 PDM_writer_option_t **options
)
{

  if (options_str == NULL) {
    *n_options = 0;
    *options = NULL;
  }

  char *_options_str  = malloc (sizeof(char) * (strlen(options_str) + 1));
  strcpy(_options_str, options_str);

  *n_options = 0;
  char *str2 = _options_str;
  char *pch;

  do {
    pch = strtok (str2,"=");
    str2 = NULL;
    if (pch != NULL) {
      pch = strtok (str2, ":");
      if (pch == NULL) {
        PDM_error(__FILE__, __LINE__, 0, "CS_cree : Erreur dans le parsing des options specifiques :"
                 "verifier les separateurs dans la chaine 'options'\n");
        exit(1);
      }
      *n_options += 1;
    }
  } while (pch != NULL);

  strcpy(_options_str, options_str);
  str2 = _options_str;
  *options = malloc (sizeof(PDM_writer_option_t) * (*n_options));
  PDM_writer_option_t *_curr = *options;

  do {
    pch = strtok (str2,"=");
    str2 = NULL;
    if (pch != NULL) {
      _curr->nom = PDM_remove_blank (pch);
      pch = strtok (str2, ":");
      if (pch != NULL) {
        _curr->val = PDM_remove_blank (pch);
      }
    }
    _curr += 1;
  } while (pch != NULL);


  free (_options_str);
}

/**
 *
 * \brief Load built-in formats
 *
 */

static void
_load_intern_fmt (void)
{
  if (fmt_tab != NULL) {
    return;
  }

  s_fmt_tab = 2 * n_intern_fmt;
  n_fmt_tab = 0;
  fmt_tab = (PDM_writer_fmt_t **) malloc (sizeof(PDM_writer_fmt_t *) * s_fmt_tab);

  /* Ensight */

  PDM_writer_fmt_t *fmt = malloc (sizeof(PDM_writer_fmt_t));
  fmt->name = malloc (sizeof(int) * 8);
  strcpy (fmt->name, "Ensight");
  fmt->create_fct       = PDM_writer_ensight_create;
  fmt->free_fct         = PDM_writer_ensight_free;
  fmt->beg_step_fct     = PDM_writer_ensight_step_beg;
  fmt->end_step_fct     = PDM_writer_ensight_step_end;
  fmt->geom_create_fct  = PDM_writer_ensight_geom_create;
  fmt->geom_write_fct   = PDM_writer_ensight_geom_write;
  fmt->geom_free_fct    = PDM_writer_ensight_geom_free;
  fmt->var_create_fct   = PDM_writer_ensight_var_create;
  fmt->var_write_fct    = PDM_writer_ensight_var_write;
  fmt->var_free_fct     = PDM_writer_ensight_var_free;

  fmt_tab[n_fmt_tab++] = fmt;
}

/*============================================================================
 * Definition des fonctions publiques
 *============================================================================*/

/**
 *
 * \brief Cree un objet CS (Cedre Sortie) et retoure un pointeur sur cet objet
 *
 * \param [in] fmt                  Format de sortie
 * \param [in] fmt_fic              Binary or ASCII
 * \param [in] topologie            Indique le maillage est mobile ou non
 * \param [in] st_reprise           Complete les sorties des calculs precedents en reprise
 * \param [in] rep_sortie           Repertoire de sortie
 * \param [in] nom_sortie           Nom de la sortie
 * \param [in] pdm_mpi_com          Communicateur MSG
 * \param [in] acces                Type d'acces
 * \param [in] prop_noeuds_actifs   Proportion des noeuds actifs dans les acces au fichier
 *                                    *  -1 : tous les processus actifs
 *                                    *   1 : un processus par noeud
 *                                    * 0 < val < 1 : un processus par noeud actif
 * \param [in] options              Options complementaires propres au format sous
 *                                 la forme ("nom_1 = val_1 : ... : nom_n = val_n")
 *
 * \return   Pointer to \ref PDM_writer object
 *
 */

PDM_writer_t *
PDM_writer_create
(
const char                   *fmt,
const PDM_writer_fmt_fic_t    fmt_fic,
const PDM_writer_topology_t  topologie,
const PDM_writer_status_t     st_reprise,
const char                   *rep_sortie,
const char                   *nom_sortie,
const PDM_MPI_Comm            pdm_mpi_comm,
const PDM_io_kind_t          acces,
const double                  prop_noeuds_actifs,
const char                   *options
)
{

  if (fmt_tab == NULL) {
    _load_intern_fmt();
  }

  /* Look for fmt */

  int fmt_id = -1;

  for (int i = 0; i < n_fmt_tab; i++) {
    PDM_writer_fmt_t *fmt_ptr = fmt_tab[i];
    if (!strcmp(fmt, fmt_ptr->name)) {
      fmt_id = i;
      break;
    }
  }

  if (fmt_id == -1) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_create : unknown format '%s'", fmt);
    abort();
  }

  /* Mise a jour du tableau de stockage */

  PDM_io_mkdir(rep_sortie);


  /* Creation du repertoire de sortie si non cree */

#ifdef _WIN32
  mkdir(rep_sortie);
#else
  mkdir(rep_sortie, 0775);
#endif

  /* Allocation de la structure PDM_writer_t */

  PDM_writer_t *cs = (PDM_writer_t *) malloc(sizeof(PDM_writer_t));


  /* Initialisation de la structure PDM_writer_t */

  cs->fmt_id      = fmt_id;     /* Format de sortie */
  cs->fmt_fic     = fmt_fic;    /* Format du fichier de sortie */
  cs->topologie   = topologie;  /* Type de toplogie du maillage */
  cs->st_reprise  = st_reprise; /* Reprise d'une sortie existante */

  size_t l_rep_sortie = strlen(rep_sortie);
  cs->rep_sortie = (char *) malloc(sizeof(char) * (l_rep_sortie + 1));
  strcpy(cs->rep_sortie, rep_sortie);  /* Nom du repertoire de sortie */
  // Gestion des options

  cs->n_options = 0;
  cs->options = NULL;
  if (options != NULL) {
    _parse_options (options, &(cs->n_options), &(cs->options));
  }

  cs->cst_global_var_tab.var = NULL;
  cs->cst_global_var_tab.n_var = 0;
  cs->cst_global_var_tab.s_var = 0;

  size_t l_nom_sortie = strlen(nom_sortie);
  cs->nom_sortie = (char *) malloc(sizeof(char) * (l_nom_sortie + 1));
  strcpy(cs->nom_sortie, nom_sortie);  /* Nom de la sortie */

  cs->pdm_mpi_comm    = pdm_mpi_comm;  /* Communicateur MPI */
  cs->sortie_fmt  = NULL;      /* Pointeur sur l'objet sortie de format fmt */

  cs->var_tab     = NULL;      /* Tableau des variables */
  cs->geom_tab    = NULL;      /* Tableau des geometries */
  cs->physical_time = 0;       /* Temps physique de simulation */
  cs->acces       = acces;
  cs->prop_noeuds_actifs = prop_noeuds_actifs;
  cs->name_map_tab = NULL;

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->create_fct != NULL) {
    (fmt_ptr->create_fct) (cs);
  }

  cs->is_there_open_step = 0;

  return cs;

}

/**
 * \brief Libere un objet CS (Cedre Sortie) et retourne un pointeur NULL si pas d'erreur
 *
 * \param [in] cs    Pointer to \ref PDM_writer object
 *
 */

void
PDM_writer_free
(
 PDM_writer_t *cs
)
{
  if (cs == NULL) {
    return;
  }

  PDM_writer_step_end (cs);
  
  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->free_fct != NULL) {
    (fmt_ptr->free_fct) (cs);
  }

  /* Liberation des differents elements de la structure */

  free(cs->rep_sortie);
  free(cs->nom_sortie);

  /* Liberation des variables */

  if (cs->var_tab != NULL) {
    if (cs->var_tab->var != NULL) {
      for (int i = 0; i < cs->var_tab->n_var; i++) {
        if (cs->var_tab->var[i] != NULL) {
          PDM_writer_var_free(cs, i);
        }
      }

      free(cs->var_tab->var);
      cs->var_tab->var = NULL;
    }

    free(cs->var_tab);
    cs->var_tab = NULL;
  }

  if (cs->cst_global_var_tab.var != NULL) {
    for (int i = 0; i < cs->cst_global_var_tab.s_var; i++) {
      if (cs->cst_global_var_tab.var[i] != NULL) {
        if (cs->cst_global_var_tab.var[i]->nom_var != NULL) {
          free (cs->cst_global_var_tab.var[i]->nom_var);
          cs->cst_global_var_tab.var[i]->nom_var = NULL;
        } 
        free (cs->cst_global_var_tab.var[i]);
        cs->cst_global_var_tab.var[i] = NULL;
      }
    }

    free(cs->cst_global_var_tab.var);
    cs->cst_global_var_tab.var = NULL;
  }


  if (cs->options != NULL) {
    for (int i = 0; i < cs->n_options; i++) {
      if ((cs->options[i]).nom != NULL) {
        free ((cs->options[i]).nom);
      }
      if ((cs->options[i]).val != NULL) {
        free ((cs->options[i]).val);
      }
    }
    free (cs->options);
  }


  /* Liberation de la geometrie */

  if (cs->geom_tab != NULL) {
    if (cs->geom_tab->geom != NULL) {
      for (int i = 0; i < cs->geom_tab->n_geom; i++) {
        if (cs->geom_tab->geom[i] != NULL) {
          PDM_writer_geom_free(cs, i);
        }
      }

      free(cs->geom_tab->geom);
      cs->geom_tab->geom = NULL;
    }

    free(cs->geom_tab);
    cs->geom_tab = NULL;
  }

  if (cs->name_map_tab != NULL) {
    if (cs->name_map_tab->name_map != NULL) {
      for (int i = 0; i < cs->name_map_tab->n_name_map; i++) {
        if (cs->name_map_tab->name_map[i] != NULL) {
          free(cs->name_map_tab->name_map[i]->public_name);
          free(cs->name_map_tab->name_map[i]->private_name);
        }
      }

      free(cs->name_map_tab->name_map);
      cs->name_map_tab->name_map = NULL;
    }

    free(cs->name_map_tab);
    cs->name_map_tab = NULL;
  }

  /* Liberation de la structure */

  free(cs);
  cs = NULL;
}


/**
 * \brief Debut d'increment
 *
 * \param [in] cs             Pointer to \ref PDM_writer object
 * \param [in] physical_time  Temps
 */

void
PDM_writer_step_beg
(
 PDM_writer_t  *cs,
 const double   physical_time
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }
 
  if (cs->is_there_open_step) {
    PDM_error (__FILE__, __LINE__, 0, "Error PDM_writer_step_beg : A step is already open\n");
  }

  cs->physical_time = physical_time;

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->beg_step_fct != NULL) {
    (fmt_ptr->beg_step_fct) (cs);
  }

  cs->is_there_open_step = 1;

}



/**
 * \brief Is there a open step
 *
 * \param [in] cs             Pointer to \ref PDM_writer object
 *
 * \return   Flag which indicates if a step is open
 */

int 
PDM_writer_is_open_step
(
 PDM_writer_t  *cs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  return cs->is_there_open_step;  
}



/**
 * \brief Fin d'increment
 *
 * \param [in] cs             Pointer to \ref PDM_writer object
 *
 */

void
PDM_writer_step_end
(
 PDM_writer_t  *cs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Appel de la fonction complementaire propre au format */

  if (cs->is_there_open_step) {

    PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

    if (fmt_ptr->end_step_fct != NULL) {
      (fmt_ptr->end_step_fct) (cs);
    }

  }
  
  cs->is_there_open_step = 0;

}

/**
 * \brief Cree une nouvelle geometrie dans l'objet CS (Cedre Sortie)
 *
 * \param [in]  cs                Pointer to \ref PDM_writer object
 * \param [in]  nom_geom          Nom de l'objet geometrique
 *
 * \return   Identificateur de l'objet geom dans cs
 *
 */

int
PDM_writer_geom_create
(
 PDM_writer_t               *cs,
 const char                 *nom_geom,
 const int                   n_part
)
{

  if (n_part <= 0) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_geom_create : Le nombre de partition doit etre >\n"
                    "                      Ajuster le communicateur MPI ou\n"
                    "                      Creer un sous-domaine avec 0 element\n");
  }

  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->geom_tab == NULL) {
    cs->geom_tab = _pdm_writer_geom_tab_create(4);
  }

  /* Allocation de la structure PDM_writer_geom_t */

  PDM_writer_geom_t *geom = (PDM_writer_geom_t *) malloc(sizeof(PDM_writer_geom_t));

  int id_geom = _pdm_writer_geom_tab_add(cs->geom_tab, geom);


  /* Initialisation de la structure PDM_writer_geom_t */

  _geom_init(geom, n_part, cs->pdm_mpi_comm);

  geom->_cs = cs;
  geom->pdm_mpi_comm = cs->pdm_mpi_comm;
  size_t l_nom_geom = strlen(nom_geom);
  geom->nom_geom = (char *) malloc(sizeof(char) * (l_nom_geom + 1));
  strcpy(geom->nom_geom, nom_geom);  /* Nom de la geometrie */

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->geom_create_fct != NULL) {
    (fmt_ptr->geom_create_fct) (geom);
  }

  return id_geom;
}

//-->>
int
PDM_writer_geom_create_from_mesh_nodal
(
 PDM_writer_t              *cs,
 const char                *nom_geom,
 PDM_part_mesh_nodal_t     *mesh
)
{
  /* Erreur si le decoupage des polygones ou polyedres est choisi */

  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->geom_tab == NULL) {
    cs->geom_tab = _pdm_writer_geom_tab_create(4);
  }

  /* Allocation de la structure PDM_writer_geom_t */

  PDM_writer_geom_t *geom = (PDM_writer_geom_t *) malloc(sizeof(PDM_writer_geom_t));

  int id_geom = _pdm_writer_geom_tab_add(cs->geom_tab, geom);

  /* Initialisation de la structure PDM_writer_geom_t */

  //_geom_init(geom, n_part, cs->pdm_mpi_comm);
  geom->nom_geom    = NULL;
  geom->_mesh_nodal = NULL;
  geom->mesh_nodal  = mesh;
  geom->geom_fmt    = NULL;

  geom->_cs = cs;
  geom->pdm_mpi_comm = cs->pdm_mpi_comm;
  size_t l_nom_geom = strlen(nom_geom);
  geom->nom_geom = (char *) malloc(sizeof(char) * (l_nom_geom + 1));
  strcpy(geom->nom_geom, nom_geom);  /* Nom de la geometrie */

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->geom_create_fct != NULL) {
    (fmt_ptr->geom_create_fct) (geom);
  }

  geom->s_section = 10;
  geom->section_owner = malloc(sizeof(PDM_ownership_t) * geom->s_section);

  geom->n_part = PDM_part_mesh_nodal_n_part_get(mesh);
  geom->_face_vtx_idx  = NULL;
  geom->_cell_face_idx = NULL;

  /* Infer geometry kind from part mesh_nodal */
  PDM_geometry_kind_t geom_kind_min = PDM_GEOMETRY_KIND_MAX;
  switch (mesh->mesh_dimension) {
  case 3:
    geom_kind_min = PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  case 2:
    geom_kind_min = PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 1:
    geom_kind_min = PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 0:
    geom_kind_min = PDM_GEOMETRY_KIND_CORNER;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid mesh_dimension %d\n", mesh->mesh_dimension);
  }

  for (PDM_geometry_kind_t geom_kind = geom_kind_min; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    int n_section = PDM_part_mesh_nodal_n_section_in_geom_kind_get(geom->mesh_nodal,
                                                                   geom_kind);

    if (n_section > 0) {
      geom->geom_kind = geom_kind;
      break;
    }
  }

  return id_geom;
}
//<<--


void
PDM_writer_geom_set_from_mesh_nodal
(
 PDM_writer_t              *cs,
 const int                  id_geom,
 PDM_part_mesh_nodal_t     *mesh
)
{
  /* Erreur si le decoupage des polygones ou polyedres est choisi */

  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  /* Initialisation de la structure PDM_writer_geom_t */

  //_geom_init(geom, n_part, cs->pdm_mpi_comm);
  geom->mesh_nodal  = mesh;

  geom->s_section = 10;
  if (geom->section_owner != NULL) {
    free(geom->section_owner);
    geom->section_owner = NULL;
  }
  geom->section_owner = malloc(sizeof(PDM_ownership_t) * geom->s_section);

  geom->n_part = PDM_part_mesh_nodal_n_part_get(mesh);
  geom->_face_vtx_idx  = NULL;
  geom->_cell_face_idx = NULL;

  /* Infer geometry kind from part mesh_nodal */
  PDM_geometry_kind_t geom_kind_min = PDM_GEOMETRY_KIND_MAX;
  switch (mesh->mesh_dimension) {
  case 3:
    geom_kind_min = PDM_GEOMETRY_KIND_VOLUMIC;
    break;
  case 2:
    geom_kind_min = PDM_GEOMETRY_KIND_SURFACIC;
    break;
  case 1:
    geom_kind_min = PDM_GEOMETRY_KIND_RIDGE;
    break;
  case 0:
    geom_kind_min = PDM_GEOMETRY_KIND_CORNER;
    break;
  default:
    PDM_error(__FILE__, __LINE__, 0, "Invalid mesh_dimension %d\n", mesh->mesh_dimension);
  }

  for (PDM_geometry_kind_t geom_kind = geom_kind_min; geom_kind < PDM_GEOMETRY_KIND_MAX; geom_kind++) {
    int n_section = PDM_part_mesh_nodal_n_section_in_geom_kind_get(geom->mesh_nodal,
                                                                   geom_kind);

    if (n_section > 0) {
      geom->geom_kind = geom_kind;
      break;
    }
  }

}



/**
 * \brief Definition des coordonnees de la partition courante
 *
 * \param [in] cs        Pointer to \ref PDM_writer object
 * \param [in] id_geom   Identificateur de l'objet geometrique
 * \param [in] id_part   Indice de partition
 * \param [in] n_som     Nombre de sommets de la partition
 * \param [in] coords    Coordonnes des sommets
 * \param [in] numabs    Numerotation absolue des sommets
 *
 */

void
PDM_writer_geom_coord_set
(
 PDM_writer_t      *cs,
 const int          id_geom,
 const int          id_part,
 const int          n_som,
 const PDM_real_t  *coords,
 const PDM_g_num_t *numabs,
 const PDM_ownership_t owner
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_part_mesh_nodal_coord_set(geom->_mesh_nodal, id_part, n_som, coords, numabs, owner);

  if (0 == 1) {
    printf("n_vtx : %d\n", n_som);
    for (int i = 0; i < n_som; i++) {
      printf ("%d "PDM_FMT_G_NUM" : %12.5e %12.5e %12.5e\n", i+1, numabs[i],
              coords[3*i], coords[3*i+1], coords[3*i+2]);
    }
  }

}


/**
 * \brief Definition des coordonnees des sommets de la partition courante
 * a partir d'un ensemble parent
 *
 *
 * \param [in] cs               Pointer to \ref PDM_writer object
 * \param [in] id_geom          Identificateur de l'objet geometrique
 * \param [in] id_part          Indice de partition
 * \param [in] n_som            Nombre de sommets de la partition
 * \param [in] n_som_parent     Nombre de sommets parent
 * \param [in] numabs           Numerotation absolue des sommets (size = n_som)
 * \param [in] num_parent       Numerotation des sommets dans la numerotation parente (size = n_som)
 * \param [in] coords_parent    Coordonnes des sommets parents (size = 3 * n_som_parent)
 * \param [in] numabs_parent    Numerotation absolue des sommets parents (size = n_som_parent)
 *
 */

void
PDM_writer_geom_coord_from_parent_set
(
 PDM_writer_t      *cs,
 const int          id_geom,
 const int          id_part,
 const int          n_som,
 const int          n_som_parent,
 const PDM_g_num_t *numabs,
 const int         *num_parent,
 const PDM_real_t  *coords_parent,
 const PDM_g_num_t *numabs_parent,
 const PDM_ownership_t ownership
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
  }

  PDM_part_mesh_nodal_coord_from_parent_set(geom->_mesh_nodal,
                                            id_part,
                                            n_som,
                                            n_som_parent,
                                            numabs,
                                            num_parent,
                                            coords_parent,
                                            numabs_parent,
                                            ownership);
}

/**
 * \brief Ajout d'un bloc d'elements d'un type donne
 *
 * \param [in] cs             Pointer to \ref PDM_writer object
 * \param [in] id_geom        Identificateur de l'objet geometrique
 * \param [in] t_elt          Type d'element
 *
 * \return   Identificateur du bloc
 *
 */

int
PDM_writer_geom_bloc_add
(
 PDM_writer_t                *cs,
 const int                    id_geom,
 const PDM_writer_elt_geom_t  t_elt,
 const PDM_ownership_t        owner
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  /* Check geom_kind coherence */
  PDM_geometry_kind_t geom_kind = PDM_Mesh_nodal_geom_kind_from_elt_type((PDM_Mesh_nodal_elt_t) t_elt);

  if (geom->geom_kind == PDM_GEOMETRY_KIND_MAX) {
    geom->geom_kind = geom_kind;
  }
  else {
    if (geom_kind != geom->geom_kind) {
      PDM_error(__FILE__, __LINE__, 0,
                "Current geometry kind (%d) cannot contain element of type %d\n",
                (int) geom->geom_kind, (int) t_elt);
    }
  }

  int id_block = PDM_part_mesh_nodal_section_add(geom->_mesh_nodal,
                                                 (PDM_Mesh_nodal_elt_t) t_elt);


  if (id_block >= geom->s_section) {
    geom->s_section = PDM_MAX(geom->s_section, id_block);
    geom->section_owner = realloc(geom->section_owner, sizeof(PDM_ownership_t) * geom->s_section);
  }
  geom->section_owner[id_block] = owner;

  return id_block;

}


/**
 * \brief Ajout d'un bloc d'elements d'un type donne dans la partition courante
 *
 *  - PDM_writer_POINT :
 *
 *   1 x
 *
 *  - PDM_writer_BAR2 :
 *
 *   1 x-------x 2
 *
 *  - PDM_writer_TRIA3 :
 *
 *   1 x-------x 3
 *      \     /
 *       \   /
 *        \ /
 *         x 2
 *
 *  - PDM_writer_QUAD4 :
 *
 *      4 x-------x 3
 *       /       /
 *      /       /
 *   1 x-------x2
 *
 *   - PDM_writer_TETRA4 :
 *
 *         x 4
 *        /|\
 *       / | \
 *      /  |  \
 *   1 x- -|- -x 3
 *      \  |  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *   - PDM_writer_PYRAMID5 :
 *
 *          5 x
 *           /|\
 *          //| \
 *         // |  \
 *      4 x/--|---x 3
 *       //   |  /
 *      //    | /
 *   1 x-------x 2
 *
 *  - PDM_writer_PRSIM6 :
 *
 *   4 x-------x 6
 *     |\     /|
 *     | \   / |
 *   1 x- \-/ -x 3
 *      \ 5x  /
 *       \ | /
 *        \|/
 *         x 2
 *
 *  - PDM_writer_HEXA8 :
 *
 *      8 x-------x 7
 *       /|      /|
 *      / |     / |
 *   5 x-------x6 |
 *     | 4x----|--x 3
 *     | /     | /
 *     |/      |/
 *   1 x-------x 2
 *
 * \param [in] cs                  Pointer to \ref PDM_writer object
 * \param [in] id_geom             Identificateur de l'objet geometrique
 * \param [in] id_bloc             Identificateur du bloc
 * \param [in] id_part             Indice de partition
 * \param [in] t_elt               Type d'element
 * \param [in] n_elt               Nombre d'elements dans le bloc
 * \param [in] connec              Table de connectivite des elements
 * \param [in] num_part            Numerotation dans la partition
 *
 */

void
PDM_writer_geom_bloc_std_set
(
 PDM_writer_t  *cs,
 const int      id_geom,
 const int      id_bloc,
 const int      id_part,
 const int      n_elt,
 PDM_l_num_t   *connec,
 PDM_g_num_t   *numabs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_part_mesh_nodal_section_std_set(geom->_mesh_nodal, id_bloc, id_part,
                                      n_elt, connec, numabs, NULL,
                                      NULL, geom->section_owner[id_bloc]);

}


/**
 * \brief Ajout d'un bloc de polygones dans la partition courante
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] id_part         Indice de partition
 * \param [in] n_elt           Nombre d'elements dans le bloc
 * \param [in] connec_idx      Index dans la table de connectivite (dim = n_elt+1)
 * \param [in] connec          Table de connectivite des elements (dim = connec_idx[n_elt])
 * \param [in] numabs          Numerotation absolue des elements
 *
 */

void
PDM_writer_geom_bloc_poly2d_set
(
PDM_writer_t        *cs,
const int            id_geom,
const int            id_bloc,
const int            id_part,
const PDM_l_num_t    n_elt,
      PDM_l_num_t   *connec_idx,
      PDM_l_num_t   *connec,
      PDM_g_num_t   *numabs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_part_mesh_nodal_section_poly2d_set(geom->_mesh_nodal, id_bloc, id_part,
                                         n_elt, connec_idx, connec, numabs, NULL,
                                         geom->section_owner[id_bloc]);

}


/**
 * \brief Ajout d'un bloc de polyedres dans la partition courante
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] id_part         Indice de partition
 * \param [in] n_elt           Nombre d'elements dans le bloc
 * \param [in] n_face          Nombre de faces de chaque element (dim = n_elt)
 * \param [in] facsom_idx      Index dans la table de connectivite des faces (dim = n_face_total+1)
 * \param [in] facsom          Table de connectivite des faces (dim = facsom_idx[n_face_total}
 * \param [in] cellfac_idx     Index dans la table de connectivite des cellules (dim = n_elt+1)
 * \param [in] cellfac         Table de connectivite des elements (dim = cellfac_idx[n_elt])
 * \param [in] numabs          Numerotation absolue des elements
 *
 */

void
PDM_writer_geom_bloc_poly3d_set
(
PDM_writer_t        *cs,
const int            id_geom,
const int            id_bloc,
const int            id_part,
const PDM_l_num_t    n_elt,
const PDM_l_num_t    n_face,
      PDM_l_num_t   *facsom_idx,
      PDM_l_num_t   *facsom,
      PDM_l_num_t   *cellfac_idx,
      PDM_l_num_t   *cellfac,
      PDM_g_num_t   *numabs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_part_mesh_nodal_section_poly3d_set(geom->_mesh_nodal,
                                         id_bloc,
                                         id_part,
                                         n_elt,
                                         n_face,
                                         facsom_idx,
                                         facsom,
                                         NULL,
                                         cellfac_idx,
                                         cellfac,
                                         numabs,
                                         NULL,
                                         NULL,
                                         geom->section_owner[id_bloc]);
}

/**
 *
 * \brief Ajout de cellules 3D decrites en fonctions des faces.
 *
 * Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * \param [in]  cs              Pointer to \ref PDM_writer object
 * \param [in]  id_geom         Identificateur de l'objet geometrique
 * \param [in]  id_part         Identificateur de partition
 * \param [in]  n_cell          Nombre de cellules 3D ajoutées
 * \param [in]  n_face          Nombre de faces décrites
 * \param [in]  face_som_idx    Index de connectivite faces -> sommets
 * \param [in]  face_som        Connectivite faces -> sommets
 * \param [in]  cell_face_idx   Index de connectivite cellules -> faces
 * \param [in]  cell_face       Connectivite cellules -> faces
 * \param [in]  numabs          Numerotation absolue des cellules
 *
 */

void
PDM_writer_geom_cell3d_cellface_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_cell,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_l_num_t  *cell_face_idx,
 PDM_l_num_t  *cell_face_nb,
 PDM_l_num_t  *cell_face,
 PDM_g_num_t  *numabs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  /* Check geom_kind coherence */
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_VOLUMIC;

  if (geom->geom_kind == PDM_GEOMETRY_KIND_MAX) {
    geom->geom_kind = geom_kind;
  }
  else {
    if (geom_kind != geom->geom_kind) {
      PDM_error(__FILE__, __LINE__, 0,
                "Current geometry kind (%d) cannot contain element of dimension 3\n",
                (int) geom->geom_kind);
    }
  }

  int *_face_som_idx = face_som_idx;
  if (face_som_nb != NULL) {
    _face_som_idx = PDM_array_new_idx_from_sizes_int(face_som_nb, n_face);
  }

  int *_cell_face_idx = cell_face_idx;
  if (cell_face_nb != NULL) {
    _cell_face_idx = PDM_array_new_idx_from_sizes_int(cell_face_nb, n_cell);
  }

  PDM_part_mesh_nodal_cell3d_cellface_add(geom->_mesh_nodal,
                                          id_part,
                                          n_cell,
                                          n_face,
                                          _face_som_idx,
                                          face_som,
                                          NULL,
                                          _cell_face_idx,
                                          cell_face,
                                          numabs,
                                          PDM_OWNERSHIP_KEEP);

  if (face_som_nb != NULL) {
    if (geom->_face_vtx_idx[id_part] != NULL) {
      free(geom->_face_vtx_idx[id_part]);
      geom->_face_vtx_idx[id_part] = NULL;
    }
    geom->_face_vtx_idx[id_part] = _face_som_idx;
  }
  if (cell_face_nb != NULL) {
    if (geom->_cell_face_idx[id_part] != NULL) {
      free(geom->_cell_face_idx[id_part]);
      geom->_cell_face_idx[id_part] = NULL;
    }
    geom->_cell_face_idx[id_part] = _cell_face_idx;
  }

  if (0 == 1) {
    printf("n_cell : %d\n", n_cell);
    for (int i = 0; i < n_cell; i++) {
      printf ("%d "PDM_FMT_G_NUM" : \n", i+1, numabs[i]);
      for (int j = cell_face_idx[i]; j < cell_face_idx[i+1]; j++) {
        printf (" %d", cell_face[j]);
      }
      printf ("\n");
    }
    printf("n_face : %d\n", n_face);
    for (int i = 0; i < n_face; i++) {
      printf ("%d: \n", i+1);
      for (int j = face_som_idx[i]; j < face_som_idx[i+1]; j++) {
        printf (" %d", face_som[j]);
      }
      printf ("\n");
    }
  }


}


/**
 *
 * \brief Ajout de cellules 2D decrites en fonctions des faces.
 *
 * Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] n_cell          Nombre de cellules 3D ajoutées
 * \param [in] n_face          Nombre de faces décrites
 * \param [in] face_som_idx    Index de connectivite faces -> sommets
 * \param [in] face_som        Connectivite faces -> sommets
 * \param [in] cell_face_idx   Index de connectivite cellules -> faces
 * \param [in] cell_face       Connectivite cellules -> faces
 * \param [in] numabs          Numerotation absolue des cellules
 *
 */

void
PDM_writer_geom_cell2d_cellface_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_cell,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_l_num_t  *cell_face_idx,
 PDM_l_num_t  *cell_face_nb,
 PDM_l_num_t  *cell_face,
 PDM_g_num_t  *numabs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  /* Check geom_kind coherence */
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_SURFACIC;

  if (geom->geom_kind == PDM_GEOMETRY_KIND_MAX) {
    geom->geom_kind = geom_kind;
  }
  else {
    if (geom_kind != geom->geom_kind) {
      PDM_error(__FILE__, __LINE__, 0,
                "Current geometry kind (%d) cannot contain element of dimension 2\n",
                (int) geom->geom_kind);
    }
  }

  int *_cell_face_idx = cell_face_idx;
  if (cell_face_nb != NULL) {
    _cell_face_idx = PDM_array_new_idx_from_sizes_int(cell_face_nb, n_cell);
  }

  PDM_UNUSED(face_som_idx);
  PDM_UNUSED(face_som_nb );

  PDM_part_mesh_nodal_face2d_faceedge_add(geom->_mesh_nodal,
                                          id_part,
                                          n_cell,
                                          n_face,
                                          face_som,
                                          _cell_face_idx,
                                          cell_face,
                                          numabs,
                                          PDM_OWNERSHIP_KEEP);
  if (cell_face_nb != NULL) {
    geom->_cell_face_idx[id_part] = _cell_face_idx;
  }
}


/**
 *
 * \brief Ajout de faces decrites en fonctions des sommets.
 *
 * Cette fonction détermine les types des éléments et crée des blocs regrouppant les éléments
 * de même type. Elle retourne l'indirection vers le nouvel ordre de rangement
 * des cellules.
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] n_elt           Nombre de cellules 3D ajoutées
 * \param [in] n_face          Nombre de faces décrites
 * \param [in] face_som_idx    Index de connectivite faces -> sommets
 * \param [in] face_som        Connectivite faces -> sommets
 * \param [in] numabs          Numerotation absolue des faces
 *
 */

void
PDM_writer_geom_faces_facesom_add
(
 PDM_writer_t *cs,
 const int     id_geom,
 const int     id_part,
 const int     n_face,
 PDM_l_num_t  *face_som_idx,
 PDM_l_num_t  *face_som_nb,
 PDM_l_num_t  *face_som,
 PDM_g_num_t  *numabs
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  /* Check geom_kind coherence */
  PDM_geometry_kind_t geom_kind = PDM_GEOMETRY_KIND_SURFACIC;

  if (geom->geom_kind == PDM_GEOMETRY_KIND_MAX) {
    geom->geom_kind = geom_kind;
  }
  else {
    if (geom_kind != geom->geom_kind) {
      PDM_error(__FILE__, __LINE__, 0,
                "Current geometry kind (%d) cannot contain element of dimension 2\n",
                (int) geom->geom_kind);
    }
  }

  int *_face_som_idx = face_som_idx;
  if (face_som_nb != NULL) {
    _face_som_idx = PDM_array_new_idx_from_sizes_int(face_som_nb, n_face);
  }

  PDM_part_mesh_nodal_faces_facevtx_add(geom->_mesh_nodal,
                                        id_part,
                                        n_face,
                                        _face_som_idx,
                                        face_som,
                                        numabs,
                                        PDM_OWNERSHIP_KEEP);
  if (face_som_nb != NULL) {
    free(_face_som_idx);
  }
}

/**
 * \brief Ecriture du maillage courant
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */

void
PDM_writer_geom_write
(
 PDM_writer_t *cs,
 const int     id_geom
 )
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  //TODO  faire un retour si geometrie n'est pas dependante du temps
  //       et si on n'est pas au premier increment
  /* Mise a jour du nombre total d'elements */


  /* Determination de la numerotation absolue interne des elements
     Independante du parallelisme */
  const int n_blocks   = PDM_part_mesh_nodal_n_section_get  (geom->mesh_nodal);
  // const int *blocks_id = PDM_part_mesh_nodal_sections_id_get(geom->mesh_nodal);

  for (int i = 0; i < n_blocks; i++) {
    PDM_part_mesh_nodal_g_num_in_section_compute(geom->mesh_nodal, i,
                                                 PDM_OWNERSHIP_KEEP);
  }

  /* Ecriture au format */

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->geom_write_fct != NULL) {
    (fmt_ptr->geom_write_fct) (geom);
  }

}


/**
 * \brief Liberation des donnees decrivant le maillage courant
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */

void
PDM_writer_geom_free
(
 PDM_writer_t *cs,
 const int     id_geom
)
{
  PDM_writer_geom_data_free(cs,
                            id_geom);

  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (cs->geom_tab == NULL) {
    return;
  }

  if (cs->geom_tab->geom == NULL) {
    return;
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom != NULL) {

    if (geom->_mesh_nodal != NULL) {
      PDM_part_mesh_nodal_free(geom->_mesh_nodal);
      geom->_mesh_nodal = NULL;
      geom->mesh_nodal = NULL;      
    }
    if (geom->nom_geom != NULL) {
      free(geom->nom_geom);
      geom->nom_geom = NULL;
    }

    /* Liberation specifique au format */

    /* Appel de la fonction complementaire propre au format */

    PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

    if (fmt_ptr->geom_free_fct != NULL) {
      (fmt_ptr->geom_free_fct) (geom);
    }

    if (geom->section_owner != NULL) {
      free(geom->section_owner);
    }

    if (geom->_face_vtx_idx != NULL) {
      for (int i = 0; i < geom->n_part; i++) {
        if (geom->_face_vtx_idx[i] != NULL) {
          free(geom->_face_vtx_idx[i]);
        }
      }
      free(geom->_face_vtx_idx);
    }

    if (geom->_cell_face_idx != NULL) {
      for (int i = 0; i < geom->n_part; i++) {
        if (geom->_cell_face_idx[i] != NULL) {
          free(geom->_cell_face_idx[i]);
        }
      }
      free(geom->_cell_face_idx);
    }

    free(geom);
    cs->geom_tab->geom[id_geom] = NULL;
  }
}


/**
 * \brief Liberation des donnees decrivant le maillage courant
 * les indirections sur les numérotation absolues sont conservées
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */

void
PDM_writer_geom_data_free
(
 PDM_writer_t *cs,
 const int     id_geom
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (cs->geom_tab == NULL) {
    return;
  }

  if (cs->geom_tab->geom == NULL) {
    return;
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom != NULL) {

    if (geom->_mesh_nodal != NULL) {
      PDM_part_mesh_nodal_partial_free(geom->_mesh_nodal);
    }

  }
}

/**
 * \brief Mapping des noms de variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] public_name     Nom Public de la variable
 * \param [in] pivate_name     Nom privé de la variable
 *
 */

void
PDM_writer_name_map_add
(
 PDM_writer_t *cs,
 const char   *public_name,
 const char   *private_name
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->name_map_tab == NULL) {
    cs->name_map_tab = _pdm_writer_name_map_tab_create(3);
  }


  PDM_writer_name_map_t *name_map = (PDM_writer_name_map_t *) malloc (sizeof(PDM_writer_name_map_t));

  _pdm_writer_name_map_tab_add(cs->name_map_tab, name_map);


  name_map->public_name = malloc ((strlen(public_name) + 1) * sizeof(char));
  name_map->private_name = malloc ((strlen(private_name) + 1) * sizeof(char));

  strcpy(name_map->public_name, public_name);
  strcpy(name_map->private_name, private_name);

}


/**
 * \brief Create a global constant variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] nom_var         Nom de la variable
 * \param [in] val_var         Valeur de la variable
 *
 * \return  Identificateur de l'objet variable
 *
 */

int
PDM_writer_cst_global_var_create
(
 PDM_writer_t               *cs,
 const char                 *nom_var,
 const double                val_var
)
{
 
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->cst_global_var_tab.s_var == 0) {

    cs->cst_global_var_tab.n_var = 0;
    cs->cst_global_var_tab.s_var = 4;

    cs->cst_global_var_tab.var = (PDM_writer_cst_global_var_t **) malloc(sizeof(PDM_writer_cst_global_var_t *) * cs->cst_global_var_tab.s_var);
    for (int i = 0; i < cs->cst_global_var_tab.s_var; i++) {
      cs->cst_global_var_tab.var[i] = NULL;
    }
  }

  /* Allocation de la structure PDM_writer_var_t */

  PDM_writer_cst_global_var_t *var = (PDM_writer_cst_global_var_t *) malloc(sizeof(PDM_writer_cst_global_var_t));

  var->nom_var = malloc(sizeof(char) * (1 + strlen (nom_var)));
  strcpy (var->nom_var, nom_var);

  var->_val = val_var;  

  if (cs->cst_global_var_tab.n_var >= cs->cst_global_var_tab.s_var) {
    cs->cst_global_var_tab.s_var = PDM_MAX(2*cs->cst_global_var_tab.s_var, cs->cst_global_var_tab.n_var+1);
    cs->cst_global_var_tab.var = 
    (PDM_writer_cst_global_var_t **) realloc(cs->cst_global_var_tab.var, sizeof(PDM_writer_cst_global_var_t *) * cs->cst_global_var_tab.s_var);

    for (int i = cs->cst_global_var_tab.n_var+1; i < cs->cst_global_var_tab.s_var; i++) {
      cs->cst_global_var_tab.var[i] = NULL;
    }
  }

  int id_var = cs->cst_global_var_tab.n_var;
  cs->cst_global_var_tab.n_var++;

  cs->cst_global_var_tab.var[id_var] = var;

  return id_var;

}




/**
 * \brief Set a global constant variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_var          Variable id
 * \param [in] val_var         Valeur de la variable
 *
 * \return  Identificateur de l'objet variable
 *
 */

void
PDM_writer_cst_global_var_set
(
 PDM_writer_t               *cs,
 const int                   id_var,
 const double                val_var
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  assert (id_var < cs->cst_global_var_tab.s_var);
  assert (cs->cst_global_var_tab.var[id_var] != NULL);

  cs->cst_global_var_tab.var[id_var]->_val = val_var;

}



/**
 * \brief Creation d'une variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] st_dep_temps    Indique si la variable est dependante du temps
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] dim             Dimension de la variable
 * \param [in] loc             Localisation de la variable
 * \param [in] nom_var         Nom de la variable
 *
 * \return  Identificateur de l'objet variable
 *
 */

int
PDM_writer_var_create
(
 PDM_writer_t               *cs,
 const PDM_writer_status_t   st_dep_tps,
 const PDM_writer_var_dim_t  dim,
 const PDM_writer_var_loc_t  loc,
 const char                 *nom_var
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  /* Mise a jour du tableau de stockage */

  if (cs->var_tab == NULL) {
    cs->var_tab = _pdm_writer_var_tab_create(4);
  }

  /* Allocation de la structure PDM_writer_var_t */

  PDM_writer_var_t *var = (PDM_writer_var_t *) malloc(sizeof(PDM_writer_var_t));
  int id_var = _pdm_writer_var_tab_add(cs->var_tab, var);

  /* Initialisation de la structure PDM_writer_var_t */

  _var_init (var);

  size_t l_nom_var = strlen(nom_var);
  var->nom_var = (char *) malloc(sizeof(char) * (l_nom_var + 1));
  strcpy(var->nom_var, nom_var);   /* Nom de la variable */

  var->st_dep_tps = st_dep_tps;    /* Variable en temps */
  var->dim        = dim;           /* Dimension de la variable */
  var->loc        = loc;           /* Dimension de la variable */
  var->_cs        = cs;
  var->private_name = NULL;

  if (cs->name_map_tab != NULL) {
    const int n_map = cs->name_map_tab->n_name_map;

    for (int i = 0; i < n_map; i++) {
      PDM_writer_name_map_t *map = cs->name_map_tab->name_map[i];
      if (!strcmp(nom_var, map->public_name)) {
        var->private_name = map->private_name;
      }
    }
  }

  /* Appel de la fonction complementaire propre au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->var_create_fct != NULL) {
    (fmt_ptr->var_create_fct) (var);
  }

  return id_var;
}


/**
 * \brief Ecriture des valeurs de la variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_var          Identificateur de la variable a ecrire
 *
 */

void
PDM_writer_var_write
(
 PDM_writer_t *cs,
 const int     id_var
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_var >= cs->var_tab->n_var) {
    PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
    abort();
  }

  PDM_writer_var_t *var = cs->var_tab->var[id_var];

  if (var == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad var identifier\n");
  }

  /* Ecriture au format */

  PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

  if (fmt_ptr->var_write_fct != NULL) {
    (fmt_ptr->var_write_fct) (var);
  }

}



/**
 * \brief Mise a jour des valeurs de la variable.
 *
 * Attention, les valeurs définies aux elements doivent être définies suivant l'ordre de définition des blocs !
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 * \param [in] id_part         Identificateur de la partition dans l'objet geometrique
 * \param [in] id_var          Identificateur de la variable mise à jour
 * \param [in] val             Valeurs
 *
 */

void
PDM_writer_var_set
(
 PDM_writer_t     *cs,
 const int         id_var,
 const int         id_geom,
 const int         id_part,
 const PDM_real_t *val
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_var >= cs->var_tab->n_var) {
    PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
    abort();
  }

  PDM_writer_var_t *var = cs->var_tab->var[id_var];

  if (var == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad var identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }

  const int n_ind = cs->geom_tab->n_geom;

  if (var->_val == NULL) {
    var->_val = (double ***) malloc(sizeof(double **) * n_ind);
    for (int i = 0; i < n_ind; i++) {
      var->_val[i] = NULL;
    }
  }

  if (n_ind <= id_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_var_set    : Indice de geometrie incorrect\n");
    abort();
  }

  int n_part = geom->n_part;

  if (var->_val[id_geom] == NULL) {
    var->_val[id_geom] = (double **) malloc(sizeof(double *) * n_part);
    for (int i = 0; i < n_part; i++)
      var->_val[id_geom][i] = NULL;
  }

  double **val_geom = var->_val[id_geom];

  if (n_part <= id_part) {
    PDM_error(__FILE__, __LINE__, 0, "Erreur cs_var_set    : Indice de partition incorrect\n");
    abort();
  }

  if (geom->geom_kind == PDM_GEOMETRY_KIND_MAX) {
    PDM_error(__FILE__, __LINE__, 0, "Undefined geometry kind\n");
  }
  int n_cell = PDM_part_mesh_nodal_n_elmts_get(geom->mesh_nodal,
                                               geom->geom_kind,
                                               id_part);
  int *num_cell_parent_to_local = PDM_part_mesh_nodal_num_elmt_parent_to_local_get(geom->mesh_nodal,
                                                                                   geom->geom_kind,
                                                                                   id_part);
  int n_som = PDM_part_mesh_nodal_n_vtx_get(geom->mesh_nodal, id_part);

  if (var->loc == PDM_WRITER_VAR_ELEMENTS) {
    val_geom[id_part] = (double *) malloc(sizeof(double) * var->dim * n_cell);
    if (num_cell_parent_to_local != NULL) {
      for (int i = 0; i < n_cell; i++) {
        for (int j = 0; j < (int) var->dim; j++)
          val_geom[id_part][var->dim * num_cell_parent_to_local[i]+j] = val[i*var->dim + j];
      }
    }
    else {
      for (int i = 0; i < n_cell; i++) {
        for (int j = 0; j < (int) var->dim; j++)
          val_geom[id_part][var->dim * i+j] = val[i*var->dim + j];
      }
    }
  }
  else {
    val_geom[id_part] = (double *) malloc(sizeof(double) * var->dim * n_som);
    for (int i = 0; i < n_som * (int) var->dim; i++) {
      val_geom[id_part][i] = val[i];
    }
  }

}


/**
 * \brief Liberation du tableau de donnees des variables
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_var          Identificateur de la variable
 *
 */

void
PDM_writer_var_data_free
(
 PDM_writer_t *cs,
 const int     id_var
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_var >= cs->var_tab->n_var) {
    PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
    abort();
  }

  PDM_writer_var_t *var = cs->var_tab->var[id_var];

  if (var != NULL) {

    if (var->_val != NULL) {

      const int n_ind = cs->geom_tab->n_geom;

      for (int i = 0; i < n_ind; i++) {
        int idx = i;
        PDM_writer_geom_t *geom = cs->geom_tab->geom[idx];

        if (geom == NULL) {
          PDM_error(__FILE__, __LINE__, 0, 
            "PDM_writer_var_data_free - Bad geom identifier : An associated geom of var '%s' is free before the var\n", var->nom_var);
          abort();
        }

        int n_part = geom->n_part;

        if ((geom != NULL) && (var->_val[idx] != NULL)) {
          for (int j = 0; j < n_part; j++) {
            if (var->_val[idx][j] != NULL)
              free(var->_val[idx][j]);
            var->_val[idx][j] = NULL;
          }
        }
      }
    }
  }
}


/**
 * \brief Liberation d'une variable
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_var          Identificateur de la variable
 *
 */

void
PDM_writer_var_free
(
 PDM_writer_t *cs,
 const int     id_var
 )
{

  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (cs->var_tab != NULL) {

    PDM_writer_var_data_free(cs, id_var);

    /* Acces a l'objet de geometrie courant */

    if (id_var >= cs->var_tab->n_var) {
      PDM_error(__FILE__, __LINE__, 0, "Bad var identifier\n");
      abort();
    }

    PDM_writer_var_t *var = cs->var_tab->var[id_var];

    if (var != NULL) {

      free(var->nom_var);

      const int n_ind = cs->geom_tab->n_geom;

      if (var->_val != NULL) {
        for (int i = 0; i < n_ind; i++) {
          int idx = i;
          if (var->_val[idx] != NULL)
            free(var->_val[idx]);
          var->_val[idx] = NULL;
        }

        free(var->_val);
        var->_val = NULL;
      }

      /* Liberation specifique au format */

      PDM_writer_fmt_t *fmt_ptr = fmt_tab[cs->fmt_id];

      if (fmt_ptr->var_free_fct != NULL) {
        (fmt_ptr->var_free_fct) (var);
      }

      free (var);
      cs->var_tab->var[id_var] = NULL;

    }
  }
}

/**
 * \brief Add a writer format
 *
 * Define a new format writer
 *
 * \param [in] name            Name
 * \param [in] create_fct      Customize \ref PDM_writer_create function for the new format  (or NULL)
 * \param [in] free_fct        Customize \ref PDM_writer_free function for the new format (or NULL)
 * \param [in] beg_step_fct    Customize \ref PDM_writer_step_beg function for the new format (or NULL)
 * \param [in] end_step_fct    Customize \ref PDM_writer_step_end function for the new format (or NULL)
 * \param [in] geom_create_fct Customize \ref PDM_writer_geom_create function for the new format (or NULL)
 * \param [in] geom_write_fct  Customize \ref PDM_writer_geom_write function for the new format
 * \param [in] geom_free_fct   Customize \ref PDM_writer_geom_free function for the new format (or NULL)
 * \param [in] var_create_fct  Customize \ref PDM_writer_var_create function for the new format (or NULL)
 * \param [in] var_write_fct   Customize \ref PDM_writer_var_write function for the new format
 * \param [in] var_free_fct    Customize \ref PDM_writer_var_free function for the new format (or NULL)
 *
 */

void
PDM_writer_fmt_add
(
 const char                  *name,            /*!< Name                                                     */
 const PDM_writer_fct_t       create_fct,      /*!< Customize \ref PDM_writer_create function for the format */
 const PDM_writer_fct_t       free_fct,        /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t       beg_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_fct_t       end_step_fct,    /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t  geom_create_fct, /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t  geom_write_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_geom_fct_t  geom_free_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t   var_create_fct,  /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t   var_write_fct,   /*!< Customize \ref PDM_writer_free function for the format   */
 const PDM_writer_var_fct_t   var_free_fct     /*!< Customize \ref PDM_writer_free function for the format   */
)
{
  _load_intern_fmt();

  if (geom_write_fct == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_fmt_add : Undefined geom write function\n");
    abort ();
  }

  if (var_write_fct == NULL) {
    PDM_error(__FILE__, __LINE__, 0, "Error PDM_writer_fmt_add : Undefined var write function\n");
    abort ();
  }

  if (n_fmt_tab >= s_fmt_tab) {
    s_fmt_tab = PDM_MAX (2*s_fmt_tab, n_fmt_tab+1);
    fmt_tab = realloc(fmt_tab, sizeof(PDM_writer_fmt_t *) * s_fmt_tab);
  }

  PDM_writer_fmt_t *fmt_ptr = malloc (sizeof(PDM_writer_fmt_t));
  fmt_tab[n_fmt_tab++] = fmt_ptr;

  fmt_ptr->name            = malloc(sizeof(char) * (strlen(name) + 1));
  strcpy (fmt_ptr->name, name);

  fmt_ptr->create_fct      = create_fct;
  fmt_ptr->free_fct        = free_fct;
  fmt_ptr->beg_step_fct    = beg_step_fct;
  fmt_ptr->end_step_fct    = end_step_fct;
  fmt_ptr->geom_create_fct = geom_create_fct;
  fmt_ptr->geom_write_fct  = geom_write_fct;
  fmt_ptr->geom_free_fct   = geom_free_fct;
  fmt_ptr->var_create_fct  = var_create_fct;
  fmt_ptr->var_write_fct   = var_write_fct;
  fmt_ptr->var_free_fct    = var_free_fct;

}


/**
 * \brief Free formats
 *
 */

void
PDM_writer_fmt_free
(
 void
)
{
  if (fmt_tab != NULL) {

    for (int i = 0; i < n_fmt_tab; i++) {
      if (fmt_tab[i] != NULL) {
        free(fmt_tab[i]->name);
        free(fmt_tab[i]);
        fmt_tab[i] = NULL;
      }
    }

    n_fmt_tab = 0;

    free(fmt_tab);
    fmt_tab = NULL;

  }
}

/**
 * \brief Réinitialisation des donnees decrivant le maillage courant
 *
 * \param [in] cs              Pointer to \ref PDM_writer object
 * \param [in] id_geom         Identificateur de l'objet geometrique
 *
 */

void
PDM_writer_geom_data_reset
(
 PDM_writer_t *cs,
 const int     id_geom
)
{
  if (cs == NULL) {
    PDM_error (__FILE__, __LINE__, 0, "Bad writer identifier\n");
  }

  if (id_geom >= cs->geom_tab->n_geom) {
    PDM_error(__FILE__, __LINE__, 0, "Bad geom identifier\n");
    abort();
  }
  PDM_writer_geom_t *geom = cs->geom_tab->geom[id_geom];

  if (geom != NULL) {

    if (geom->_mesh_nodal != NULL) {
      PDM_part_mesh_nodal_reset(geom->_mesh_nodal);
    }

  }
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
