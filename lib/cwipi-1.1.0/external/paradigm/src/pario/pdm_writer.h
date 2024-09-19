/*
 * \file
 */

#ifndef __PDM_WRITER_H__
#define __PDM_WRITER_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm.h"
#include "pdm_mpi.h"
#include "pdm_io.h"
#include "pdm_mesh_nodal.h"
#include "pdm_part_mesh_nodal.h"

/*=============================================================================
 * Definitions des macro
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Definition des types
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Statut
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_OFF,
  PDM_WRITER_ON

} PDM_writer_status_t;

/*----------------------------------------------------------------------------
 * Type de topologie
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_TOPO_CST ,
  PDM_WRITER_TOPO_DEFORMABLE,
  PDM_WRITER_TOPO_VARIABLE

} PDM_writer_topology_t;

/*----------------------------------------------------------------------------
 * Type d'elements géometriques (It's the same than the type defined into PDM_Mesh_nodal)
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_POINT = PDM_MESH_NODAL_POINT,
  PDM_WRITER_BAR2 = PDM_MESH_NODAL_BAR2,
  PDM_WRITER_TRIA3 = PDM_MESH_NODAL_TRIA3,
  PDM_WRITER_QUAD4 = PDM_MESH_NODAL_QUAD4,
  PDM_WRITER_POLY_2D = PDM_MESH_NODAL_POLY_2D,
  PDM_WRITER_TETRA4 = PDM_MESH_NODAL_TETRA4,
  PDM_WRITER_PYRAMID5 = PDM_MESH_NODAL_PYRAMID5,
  PDM_WRITER_PRISM6 = PDM_MESH_NODAL_PRISM6,
  PDM_WRITER_HEXA8 = PDM_MESH_NODAL_HEXA8,
  PDM_WRITER_POLY_3D  = PDM_MESH_NODAL_POLY_3D

} PDM_writer_elt_geom_t;


typedef struct _PDM_writer_t PDM_writer_t;

typedef struct _PDM_writer_geom_t PDM_writer_geom_t;

typedef struct _PDM_writer_var_t PDM_writer_var_t;

typedef struct _PDM_writer_cst_global_var_t PDM_writer_cst_global_var_t;


/**
 * \struct PDM_writer_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_t instance
 *
 */

typedef void (*PDM_writer_fct_t) (PDM_writer_t *pw);


/**
 * \struct PDM_writer_geom_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_geom_t instance
 *
 */

typedef void (*PDM_writer_geom_fct_t) (PDM_writer_geom_t *geom);


/**
 * \struct PDM_writer_var_fct_t
 *
 * \brief  Function pointer used to define an user format with \ref PDM_writer_var_t instance
 *
 */

typedef void (*PDM_writer_var_fct_t) (PDM_writer_var_t *var);


/*----------------------------------------------------------------------------
 * Format du fichie de sortie
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_FMT_BIN,
  PDM_WRITER_FMT_ASCII

} PDM_writer_fmt_fic_t;

/*----------------------------------------------------------------------------
 * Types des variables
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_VAR_CST           = 0,
  PDM_WRITER_VAR_SCALAR        = 1,
  PDM_WRITER_VAR_VECTOR        = 3,
  PDM_WRITER_VAR_TENSOR_SYM    = 6,
  PDM_WRITER_VAR_TENSOR_ASYM   = 9

} PDM_writer_var_dim_t;

/*----------------------------------------------------------------------------
 * Localisation des variables
 *----------------------------------------------------------------------------*/

typedef enum {

  PDM_WRITER_VAR_VERTICES,
  PDM_WRITER_VAR_ELEMENTS,
  PDM_WRITER_VAR_PARTICLES

} PDM_writer_var_loc_t;


/*=============================================================================
 * Variables globales
 *============================================================================*/

/*=============================================================================
 * Prototypes des fonctions publiques
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
);


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
);


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
);


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
);


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
);


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
);


int
PDM_writer_geom_create_from_mesh_nodal
(
 PDM_writer_t              *cs,
 const char                *nom_geom,
 PDM_part_mesh_nodal_t     *mesh
);


void
PDM_writer_geom_set_from_mesh_nodal
(
 PDM_writer_t              *cs,
 const int                  id_geom,
 PDM_part_mesh_nodal_t     *mesh
);


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
 const PDM_ownership_t ownership
);


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
);


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
);


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
);


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
);


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
);


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
 );


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
);


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
);


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
 );


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
 );


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
);


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
);


/**
 * \brief Create a global constant variable
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
);



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
);


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
);


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
);


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
);


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
);


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
);


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
);


/**
 * \brief Free format
 *
 */

void
PDM_writer_fmt_free
(
 void
);


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
 );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_WRITER_H__ */
