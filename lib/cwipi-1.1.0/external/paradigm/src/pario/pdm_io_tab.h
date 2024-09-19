/*
 * \file
 */

#ifndef __PDM_IO_TAB_H__
#define __PDM_IO_TAB_H__

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_io.h"

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
 * Types decrivant un tableau
 *----------------------------------------------------------------------------*/

typedef struct _PDM_io_array_t PDM_io_array_t;

/*=============================================================================
 * Variables globales
 *============================================================================*/


/*=============================================================================
 * Prototypes des fonctions publiques
 *============================================================================*/


/**
 * \brief Initialise une phase d'écriture parallèle de tableaux de données associées
 * aux numéros de variables PDM
 * Chaque tableau a ses propres caractéristiques :
 *         - taille de données
 *         - nombre de donnée
 *         - indirection (numérotation absolue)
 *
 * \param [in] unite              Unite du fichier
 * \param [in] t_rangement        Type de rangement
 * \param [in] num_var_cedre_max  Numéro max de variable PDM
 * \param [in] n_partition_local  Nombre de partitions locales
 *
 */

void PDM_io_array_write_beg
(
 PDM_io_file_t         *unite,
 const PDM_stride_t  t_rangement,
 const PDM_l_num_t         num_var_cedre_max,
 const PDM_l_num_t         n_partition_local
);


/**
 * \brief Ajoute une partie des donnees dans un tableau associés à une variable
 * PDM
 *
 * \param [in] num_var_cedre          Numéro de variable PDM
 * \param [in] i_part                 indice de partition
 * \param [in] n_composantes          Nombre de composantes pour chaque donnee
 * \param [in] n_donnees              Nombre de donnees a lire
 * \param [in] indirection            Indirection de redistribition des donnees
 * \param [in] donnees                Donnees a écrire
 *
 */


void PDM_io_array_write_data_append
(
 const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            i_part,
 const PDM_l_num_t           *n_composantes,
 const PDM_l_num_t            n_donnees,
 const PDM_g_num_t           *indirection,
 void                        *donnees
 );

/**
 * \brief Definition d'une variable en ecriture
 *
 * \param [in] num_var_cedre          Numéro de variable PDM
 * \param [in] num_indirection_cedre  Numéro d'indirection PDM
 * \param [in] t_n_composantes        Type de tailles composantes (PDM_STRIDE_CST_INTERLACED ou PDM_STRIDE_VAR_INTERLACED)
 * \param [in] n_composantes          Nombre de composantes pour chaque donnee
 * \param [in] taille_donnee          Taille unitaire de la donnnee
 *
 */

void PDM_io_array_write_var_def
(
 const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            num_indirection_cedre,
 const PDM_stride_t t_n_composantes,
 const PDM_l_num_t            n_composantes,
 const PDM_l_num_t            taille_donnee
 );

/**
 * \brief Finalise une phase d'écriture parallèle de tableaux de données associées
 * aux numéros de variables PDM. Cette fonction déclenche réellement
 * les écritures
 *
 */

void PDM_io_array_write_end
(void);


/**
 * \brief Initialise une phase de lecture parallèle de tableaux de données associées
 * aux numéros de variables PDM
 * Chaque tableau a ses propres caractéristiques :
 *         - taille de données
 *         - nombre de donnée
 *         - indirection (numérotation absolue)
 *
 * \param [in] unite              Unite du fichier
 * \param [in] t_rangement        Type de rangement
 * \param [in] num_var_cedre_max  Numéro max de variable PDM
 * \param [in] n_partition_local  Nombre de partitions locales
 *
 */

void PDM_io_array_read_beg
(
 PDM_io_file_t         *unite,
 const PDM_stride_t  t_rangement,
 const PDM_l_num_t         num_var_cedre_max,
 const PDM_l_num_t         n_partition_local
);


/**
 * \brief Ajoute une partie des donnees dans un tableau associés à une variable PDM
 *
 * \param [in] num_var_cedre          Numéro de variable PDM
 * \param [in] i_part                 indice de partition
 * \param [in] n_composantes          Nombre de composantes pour chaque donnee
 * \param [in] n_donnees              Nombre de donnees a lire
 * \param [in] indirection            Indirection de redistribition des donnees
 * \param [in] donnees                Donnees a écrire
 *
 */

void PDM_io_array_read_data_append
(
 const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            i_part,
 const PDM_l_num_t           *n_composantes,
 const PDM_l_num_t            n_donnees,
 const PDM_g_num_t           *indirection,
 void                        *donnees
 );

/**
 * \brief Definition d'une variable en ecriture
 *
 * \param [in] num_var_cedre          Numéro de variable PDM
 * \param [in] num_indirection_cedre  Numéro d'indirection PDM
 * \param [in] t_n_composantes        Type de tailles composantes (PDM_STRIDE_CST_INTERLACED ou PDM_STRIDE_VAR_INTERLACED)
 * \param [in] n_composantes          Nombre de composantes pour chaque donnee
 * \param [in] taille_donnee          Taille unitaire de la donnnee
 *
 */

void PDM_io_array_read_var_def
(
 const PDM_l_num_t            num_var_cedre,
 const PDM_l_num_t            num_indirection_cedre,
 const PDM_stride_t t_n_composantes,
 const PDM_l_num_t            n_composantes,
 const PDM_l_num_t            taille_donnee
 );

/**
 * \brief Finalise une phase de lecture parallèle de tableaux de données associées
 * aux numéros de variables PDM. Cette fonction déclenche réellement
 * les écritures
 *
 */

void PDM_io_array_read_end
(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_IO_TAB_H__ */
