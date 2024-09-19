#ifndef __CWIPI_TIMER_H__
#define __CWIPI_TIMER_H__

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

/*============================================================================
 * Description d'un fichier sequentiel
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
 * Structure de mesure des temps d'execution
 *----------------------------------------------------------------------------*/

typedef struct _cwipi_timer_t CWIPI_timer_t;

/*============================================================================
 * Interfaces des fonctions publiques
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Creation d'un objet timer
 *
 * return
 *   timer
 *
 *----------------------------------------------------------------------------*/

CWIPI_timer_t *CWIPI_timer_cree(void);

/*----------------------------------------------------------------------------
 * Reinitialisation des compteurs de temps
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void CWIPI_timer_init(CWIPI_timer_t *timer);

/*----------------------------------------------------------------------------
 * Reprend la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void CWIPI_timer_reprise(CWIPI_timer_t *timer);

/*----------------------------------------------------------------------------
 * Suspend la mesure du temps ecoule
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

void CWIPI_timer_suspend(CWIPI_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps CPU en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double CWIPI_timer_tps_cpu(CWIPI_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps CPU utilisateur en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double CWIPI_timer_tps_cpu_user(CWIPI_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps CPU systeme en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double CWIPI_timer_tps_cpu_sys(CWIPI_timer_t *timer);

/*----------------------------------------------------------------------------
 * Retourne le temps elaps en secondes
 *
 * parameters :
 *   timer            <-- Timer
 * return
 *----------------------------------------------------------------------------*/

double CWIPI_timer_tps_elapsed(CWIPI_timer_t *timer);

/*----------------------------------------------------------------------------
 * Destruction d'un objet timer
 *
 * parameters :
 *   timer            <-- Timer
 *
 *----------------------------------------------------------------------------*/

void CWIPI_timer_detruit(CWIPI_timer_t *timer);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __FICHIER_SEQ_H__ */
