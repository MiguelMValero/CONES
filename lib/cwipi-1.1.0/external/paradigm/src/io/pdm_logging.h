/*
 * \file
 */

#ifndef __PDM_LOGGING_H__
#define __PDM_LOGGING_H__

/*-----------------------------------------------------------------------------*/

/* Standard C library headers */
#include <stdio.h>
#include <stdarg.h>

/* BFT library headers */

#include "pdm.h"
/*-----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/*============================================================================
 * Public types
 *============================================================================*/

/*============================================================================
 * Public function prototypes
 *============================================================================*/

typedef void (*log_lock_fn)(void *udata, int lock);

enum { LOG_TRACE, LOG_DEBUG, LOG_INFO, LOG_WARN, LOG_ERROR, LOG_FATAL };

#define log_trace(...) log_log(LOG_TRACE, __FILE__, __LINE__, __VA_ARGS__)
#define log_debug(...) log_log(LOG_DEBUG, __FILE__, __LINE__, __VA_ARGS__)
#define log_info(...)  log_log(LOG_INFO,  __FILE__, __LINE__, __VA_ARGS__)
#define log_warn(...)  log_log(LOG_WARN,  __FILE__, __LINE__, __VA_ARGS__)
#define log_error(...) log_log(LOG_ERROR, __FILE__, __LINE__, __VA_ARGS__)
#define log_fatal(...) log_log(LOG_FATAL, __FILE__, __LINE__, __VA_ARGS__)

void log_set_udata(void *udata);
void log_set_lock(log_lock_fn fn);
void log_set_fp(FILE *fp);
void log_set_level(int level);
void log_set_quiet(int enable);

void log_log(int level, const char *file, int line, const char *fmt, ...);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    larray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_int
(
 const int*  array,
 const int   larray,
 const char* header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    larray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_long
(
 const PDM_g_num_t* array,
 const int          larray,
 const char*        header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    larray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_double
(
 const double*      array,
 const int          larray,
 const char*        header
);


/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array        Array to print
 * \param [in]    larray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_array_size_t
(
 const size_t *array,
 const int     larray,
 const char   *header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array_idx    Index of array to print
 * \param [inout] array        Array to print
 * \param [in]    larray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_connectivity_long
(
 const int         *array_idx,
 const PDM_g_num_t *array,
 const int          larray,
 const char*        header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [in] entitiy1_entity2_idx Index of connectivity of entity1->entity2
 * \param [in] entitiy1_entity2     Connectivity of entity1->entity2
 * \param [in] entitiy1_ln_to_gn    Global number of 1-entities
 * \param [in] entitiy2_ln_to_gn    Global number of 2-entities
 * \param [in] n_entity1            Number of 1-entities
 * \param [in] header               Print header
 *
 */
void
PDM_log_trace_part_connectivity_gnum
(
 const int         *entitiy1_entity2_idx,
 const int         *entitiy1_entity2,
 const PDM_g_num_t *entitiy1_ln_to_gn,
 const PDM_g_num_t *entitiy2_ln_to_gn,
 const int          n_entity1,
 const char*        header
);


/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array_idx    Index of array to print
 * \param [inout] array        Array to print
 * \param [in]    larray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_connectivity_int
(
 const int         *array_idx,
 const int         *array,
 const int          larray,
 const char*        header
);

/**
 *
 * \brief Pretty print of array in trace_log
 *
 * \param [inout] array_idx    Index of array to print
 * \param [inout] array        Array to print
 * \param [in]    g_num        Associated global number
 * \param [in]    larray       Array length
 * \param [inout] header       First line of log
 *
 */
void
PDM_log_trace_connectivity_int2
(
 const int         *array_idx,
 const int         *array,
 const PDM_g_num_t *g_num,
 const int          larray,
 const char*        header
);
/*----------------------------------------------------------------------------*/


void
PDM_log_trace_graph_nuplet_int
(
 const int         *array_idx,
 const int         *array_desc,
 const int          strid,
 const int          larray,
 const char*        header
);



#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __PDM_PRINTF_H__ */
