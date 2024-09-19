/*============================================================================
 * Mesure des temps CPU et elapsed
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include "pdm_config.h"


/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "pdm_memory_stats.h"
#include "pdm_memory_stats_priv.h"
#include "pdm_printf.h"
#include "pdm_error.h"
#include "pdm_logging.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force back Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */



static
void
readable_fs
(
  long   size,
  char  *buf
)
{
  double dbl_size = (double) size;
  int i = 0;
  const char* units[] = {"B", "kB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"};
  while (dbl_size > 1024) {
    dbl_size /= 1024;
    i++;
  }
  sprintf(buf, "%12.5e %s", dbl_size, units[i]);
}


/*
 * Measures the current (and peak) resident and virtual memories
 * usage of your linux C process, in kB --> https://stackoverflow.com/questions/1558402/memory-usage-of-current-process-in-c
 */
static
void
PDM_get_current_memory
(
  long* curr_real_mem,
  long* peak_real_mem,
  long* curr_virt_mem,
  long* peak_virt_mem
)
{
  // stores each word in status file
  char buffer[1024] = "";

  // linux file contains this-process info
  FILE* file = fopen("/proc/self/status", "r");

  // read the entire file
  while (fscanf(file, " %1023s", buffer) == 1) {

    if (strcmp(buffer, "VmRSS:") == 0) {
        fscanf(file, " %ld", curr_real_mem);
    }
    if (strcmp(buffer, "VmHWM:") == 0) {
        fscanf(file, " %ld", peak_real_mem);
    }
    if (strcmp(buffer, "VmSize:") == 0) {
        fscanf(file, " %ld", curr_virt_mem);
    }
    if (strcmp(buffer, "VmPeak:") == 0) {
        fscanf(file, " %ld", peak_virt_mem);
    }
  }
  fclose(file);
}

PDM_memory_stats_t*
PDM_memory_stats_create
(
 int          n_memory_snapshot,
 PDM_MPI_Comm comm
)
{
  PDM_memory_stats_t* ms = malloc(sizeof(PDM_memory_stats_t));

  ms->comm              = comm;
  ms->n_memory_snapshot = n_memory_snapshot;

  ms->snapshot_name = malloc( n_memory_snapshot * sizeof(char *));
  ms->curr_real_mem = malloc( n_memory_snapshot * sizeof(long  ));
  ms->peak_real_mem = malloc( n_memory_snapshot * sizeof(long  ));
  ms->curr_virt_mem = malloc( n_memory_snapshot * sizeof(long  ));
  ms->peak_virt_mem = malloc( n_memory_snapshot * sizeof(long  ));

  for(int i = 0; i < ms->n_memory_snapshot; ++i) {
    ms->snapshot_name[i] = NULL;
  }

  return ms;
}

void
PDM_memory_stats_add
(
 PDM_memory_stats_t *ms,
 int                 i_snapshot,
 const char         *name
)
{
  assert(ms->snapshot_name[i_snapshot] == NULL);
  ms->snapshot_name[i_snapshot] = malloc( (strlen(name)+1) * sizeof(char));
  strcpy(ms->snapshot_name[i_snapshot], name);

  PDM_get_current_memory(&ms->curr_real_mem[i_snapshot],
                         &ms->peak_real_mem[i_snapshot],
                         &ms->curr_virt_mem[i_snapshot],
                         &ms->peak_virt_mem[i_snapshot]);

  ms->curr_real_mem[i_snapshot] *= 1024; // Because orinaly in kB
  ms->peak_real_mem[i_snapshot] *= 1024;
  ms->curr_virt_mem[i_snapshot] *= 1024;
  ms->peak_virt_mem[i_snapshot] *= 1024;

}

void
PDM_memory_stats_log
(
 PDM_memory_stats_t* ms
)
{
  char curr_real[99];
  char peak_real[99];
  char curr_virt[99];
  char peak_virt[99];


  for(int i = 0; i < ms->n_memory_snapshot; ++i) {

    if(i == 0) {
      readable_fs(ms->curr_real_mem[i], curr_real);
      readable_fs(ms->peak_real_mem[i], peak_real);
      readable_fs(ms->curr_virt_mem[i], curr_virt);
      readable_fs(ms->peak_virt_mem[i], peak_virt);
    } else {
      readable_fs(ms->curr_real_mem[i] - ms->curr_real_mem[0], curr_real);
      readable_fs(ms->peak_real_mem[i] - ms->peak_real_mem[0], peak_real);
      readable_fs(ms->curr_virt_mem[i] - ms->curr_virt_mem[0], curr_virt);
      readable_fs(ms->peak_virt_mem[i] - ms->peak_virt_mem[0], peak_virt);

    }

    if(i == 0) {
      log_trace("%s : [curr_real = %s | peak_real = %s | curr_virt = %s | peak_virt = %s]\n", ms->snapshot_name[i],
                curr_real,
                peak_real,
                curr_virt,
                peak_virt);
    } else {
      log_trace("%s : delta [curr_real = %s | peak_real = %s | curr_virt = %s | peak_virt = %s]\n", ms->snapshot_name[i],
                curr_real,
                peak_real,
                curr_virt,
                peak_virt);

    }


  }


}


void
PDM_memory_stats_free
(
 PDM_memory_stats_t* ms
)
{
  for(int i = 0; i < ms->n_memory_snapshot; ++i) {
    if(ms->snapshot_name[i] != NULL) {
      free(ms->snapshot_name[i]);
    }
  }

  free(ms->snapshot_name);
  free(ms->curr_real_mem);
  free(ms->peak_real_mem);
  free(ms->curr_virt_mem);
  free(ms->peak_virt_mem);

  free(ms);
}


// TO BE CONSIDERED

// malloc_stats();

// struct mallinfo stat = mallinfo();
// printf("Non-mmapped space allocated (bytes)         : %i \n", stat.arena   );
// printf("Number of free chunks                       : %i \n", stat.ordblks );
// printf("Number of free fastbin blocks               : %i \n", stat.smblks  );
// printf("Number of mmapped regions                   : %i \n", stat.hblks   );
// printf("Space allocated in mmapped regions (bytes)  : %i \n", stat.hblkhd  );
// printf("See below                                   : %i \n", stat.usmblks );
// printf("Space in freed fastbin blocks (bytes)       : %i \n", stat.fsmblks );
// printf("Total allocated space (bytes)               : %i \n", stat.uordblks);
// printf("Total free space (bytes)                    : %i \n", stat.fordblks);
// printf("Top-most, releasable space (bytes)          : %i \n", stat.keepcost);


#ifdef __cplusplus
}
#endif /* __cplusplus */
