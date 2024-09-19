#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "pdm.h"
#include "pdm_priv.h"

#include "pdm_mpi.h"
#include "pdm_config.h"

#include "pdm_mesh_adapt.h"


/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private function definitions
 *============================================================================*/

/**
 *
 * \brief  Main
 *
 */

int
main
(
 int argc,
 char *argv[]
)
{
  // check pdm_mesh_adapt.h compilation only
  PDM_UNUSED(argc);
  PDM_UNUSED(argv);

  return 0;
}
