#include <stdio.h>
#include "dlclose.h"
#include "pdm.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

int dlclose(void *handle) {
	PDM_UNUSED(handle);
	return 0;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
