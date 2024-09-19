#include <stdio.h>
#include "dlclose.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

int dlclose(void *handle) {
	(void) (handle);
	return 0;
}

#ifdef __cplusplus
}
#endif /* __cplusplus */
