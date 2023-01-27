/* ID: base.c, last updated 2021-02-02, F.Osorio */

#include "base.h"

/* 'dims' functions */

DIMS
dims(int *pdims)
{ /* dims object */
  DIMS ans;

  ans = (DIMS) Calloc(1, DIMS_struct);
  ans->n = (int) pdims[0];
  ans->p = (int) pdims[1];
  return ans;
}

void
dims_free(DIMS this)
{ /* destructor for a dims object */
  Free(this);
}
