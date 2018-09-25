/*
 * $Id: ncio.c 5717 2016-01-12 11:35:24Z mourits $
 */

#if defined(_CRAY)
#   include "ffio.c"
#else
#   include "posixio.c"
#endif
