/*!
** \file rfile.cc
** \ingroup (CIOMR)
*/

#include "includes.h"
#include "pointers.h"
#include <libciomr/libciomr.h>

extern "C" {

/*!
** rfile: open a binary file
**
** \param unit = file number
**
** \ingroup (CIOMR)
*/
void rfile(int unit)
{
       if (ptr.wptr == NULL) init_ptrs();

       sector = 1024;
       ptr.wptr[unit] = 0;
       ioopen_(&unit);
    }

} /* extern "C" */
