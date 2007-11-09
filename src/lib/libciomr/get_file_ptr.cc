/*!
** \file get_file_ptr.cc
** \ingroup (CIOMR)
*/

#include "includes.h"
#include "pointers.h"
 
extern "C" {

/*!
** GET_FILE_PTR():  Simply return the bytewise global file pointer for the
** specified unit.  Allows one to get this pointer without having to go to
** the trouble of declaring things like ptr, etc.
**
**   \param unit = file unit number
**
** Returns:
**   bytewise file pointer
**
** Daniel, January 1996
** \ingroup (CIOMR)
*/
PSI_FPTR get_file_ptr(int unit)
{
  return(ptr.wptr[unit]);
}

} /* extern "C" */
