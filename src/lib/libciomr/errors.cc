/*!
  \file errors.cc
  \brief Print error messages and abort for various errors
  \ingroup (CIOMR)
*/

#include "iomrparam.h"
#include "includes.h"
#include "types.h"

extern "C" {

/*!
** no_path_given(): Print error message for no path given and abort
**
** \param name = name of calling routine
**
** \ingroup (CIOMR)
*/ 
void no_path_given(char *name)
{
  fprintf(stderr,"%s: no path given\n",name);
  ioabort();
}

/*!
** malloc_check(): Check to see if malloc succeeded or failed.  If failure,
** print error and abort.
**
** \param caller = name of calling routine
** \param data = pointer to new data (supposed to be char *, not very
**               useful anymore...) 
** \ingroup (CIOMR)
*/
void malloc_check(char *caller, char *data)
{
  if (!data) {
    fprintf(stderr,"%s: malloc failed\n",caller);
    perror("malloc");
    ioabort();
    }
}

/*!
** fopen_check(): See if fopen worked; if not, print error and abort
** \param caller = name of calling routine
** \param path = path for fopen
** \param data = pointer for output stream (probably shouldn't really 
**               be char *)
** \ingroup (CIOMR)
*/
void fopen_check(char *caller, char *path, char *data)
{
  if (!data) {
    fprintf(stderr,"%s: fopen failed for %s\n",caller,path);
    perror("fopen");
    ioabort();
    }
}

/*!
** fread_error(): If error in fread, print error and abort
** \ingroup (CIOMR)
*/
void fread_error(char *caller)
{
  fprintf(stderr,"%s: fread failed\n",caller);
  perror("fread");
  ioabort();
}

/*!
** fwrite_error(): If error in fwrite, print error and abort
** \ingroup (CIOMR)
*/
void fwrite_error(char *caller)
{
  fprintf(stderr,"%s: fwrite failed\n",caller);
  perror("fwrite");
  ioabort();
}

} /* extern "C" */
