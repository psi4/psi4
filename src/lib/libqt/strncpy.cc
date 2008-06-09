/*!
  \file
  \brief Override strncpy to ensure strings are terminated
  \ingroup QT
  By Edward Valeev
*/

#include <cstring>
#include <libqt/qt.h>

namespace psi {

/*!
** strncpy(): Override default strncpy to ensure last byte is a string
**   terminating character.
**
** \param dest   = destination string
** \param source = source string
** \param n      = number of characters to copy
**
** Returns: pointer to destination string
*/
char* strncpy(char* dest, const char* source, size_t n) {
  if (n > 0) {
    ::strncpy(dest,source,n);
  }
  dest[n-1] = '\0';
  return dest;
}

}

