/*! \defgroup CIOMR libciomr: The PSI I/O and Math Library */

/*!
  \file
  \brief Add two arrays
  \ingroup CIOMR
*/

namespace psi {
/*!
** add_arr(): Add arrays a and b and put the result in array c.  Adds
** the first n elements
**
** \param a = first array to add
** \param b = second array to add
** \param c = array to hold the result of a+b
** \param n = number of elements to add
**
** Returns: none
**
** \ingroup CIOMR
*/
void add_arr(double *a, double *b, double *c, int n)
{
  register int i;

  for (i=0; i < n; i++) {
    c[i] = a[i]+b[i];
  }
}

}

