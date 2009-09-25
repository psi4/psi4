/*! \file
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/

namespace psi { namespace response {

double dot(double *A, double *B, int n)
{
  register int i;
  double tval;

  tval = 0.0;
  for(i=0; i < n; i++,A++,B++)
    tval += (*A) * (*B);

  return tval;
}

}} // namespace psi::response
