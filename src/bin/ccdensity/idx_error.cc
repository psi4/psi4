/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <exception.h>
#include <psifiles.h>

namespace psi {
extern FILE* outfile;
namespace ccdensity {

void idx_error(const char *message, int p, int q, int r, int s, int pq, int rs,
               int pq_sym, int rs_sym, FILE *outfile)
{

  fprintf(outfile, "\n\tDPD Parameter Error in %s\n", message);
  fprintf(outfile,
          "\t-------------------------------------------------\n");
  fprintf(outfile,
          "\t    p      q      r      s  [   pq]  [   rs] pq_symm rs_symm\n");
  fprintf(outfile,"\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d   %1d\n", p,q,r,s,
          pq,rs,pq_sym,rs_sym);
  throw PsiException("ccdensity: error", __FILE__, __LINE__);
}


}} // namespace psi::ccdensity
