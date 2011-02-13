#include "mp2.h"

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

MP2::MP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
}

MP2::~MP2()
{
}

double MP2::compute_energy()
{
  return(0.0);
}

void MP2::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                        DF-MP2                        *\n");
    fprintf(outfile, "\t*    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t*        with Laplace and Pseudospectral Grids         *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");

    //CC::print_header();
            
}

}}
