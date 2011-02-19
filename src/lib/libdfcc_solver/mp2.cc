#include "mp2.h"
#include <libmints/mints.h>
#include <lib3index/3index.h>

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

MP2::MP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  CC::print_header();

  shared_ptr<PseudoGrid> g(new PseudoGrid(basisset_->molecule(), "File"));
  g->parse(options_.get_str("PS_GRID_FILE"));

  shared_ptr<Pseudospectral> ps(new Pseudospectral(psio_, basisset_, dealias_, g));
  shared_ptr<Matrix> I = ps->form_I();
  I->print();

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
    fprintf(outfile, "\n");
    //CC::print_header();
            
}

}}
