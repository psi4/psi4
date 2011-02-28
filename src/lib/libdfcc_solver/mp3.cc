#include "mp3.h"

#include <time.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libmints/mints.h>

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

MP3::MP3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  shared_ptr<DFTensor> df(new DFTensor(psio_, basisset_, ribasis_));
//df->form_MO_integrals((ULI)(0.9*(double)doubles_), C_aocc_, C_avir_, false, 
//  fitting_algorithm_, fitting_condition_, schwarz_cutoff_);
}

MP3::~MP3()
{
}

double MP3::compute_energy()
{
}

void MP3::print_header()
{
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*                       DF-MP3                         *\n");
    fprintf(outfile, "\t*    Third-Order Moller-Plesset Perturbation Theory    *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t*            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t*                                                      *\n");
    fprintf(outfile, "\t********************************************************\n");
    fprintf(outfile, "\n");
    CC::print_header();

}

}}
