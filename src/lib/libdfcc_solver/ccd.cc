#include "ccd.h"

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

CCD::CCD(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
}

CCD::~CCD()
{
}

double CCD::compute_energy()
{
  return(0.0);
}

void CCD::print_header()
{
     fprintf(outfile,"           CCD   \n");
     fprintf(outfile,"      Ed Hohenstein\n") ;
     fprintf(outfile,"     12 February 2011\n") ;
     fprintf(outfile,"\n");
     fprintf(outfile,"    Orbital Information\n");
     fprintf(outfile,"  -----------------------\n");
     fprintf(outfile,"    NSO      = %8d\n",nso_);
     fprintf(outfile,"    NMO      = %8d\n",nmo_);
     fprintf(outfile,"    NOCC Tot = %8d\n",nocc_);
     fprintf(outfile,"    NOCC Frz = %8d\n",nfocc_);
     fprintf(outfile,"    NOCC Act = %8d\n",naocc_);
     fprintf(outfile,"    NVIR Tot = %8d\n",nvir_);
     fprintf(outfile,"    NVIR Frz = %8d\n",nfvir_);
     fprintf(outfile,"    NVIR Act = %8d\n",navir_);
     fprintf(outfile,"    NDF      = %8d\n",ndf_);
     fprintf(outfile,"\n");
     fflush(outfile);
}

}}
