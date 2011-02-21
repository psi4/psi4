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

  common_init();

  //shared_ptr<PseudoGrid> g(new PseudoGrid(basisset_->molecule(), "File"));
  //g->parse(options_.get_str("PS_GRID_FILE"));

  //shared_ptr<Pseudospectral> ps(new Pseudospectral(psio_, basisset_, dealias_, g));
  //shared_ptr<Matrix> I = ps->form_I();
  //I->print();

  //shared_ptr<LaplaceDenominator> denom(new LaplaceDenominator(evals_aocc_,evals_avir_,1.0E-6)); 
  //denom->debug();
  //shared_ptr<CholeskyDenominator> cdenom(new CholeskyDenominator(evals_aocc_,evals_avir_,1.0E-6)); 
  //cdenom->debug();

}
void MP2::common_init()
{

    mp2_algorithm_ = options_.get_str("MP2_ALGORITHM");

}

MP2::~MP2()
{
}

double MP2::compute_energy()
{
  if (mp2_algorithm_ == "DF")
    compute_DF_MP2();
  else if (mp2_algorithm_ == "SOS")
    compute_OS_MP2();
  else if (mp2_algorithm_ == "MOS")
    compute_OS_MP2();
  else if (mp2_algorithm_ == "PS")
    compute_PS_MP2();
  else if (mp2_algorithm_ == "PS2")
    compute_PS2_MP2();
  else if (mp2_algorithm_ == "PS3")
    compute_PS3_MP2();

  
  energies_["Opposite-Spin Energy"] = 0.5*energies_["MP2J Energy"]; 
  energies_["Same-Spin Energy"] = 0.5*energies_["MP2J Energy"] +  energies_["MP2K Energy"]; 
  energies_["Correlation Energy"] = energies_["MP2J Energy"] + energies_["MP2K Energy"];
  energies_["Total Energy"] = energies_["Reference Energy"] + energies_["Correlation Energy"];

  energies_["SCS Opposite-Spin Energy"] = 0.5*oss_*energies_["MP2J Energy"]; 
  energies_["SCS Same-Spin Energy"] = 0.5*sss_*energies_["MP2J Energy"] +  sss_*energies_["MP2K Energy"]; 
  energies_["SCS Correlation Energy"] = energies_["SCS Opposite-Spin Energy"] + energies_["SCS Same-Spin Energy"];
  energies_["SCS Total Energy"] = energies_["Reference Energy"] + energies_["SCS Correlation Energy"];

  return energies_["Total Energy"];
}
void MP2::compute_DF_MP2()
{
    shared_ptr<DFTensor> df(new DFTensor(psio_, basisset_, ribasis_));
    df->form_OV_integrals((ULI)(0.9*(double)doubles_), C_aocc_, C_avir_, true, fitting_algorithm_, fitting_condition_, schwarz_cutoff_);

    int nocc = df->nocc();
    int nvir = df->nvir();
    int naux = df->naux();   

}
void MP2::compute_OS_MP2()
{
    throw FeatureNotImplemented("libdfcc_solver", "psi::dfcc::MP2::compute_OS_MP2", __FILE__, __LINE__); 
}
void MP2::compute_PS_MP2()
{
    throw FeatureNotImplemented("libdfcc_solver", "psi::dfcc::MP2::compute_PS_MP2", __FILE__, __LINE__); 
}
void MP2::compute_PS2_MP2()
{
    throw FeatureNotImplemented("libdfcc_solver", "psi::dfcc::MP2::compute_PS2_MP2", __FILE__, __LINE__); 
}
void MP2::compute_PS3_MP2()
{
    throw FeatureNotImplemented("libdfcc_solver", "psi::dfcc::MP2::compute_PS3_MP2", __FILE__, __LINE__); 
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
