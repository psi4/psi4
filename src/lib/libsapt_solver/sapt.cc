/* 
 *  SAPT.CC 
 *
 */
#include "sapt.h"
#include "structs.h"

#ifdef _MKL
#include <mkl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/basisset.h>
#include <libmints/basisset_parser.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT::SAPT(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt)
{
    get_params();
    get_ribasis();
}

SAPT::~SAPT()
{
}

void SAPT::get_ribasis()
{
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    ribasis_ = BasisSet::construct(parser, molecule_, "RI_BASIS_SAPT");
    zero_ = BasisSet::zero_ao_basis_set();
}

void SAPT::get_params()
{
    //CPHF convergence parameters
    params_.e_conv = pow(10.0,-options_.get_int("E_CONVERGE"));
    params_.d_conv = pow(10.0,-options_.get_int("D_CONVERGE"));
    params_.maxiter = options_.get_int("MAXITER");
    params_.diisvec = options_.get_int("DIISVECS");

    //Print 
    params_.print = options_.get_int("PRINT");

    //Schwarz cutoff
    params_.schwarz = options_.get_double("SCHWARZ_CUTOFF");

    //Memory
    params_.memory = (long int) ((double) memory_ *
      options_.get_double("SAPT_MEM_SAFETY"));

  //Get Frozen Orbital Info
  if (options_["NFRZ_A"].has_changed() || options_["NFRZ_B"].has_changed() ||
    options_["NFRZ_B"].has_changed()) {
    params_.foccA = options_.get_int("NFRZ_A");
    params_.foccB = options_.get_int("NFRZ_B");
  }
  else {
    std::vector<int> realsA;
    realsA.push_back(0);
    std::vector<int> ghostsA;
    ghostsA.push_back(1);
    shared_ptr<Molecule> monomerA = molecule_->extract_subsets(realsA,
      ghostsA);
    params_.foccA = monomerA->nfrozen_core(options_.get_str("FREEZE_CORE"));

    std::vector<int> realsB;
    realsB.push_back(1);
    std::vector<int> ghostsB;
    ghostsB.push_back(0);
    shared_ptr<Molecule> monomerB = molecule_->extract_subsets(realsB,
      ghostsB);
    params_.foccB = monomerB->nfrozen_core(options_.get_str("FREEZE_CORE"));
  }  

    //Natural Orbital Stuff
    params_.nat_orbs = options_.get_bool("NAT_ORBS");
    params_.nat_orbs_t2 = options_.get_bool("NAT_ORBS_T2");
    params_.occ_cutoff = options_.get_double("OCC_CUTOFF");
}

}}
