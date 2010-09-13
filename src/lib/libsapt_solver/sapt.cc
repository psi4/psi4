/* 
 *  SAPT.CC 
 *
 */
#ifdef HAVE_MKL
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
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "sapt.h"
#include "structs.h"

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
    get_calc_info();
}

SAPT::~SAPT()
{
    cleanup_calc_info();
    psio_->close(PSIF_SAPT_AA_DF_INTS,1);
    psio_->close(PSIF_SAPT_BB_DF_INTS,1);
    psio_->close(PSIF_SAPT_AB_DF_INTS,1);
    psio_->close(PSIF_SAPT_AMPS,1);
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
    std::vector<int> realsA;
    realsA.push_back(0); 
    std::vector<int> ghostsA; 
    ghostsA.push_back(1); 
    std::vector<int> realsB;
    realsB.push_back(1); 
    std::vector<int> ghostsB; 
    ghostsB.push_back(0);
    shared_ptr<Molecule> monomerA = molecule_->extract_subsets(realsA,ghostsA);
    shared_ptr<Molecule> monomerB = molecule_->extract_subsets(realsB,ghostsB);
 
    params_.foccA = monomerA->nfrozen_core(options_.get_str("FREEZE_CORE"));
    params_.foccB = monomerB->nfrozen_core(options_.get_str("FREEZE_CORE"));
    //params_.foccA = options_.get_int("NFRZ_A");
    //params_.foccB = options_.get_int("NFRZ_B");

    //Natural Orbital Stuff
    params_.nat_orbs = options_.get_bool("NAT_ORBS");
    params_.occ_cutoff = options_.get_double("OCC_CUTOFF");

    if (options_.get_str("SAPT_LEVEL") == "SAPT0") {
      workflow_.W_ov = 1;
      workflow_.W_oo = 0;
      workflow_.W_vv = 0;
      workflow_.save_s = 0;
      workflow_.save_chf = 1;
      workflow_.save_Jmhalf = 0;
      workflow_.t_arar = 0;
      workflow_.theta_arar = 0;
      workflow_.t_bsbs = 0;
      workflow_.theta_bsbs = 0;
      workflow_.Y2_ar = 0;
      workflow_.Y2_bs = 0;
      workflow_.t_ar = 0;
      workflow_.t_bs = 0;
      workflow_.mp2_opdm = 0;
      workflow_.theta_ar_ar = 0;
      workflow_.theta_bs_bs = 0;
      workflow_.t2_arar = 0;
      workflow_.theta2_arar = 0;
      workflow_.t2_bsbs = 0;
      workflow_.theta2_bsbs = 0;
      workflow_.theta2_ar_ar = 0;
      workflow_.theta2_bs_bs = 0;
      workflow_.g_arar = 0;
      workflow_.g_bsbs = 0;
      workflow_.t_arbs = 0;
      workflow_.t_bsar = 0;
      workflow_.gt_ar_arbs = 0;
      workflow_.gt_bs_arbs = 0;
    }
    else if (options_.get_str("SAPT_LEVEL") == "SAPT_DFT") {
      workflow_.W_ov = 1;
      workflow_.W_oo = 0;
      workflow_.W_vv = 0;
      workflow_.save_s = 0;
      workflow_.save_chf = 1;
      workflow_.save_Jmhalf = 1;
      workflow_.t_arar = 0;
      workflow_.theta_arar = 0;
      workflow_.t_bsbs = 0;
      workflow_.theta_bsbs = 0;
      workflow_.Y2_ar = 0;
      workflow_.Y2_bs = 0;
      workflow_.t_ar = 0;
      workflow_.t_bs = 0;
      workflow_.mp2_opdm = 0;
      workflow_.theta_ar_ar = 0;
      workflow_.theta_bs_bs = 0;
      workflow_.t2_arar = 0;
      workflow_.theta2_arar = 0;
      workflow_.t2_bsbs = 0;
      workflow_.theta2_bsbs = 0;
      workflow_.theta2_ar_ar = 0;
      workflow_.theta2_bs_bs = 0;
      workflow_.g_arar = 0;
      workflow_.g_bsbs = 0;
      workflow_.t_arbs = 0;
      workflow_.t_bsar = 0;
      workflow_.gt_ar_arbs = 0;
      workflow_.gt_bs_arbs = 0;
    }
    else if (options_.get_str("SAPT_LEVEL") == "SAPT2") {
      workflow_.W_ov = 1;
      workflow_.W_oo = 1;
      workflow_.W_vv = 1;
      workflow_.save_s = 1;
      workflow_.save_chf = 1;
      workflow_.save_Jmhalf = 0;
      workflow_.t_arar = 1;
      workflow_.theta_arar = 1;
      workflow_.t_bsbs = 1;
      workflow_.theta_bsbs = 1;
      workflow_.Y2_ar = 1;
      workflow_.Y2_bs = 1;
      workflow_.t_ar = 1;
      workflow_.t_bs = 1;
      workflow_.mp2_opdm = 1;
      workflow_.theta_ar_ar = 1;
      workflow_.theta_bs_bs = 1;
      workflow_.t2_arar = 1;
      workflow_.theta2_arar = 1;
      workflow_.t2_bsbs = 1;
      workflow_.theta2_bsbs = 1;
      workflow_.theta2_ar_ar = 1;
      workflow_.theta2_bs_bs = 1;
      workflow_.g_arar = 0;
      workflow_.g_bsbs = 0;
      workflow_.t_arbs = 0;
      workflow_.t_bsar = 0;
      workflow_.gt_ar_arbs = 0;
      workflow_.gt_bs_arbs = 0;
    }
    else if (options_.get_str("SAPT_LEVEL") == "SAPT2+") {
      workflow_.W_ov = 1;
      workflow_.W_oo = 1;
      workflow_.W_vv = 1;
      workflow_.save_s = 1;
      workflow_.save_chf = 1;
      workflow_.save_Jmhalf = 0;
      workflow_.t_arar = 1;
      workflow_.theta_arar = 1;
      workflow_.t_bsbs = 1;
      workflow_.theta_bsbs = 1;
      workflow_.Y2_ar = 1;
      workflow_.Y2_bs = 1;
      workflow_.t_ar = 1;
      workflow_.t_bs = 1;
      workflow_.mp2_opdm = 1;
      workflow_.theta_ar_ar = 1;
      workflow_.theta_bs_bs = 1;
      workflow_.t2_arar = 1; 
      workflow_.theta2_arar = 1;
      workflow_.t2_bsbs = 1; 
      workflow_.theta2_bsbs = 1;
      workflow_.theta2_ar_ar = 1;
      workflow_.theta2_bs_bs = 1;
      workflow_.g_arar = 1;
      workflow_.g_bsbs = 1;
      workflow_.t_arbs = 1;
      workflow_.t_bsar = 1;
      workflow_.gt_ar_arbs = 1;
      workflow_.gt_bs_arbs = 1;
    }
    else if (options_.get_str("SAPT_LEVEL") == "SAPT2+3") {
    }
    else {
      fprintf(outfile,"Invalid SAPT level\n"); fflush(outfile);
      abort();
    }

    if (params_.nat_orbs) {
      workflow_.t_arar = 1;
      workflow_.theta_arar = 1;
      workflow_.t_bsbs = 1;
      workflow_.theta_bsbs = 1;
      workflow_.mp2_opdm = 1;
    }
}

void SAPT::get_ribasis()
{ 
    if (!options_.get_bool("NO_INPUT")) {
      ribasis_ = shared_ptr<BasisSet>(new BasisSet(chkpt_, "DF_BASIS_SAPT"));
    } 
    else {
      shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(
        options_.get_str("BASIS_PATH")));
      ribasis_ = BasisSet::construct(parser, molecule_, options_.get_str(
        "RI_BASIS_SAPT"));
    }
    zero_ = BasisSet::zero_basis_set();
    calc_info_.nri = ribasis_->nbf();
    calc_info_.nrio = ribasis_->nbf() + 3;
}

void SAPT::get_calc_info()
{
    psio_->open(PSIF_SAPT_DIMER,PSIO_OPEN_OLD);

    int errcod = 0;
    psio_->read_entry(PSIF_SAPT_DIMER,"Dimer NSO",(char *) &calc_info_.nso, 
      sizeof(int));
    psio_->read_entry(PSIF_SAPT_DIMER,"Dimer NMO",(char *) &calc_info_.nmo, 
      sizeof(int));
    psio_->read_entry(PSIF_SAPT_DIMER,"Dimer HF Energy",(char *) 
      &calc_info_.eHF_D, sizeof(double));
    psio_->read_entry(PSIF_SAPT_DIMER,"Dimer Nuclear Repulsion Energy",(char *) 
      &calc_info_.enuc_D, sizeof(double));
  
    calc_info_.nsotri = calc_info_.nso*(calc_info_.nso+1)/2;
    calc_info_.nmotri = calc_info_.nmo*(calc_info_.nmo+1)/2;
    calc_info_.ntei = calc_info_.nsotri*(calc_info_.nsotri+1)/2;
  
    /* Store overlap integrals */
    calc_info_.S = init_array(calc_info_.nsotri);
    psio_->read_entry(PSIF_SAPT_DIMER,"Dimer Overlap Integrals",(char *) 
      &calc_info_.S[0], sizeof(double)*calc_info_.nsotri);
  
    psio_->close(PSIF_SAPT_DIMER,1);
  
    calc_info_.ioff = (int *) malloc (calc_info_.nsotri * sizeof(int));
    calc_info_.index2i = (int *) malloc (calc_info_.nsotri * sizeof(int));
    calc_info_.index2j = (int *) malloc (calc_info_.nsotri * sizeof(int));
  
    calc_info_.ioff[0] = 0; // Create ioff array
  
    for (int i=1; i < calc_info_.nsotri; i++) {
      calc_info_.ioff[i] = calc_info_.ioff[i-1] + i;
    }
  
    for (int i=0; i<calc_info_.nso; i++) {
      for (int j=0; j<=i; j++) {
        calc_info_.index2i[INDEX(i,j)] = i;
        calc_info_.index2j[INDEX(i,j)] = j;
    }}
  
    psio_->open(PSIF_SAPT_MONOMERA,PSIO_OPEN_OLD);
  
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NSO",(char *) 
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NMO",(char *) 
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NOCC",(char *) 
      &calc_info_.noccA, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer NVIR",(char *) 
      &calc_info_.nvirA, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer Number of Electrons",(char *) 
      &calc_info_.NA, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Energy",(char *) 
      &calc_info_.eHF_A, sizeof(double));
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer Nuclear Repulsion Energy",
      (char *) &calc_info_.enuc_A, sizeof(double));
  
    calc_info_.evalsA = init_array(calc_info_.nmo);
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Eigenvalues",(char *) 
      &(calc_info_.evalsA[0]), sizeof(double)*calc_info_.nmo);
  
    calc_info_.VA = init_array(calc_info_.nsotri);
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer Nuclear Attraction Integrals",
      (char *) &(calc_info_.VA[0]), sizeof(double)*calc_info_.nsotri);
  
    calc_info_.CA = block_matrix(calc_info_.nso,calc_info_.nmo);
    psio_->read_entry(PSIF_SAPT_MONOMERA,"Monomer HF Coefficients",(char *) 
      &(calc_info_.CA[0][0]), sizeof(double)*calc_info_.nmo*calc_info_.nso);
  
    psio_->close(PSIF_SAPT_MONOMERA,1);
  
    psio_->open(PSIF_SAPT_MONOMERB,PSIO_OPEN_OLD);
  
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer NOCC",(char *) 
      &calc_info_.noccB, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer NVIR",(char *) 
      &calc_info_.nvirB, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer Number of Electrons",(char *) 
      &calc_info_.NB, sizeof(int));
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Energy",(char *) 
      &calc_info_.eHF_B, sizeof(double));
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer Nuclear Repulsion Energy",
      (char *) &calc_info_.enuc_B, sizeof(double));
  
    calc_info_.evalsB = init_array(calc_info_.nmo);
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Eigenvalues",(char *) 
      &(calc_info_.evalsB[0]), sizeof(double)*calc_info_.nmo);
  
    calc_info_.VB = init_array(calc_info_.nsotri);
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer Nuclear Attraction Integrals",
      (char *) &(calc_info_.VB[0]), sizeof(double)*calc_info_.nsotri);
  
    calc_info_.CB = block_matrix(calc_info_.nso,calc_info_.nmo);
    psio_->read_entry(PSIF_SAPT_MONOMERB,"Monomer HF Coefficients",(char *) 
      &(calc_info_.CB[0][0]), sizeof(double)*calc_info_.nmo*calc_info_.nso);
  
    psio_->close(PSIF_SAPT_MONOMERB,1);
}

void SAPT::cleanup_calc_info()
{ 
    free_block(calc_info_.CA);
    free_block(calc_info_.CB);
    free_block(calc_info_.S_AB);
    free_block(calc_info_.VABB);
    free_block(calc_info_.VAAB);
    free_block(calc_info_.VBAA);
    free_block(calc_info_.VBAB);
    if (workflow_.save_s) {
      free_block(calc_info_.sA);
      free_block(calc_info_.sB);
    }
    if (workflow_.save_chf) {
      free_block(calc_info_.CHFA);
      free_block(calc_info_.CHFB);
    }
    if (workflow_.W_oo) {
      free_block(calc_info_.WABB);
      free_block(calc_info_.WBAA);
    }
    if (workflow_.W_ov) {
      free_block(calc_info_.WABS);
      free_block(calc_info_.WBAR);
    }
    if (workflow_.W_vv) {
      free_block(calc_info_.WASS);
      free_block(calc_info_.WBRR);
    }
    if (workflow_.save_Jmhalf) {
      free_block(calc_info_.Jmhalf);
    }

    free(calc_info_.diagAA);
    free(calc_info_.diagBB);
    free(calc_info_.evalsA);
    free(calc_info_.evalsB);
    free(calc_info_.ioff);
    free(calc_info_.index2i);
    free(calc_info_.index2j);    
} 

}}
