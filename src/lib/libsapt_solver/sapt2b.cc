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

#include "sapt2b.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT2B::SAPT2B(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : SAPT(options, psio, chkpt)
{
    get_workflow();
    get_calc_info();
}

SAPT2B::~SAPT2B()
{
    cleanup_calc_info();
    psio_->close(PSIF_SAPT_AA_DF_INTS,1);
    psio_->close(PSIF_SAPT_BB_DF_INTS,1);
    psio_->close(PSIF_SAPT_AB_DF_INTS,1);
    psio_->close(PSIF_SAPT_AMPS,1);
}

void SAPT2B::get_workflow()
{
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
      workflow_.Y3_ar = 0;
      workflow_.Y3_bs = 0;
    }
    else if (options_.get_str("SAPT_LEVEL") == "SCS_SAPT") {
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
      workflow_.Y3_ar = 0;
      workflow_.Y3_bs = 0;
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
      workflow_.Y3_ar = 0;
      workflow_.Y3_bs = 0;
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
      workflow_.Y3_ar = 0;
      workflow_.Y3_bs = 0;
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
      workflow_.Y3_ar = 0;
      workflow_.Y3_bs = 0;
    }
    else if (options_.get_str("SAPT_LEVEL") == "SAPT2+3") {
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
      workflow_.Y3_ar = 1;
      workflow_.Y3_bs = 1;
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

void SAPT2B::get_calc_info()
{
    calc_info_.nri = ribasis_->nbf();
    calc_info_.nrio = ribasis_->nbf() + 3;

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

    calc_info_.C = block_matrix(calc_info_.nso,calc_info_.nmo);
    psio_->read_entry(PSIF_SAPT_DIMER,"Dimer HF Coefficients",(char *)
      &(calc_info_.C[0][0]), sizeof(double)*calc_info_.nmo*calc_info_.nso);

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

void SAPT2B::cleanup_calc_info()
{
    free_block(calc_info_.C);
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
