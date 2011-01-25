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

#include "sapt3b.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace sapt {

SAPT3B::SAPT3B(Options& options, shared_ptr<PSIO> psio,
    shared_ptr<Chkpt> chkpt) : SAPT(options, psio, chkpt)
{
    get_calc_info();
}

SAPT3B::~SAPT3B()
{
    psio_->close(PSIF_3B_SAPT_AA_DF_INTS,1);
    psio_->close(PSIF_3B_SAPT_BB_DF_INTS,1);
    psio_->close(PSIF_3B_SAPT_CC_DF_INTS,1);
    psio_->close(PSIF_3B_SAPT_AMPS,1);
}

void SAPT3B::get_calc_info()
{
    calc_info_.nrio = ribasis_->nbf() + 3;

    /* Get trimer info */
    psio_->open(PSIF_3B_SAPT_TRIMER,PSIO_OPEN_OLD);

    psio_->read_entry(PSIF_3B_SAPT_TRIMER,"Trimer NSO",(char *)
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_TRIMER,"Trimer NMO",(char *)
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_TRIMER,"Trimer HF Energy",(char *)
      &calc_info_.eHF_ABC, sizeof(double));

    calc_info_.nsotri = calc_info_.nso*(calc_info_.nso+1)/2;
    calc_info_.nmotri = calc_info_.nmo*(calc_info_.nmo+1)/2;

    /* Store overlap integrals */
    calc_info_.S = init_array(calc_info_.nsotri);
    psio_->read_entry(PSIF_3B_SAPT_TRIMER,"Trimer Overlap Integrals", (char *)
      &calc_info_.S[0], sizeof(double)*calc_info_.nsotri);

    psio_->close(PSIF_3B_SAPT_TRIMER,1);

    calc_info_.ioff = init_int_array(calc_info_.nsotri);
    calc_info_.index2i = init_int_array(calc_info_.nsotri);
    calc_info_.index2j = init_int_array(calc_info_.nsotri);

    calc_info_.ioff[0] = 0; // Create ioff array

    for (int i=1; i < calc_info_.nsotri; i++) {
      calc_info_.ioff[i] = calc_info_.ioff[i-1] + i;
    }

    for (int i=0; i<calc_info_.nso; i++) {
      for (int j=0; j<=i; j++) {
        calc_info_.index2i[INDEX(i,j)] = i;
        calc_info_.index2j[INDEX(i,j)] = j;
      }}

    /* Get dimer AB info */
    psio_->open(PSIF_3B_SAPT_DIMER_AB,PSIO_OPEN_OLD);

    psio_->read_entry(PSIF_3B_SAPT_DIMER_AB,"Dimer NSO",(char *)
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_DIMER_AB,"Dimer NMO",(char *)
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_DIMER_AB,"Dimer HF Energy",(char *)
      &calc_info_.eHF_AB, sizeof(double));

    psio_->close(PSIF_3B_SAPT_DIMER_AB,1);

    /* Get dimer AC info */
    psio_->open(PSIF_3B_SAPT_DIMER_AC,PSIO_OPEN_OLD);

    psio_->read_entry(PSIF_3B_SAPT_DIMER_AC,"Dimer NSO",(char *)
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_DIMER_AC,"Dimer NMO",(char *)
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_DIMER_AC,"Dimer HF Energy",(char *)
      &calc_info_.eHF_AC, sizeof(double));

    psio_->close(PSIF_3B_SAPT_DIMER_AC,1);

    /* Get dimer BC info */
    psio_->open(PSIF_3B_SAPT_DIMER_BC,PSIO_OPEN_OLD);

    psio_->read_entry(PSIF_3B_SAPT_DIMER_BC,"Dimer NSO",(char *)
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_DIMER_BC,"Dimer NMO",(char *)
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_DIMER_BC,"Dimer HF Energy",(char *)
      &calc_info_.eHF_BC, sizeof(double));

    psio_->close(PSIF_3B_SAPT_DIMER_BC,1);

    /* Get monomer A info */
    psio_->open(PSIF_3B_SAPT_MONOMER_A,PSIO_OPEN_OLD);

    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,"Monomer NSO",(char *)
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,"Monomer NMO",(char *)
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,"Monomer NOCC",(char *)
      &calc_info_.noccA, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,"Monomer NVIR",(char *)
      &calc_info_.nvirA, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,"Monomer HF Energy",(char *)
      &calc_info_.eHF_A,sizeof(double));

    calc_info_.evalsA = init_array(calc_info_.nmo);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,"Monomer HF Eigenvalues",(char *)
      &(calc_info_.evalsA[0]),sizeof(double)*calc_info_.nmo);

    calc_info_.VA = init_array(calc_info_.nsotri);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,
      "Monomer Nuclear Attraction Integrals",(char *) &(calc_info_.VA[0]),
      sizeof(double)*calc_info_.nsotri);

    calc_info_.CA = block_matrix(calc_info_.nso,calc_info_.nmo);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_A,"Monomer HF Coefficients",(char *)
      &(calc_info_.CA[0][0]), sizeof(double)*calc_info_.nmo*calc_info_.nso);

    psio_->close(PSIF_3B_SAPT_MONOMER_A,1);

    /* Get monomer B info */
    psio_->open(PSIF_3B_SAPT_MONOMER_B,PSIO_OPEN_OLD);

    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,"Monomer NSO",(char *)
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,"Monomer NMO",(char *)
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,"Monomer NOCC",(char *)
      &calc_info_.noccB, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,"Monomer NVIR",(char *)
      &calc_info_.nvirB, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,"Monomer HF Energy",(char *)
      &calc_info_.eHF_B,sizeof(double));

    calc_info_.evalsB = init_array(calc_info_.nmo);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,"Monomer HF Eigenvalues",(char *)
      &(calc_info_.evalsB[0]),sizeof(double)*calc_info_.nmo);

    calc_info_.VB = init_array(calc_info_.nsotri);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,
      "Monomer Nuclear Attraction Integrals",(char *) &(calc_info_.VB[0]),
      sizeof(double)*calc_info_.nsotri);

    calc_info_.CB = block_matrix(calc_info_.nso,calc_info_.nmo);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_B,"Monomer HF Coefficients",(char *)
      &(calc_info_.CB[0][0]), sizeof(double)*calc_info_.nmo*calc_info_.nso);

    psio_->close(PSIF_3B_SAPT_MONOMER_B,1);

    /* Get monomer C info */
    psio_->open(PSIF_3B_SAPT_MONOMER_C,PSIO_OPEN_OLD);

    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,"Monomer NSO",(char *)
      &calc_info_.nso, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,"Monomer NMO",(char *)
      &calc_info_.nmo, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,"Monomer NOCC",(char *)
      &calc_info_.noccC, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,"Monomer NVIR",(char *)
      &calc_info_.nvirC, sizeof(int));
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,"Monomer HF Energy",(char *)
      &calc_info_.eHF_C,sizeof(double));

    calc_info_.evalsC = init_array(calc_info_.nmo);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,"Monomer HF Eigenvalues",(char *)
      &(calc_info_.evalsC[0]),sizeof(double)*calc_info_.nmo);

    calc_info_.VC = init_array(calc_info_.nsotri);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,
      "Monomer Nuclear Attraction Integrals",(char *) &(calc_info_.VC[0]),
      sizeof(double)*calc_info_.nsotri);

    calc_info_.CC = block_matrix(calc_info_.nso,calc_info_.nmo);
    psio_->read_entry(PSIF_3B_SAPT_MONOMER_C,"Monomer HF Coefficients",(char *)
      &(calc_info_.CC[0][0]), sizeof(double)*calc_info_.nmo*calc_info_.nso);

    psio_->close(PSIF_3B_SAPT_MONOMER_C,1);

    double e_abc = calc_info_.eHF_ABC - calc_info_.eHF_A - calc_info_.eHF_B -
      calc_info_.eHF_C;
    double e_ab = calc_info_.eHF_AB - calc_info_.eHF_A - calc_info_.eHF_B;
    double e_ac = calc_info_.eHF_AC - calc_info_.eHF_A - calc_info_.eHF_C;
    double e_bc = calc_info_.eHF_BC - calc_info_.eHF_B - calc_info_.eHF_C;
    double non_add = e_abc - e_ab - e_ac - e_bc;

    fprintf(outfile,"    Hartree-Fock Interaction Energies\n");
    fprintf(outfile,"   ***********************************\n");
    fprintf(outfile,"  E_ABC         %16.8lf mH %16.8lf kcal mol^-1\n",
      e_abc*1000.0,e_abc*627.5095);
    fprintf(outfile,"  E_AB          %16.8lf mH %16.8lf kcal mol^-1\n",
      e_ab*1000.0,e_ab*627.5095);
    fprintf(outfile,"  E_AC          %16.8lf mH %16.8lf kcal mol^-1\n",
      e_ac*1000.0,e_ac*627.5095);
    fprintf(outfile,"  E_BC          %16.8lf mH %16.8lf kcal mol^-1\n",
      e_bc*1000.0,e_bc*627.5095);
    fprintf(outfile,"  E_(ABC)       %16.8lf mH %16.8lf kcal mol^-1\n",
      non_add*1000.0,non_add*627.5095);
    fprintf(outfile,"\n");

    results_.e_HF = non_add;
}

void SAPT3B::cleanup_calc_info()
{
}

}}
