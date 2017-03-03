/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "usapt0.h"
#include "psi4/libfock/jk.h"
#include "psi4/libthce/thce.h"
#include "psi4/libthce/lreri.h"
#include "psi4/physconst.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/integral.h"

using namespace psi;
using namespace std;

namespace psi{ 

namespace sapt {

// TODO: ROHF orbitals for coupled induction
// TODO: SAPT charge transfer energy
USAPT0::USAPT0(SharedWavefunction d,
               SharedWavefunction mA,
               SharedWavefunction mB, Options& options, std::shared_ptr<PSIO> psio) : options_(options)
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
    bench_ = options_.get_int("BENCH");
    coupled_ind_ = options_.get_bool("COUPLED_INDUCTION");
// Get memory and convert it into words.
    memory_ = (unsigned long int)(Process::environment.get_memory() * options_.get_double("SAPT_MEM_FACTOR") * 0.125);

    cpks_maxiter_ = options_.get_int("MAXITER");
    cpks_delta_ = options_.get_double("D_CONVERGENCE");

    dimer_     = d->molecule();
    monomer_A_ = mA->molecule();
    monomer_B_ = mB->molecule();

    E_dimer_     = d->reference_energy();
    E_monomer_A_ = mA->reference_energy();
    E_monomer_B_ = mB->reference_energy();

    primary_   = d->basisset();
    primary_A_ = mA->basisset();
    primary_B_ = mB->basisset();

    mp2fit_ = d->get_basisset("DF_BASIS_SAPT");
    jkfit_ = d->get_basisset("DF_BASIS_SCF");

    if(options_.get_str("EXCH_SCALE_ALPHA") == "FALSE") {
        exch_scale_alpha_ = 0.0;
    } else if (options_.get_str("EXCH_SCALE_ALPHA") == "TRUE") {
        exch_scale_alpha_ = 1.0;                // Default value for alpha
    } else {
        exch_scale_alpha_ = std::atof(options_.get_str("EXCH_SCALE_ALPHA").c_str());
    }
    Process::environment.globals["SAPT ALPHA"] = exch_scale_alpha_;

//TODO: Modify this if all formulae are finally adapted for all bases.
    if (primary_A_->nbf() != primary_B_->nbf() || primary_->nbf() != primary_A_->nbf()) {
        throw PSIEXCEPTION("Monomer-centered bases not allowed in open-shell SAPT0");
    }

    initialize(mA,mB);
}

void USAPT0::initialize(SharedWavefunction mA, SharedWavefunction mB) {

    type_  = "USAPT0";

//    freq_points_ = options.get_int("FREQ_POINTS");
//    freq_scale_  = options.get_double("FREQ_SCALE");
//    freq_max_k_  = options.get_int("FREQ_MAX_K");

    Cocca_A_     = mA->Ca_subset("AO","OCC");
    Coccb_A_     = mA->Cb_subset("AO","OCC");
    Cvira_A_     = mA->Ca_subset("AO","VIR");
    Cvirb_A_     = mA->Cb_subset("AO","VIR");
    eps_occa_A_  = mA->epsilon_a_subset("AO","OCC");
    eps_occb_A_  = mA->epsilon_b_subset("AO","OCC");
    eps_vira_A_  = mA->epsilon_a_subset("AO","VIR");
    eps_virb_A_  = mA->epsilon_b_subset("AO","VIR");

    Cfocca_A_    = mA->Ca_subset("AO","FROZEN_OCC");
    Cfoccb_A_    = mA->Cb_subset("AO","FROZEN_OCC");
    Caocca_A_    = mA->Ca_subset("AO","ACTIVE_OCC");
    Caoccb_A_    = mA->Cb_subset("AO","ACTIVE_OCC");
    Cavira_A_    = mA->Ca_subset("AO","ACTIVE_VIR");
    Cavirb_A_    = mA->Cb_subset("AO","ACTIVE_VIR");
    Cfvira_A_    = mA->Ca_subset("AO","FROZEN_VIR");
    Cfvirb_A_    = mA->Cb_subset("AO","FROZEN_VIR");

    eps_focca_A_ = mA->epsilon_a_subset("AO","FROZEN_OCC");
    eps_foccb_A_ = mA->epsilon_b_subset("AO","FROZEN_OCC");
    eps_aocca_A_ = mA->epsilon_a_subset("AO","ACTIVE_OCC");
    eps_aoccb_A_ = mA->epsilon_b_subset("AO","ACTIVE_OCC");
    eps_avira_A_ = mA->epsilon_a_subset("AO","ACTIVE_VIR");
    eps_avirb_A_ = mA->epsilon_b_subset("AO","ACTIVE_VIR");
    eps_fvira_A_ = mA->epsilon_a_subset("AO","FROZEN_VIR");
    eps_fvirb_A_ = mA->epsilon_b_subset("AO","FROZEN_VIR");

    Cocca_B_     = mB->Ca_subset("AO","OCC");
    Coccb_B_     = mB->Cb_subset("AO","OCC");
    Cvira_B_     = mB->Ca_subset("AO","VIR");
    Cvirb_B_     = mB->Cb_subset("AO","VIR");
    eps_occa_B_  = mB->epsilon_a_subset("AO","OCC");
    eps_occb_B_  = mB->epsilon_b_subset("AO","OCC");
    eps_vira_B_  = mB->epsilon_a_subset("AO","VIR");
    eps_virb_B_  = mB->epsilon_b_subset("AO","VIR");

    Cfocca_B_    = mB->Ca_subset("AO","FROZEN_OCC");
    Cfoccb_B_    = mB->Cb_subset("AO","FROZEN_OCC");
    Caocca_B_    = mB->Ca_subset("AO","ACTIVE_OCC");
    Caoccb_B_    = mB->Cb_subset("AO","ACTIVE_OCC");
    Cavira_B_    = mB->Ca_subset("AO","ACTIVE_VIR");
    Cavirb_B_    = mB->Cb_subset("AO","ACTIVE_VIR");
    Cfvira_B_    = mB->Ca_subset("AO","FROZEN_VIR");
    Cfvirb_B_    = mB->Cb_subset("AO","FROZEN_VIR");

    eps_focca_B_ = mB->epsilon_a_subset("AO","FROZEN_OCC");
    eps_foccb_B_ = mB->epsilon_b_subset("AO","FROZEN_OCC");
    eps_aocca_B_ = mB->epsilon_a_subset("AO","ACTIVE_OCC");
    eps_aoccb_B_ = mB->epsilon_b_subset("AO","ACTIVE_OCC");
    eps_avira_B_ = mB->epsilon_a_subset("AO","ACTIVE_VIR");
    eps_avirb_B_ = mB->epsilon_b_subset("AO","ACTIVE_VIR");
    eps_fvira_B_ = mB->epsilon_a_subset("AO","FROZEN_VIR");
    eps_fvirb_B_ = mB->epsilon_b_subset("AO","FROZEN_VIR");

}
USAPT0::~USAPT0()
{
}

double USAPT0::compute_energy()
{
    energies_["HF"] = E_dimer_ - E_monomer_A_ - E_monomer_B_; // TODO: get dHF loaded correctly

    print_header();

    if (type_ == "USAPT0") {
        fock_terms();
        mp2_terms();
        print_trailer();
    } else {
        throw PSIEXCEPTION("USAPT: Unrecognized type");
    }

    return 0.0;
}
void USAPT0::print_header() const
{
    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\t                         SAPT                      \n");
    outfile->Printf( "\t               Rob Parrish and Ed Hohenstein             \n");
    outfile->Printf( "\t                Open-shell: Jérôme Gonthier              \n");
    outfile->Printf( "\t --------------------------------------------------------\n");
    outfile->Printf( "\n");

    outfile->Printf( "  ==> Sizes <==\n");
    outfile->Printf( "\n");

    outfile->Printf( "   => Resources <=\n\n");

    outfile->Printf( "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
    outfile->Printf( "\n");

    outfile->Printf( "   => Orbital Ranges <=\n\n");

    int nmoa_A = eps_focca_A_->dim() + eps_aocca_A_->dim() + eps_avira_A_->dim() + eps_fvira_A_->dim();
    int nmob_A = eps_foccb_A_->dim() + eps_aoccb_A_->dim() + eps_avirb_A_->dim() + eps_fvirb_A_->dim();
    int nmoa_B = eps_focca_B_->dim() + eps_aocca_B_->dim() + eps_avira_B_->dim() + eps_fvira_B_->dim();
    int nmob_B = eps_foccb_B_->dim() + eps_aoccb_B_->dim() + eps_avirb_B_->dim() + eps_fvirb_B_->dim();

    int nA = 0;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) != 0.0) nA++;
    }

    int nB = 0;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) != 0.0) nB++;
    }

    outfile->Printf( "    ------------------\n");
    outfile->Printf( "    %-6s %5s %5s\n", "Range", "M_A", "M_B");
    outfile->Printf( "    ------------------\n");
    outfile->Printf( "    %-6s %5d %5d\n", "natom", nA, nB);
    outfile->Printf( "    %-6s %5d %5d\n", "nso", primary_A_->nbf(), primary_B_->nbf());
    outfile->Printf( "    ------------------\n");
    outfile->Printf( "      Alpha orbitals  \n");
    outfile->Printf( "    ------------------\n");
    outfile->Printf( "    %-6s %5d %5d\n", "nmo", nmoa_A, nmoa_B);
    outfile->Printf( "    %-6s %5d %5d\n", "nocc", eps_aocca_A_->dim() + eps_focca_A_->dim(), eps_aocca_B_->dim() + eps_focca_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "nvir", eps_avira_A_->dim() + eps_fvira_A_->dim(), eps_avira_B_->dim() + eps_fvira_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "nfocc", eps_focca_A_->dim(), eps_focca_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "naocc", eps_aocca_A_->dim(), eps_aocca_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "navir", eps_avira_A_->dim(), eps_avira_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "nfvir", eps_fvira_A_->dim(), eps_fvira_B_->dim());
    outfile->Printf( "    ------------------\n");
    outfile->Printf( "      Beta orbitals  \n");
    outfile->Printf( "    ------------------\n");
    outfile->Printf( "    %-6s %5d %5d\n", "nmo", nmob_A, nmob_B);
    outfile->Printf( "    %-6s %5d %5d\n", "nocc", eps_aoccb_A_->dim() + eps_foccb_A_->dim(), eps_aoccb_B_->dim() + eps_foccb_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "nvir", eps_avirb_A_->dim() + eps_fvirb_A_->dim(), eps_avirb_B_->dim() + eps_fvirb_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "nfocc", eps_foccb_A_->dim(), eps_foccb_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "naocc", eps_aoccb_A_->dim(), eps_aoccb_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "navir", eps_avirb_A_->dim(), eps_avirb_B_->dim());
    outfile->Printf( "    %-6s %5d %5d\n", "nfvir", eps_fvirb_A_->dim(), eps_fvirb_B_->dim());
    outfile->Printf( "    ------------------\n");
    outfile->Printf( "\n");

    outfile->Printf( "   => Primary Basis Set <=\n\n");
    primary_->print_by_level("outfile", print_);

}

void USAPT0::print_trailer()
{
  // Note: sSAPT0 scaling is not implemented here since it was
  // never tested on open-shell systems.

  // The tolerance to scale exchange energies, i.e. if E_exch10 is
  // less than the scaling tolerance, we do not scale.
  double scaling_tol = 1.0e-5;

  double sapt_Xscal = ( energies_["Exch10"] < scaling_tol ? 1.0 : energies_["Exch10"] / energies_["Exch10(S^2)"] );
  if(exch_scale_alpha_ != 0.0) {
      sapt_Xscal = pow(sapt_Xscal, exch_scale_alpha_);
  }

  // Now we compute everything once without scaling, and then with scaling
  // if requested.
  std::vector<double> Xscal;
  Xscal.push_back(1.0);
  if(exch_scale_alpha_ != 0.0)
      Xscal.push_back(sapt_Xscal);

  // The main loop, computes everything with all scaling factors in
  // the Xscal vector. Only exports variables once, for the scaling factor
  // of 1.0.
  std::vector<double>::iterator scal_it;
  for(scal_it = Xscal.begin(); scal_it != Xscal.end(); ++scal_it) {

    energies_["delta HF,r (2)"] = 0.0;
    energies_["delta HF,u (2)"] = 0.0;
    if (energies_["HF"] != 0.0) {
        if (coupled_ind_) {
          energies_["delta HF,r (2)"] = energies_["HF"] - energies_["Elst10,r"] - energies_["Exch10"] - energies_["Ind20,r"] - *scal_it * energies_["Exch-Ind20,r"];
        } else {
          energies_["delta HF,u (2)"] = energies_["HF"] - energies_["Elst10,r"] - energies_["Exch10"] - energies_["Ind20,u"] - *scal_it * energies_["Exch-Ind20,u"];
        }
    }
    energies_["Electrostatics"] = energies_["Elst10,r"];
    energies_["Exchange"]       = energies_["Exch10"];
    if (coupled_ind_) {
        energies_["Induction"]      = energies_["Ind20,r"] + *scal_it * energies_["Exch-Ind20,r"] + energies_["delta HF,r (2)"];
    } else {
        energies_["Induction"]      = energies_["Ind20,u"] + *scal_it * energies_["Exch-Ind20,u"] + energies_["delta HF,u (2)"];
    }
    energies_["Dispersion"]     = energies_["Disp20"] + *scal_it * energies_["Exch-Disp20"];
    energies_["SAPT"]           = energies_["Electrostatics"] + energies_["Exchange"] + energies_["Induction"] + energies_["Dispersion"];

    if(scal_it == Xscal.begin()) {
        outfile->Printf("\n    SAPT Results \n");
    } else {
        outfile->Printf("\n    SAPT Results ==> ALL S2 TERMS SCALED (see Manual) <== \n");
        outfile->Printf("\n    Scaling factor (Exch10/Exch10(S^2))^{Alpha} = %12.6f\n", *scal_it);
        outfile->Printf("    with Alpha = %12.6f \n", exch_scale_alpha_);
    }
    std::string scaled = (scal_it != Xscal.begin() ? "sc." : "   ");
    outfile->Printf("  --------------------------------------------------------------------------------------------------------\n");
    outfile->Printf("    Electrostatics            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      energies_["Electrostatics"] * 1000.0,
      energies_["Electrostatics"] * pc_hartree2kcalmol,
      energies_["Electrostatics"] * pc_hartree2kJmol);
    outfile->Printf("      Elst10,r                %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      energies_["Elst10,r"] * 1000.0,
      energies_["Elst10,r"] * pc_hartree2kcalmol,
      energies_["Elst10,r"] * pc_hartree2kJmol);
    outfile->Printf("\n    Exchange %3s              %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      scaled.c_str(),
      energies_["Exchange"] * 1000.0,
      energies_["Exchange"] * pc_hartree2kcalmol,
      energies_["Exchange"] * pc_hartree2kJmol);
    outfile->Printf("      Exch10                  %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      energies_["Exch10"] * 1000.0,
      energies_["Exch10"] * pc_hartree2kcalmol,
      energies_["Exch10"] * pc_hartree2kJmol);
    outfile->Printf("      Exch10(S^2)             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      energies_["Exch10(S^2)"] * 1000.0,
      energies_["Exch10(S^2)"] * pc_hartree2kcalmol,
      energies_["Exch10(S^2)"] * pc_hartree2kJmol);
    outfile->Printf("\n    Induction %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      scaled.c_str(),
      energies_["Induction"] * 1000.0,
      energies_["Induction"] * pc_hartree2kcalmol,
      energies_["Induction"] * pc_hartree2kJmol);
    if (coupled_ind_) {
        outfile->Printf("      Ind20,r                 %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
          energies_["Ind20,r"] * 1000.0,
          energies_["Ind20,r"] * pc_hartree2kcalmol,
          energies_["Ind20,r"] * pc_hartree2kJmol);
        outfile->Printf("      Exch-Ind20,r %3s        %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
          scaled.c_str(),
          *scal_it * energies_["Exch-Ind20,r"] * 1000.0,
          *scal_it * energies_["Exch-Ind20,r"] * pc_hartree2kcalmol,
          *scal_it * energies_["Exch-Ind20,r"] * pc_hartree2kJmol);
        outfile->Printf("      delta HF,r (2) %3s      %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
          scaled.c_str(),
          energies_["delta HF,r (2)"] * 1000.0,
          energies_["delta HF,r (2)"] * pc_hartree2kcalmol,
          energies_["delta HF,r (2)"] * pc_hartree2kJmol);
    } else {
        outfile->Printf("      Ind20,u                 %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
          energies_["Ind20,u"] * 1000.0,
          energies_["Ind20,u"] * pc_hartree2kcalmol,
          energies_["Ind20,u"] * pc_hartree2kJmol);
        outfile->Printf("      Exch-Ind20,u %3s        %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
          scaled.c_str(),
          *scal_it * energies_["Exch-Ind20,u"] * 1000.0,
          *scal_it * energies_["Exch-Ind20,u"] * pc_hartree2kcalmol,
          *scal_it * energies_["Exch-Ind20,u"] * pc_hartree2kJmol);
        outfile->Printf("      delta HF,u (2) %3s      %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
          scaled.c_str(),
          energies_["delta HF,u (2)"] * 1000.0,
          energies_["delta HF,u (2)"] * pc_hartree2kcalmol,
          energies_["delta HF,u (2)"] * pc_hartree2kJmol);
    }
    outfile->Printf("\n    Dispersion %3s            %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      scaled.c_str(),
      energies_["Dispersion"] * 1000.0,
      energies_["Dispersion"] * pc_hartree2kcalmol,
      energies_["Dispersion"] * pc_hartree2kJmol);
    outfile->Printf("      Disp20                  %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      energies_["Disp20"] * 1000.0,
      energies_["Disp20"] * pc_hartree2kcalmol,
      energies_["Disp20"] * pc_hartree2kJmol);
    outfile->Printf("      Exch-Disp20 %3s         %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      scaled.c_str(),
      *scal_it * energies_["Exch-Disp20"] * 1000.0,
      *scal_it * energies_["Exch-Disp20"] * pc_hartree2kcalmol,
      *scal_it * energies_["Exch-Disp20"] * pc_hartree2kJmol);

    outfile->Printf("\n  Total HF                    %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      energies_["HF"] * 1000.0,
      energies_["HF"] * pc_hartree2kcalmol,
      energies_["HF"] * pc_hartree2kJmol);
    outfile->Printf("  Total SAPT0 %3s             %16.8lf [mEh] %16.8lf [kcal/mol] %16.8lf [kJ/mol]\n",
      scaled.c_str(),
      energies_["SAPT"] * 1000.0,
      energies_["SAPT"] * pc_hartree2kcalmol,
      energies_["SAPT"] * pc_hartree2kJmol);
    outfile->Printf("  --------------------------------------------------------------------------------------------------------\n");

    // Only export if not scaled.
    if(scal_it == Xscal.begin()) {

        Process::environment.globals["SAPT ELST ENERGY"] = energies_["Electrostatics"];
        Process::environment.globals["SAPT ELST10,R ENERGY"] = energies_["Elst10,r"];

        Process::environment.globals["SAPT EXCH ENERGY"] = energies_["Exchange"];
        Process::environment.globals["SAPT EXCH10 ENERGY"] = energies_["Exch10"];
        Process::environment.globals["SAPT EXCH10(S^2) ENERGY"] = energies_["Exch10(S^2)"];

        Process::environment.globals["SAPT IND ENERGY"] = energies_["Induction"];
        if (coupled_ind_) {
            Process::environment.globals["SAPT IND20,R ENERGY"] = energies_["Ind20,r"];
            Process::environment.globals["SAPT EXCH-IND20,R ENERGY"] = energies_["Exch-Ind20,r"];
        } else {
            // We still store in the R variants so the PsiVars machinery works.
            outfile->Printf("    WARNING: **Uncoupled** SAPT induction stored in SAPT IND20,R ENERGY \n");
            outfile->Printf("             and in SAPT EXCH-IND20,R ENERGY \n");
            Process::environment.globals["SAPT IND20,R ENERGY"] = energies_["Ind20,u"];
            Process::environment.globals["SAPT EXCH-IND20,R ENERGY"] = energies_["Exch-Ind20,u"];
            Process::environment.globals["SAPT IND20,U ENERGY"] = energies_["Ind20,u"];
            Process::environment.globals["SAPT EXCH-IND20,U ENERGY"] = energies_["Exch-Ind20,u"];
        }
        Process::environment.globals["SAPT HF TOTAL ENERGY"] = energies_["HF"];
        // Process::environment.globals["SAPT CT ENERGY"] = e_ind20_ + e_exch_ind20_;

        Process::environment.globals["SAPT DISP ENERGY"] = energies_["Dispersion"];
        Process::environment.globals["SAPT DISP20 ENERGY"] = energies_["Disp20"];
        Process::environment.globals["SAPT EXCH-DISP20 ENERGY"] = energies_["Exch-Disp20"];

        Process::environment.globals["SAPT0 TOTAL ENERGY"] = energies_["SAPT"];
        Process::environment.globals["SAPT ENERGY"] = energies_["SAPT"];
        Process::environment.globals["CURRENT ENERGY"] = Process::environment.globals["SAPT ENERGY"];

    }
  }

}

void USAPT0::fock_terms()
{
    outfile->Printf( "  SCF TERMS:\n\n");

    // ==> Setup <== //

    // => Compute the D matrices <= //

    std::shared_ptr<Matrix> Da_A    = Matrix::doublet(Cocca_A_, Cocca_A_,false,true);
    std::shared_ptr<Matrix> Db_A    = Matrix::doublet(Coccb_A_, Coccb_A_,false,true);
    std::shared_ptr<Matrix> Da_B    = Matrix::doublet(Cocca_B_, Cocca_B_,false,true);
    std::shared_ptr<Matrix> Db_B    = Matrix::doublet(Coccb_B_, Coccb_B_,false,true);

    // => Compute the P matrices <= //

    std::shared_ptr<Matrix> Pa_A    = Matrix::doublet(Cvira_A_, Cvira_A_,false,true);
    std::shared_ptr<Matrix> Pb_A    = Matrix::doublet(Cvirb_A_, Cvirb_A_,false,true);
    std::shared_ptr<Matrix> Pa_B    = Matrix::doublet(Cvira_B_, Cvira_B_,false,true);
    std::shared_ptr<Matrix> Pb_B    = Matrix::doublet(Cvirb_B_, Cvirb_B_,false,true);

    // => Compute the S matrix <= //

    std::shared_ptr<Matrix> S = build_S(primary_);

    // => Compute the V matrices <= //

    std::shared_ptr<Matrix> V_A = build_V(primary_A_);
    std::shared_ptr<Matrix> V_B = build_V(primary_B_);

    // => JK Object <= //

    std::shared_ptr<JK> jk = JK::build_JK(primary_, jkfit_, options_);

    // TODO: Recompute exactly how much memory is needed
    int naA = Cocca_A_->ncol();
    int nbA = Coccb_A_->ncol();
    int naB = Cocca_B_->ncol();
    int nbB = Coccb_B_->ncol();
    int nso = Cocca_A_->nrow();
    long int jk_memory = (long int)memory_;
    jk_memory -= 24 * nso * nso;
// Not sure why it should be 4, just taken from DFTSAPT code.
    jk_memory -=  4 * naA * nso;
    jk_memory -=  4 * nbA * nso;
    jk_memory -=  4 * naB * nso;
    jk_memory -=  4 * nbB * nso;
    if (jk_memory < 0L) {
        throw PSIEXCEPTION("Too little static memory for USAPT::fock_terms");
    }
    jk->set_memory((unsigned long int )jk_memory);

    // ==> Generalized Fock Source Terms [Elst/Exch] <== //

    // => Steric Interaction Density Terms (T) <= //
    std::shared_ptr<Matrix> Sija = build_Sija(S);
    std::shared_ptr<Matrix> Sijb = build_Sijb(S);
    std::shared_ptr<Matrix> Sija_n = build_Sij_n(Sija);
    Sija_n->set_name("Sija^inf (MO)");
    std::shared_ptr<Matrix> Sijb_n = build_Sij_n(Sijb);
    Sijb_n->set_name("Sijb^inf (MO)");

/* Build all of these matrices at once.
   I like that C_T_BA is actually CB T[BA] , so my naming convention
   is reversed compared to the DFTSAPT code.
   And AB quantities are not needed so we never compute them. */
    std::map<std::string, std::shared_ptr<Matrix> > Cbar_n = build_Cbar(Sija_n,Sijb_n);
    std::shared_ptr<Matrix> C_Ta_A_n  = Cbar_n["C_Ta_A"];
    std::shared_ptr<Matrix> C_Ta_B_n  = Cbar_n["C_Ta_B"];
    std::shared_ptr<Matrix> C_Ta_BA_n = Cbar_n["C_Ta_BA"];
    std::shared_ptr<Matrix> C_Tb_A_n  = Cbar_n["C_Tb_A"];
    std::shared_ptr<Matrix> C_Tb_B_n  = Cbar_n["C_Tb_B"];
    std::shared_ptr<Matrix> C_Tb_BA_n = Cbar_n["C_Tb_BA"];

    Sija.reset();
    Sijb.reset();
    Sija_n.reset();
    Sijb_n.reset();

//    Create matrices needed for the S^{2} formula as well

    std::shared_ptr<Matrix> S_CA = Matrix::doublet(S,Cocca_A_);
    std::shared_ptr<Matrix> Ca_AS = Matrix::doublet(Pa_B,S_CA);
//    Next one for DM-based formula (debug purposes for now)
    std::shared_ptr<Matrix> Ca_AB = Matrix::doublet(Da_B,S_CA);
    S_CA = Matrix::doublet(S,Coccb_A_);
    std::shared_ptr<Matrix> Cb_AS = Matrix::doublet(Pb_B,S_CA);
//    Next one for DM-based formula (debug purposes for now)
    std::shared_ptr<Matrix> Cb_AB = Matrix::doublet(Db_B,S_CA);
    S_CA.reset();

    // => Load the JK Object <= //

    std::vector<SharedMatrix>& Cl = jk->C_left();
    std::vector<SharedMatrix>& Cr = jk->C_right();
    Cl.clear();
    Cr.clear();

    // J/K[P^A]
    Cl.push_back(Cocca_A_);
    Cr.push_back(Cocca_A_);
    Cl.push_back(Coccb_A_);
    Cr.push_back(Coccb_A_);
    // J/K[T^B, S^\infty]
    Cl.push_back(C_Ta_B_n);
    Cr.push_back(Cocca_B_);
    Cl.push_back(C_Tb_B_n);
    Cr.push_back(Coccb_B_);
    // J/K[T^BA, S^\infty]
    Cl.push_back(C_Ta_BA_n);
    Cr.push_back(Cocca_A_);
    Cl.push_back(C_Tb_BA_n);
    Cr.push_back(Coccb_A_);

    // J/K[S_as]
    Cl.push_back(Cocca_A_);
    Cr.push_back(Ca_AS);
    Cl.push_back(Coccb_A_);
    Cr.push_back(Cb_AS);
    // J/K[S_ab] These next quantities relate to the DM formula for
    // E_{exch}^{(10)}(S^{2}), essentially for debug purposes
    Cl.push_back(Cocca_A_);
    Cr.push_back(Ca_AB);
    Cl.push_back(Coccb_A_);
    Cr.push_back(Cb_AB);
    // J/K[P^B]
    Cl.push_back(Cocca_B_);
    Cr.push_back(Cocca_B_);
    Cl.push_back(Coccb_B_);
    Cr.push_back(Coccb_B_);

    // => Initialize the JK object <= //

    jk->set_do_J(true);
    jk->set_do_K(true);
    jk->initialize();
    jk->print_header();

    // => Compute the JK matrices <= //

    jk->compute();

    // => Unload the JK Object <= //

    const std::vector<SharedMatrix>& J = jk->J();
    const std::vector<SharedMatrix>& K = jk->K();

    SharedMatrix Ja_A      = J[0];
    SharedMatrix Jb_A      = J[1];
    SharedMatrix J_Ta_B_n  = J[2];
    SharedMatrix J_Tb_B_n  = J[3];
    SharedMatrix J_Ta_BA_n = J[4];
    SharedMatrix J_Tb_BA_n = J[5];
    SharedMatrix Ja_B      = J[10];
    SharedMatrix Jb_B      = J[11];

    SharedMatrix Ka_A      = K[0];
    SharedMatrix Kb_A      = K[1];
    SharedMatrix K_Ta_B_n  = K[2];
    SharedMatrix K_Tb_B_n  = K[3];
    SharedMatrix K_Ta_BA_n = K[4];
    SharedMatrix K_Tb_BA_n = K[5];
    SharedMatrix Ka_AS     = K[6];
    SharedMatrix Kb_AS     = K[7];
    SharedMatrix Ka_AB     = K[8];
    SharedMatrix Kb_AB     = K[9];
    SharedMatrix Ka_B      = K[10];
    SharedMatrix Kb_B      = K[11];

    // ==> Electrostatic Terms <== //

    // Classical physics (watch for cancellation!)

    double Enuc = 0.0;
    Enuc += dimer_->nuclear_repulsion_energy();
    Enuc -= monomer_A_->nuclear_repulsion_energy();
    Enuc -= monomer_B_->nuclear_repulsion_energy();

    SharedMatrix D_A(Da_A->clone());
    D_A->add(Db_A);
    SharedMatrix D_B(Da_B->clone());
    D_B->add(Db_B);
    SharedMatrix El_pot_A(V_A->clone());
    El_pot_A->add(Ja_A);
    El_pot_A->add(Jb_A);
    double Elst10 = 0.0;
    std::vector<double> Elst10_terms;
    Elst10_terms.resize(3);
    Elst10_terms[0] += D_A->vector_dot(V_B);
    Elst10_terms[1] += D_B->vector_dot(El_pot_A);
    Elst10_terms[2] += Enuc;
    for (int k = 0; k < Elst10_terms.size(); k++) {
        Elst10 += Elst10_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Elst10_terms.size(); k++) {
            outfile->Printf("    Elst10,r (%1d)        = %18.12lf H\n",k+1,Elst10_terms[k]);
        }
    }
    D_A.reset();
    D_B.reset();
    energies_["Elst10,r"] = Elst10;
    outfile->Printf("    Elst10,r            = %18.12lf H\n",Elst10);


    // ==> Exchange Terms (S^\infty) <== //

    // => Compute the T matrices <= //

    std::shared_ptr<Matrix> Ta_A_n  = Matrix::doublet(C_Ta_A_n, Cocca_A_,false,true);
    std::shared_ptr<Matrix> Ta_B_n  = Matrix::doublet(C_Ta_B_n, Cocca_B_,false,true);
    std::shared_ptr<Matrix> Ta_BA_n = Matrix::doublet(C_Ta_BA_n,Cocca_A_,false,true);

    std::shared_ptr<Matrix> Tb_A_n  = Matrix::doublet(C_Tb_A_n, Coccb_A_,false,true);
    std::shared_ptr<Matrix> Tb_B_n  = Matrix::doublet(C_Tb_B_n, Coccb_B_,false,true);
    std::shared_ptr<Matrix> Tb_BA_n = Matrix::doublet(C_Tb_BA_n,Coccb_A_,false,true);

    C_Ta_A_n.reset();
    C_Ta_B_n.reset();
    C_Ta_BA_n.reset();

    C_Tb_A_n.reset();
    C_Tb_B_n.reset();
    C_Tb_BA_n.reset();

//    Here we compute some intermediates to reduce the number of matrix
//    multiplications
    std::shared_ptr<Matrix> El_pot_B(V_B->clone());
    El_pot_B->add(Ja_B);
    El_pot_B->add(Jb_B);
    std::shared_ptr<Matrix> ha_B(El_pot_B->clone());
    std::shared_ptr<Matrix> hb_B(El_pot_B->clone());
    ha_B->subtract(Ka_B);
    hb_B->subtract(Kb_B);
    std::shared_ptr<Matrix> ha_A(El_pot_A->clone());
    std::shared_ptr<Matrix> hb_A(El_pot_A->clone());
    ha_A->subtract(Ka_A);
    hb_A->subtract(Kb_A);
    std::shared_ptr<Matrix> Exch_pota_B(J_Ta_B_n->clone());
    Exch_pota_B->add(J_Tb_B_n);
    Exch_pota_B->add(J_Ta_BA_n);
    Exch_pota_B->add(J_Tb_BA_n);
    std::shared_ptr<Matrix> Exch_potb_B(Exch_pota_B->clone());
    Exch_pota_B->subtract(K_Ta_B_n);
    Exch_pota_B->subtract(K_Ta_BA_n);
    Exch_potb_B->subtract(K_Tb_B_n);
    Exch_potb_B->subtract(K_Tb_BA_n);

    std::shared_ptr<Matrix> inter_a(ha_B->clone());
    std::shared_ptr<Matrix> inter_b(hb_B->clone());
    inter_a->add(Exch_pota_B);
    inter_b->add(Exch_potb_B);
    double Exch10_n = 0.0;
    std::vector<double> Exch10_n_terms;
    Exch10_n_terms.resize(8);
    Exch10_n_terms[0] -= Da_A->vector_dot(Ka_B);
    Exch10_n_terms[1] -= Db_A->vector_dot(Kb_B);
    Exch10_n_terms[2] += Ta_B_n->vector_dot(ha_A);
    Exch10_n_terms[3] += Tb_B_n->vector_dot(hb_A);
    Exch10_n_terms[4] += Ta_A_n->vector_dot(inter_a);
    Exch10_n_terms[5] += Tb_A_n->vector_dot(inter_b);
    inter_a->add(ha_A);
    inter_b->add(hb_A);
    Exch10_n_terms[6] += Ta_BA_n->vector_dot(inter_a);
    Exch10_n_terms[7] += Tb_BA_n->vector_dot(inter_b);

    for (int k = 0; k < Exch10_n_terms.size(); k++) {
        Exch10_n += Exch10_n_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Exch10_n_terms.size(); k++) {
            outfile->Printf("    Exch10 (%1d)          = %18.12lf H\n",k,Exch10_n_terms[k]);
        }
    }
    energies_["Exch10"] = Exch10_n;
    outfile->Printf("    Exch10              = %18.12lf H\n",Exch10_n);

    inter_a.reset();
    inter_b.reset();
    Ta_A_n.reset();
    Ta_B_n.reset();
    Ta_BA_n.reset();
    Tb_A_n.reset();
    Tb_B_n.reset();
    Tb_BA_n.reset();
    Exch_pota_B.reset();
    Exch_potb_B.reset();


    // ==> Exchange Terms (S^2) <== //
//    First, in the second quantized formalism

    inter_a = std::shared_ptr<Matrix>(Matrix::triplet(El_pot_B,Da_A,S));
    inter_a->add(Ka_AS);
    inter_b = std::shared_ptr<Matrix>(Matrix::triplet(El_pot_B,Db_A,S));
    inter_b->add(Kb_AS);

    double Exch10_2 = 0.0;
    std::vector<double> Exch10_2_terms;
/*    Exch10_2_terms.resize(3);
    Exch10_2_terms[0] -= Matrix::triplet(Matrix::triplet(Da_A,S,Da_B),S,Pa_A)->vector_dot(V_B);
    Exch10_2_terms[0] -= Matrix::triplet(Matrix::triplet(Da_A,S,Da_B),S,Pa_A)->vector_dot(Ja_B);
    Exch10_2_terms[0] -= Matrix::triplet(Matrix::triplet(Da_A,S,Da_B),S,Pa_A)->vector_dot(Jb_B);
    Exch10_2_terms[1] -= Matrix::triplet(Matrix::triplet(Da_B,S,Da_A),S,Pa_B)->vector_dot(V_A);
    Exch10_2_terms[1] -= Matrix::triplet(Matrix::triplet(Da_B,S,Da_A),S,Pa_B)->vector_dot(Ja_A);
    Exch10_2_terms[1] -= Matrix::triplet(Matrix::triplet(Da_B,S,Da_A),S,Pa_B)->vector_dot(Jb_A);
    Exch10_2_terms[2] -= Matrix::triplet(Pa_A,S,Da_B)->vector_dot(Ka_AS);
    Exch10_2_terms[0] -= Matrix::triplet(Matrix::triplet(Db_A,S,Db_B),S,Pb_A)->vector_dot(V_B);
    Exch10_2_terms[0] -= Matrix::triplet(Matrix::triplet(Db_A,S,Db_B),S,Pb_A)->vector_dot(Ja_B);
    Exch10_2_terms[0] -= Matrix::triplet(Matrix::triplet(Db_A,S,Db_B),S,Pb_A)->vector_dot(Jb_B);
    Exch10_2_terms[1] -= Matrix::triplet(Matrix::triplet(Db_B,S,Db_A),S,Pb_B)->vector_dot(V_A);
    Exch10_2_terms[1] -= Matrix::triplet(Matrix::triplet(Db_B,S,Db_A),S,Pb_B)->vector_dot(Ja_A);
    Exch10_2_terms[1] -= Matrix::triplet(Matrix::triplet(Db_B,S,Db_A),S,Pb_B)->vector_dot(Jb_A);
    Exch10_2_terms[2] -= Matrix::triplet(Pb_A,S,Db_B)->vector_dot(Kb_AS);
*/    Exch10_2_terms.resize(4);
    Exch10_2_terms[0] -= Matrix::triplet(Pa_A,S,Da_B)->vector_dot(inter_a);
    Exch10_2_terms[1] -= Matrix::triplet(Pa_B,S,Da_A)->vector_dot(Matrix::triplet(El_pot_A,Da_B,S));
    Exch10_2_terms[2] -= Matrix::triplet(Pb_A,S,Db_B)->vector_dot(inter_b);
    Exch10_2_terms[3] -= Matrix::triplet(Pb_B,S,Db_A)->vector_dot(Matrix::triplet(El_pot_A,Db_B,S));
    for (int k = 0; k < Exch10_2_terms.size(); k++) {
        Exch10_2 += Exch10_2_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Exch10_2_terms.size(); k++) {
            outfile->Printf("    Exch10(S^2) (%1d)     = %18.12lf H\n",k+1,Exch10_2_terms[k]);
        }
    }
    energies_["Exch10(S^2)"] = Exch10_2;
    outfile->Printf("    Exch10(S^2)         = %18.12lf H\n",Exch10_2);
    outfile->Printf( "\n");

//    Now compute the same quantity but in the density matrix formalism,
//    valid for both dimer- and monomer-centered basis sets.

    inter_a = std::shared_ptr<Matrix>(Matrix::triplet(S,Da_B,El_pot_A));
    inter_a->scale(-1.0);
    inter_a->subtract(Matrix::triplet(El_pot_B,Da_A,S));
    inter_a->add(Ka_AB);
    inter_a->add(ha_A);
    inter_a->add(ha_B);

    inter_b = std::shared_ptr<Matrix>(Matrix::triplet(S,Db_B,El_pot_A));
    inter_b->scale(-1.0);
    inter_b->subtract(Matrix::triplet(El_pot_B,Db_A,S));
    inter_b->add(Kb_AB);
    inter_b->add(hb_A);
    inter_b->add(hb_B);

    double Exch10_2dm = 0.0;
    std::vector<double> Exch10_2dm_terms;
    Exch10_2dm_terms.resize(4);
    Exch10_2dm_terms[0] -= Da_A->vector_dot(Ka_B);
    Exch10_2dm_terms[1] -= Db_A->vector_dot(Kb_B);
    Exch10_2dm_terms[2] -= Matrix::triplet(Da_A,S,Da_B)->vector_dot(inter_a);
    Exch10_2dm_terms[3] -= Matrix::triplet(Db_A,S,Db_B)->vector_dot(inter_b);

/*    double Exch10_2dm = 0.0;
    std::vector<double> Exch10_2dm_terms;
    Exch10_2dm_terms.resize(6);
    Exch10_2dm_terms[0] += Matrix::triplet(Da_A,S,Da_B)->vector_dot(Ka_A);
    Exch10_2dm_terms[0] -= Matrix::triplet(Da_A,S,Da_B)->vector_dot(Ja_A);
    Exch10_2dm_terms[0] -= Matrix::triplet(Da_A,S,Da_B)->vector_dot(Jb_A);
    Exch10_2dm_terms[0] -= Matrix::triplet(Da_A,S,Da_B)->vector_dot(V_A);
    Exch10_2dm_terms[0] -= Matrix::triplet(Da_A,S,Da_B)->vector_dot(Ka_AB);
    Exch10_2dm_terms[0] += Matrix::triplet(Da_A,S,Matrix::triplet(Da_B,Ja_A,Da_B))->vector_dot(S);
    Exch10_2dm_terms[0] += Matrix::triplet(Da_A,S,Matrix::triplet(Da_B,Jb_A,Da_B))->vector_dot(S);
    Exch10_2dm_terms[0] += Matrix::triplet(Da_A,S,Matrix::triplet(Da_B,V_A,Da_B))->vector_dot(S);
    Exch10_2dm_terms[0] += Matrix::triplet(Da_A,S,Matrix::triplet(Da_B,S,Da_A))->vector_dot(Ja_B);
    Exch10_2dm_terms[0] += Matrix::triplet(Da_A,S,Matrix::triplet(Da_B,S,Da_A))->vector_dot(Jb_B);
    Exch10_2dm_terms[0] += Matrix::triplet(Da_A,S,Matrix::triplet(Da_B,S,Da_A))->vector_dot(V_B);
    Exch10_2dm_terms[1] += Matrix::triplet(Db_A,S,Db_B)->vector_dot(Kb_A);
    Exch10_2dm_terms[1] -= Matrix::triplet(Db_A,S,Db_B)->vector_dot(Jb_A);
    Exch10_2dm_terms[1] -= Matrix::triplet(Db_A,S,Db_B)->vector_dot(Ja_A);
    Exch10_2dm_terms[1] -= Matrix::triplet(Db_A,S,Db_B)->vector_dot(V_A);
    Exch10_2dm_terms[1] -= Matrix::triplet(Db_A,S,Db_B)->vector_dot(Kb_AB);
    Exch10_2dm_terms[1] += Matrix::triplet(Db_A,S,Matrix::triplet(Db_B,Jb_A,Db_B))->vector_dot(S);
    Exch10_2dm_terms[1] += Matrix::triplet(Db_A,S,Matrix::triplet(Db_B,Ja_A,Db_B))->vector_dot(S);
    Exch10_2dm_terms[1] += Matrix::triplet(Db_A,S,Matrix::triplet(Db_B,V_A,Db_B))->vector_dot(S);
    Exch10_2dm_terms[1] += Matrix::triplet(Db_A,S,Matrix::triplet(Db_B,S,Db_A))->vector_dot(Jb_B);
    Exch10_2dm_terms[1] += Matrix::triplet(Db_A,S,Matrix::triplet(Db_B,S,Db_A))->vector_dot(Ja_B);
    Exch10_2dm_terms[1] += Matrix::triplet(Db_A,S,Matrix::triplet(Db_B,S,Db_A))->vector_dot(V_B);
    Exch10_2dm_terms[2] += Matrix::triplet(Da_B,S,Da_A)->vector_dot(Ka_B);
    Exch10_2dm_terms[2] -= Matrix::triplet(Da_B,S,Da_A)->vector_dot(Jb_B);
    Exch10_2dm_terms[2] -= Matrix::triplet(Da_B,S,Da_A)->vector_dot(Ja_B);
    Exch10_2dm_terms[2] -= Matrix::triplet(Da_B,S,Da_A)->vector_dot(V_B);
    Exch10_2dm_terms[3] += Matrix::triplet(Db_B,S,Db_A)->vector_dot(Kb_B);
    Exch10_2dm_terms[3] -= Matrix::triplet(Db_B,S,Db_A)->vector_dot(Ja_B);
    Exch10_2dm_terms[3] -= Matrix::triplet(Db_B,S,Db_A)->vector_dot(Jb_B);
    Exch10_2dm_terms[3] -= Matrix::triplet(Db_B,S,Db_A)->vector_dot(V_B);
    Exch10_2dm_terms[4] -= Da_A->vector_dot(Ka_B);
    Exch10_2dm_terms[5] -= Db_A->vector_dot(Kb_B);
*/    for (int k = 0; k < Exch10_2dm_terms.size(); k++) {
        Exch10_2dm += Exch10_2dm_terms[k];
    }
    if (debug_) {
        for (int k = 0; k < Exch10_2dm_terms.size(); k++) {
            outfile->Printf("    Exch10(S^2) DM (%1d)     = %18.12lf H\n",k+1,Exch10_2dm_terms[k]);
        }
    }
    energies_["Exch10(S^2) DM"] = Exch10_2dm;
    outfile->Printf("    Exch10(S^2) DM      = %18.12lf H\n",Exch10_2dm);
    outfile->Printf( "\n");

    // Clear up some memory
    inter_a.reset();
    inter_b.reset();
    Ka_AS.reset();
    Kb_AS.reset();

    // ==> Uncorrelated Second-Order Response Terms [Induction] <== //

    // => ExchInd pertubations <= //

    std::shared_ptr<Matrix> C_Oa_B = Matrix::triplet(Da_A,S,Cocca_B_);
    std::shared_ptr<Matrix> C_Ob_B = Matrix::triplet(Db_A,S,Coccb_B_);
    std::shared_ptr<Matrix> C_Pa_B = Matrix::triplet(Matrix::triplet(Da_B,S,Da_A),S,Cocca_B_);
    std::shared_ptr<Matrix> C_Pb_B = Matrix::triplet(Matrix::triplet(Db_B,S,Db_A),S,Coccb_B_);
    std::shared_ptr<Matrix> C_Pa_A = Matrix::triplet(Matrix::triplet(Da_A,S,Da_B),S,Cocca_A_);
    std::shared_ptr<Matrix> C_Pb_A = Matrix::triplet(Matrix::triplet(Db_A,S,Db_B),S,Coccb_A_);

    Cl.clear();
    Cr.clear();

    // J/K[O]
    Cl.push_back(C_Oa_B);
    Cr.push_back(Cocca_B_);
    Cl.push_back(C_Ob_B);
    Cr.push_back(Coccb_B_);
    // J/K[P_B]
    Cl.push_back(C_Pa_B);
    Cr.push_back(Cocca_B_);
    Cl.push_back(C_Pb_B);
    Cr.push_back(Coccb_B_);
    // J/K[P_A]
    Cl.push_back(C_Pa_A);
    Cr.push_back(Cocca_A_);
    Cl.push_back(C_Pb_A);
    Cr.push_back(Coccb_A_);

    // => Compute the JK matrices <= //

    jk->compute();

    // => Unload the JK Object <= //

    std::shared_ptr<Matrix> Ja_O     = J[0];
    std::shared_ptr<Matrix> Jb_O     = J[1];
    std::shared_ptr<Matrix> J_Pa_B   = J[2];
    std::shared_ptr<Matrix> J_Pb_B   = J[3];
    std::shared_ptr<Matrix> J_Pa_A   = J[4];
    std::shared_ptr<Matrix> J_Pb_A   = J[5];

    std::shared_ptr<Matrix> Ka_O     = K[0];
    std::shared_ptr<Matrix> Kb_O     = K[1];

    // ==> Generalized ESP (Flat and Exchange) <== //

    std::map<std::string, std::shared_ptr<Matrix> > mapA;
    mapA["Cocc_A"] = Cocca_A_;
    mapA["Cvir_A"] = Cvira_A_;
    mapA["S"] = S;
    mapA["D_A"] = Da_A;
    mapA["El_pot_A"] = El_pot_A;
    mapA["h_A"] = ha_A;
    mapA["K_B"] = Ka_B;
    mapA["D_B"] = Da_B;
    mapA["El_pot_B"] = El_pot_B;
    mapA["h_B"] = ha_B;
    mapA["Ja_O"] = Ja_O;
    mapA["Jb_O"] = Jb_O;
    mapA["K_O"] = Ka_O;
    mapA["Ja_P"] = J_Pa_B;
    mapA["Jb_P"] = J_Pb_B;

    std::shared_ptr<Matrix> waB = build_ind_pot(mapA);
    std::shared_ptr<Matrix> uaB = build_exch_ind_pot(mapA);

    mapA["Cocc_A"] = Coccb_A_;
    mapA["Cvir_A"] = Cvirb_A_;
    mapA["S"] = S;
    mapA["D_A"] = Db_A;
    mapA["El_pot_A"] = El_pot_A;
    mapA["h_A"] = hb_A;
    mapA["K_B"] = Kb_B;
    mapA["D_B"] = Db_B;
    mapA["El_pot_B"] = El_pot_B;
    mapA["h_B"] = hb_B;
    mapA["Ja_O"] = Jb_O;
    mapA["Jb_O"] = Ja_O;
    mapA["K_O"] = Kb_O;
    mapA["Ja_P"] = J_Pb_B;
    mapA["Jb_P"] = J_Pa_B;

    std::shared_ptr<Matrix> wbB = build_ind_pot(mapA);
    std::shared_ptr<Matrix> ubB = build_exch_ind_pot(mapA);

    Ka_O->transpose_this();
    Kb_O->transpose_this();

    std::map<std::string, std::shared_ptr<Matrix> > mapB;
    mapB["Cocc_A"] = Cocca_B_;
    mapB["Cvir_A"] = Cvira_B_;
    mapB["S"] = S;
    mapB["D_A"] = Da_B;
    mapB["El_pot_A"] = El_pot_B;
    mapB["h_A"] = ha_B;
    mapB["K_B"] = Ka_A;
    mapB["D_B"] = Da_A;
    mapB["El_pot_B"] = El_pot_A;
    mapB["h_B"] = ha_A;
    mapB["Ja_O"] = Ja_O;
    mapB["Jb_O"] = Jb_O;
    mapB["K_O"] = Ka_O;
    mapB["Ja_P"] = J_Pa_A;
    mapB["Jb_P"] = J_Pb_A;

    std::shared_ptr<Matrix> waA = build_ind_pot(mapB);
    std::shared_ptr<Matrix> uaA = build_exch_ind_pot(mapB);

    mapB["Cocc_A"] = Coccb_B_;
    mapB["Cvir_A"] = Cvirb_B_;
    mapB["S"] = S;
    mapB["D_A"] = Db_B;
    mapB["El_pot_A"] = El_pot_B;
    mapB["h_A"] = hb_B;
    mapB["K_B"] = Kb_A;
    mapB["D_B"] = Db_A;
    mapB["El_pot_B"] = El_pot_A;
    mapB["h_B"] = hb_A;
    mapB["Ja_O"] = Jb_O;
    mapB["Jb_O"] = Ja_O;
    mapB["K_O"] = Kb_O;
    mapB["Ja_P"] = J_Pb_A;
    mapB["Jb_P"] = J_Pa_A;

    std::shared_ptr<Matrix> wbA = build_ind_pot(mapB);
    std::shared_ptr<Matrix> ubA = build_exch_ind_pot(mapB);

    Ka_O->transpose_this();
    Kb_O->transpose_this();

    // ==> Uncoupled Induction <== //

    std::shared_ptr<Matrix> xuaA(waB->clone());
    std::shared_ptr<Matrix> xuaB(waA->clone());
    std::shared_ptr<Matrix> xubA(wbB->clone());
    std::shared_ptr<Matrix> xubB(wbA->clone());

    {
        int naa = eps_occa_A_->dimpi()[0];
        int nab = eps_occa_B_->dimpi()[0];
        int nar = eps_vira_A_->dimpi()[0];
        int nas = eps_vira_B_->dimpi()[0];
        int nba = eps_occb_A_->dimpi()[0];
        int nbb = eps_occb_B_->dimpi()[0];
        int nbr = eps_virb_A_->dimpi()[0];
        int nbs = eps_virb_B_->dimpi()[0];

        double** xuaAp = xuaA->pointer();
        double** xuaBp = xuaB->pointer();
        double*  eaap = eps_occa_A_->pointer();
        double*  erap = eps_vira_A_->pointer();
        double*  ebap = eps_occa_B_->pointer();
        double*  esap = eps_vira_B_->pointer();
        double** xubAp = xubA->pointer();
        double** xubBp = xubB->pointer();
        double*  eabp = eps_occb_A_->pointer();
        double*  erbp = eps_virb_A_->pointer();
        double*  ebbp = eps_occb_B_->pointer();
        double*  esbp = eps_virb_B_->pointer();

        for (int a = 0; a < naa; a++) {
            for (int r = 0; r < nar; r++) {
                xuaAp[a][r] = xuaAp[a][r] / (eaap[a] - erap[r]);
            }
        }
        for (int a = 0; a < nba; a++) {
            for (int r = 0; r < nbr; r++) {
                xubAp[a][r] = xubAp[a][r] / (eabp[a] - erbp[r]);
            }
        }

        for (int b = 0; b < nab; b++) {
            for (int s = 0; s < nas; s++) {
                xuaBp[b][s] = xuaBp[b][s] / (ebap[b] - esap[s]);
            }
        }
        for (int b = 0; b < nbb; b++) {
            for (int s = 0; s < nbs; s++) {
                xubBp[b][s] = xubBp[b][s] / (ebbp[b] - esbp[s]);
            }
        }
    }

    // ==> Induction <== //

    double Ind20u_AB = xuaA->vector_dot(waB);
    Ind20u_AB += xubA->vector_dot(wbB);
    double Ind20u_BA = xuaB->vector_dot(waA);
    Ind20u_BA += xubB->vector_dot(wbA);
    double Ind20u = Ind20u_AB + Ind20u_BA;
    energies_["Ind20,u (A<-B)"] = Ind20u_AB;
    energies_["Ind20,u (B<-A)"] = Ind20u_BA;
    energies_["Ind20,u"] = Ind20u;
    outfile->Printf("    Ind20,u (A<-B)      = %18.12lf H\n",Ind20u_AB);
    outfile->Printf("    Ind20,u (B<-A)      = %18.12lf H\n",Ind20u_BA);
    outfile->Printf("    Ind20,u             = %18.12lf H\n",Ind20u);

    // => Exchange-Induction <= //

    double ExchInd20u_AB = xuaA->vector_dot(uaB);
    ExchInd20u_AB += xubA->vector_dot(ubB);
    double ExchInd20u_BA = xuaB->vector_dot(uaA);
    ExchInd20u_BA += xubB->vector_dot(ubA);
    double ExchInd20u = ExchInd20u_AB + ExchInd20u_BA;
    outfile->Printf("    Exch-Ind20,u (A<-B) = %18.12lf H\n",ExchInd20u_AB);
    outfile->Printf("    Exch-Ind20,u (B<-A) = %18.12lf H\n",ExchInd20u_BA);
    outfile->Printf("    Exch-Ind20,u        = %18.12lf H\n",ExchInd20u);
    outfile->Printf("\n");

    energies_["Exch-Ind20,u (A<-B)"] = ExchInd20u_AB;
    energies_["Exch-Ind20,u (B<-A)"] = ExchInd20u_BA;
    energies_["Exch-Ind20,u"] = ExchInd20u_AB + ExchInd20u_BA;

    if (coupled_ind_) {
        // => Coupled Induction <= //
        // TODO: Write an RO CPHF solver for ROHF reference orbitals
    
        // Compute CPKS for both alpha and beta electrons.
        std::map < std::string, std::shared_ptr<Matrix> > x_sol = compute_x(jk,waB,wbB,waA,wbA);
        std::shared_ptr<Matrix> xaA = x_sol["Aa"];
        std::shared_ptr<Matrix> xbA = x_sol["Ab"];
        std::shared_ptr<Matrix> xaB = x_sol["Ba"];
        std::shared_ptr<Matrix> xbB = x_sol["Bb"];
    
        // Backward in Ed's convention
        xaA->scale(-1.0);
        xaB->scale(-1.0);
        xbA->scale(-1.0);
        xbB->scale(-1.0);
    
        // => Induction <= //
    
        double Ind20r_AB = xaA->vector_dot(waB);
        Ind20r_AB += xbA->vector_dot(wbB);
        double Ind20r_BA = xaB->vector_dot(waA);
        Ind20r_BA += xbB->vector_dot(wbA);
        double Ind20r = Ind20r_AB + Ind20r_BA;
        energies_["Ind20,r (A<-B)"] = Ind20r_AB;
        energies_["Ind20,r (B<-A)"] = Ind20r_BA;
        energies_["Ind20,r"] = Ind20r;
        outfile->Printf("    Ind20,r (A<-B)      = %18.12lf H\n",Ind20r_AB);
        outfile->Printf("    Ind20,r (B<-A)      = %18.12lf H\n",Ind20r_BA);
        outfile->Printf("    Ind20,r             = %18.12lf H\n",Ind20r);
    
        // => Exchange-Induction <= //
    
        double ExchInd20r_AB = xaA->vector_dot(uaB);
        ExchInd20r_AB += xbA->vector_dot(ubB);
        double ExchInd20r_BA = xaB->vector_dot(uaA);
        ExchInd20r_BA += xbB->vector_dot(ubA);
        double ExchInd20r = ExchInd20r_AB + ExchInd20r_BA;
        outfile->Printf("    Exch-Ind20,r (A<-B) = %18.12lf H\n",ExchInd20r_AB);
        outfile->Printf("    Exch-Ind20,r (B<-A) = %18.12lf H\n",ExchInd20r_BA);
        outfile->Printf("    Exch-Ind20,r        = %18.12lf H\n",ExchInd20r);
        outfile->Printf("\n");
    
        energies_["Exch-Ind20,r (A<-B)"] = ExchInd20r_AB;
        energies_["Exch-Ind20,r (B<-A)"] = ExchInd20r_BA;
        energies_["Exch-Ind20,r"] = ExchInd20r_AB + ExchInd20r_BA;

    }

    vars_["S"]   = S;
    vars_["Da_A"] = Da_A;
    vars_["Db_A"] = Db_A;
    vars_["El_pot_A"] = El_pot_A;
    vars_["ha_A"] = ha_A;
    vars_["hb_A"] = hb_A;
    vars_["Da_B"] = Da_B;
    vars_["Db_B"] = Db_B;
    vars_["El_pot_B"] = El_pot_B;
    vars_["ha_B"] = ha_B;
    vars_["hb_B"] = hb_B;
    vars_["Ka_O"] = Ka_O;
    vars_["Kb_O"] = Kb_O;
}
std::shared_ptr<Matrix> USAPT0::build_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars)
{
    std::shared_ptr<Matrix> Ca = vars["Cocc_A"];
    std::shared_ptr<Matrix> Cr = vars["Cvir_A"];
    std::shared_ptr<Matrix> El_pot = vars["El_pot_B"];

    return Matrix::triplet(Ca,El_pot,Cr,true,false,false);
}
std::shared_ptr<Matrix> USAPT0::build_exch_ind_pot(std::map<std::string, std::shared_ptr<Matrix> >& vars)
{
    // By convention, a denotes the spin of the D_A density matrix, and
    // b the opposite spin, since this routine applies to both spins.
    std::shared_ptr<Matrix> Ca = vars["Cocc_A"];
    std::shared_ptr<Matrix> Cr = vars["Cvir_A"];

    std::shared_ptr<Matrix> S = vars["S"];

    std::shared_ptr<Matrix> D_A = vars["D_A"];
    std::shared_ptr<Matrix> El_pot_A = vars["El_pot_A"];
    std::shared_ptr<Matrix> Jb_A = vars["Jb_A"];
    std::shared_ptr<Matrix> h_A = vars["h_A"];
    std::shared_ptr<Matrix> K_B = vars["K_B"];
    std::shared_ptr<Matrix> D_B = vars["D_B"];
    std::shared_ptr<Matrix> El_pot_B = vars["El_pot_B"];
    std::shared_ptr<Matrix> h_B = vars["h_B"];

    std::shared_ptr<Matrix> Ja_O = vars["Ja_O"];
    std::shared_ptr<Matrix> Jb_O = vars["Jb_O"]; // J[D^A S D^B]
    std::shared_ptr<Matrix> K_O = vars["K_O"]; // K[D^A S D^B]
    std::shared_ptr<Matrix> Ja_P = vars["Ja_P"]; // J[D^B S D^A S D^B]
    std::shared_ptr<Matrix> Jb_P = vars["Jb_P"]; // J[D^B S D^A S D^B]

    std::shared_ptr<Matrix> W(K_B->clone());
    std::shared_ptr<Matrix> T;

    // 1
    W->scale(-1.0);

    // 2
    W->subtract(Ja_O);

    //3
    W->subtract(Jb_O);

    //4
    W->add(K_O);

    //5
    W->add(Ja_P);

    //6
    W->add(Jb_P);

    //7 Use T to compute intermediate

    T = h_A->clone();
    T->scale(-1.0);
    T->add(Matrix::triplet(S,D_A,El_pot_B));
    T->add(Matrix::triplet(El_pot_A,D_B,S));
    T->subtract(K_O->transpose());

    W->add(Matrix::triplet(S,D_B,T));

    //8 Use again T for intermediate

    T = h_B->clone();
    T->scale(-1.0);
    T->add(Matrix::triplet(El_pot_B,D_A,S));
    T->subtract(K_O);

    W->add(Matrix::triplet(T,D_B,S));

    return Matrix::triplet(Ca,W,Cr,true,false,false);
}
std::shared_ptr<Matrix> USAPT0::build_S(std::shared_ptr<BasisSet> basis)
{
    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(basis));
    std::shared_ptr<OneBodyAOInt> Sint(factory->ao_overlap());
    std::shared_ptr<Matrix> S(new Matrix("S (AO)", basis->nbf(), basis->nbf()));
    Sint->compute(S);
    return S;
}
std::shared_ptr<Matrix> USAPT0::build_V(std::shared_ptr<BasisSet> basis)
{
    std::shared_ptr<IntegralFactory> factory(new IntegralFactory(basis));
    std::shared_ptr<OneBodyAOInt> Sint(factory->ao_potential());
    std::shared_ptr<Matrix> S(new Matrix("V (AO)", basis->nbf(), basis->nbf()));
    Sint->compute(S);
    return S;
}
std::shared_ptr<Matrix> USAPT0::build_Sija(std::shared_ptr<Matrix> S)
{
    int nso = Cocca_A_->nrow();
    int nocc_A = Cocca_A_->ncol();
    int nocc_B = Cocca_B_->ncol();
    int nocc = nocc_A + nocc_B;

    std::shared_ptr<Matrix> Sij(new Matrix("Sija (MO)", nocc, nocc));
    std::shared_ptr<Matrix> T(new Matrix("T", nso, nocc_B));

    double** Sp = S->pointer();
    double** Tp = T->pointer();
    double** Sijp = Sij->pointer();
    double** CAp = Cocca_A_->pointer();
    double** CBp = Cocca_B_->pointer();

    C_DGEMM('N','N',nso,nocc_B,nso,1.0,Sp[0],nso,CBp[0],nocc_B,0.0,Tp[0],nocc_B);
    C_DGEMM('T','N',nocc_A,nocc_B,nso,1.0,CAp[0],nocc_A,Tp[0],nocc_B,0.0,&Sijp[0][nocc_A],nocc);

    Sij->copy_upper_to_lower();

    return Sij;
}
std::shared_ptr<Matrix> USAPT0::build_Sijb(std::shared_ptr<Matrix> S)
{
    int nso = Coccb_A_->nrow();
    int nocc_A = Coccb_A_->ncol();
    int nocc_B = Coccb_B_->ncol();
    int nocc = nocc_A + nocc_B;

    std::shared_ptr<Matrix> Sij(new Matrix("Sijb (MO)", nocc, nocc));
    std::shared_ptr<Matrix> T(new Matrix("T", nso, nocc_B));

    double** Sp = S->pointer();
    double** Tp = T->pointer();
    double** Sijp = Sij->pointer();
    double** CAp = Coccb_A_->pointer();
    double** CBp = Coccb_B_->pointer();

    C_DGEMM('N','N',nso,nocc_B,nso,1.0,Sp[0],nso,CBp[0],nocc_B,0.0,Tp[0],nocc_B);
    C_DGEMM('T','N',nocc_A,nocc_B,nso,1.0,CAp[0],nocc_A,Tp[0],nocc_B,0.0,&Sijp[0][nocc_A],nocc);

    Sij->copy_upper_to_lower();

    return Sij;
}
// TOMODIF - spin
std::shared_ptr<Matrix> USAPT0::build_Sij_n(std::shared_ptr<Matrix> Sij)
{
    int nocc = Sij->nrow();

    std::shared_ptr<Matrix> Sij2(new Matrix("Sij^inf (MO)", nocc, nocc));

    double** Sijp = Sij->pointer();
    double** Sij2p = Sij2->pointer();

    Sij2->copy(Sij);
    for (int i = 0; i < nocc; i++) {
        Sij2p[i][i] = 1.0;
    }

    int info;

    info = C_DPOTRF('L',nocc,Sij2p[0],nocc);
    if (info) {
        throw PSIEXCEPTION("Sij DPOTRF failed. How far up the steric wall are you?");
    }

    info = C_DPOTRI('L',nocc,Sij2p[0],nocc);
    if (info) {
        throw PSIEXCEPTION("Sij DPOTRI failed. How far up the steric wall are you?");
    }

    Sij2->copy_upper_to_lower();

    for (int i = 0; i < nocc; i++) {
        Sij2p[i][i] -= 1.0;
    }

    return Sij2;
}

std::map<std::string, std::shared_ptr<Matrix> > USAPT0::build_Cbar(std::shared_ptr<Matrix> Sa, std::shared_ptr<Matrix> Sb)
{
    std::map<std::string, std::shared_ptr<Matrix> > Cbar;

    int nso = Cocca_A_->nrow();
    int nA = Cocca_A_->ncol();
    int nB = Cocca_B_->ncol();
    int no = nA + nB;

    double** Sp = Sa->pointer();
    double** CAp = Cocca_A_->pointer();
    double** CBp = Cocca_B_->pointer();
    double** Cp;

    Cbar["C_Ta_A"] = std::shared_ptr<Matrix>(new Matrix("C_Ta_A", nso, nA));
    Cp = Cbar["C_Ta_A"]->pointer();
    C_DGEMM('N','N',nso,nA,nA,1.0,CAp[0],nA,&Sp[0][0],no,0.0,Cp[0],nA);

    Cbar["C_Ta_B"] = std::shared_ptr<Matrix>(new Matrix("C_Ta_B", nso, nB));
    Cp = Cbar["C_Ta_B"]->pointer();
    C_DGEMM('N','N',nso,nB,nB,1.0,CBp[0],nB,&Sp[nA][nA],no,0.0,Cp[0],nB);

    Cbar["C_Ta_BA"] = std::shared_ptr<Matrix>(new Matrix("C_Ta_BA", nso, nA));
    Cp = Cbar["C_Ta_BA"]->pointer();
    C_DGEMM('N','N',nso,nA,nB,1.0,CBp[0],nB,&Sp[nA][0],no,0.0,Cp[0],nA);

//    Now we switch to the beta quantities
    nso = Coccb_A_->nrow();
    nA = Coccb_A_->ncol();
    nB = Coccb_B_->ncol();
    no = nA + nB;

    Sp = Sb->pointer();
    CAp = Coccb_A_->pointer();
    CBp = Coccb_B_->pointer();

    Cbar["C_Tb_A"] = std::shared_ptr<Matrix>(new Matrix("C_Tb_A", nso, nA));
    Cp = Cbar["C_Tb_A"]->pointer();
    C_DGEMM('N','N',nso,nA,nA,1.0,CAp[0],nA,&Sp[0][0],no,0.0,Cp[0],nA);

    Cbar["C_Tb_B"] = std::shared_ptr<Matrix>(new Matrix("C_Tb_B", nso, nB));
    Cp = Cbar["C_Tb_B"]->pointer();
    C_DGEMM('N','N',nso,nB,nB,1.0,CBp[0],nB,&Sp[nA][nA],no,0.0,Cp[0],nB);

    Cbar["C_Tb_BA"] = std::shared_ptr<Matrix>(new Matrix("C_Tb_BA", nso, nA));
    Cp = Cbar["C_Tb_BA"]->pointer();
    C_DGEMM('N','N',nso,nA,nB,1.0,CBp[0],nB,&Sp[nA][0],no,0.0,Cp[0],nA);


    return Cbar;
}

std::map< std::string,  std::shared_ptr<Matrix> > USAPT0::compute_x(std::shared_ptr<JK> jk,
                                                                     std::shared_ptr<Matrix> wa_B,
                                                                     std::shared_ptr<Matrix> wb_B,
                                                                     std::shared_ptr<Matrix> wa_A,
                                                                     std::shared_ptr<Matrix> wb_A)
{
    std::shared_ptr<CPKS_USAPT0> cpks(new CPKS_USAPT0);

    // Effective constructor
    cpks->delta_ = cpks_delta_;
    cpks->maxiter_ = cpks_maxiter_;
    cpks->jk_ = jk;

    cpks->wa_A_ = wa_A; // I don't like convention reversal.
    cpks->wb_A_ = wb_A;
    cpks->Cocca_A_ = Cocca_A_;
    cpks->Coccb_A_ = Coccb_A_;
    cpks->Cvira_A_ = Cvira_A_;
    cpks->Cvirb_A_ = Cvirb_A_;
    cpks->eps_occa_A_ = eps_occa_A_;
    cpks->eps_occb_A_ = eps_occb_A_;
    cpks->eps_vira_A_ = eps_vira_A_;
    cpks->eps_virb_A_ = eps_virb_A_;

    cpks->wa_B_ = wa_B; // Still don't like convention reversal.
    cpks->wb_B_ = wb_B;
    cpks->Cocca_B_ = Cocca_B_;
    cpks->Cvira_B_ = Cvira_B_;
    cpks->eps_occa_B_ = eps_occa_B_;
    cpks->eps_vira_B_ = eps_vira_B_;
    cpks->Coccb_B_ = Coccb_B_;
    cpks->Cvirb_B_ = Cvirb_B_;
    cpks->eps_occb_B_ = eps_occb_B_;
    cpks->eps_virb_B_ = eps_virb_B_;

    // Gogo CPKS
    cpks->compute_cpks();

    // Unpack
    std::map < std::string, std::shared_ptr<Matrix> > x_sol;
    x_sol["Aa"] = cpks->xa_A_;
    x_sol["Ab"] = cpks->xb_A_;
    x_sol["Ba"] = cpks->xa_B_;
    x_sol["Bb"] = cpks->xb_B_;

    return x_sol;
}

CPKS_USAPT0::CPKS_USAPT0()
{
}
CPKS_USAPT0::~CPKS_USAPT0()
{
}
void CPKS_USAPT0::compute_cpks()
{
    // Allocate
    xa_A_ = std::shared_ptr<Matrix>(wa_B_->clone());
    xa_B_ = std::shared_ptr<Matrix>(wa_A_->clone());
    xb_A_ = std::shared_ptr<Matrix>(wb_B_->clone());
    xb_B_ = std::shared_ptr<Matrix>(wb_A_->clone());
    xa_A_->zero();
    xa_B_->zero();
    xb_A_->zero();
    xb_B_->zero();

    std::shared_ptr<Matrix> ra_A(wa_B_->clone());
    std::shared_ptr<Matrix> za_A(wa_B_->clone());
    std::shared_ptr<Matrix> pa_A(wa_B_->clone());
    std::shared_ptr<Matrix> ra_B(wa_A_->clone());
    std::shared_ptr<Matrix> za_B(wa_A_->clone());
    std::shared_ptr<Matrix> pa_B(wa_A_->clone());
    std::shared_ptr<Matrix> rb_A(wb_B_->clone());
    std::shared_ptr<Matrix> zb_A(wb_B_->clone());
    std::shared_ptr<Matrix> pb_A(wb_B_->clone());
    std::shared_ptr<Matrix> rb_B(wb_A_->clone());
    std::shared_ptr<Matrix> zb_B(wb_A_->clone());
    std::shared_ptr<Matrix> pb_B(wb_A_->clone());
    // This is using PCG for convergence acceleration

    // Initialization done above if x = 0

    preconditioner(ra_A,za_A,eps_occa_A_,eps_vira_A_);
    preconditioner(ra_B,za_B,eps_occa_B_,eps_vira_B_);
    preconditioner(rb_A,zb_A,eps_occb_A_,eps_virb_A_);
    preconditioner(rb_B,zb_B,eps_occb_B_,eps_virb_B_);

    // Uncoupled value
//    outfile->Printf( "(A<-B): %24.16E\n", - za_A->vector_dot(wa_B_) - zb_A->vector_dot(wb_B_));
//    outfile->Printf( "(B<-A): %24.16E\n", - za_B->vector_dot(wa_A_)- zb_B->vector_dot(wb_A_));

    pa_A->copy(za_A);
    pa_B->copy(za_B);
    pb_A->copy(zb_A);
    pb_B->copy(zb_B);

    double zr_old_A = za_A->vector_dot(ra_A);
    zr_old_A += zb_A->vector_dot(rb_A);
    double zr_old_B = za_B->vector_dot(ra_B);
    zr_old_B += zb_B->vector_dot(rb_B);

    double r2A = 1.0;
    double r2B = 1.0;

    double b2A = sqrt(wa_B_->vector_dot(wa_B_) + wb_B_->vector_dot(wb_B_));
    double b2B = sqrt(wa_A_->vector_dot(wa_A_) + wb_A_->vector_dot(wb_A_));

    outfile->Printf( "  ==> CPKS Iterations <==\n\n");

    outfile->Printf( "    Maxiter     = %11d\n", maxiter_);
    outfile->Printf( "    Convergence = %11.3E\n", delta_);
    outfile->Printf( "\n");

    time_t start;
    time_t stop;

    start = time(NULL);

    outfile->Printf( "    -----------------------------------------\n");
    outfile->Printf( "    %-4s %11s  %11s  %10s\n", "Iter", "Monomer A", "Monomer B", "Time [s]");
    outfile->Printf( "    -----------------------------------------\n");

    int iter;
    for (iter = 0; iter < maxiter_; iter++) {

        std::map<std::string, std::shared_ptr<Matrix> > b;
        if (r2A > delta_) {
            b["Aa"] = pa_A;
            b["Ab"] = pb_A;
        }
        if (r2B > delta_) {
            b["Ba"] = pa_B;
            b["Bb"] = pb_B;
        }

        std::map<std::string, std::shared_ptr<Matrix> > s =
            product(b);

        if (r2A > delta_) {
            std::shared_ptr<Matrix> sa_A = s["Aa"];
            std::shared_ptr<Matrix> sb_A = s["Ab"];
            double alpha = (ra_A->vector_dot(za_A) + rb_A->vector_dot(zb_A))
                           / (pa_A->vector_dot(sa_A) + pb_A->vector_dot(sb_A));
            if (alpha < 0.0) {
                throw PSIEXCEPTION("Monomer A: A Matrix is not SPD");
            }
            int no = xa_A_->nrow();
            int nv = xa_A_->ncol();
            double** xp = xa_A_->pointer();
            double** rp = ra_A->pointer();
            double** pp = pa_A->pointer();
            double** sp = sa_A->pointer();
            C_DAXPY(no*nv, alpha,pp[0],1,xp[0],1);
            C_DAXPY(no*nv,-alpha,sp[0],1,rp[0],1);
            r2A = C_DDOT(no*nv,rp[0],1,rp[0],1);
            no = xb_A_->nrow();
            nv = xb_A_->ncol();
            xp = xb_A_->pointer();
            rp = rb_A->pointer();
            pp = pb_A->pointer();
            sp = sb_A->pointer();
            C_DAXPY(no*nv, alpha,pp[0],1,xp[0],1);
            C_DAXPY(no*nv,-alpha,sp[0],1,rp[0],1);
            r2A += C_DDOT(no*nv,rp[0],1,rp[0],1);
            r2A = sqrt(r2A) / b2A;
        }

        if (r2B > delta_) {
            std::shared_ptr<Matrix> sa_B = s["Ba"];
            std::shared_ptr<Matrix> sb_B = s["Bb"];
            double alpha = (ra_B->vector_dot(za_B) + rb_B->vector_dot(zb_B))
                           / (pa_B->vector_dot(sa_B) + pb_B->vector_dot(sb_B));
            if (alpha < 0.0) {
                throw PSIEXCEPTION("Monomer B: A Matrix is not SPD");
            }
            int no = xa_B_->nrow();
            int nv = xa_B_->ncol();
            double** xp = xa_B_->pointer();
            double** rp = ra_B->pointer();
            double** pp = pa_B->pointer();
            double** sp = sa_B->pointer();
            C_DAXPY(no*nv, alpha,pp[0],1,xp[0],1);
            C_DAXPY(no*nv,-alpha,sp[0],1,rp[0],1);
            r2B = C_DDOT(no*nv,rp[0],1,rp[0],1);
            no = xb_B_->nrow();
            nv = xb_B_->ncol();
            xp = xb_B_->pointer();
            rp = rb_B->pointer();
            pp = pb_B->pointer();
            sp = sb_B->pointer();
            C_DAXPY(no*nv, alpha,pp[0],1,xp[0],1);
            C_DAXPY(no*nv,-alpha,sp[0],1,rp[0],1);
            r2B += C_DDOT(no*nv,rp[0],1,rp[0],1);
            r2B = sqrt(r2B) / b2B;
        }

        stop = time(NULL);
        outfile->Printf( "    %-4d %11.3E%1s %11.3E%1s %10ld\n", iter+1,
            r2A, (r2A < delta_ ? "*" : " "),
            r2B, (r2B < delta_ ? "*" : " "),
            stop-start
            );

        if (r2A <= delta_ && r2B <= delta_) {
            break;
        }

        if (r2A > delta_) {
            preconditioner(ra_A,za_A,eps_occa_A_,eps_vira_A_);
            preconditioner(rb_A,zb_A,eps_occb_A_,eps_virb_A_);
            double zr_new = za_A->vector_dot(ra_A);
            zr_new += zb_A->vector_dot(rb_A);
            double beta = zr_new / zr_old_A;
            zr_old_A = zr_new;

            int no = pa_A->nrow();
            int nv = pa_A->ncol();
            double** pp = pa_A->pointer();
            double** zp = za_A->pointer();
            C_DSCAL(no*nv,beta,pp[0],1);
            C_DAXPY(no*nv,1.0,zp[0],1,pp[0],1);
            no = pb_A->nrow();
            nv = pb_A->ncol();
            pp = pb_A->pointer();
            zp = zb_A->pointer();
            C_DSCAL(no*nv,beta,pp[0],1);
            C_DAXPY(no*nv,1.0,zp[0],1,pp[0],1);
        }

        if (r2B > delta_) {
            preconditioner(ra_B,za_B,eps_occa_B_,eps_vira_B_);
            preconditioner(rb_B,zb_B,eps_occb_B_,eps_virb_B_);
            double zr_new = za_B->vector_dot(ra_B);
            zr_new += zb_B->vector_dot(rb_B);
            double beta = zr_new / zr_old_B;
            zr_old_B = zr_new;

            int no = pa_B->nrow();
            int nv = pa_B->ncol();
            double** pp = pa_B->pointer();
            double** zp = za_B->pointer();
            C_DSCAL(no*nv,beta,pp[0],1);
            C_DAXPY(no*nv,1.0,zp[0],1,pp[0],1);
            no = pb_B->nrow();
            nv = pb_B->ncol();
            pp = pb_B->pointer();
            zp = zb_B->pointer();
            C_DSCAL(no*nv,beta,pp[0],1);
            C_DAXPY(no*nv,1.0,zp[0],1,pp[0],1);
        }
    }

    outfile->Printf( "    -----------------------------------------\n");
    outfile->Printf( "\n");

    if (iter == maxiter_)
        throw PSIEXCEPTION("CPKS did not converge.");
}
void CPKS_USAPT0::preconditioner(std::shared_ptr<Matrix> r,
                               std::shared_ptr<Matrix> z,
                               std::shared_ptr<Vector> o,
                               std::shared_ptr<Vector> v)
{
    int no = o->dim();
    int nv = v->dim();

    double** rp = r->pointer();
    double** zp = z->pointer();

    double* op = o->pointer();
    double* vp = v->pointer();

    for (int i = 0; i < no; i++) {
        for (int a = 0; a < nv; a++) {
            zp[i][a] = rp[i][a] / (vp[a] - op[i]);
        }
    }
}

std::map<std::string, std::shared_ptr<Matrix> > CPKS_USAPT0::product(std::map<std::string, std::shared_ptr<Matrix> >& b)
{
    std::map<std::string, std::shared_ptr<Matrix> > s;

    bool do_A = b.count("Aa") || b.count("Ab");
    bool do_B = b.count("Ba") || b.count("Bb");


    std::vector<SharedMatrix>& Cl = jk_->C_left();
    std::vector<SharedMatrix>& Cr = jk_->C_right();
    Cl.clear();
    Cr.clear();

    if (do_A) {
        Cl.push_back(Cocca_A_);
        Cl.push_back(Coccb_A_);
        int no = b["Aa"]->nrow();
        int nv = b["Aa"]->ncol();
        int nso = Cvira_A_->nrow();
        double** Cp = Cvira_A_->pointer();
        double** bp = b["Aa"]->pointer();
        std::shared_ptr<Matrix> T(new Matrix("T",nso,no));
        double** Tp = T->pointer();
        C_DGEMM('N','T',nso,no,nv,1.0,Cp[0],nv,bp[0],nv,0.0,Tp[0],no);
        Cr.push_back(T);
        no = b["Ab"]->nrow();
        nv = b["Ab"]->ncol();
        nso = Cvirb_A_->nrow();
        Cp = Cvirb_A_->pointer();
        bp = b["Ab"]->pointer();
        T = std::shared_ptr<Matrix>(new Matrix("T",nso,no));
        Tp = T->pointer();
        C_DGEMM('N','T',nso,no,nv,1.0,Cp[0],nv,bp[0],nv,0.0,Tp[0],no);
        Cr.push_back(T);
    }

    if (do_B) {
        Cl.push_back(Cocca_B_);
        Cl.push_back(Coccb_B_);
        int no = b["Ba"]->nrow();
        int nv = b["Ba"]->ncol();
        int nso = Cvira_B_->nrow();
        double** Cp = Cvira_B_->pointer();
        double** bp = b["Ba"]->pointer();
        std::shared_ptr<Matrix> T(new Matrix("T",nso,no));
        double** Tp = T->pointer();
        C_DGEMM('N','T',nso,no,nv,1.0,Cp[0],nv,bp[0],nv,0.0,Tp[0],no);
        Cr.push_back(T);
        no = b["Bb"]->nrow();
        nv = b["Bb"]->ncol();
        nso = Cvirb_B_->nrow();
        Cp = Cvirb_B_->pointer();
        bp = b["Bb"]->pointer();
        T = std::shared_ptr<Matrix>(new Matrix("T",nso,no));
        Tp = T->pointer();
        C_DGEMM('N','T',nso,no,nv,1.0,Cp[0],nv,bp[0],nv,0.0,Tp[0],no);
        Cr.push_back(T);
    }

    jk_->compute();

    const std::vector<SharedMatrix>& J = jk_->J();
    const std::vector<SharedMatrix>& K = jk_->K();

    int indA = 0;
    int indB = (do_A ? 2 : 0);

    if (do_A) {
        std::shared_ptr<Matrix> Jva = J[indA];
        std::shared_ptr<Matrix> Jvb = J[indA + 1];
        std::shared_ptr<Matrix> Kva = K[indA];
        std::shared_ptr<Matrix> Kvb = K[indA + 1];
        Jva->scale(2.0);
        Jvb->scale(2.0);
        std::shared_ptr<Matrix> T(Jva->clone());
        T->add(Jvb);
        Jva->copy(T);
        Jva->subtract(Kva);
        Jva->subtract(Kva->transpose());
        Jvb->copy(T);
        Jvb->subtract(Kvb);
        Jvb->subtract(Kvb->transpose());

        int no = b["Aa"]->nrow();
        int nv = b["Aa"]->ncol();
        int nso = Cvira_A_->nrow();
        T = std::shared_ptr<Matrix>(new Matrix("T", no, nso));
        s["Aa"] = std::shared_ptr<Matrix>(new Matrix("SAa", no, nv));
        double** Cop = Cocca_A_->pointer();
        double** Cvp = Cvira_A_->pointer();
        double** Jp = Jva->pointer();
        double** Tp = T->pointer();
        double** Sp = s["Aa"]->pointer();
        C_DGEMM('T','N',no,nso,nso,1.0,Cop[0],no,Jp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',no,nv,nso,1.0,Tp[0],nso,Cvp[0],nv,0.0,Sp[0],nv);

        double** bp = b["Aa"]->pointer();
        double* op = eps_occa_A_->pointer();
        double* vp = eps_vira_A_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
        // Reproduce the above for the beta part.

        no = b["Ab"]->nrow();
        nv = b["Ab"]->ncol();
        nso = Cvirb_A_->nrow();
        T = std::shared_ptr<Matrix>(new Matrix("T", no, nso));
        s["Ab"] = std::shared_ptr<Matrix>(new Matrix("SAb", no, nv));
        Cop = Coccb_A_->pointer();
        Cvp = Cvirb_A_->pointer();
        Jp = Jvb->pointer();
        Tp = T->pointer();
        Sp = s["Ab"]->pointer();
        C_DGEMM('T','N',no,nso,nso,1.0,Cop[0],no,Jp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',no,nv,nso,1.0,Tp[0],nso,Cvp[0],nv,0.0,Sp[0],nv);

        bp = b["Ab"]->pointer();
        op = eps_occb_A_->pointer();
        vp = eps_virb_A_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    }

    if (do_B) {
        std::shared_ptr<Matrix> Jva = J[indB];
        std::shared_ptr<Matrix> Jvb = J[indB + 1];
        std::shared_ptr<Matrix> Kva = K[indB];
        std::shared_ptr<Matrix> Kvb = K[indB + 1];
        Jva->scale(2.0);
        Jvb->scale(2.0);
        std::shared_ptr<Matrix> T(Jva->clone());
        T->add(Jvb);
        Jva->copy(T);
        Jva->subtract(Kva);
        Jva->subtract(Kva->transpose());
        Jvb->copy(T);
        Jvb->subtract(Kvb);
        Jvb->subtract(Kvb->transpose());

        int no = b["Ba"]->nrow();
        int nv = b["Ba"]->ncol();
        int nso = Cvira_B_->nrow();
        T = std::shared_ptr<Matrix>(new Matrix("T", no, nso));
        s["Ba"] = std::shared_ptr<Matrix>(new Matrix("SBa", no, nv));
        double** Cop = Cocca_B_->pointer();
        double** Cvp = Cvira_B_->pointer();
        double** Jp = Jva->pointer();
        double** Tp = T->pointer();
        double** Sp = s["Ba"]->pointer();
        C_DGEMM('T','N',no,nso,nso,1.0,Cop[0],no,Jp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',no,nv,nso,1.0,Tp[0],nso,Cvp[0],nv,0.0,Sp[0],nv);

        double** bp = b["Ba"]->pointer();
        double* op = eps_occa_B_->pointer();
        double* vp = eps_vira_B_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }

        // Reproduce the above for the beta part

        no = b["Bb"]->nrow();
        nv = b["Bb"]->ncol();
        nso = Cvirb_B_->nrow();
        T = std::shared_ptr<Matrix>(new Matrix("T", no, nso));
        s["Bb"] = std::shared_ptr<Matrix>(new Matrix("SBb", no, nv));
        Cop = Coccb_B_->pointer();
        Cvp = Cvirb_B_->pointer();
        Jp = Jvb->pointer();
        Tp = T->pointer();
        Sp = s["Bb"]->pointer();
        C_DGEMM('T','N',no,nso,nso,1.0,Cop[0],no,Jp[0],nso,0.0,Tp[0],nso);
        C_DGEMM('N','N',no,nv,nso,1.0,Tp[0],nso,Cvp[0],nv,0.0,Sp[0],nv);

        bp = b["Bb"]->pointer();
        op = eps_occb_B_->pointer();
        vp = eps_virb_B_->pointer();
        for (int i = 0; i < no; i++) {
            for (int a = 0; a < nv; a++) {
                Sp[i][a] += bp[i][a] * (vp[a] - op[i]);
            }
        }
    }

    return s;
}

void USAPT0::mp2_terms()
{
// Note: the unrestricted version implements the second-quantized formula, which was also the one
// implemented for induction. The restricted code apparently implements the density matrix
// formula for dispersion, although it uses the second-quantized one for induction.
// On looking closer at the Hesselmann paper, it may be that induction is the same in the second-quantized
// and density-matrix formalism anyway.

    outfile->Printf( "  PT2 TERMS:\n\n");

    // => Sizing <= //

    int nn = primary_->nbf();

    int naa = Caocca_A_->colspi()[0];
    int nab = Caocca_B_->colspi()[0];
    int nar = Cavira_A_->colspi()[0];
    int nas = Cavira_B_->colspi()[0];
    int nba = Caoccb_A_->colspi()[0];
    int nbb = Caoccb_B_->colspi()[0];
    int nbr = Cavirb_A_->colspi()[0];
    int nbs = Cavirb_B_->colspi()[0];
    int nQ = mp2fit_->nbf();
    size_t narQ = nar * (size_t) nQ;
    size_t nbrQ = nbr * (size_t) nQ;
    size_t nasQ = nas * (size_t) nQ;
    size_t nbsQ = nbs * (size_t) nQ;

    int nT = 1;
    #ifdef _OPENMP
        nT = Process::environment.get_n_threads();
    #endif

    // => Stashed Variables <= //

    std::shared_ptr<Matrix> S   = vars_["S"];
    std::shared_ptr<Matrix> Da_A = vars_["Da_A"];
    std::shared_ptr<Matrix> Db_A = vars_["Db_A"];
    std::shared_ptr<Matrix> El_pot_A = vars_["El_pot_A"];
    std::shared_ptr<Matrix> ha_A = vars_["ha_A"];
    std::shared_ptr<Matrix> hb_A = vars_["hb_A"];
    std::shared_ptr<Matrix> Da_B = vars_["Da_B"];
    std::shared_ptr<Matrix> Db_B = vars_["Db_B"];
    std::shared_ptr<Matrix> El_pot_B = vars_["El_pot_B"];
    std::shared_ptr<Matrix> ha_B = vars_["ha_B"];
    std::shared_ptr<Matrix> hb_B = vars_["hb_B"];
    std::shared_ptr<Matrix> Ka_O = vars_["Ka_O"];
    std::shared_ptr<Matrix> Kb_O = vars_["Kb_O"];

    // => Auxiliary C matrices <= //
//    We build them to maximize the reuse of intermediates

    std::shared_ptr<Matrix> Ca_a2 = Matrix::doublet(Da_B,S);
    std::shared_ptr<Matrix> Cb_a2 = Matrix::doublet(Db_B,S);
    std::shared_ptr<Matrix> Ca_b2 = Matrix::doublet(Da_A,S);
    std::shared_ptr<Matrix> Cb_b2 = Matrix::doublet(Db_A,S);

    std::shared_ptr<Matrix> Ca_s1 = Matrix::doublet(Ca_b2,Cavira_B_);
    std::shared_ptr<Matrix> Cb_s1 = Matrix::doublet(Cb_b2,Cavirb_B_);
    std::shared_ptr<Matrix> Ca_r1 = Matrix::doublet(Ca_a2,Cavira_A_);
    std::shared_ptr<Matrix> Cb_r1 = Matrix::doublet(Cb_a2,Cavirb_A_);

    std::shared_ptr<Matrix> Ca_s3 = Matrix::triplet(Da_B,S,Ca_s1);
    std::shared_ptr<Matrix> Cb_s3 = Matrix::triplet(Db_B,S,Cb_s1);
    std::shared_ptr<Matrix> Ca_r3 = Matrix::triplet(Da_A,S,Ca_r1);
    std::shared_ptr<Matrix> Cb_r3 = Matrix::triplet(Db_A,S,Cb_r1);

    Ca_s3->subtract(Ca_s1);
    Cb_s3->subtract(Cb_s1);
    Ca_r3->subtract(Ca_r1);
    Cb_r3->subtract(Cb_r1);

    Ca_s1->subtract(Cavira_B_);
    Cb_s1->subtract(Cavirb_B_);
    Ca_r1->subtract(Cavira_A_);
    Cb_r1->subtract(Cavirb_A_);

    Ca_r1->scale(-1.0);
    Cb_r1->scale(-1.0);

    Ca_a2 = Matrix::doublet(Ca_a2,Caocca_A_);
    Cb_a2 = Matrix::doublet(Cb_a2,Caoccb_A_);
    Ca_b2 = Matrix::doublet(Ca_b2,Caocca_B_);
    Cb_b2 = Matrix::doublet(Cb_b2,Caoccb_B_);

    std::shared_ptr<Matrix> Ca_a4 = Matrix::triplet(Da_A,S,Ca_a2);
    std::shared_ptr<Matrix> Cb_a4 = Matrix::triplet(Db_A,S,Cb_a2);
    std::shared_ptr<Matrix> Ca_b4 = Matrix::triplet(Da_B,S,Ca_b2);
    std::shared_ptr<Matrix> Cb_b4 = Matrix::triplet(Db_B,S,Cb_b2);

    // => Auxiliary Fock-derived matrices <= //
    // We build them to maximize intermediate reuse and minimize the number
    // of matrix multiplications

    std::shared_ptr<Matrix> Da_BS = Matrix::doublet(Da_B,S);
    std::shared_ptr<Matrix> Db_BS = Matrix::doublet(Db_B,S);
    std::shared_ptr<Matrix> Da_AS = Matrix::doublet(Da_A,S);
    std::shared_ptr<Matrix> Db_AS = Matrix::doublet(Db_A,S);

    std::shared_ptr<Matrix> Ta_as = Matrix::doublet(El_pot_B,Da_AS);
    std::shared_ptr<Matrix> Tb_as = Matrix::doublet(El_pot_B,Db_AS);
    std::shared_ptr<Matrix> Ta_br = Matrix::doublet(El_pot_A,Da_BS);
    std::shared_ptr<Matrix> Tb_br = Matrix::doublet(El_pot_A,Db_BS);

    std::shared_ptr<Matrix> Sa_Bar = Matrix::triplet(Caocca_A_,S,Matrix::doublet(Da_BS,Cavira_A_),true,false,false);
    std::shared_ptr<Matrix> Sb_Bar = Matrix::triplet(Caoccb_A_,S,Matrix::doublet(Db_BS,Cavirb_A_),true,false,false);
    std::shared_ptr<Matrix> Sa_Abs = Matrix::triplet(Caocca_B_,S,Matrix::doublet(Da_AS,Cavira_B_),true,false,false);
    std::shared_ptr<Matrix> Sb_Abs = Matrix::triplet(Caoccb_B_,S,Matrix::doublet(Db_AS,Cavirb_B_),true,false,false);

    Da_BS.reset();
    Db_BS.reset();
    Da_AS.reset();
    Db_AS.reset();

//  Build the other auxiliary matrices in the AO basis

    Ta_as->add(Matrix::triplet(S,Da_B,El_pot_A));
    Tb_as->add(Matrix::triplet(S,Db_B,El_pot_A));
    Ta_br->add(Matrix::triplet(S,Da_A,El_pot_B));
    Tb_br->add(Matrix::triplet(S,Db_A,El_pot_B));

    Ta_as->subtract(ha_B);
    Tb_as->subtract(hb_B);
    Ta_br->subtract(ha_A);
    Tb_br->subtract(hb_A);

    Ta_as->subtract(Ka_O);
    Tb_as->subtract(Kb_O);
    Ta_br->subtract(Ka_O->transpose());
    Tb_br->subtract(Kb_O->transpose());


    std::shared_ptr<Matrix> Sa_as = Matrix::triplet(Caocca_A_,S,Cavira_B_,true,false,false);
    std::shared_ptr<Matrix> Sb_as = Matrix::triplet(Caoccb_A_,S,Cavirb_B_,true,false,false);
    std::shared_ptr<Matrix> Sa_br = Matrix::triplet(Caocca_B_,S,Cavira_A_,true,false,false);
    std::shared_ptr<Matrix> Sb_br = Matrix::triplet(Caoccb_B_,S,Cavirb_A_,true,false,false);

    std::shared_ptr<Matrix> Qa_as = Matrix::triplet(Caocca_A_,Ta_as,Cavira_B_,true,false,false);
    std::shared_ptr<Matrix> Qb_as = Matrix::triplet(Caoccb_A_,Tb_as,Cavirb_B_,true,false,false);
    std::shared_ptr<Matrix> Qa_br = Matrix::triplet(Caocca_B_,Ta_br,Cavira_A_,true,false,false);
    std::shared_ptr<Matrix> Qb_br = Matrix::triplet(Caoccb_B_,Tb_br,Cavirb_A_,true,false,false);

    std::shared_ptr<Matrix> Qa_ar = Matrix::triplet(Caocca_A_,El_pot_B,Cavira_A_,true,false,false);
    std::shared_ptr<Matrix> Qb_ar = Matrix::triplet(Caoccb_A_,El_pot_B,Cavirb_A_,true,false,false);
    std::shared_ptr<Matrix> Qa_bs = Matrix::triplet(Caocca_B_,El_pot_A,Cavira_B_,true,false,false);
    std::shared_ptr<Matrix> Qb_bs = Matrix::triplet(Caoccb_B_,El_pot_A,Cavirb_B_,true,false,false);

    Ta_as.reset();
    Tb_as.reset();
    Ta_br.reset();
    Tb_br.reset();

    S.reset();
    Da_A.reset();
    Db_A.reset();
    El_pot_A.reset();
    ha_A.reset();
    hb_A.reset();
    Da_B.reset();
    Db_B.reset();
    El_pot_B.reset();
    ha_B.reset();
    hb_B.reset();
    Ka_O.reset();
    Kb_O.reset();

    vars_.clear();

    // => Memory <= //

    // => Integrals from the THCE <= //

    std::shared_ptr<DFERI> df = DFERI::build(primary_,mp2fit_,Process::environment.options);
    df->clear();

    std::vector<std::shared_ptr<Matrix> > Cs;
    Cs.push_back(Caocca_A_);
    Cs.push_back(Cavira_A_);
    Cs.push_back(Caocca_B_);
    Cs.push_back(Cavira_B_);
    Cs.push_back(Ca_r1);
    Cs.push_back(Ca_s1);
    Cs.push_back(Ca_a2);
    Cs.push_back(Ca_b2);
    Cs.push_back(Ca_r3);
    Cs.push_back(Ca_s3);
    Cs.push_back(Ca_a4);
    Cs.push_back(Ca_b4);

    Cs.push_back(Caoccb_A_);
    Cs.push_back(Cavirb_A_);
    Cs.push_back(Caoccb_B_);
    Cs.push_back(Cavirb_B_);
    Cs.push_back(Cb_r1);
    Cs.push_back(Cb_s1);
    Cs.push_back(Cb_a2);
    Cs.push_back(Cb_b2);
    Cs.push_back(Cb_r3);
    Cs.push_back(Cb_s3);
    Cs.push_back(Cb_a4);
    Cs.push_back(Cb_b4);
    std::shared_ptr<Matrix> Call = Matrix::horzcat(Cs);
    Cs.clear();

    df->set_C(Call);
    df->set_memory(memory_ - Call->nrow() * Call->ncol());

    int offset = 0;
    df->add_space("a_a",offset,offset+Caocca_A_->colspi()[0]); offset += Caocca_A_->colspi()[0];
    df->add_space("a_r",offset,offset+Cavira_A_->colspi()[0]); offset += Cavira_A_->colspi()[0];
    df->add_space("a_b",offset,offset+Caocca_B_->colspi()[0]); offset += Caocca_B_->colspi()[0];
    df->add_space("a_s",offset,offset+Cavira_B_->colspi()[0]); offset += Cavira_B_->colspi()[0];
    df->add_space("a_r1",offset,offset+Ca_r1->colspi()[0]); offset += Ca_r1->colspi()[0];
    df->add_space("a_s1",offset,offset+Ca_s1->colspi()[0]); offset += Ca_s1->colspi()[0];
    df->add_space("a_a2",offset,offset+Ca_a2->colspi()[0]); offset += Ca_a2->colspi()[0];
    df->add_space("a_b2",offset,offset+Ca_b2->colspi()[0]); offset += Ca_b2->colspi()[0];
    df->add_space("a_r3",offset,offset+Ca_r3->colspi()[0]); offset += Ca_r3->colspi()[0];
    df->add_space("a_s3",offset,offset+Ca_s3->colspi()[0]); offset += Ca_s3->colspi()[0];
    df->add_space("a_a4",offset,offset+Ca_a4->colspi()[0]); offset += Ca_a4->colspi()[0];
    df->add_space("a_b4",offset,offset+Ca_b4->colspi()[0]); offset += Ca_b4->colspi()[0];

    df->add_space("b_a",offset,offset+Caoccb_A_->colspi()[0]); offset += Caoccb_A_->colspi()[0];
    df->add_space("b_r",offset,offset+Cavirb_A_->colspi()[0]); offset += Cavirb_A_->colspi()[0];
    df->add_space("b_b",offset,offset+Caoccb_B_->colspi()[0]); offset += Caoccb_B_->colspi()[0];
    df->add_space("b_s",offset,offset+Cavirb_B_->colspi()[0]); offset += Cavirb_B_->colspi()[0];
    df->add_space("b_r1",offset,offset+Cb_r1->colspi()[0]); offset += Cb_r1->colspi()[0];
    df->add_space("b_s1",offset,offset+Cb_s1->colspi()[0]); offset += Cb_s1->colspi()[0];
    df->add_space("b_a2",offset,offset+Cb_a2->colspi()[0]); offset += Cb_a2->colspi()[0];
    df->add_space("b_b2",offset,offset+Cb_b2->colspi()[0]); offset += Cb_b2->colspi()[0];
    df->add_space("b_r3",offset,offset+Cb_r3->colspi()[0]); offset += Cb_r3->colspi()[0];
    df->add_space("b_s3",offset,offset+Cb_s3->colspi()[0]); offset += Cb_s3->colspi()[0];
    df->add_space("b_a4",offset,offset+Cb_a4->colspi()[0]); offset += Cb_a4->colspi()[0];
    df->add_space("b_b4",offset,offset+Cb_b4->colspi()[0]); offset += Cb_b4->colspi()[0];


    df->add_pair_space("Aa_ar", "a_a", "a_r");
    df->add_pair_space("Aa_bs", "a_b", "a_s");
    df->add_pair_space("Ba_as", "a_a", "a_s1");
    df->add_pair_space("Ba_br", "a_b", "a_r1");
    df->add_pair_space("Ca_as", "a_a2","a_s");
    df->add_pair_space("Ca_br", "a_b2","a_r");
    df->add_pair_space("Da_ar", "a_a", "a_r3");
    df->add_pair_space("Da_bs", "a_b", "a_s3");
    df->add_pair_space("Ea_ar", "a_a4","a_r");
    df->add_pair_space("Ea_bs", "a_b4","a_s");

    df->add_pair_space("Ab_ar", "b_a", "b_r");
    df->add_pair_space("Ab_bs", "b_b", "b_s");
    df->add_pair_space("Bb_as", "b_a", "b_s1");
    df->add_pair_space("Bb_br", "b_b", "b_r1");
    df->add_pair_space("Cb_as", "b_a2","b_s");
    df->add_pair_space("Cb_br", "b_b2","b_r");
    df->add_pair_space("Db_ar", "b_a", "b_r3");
    df->add_pair_space("Db_bs", "b_b", "b_s3");
    df->add_pair_space("Eb_ar", "b_a4","b_r");
    df->add_pair_space("Eb_bs", "b_b4","b_s");

    Ca_r1.reset();
    Ca_s1.reset();
    Ca_a2.reset();
    Ca_b2.reset();
    Ca_r3.reset();
    Ca_s3.reset();
    Ca_a4.reset();
    Ca_b4.reset();
    Cb_r1.reset();
    Cb_s1.reset();
    Cb_a2.reset();
    Cb_b2.reset();
    Cb_r3.reset();
    Cb_s3.reset();
    Cb_a4.reset();
    Cb_b4.reset();
    Call.reset();

    df->print_header();
    df->compute();

    std::map<std::string, std::shared_ptr<Tensor> >& ints = df->ints();

    std::shared_ptr<Tensor> Aa_arT = ints["Aa_ar"];
    std::shared_ptr<Tensor> Aa_bsT = ints["Aa_bs"];
    std::shared_ptr<Tensor> Ba_asT = ints["Ba_as"];
    std::shared_ptr<Tensor> Ba_brT = ints["Ba_br"];
    std::shared_ptr<Tensor> Ca_asT = ints["Ca_as"];
    std::shared_ptr<Tensor> Ca_brT = ints["Ca_br"];
    std::shared_ptr<Tensor> Da_arT = ints["Da_ar"];
    std::shared_ptr<Tensor> Da_bsT = ints["Da_bs"];
    std::shared_ptr<Tensor> Ea_arT = ints["Ea_ar"];
    std::shared_ptr<Tensor> Ea_bsT = ints["Ea_bs"];

    std::shared_ptr<Tensor> Ab_arT = ints["Ab_ar"];
    std::shared_ptr<Tensor> Ab_bsT = ints["Ab_bs"];
    std::shared_ptr<Tensor> Bb_asT = ints["Bb_as"];
    std::shared_ptr<Tensor> Bb_brT = ints["Bb_br"];
    std::shared_ptr<Tensor> Cb_asT = ints["Cb_as"];
    std::shared_ptr<Tensor> Cb_brT = ints["Cb_br"];
    std::shared_ptr<Tensor> Db_arT = ints["Db_ar"];
    std::shared_ptr<Tensor> Db_bsT = ints["Db_bs"];
    std::shared_ptr<Tensor> Eb_arT = ints["Eb_ar"];
    std::shared_ptr<Tensor> Eb_bsT = ints["Eb_bs"];

    df.reset();

    // => Blocking <= //

    long int overhead = 0L;
    overhead += 2L * nT * nar * nas + 2L * nT * nbr * nbs;
    overhead += 2L * nT * nbr * nas + 2L * nT * nar * nbs;
    overhead += 2L * naa * nas + 2L * nab * nar + 2L * naa * nar + 2L * nab * nas;
    overhead += 2L * nba * nbs + 2L * nbb * nbr + 2L * nba * nbr + 2L * nbb * nbs;
    long int rem = memory_ - overhead;

    if (rem < 0L) {
        throw PSIEXCEPTION("Too little static memory for USAPT0::mp2_terms");
    }

    long int cost_a = 2L * nar * nQ + 2L * nas * nQ;
    cost_a += 2L * nbr * nQ + 2L * nbs * nQ;
    long int maxa_a = rem / (2L * cost_a);
    long int maxb_a = maxa_a;
    long int maxa_b = maxa_a;
    long int maxb_b = maxb_a;
    maxa_a = (maxa_a > naa ? naa : maxa_a);
    maxa_b = (maxa_b > nab ? nab : maxa_b);
    maxb_a = (maxb_a > nba ? nba : maxb_a);
    maxb_b = (maxb_b > nbb ? nbb : maxb_b);
    if (maxa_a < 1L || maxa_b < 1L || maxb_a < 1L || maxb_b < 1L) {
        throw PSIEXCEPTION("Too little dynamic memory for USAPT0::mp2_terms");
    }

    // => Tensor Slices <= //

    std::shared_ptr<Matrix> Aa_ar(new Matrix("Aa_ar",maxa_a*nar,nQ));
    std::shared_ptr<Matrix> Aa_bs(new Matrix("Aa_bs",maxa_b*nas,nQ));
    std::shared_ptr<Matrix> Ba_as(new Matrix("Ba_as",maxa_a*nas,nQ));
    std::shared_ptr<Matrix> Ba_br(new Matrix("Ba_br",maxa_b*nar,nQ));
    std::shared_ptr<Matrix> Ca_as(new Matrix("Ca_as",maxa_a*nas,nQ));
    std::shared_ptr<Matrix> Ca_br(new Matrix("Ca_br",maxa_b*nar,nQ));
    std::shared_ptr<Matrix> Da_ar(new Matrix("Da_ar",maxa_a*nar,nQ));
    std::shared_ptr<Matrix> Da_bs(new Matrix("Da_bs",maxa_b*nas,nQ));

    std::shared_ptr<Matrix> Ab_ar(new Matrix("Ab_ar",maxb_a*nbr,nQ));
    std::shared_ptr<Matrix> Ab_bs(new Matrix("Ab_bs",maxb_b*nbs,nQ));
    std::shared_ptr<Matrix> Bb_as(new Matrix("Bb_as",maxb_a*nbs,nQ));
    std::shared_ptr<Matrix> Bb_br(new Matrix("Bb_br",maxb_b*nbr,nQ));
    std::shared_ptr<Matrix> Cb_as(new Matrix("Cb_as",maxb_a*nbs,nQ));
    std::shared_ptr<Matrix> Cb_br(new Matrix("Cb_br",maxb_b*nbr,nQ));
    std::shared_ptr<Matrix> Db_ar(new Matrix("Db_ar",maxb_a*nbr,nQ));
    std::shared_ptr<Matrix> Db_bs(new Matrix("Db_bs",maxb_b*nbs,nQ));

    // => Thread Work Arrays <= //

    std::vector<std::shared_ptr<Matrix> > Taa_rs;
    std::vector<std::shared_ptr<Matrix> > Vaa_rs;
    std::vector<std::shared_ptr<Matrix> > Tbb_rs;
    std::vector<std::shared_ptr<Matrix> > Vbb_rs;
    std::vector<std::shared_ptr<Matrix> > Tba_rs;
    std::vector<std::shared_ptr<Matrix> > Vba_rs;
    std::vector<std::shared_ptr<Matrix> > Tab_rs;
    std::vector<std::shared_ptr<Matrix> > Vab_rs;
    for (int t = 0; t < nT; t++) {
        Taa_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Taa_rs",nar,nas)));
        Vaa_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Vaa_rs",nar,nas)));
        Tbb_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Tbb_rs",nbr,nbs)));
        Vbb_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Vbb_rs",nbr,nbs)));
        Tab_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Tab_rs",nar,nbs)));
        Vab_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Vab_rs",nar,nbs)));
        Tba_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Tba_rs",nbr,nas)));
        Vba_rs.push_back(std::shared_ptr<Matrix>(new Matrix("Vba_rs",nbr,nas)));
    }

    // => Pointers <= //

    double** Aa_arp = Aa_ar->pointer();
    double** Aa_bsp = Aa_bs->pointer();
    double** Ba_asp = Ba_as->pointer();
    double** Ba_brp = Ba_br->pointer();
    double** Ca_asp = Ca_as->pointer();
    double** Ca_brp = Ca_br->pointer();
    double** Da_arp = Da_ar->pointer();
    double** Da_bsp = Da_bs->pointer();

    double** Ab_arp = Ab_ar->pointer();
    double** Ab_bsp = Ab_bs->pointer();
    double** Bb_asp = Bb_as->pointer();
    double** Bb_brp = Bb_br->pointer();
    double** Cb_asp = Cb_as->pointer();
    double** Cb_brp = Cb_br->pointer();
    double** Db_arp = Db_ar->pointer();
    double** Db_bsp = Db_bs->pointer();

    double** Sa_asp = Sa_as->pointer();
    double** Sa_brp = Sa_br->pointer();
    double** Sb_asp = Sb_as->pointer();
    double** Sb_brp = Sb_br->pointer();
    double** Sa_Barp = Sa_Bar->pointer();
    double** Sa_Absp = Sa_Abs->pointer();
    double** Sb_Barp = Sb_Bar->pointer();
    double** Sb_Absp = Sb_Abs->pointer();

    double** Qa_asp = Qa_as->pointer();
    double** Qa_brp = Qa_br->pointer();
    double** Qa_arp = Qa_ar->pointer();
    double** Qa_bsp = Qa_bs->pointer();

    double** Qb_asp = Qb_as->pointer();
    double** Qb_brp = Qb_br->pointer();
    double** Qb_arp = Qb_ar->pointer();
    double** Qb_bsp = Qb_bs->pointer();

    double*  ea_ap  = eps_aocca_A_->pointer();
    double*  ea_bp  = eps_aocca_B_->pointer();
    double*  ea_rp  = eps_avira_A_->pointer();
    double*  ea_sp  = eps_avira_B_->pointer();

    double*  eb_ap  = eps_aoccb_A_->pointer();
    double*  eb_bp  = eps_aoccb_B_->pointer();
    double*  eb_rp  = eps_avirb_A_->pointer();
    double*  eb_sp  = eps_avirb_B_->pointer();

    // => File Pointers <= //

    FILE* Aa_arf = Aa_arT->file_pointer();
    FILE* Aa_bsf = Aa_bsT->file_pointer();
    FILE* Ba_asf = Ba_asT->file_pointer();
    FILE* Ba_brf = Ba_brT->file_pointer();
    FILE* Ca_asf = Ca_asT->file_pointer();
    FILE* Ca_brf = Ca_brT->file_pointer();
    FILE* Da_arf = Da_arT->file_pointer();
    FILE* Da_bsf = Da_bsT->file_pointer();
    FILE* Ea_arf = Ea_arT->file_pointer();
    FILE* Ea_bsf = Ea_bsT->file_pointer();

    FILE* Ab_arf = Ab_arT->file_pointer();
    FILE* Ab_bsf = Ab_bsT->file_pointer();
    FILE* Bb_asf = Bb_asT->file_pointer();
    FILE* Bb_brf = Bb_brT->file_pointer();
    FILE* Cb_asf = Cb_asT->file_pointer();
    FILE* Cb_brf = Cb_brT->file_pointer();
    FILE* Db_arf = Db_arT->file_pointer();
    FILE* Db_bsf = Db_bsT->file_pointer();
    FILE* Eb_arf = Eb_arT->file_pointer();
    FILE* Eb_bsf = Eb_bsT->file_pointer();

    // => Slice D + E -> D <= //

//  Here we could read more slices in memory, or just move that to the master loop directly.
    fseek(Da_arf,0L,SEEK_SET);
    fseek(Ea_arf,0L,SEEK_SET);
    for (int astart = 0; astart < naa; astart += 2L * maxa_a) {
        int nablock = (astart + maxa_a >= naa ? naa - astart : maxa_a);
        fread(Da_arp[0],sizeof(double),nablock*narQ,Da_arf);
        fread(Aa_arp[0],sizeof(double),nablock*narQ,Ea_arf);
        C_DAXPY(nablock*narQ,1.0,Aa_arp[0],1,Da_arp[0],1);
        fseek(Da_arf,sizeof(double)*astart*narQ,SEEK_SET);
        fwrite(Da_arp[0],sizeof(double),nablock*narQ,Da_arf);
    }

    fseek(Db_arf,0L,SEEK_SET);
    fseek(Eb_arf,0L,SEEK_SET);
    for (int astart = 0; astart < nba; astart += 2L * maxb_a) {
        int nablock = (astart + maxb_a >= nba ? nba - astart : maxb_a);
        fread(Db_arp[0],sizeof(double),nablock*nbrQ,Db_arf);
        fread(Ab_arp[0],sizeof(double),nablock*nbrQ,Eb_arf);
        C_DAXPY(nablock*nbrQ,1.0,Ab_arp[0],1,Db_arp[0],1);
        fseek(Db_arf,sizeof(double)*astart*nbrQ,SEEK_SET);
        fwrite(Db_arp[0],sizeof(double),nablock*nbrQ,Db_arf);
    }

    fseek(Da_bsf,0L,SEEK_SET);
    fseek(Ea_bsf,0L,SEEK_SET);
    for (int bstart = 0; bstart < nab; bstart += 2L * maxa_b) {
        int nbblock = (bstart + maxa_b >= nab ? nab - bstart : maxa_b);
        fread(Da_bsp[0],sizeof(double),nbblock*nasQ,Da_bsf);
        fread(Aa_bsp[0],sizeof(double),nbblock*nasQ,Ea_bsf);
        C_DAXPY(nbblock*nasQ,1.0,Aa_bsp[0],1,Da_bsp[0],1);
        fseek(Da_bsf,sizeof(double)*bstart*nasQ,SEEK_SET);
        fwrite(Da_bsp[0],sizeof(double),nbblock*nasQ,Da_bsf);
    }

    fseek(Db_bsf,0L,SEEK_SET);
    fseek(Eb_bsf,0L,SEEK_SET);
    for (int bstart = 0; bstart < nbb; bstart += 2L * maxb_b) {
        int nbblock = (bstart + maxb_b >= nbb ? nbb - bstart : maxb_b);
        fread(Db_bsp[0],sizeof(double),nbblock*nbsQ,Db_bsf);
        fread(Ab_bsp[0],sizeof(double),nbblock*nbsQ,Eb_bsf);
        C_DAXPY(nbblock*nbsQ,1.0,Ab_bsp[0],1,Db_bsp[0],1);
        fseek(Db_bsf,sizeof(double)*bstart*nbsQ,SEEK_SET);
        fwrite(Db_bsp[0],sizeof(double),nbblock*nbsQ,Db_bsf);
    }

    // => Targets <= //

    double Disp20 = 0.0;
    double ExchDisp20 = 0.0;

    // ==> Master Loop <== //

    fseek(Aa_arf,0L,SEEK_SET);
    fseek(Ba_asf,0L,SEEK_SET);
    fseek(Ca_asf,0L,SEEK_SET);
    fseek(Da_arf,0L,SEEK_SET);
    fseek(Ab_arf,0L,SEEK_SET);
    fseek(Bb_asf,0L,SEEK_SET);
    fseek(Cb_asf,0L,SEEK_SET);
    fseek(Db_arf,0L,SEEK_SET);
    for (int astart = 0; astart < max(naa, nba); astart += maxa_a) {
        int na_ablock = (astart + maxa_a >= naa ? naa - astart : maxa_a);
        int nb_ablock = (astart + maxb_a >= nba ? nba - astart : maxb_a);

        if ( na_ablock > 0) {
            fread(Aa_arp[0],sizeof(double),na_ablock*narQ,Aa_arf);
            fread(Ba_asp[0],sizeof(double),na_ablock*nasQ,Ba_asf);
            fread(Ca_asp[0],sizeof(double),na_ablock*nasQ,Ca_asf);
            fread(Da_arp[0],sizeof(double),na_ablock*narQ,Da_arf);
        }

        if ( nb_ablock > 0) {
            fread(Ab_arp[0],sizeof(double),nb_ablock*nbrQ,Ab_arf);
            fread(Bb_asp[0],sizeof(double),nb_ablock*nbsQ,Bb_asf);
            fread(Cb_asp[0],sizeof(double),nb_ablock*nbsQ,Cb_asf);
            fread(Db_arp[0],sizeof(double),nb_ablock*nbrQ,Db_arf);
        }


        fseek(Aa_bsf,0L,SEEK_SET);
        fseek(Ba_brf,0L,SEEK_SET);
        fseek(Ca_brf,0L,SEEK_SET);
        fseek(Da_bsf,0L,SEEK_SET);
        fseek(Ab_bsf,0L,SEEK_SET);
        fseek(Bb_brf,0L,SEEK_SET);
        fseek(Cb_brf,0L,SEEK_SET);
        fseek(Db_bsf,0L,SEEK_SET);
        for (int bstart = 0; bstart < max(nab, nbb); bstart += maxa_b) {
            int na_bblock = (bstart + maxa_b >= nab ? nab - bstart : maxa_b);
            int nb_bblock = (bstart + maxb_b >= nbb ? nbb - bstart : maxb_b);

            if ( na_bblock > 0) {
                fread(Aa_bsp[0],sizeof(double),na_bblock*nasQ,Aa_bsf);
                fread(Ba_brp[0],sizeof(double),na_bblock*narQ,Ba_brf);
                fread(Ca_brp[0],sizeof(double),na_bblock*narQ,Ca_brf);
                fread(Da_bsp[0],sizeof(double),na_bblock*nasQ,Da_bsf);
            }

            if ( nb_bblock > 0) {
                fread(Ab_bsp[0],sizeof(double),nb_bblock*nbsQ,Ab_bsf);
                fread(Bb_brp[0],sizeof(double),nb_bblock*nbrQ,Bb_brf);
                fread(Cb_brp[0],sizeof(double),nb_bblock*nbrQ,Cb_brf);
                fread(Db_bsp[0],sizeof(double),nb_bblock*nbsQ,Db_bsf);
            }

            long int nab = (na_ablock + nb_ablock) * (na_bblock + nb_bblock);

            #pragma omp parallel for schedule(dynamic) reduction(+: Disp20, ExchDisp20)
            for (long int ab = 0L; ab < nab; ab++) {
                int a = ab / (na_bblock + nb_bblock);
                int b = ab % (na_bblock + nb_bblock);

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

//                We write one condition handling the different cases here

                double** Trsp;
                double** Vrsp;
                double** Aarp;
                double** Absp;
                double** Bbrp = NULL;
                double** Basp = NULL;
                double** Casp = NULL;
                double** Cbrp = NULL;
                double** Darp;
                double** Dbsp;
                double** Qbrp = NULL;
                double** Qasp = NULL;
                double** Sbrp = NULL;
                double** Sasp = NULL;
                double** Qarp;
                double** Qbsp;
                double** SAbsp;
                double** SBarp;
                double*  eap;
                double*  ebp;
                double*  erp;
                double*  esp;
                int nr;
                int ns;

                if ( a >= na_ablock && b >= na_bblock ) {
                    a = a - na_ablock;
                    b = b - na_bblock;
                    Trsp = Tbb_rs[thread]->pointer();
                    Vrsp = Vbb_rs[thread]->pointer();
                    Aarp = Ab_arp;
                    Absp = Ab_bsp;
                    Bbrp = Bb_brp;
                    Basp = Bb_asp;
                    Cbrp = Cb_brp;
                    Casp = Cb_asp;
                    Darp = Db_arp;
                    Dbsp = Db_bsp;
                    Qbrp = Qb_brp;
                    Qasp = Qb_asp;
                    Sbrp = Sb_brp;
                    Sasp = Sb_asp;
                    Qarp = Qb_arp;
                    Qbsp = Qb_bsp;
                    SAbsp = Sb_Absp;
                    SBarp = Sb_Barp;
                    eap = eb_ap;
                    ebp = eb_bp;
                    erp = eb_rp;
                    esp = eb_sp;
                    nr = nbr;
                    ns = nbs;
                } else if ( a >= na_ablock && b < na_bblock) {
                    a = a - na_ablock;
                    Trsp = Tba_rs[thread]->pointer();
                    Vrsp = Vba_rs[thread]->pointer();
                    Aarp = Ab_arp;
                    Absp = Aa_bsp;
                    Darp = Db_arp;
                    Dbsp = Da_bsp;
                    Qarp = Qb_arp;
                    Qbsp = Qa_bsp;
                    SAbsp = Sa_Absp;
                    SBarp = Sb_Barp;
                    eap = eb_ap;
                    ebp = ea_bp;
                    erp = eb_rp;
                    esp = ea_sp;
                    nr = nbr;
                    ns = nas;
                } else if ( a < na_ablock && b >= na_bblock ) {
                    b = b - na_bblock;
                    Trsp = Tab_rs[thread]->pointer();
                    Vrsp = Vab_rs[thread]->pointer();
                    Aarp = Aa_arp;
                    Absp = Ab_bsp;
                    Darp = Da_arp;
                    Dbsp = Db_bsp;
                    Qarp = Qa_arp;
                    Qbsp = Qb_bsp;
                    SAbsp = Sb_Absp;
                    SBarp = Sa_Barp;
                    eap = ea_ap;
                    ebp = eb_bp;
                    erp = ea_rp;
                    esp = eb_sp;
                    nr = nar;
                    ns = nbs;
                } else {
                    Trsp = Taa_rs[thread]->pointer();
                    Vrsp = Vaa_rs[thread]->pointer();
                    Aarp = Aa_arp;
                    Absp = Aa_bsp;
                    Bbrp = Ba_brp;
                    Basp = Ba_asp;
                    Cbrp = Ca_brp;
                    Casp = Ca_asp;
                    Darp = Da_arp;
                    Dbsp = Da_bsp;
                    Qbrp = Qa_brp;
                    Qasp = Qa_asp;
                    Sbrp = Sa_brp;
                    Sasp = Sa_asp;
                    Qarp = Qa_arp;
                    Qbsp = Qa_bsp;
                    SAbsp = Sa_Absp;
                    SBarp = Sa_Barp;
                    eap = ea_ap;
                    ebp = ea_bp;
                    erp = ea_rp;
                    esp = ea_sp;
                    nr = nar;
                    ns = nas;
                }


                // => Amplitudes, Disp20 <= //

                C_DGEMM('N','T',nr,ns,nQ,1.0,Aarp[(a)*nr],nQ,Absp[(b)*ns],nQ,0.0,Vrsp[0],ns);

                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        Trsp[r][s] = Vrsp[r][s] / (eap[a + astart] + ebp[b + bstart] - erp[r] - esp[s]);
                        Disp20 += Trsp[r][s] * Vrsp[r][s];
                    }
                }

                // => Exch-Disp20 <= //

                // > Q1-Q3 < //

                C_DGEMM('N','T',nr,ns,nQ,1.0,Aarp[(a)*nr],nQ,Dbsp[(b)*ns],nQ,0.0,Vrsp[0],ns);
                C_DGEMM('N','T',nr,ns,nQ,1.0,Darp[(a)*nr],nQ,Absp[(b)*ns],nQ,1.0,Vrsp[0],ns);
                if ( Bbrp != NULL ) {
                   C_DGEMM('N','T',nr,ns,nQ,1.0,Bbrp[(b)*nr],nQ,Basp[(a)*ns],nQ,1.0,Vrsp[0],ns);
                   C_DGEMM('N','T',nr,ns,nQ,-1.0,Cbrp[(b)*nr],nQ,Casp[(a)*ns],nQ,1.0,Vrsp[0],ns);

                   // > V,J,K < //

                   C_DGER(nr,ns,1.0,Qbrp[b + bstart],1,Sasp[a + astart],1,Vrsp[0],ns);
                   C_DGER(nr,ns,1.0,Sbrp[b + bstart],1,Qasp[a + astart],1,Vrsp[0],ns);
                }


                C_DGER(nr,ns,-1.0,Qarp[a + astart],1,SAbsp[b + bstart],1,Vrsp[0],ns);
                C_DGER(nr,ns,-1.0,SBarp[a + astart],1,Qbsp[b + bstart],1,Vrsp[0],ns);

                for (int r = 0; r < nr; r++) {
                    for (int s = 0; s < ns; s++) {
                        ExchDisp20 += Trsp[r][s] * Vrsp[r][s];
                    }
                }
            }
        }
    }

    energies_["Disp20"] = Disp20;
    energies_["Exch-Disp20"] = ExchDisp20;
    outfile->Printf("    Disp20              = %18.12lf H\n",Disp20);
    outfile->Printf("    Exch-Disp20         = %18.12lf H\n",ExchDisp20);
    outfile->Printf("\n");
}


}} // End namespaces
