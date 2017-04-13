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

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

#include "psi4/psifiles.h"
#include "psi4/physconst.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libqt/qt.h"
#include "psi4/liboptions/liboptions_python.h"
#include "psi4/psifiles.h"
#include "psi4/libfock/jk.h"
#include "psi4/libfock/v.h"
#include "psi4/libfunctional/superfunctional.h"

#ifdef USING_PCMSolver
#include "psi4/libpsipcm/psipcm.h"
#endif

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/basisset_parser.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/extern.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/pointgrp.h"
#include "psi4/libmints/oeprop.h"
#include "hf.h"

#include "psi4/psi4-dec.h"
#include "psi4/libefp_solver/efp_solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi { namespace scf {

HF::HF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> func,
       Options &options, std::shared_ptr<PSIO> psio)
    : Wavefunction(options),
      functional_(func),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6) {
    shallow_copy(ref_wfn);
    psio_ = psio;
    common_init();
}


HF::~HF()
{
}

void HF::common_init()
{

    attempt_number_ = 1;
    ref_C_ = false;
    reset_occ_ = false;

    max_attempts_ = options_.get_int("MAX_ATTEMPTS");

    // This quantity is needed fairly soon
    nirrep_ = factory_->nirrep();

    integral_threshold_ = options_.get_double("INTS_TOLERANCE");

    scf_type_ = options_.get_str("SCF_TYPE");

    H_.reset(factory_->create_matrix("One-electron Hamiltonian"));
    X_.reset(factory_->create_matrix("X"));

    nmo_ = 0;
    nso_ = 0;
    const Dimension& dimpi = factory_->colspi();
    for (int h = 0; h< factory_->nirrep(); h++){
        nsopi_[h] = dimpi[h];
        nmopi_[h] = nsopi_[h]; //For now, may change in S^-1/2
        nso_ += nsopi_[h];
        nmo_ += nmopi_[h]; //For now, may change in S^-1/2
    }

    // Read in energy convergence threshold
    energy_threshold_ = options_.get_double("E_CONVERGENCE");

    // Read in density convergence threshold
    density_threshold_ = options_.get_double("D_CONVERGENCE");;

    density_fitted_ = false;

    Eold_    = 0.0;
    E_       = 0.0;
    maxiter_ = 40;

    // Read information from input file
    maxiter_ = options_.get_int("MAXITER");

    // Should we continue if we fail to converge?
    fail_on_maxiter_ = options_.get_bool("FAIL_ON_MAXITER");

    // Set name
    if(options_.get_str("REFERENCE") == "RKS" || options_.get_str("REFERENCE") == "UKS")
        name_ = "DFT";
    else
        name_ = "SCF";

    // Read in DOCC and SOCC from memory
    int nirreps = factory_->nirrep();
    input_docc_ = false;
    if (options_["DOCC"].has_changed()) {
        input_docc_ = true;
        // Map the symmetry of the input DOCC, to account for displacements
        std::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
        if(old_pg){
            // This is one of a series of displacements;  check the dimension against the parent point group
            size_t full_nirreps = old_pg->char_table().nirrep();
	    if(options_["DOCC"].size() != full_nirreps)
	        throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            Dimension temp_docc(full_nirreps);
            for(int h = 0; h < full_nirreps; ++h) {
                temp_docc[h] = options_["DOCC"][h].to_integer();
	    }
            doccpi_ = map_irreps(temp_docc);
        }else{
            // This is a normal calculation; check the dimension against the current point group then read
            if(options_["DOCC"].size() != nirreps)
                throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            for(int h = 0; h < nirreps; ++h) {
	      doccpi_[h] = options_["DOCC"][h].to_integer();
	    }
        }
    } // else take the reference wavefunctions doccpi

    input_socc_ = false;
    if (options_["SOCC"].has_changed()) {
        input_socc_ = true;
        // Map the symmetry of the input SOCC, to account for displacements
        std::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
        if(old_pg){
            // This is one of a series of displacements;  check the dimension against the parent point group
            size_t full_nirreps = old_pg->char_table().nirrep();
            if(options_["SOCC"].size() != full_nirreps)
                throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            Dimension temp_socc(full_nirreps);
            for(int h = 0; h < full_nirreps; ++h) {
                temp_socc[h] = options_["SOCC"][h].to_integer();
	    }
            soccpi_ = map_irreps(temp_socc);
        }else{
            // This is a normal calculation; check the dimension against the current point group then read
            if(options_["SOCC"].size() != nirreps)
                throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            for(int h = 0; h < nirreps; ++h) {
                soccpi_[h] = options_["SOCC"][h].to_integer();
	    }
        }
    } // else take the reference wavefunctions soccpi

    // Check that we have enough basis functions
    for(int h = 0; h < nirreps; ++h) {
      if(doccpi_[h]+soccpi_[h] > nmopi_[h]) {
	throw PSIEXCEPTION("Not enough basis functions to satisfy requested occupancies");
      }
    }

    if (input_socc_ || input_docc_) {
        for (int h = 0; h < nirrep_; h++) {
            nalphapi_[h] = doccpi_[h] + soccpi_[h];
            nbetapi_[h]  = doccpi_[h];
        }
    }


    // Set additional information
    nuclearrep_ = molecule_->nuclear_repulsion_energy();
    charge_ = molecule_->molecular_charge();
    multiplicity_ = molecule_->multiplicity();
    nelectron_ = nbeta_ + nalpha_;

    // Copy data for storage
    original_doccpi_ = doccpi_;
    original_soccpi_ = soccpi_;
    original_nalpha_ = nalpha_;
    original_nbeta_ = nbeta_;

    // Check if it is a broken symmetry solution
    // broken_symmetry_ = false;
    // int socc = 0;
    // for (int h = 0; h < nirrep_; h++) {
    //     socc += soccpi_[h];
    // }

    // if (multiplicity_ == 1 && socc == 2) {
    //     // Set up occupation for the broken symmetry solution
    //     outfile->Printf( "  Broken symmetry solution detected... \n"); //TEST
    //     broken_symmetry_ = true;
    //     int socc_count = 0;
    //     nalphapi_ = doccpi_;
    //     nbetapi_  = doccpi_;

    //     for (int h = 0; h < nirrep_; h++) {
    //         for (int i = 0; i<soccpi_[h]; i++) {
    //             socc_count++;
    //             if (socc_count == 1) {
    //                 nalphapi_[h]++;
    //             }
    //             if (socc_count == 2) {
    //                 nbetapi_[h]++;
    //             }
    //         }
    //     }
    //     if (print_ > 2) {
    //         nalphapi_.print();
    //         nbetapi_.print();
    //     }
    // }

    perturb_h_ = false;
    perturb_h_ = options_.get_bool("PERTURB_H");
    perturb_ = nothing;
    if (perturb_h_) {
        std::string perturb_with;


        if (options_["PERTURB_WITH"].has_changed()) {
            perturb_with = options_.get_str("PERTURB_WITH");
            // Do checks to see what perturb_with is.
            if (perturb_with == "DIPOLE_X") {
                perturb_ = dipole_x;
                perturb_dipoles_[0] = options_.get_double("PERTURB_MAGNITUDE");
                nuclearrep_ += perturb_dipoles_[0]*molecule_->nuclear_dipole()[0];
                outfile->Printf(" WARNING: the DIPOLE_X and PERTURB_MAGNITUDE keywords are deprecated."
                                "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE_Y") {
                perturb_ = dipole_y;
                perturb_dipoles_[1] = options_.get_double("PERTURB_MAGNITUDE");
                nuclearrep_ += perturb_dipoles_[1]*molecule_->nuclear_dipole()[1];
                outfile->Printf(" WARNING: the DIPOLE_Y and PERTURB_MAGNITUDE keywords are deprecated."
                                "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE_Z") {
                perturb_ = dipole_z;
                perturb_dipoles_[2] = options_.get_double("PERTURB_MAGNITUDE");
                nuclearrep_ += perturb_dipoles_[2]*molecule_->nuclear_dipole()[2];
                outfile->Printf(" WARNING: the DIPOLE_Z and PERTURB_MAGNITUDE keywords are deprecated."
                                "  Use DIPOLE and the PERTURB_DIPOLE array instead.");
            } else if (perturb_with == "DIPOLE") {
                perturb_ = dipole;
                if(options_["PERTURB_DIPOLE"].size() !=3)
                    throw PSIEXCEPTION("The PERTURB dipole should have exactly three floating point numbers.");
                for(int n = 0; n < 3; ++n)
                    perturb_dipoles_[n] = options_["PERTURB_DIPOLE"][n].to_double();
                nuclearrep_ += perturb_dipoles_.dot(molecule_->nuclear_dipole());
            } else if (perturb_with == "EMBPOT") {
                perturb_ = embpot;
                perturb_dipoles_[0] = 1.0;
            }
            else if (perturb_with == "DX") {
                perturb_ = dx;
                perturb_dipoles_[0] = 1.0;
            }
            else if (perturb_with == "SPHERE") {
                perturb_ = sphere;
                perturb_dipoles_[0] = 1.0;
            }
            else {

                    outfile->Printf( "Unknown PERTURB_WITH. Applying no perturbation.\n");

            }
        } else {

                outfile->Printf( "PERTURB_H is true, but PERTURB_WITH not found, applying no perturbation.\n");

        }
    }

    // How much stuff shall we echo to the user?
    if(options_["PRINT"].has_changed())
        print_ = options_.get_int("PRINT");

    if(options_["DAMPING_PERCENTAGE"].has_changed()){
        // The user has asked for damping to be turned on
        damping_enabled_ = true;
        damping_percentage_ = options_.get_double("DAMPING_PERCENTAGE") / 100.0;
        if(damping_percentage_ < 0.0 || damping_percentage_ > 1.0)
            throw PSIEXCEPTION("DAMPING_PERCENTAGE must be between 0 and 100.");
        damping_convergence_ = options_.get_double("DAMPING_CONVERGENCE");
    }else{
        damping_enabled_ = false;
    }

    // Handle common diis info
    diis_enabled_ = true;
    min_diis_vectors_ = 4;

    // Allocate memory for DIISmin_diis_vectors_
    //  First, did the user request a different number of diis vectors?
    min_diis_vectors_ = options_.get_int("DIIS_MIN_VECS");
    max_diis_vectors_ = options_.get_int("DIIS_MAX_VECS");
    diis_start_ = options_.get_int("DIIS_START");
    diis_enabled_ = options_.get_bool("DIIS");

    // Don't perform DIIS if less than 2 vectors requested, or user requested a negative number
    if (min_diis_vectors_ < 2) {
        // disable diis
        diis_enabled_ = false;
    }

    initialized_diis_manager_ = false;

    // Second-order convergence acceleration
    soscf_enabled_ = options_.get_bool("SOSCF");
    soscf_r_start_ = options_.get_double("SOSCF_START_CONVERGENCE");
    soscf_min_iter_ = options_.get_int("SOSCF_MIN_ITER");
    soscf_max_iter_ = options_.get_int("SOSCF_MAX_ITER");
    soscf_conv_ = options_.get_double("SOSCF_CONV");
    soscf_print_ = options_.get_bool("SOSCF_PRINT");

    // MOM convergence acceleration
    MOM_enabled_ = (options_.get_int("MOM_START") != 0);
    MOM_excited_ = (options_["MOM_OCC"].size() != 0 && MOM_enabled_);
    MOM_started_ = false;
    MOM_performed_ = false;

    frac_enabled_ = (options_.get_int("FRAC_START") != 0);
    frac_performed_ = false;
    print_header();

    // DFT stuff
    if (functional_->needs_xc()){
        potential_ = VBase::build_V(basisset_, functional_, options_, (options_.get_str("REFERENCE") == "RKS" ? "RV" : "UV"));
        potential_->initialize();

        // Print the KS-specific stuff
        potential_->print_header();
    } else {
        potential_ = nullptr;
    }

    // -D is zero by default
    variables_["-D Energy"] = 0.0;
    energies_["-D"] = 0.0;

    // Initialize PCM object, if requested
#ifdef USING_PCMSolver
    if(pcm_enabled_ = (options_.get_bool("PCM")))
      hf_pcm_ = static_cast<SharedPCM>(new PCM(options_, psio_, nirrep_, basisset_));
#endif
}

void HF::damp_update()
{
    throw PSIEXCEPTION("Sorry, damping has not been implemented for this "
                       "type of SCF wavefunction yet.");
}

int HF::soscf_update()
{
    throw PSIEXCEPTION("Sorry, second-order convergence has not been implemented for this "
                       "type of SCF wavefunction yet.");
}
void HF::form_V()
{
    throw PSIEXCEPTION("Sorry, DFT functionals are not suppored for this type of SCF wavefunction.");
}
void HF::form_C()
{
    throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot construct orbitals.");
}
void HF::form_D()
{
    throw PSIEXCEPTION("Sorry, the base HF wavefunction cannot construct densities.");
}
void HF::rotate_orbitals(SharedMatrix C, const SharedMatrix x)
{
    // => Rotate orbitals <= //
    SharedMatrix U(new Matrix("Ck", nirrep_, nmopi_, nmopi_));
    std::string reference = options_.get_str("REFERENCE");

    // We guess occ x vir block size by the size of x to make this method easy to use
    Dimension tsize = x->colspi() + x->rowspi();
    if ((reference != "ROHF") && (tsize != nmopi_)){
        throw PSIEXCEPTION("HF::rotate_orbitals: x dimensions do not match nmo_ dimension.");
    }
    tsize = x->colspi() + x->rowspi() - soccpi_;
    if ((reference == "ROHF") && (tsize != nmopi_)){
        throw PSIEXCEPTION("HF::rotate_orbitals: x dimensions do not match nmo_ dimension.");
    }

    // Form full antisymmetric matrix
    for (size_t h=0; h<nirrep_; h++){

        // Whatever the dimension are, we set top right/bot left
        size_t doccpih = (size_t)x->rowspi()[h];
        size_t virpih = (size_t)x->colspi()[h];
        if (!doccpih || !virpih) continue;
        double** up = U->pointer(h);
        double*  xp = x->pointer(h)[0];

        // Matrix::schmidt orthogonalizes rows not columns so we need to transpose
        for (size_t i=0, target=0; i<doccpih; i++){
            for (size_t a=(nmopi_[h] - virpih); a < nmopi_[h]; a++){
                up[a][i] = xp[target];
                up[i][a] = -1.0 * xp[target++];
            }
        }

    }
    U->expm(4, true);

    // Need to build a new one here incase nmo != nso
    SharedMatrix tmp = Matrix::doublet(C, U, false, false);
    C->copy(tmp);

    U.reset();
    tmp.reset();
}
void HF::integrals()
{
    if (print_ )
        outfile->Printf( "  ==> Integral Setup <==\n\n");

    // Build the JK from options, symmetric type
    // try {
    if(options_.get_str("SCF_TYPE") == "GTFOCK") {
      #ifdef HAVE_JK_FACTORY
        //DGAS is adding to the ghetto, this Python -> C++ -> C -> C++ -> back to C is FUBAR
        std::shared_ptr<Molecule> other_legacy = Process::environment.legacy_molecule();
        Process::environment.set_legacy_molecule(molecule_);
        if(options_.get_bool("SOSCF"))
            jk_ = std::shared_ptr<JK>(new GTFockJK(basisset_,2,false));
        else
            jk_ = std::shared_ptr<JK>(new GTFockJK(basisset_,2,true));
        Process::environment.set_legacy_molecule(other_legacy);
      #else
        throw PSIEXCEPTION("GTFock was not compiled in this version.\n");
      #endif
    } else {
        if (options_.get_str("SCF_TYPE") == "DF"){
            jk_ = JK::build_JK(get_basisset("ORBITAL"), get_basisset("DF_BASIS_SCF"), options_);
        } else {
            jk_ = JK::build_JK(get_basisset("ORBITAL"), BasisSet::zero_ao_basis_set(), options_);

        }
    }

    // Tell the JK to print
    jk_->set_print(print_);
    // Give the JK 75% of the memory
    jk_->set_memory((ULI)(options_.get_double("SCF_MEM_SAFETY_FACTOR")*(Process::environment.get_memory() / 8L)));

    // DFT sometimes needs custom stuff
    // K matrices
    jk_->set_do_K(functional_->is_x_hybrid());
    // wK matrices
    jk_->set_do_wK(functional_->is_x_lrc());
    // w Value
    jk_->set_omega(functional_->x_omega());

    // Initialize
    jk_->initialize();
    // Print the header
    jk_->print_header();
}

double HF::compute_energy()
{
  initialize();
  iterations();
  return finalize_E();
}

double HF::finalize_E()
{
    // Perform wavefunction stability analysis before doing
    // anything on a wavefunction that may not be truly converged.
    if(options_.get_str("STABILITY_ANALYSIS") != "NONE") {
        // We need the integral file, make sure it is written and
        // compute it if needed
        if(options_.get_str("REFERENCE") != "UHF") {
            psio_->open(PSIF_SO_TEI, PSIO_OPEN_OLD);
            if (psio_->tocscan(PSIF_SO_TEI, IWL_KEY_BUF) == NULL) {
                psio_->close(PSIF_SO_TEI,1);
                outfile->Printf("    SO Integrals not on disk, computing...");
                std::shared_ptr<MintsHelper> mints(new MintsHelper(basisset_, options_, 0));
                mints->integrals();
                outfile->Printf("done.\n");
            } else {
                psio_->close(PSIF_SO_TEI,1);
            }

        }
        bool follow = stability_analysis();

        while ( follow && !(attempt_number_ > max_attempts_) ) {

          attempt_number_++;
          outfile->Printf("    Running SCF again with the rotated orbitals.\n");

          if(initialized_diis_manager_) diis_manager_->reset_subspace();
          // Reading the rotated orbitals in before starting iterations
          form_D();
          E_ = compute_initial_E();
          iterations();
          follow = stability_analysis();
        }
        if ( follow && (attempt_number_ > max_attempts_) ) {
          outfile->Printf( "    There's still a negative eigenvalue. Try modifying FOLLOW_STEP_SCALE\n");
          outfile->Printf("    or increasing MAX_ATTEMPTS (not available for PK integrals).\n");
        }
    }

    // At this point, we are not doing any more SCF cycles
    // and we can compute and print final quantities.
#ifdef USING_libefp
    if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        Process::environment.get_efp()->compute();

        double efp_wfn_independent_energy = Process::environment.globals["EFP TOTAL ENERGY"] -
                                            Process::environment.globals["EFP IND ENERGY"];
        energies_["EFP"] = Process::environment.globals["EFP TOTAL ENERGY"];

        outfile->Printf("    EFP excluding EFP Induction   %20.12f [Eh]\n", efp_wfn_independent_energy);
        outfile->Printf("    SCF including EFP Induction   %20.12f [Eh]\n", E_);

        E_ += efp_wfn_independent_energy;

        outfile->Printf("    Total SCF including Total EFP %20.12f [Eh]\n", E_);
    }
#endif

    outfile->Printf( "\n  ==> Post-Iterations <==\n\n");

    check_phases();
    compute_spin_contamination();
    frac_renormalize();
    std::string reference = options_.get_str("REFERENCE");

    if (converged_ || !fail_on_maxiter_) {

        // Print the orbitals
        if(print_)
            print_orbitals();

        if (converged_) {
            outfile->Printf( "  Energy converged.\n\n");
        }
        if (!converged_) {
            outfile->Printf( "  Energy did not converge, but proceeding anyway.\n\n");
        }

        bool df = (options_.get_str("SCF_TYPE") == "DF");

        outfile->Printf( "  @%s%s Final Energy: %20.14f", df ? "DF-" : "", reference.c_str(), E_);
        if (perturb_h_) {
            outfile->Printf( " with %f %f %f perturbation", perturb_dipoles_[0], perturb_dipoles_[1], perturb_dipoles_[2]);
        }
        outfile->Printf( "\n\n");
        print_energies();

        // Need to recompute the Fock matrices, as they are modified during the SCF iteration
        // and might need to be dumped to checkpoint later
        form_F();
#ifdef USING_PCMSolver
        if(pcm_enabled_) {
            // Prepare the density
            SharedMatrix D_pcm(Da_->clone());
            if(same_a_b_orbs()) {
              D_pcm->scale(2.0); // PSI4's density doesn't include the occupation
            }
            else {
              D_pcm->add(Db_);
            }

            // Add the PCM potential to the Fock matrix
            SharedMatrix V_pcm;
            V_pcm = hf_pcm_->compute_V();
            if(same_a_b_orbs()) Fa_->add(V_pcm);
            else {
              Fa_->add(V_pcm);
              Fb_->add(V_pcm);
            }
        }
#endif

        // Properties
//  Comments so that autodoc utility will find these PSI variables
//
//  Process::environment.globals["SCF DIPOLE X"] =
//  Process::environment.globals["SCF DIPOLE Y"] =
//  Process::environment.globals["SCF DIPOLE Z"] =
//  Process::environment.globals["SCF QUADRUPOLE XX"] =
//  Process::environment.globals["SCF QUADRUPOLE XY"] =
//  Process::environment.globals["SCF QUADRUPOLE XZ"] =
//  Process::environment.globals["SCF QUADRUPOLE YY"] =
//  Process::environment.globals["SCF QUADRUPOLE YZ"] =
//  Process::environment.globals["SCF QUADRUPOLE ZZ"] =

        save_information();
    } else {
            outfile->Printf( "  Failed to converge.\n");
        E_ = 0.0;
        if(psio_->open_check(PSIF_CHKPT))
            psio_->close(PSIF_CHKPT, 1);

        // Throw if we didn't converge?
        die_if_not_converged();
    }

    // Orbitals are always saved, in case an MO guess is requested later
    // save_orbitals();

    finalize();

    //outfile->Printf("\nComputation Completed\n");

    return E_;

}

void HF::finalize()
{
    // Clean memory off, handle diis closeout, etc

    // This will be the only one
    if (!options_.get_bool("SAVE_JK")) {
        jk_.reset();
    }

    // Clean up after DIIS
    if(initialized_diis_manager_)
        diis_manager_->delete_diis_file();
    diis_manager_.reset();
    initialized_diis_manager_ = false;

    // Figure out how many frozen virtual and frozen core per irrep
    compute_fcpi();
    compute_fvpi();
    energy_ = E_;

    //Sphalf_.reset();
    X_.reset();
    T_.reset();
    diag_temp_.reset();
    diag_F_temp_.reset();
    diag_C_temp_.reset();


}

void HF::semicanonicalize()
{
    throw PSIEXCEPTION("This type of wavefunction cannot be semicanonicalized!");
}

void HF::find_occupation()
{
    // Don't mess with the occ, MOM's got it!
    if (MOM_started_) {
        MOM();
    } else {
        std::vector<std::pair<double, int> > pairs_a;
        std::vector<std::pair<double, int> > pairs_b;
        for (int h=0; h<epsilon_a_->nirrep(); ++h) {
            for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
                pairs_a.push_back(std::make_pair(epsilon_a_->get(h, i), h));
        }
        for (int h=0; h<epsilon_b_->nirrep(); ++h) {
            for (int i=0; i<epsilon_b_->dimpi()[h]; ++i)
                pairs_b.push_back(std::make_pair(epsilon_b_->get(h, i), h));
        }
        sort(pairs_a.begin(),pairs_a.end());
        sort(pairs_b.begin(),pairs_b.end());

        if(!input_docc_ && !input_socc_){
            memset(nalphapi_, 0, sizeof(int) * epsilon_a_->nirrep());
            for (int i=0; i<nalpha_; ++i)
                nalphapi_[pairs_a[i].second]++;
        }
        if(!input_docc_ && !input_socc_){
            memset(nbetapi_, 0, sizeof(int) * epsilon_b_->nirrep());
            for (int i=0; i<nbeta_; ++i)
                nbetapi_[pairs_b[i].second]++;
        }

        int old_socc[8];
        int old_docc[8];
        for(int h = 0; h < nirrep_; ++h){
            old_socc[h] = soccpi_[h];
            old_docc[h] = doccpi_[h];
        }

        if(!input_docc_ && !input_socc_){
            for (int h = 0; h < nirrep_; ++h) {
                soccpi_[h] = std::abs(nalphapi_[h] - nbetapi_[h]);
                doccpi_[h] = std::min(nalphapi_[h] , nbetapi_[h]);
            }
        }

        bool occ_changed = false;
        for(int h = 0; h < nirrep_; ++h){
            if( old_socc[h] != soccpi_[h] || old_docc[h] != doccpi_[h]){
                occ_changed = true;
                break;
            }
        }

        // If print > 2 (diagnostics), print always
        if((print_ > 2 || (print_ && occ_changed)) && iteration_ > 0){

                outfile->Printf( "    Occupation by irrep:\n");
            print_occupation();
        }
        // Start MOM if needed (called here because we need the nocc
        // to be decided by Aufbau ordering prior to MOM_start)
        MOM_start();
    }
    // Do fractional orbital normalization here.
    frac();
}

void HF::print_header()
{
    int nthread = 1;
    #ifdef _OPENMP
        nthread = Process::environment.get_n_threads();
    #endif


    outfile->Printf( "\n");
    outfile->Printf( "         ---------------------------------------------------------\n");
    outfile->Printf( "                                   SCF\n");
    outfile->Printf( "            by Justin Turney, Rob Parrish, and Andy Simmonett\n");
    outfile->Printf( "                             %4s Reference\n", options_.get_str("REFERENCE").c_str());
    outfile->Printf( "                      %3d Threads, %6ld MiB Core\n", nthread, memory_ / 1048576L);
    outfile->Printf( "         ---------------------------------------------------------\n");
    outfile->Printf( "\n");
    outfile->Printf( "  ==> Geometry <==\n\n");


    molecule_->print();


    outfile->Printf( "  Running in %s symmetry.\n\n", molecule_->point_group()->symbol().c_str());

    molecule_->print_rotational_constants();

    outfile->Printf( "  Nuclear repulsion = %20.15f\n\n", nuclearrep_);
    outfile->Printf( "  Charge       = %d\n", charge_);
    outfile->Printf( "  Multiplicity = %d\n", multiplicity_);
    outfile->Printf( "  Electrons    = %d\n", nelectron_);
    outfile->Printf( "  Nalpha       = %d\n", nalpha_);
    outfile->Printf( "  Nbeta        = %d\n\n", nbeta_);

    outfile->Printf( "  ==> Algorithm <==\n\n");
    outfile->Printf( "  SCF Algorithm Type is %s.\n", options_.get_str("SCF_TYPE").c_str());
    outfile->Printf( "  DIIS %s.\n", diis_enabled_ ? "enabled" : "disabled");
    if (MOM_excited_)
        outfile->Printf( "  Excited-state MOM enabled.\n");
    else
        outfile->Printf( "  MOM %s.\n", MOM_enabled_ ? "enabled" : "disabled");
    outfile->Printf( "  Fractional occupation %s.\n", frac_enabled_ ? "enabled" : "disabled");
    outfile->Printf( "  Guess Type is %s.\n", options_.get_str("GUESS").c_str());
    outfile->Printf( "  Energy threshold   = %3.2e\n", energy_threshold_);
    outfile->Printf( "  Density threshold  = %3.2e\n", density_threshold_);
    outfile->Printf( "  Integral threshold = %3.2e\n\n", integral_threshold_);


    outfile->Printf( "  ==> Primary Basis <==\n\n");

    basisset_->print_by_level("outfile", print_);

}
void HF::print_preiterations()
{
    CharacterTable ct = molecule_->point_group()->char_table();


    outfile->Printf( "   -------------------------------------------------------\n");
    outfile->Printf( "    Irrep   Nso     Nmo     Nalpha   Nbeta   Ndocc  Nsocc\n");
    outfile->Printf( "   -------------------------------------------------------\n");
    for (int h= 0; h < nirrep_; h++) {
        outfile->Printf("     %-3s   %6d  %6d  %6d  %6d  %6d  %6d\n",
                        ct.gamma(h).symbol(), nsopi_[h], nmopi_[h],
                        nalphapi_[h], nbetapi_[h], doccpi_[h], soccpi_[h]);
    }
    outfile->Printf( "   -------------------------------------------------------\n");
    outfile->Printf("    Total  %6d  %6d  %6d  %6d  %6d  %6d\n", nso_, nmo_,
                    nalpha_, nbeta_, nbeta_, nalpha_ - nbeta_);
    outfile->Printf("   -------------------------------------------------------\n\n");
}

void HF::form_H()
{
    T_ = SharedMatrix(factory_->create_matrix(PSIF_SO_T));
    V_ = SharedMatrix(factory_->create_matrix(PSIF_SO_V));

    // Assumes these have already been created and stored
    T_->load(psio_, PSIF_OEI);
    V_->load(psio_, PSIF_OEI);

    if (debug_ > 2)
        T_->print("outfile");

    if (debug_ > 2)
        V_->print("outfile");

    if (perturb_h_) {
      if(perturb_ == embpot || perturb_ == sphere || perturb_ == dx) { // embedding potential read from file
        if(nirrep_ > 1)
          throw PSIEXCEPTION("RHF_embed: embedding, dx, and spherical potentials require 'symmetry c1'.");
        int nso = 0;
        for(int h=0; h < nirrep_; h++) nso += nsopi_[h];
        int nao = basisset_->nao();

        // Set up AO->SO transformation matrix (u)
        MintsHelper helper(basisset_, options_, 0);
        SharedMatrix aotoso = helper.petite_list(true)->aotoso();
        int *col_offset = new int[nirrep_];
        col_offset[0] = 0;
        for(int h=1; h < nirrep_; h++)
          col_offset[h] = col_offset[h-1] + aotoso->coldim(h-1);

        double **u = block_matrix(nao, nso);
        for(int h=0; h < nirrep_; h++)
          for(int j=0; j < aotoso->coldim(h); j++)
            for(int i=0; i < nao; i++)
              u[i][j+col_offset[h]] = aotoso->get(h, i, j);
        delete[] col_offset;

        double *phi_ao, *phi_so, **V_eff;
        phi_ao = init_array(nao);
        phi_so = init_array(nso);
        V_eff = block_matrix(nso, nso);

        if(perturb_ == embpot) {

          FILE* input = fopen("EMBPOT", "r");
          int npoints;
          int statusvalue=fscanf(input, "%d", &npoints);
          outfile->Printf( "  npoints = %d\n", npoints);
          double x, y, z, w, v;
          double max = 0;
          for(int k=0; k < npoints; k++) {
            statusvalue=fscanf(input, "%lf %lf %lf %lf %lf", &x, &y, &z, &w, &v);
            if(fabs(v) > max) max = fabs(v);

            basisset_->compute_phi(phi_ao, x, y, z);
            // Transform phi_ao to SO basis
            C_DGEMV('t', nao, nso, 1.0, &(u[0][0]), nso, &(phi_ao[0]), 1, 0.0, &(phi_so[0]), 1);
            for(int i=0; i < nso; i++)
              for(int j=0; j < nso; j++)
                V_eff[i][j] += w * v * phi_so[i] * phi_so[j];
          } // npoints

          outfile->Printf( "  Max. embpot value = %20.10f\n", max);
          fclose(input);

        } // embpot
        else if(perturb_ == dx) {
          dx_read(V_eff, phi_ao, phi_so, nao, nso, u);

        } // dx file
        else if(perturb_ == sphere) {
          radius_ = options_.get_double("RADIUS");
          thickness_ = options_.get_double("THICKNESS");
          r_points_ = options_.get_int("R_POINTS");
          theta_points_ = options_.get_int("THETA_POINTS");
          phi_points_ = options_.get_int("PHI_POINTS");
          outfile->Printf( "  Hard spherical potential radius         = %3.2f bohr\n", radius_);
          outfile->Printf( "  Spherical potential thickness           = %3.2f bohr\n", thickness_);
          outfile->Printf( "  Number of radial integration points     = %d\n", r_points_);
          outfile->Printf( "  Number of colatitude integration points = %d\n", theta_points_);
          outfile->Printf( "  Number of azimuthal integration points  = %d\n", phi_points_);

          double r_step = thickness_/r_points_; // bohr
          double theta_step = 2*pc_pi/theta_points_; // 1 degree in radians
          double phi_step = 2*pc_pi/phi_points_; // 1 degree in radians
          double weight = r_step * theta_step * phi_step;
          for(double r=radius_; r < radius_+thickness_; r += r_step) {
            for(double theta=0.0; theta < pc_pi; theta += theta_step) {  /* colatitude */
              for(double phi=0.0; phi < 2*pc_pi; phi += phi_step) { /* azimuthal */

                double x = r * sin(theta) * cos(phi);
                double y = r * sin(theta) * sin(phi);
                double z = r * cos(theta);

                double jacobian = weight * r * r * sin(theta);

                basisset_->compute_phi(phi_ao, x, y, z);

                C_DGEMV('t', nao, nso, 1.0, &(u[0][0]), nso, &(phi_ao[0]), 1,
                        0.0, &(phi_so[0]), 1);

                for(int i=0; i < nso; i++)
                  for(int j=0; j < nso; j++)
                    V_eff[i][j] += jacobian * (-1.0e6) * phi_so[i] * phi_so[j];
              }
            }
          }
        } // sphere


          outfile->Printf( "  Perturbing H by %f %f %f V_eff.\n", perturb_dipoles_[0], perturb_dipoles_[1], perturb_dipoles_[2]);
          if(options_.get_int("PRINT") > 3) mat_print(V_eff, nso, nso, "outfile");


        if(perturb_ == dx) {
          for(int i=0; i < nso; i++)
            for(int j=0; j < nso; j++)
              V_->set(i, j, V_eff[i][j]); // ignore nuclear potential
        }
        else {
          for(int i=0; i < nso; i++)
            for(int j=0; j < nso; j++)
              V_->set(i, j, (V_eff[i][j] + V_->get(i,j)));
        }

        free(phi_ao);
        free(phi_so);
        free_block(V_eff);
      }  // embpot or sphere
    } // end perturb_h_

    // If an external field exists, add it to the one-electron Hamiltonian
    py::object pyExtern = dynamic_cast<PythonDataType*>(options_["EXTERN"].get())->to_python();
    std::shared_ptr<ExternalPotential> external;
    if (pyExtern)
        external = pyExtern.cast<std::shared_ptr<ExternalPotential>>();
    if (external) {
        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY") == false && H_->nirrep() != 1)
            throw PSIEXCEPTION("SCF: External Fields are not consistent with symmetry. Set symmetry c1.");

        SharedMatrix Vprime = external->computePotentialMatrix(basisset_);

        if (options_.get_bool("EXTERNAL_POTENTIAL_SYMMETRY")) {
            // Attempt to apply symmetry. No error checking is performed.
            SharedMatrix Vprimesym = factory_->create_shared_matrix("External Potential");
            Vprimesym->apply_symmetry(Vprime, AO2SO_);
            Vprime = Vprimesym;
        }

        if (print_) {
            external->set_print(print_);
            external->print();
        }
        if (print_ > 3)
            Vprime->print();
        V_->add(Vprime);


        // Extra nuclear repulsion
        double enuc2 = external->computeNuclearEnergy(molecule_);
        if (print_) {
               outfile->Printf( "  Old nuclear repulsion        = %20.15f\n", nuclearrep_);
               outfile->Printf( "  Additional nuclear repulsion = %20.15f\n", enuc2);
               outfile->Printf( "  Total nuclear repulsion      = %20.15f\n\n", nuclearrep_ + enuc2);
        }
        nuclearrep_ += enuc2;

    }  // end external

    // Save perturbed V_ for future (e.g. correlated) calcs
    V_->save(psio_, PSIF_OEI);

    H_->copy(T_);
    H_->add(V_);

    if (print_ > 3)
        H_->print("outfile");
}

void HF::form_Shalf()
{
    // ==> SYMMETRIC ORTHOGONALIZATION <== //

    // S_ is computed by wavefunction

    SharedMatrix eigvec= factory_->create_shared_matrix("L");
    SharedMatrix eigtemp= factory_->create_shared_matrix("Temp");
    SharedMatrix eigtemp2= factory_->create_shared_matrix("Temp2");
    SharedMatrix eigvec_store= factory_->create_shared_matrix("TempStore");
    SharedVector eigval(factory_->create_vector());
    SharedVector eigval_store(factory_->create_vector());

    //Used to do this 3 times, now only once
    S_->diagonalize(eigvec, eigval);
    eigvec_store->copy(eigvec);
    eigval_store->copy(eigval.get());

    // Convert the eigenvales to 1/sqrt(eigenvalues)
    const Dimension& dimpi = eigval->dimpi();
    double min_S = fabs(eigval->get(0,0));
    for (int h=0; h<nirrep_; ++h) {
        for (int i=0; i<dimpi[h]; ++i) {
            if (min_S > eigval->get(h,i))
                min_S = eigval->get(h,i);
            double scale = 1.0 / sqrt(eigval->get(h, i));
            eigval->set(h, i, scale);
        }
    }
    if (print_ )
        outfile->Printf("  Minimum eigenvalue in the overlap matrix is %14.10E.\n",min_S);
    // Create a vector matrix from the converted eigenvalues
    eigtemp2->set_diagonal(eigval);

    eigtemp->gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    X_->gemm(false, false, 1.0, eigvec, eigtemp, 0.0);

    // ==> CANONICAL ORTHOGONALIZATION <== //

    // Decide symmetric or canonical
    double S_cutoff = options_.get_double("S_TOLERANCE");
    if (min_S > S_cutoff && options_.get_str("S_ORTHOGONALIZATION") == "SYMMETRIC") {

        if (print_)
            outfile->Printf("  Using Symmetric Orthogonalization.\n\n");

    } else {

        if (print_)
            outfile->Printf("  Using Canonical Orthogonalization with cutoff of %14.10E.\n",S_cutoff);

        //Diagonalize S (or just get a fresh copy)
        eigvec->copy(eigvec_store.get());
        eigval->copy(eigval_store.get());
        int delta_mos = 0;
        for (int h=0; h<nirrep_; ++h) {
            //in each irrep, scale significant cols i  by 1.0/sqrt(s_i)
            int start_index = 0;
            for (int i=0; i<dimpi[h]; ++i) {
                if (S_cutoff  < eigval->get(h,i)) {
                    double scale = 1.0 / sqrt(eigval->get(h, i));
                    eigvec->scale_column(h, i, scale);
                } else {
                    start_index++;
                    nmopi_[h]--;
                    nmo_--;
                }
            }
            if (print_>2)
                outfile->Printf("  Irrep %d, %d of %d possible MOs eliminated.\n",h,start_index,nsopi_[h]);

            delta_mos += start_index;
        }

        X_->init(nirrep_,nsopi_,nmopi_,"X (Canonical Orthogonalization)");
        for (int h=0; h<eigval->nirrep(); ++h) {
            //Copy significant columns of eigvec into X in
            //descending order
            int start_index = 0;
            for (int i=0; i<dimpi[h]; ++i) {
                if (S_cutoff  < eigval->get(h,i)) {
                } else {
                    start_index++;
                }
            }
            for (int i=0; i<dimpi[h]-start_index; ++i) {
                for (int m = 0; m < dimpi[h]; m++)
                    X_->set(h,m,i,eigvec->get(h,m,dimpi[h]-i-1));
            }
        }

        if (print_)
            outfile->Printf("  Overall, %d of %d possible MOs eliminated.\n\n",delta_mos,nso_);

	// Double check occupation vectors
	for(int h = 0; h < eigval->nirrep(); ++h) {
	  if(doccpi_[h]+soccpi_[h] > nmopi_[h]) {
	    throw PSIEXCEPTION("Not enough molecular orbitals to satisfy requested occupancies");
	  }
	}

        // Refreshes twice in RHF, no big deal
        epsilon_a_->init(nmopi_);
        Ca_->init(nirrep_,nsopi_,nmopi_,"Alpha MO coefficients");
        epsilon_b_->init(nmopi_);
        Cb_->init(nirrep_,nsopi_,nmopi_,"Beta MO coefficients");

        // Extra matrix dimension changes for specific derived classes
        prepare_canonical_orthogonalization();

    }

    // Temporary variables needed by diagonalize_F
    diag_temp_   = SharedMatrix(new Matrix(nirrep_, nmopi_, nsopi_));
    diag_F_temp_ = SharedMatrix(new Matrix(nirrep_, nmopi_, nmopi_));
    diag_C_temp_ = SharedMatrix(new Matrix(nirrep_, nmopi_, nmopi_));

    if (print_ > 3) {
        S_->print("outfile");
        X_->print("outfile");
    }

}

void HF::compute_fcpi()
{
    // FROZEN_DOCC takes precedence, FREEZE_CORE directive has second priority
    if (options_["FROZEN_DOCC"].has_changed()) {
        if (options_["FROZEN_DOCC"].size() != epsilon_a_->nirrep()) {
            throw PSIEXCEPTION("The FROZEN_DOCC array has the wrong dimensions");
        }
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            frzcpi_[h] = options_["FROZEN_DOCC"][h].to_integer();
        }
    } else {

        int nfzc = 0;
        if (options_.get_int("NUM_FROZEN_DOCC") != 0) {
            nfzc = options_.get_int("NUM_FROZEN_DOCC");
        } else {
            nfzc = molecule_->nfrozen_core(options_.get_str("FREEZE_CORE"));
        }
        // Print out orbital energies.
        std::vector<std::pair<double, int> > pairs;
        for (int h=0; h<epsilon_a_->nirrep(); ++h) {
            for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
                pairs.push_back(std::make_pair(epsilon_a_->get(h, i), h));
            frzcpi_[h] = 0;
        }
        sort(pairs.begin(),pairs.end());

        for (int i=0; i<nfzc; ++i)
            frzcpi_[pairs[i].second]++;
    }
    // total frozen core
    nfrzc_ = 0;
    for (int h = 0; h < epsilon_a_->nirrep(); h++) nfrzc_ += frzcpi_[h];
}

void HF::compute_fvpi()
{
    // FROZEN_UOCC takes precedence, FREEZE_UOCC directive has second priority
    if (options_["FROZEN_UOCC"].has_changed()) {
        if (options_["FROZEN_UOCC"].size() != epsilon_a_->nirrep()) {
            throw PSIEXCEPTION("The FROZEN_UOCC array has the wrong dimensions");
        }
        for (int h = 0; h < epsilon_a_->nirrep(); h++) {
            frzvpi_[h] = options_["FROZEN_UOCC"][h].to_integer();
        }
    } else {
        int nfzv = options_.get_int("NUM_FROZEN_UOCC");
        // Print out orbital energies.
        std::vector<std::pair<double, int> > pairs;
        for (int h=0; h<epsilon_a_->nirrep(); ++h) {
            for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
                pairs.push_back(std::make_pair(epsilon_a_->get(h, i), h));
            frzvpi_[h] = 0;
        }
        sort(pairs.begin(),pairs.end(), std::greater<std::pair<double, int> >());

        for (int i=0; i<nfzv; ++i)
            frzvpi_[pairs[i].second]++;
    }
}

void HF::print_orbitals(const char* header, std::vector<std::pair<double, std::pair<const char*, int> > > orbs)
{

        outfile->Printf( "    %-70s\n\n    ", header);
        int count = 0;
        for (int i = 0; i < orbs.size(); i++) {
            outfile->Printf( "%4d%-4s%11.6f  ", orbs[i].second.second, orbs[i].second.first, orbs[i].first);
            if (count++ % 3 == 2 && count != orbs.size())
                outfile->Printf( "\n    ");
        }
        outfile->Printf( "\n\n");

}

void HF::print_orbitals()
{
    char **labels = molecule_->irrep_labels();

        outfile->Printf( "    Orbital Energies (a.u.)\n    -----------------------\n\n");

    std::string reference = options_.get_str("REFERENCE");
    if((reference == "RHF") || (reference == "RKS")){

        std::vector<std::pair<double, std::pair<const char*, int> > > occ;
        std::vector<std::pair<double, std::pair<const char*, int> > > vir;

        for (int h = 0; h < nirrep_; h++) {

            std::vector<std::pair<double, int> > orb_e;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_e.push_back(std::make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occ.push_back(std::make_pair(epsilon_a_->get(h,a), std::make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                vir.push_back(std::make_pair(epsilon_a_->get(h,a), std::make_pair(labels[h],orb_order[a] + 1)));

        }
        std::sort(occ.begin(), occ.end());
        std::sort(vir.begin(), vir.end());

        print_orbitals("Doubly Occupied:", occ);
        print_orbitals("Virtual:", vir);

    }else if((reference == "UHF") || (reference == "UKS") ||
        (reference == "CUHF")){

        std::vector<std::pair<double, std::pair<const char*, int> > > occA;
        std::vector<std::pair<double, std::pair<const char*, int> > > virA;
        std::vector<std::pair<double, std::pair<const char*, int> > > occB;
        std::vector<std::pair<double, std::pair<const char*, int> > > virB;

        for (int h = 0; h < nirrep_; h++) {

            std::vector<std::pair<double, int> > orb_eA;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_eA.push_back(std::make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_eA.begin(), orb_eA.end());

            std::vector<int> orb_orderA(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_orderA[orb_eA[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occA.push_back(std::make_pair(epsilon_a_->get(h,a), std::make_pair(labels[h],orb_orderA[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                virA.push_back(std::make_pair(epsilon_a_->get(h,a), std::make_pair(labels[h],orb_orderA[a] + 1)));

            std::vector<std::pair<double, int> > orb_eB;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_eB.push_back(std::make_pair(epsilon_b_->get(h,a), a));
            std::sort(orb_eB.begin(), orb_eB.end());

            std::vector<int> orb_orderB(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_orderB[orb_eB[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                occB.push_back(std::make_pair(epsilon_b_->get(h,a), std::make_pair(labels[h],orb_orderB[a] + 1)));
            for (int a = nbetapi_[h]; a < nmopi_[h]; a++)
                virB.push_back(std::make_pair(epsilon_b_->get(h,a), std::make_pair(labels[h],orb_orderB[a] + 1)));

        }
        std::sort(occA.begin(), occA.end());
        std::sort(virA.begin(), virA.end());
        std::sort(occB.begin(), occB.end());
        std::sort(virB.begin(), virB.end());

        print_orbitals("Alpha Occupied:", occA);
        print_orbitals("Alpha Virtual:", virA);
        print_orbitals("Beta Occupied:", occB);
        print_orbitals("Beta Virtual:", virB);

    }else if(reference == "ROHF"){

        std::vector<std::pair<double, std::pair<const char*, int> > > docc;
        std::vector<std::pair<double, std::pair<const char*, int> > > socc;
        std::vector<std::pair<double, std::pair<const char*, int> > > vir;

        for (int h = 0; h < nirrep_; h++) {

            std::vector<std::pair<double, int> > orb_e;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_e.push_back(std::make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                docc.push_back(std::make_pair(epsilon_a_->get(h,a), std::make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nbetapi_[h] ; a < nalphapi_[h]; a++)
                socc.push_back(std::make_pair(epsilon_a_->get(h,a), std::make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nalphapi_[h] ; a < nmopi_[h]; a++)
                vir.push_back(std::make_pair(epsilon_a_->get(h,a), std::make_pair(labels[h],orb_order[a] + 1)));

        }
        std::sort(docc.begin(), docc.end());
        std::sort(socc.begin(), socc.end());
        std::sort(vir.begin(), vir.end());

        print_orbitals("Doubly Occupied:", docc);
        print_orbitals("Singly Occupied:", socc);
        print_orbitals("Virtual:", vir);

    }else{
        throw PSIEXCEPTION("Unknown reference in HF::print_orbitals");
    }

    for(int h = 0; h < nirrep_; ++h)
        free(labels[h]);
    free(labels);


    outfile->Printf( "    Final Occupation by Irrep:\n");
    print_occupation();
}

void HF::guess()
{
    // don't save guess energy as "the" energy because we need to avoid
    // a false positive test for convergence on the first iteration (that
    // was happening before in tests/scf-guess-read before I removed
    // the statements putting this into E_).  -CDS 3/25/13
    double guess_E;

    //What does the user want?
    //Options will be:
    // ref_C_-C matrices were detected in the incoming wavefunction
    // "CORE"-CORE Hamiltonain
    // "GWH"-Generalized Wolfsberg-Helmholtz
    // "SAD"-Superposition of Atomic Denisties
    std::string guess_type = options_.get_str("GUESS");

    // DGAS broke SAD
    // if (guess_type == "SAD"){
    //     outfile->Printf("\nWarning! SAD is temporarily broken, switching to CORE!\n\n");
    //     guess_type = "CORE";
    // }
    // Take care of options that should be overridden
    if (guess_type == "AUTO"){
        outfile->Printf("\nWarning! Guess was AUTO, switching to CORE!\n\n");
        outfile->Printf("           This option should have been configured at the driver level.\n\n");
        guess_type = "CORE";
    }

    if ((guess_type == "READ") && !guess_Ca_){
        outfile->Printf("\nWarning! Guess was READ without Ca set, switching to CORE!\n");
        outfile->Printf("           This option should have been configured at the driver level.\n\n");
        guess_type = "CORE";
    }

    if (guess_Ca_){
        if (print_)
            outfile->Printf( "  SCF Guess: Orbitals guess was supplied from a previous computation.\n\n");

        std::string reference = options_.get_str("REFERENCE");
        bool single_orb = (reference == "RHF");

        if (single_orb){
            guess_Cb_ = guess_Ca_;
        } else {
            if (!guess_Cb_){
                throw PSIEXCEPTION("Guess Ca was set, but did not find a matching Cb!\n");
            }
        }


        if ((guess_Ca_->nirrep() != nirrep_) or (guess_Cb_->nirrep() != nirrep_)) {
            throw PSIEXCEPTION("Number of guess of the input orbitals do not match number of irreps of the wavefunction.");
        }
        if ((guess_Ca_->rowspi() != nsopi_) or (guess_Cb_->rowspi() != nsopi_)) {
            throw PSIEXCEPTION("Nso of the guess orbitals do not match Nso of the wavefunction.");
        }

       for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < guess_Ca_->colspi()[h]; i++) {
                C_DCOPY(nsopi_[h], &guess_Ca_->pointer(h)[0][i], guess_Ca_->colspi()[h],
                        &Ca_->pointer(h)[0][i], nmopi_[h]);
            }
        }

        if (!single_orb){
           for (int h = 0; h < nirrep_; h++) {
                for (int i = 0; i < guess_Cb_->colspi()[h]; i++) {
                    C_DCOPY(nsopi_[h], &guess_Cb_->pointer(h)[0][i], guess_Cb_->colspi()[h],
                            &Cb_->pointer(h)[0][i], nmopi_[h]);
                }
            }
        } else {
            Cb_ = Ca_;
        }

        // Figure out occupations from given input
        if (!(input_socc_ || input_docc_)) {
            nalphapi_ = guess_Ca_->colspi();
            nbetapi_ = guess_Cb_->colspi();
            nalpha_ = nalphapi_.sum();
            nbeta_ = nbetapi_.sum();
            soccpi_ = nalphapi_ - nbetapi_;
            doccpi_ = nalphapi_ - soccpi_;
        }

        format_guess();
        form_D();

        // This is a guess iteration similar to SAD
        iteration_ = -1;
        guess_E = compute_initial_E();

    } else if (guess_type == "SAD") {

        if (print_)
            outfile->Printf( "  SCF Guess: Superposition of Atomic Densities via on-the-fly atomic UHF.\n\n");

        //Superposition of Atomic Density
        iteration_ = -1;
        reset_occ_ = true;
        compute_SAD_guess();
        guess_E = compute_initial_E();

    } else if (guess_type == "GWH") {
        //Generalized Wolfsberg Helmholtz (Sounds cool, easy to code)
        if (print_)
            outfile->Printf( "  SCF Guess: Generalized Wolfsberg-Helmholtz.\n\n");

        Fa_->zero(); //Try Fa_{mn} = S_{mn} (H_{mm} + H_{nn})/2
        int h, i, j;
        const int *opi = S_->rowspi();
        int nirreps = S_->nirrep();
        for (h=0; h<nirreps; ++h) {
            for (i=0; i<opi[h]; ++i) {
                Fa_->set(h,i,i,H_->get(h,i,i));
                for (j=0; j<i; ++j) {
                    Fa_->set(h,i,j,0.875*S_->get(h,i,j)*(H_->get(h,i,i)+H_->get(h,j,j)));
                    Fa_->set(h,j,i,Fa_->get(h,i,j));
                }
            }
        }
        Fb_->copy(Fa_);
        form_initial_C();
        find_occupation();
        form_D();
        guess_E = compute_initial_E();

    } else if (guess_type == "CORE") {

        if (print_)
            outfile->Printf( "  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");

        Fa_->copy(H_); //Try the core Hamiltonian as the Fock Matrix
        Fb_->copy(H_);

        form_initial_C();
        find_occupation();
        form_D();
        guess_E = compute_initial_E();
    } else {
        throw PSIEXCEPTION("  SCF Guess: No guess was found!");

    }

    if (print_ > 3) {
        Ca_->print();
        Cb_->print();
        Da_->print();
        Db_->print();
        Fa_->print();
        Fb_->print();
    }


    E_ = 0.0; // don't use this guess in our convergence checks
}

void HF::format_guess()
{
    // Nothing to do, only for special cases
}


void HF::check_phases()
{
    for (int h=0; h<nirrep_; ++h) {
        for (int p = 0; p < Ca_->colspi(h); ++p) {
            for (int mu = 0; mu < Ca_->rowspi(h); ++mu) {
                if (fabs(Ca_->get(h, mu, p)) > 1.0E-3) {
                    if (Ca_->get(h, mu, p) < 1.0E-3) {
                        Ca_->scale_column(h, p, -1.0);
                    }
                    break;
                }
            }
        }
    }

    if (Ca_ != Cb_) {
        for (int h=0; h<nirrep_; ++h) {
            for (int p = 0; p < Cb_->colspi(h); ++p) {
                for (int mu = 0; mu < Cb_->rowspi(h); ++mu) {
                    if (fabs(Cb_->get(h, mu, p)) > 1.0E-3) {
                        if (Cb_->get(h, mu, p) < 1.0E-3) {
                            Cb_->scale_column(h, p, -1.0);
                        }
                        break;
                    }
                }
            }
        }
    }
}


void HF::initialize()
{
    converged_ = false;

    iteration_ = 0;

    if (print_)
        outfile->Printf( "  ==> Pre-Iterations <==\n\n");

    if (print_)
        print_preiterations();

    // Andy trick 2.0
    old_scf_type_ = options_.get_str("SCF_TYPE");
    if (options_.get_bool("DF_SCF_GUESS") && (old_scf_type_ == "DIRECT") ) {
         outfile->Printf( "  Starting with a DF guess...\n\n");
         if(!options_["DF_BASIS_SCF"].has_changed()) {
             // TODO: Match Dunning basis sets
             molecule_->set_basis_all_atoms("CC-PVDZ-JKFIT", "DF_BASIS_SCF");
         }
         scf_type_ = "DF";
         options_.set_str("SCF","SCF_TYPE","DF"); // Scope is reset in proc.py. This is not pretty, but it works
    }

    if(attempt_number_ == 1){
        std::shared_ptr<MintsHelper> mints (new MintsHelper(basisset_, options_, 0, ecpbasisset_));
        if ((options_.get_str("RELATIVISTIC") == "X2C") ||
            (options_.get_str("RELATIVISTIC") == "DKH")) {
            mints->set_rel_basisset(get_basisset("BASIS_RELATIVISTIC"));
        }

        mints->one_electron_integrals();

        integrals();

        timer_on("HF: Form H");
        form_H(); //Core Hamiltonian
        timer_off("HF: Form H");

#ifdef USING_libefp
        // EFP: Add in permanent moment contribution and cache
        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
            std::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_permanent();
            H_->add(Vefp);
            Horig_ = SharedMatrix(new Matrix("H orig Matrix", basisset_->nbf(), basisset_->nbf()));
            Horig_->copy(H_);
            outfile->Printf( "  QM/EFP: iterating Total Energy including QM/EFP Induction\n");
        }
#endif

        timer_on("HF: Form S/X");
        form_Shalf(); //S and X Matrix
        timer_off("HF: Form S/X");

        timer_on("HF: Guess");
        guess(); // Guess
        timer_off("HF: Guess");

    }else{
        // We're reading the orbitals from the previous set of iterations.
        form_D();
        E_ = compute_initial_E();
    }

#ifdef USING_libefp
    if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        Process::environment.get_efp()->set_qm_atoms();
    }
#endif
}

void HF::iterations()
{
    std::string reference = options_.get_str("REFERENCE");

    MOM_performed_ = false;
    diis_performed_ = false;

    bool df = (options_.get_str("SCF_TYPE") == "DF");

    outfile->Printf( "  ==> Iterations <==\n\n");
    outfile->Printf( "%s                        Total Energy        Delta E     RMS |[F,P]|\n\n", df ? "   " : "");


    // SCF iterations
    do {
        iteration_++;

        save_density_and_energy();

        // Call any preiteration callbacks
//        call_preiteration_callbacks();

#ifdef USING_libefp
        // add efp contribution to Fock matrix
        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
            H_->copy(Horig_);
            std::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_induced();
            H_->add(Vefp);
        }
#endif

        E_ = 0.0;

        timer_on("HF: Form G");
        form_G();
        timer_off("HF: Form G");

        // Reset fractional SAD occupation
        if ((iteration_ == 0) && reset_occ_){
            reset_occupation();
        }

        timer_on("HF: Form F");
        form_F();
        timer_off("HF: Form F");

        if (print_>3) {
            Fa_->print("outfile");
            Fb_->print("outfile");
        }

        E_ += compute_E();

#ifdef USING_libefp
        // add efp contribution to energy
        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
            double efp_wfn_dependent_energy = Process::environment.get_efp()->scf_energy_update();
            E_ += efp_wfn_dependent_energy;
        }
#endif

#ifdef USING_PCMSolver
        // The PCM potential must be added to the Fock operator *after* the
        // energy computation, not in form_F()
        if(pcm_enabled_) {
          // Prepare the density
          SharedMatrix D_pcm(Da_->clone());
          if(same_a_b_orbs()) {
            D_pcm->scale(2.0); // PSI4's density doesn't include the occupation
          }
          else {
            D_pcm->add(Db_);
          }

          // Compute the PCM charges and polarization energy
          double Epcm = 0.0;
      if (options_.get_str("PCM_SCF_TYPE") == "TOTAL")
      {
            Epcm = hf_pcm_->compute_E(D_pcm, PCM::Total);
      }
      else
      {
            Epcm = hf_pcm_->compute_E(D_pcm, PCM::NucAndEle);
      }
          energies_["PCM Polarization"] = Epcm;
      Process::environment.globals["PCM POLARIZATION ENERGY"] = Epcm;
          E_ += Epcm;

          // Add the PCM potential to the Fock matrix
          SharedMatrix V_pcm;
          V_pcm = hf_pcm_->compute_V();
          if(same_a_b_orbs()) Fa_->add(V_pcm);
          else {
            Fa_->add(V_pcm);
            Fb_->add(V_pcm);
          }
        }
#endif
        std::string status = "";

        // We either do SOSCF or DIIS
        bool did_soscf = false;
        if (soscf_enabled_ && (Drms_ < soscf_r_start_) && (iteration_ > 3)){
            compute_orbital_gradient(false);
            diis_performed_ = false;

            if (!test_convergency()){
                int nmicro = soscf_update();
                if (nmicro > 0){ // If zero the soscf call bounced for some reason
                    find_occupation();
                    status += "SOSCF, nmicro = ";
                    status += psi::to_string(nmicro);
                    did_soscf = true; // Stops DIIS
                }
                else{
                    if (print_){
                        outfile->Printf("Did not take a SOSCF step, using normal convergence methods\n");
                    }
                    did_soscf = false; // Back to DIIS
                }
            }
            else{
                // We need to ensure orthogonal orbitals and set epsilon
                status += "SOSCF, conv";
                timer_on("HF: Form C");
                form_C();
                timer_off("HF: Form C");
                did_soscf = true; // Stops DIIS
            }
        } // End SOSCF block

        if (!did_soscf){ // Normal convergence procedures if we do not do SOSCF

            timer_on("HF: DIIS");
            bool add_to_diis_subspace = false;
            if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ )
                add_to_diis_subspace = true;

            compute_orbital_gradient(add_to_diis_subspace);

            if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
                diis_performed_ = diis();
            } else {
                diis_performed_ = false;
            }
            timer_off("HF: DIIS");

            if (print_>4 && diis_performed_) {
                outfile->Printf("  After DIIS:\n");
                Fa_->print("outfile");
                Fb_->print("outfile");
            }

            timer_on("HF: Form C");
            form_C();
            timer_off("HF: Form C");
        }

        // If we're too well converged, or damping wasn't enabled, do DIIS
        damping_performed_ = (damping_enabled_ && iteration_ > 1 && Drms_ > damping_convergence_);

        if(diis_performed_){
            if(status != "") status += "/";
            status += "DIIS";
        }
        if(MOM_performed_){
            if(status != "") status += "/";
            status += "MOM";
        }
        if(damping_performed_){
            if(status != "") status += "/";
            status += "DAMP";
        }
        if(frac_performed_){
            if(status != "") status += "/";
            status += "FRAC";
        }

        timer_on("HF: Form D");
        form_D();
        timer_off("HF: Form D");

        Process::environment.globals["SCF ITERATION ENERGY"] = E_;

        // After we've built the new D, damp the update if
        if(damping_performed_) damp_update();

        if (print_ > 3){
            Ca_->print("outfile");
            Cb_->print("outfile");
            Da_->print("outfile");
            Db_->print("outfile");
        }

        converged_ = test_convergency();

        df = (options_.get_str("SCF_TYPE") == "DF");


        outfile->Printf( "   @%s%s iter %3d: %20.14f   %12.5e   %-11.5e %s\n", df ? "DF-" : "",
                          reference.c_str(), iteration_, E_, E_ - Eold_, Drms_, status.c_str());


        // If a an excited MOM is requested but not started, don't stop yet
        if (MOM_excited_ && !MOM_started_) converged_ = false;

        // If a fractional occupation is requested but not started, don't stop yet
        if (frac_enabled_ && !frac_performed_) converged_ = false;

        // If a DF Guess environment, reset the JK object, and keep running
        if (converged_ && options_.get_bool("DF_SCF_GUESS") && (old_scf_type_ == "DIRECT")) {
            outfile->Printf( "\n  DF guess converged.\n\n"); // Be cool dude.
            converged_ = false;
            if(initialized_diis_manager_)
                diis_manager_->reset_subspace();
            scf_type_ = old_scf_type_;
            options_.set_str("SCF","SCF_TYPE",old_scf_type_);
            old_scf_type_ = "DF";
            integrals();
        }

        // Call any postiteration callbacks
//        call_postiteration_callbacks();

    } while (!converged_ && iteration_ < maxiter_ );

}

void HF::print_energies()
{
    outfile->Printf("   => Energetics <=\n\n");
    outfile->Printf("    Nuclear Repulsion Energy =        %24.16f\n", energies_["Nuclear"]);
    outfile->Printf("    One-Electron Energy =             %24.16f\n", energies_["One-Electron"]);
    outfile->Printf("    Two-Electron Energy =             %24.16f\n", energies_["Two-Electron"]);
    outfile->Printf("    DFT Exchange-Correlation Energy = %24.16f\n", energies_["XC"]);
    outfile->Printf("    Empirical Dispersion Energy =     %24.16f\n", energies_["-D"]);
    if (!pcm_enabled_)
        energies_["PCM Polarization"] = 0.0;
    outfile->Printf("    PCM Polarization Energy =         %24.16f\n", energies_["PCM Polarization"]);
    outfile->Printf("    EFP Energy =                      %24.16f\n", energies_["EFP"]);
    outfile->Printf("    Total Energy =                    %24.16f\n", energies_["Nuclear"] +
        energies_["One-Electron"] + energies_["Two-Electron"] + energies_["XC"] +
        energies_["-D"] + energies_["EFP"] + energies_["PCM Polarization"]);
    outfile->Printf( "\n");

    Process::environment.globals["NUCLEAR REPULSION ENERGY"] = energies_["Nuclear"];
    Process::environment.globals["ONE-ELECTRON ENERGY"] = energies_["One-Electron"];
    Process::environment.globals["TWO-ELECTRON ENERGY"] = energies_["Two-Electron"];
    if (fabs(energies_["XC"]) > 1.0e-14) {
        Process::environment.globals["DFT XC ENERGY"] = energies_["XC"];
        Process::environment.globals["DFT FUNCTIONAL TOTAL ENERGY"] = energies_["Nuclear"] +
            energies_["One-Electron"] + energies_["Two-Electron"] + energies_["XC"];
        Process::environment.globals["DFT TOTAL ENERGY"] = energies_["Nuclear"] +
            energies_["One-Electron"] + energies_["Two-Electron"] + energies_["XC"] + energies_["-D"];
    } else {
        Process::environment.globals["HF TOTAL ENERGY"] = energies_["Nuclear"] +
            energies_["One-Electron"] + energies_["Two-Electron"];
    }
    if (fabs(energies_["-D"]) > 1.0e-14) {
        Process::environment.globals["DISPERSION CORRECTION ENERGY"] = energies_["-D"];
    }

    Process::environment.globals["SCF ITERATIONS"] = iteration_;

    // Only print this alert if we are actually doing EFP or PCM
    if(pcm_enabled_ || ( Process::environment.get_efp()->get_frag_count() > 0 ) ) {
        outfile->Printf("    Alert: EFP and PCM quantities not currently incorporated into SCF psivars.");
    }

//  Comment so that autodoc utility will find this PSI variable
//     It doesn't really belong here but needs to be linked somewhere
//  Process::environment.globals["DOUBLE-HYBRID CORRECTION ENERGY"]
}

void HF::print_occupation()
{

        char **labels = molecule_->irrep_labels();
        std::string reference = options_.get_str("REFERENCE");
        outfile->Printf( "          ");
        for(int h = 0; h < nirrep_; ++h) outfile->Printf( " %4s ", labels[h]); outfile->Printf( "\n");
        outfile->Printf( "    DOCC [ ");
        for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", doccpi_[h]);
        outfile->Printf( " %4d ]\n", doccpi_[nirrep_-1]);
        if(reference != "RHF" && reference != "RKS"){
            outfile->Printf( "    SOCC [ ");
            for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", soccpi_[h]);
            outfile->Printf( " %4d ]\n", soccpi_[nirrep_-1]);
        }
        if (MOM_excited_) {
            // Also print nalpha and nbeta per irrep, which are more physically meaningful
            outfile->Printf( "    NA   [ ");
            for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", nalphapi_[h]);
            outfile->Printf( " %4d ]\n", nalphapi_[nirrep_-1]);
            outfile->Printf( "    NB   [ ");
            for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", nbetapi_[h]);
            outfile->Printf( " %4d ]\n", nbetapi_[nirrep_-1]);
        }

        for(int h = 0; h < nirrep_; ++h) free(labels[h]); free(labels);
        outfile->Printf("\n");

}

//  Returns a vector of the occupation of the a orbitals
std::shared_ptr<Vector> HF::occupation_a() const
{
  SharedVector occA = SharedVector(new Vector(nmopi_));
  for(int h=0; h < nirrep_;++h)
    for(int n=0; n < nalphapi()[h]; n++)
      occA->set(h, n, 1.0);

  return occA;
}

//  Returns a vector of the occupation of the b orbitals
std::shared_ptr<Vector> HF::occupation_b() const
{
  SharedVector occB = SharedVector(new Vector(nmopi_));
  for(int h=0; h < nirrep_;++h)
    for(int n=0; n < nbetapi()[h]; n++)
      occB->set(h, n, 1.0);

  return occB;
}

void HF::diagonalize_F(const SharedMatrix& Fm, SharedMatrix& Cm, std::shared_ptr<Vector>& epsm)
{
    //Form F' = X'FX for canonical orthogonalization
    diag_temp_->gemm(true, false, 1.0, X_, Fm, 0.0);
    diag_F_temp_->gemm(false, false, 1.0, diag_temp_, X_, 0.0);

    //Form C' = eig(F')
    diag_F_temp_->diagonalize(diag_C_temp_, epsm);

    //Form C = XC'
    Cm->gemm(false, false, 1.0, X_, diag_C_temp_, 0.0);
}

void HF::reset_occupation()
{
    // RHF style for now
    doccpi_ = original_doccpi_;
    soccpi_ = original_soccpi_;
    nalphapi_ = doccpi_ + soccpi_;
    nbetapi_ = doccpi_;

    // These may not match the per irrep. Will remap correctly next find_occupation call
    nalpha_ = original_nalpha_;
    nbeta_ = original_nbeta_;

}
SharedMatrix HF::form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi)
{
    const int* nsopi = Cso->rowspi();
    const int* nmopi = Cso->colspi();
    int* nvirpi = new int[nirrep_];

    for (int h = 0; h < nirrep_; h++)
        nvirpi[h] = nmopi[h] - noccpi[h];

    SharedMatrix Fia(new Matrix("Fia (Some Basis)", nirrep_, noccpi, nvirpi));

    // Hack to get orbital e for this Fock
    SharedMatrix C2(new Matrix("C2", Cso->rowspi(), Cso->colspi()));
    std::shared_ptr<Vector> E2(new Vector("E2", Cso->colspi()));
    diagonalize_F(Fso, C2, E2);

    for (int h = 0; h < nirrep_; h++) {
        int nmo = nmopi[h];
        int nso = nsopi[h];
        int nvir = nvirpi[h];
        int nocc = noccpi[h];

        if (nmo == 0 || nso == 0 || nvir == 0 || nocc == 0) continue;

        //double** C = Cso->pointer(h);
        double** C = C2->pointer(h);
        double** F = Fso->pointer(h);
        double** Fiap = Fia->pointer(h);

        double** Temp = block_matrix(nocc, nso);

        C_DGEMM('T','N',nocc,nso,nso,1.0,C[0],nmo,F[0],nso,0.0,Temp[0],nso);
        C_DGEMM('N','N',nocc,nvir,nso,1.0,Temp[0],nso,&C[0][nocc],nmo,0.0,Fiap[0],nvir);

        free_block(Temp);

        //double* eps = E2->pointer(h);
        //for (int i = 0; i < nocc; i++)
        //    for (int a = 0; a < nvir; a++)
        //        Fiap[i][a] /= eps[a + nocc] - eps[i];

    }

    //Fia->print();

    delete[] nvirpi;

    return Fia;
}
SharedMatrix HF::form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso)
{
    SharedMatrix FDSmSDF(new Matrix("FDS-SDF", nirrep_, nsopi_, nsopi_));
    SharedMatrix DS(new Matrix("DS", nirrep_, nsopi_, nsopi_));

    DS->gemm(false,false,1.0,Dso,S_,0.0);
    FDSmSDF->gemm(false,false,1.0,Fso,DS,0.0);

    SharedMatrix SDF(FDSmSDF->transpose());
    FDSmSDF->subtract(SDF);

    DS.reset();
    SDF.reset();

    SharedMatrix XP(new Matrix("X'(FDS - SDF)", nirrep_, nmopi_, nsopi_));
    SharedMatrix XPX(new Matrix("X'(FDS - SDF)X", nirrep_, nmopi_, nmopi_));
    XP->gemm(true,false,1.0,X_,FDSmSDF,0.0);
    XPX->gemm(false,false,1.0,XP,X_,0.0);

    //XPX->print();

    return XPX;
}

void HF::print_stability_analysis(std::vector<std::pair<double, int> > &vec)
{
    std::sort(vec.begin(), vec.end());
    std::vector<std::pair<double, int> >::const_iterator iter = vec.begin();
    outfile->Printf( "    ");
    char** irrep_labels = molecule_->irrep_labels();
    int count = 0;
    for(; iter != vec.end(); ++iter){
        ++count;
        outfile->Printf( "%4s %-10.6f", irrep_labels[iter->second], iter->first);
        if(count == 4){
            outfile->Printf( "\n    ");
            count = 0;
        }else{
            outfile->Printf( "    ");
        }
    }
    if(count)
        outfile->Printf( "\n\n");
    else
        outfile->Printf( "\n");

    for(int h = 0; h < nirrep_; ++h)
        free(irrep_labels[h]);
    free(irrep_labels);
}
bool HF::stability_analysis()
{
    throw PSIEXCEPTION("Stability analysis hasn't been implemented yet for this wfn type.");
    return false;
}
}}
