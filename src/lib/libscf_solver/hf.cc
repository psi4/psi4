/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
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

#include <libmints/mints.h>

#include <libfunctional/superfunctional.h>
#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <liboptions/liboptions_python.h>
#include <psifiles.h>
#include <libfock/jk.h>

#include "hf.h"

#include <psi4-dec.h>
#include <libefp_solver/efp_solver.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

HF::HF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6)
{
    common_init();
}

HF::HF(Options& options, boost::shared_ptr<PSIO> psio)
    : Wavefunction(options, psio),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6)
{
    common_init();
}

HF::~HF()
{
}

void HF::common_init()
{

    attempt_number_ = 1;

    // This quantity is needed fairly soon
    nirrep_ = factory_->nirrep();

    integral_threshold_ = options_.get_double("INTS_TOLERANCE");

    scf_type_ = options_.get_str("SCF_TYPE");

    H_.reset(factory_->create_matrix("One-electron Hamiltonion"));
    X_.reset(factory_->create_matrix("X"));

    nmo_ = 0;
    nso_ = 0;
    int* dimpi = factory_->colspi();
    for (int h = 0; h< factory_->nirrep(); h++){
        nsopi_[h] = dimpi[h];
        nmopi_[h] = nsopi_[h]; //For now
        nso_ += nsopi_[h];
        nmo_ += nmopi_[h]; //For now
    }

    Eold_    = 0.0;
    E_       = 0.0;
    maxiter_ = 40;

    // Read information from input file
    maxiter_ = options_.get_int("MAXITER");

    // Should we continue if we fail to converge?
    fail_on_maxiter_ = options_.get_bool("FAIL_ON_MAXITER");

    // Read in DOCC and SOCC from memory
    int nirreps = factory_->nirrep();
    int ndocc = 0, nsocc = 0;
    input_docc_ = false;
    if (options_["DOCC"].has_changed()) {
        input_docc_ = true;
        // Map the symmetry of the input DOCC, to account for displacements
        boost::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
        if(old_pg){
            // This is one of a series of displacements;  check the dimension against the parent point group
            int full_nirreps = old_pg->char_table().nirrep();
            if(options_["DOCC"].size() != full_nirreps)
                throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            int *temp_docc = new int[full_nirreps];
            for(int h = 0; h < full_nirreps; ++h)
                temp_docc[h] = options_["DOCC"][h].to_integer();
            map_irreps(temp_docc);
            doccpi_ = temp_docc;
            delete[] temp_docc;
        }else{
            // This is a normal calculation; check the dimension against the current point group then read
            if(options_["DOCC"].size() != nirreps)
                throw PSIEXCEPTION("Input DOCC array has the wrong dimensions");
            for(int h = 0; h < nirreps; ++h)
                doccpi_[h] = options_["DOCC"][h].to_integer();
        }
        for (int i=0; i<nirreps; ++i)
            ndocc += 2*doccpi_[i];
    } else {
        for (int i=0; i<nirreps; ++i)
            doccpi_[i] = 0;
    }

    if(options_.get_str("REFERENCE") == "RKS" || options_.get_str("REFERENCE") == "UKS")
        name_ = "DFT";
    else
        name_ = "SCF";

    input_socc_ = false;
    if (options_["SOCC"].has_changed()) {
        input_socc_ = true;
        // Map the symmetry of the input SOCC, to account for displacements
        boost::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
        if(old_pg){
            // This is one of a series of displacements;  check the dimension against the parent point group
            int full_nirreps = old_pg->char_table().nirrep();
            if(options_["SOCC"].size() != full_nirreps)
                throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            int *temp_socc = new int[full_nirreps];
            for(int h = 0; h < full_nirreps; ++h)
                temp_socc[h] = options_["SOCC"][h].to_integer();
            map_irreps(temp_socc);
            soccpi_ = temp_socc;
            delete[] temp_socc;
        }else{
            // This is a normal calculation; check the dimension against the current point group then read
            if(options_["SOCC"].size() != nirreps)
                throw PSIEXCEPTION("Input SOCC array has the wrong dimensions");
            for(int h = 0; h < nirreps; ++h)
                soccpi_[h] = options_["SOCC"][h].to_integer();
        }
        for (int i=0; i<nirreps; ++i)
            nsocc += soccpi_[i];
    } else {
        for (int i=0; i<nirreps; ++i)
            soccpi_[i] = 0;
    }



    // Read information from checkpoint
    nuclearrep_ = molecule_->nuclear_repulsion_energy();

    // Determine the number of electrons in the system
    charge_ = molecule_->molecular_charge();
    nelectron_  = 0;
    for (int i=0; i<molecule_->natom(); ++i)
        nelectron_ += (int)molecule_->Z(i);
    nelectron_ -= charge_;

    // If the user told us the multiplicity, read it from the input
    if(molecule_->multiplicity_specified()){
        multiplicity_ = molecule_->multiplicity();
    }else{
        if(nelectron_%2){
            multiplicity_ = 2;
            molecule_->set_multiplicity(2);

            // There are an odd number of electrons
                outfile->Printf("\tThere are an odd number of electrons - assuming doublet.\n"
                            "\tSpecify the multiplicity with the MULTP option in the\n"
                            "\tinput if this is incorrect\n\n");

        }else{
            multiplicity_ = 1;
            // There are an even number of electrons

                outfile->Printf("\tThere are an even number of electrons - assuming singlet.\n"
                            "\tSpecify the multiplicity with the MULTP option in the\n"
                            "\tinput if this is incorrect\n\n");

        }
    }

    // Make sure that the multiplicity is reasonable
    if(multiplicity_ - 1 > nelectron_){
        char *str = new char[100];
        sprintf(str, "There are not enough electrons for multiplicity = %d, \n"
                     "please check your input and use the MULTP keyword", multiplicity_);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }
    if(multiplicity_%2 == nelectron_%2){
        char *str = new char[100];
        sprintf(str, "A multiplicity of %d with %d electrons is impossible.\n"
                     "Please check your input and use the MULTP and/or CHARGE keywords",
                     multiplicity_, nelectron_);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }

    nbeta_  = (nelectron_ - multiplicity_ + 1)/2;
    nalpha_ = nbeta_ + multiplicity_ - 1;

    if (input_socc_ || input_docc_) {
        for (int h = 0; h < nirrep_; h++) {
            nalphapi_[h] = doccpi_[h] + soccpi_[h];
            nbetapi_[h]  = doccpi_[h];
        }
    }

    // Check if it is a broken symmetry solution
    broken_symmetry_ = false;
    int socc = 0;
    for (int h = 0; h < nirrep_; h++) {
        socc += soccpi_[h];
    }

    if (multiplicity_ == 1 && socc == 2) {
        // Set up occupation for the broken symmetry solution
        outfile->Printf( "  Broken symmetry solution detected... \n"); //TEST
        broken_symmetry_ = true;
        int socc_count = 0;
        nalphapi_ = doccpi_;
        nbetapi_  = doccpi_;

        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i<soccpi_[h]; i++) {
                socc_count++;
                if (socc_count == 1) {
                    nalphapi_[h]++;
                }
                if (socc_count == 2) {
                    nbetapi_[h]++;
                }
            }
        }
        if (print_ > 2) {
            nalphapi_.print();
            nbetapi_.print();
        }
    }

    perturb_h_ = false;
    perturb_h_ = options_.get_bool("PERTURB_H");
    perturb_ = nothing;
    lambda_ = 0.0;
    if (perturb_h_) {
        string perturb_with;

        lambda_ = options_.get_double("PERTURB_MAGNITUDE");

        if (options_["PERTURB_WITH"].has_changed()) {
            perturb_with = options_.get_str("PERTURB_WITH");
            // Do checks to see what perturb_with is.
            if (perturb_with == "DIPOLE_X")
                perturb_ = dipole_x;
            else if (perturb_with == "DIPOLE_Y")
                perturb_ = dipole_y;
            else if (perturb_with == "DIPOLE_Z")
                perturb_ = dipole_z;
            else if (perturb_with == "EMBPOT") {
                perturb_ = embpot;
                lambda_ = 1.0;
            }
            else if (perturb_with == "DX") {
                perturb_ = dx;
                lambda_ = 1.0;
            }
            else if (perturb_with == "SPHERE") {
                perturb_ = sphere;
                lambda_ = 1.0;
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

    MOM_enabled_ = (options_.get_int("MOM_START") != 0);
    MOM_excited_ = (options_["MOM_OCC"].size() != 0 && MOM_enabled_);
    MOM_started_ = false;
    MOM_performed_ = false;

    frac_enabled_ = (options_.get_int("FRAC_START") != 0);
    frac_performed_ = false;

    print_header();
}

void HF::damp_update()
{
    throw PSIEXCEPTION("Sorry, damping has not been implemented for this "
                       "type of SCF wavefunction yet.");
}

void HF::integrals()
{
    if (print_ )
        outfile->Printf( "  ==> Integral Setup <==\n\n");

    // Build the JK from options, symmetric type
    try {
        jk_ = JK::build_JK();
    }
    catch(const BasisSetNotFound& e) {
        if (options_.get_str("SCF_TYPE") == "DF" || options_.get_int("DF_SCF_GUESS") == 1) {
            outfile->Printf( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            outfile->Printf( "%s\n", e.what());
            outfile->Printf( "   Turning off DF and switching to PK method.\n");
            outfile->Printf( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
            options_.set_str("SCF", "SCF_TYPE", "PK");
            options_.set_bool("SCF", "DF_SCF_GUESS", false);
            jk_ = JK::build_JK();
        }
        else
            throw; // rethrow the error
    }

    // Tell the JK to print
    jk_->set_print(print_);
    // Give the JK 75% of the memory
    jk_->set_memory((ULI)(options_.get_double("SCF_MEM_SAFETY_FACTOR")*(Process::environment.get_memory() / 8L)));

    // DFT sometimes needs custom stuff
    if ((options_.get_str("REFERENCE") == "UKS" || options_.get_str("REFERENCE") == "RKS")) {

        // Need a temporary functional
        boost::shared_ptr<SuperFunctional> functional =
            SuperFunctional::current(options_);

        // K matrices
        jk_->set_do_K(functional->is_x_hybrid());
        // wK matrices
        jk_->set_do_wK(functional->is_x_lrc());
        // w Value
        jk_->set_omega(functional->x_omega());
    }

    // Initialize
    jk_->initialize();
    // Print the header
    jk_->print_header();
}

void HF::finalize()
{
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

    dump_to_checkpoint();

    //Sphalf_.reset();
    X_.reset();
    T_.reset();
    V_.reset();
    diag_temp_.reset();
    diag_F_temp_.reset();
    diag_C_temp_.reset();

    // Close the chkpt
    if(psio_->open_check(PSIF_CHKPT))
        psio_->close(PSIF_CHKPT, 1);
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
                pairs_a.push_back(make_pair(epsilon_a_->get(h, i), h));
        }
        for (int h=0; h<epsilon_b_->nirrep(); ++h) {
            for (int i=0; i<epsilon_b_->dimpi()[h]; ++i)
                pairs_b.push_back(make_pair(epsilon_b_->get(h, i), h));
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

                outfile->Printf( "\tOccupation by irrep:\n");
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
        nthread = omp_get_max_threads();
    #endif


        outfile->Printf( "\n");
        outfile->Printf( "         ---------------------------------------------------------\n");
        outfile->Printf( "                                   SCF\n");
        outfile->Printf( "            by Justin Turney, Rob Parrish, and Andy Simmonett\n");
        outfile->Printf( "                             %4s Reference\n", options_.get_str("REFERENCE").c_str());
        outfile->Printf( "                      %3d Threads, %6ld MiB Core\n", nthread, memory_ / 1000000L);
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
            outfile->Printf( "     %-3s   %6d  %6d  %6d  %6d  %6d  %6d\n", ct.gamma(h).symbol(), nsopi_[h], nmopi_[h], nalphapi_[h], nbetapi_[h], doccpi_[h], soccpi_[h]);
        }
        outfile->Printf( "   -------------------------------------------------------\n");
        outfile->Printf( "    Total  %6d  %6d  %6d  %6d  %6d  %6d\n", nso_, nmo_, nalpha_, nbeta_, nbeta_, nalpha_-nbeta_);
        outfile->Printf( "   -------------------------------------------------------\n\n");

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
        MintsHelper helper(options_, 0);
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
          fscanf(input, "%d", &npoints);
          outfile->Printf( "  npoints = %d\n", npoints);
          double x, y, z, w, v;
          double max = 0;
          for(int k=0; k < npoints; k++) {
            fscanf(input, "%lf %lf %lf %lf %lf", &x, &y, &z, &w, &v);
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


          outfile->Printf( "  Perturbing H by %f V_eff.\n", lambda_);
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
      else {
          // The following perturbations are handled by MintsHelper.
#if 0
        OperatorSymmetry msymm(1, molecule_, integral_, factory_);
        vector<SharedMatrix> dipoles = msymm.create_matrices("Dipole");
        OneBodySOInt *so_dipole = integral_->so_dipole();
        so_dipole->compute(dipoles);

        if (perturb_ == dipole_x ) {
            if (msymm.component_symmetry(0) != 0){
                outfile->Printf( "  WARNING: You requested mu(x) perturbation, but mu(x) is not symmetric.\n");
            }
            else {
                    outfile->Printf( "  Perturbing H by %f mu(x).\n", lambda_);
                dipoles[0]->scale(lambda_);
                V_->add(dipoles[0]);
            }
        } else if (perturb_ == dipole_y) {
            if (msymm.component_symmetry(1) != 0){
                    outfile->Printf( "  WARNING: You requested mu(y) perturbation, but mu(y) is not symmetric.\n");
            }
            else {
                    outfile->Printf( "  Perturbing H by %f mu(y).\n", lambda_);
                dipoles[1]->scale(lambda_);
                V_->add(dipoles[1]);
            }
        } else if (perturb_ == dipole_z) {
            if (msymm.component_symmetry(2) != 0){
                    outfile->Printf( "  WARNING: You requested mu(z) perturbation, but mu(z) is not symmetric.\n");
            }
            else {
                    outfile->Printf( "  Perturbing H by %f mu(z).\n", lambda_);
                dipoles[2]->scale(lambda_);
                V_->add(dipoles[2]);
            }
        }
#endif
      } // end dipole perturbations
    } // end perturb_h_

    // If an external field exists, add it to the one-electron Hamiltonian
    boost::python::object pyExtern = dynamic_cast<PythonDataType*>(options_["EXTERN"].get())->to_python();
    boost::shared_ptr<ExternalPotential> external = boost::python::extract<boost::shared_ptr<ExternalPotential> >(pyExtern);
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
    int *dimpi = eigval->dimpi();
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
            outfile->Printf("  Using Symmetric Orthogonalization.\n");

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

        // Refreshes twice in RHF, no big deal
        epsilon_a_->init(nmopi_);
        Ca_->init(nirrep_,nsopi_,nmopi_,"MO coefficients");
        epsilon_b_->init(nmopi_);
        Cb_->init(nirrep_,nsopi_,nmopi_,"MO coefficients");
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
                pairs.push_back(make_pair(epsilon_a_->get(h, i), h));
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
                pairs.push_back(make_pair(epsilon_a_->get(h, i), h));
            frzvpi_[h] = 0;
        }
        sort(pairs.begin(),pairs.end(), greater<std::pair<double, int> >());

        for (int i=0; i<nfzv; ++i)
            frzvpi_[pairs[i].second]++;
    }
}

void HF::print_orbitals(const char* header, std::vector<std::pair<double, std::pair<const char*, int> > > orbs)
{

        outfile->Printf( "\t%-70s\n\n\t", header);
        int count = 0;
        for (int i = 0; i < orbs.size(); i++) {
            outfile->Printf( "%4d%-4s%11.6f  ", orbs[i].second.second, orbs[i].second.first, orbs[i].first);
            if (count++ % 3 == 2 && count != orbs.size())
                outfile->Printf( "\n\t");
        }
        outfile->Printf( "\n\n");

}

void HF::print_orbitals()
{
    char **labels = molecule_->irrep_labels();

        outfile->Printf( "\tOrbital Energies (a.u.)\n\t-----------------------\n\n");

    std::string reference = options_.get_str("REFERENCE");
    if((reference == "RHF") || (reference == "RKS")){

        std::vector<std::pair<double, std::pair<const char*, int> > > occ;
        std::vector<std::pair<double, std::pair<const char*, int> > > vir;

        for (int h = 0; h < nirrep_; h++) {

            std::vector<std::pair<double, int> > orb_e;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_e.push_back(make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occ.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                vir.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));

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
                orb_eA.push_back(make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_eA.begin(), orb_eA.end());

            std::vector<int> orb_orderA(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_orderA[orb_eA[a].second] = a;

            for (int a = 0; a < nalphapi_[h]; a++)
                occA.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_orderA[a] + 1)));
            for (int a = nalphapi_[h]; a < nmopi_[h]; a++)
                virA.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_orderA[a] + 1)));

            std::vector<std::pair<double, int> > orb_eB;
            for (int a = 0; a < nmopi_[h]; a++)
                orb_eB.push_back(make_pair(epsilon_b_->get(h,a), a));
            std::sort(orb_eB.begin(), orb_eB.end());

            std::vector<int> orb_orderB(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_orderB[orb_eB[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                occB.push_back(make_pair(epsilon_b_->get(h,a), make_pair(labels[h],orb_orderB[a] + 1)));
            for (int a = nbetapi_[h]; a < nmopi_[h]; a++)
                virB.push_back(make_pair(epsilon_b_->get(h,a), make_pair(labels[h],orb_orderB[a] + 1)));

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
                orb_e.push_back(make_pair(epsilon_a_->get(h,a), a));
            std::sort(orb_e.begin(), orb_e.end());

            std::vector<int> orb_order(nmopi_[h]);
            for (int a = 0; a < nmopi_[h]; a++)
                orb_order[orb_e[a].second] = a;

            for (int a = 0; a < nbetapi_[h]; a++)
                docc.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nbetapi_[h] ; a < nalphapi_[h]; a++)
                socc.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));
            for (int a = nalphapi_[h] ; a < nmopi_[h]; a++)
                vir.push_back(make_pair(epsilon_a_->get(h,a), make_pair(labels[h],orb_order[a] + 1)));

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


        outfile->Printf( "\tFinal Occupation by Irrep:\n");
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
    // "READ"-try to read MOs from guess file, projecting if needed
    // "CORE"-CORE Hamiltonain
    // "GWH"-Generalized Wolfsberg-Helmholtz
    // "SAD"-Superposition of Atomic Denisties
    string guess_type = options_.get_str("GUESS");
    if (guess_type == "READ" && !psio_->exists(PSIF_SCF_MOS)) {
        outfile->Printf( "  SCF Guess was Projection but file not found.\n");
        outfile->Printf( "  Switching over to SAD guess.\n\n");
        guess_type = "SAD";
    }

    if (guess_type == "READ") {
        if (do_use_fock_guess()) {
            outfile->Printf( "  SCF Guess: Guess MOs from previously saved Fock matrix.\n\n");
            load_fock(); // won't save the energy from here
            form_C();
            form_D();
        }
        else {
            outfile->Printf( "  SCF Guess: Reading in previously saved MOs, projecting if necessary.\n\n");
            load_orbitals(); // won't save the energy from here
            form_D();
        }

    } else if (guess_type == "SAD") {

        if (print_)
            outfile->Printf( "  SCF Guess: Superposition of Atomic Densities via on-the-fly atomic UHF.\n\n");

        //Superposition of Atomic Density (RHF only at present)
        compute_SAD_guess();
        guess_E = compute_initial_E();

    } else if (guess_type == "GWH") {
        //Generalized Wolfsberg Helmholtz (Sounds cool, easy to code)
        if (print_)
            outfile->Printf( "  SCF Guess: Generalized Wolfsberg-Helmholtz.\n\n");

        Fa_->zero(); //Try Fa_{mn} = S_{mn} (H_{mm} + H_{nn})/2
        int h, i, j;
        int *opi = S_->rowspi();
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

void HF::save_orbitals()
{
    psio_->open(PSIF_SCF_MOS,PSIO_OPEN_NEW);

    if (print_)
        outfile->Printf("\n  Saving occupied orbitals to File %d.\n", PSIF_SCF_MOS);

    psio_->write_entry(PSIF_SCF_MOS,"SCF ENERGY",(char *) &(E_),sizeof(double));
    psio_->write_entry(PSIF_SCF_MOS,"NIRREP",(char *) &(nirrep_),sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"NSOPI",(char *) &(nsopi_[0]),nirrep_*sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"NALPHAPI",(char *) &(nalphapi_[0]),nirrep_*sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"NBETAPI",(char *) &(nbetapi_[0]),nirrep_*sizeof(int));

    char *basisname = strdup(options_.get_str("BASIS").c_str());
    int basislength = strlen(options_.get_str("BASIS").c_str()) + 1;
    int nbf = basisset_->nbf();

    psio_->write_entry(PSIF_SCF_MOS,"BASIS NAME LENGTH",(char *)(&basislength),sizeof(int));
    psio_->write_entry(PSIF_SCF_MOS,"BASIS NAME",basisname,basislength*sizeof(char));
    psio_->write_entry(PSIF_SCF_MOS,"NUMBER OF BASIS FUNCTIONS",(char *)(&nbf),sizeof(int));

    // upon loading, need to know what value of puream was used
    int old_puream = (basisset_->has_puream() ? 1 : 0);
    psio_->write_entry(PSIF_SCF_MOS,"PUREAM",(char *)(&old_puream),sizeof(int));

    SharedMatrix Ctemp_a(new Matrix("ALPHA MOS", nirrep_, nsopi_, nalphapi_));
    for (int h = 0; h < nirrep_; h++)
        for (int m = 0; m<nsopi_[h]; m++)
            for (int i = 0; i<nalphapi_[h]; i++)
                Ctemp_a->set(h,m,i,Ca_->get(h,m,i));
    Ctemp_a->save(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);

    SharedMatrix Ctemp_b(new Matrix("BETA MOS", nirrep_, nsopi_, nbetapi_));
    for (int h = 0; h < nirrep_; h++)
        for (int m = 0; m<nsopi_[h]; m++)
            for (int i = 0; i<nbetapi_[h]; i++)
                Ctemp_b->set(h,m,i,Cb_->get(h,m,i));
    Ctemp_b->save(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
    // Write Fock matrix to file 280 after removing symmetry
    MintsHelper helper(options_, 0);
    SharedMatrix sotoao = helper.petite_list()->sotoao();
    SharedMatrix Fa(new Matrix(nbf,nbf));
    SharedMatrix Fb(new Matrix(nbf,nbf));
    Fa->remove_symmetry(Fa_,sotoao);
    Fb->remove_symmetry(Fb_,sotoao);
    Fa->set_name("ALPHA FOCK C1");
    Fb->set_name("BETA FOCK C1");
    Fa->save(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
    Fb->save(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);

    psio_->close(PSIF_SCF_MOS,1);
    free(basisname);
}

bool HF::do_use_fock_guess() // only use this approach when symmetry changes are causing issues
{
  psio_->open(PSIF_SCF_MOS,PSIO_OPEN_OLD);
  // Compare current basis with old basis to see if fock guess will work
  bool is_same_basis = false;
  int nbf = basisset_->nbf(), basislength, old_nbf;
  psio_->read_entry(PSIF_SCF_MOS,"BASIS NAME LENGTH",(char *)(&basislength),sizeof(int));
  char *basisnamec = new char[basislength];
  psio_->read_entry(PSIF_SCF_MOS,"BASIS NAME",basisnamec,basislength*sizeof(char));
  psio_->read_entry(PSIF_SCF_MOS,"NUMBER OF BASIS FUNCTIONS",(char *)(&old_nbf),sizeof(int));
  std::string old_basisname(basisnamec); delete[] basisnamec;
  is_same_basis = (options_.get_str("BASIS") == old_basisname && nbf == old_nbf);
  // Compare number of irreps with old number to see if fock guess is necessary
  bool is_different_symmetry = false;
  int old_nirrep;
  psio_->read_entry(PSIF_SCF_MOS,"NIRREP",(char *) &(old_nirrep),sizeof(int));
  is_different_symmetry = (nirrep_ != old_nirrep);
  psio_->close(PSIF_SCF_MOS,1);

  // Final comparison: Use a Fock guess if you are using the same basis as before and
  //                   the symmetry has changed
  return (is_same_basis && is_different_symmetry);
}

void HF::load_fock()
{
  int nbf = basisset_->nbf();
  psio_->open(PSIF_SCF_MOS,PSIO_OPEN_OLD);
  // Read Fock matrix from file 280, applying current symmetry
  MintsHelper helper(options_, 0);
  SharedMatrix aotoso = helper.petite_list()->aotoso();
  SharedMatrix Fa(new Matrix("ALPHA FOCK C1",nbf,nbf));
  SharedMatrix Fb(new Matrix("BETA FOCK C1",nbf,nbf));
  Fa->load(psio_,PSIF_SCF_MOS,Matrix::SubBlocks);
  Fb->load(psio_,PSIF_SCF_MOS,Matrix::SubBlocks);
  Fa_->apply_symmetry(Fa,aotoso);
  Fb_->apply_symmetry(Fb,aotoso);
  psio_->close(PSIF_SCF_MOS,1);
}

void HF::load_orbitals()
{
    psio_->open(PSIF_SCF_MOS,PSIO_OPEN_OLD);

    int basislength, old_puream;
    psio_->read_entry(PSIF_SCF_MOS,"BASIS NAME LENGTH",
        (char *)(&basislength),sizeof(int));
    char *basisnamec = new char[basislength];
    psio_->read_entry(PSIF_SCF_MOS,"BASIS NAME",basisnamec,
        basislength*sizeof(char));
    psio_->read_entry(PSIF_SCF_MOS,"PUREAM",(char *)(&old_puream),
        sizeof(int));
    bool old_forced_puream = (old_puream) ? true : false;
    std::string basisname(basisnamec);

    if (basisname == "")
        throw PSIEXCEPTION("SCF::load_orbitals: Custom basis sets not allowed for projection from a previous SCF");

    if (print_) {
        if (basisname != options_.get_str("BASIS")) {
                outfile->Printf("  Computing basis set projection from %s to %s.\n", \
                    basisname.c_str(),options_.get_str("BASIS").c_str());
        } else {
                outfile->Printf("  Using orbitals from previous SCF, no projection.\n");
        }
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(old_forced_puream));
    molecule_->set_basis_all_atoms(basisname, "DUAL_BASIS_SCF");
    boost::shared_ptr<BasisSet> dual_basis = BasisSet::construct(parser, molecule_, "DUAL_BASIS_SCF");
    // TODO: oh my, forced_puream!
    // TODO: oh my, a basis for which a fn hasn't been set in the input translation
    // TODO: oh my, a non-fitting basis to be looked up (in Mol) not under BASIS
    //boost::shared_ptr<BasisSet> dual_basis = BasisSet::pyconstruct(molecule_, basisname,
    //            "DUAL_BASIS_SCF");
    // TODO: I think Rob was planning to rework this projection bit anyways

    psio_->read_entry(PSIF_SCF_MOS,"SCF ENERGY",(char *) &(E_),sizeof(double));

    int old_nirrep, old_nsopi[8];
    psio_->read_entry(PSIF_SCF_MOS,"NIRREP",(char *) &(old_nirrep),sizeof(int));

    if (old_nirrep != nirrep_)
        throw PSIEXCEPTION("SCF::load_orbitals: Projection of orbitals between different symmetries is not currently supported");

    psio_->read_entry(PSIF_SCF_MOS,"NSOPI",(char *) (old_nsopi),nirrep_*sizeof(int));

    // Save current alpha and beta occupation vectors
    Dimension nalphapi_current (nirrep_, "Current number of alpha electrons per irrep");
    Dimension nbetapi_current (nirrep_, "Current number of beta electrons per irrep");
    nalphapi_current = nalphapi_;
    nbetapi_current = nbetapi_;

    // Read in guess alpha and beta occupation vectors
    psio_->read_entry(PSIF_SCF_MOS,"NALPHAPI",(char *) &(nalphapi_[0]),nirrep_*sizeof(int));
    psio_->read_entry(PSIF_SCF_MOS,"NBETAPI",(char *) &(nbetapi_[0]),nirrep_*sizeof(int));

    // Check if guess is the broken symmetry solution
    int nalpha_guess = 0;
    int nbeta_guess = 0;
    for (int h = 0; h < nirrep_; h++) {
        nalpha_guess += nalphapi_[h];
        nbeta_guess += nbetapi_[h];
    }

    bool guess_broken_symmetry = false;
    if (nalpha_guess == nbeta_guess && multiplicity_ == 3) guess_broken_symmetry = true;

    if (!broken_symmetry_) {
        if (!guess_broken_symmetry) {
            outfile->Printf( "  Recomputing DOCC and SOCC from number of alpha and beta electrons from previous calculation.\n");
            for (int h = 0; h < nirrep_; h++) {
                doccpi_[h] = std::min(nalphapi_[h] , nbetapi_[h]);
                soccpi_[h] = std::abs(nalphapi_[h] - nbetapi_[h]);
            }
            print_occupation();

            SharedMatrix Ctemp_a(new Matrix("ALPHA MOS", nirrep_, old_nsopi, nalphapi_));
            Ctemp_a->load(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
            SharedMatrix Ca;
            if (basisname != options_.get_str("BASIS")) {
                Ca = BasisProjection(Ctemp_a, nalphapi_, dual_basis, basisset_);
            } else {
                Ca = Ctemp_a;
            }
            for (int h = 0; h < nirrep_; h++)
                for (int m = 0; m<nsopi_[h]; m++)
                    for (int i = 0; i<nalphapi_[h]; i++)
                        Ca_->set(h,m,i,Ca->get(h,m,i));

            SharedMatrix Ctemp_b(new Matrix("BETA MOS", nirrep_, old_nsopi, nbetapi_));
            Ctemp_b->load(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
            SharedMatrix Cb;
            if (basisname != options_.get_str("BASIS")) {
                Cb = BasisProjection(Ctemp_b, nbetapi_, dual_basis, basisset_);
            } else {
                Cb = Ctemp_b;
            }
            for (int h = 0; h < nirrep_; h++)
                for (int m = 0; m<nsopi_[h]; m++)
                    for (int i = 0; i<nbetapi_[h]; i++)
                        Cb_->set(h,m,i,Cb->get(h,m,i));
        }
        // If it's a triplet and there is a broken symmetry singlet guess - ignore the guess
        else{
            // Restore nalphapi and nbetapi requested by the user
            nalphapi_ = nalphapi_current;
            nbetapi_ = nbetapi_current;
        }
    }
    // if it's a broken symmetry solution
    else {

        SharedMatrix Ctemp_a(new Matrix("ALPHA MOS", nirrep_, old_nsopi, nalphapi_));
        Ctemp_a->load(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
        SharedMatrix Ca;
        if (basisname != options_.get_str("BASIS")) {
            Ca = BasisProjection(Ctemp_a, nalphapi_, dual_basis, basisset_);
        } else {
            Ca = Ctemp_a;
        }

        SharedMatrix Ctemp_b(new Matrix("BETA MOS", nirrep_, old_nsopi, nbetapi_));
        Ctemp_b->load(psio_, PSIF_SCF_MOS, Matrix::SubBlocks);
        SharedMatrix Cb;
        if (basisname != options_.get_str("BASIS")) {
            Cb = BasisProjection(Ctemp_b, nbetapi_, dual_basis, basisset_);
        } else {
            Cb = Ctemp_b;
        }

        // Restore nalphapi and nbetapi requested by the user
        nalphapi_ = nalphapi_current;
        nbetapi_ = nbetapi_current;

        int socc_count = 0;
        for (int h = 0; h < nirrep_; h++) {
            // Copy doubly occupied orbitals into the orbital space for alpha and beta
            for (int i = 0; i<doccpi_[h]; i++) {
                for (int m = 0; m<nsopi_[h]; m++) {
                    Ca_->set(h,m,i,Ca->get(h,m,i));
                    Cb_->set(h,m,i,Cb->get(h,m,i));
                }
            }
            // Copy singly occupied orbitals into the appropriate alpha and beta orbital spaces
            for (int i = doccpi_[h]; i<(doccpi_[h]+soccpi_[h]); i++) {
                socc_count++;
                if (socc_count == 1) {
                    for (int m = 0; m<nsopi_[h]; m++) {
                        Ca_->set(h,m,doccpi_[h],Ca->get(h,m,i));
                    }
                }
                if (socc_count == 2) {
                    for (int m = 0; m<nsopi_[h]; m++) {
                        Cb_->set(h,m,doccpi_[h],Ca->get(h,m,i));
                    }
                }
            }
        }
    }
    psio_->close(PSIF_SCF_MOS,1);    
    delete[] basisnamec;
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

void HF::dump_to_checkpoint()
{
    // avoid overwriting TOC entries if nmo has changed by writing
    // an all-new checkpoint file every time!  (Will delete any
    // post-HF stuff, hopefully don't need that in a future iteration)
    // CDS 3/19/14
    /*
    if(!psio_->open_check(PSIF_CHKPT))
        psio_->open(PSIF_CHKPT, PSIO_OPEN_OLD);
    */
    if(psio_->open_check(PSIF_CHKPT)) psio_->close(PSIF_CHKPT, 0);
    psio_->open(PSIF_CHKPT, PSIO_OPEN_NEW);

    chkpt_->wt_nirreps(nirrep_);
    char **labels = molecule_->irrep_labels();
    chkpt_->wt_irr_labs(labels);
    for(int h = 0; h < nirrep_; ++h)
        free(labels[h]);
    free(labels);
    chkpt_->wt_nmo(nmo_);
    chkpt_->wt_nso(nso_);
    chkpt_->wt_nao(basisset_->nao());
    chkpt_->wt_ref(0);
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_enuc(nuclearrep_);
    chkpt_->wt_orbspi(nmopi_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_openpi(soccpi_);
    chkpt_->wt_phase_check(1); // Jet's phase check supposed to always work?
    chkpt_->wt_sopi(nsopi_);
    // Figure out total number of frozen docc/uocc orbitals
    int nfzc = 0;
    int nfzv = 0;
    for (int h = 0; h < nirrep_; h++) {
        nfzc += frzcpi_[h];
        nfzv += frzvpi_[h];
    }
    chkpt_->wt_nfzc(nfzc);
    chkpt_->wt_nfzv(nfzv);
    // These were computed by HF::finalize()
    chkpt_->wt_frzcpi(frzcpi_);
    chkpt_->wt_frzvpi(frzvpi_);

    int m = 0;
    for(int h = 0; h < nirrep_; ++h)
        if(soccpi_[h]) ++m;
    chkpt_->wt_iopen(m*(m+1)/2);

    if(options_.get_str("REFERENCE") == "UHF" ||
        options_.get_str("REFERENCE") == "CUHF"){

        double* values = epsilon_a_->to_block_vector();
        chkpt_->wt_alpha_evals(values);
        delete[] values;
        values = epsilon_b_->to_block_vector();
        chkpt_->wt_beta_evals(values);
        delete[] values;
        double** vectors = Ca_->to_block_matrix();
        chkpt_->wt_alpha_scf(vectors);
        delete[] vectors[0];
        delete[] vectors;
        vectors = Cb_->to_block_matrix();
        chkpt_->wt_beta_scf(vectors);
        delete[] vectors[0];
        delete[] vectors;
    }else{
        // All other reference type yield restricted orbitals
        double* values = epsilon_a_->to_block_vector();
        chkpt_->wt_evals(values);
        delete[] values;
        double** vectors = Ca_->to_block_matrix();
        chkpt_->wt_scf(vectors);
        delete[] vectors[0];
        delete[] vectors;
        double *ftmp = Fa_->to_lower_triangle();
        chkpt_->wt_fock(ftmp);
        delete[] ftmp;
    }
    psio_->close(PSIF_CHKPT, 1);
}

double HF::compute_energy()
{
    std::string reference = options_.get_str("REFERENCE");

    bool converged = false;
    MOM_performed_ = false;
    diis_performed_ = false;
    // Neither of these are idempotent
    if (options_.get_str("GUESS") == "SAD" || options_.get_str("GUESS") == "READ")
        iteration_ = -1;
    else
        iteration_ = 0;

    if (print_)
        outfile->Printf( "  ==> Pre-Iterations <==\n\n");

    if (print_)
        print_preiterations();

    // Andy trick 2.0
    std::string old_scf_type = options_.get_str("SCF_TYPE");
    if (options_.get_bool("DF_SCF_GUESS") && !(old_scf_type == "DF" || old_scf_type == "CD")) {
         outfile->Printf( "  Starting with a DF guess...\n\n");
         if(!options_["DF_BASIS_SCF"].has_changed()) {
             // TODO: Match Dunning basis sets
             molecule_->set_basis_all_atoms("CC-PVDZ-JKFIT", "DF_BASIS_SCF");
         }
         scf_type_ = "DF";
         options_.set_str("SCF","SCF_TYPE","DF"); // Scope is reset in proc.py. This is not pretty, but it works
    }

    if(attempt_number_ == 1){
        boost::shared_ptr<MintsHelper> mints (new MintsHelper(options_, 0));
        mints->one_electron_integrals();

        integrals();

        timer_on("Form H");
        form_H(); //Core Hamiltonian
        timer_off("Form H");

        // EFP: Add in permenent moment contribution and cache
        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
    	    boost::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_permanent();
    	    H_->add(Vefp);
	        Horig_ = SharedMatrix(new Matrix("H orig Matrix", basisset_->nbf(), basisset_->nbf()));
	        Horig_->copy(H_);
        }

        timer_on("Form S/X");
        form_Shalf(); //S and X Matrix
        timer_off("Form S/X");

        timer_on("Guess");
        guess(); // Guess
        timer_off("Guess");

    }else{
        // We're reading the orbitals from the previous set of iterations.
        form_D();
        E_ = compute_initial_E();
    }

    bool df = (options_.get_str("SCF_TYPE") == "DF");

        outfile->Printf( "  ==> Iterations <==\n\n");
        outfile->Printf( "%s                        Total Energy        Delta E     RMS |[F,P]|\n\n", df ? "   " : "");
    

    if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        Process::environment.get_efp()->set_qm_atoms();
        double efp_wfn_dependent_energy = Process::environment.get_efp()->scf_energy_update();
    }

    // SCF iterations
    do {
        iteration_++;

        save_density_and_energy();

        // Call any preiteration callbacks
        call_preiteration_callbacks();


        // add efp contribution to Fock matrix
        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
            H_->copy(Horig_);
    	    boost::shared_ptr<Matrix> Vefp = Process::environment.get_efp()->modify_Fock_induced();
    	    H_->add(Vefp);
        }

        E_ = 0.0;

        timer_on("Form G");
        form_G();
        timer_off("Form G");

        // Reset fractional SAD occupation
        if (iteration_ == 0 && options_.get_str("GUESS") == "SAD")
            reset_SAD_occupation();

        timer_on("Form F");
        form_F();
        timer_off("Form F");

        if (print_>3) {
            Fa_->print("outfile");
            Fb_->print("outfile");
        }

        E_ += compute_E();

        // add efp contribuation to energy
        if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
            double efp_wfn_dependent_energy =
		    Process::environment.get_efp()->scf_energy_update();
            //fprintf(outfile, "  Wfn dependent Energy =  %24.16f [H]\n",
			//    efp_wfn_dependent_energy);

            E_ += efp_wfn_dependent_energy;
        }   

        timer_on("DIIS");
        bool add_to_diis_subspace = false;
        if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ )
            add_to_diis_subspace = true;

        compute_orbital_gradient(add_to_diis_subspace);

        if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
            diis_performed_ = diis();
        } else {
            diis_performed_ = false;
        }
        timer_off("DIIS");

        if (print_>4 && diis_performed_) {
            outfile->Printf("  After DIIS:\n");
            Fa_->print("outfile");
            Fb_->print("outfile");
        }

        // If we're too well converged, or damping wasn't enabled, do DIIS
        damping_performed_ = (damping_enabled_ && iteration_ > 1 && Drms_ > damping_convergence_);

        std::string status = "";
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



        timer_on("Form C");
        form_C();
        timer_off("Form C");
        timer_on("Form D");
        form_D();
        timer_off("Form D");

        Process::environment.globals["SCF ITERATION ENERGY"] = E_;

        // After we've built the new D, damp the update if
        if(damping_performed_) damp_update();

        if (print_ > 3){
            Ca_->print("outfile");
            Cb_->print("outfile");
            Da_->print("outfile");
            Db_->print("outfile");
        }

        converged = test_convergency();

        df = (options_.get_str("SCF_TYPE") == "DF");


            outfile->Printf( "   @%s%s iter %3d: %20.14f   %12.5e   %-11.5e %s\n", df ? "DF-" : "",
                              reference.c_str(), iteration_, E_, E_ - Eold_, Drms_, status.c_str());
            


        // If a an excited MOM is requested but not started, don't stop yet
        if (MOM_excited_ && !MOM_started_) converged = false;

        // If a fractional occupation is requested but not started, don't stop yet
        if (frac_enabled_ && !frac_performed_) converged = false;

        // If a DF Guess environment, reset the JK object, and keep running
        if (converged && options_.get_bool("DF_SCF_GUESS") && !(old_scf_type == "DF" || old_scf_type == "CD")) {
            outfile->Printf( "\n  DF guess converged.\n\n"); // Be cool dude.
            converged = false;
            if(initialized_diis_manager_)
                diis_manager_->reset_subspace();
            scf_type_ = old_scf_type;
            options_.set_str("SCF","SCF_TYPE",old_scf_type);
            old_scf_type = "DF";
            integrals();
        }

        // Call any postiteration callbacks
        call_postiteration_callbacks();

    } while (!converged && iteration_ < maxiter_ );

    if ( Process::environment.get_efp()->get_frag_count() > 0 ) {
        Process::environment.get_efp()->Compute();

        double efp_elst         = Process::environment.globals["EFP ELST ENERGY"];
        double efp_exch         = Process::environment.globals["EFP EXCH ENERGY"];
        double efp_pol          = Process::environment.globals["EFP POL ENERGY"];
        double efp_disp         = Process::environment.globals["EFP DISP ENERGY"];
        double efp_total_energy = efp_elst + efp_exch + efp_pol + efp_disp;
        double efp_wfn_dependent_energy = efp_pol;

	E_ += efp_total_energy - efp_wfn_dependent_energy;

        fprintf(outfile, "  EFP Electrostatics Energy = %24.16f [H]\n", efp_elst);
        fprintf(outfile, "  EFP Polarization Energy =   %24.16f [H]\n", efp_pol);
        fprintf(outfile, "  EFP Dispersion Energy =     %24.16f [H]\n", efp_disp);
        fprintf(outfile, "  EFP Exchange Energy =       %24.16f [H]\n", efp_exch);
        fprintf(outfile, "  EFP Wfn dependent Energy =  %24.16f [H]\n", efp_wfn_dependent_energy);
        fprintf(outfile, "  EFP Total Energy =          %24.16f [H]\n", efp_total_energy);
        fprintf(outfile, "  Total SCF Energy =          %24.16f [H]\n", E_);
    }

    if (WorldComm->me() == 0)
        outfile->Printf( "\n  ==> Post-Iterations <==\n\n");

    check_phases();
    compute_spin_contamination();
    frac_renormalize();

    if (converged || !fail_on_maxiter_) {
        // Need to recompute the Fock matrices, as they are modified during the SCF interation
        // and might need to be dumped to checkpoint later
        form_F();

        // Print the orbitals
        if(print_)
            print_orbitals();

        if (converged) {
            outfile->Printf( "  Energy converged.\n\n");
        }
        if (!converged) {
            outfile->Printf( "  Energy did not converge, but proceeding anyway.\n\n");
        }
            outfile->Printf( "  @%s%s Final Energy: %20.14f", df ? "DF-" : "", reference.c_str(), E_);
            if (perturb_h_) {
                outfile->Printf( " with %f perturbation", lambda_);
            }
            outfile->Printf( "\n\n");
            print_energies();


        // Properties
        if (print_) {
            boost::shared_ptr<OEProp> oe(new OEProp());
            oe->set_title("SCF");
            oe->add("DIPOLE");

            if (print_ >= 2) {
                oe->add("QUADRUPOLE");
                oe->add("MULLIKEN_CHARGES");
            }

            if (print_ >= 3) {
                oe->add("LOWDIN_CHARGES");
                oe->add("MAYER_INDICES");
                oe->add("WIBERG_LOWDIN_INDICES");
            }

                outfile->Printf( "  ==> Properties <==\n\n");
            oe->compute();

            // TODO: Hack to test CubicScalarGrid
            /*
            boost::shared_ptr<CubicScalarGrid> grid(new CubicScalarGrid(basisset_));
            grid->print_header();
            grid->compute_density_cube(Da_, "Da");
            */

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

            Process::environment.globals["CURRENT DIPOLE X"] = Process::environment.globals["SCF DIPOLE X"];
            Process::environment.globals["CURRENT DIPOLE Y"] = Process::environment.globals["SCF DIPOLE Y"];
            Process::environment.globals["CURRENT DIPOLE Z"] = Process::environment.globals["SCF DIPOLE Z"];
        }

        save_information();
    } else {
            outfile->Printf( "  Failed to converged.\n");
            outfile->Printf( "    NOTE: MO Coefficients will not be saved to Checkpoint.\n");
        E_ = 0.0;
        if(psio_->open_check(PSIF_CHKPT))
            psio_->close(PSIF_CHKPT, 1);

        // Throw if we didn't converge?
        die_if_not_converged();
    }

    // Orbitals are always saved, in case an MO guess is requested later
    save_orbitals();
    if (options_.get_str("SAPT") != "FALSE") //not a bool because it has types
        save_sapt_info();

    // Perform wavefunction stability analysis
    if(options_.get_str("STABILITY_ANALYSIS") != "NONE")
        stability_analysis();

    // Clean memory off, handle diis closeout, etc
    finalize();


    //outfile->Printf("\nComputation Completed\n");
    
    return E_;
}

void HF::print_energies()
{
    outfile->Printf( "   => Energetics <=\n\n");
    outfile->Printf( "    Nuclear Repulsion Energy =        %24.16f\n", energies_["Nuclear"]);
    outfile->Printf( "    One-Electron Energy =             %24.16f\n", energies_["One-Electron"]);
    outfile->Printf( "    Two-Electron Energy =             %24.16f\n", energies_["Two-Electron"]);
    outfile->Printf( "    DFT Exchange-Correlation Energy = %24.16f\n", energies_["XC"]);
    outfile->Printf( "    Empirical Dispersion Energy =     %24.16f\n", energies_["-D"]);
    outfile->Printf( "    Total Energy =                    %24.16f\n", energies_["Nuclear"] +
        energies_["One-Electron"] + energies_["Two-Electron"] + energies_["XC"] + energies_["-D"]);
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
//  Comment so that autodoc utility will find this PSI variable
//     It doesn't really belong here but needs to be linked somewhere
//  Process::environment.globals["DOUBLE-HYBRID CORRECTION ENERGY"]
}

void HF::print_occupation()
{

        char **labels = molecule_->irrep_labels();
        std::string reference = options_.get_str("REFERENCE");
        outfile->Printf( "\t      ");
        for(int h = 0; h < nirrep_; ++h) outfile->Printf( " %4s ", labels[h]); outfile->Printf( "\n");
        outfile->Printf( "\tDOCC [ ");
        for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", doccpi_[h]);
        outfile->Printf( " %4d ]\n", doccpi_[nirrep_-1]);
        if(reference != "RHF" && reference != "RKS"){
            outfile->Printf( "\tSOCC [ ");
            for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", soccpi_[h]);
            outfile->Printf( " %4d ]\n", soccpi_[nirrep_-1]);
        }
        if (MOM_excited_) {
            // Also print nalpha and nbeta per irrep, which are more physically meaningful
            outfile->Printf( "\tNA   [ ");
            for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", nalphapi_[h]);
            outfile->Printf( " %4d ]\n", nalphapi_[nirrep_-1]);
            outfile->Printf( "\tNB   [ ");
            for(int h = 0; h < nirrep_-1; ++h) outfile->Printf( " %4d,", nbetapi_[h]);
            outfile->Printf( " %4d ]\n", nbetapi_[nirrep_-1]);
        }

        for(int h = 0; h < nirrep_; ++h) free(labels[h]); free(labels);
        outfile->Printf("\n");

}

//  Returns a vector of the occupation of the a orbitals
boost::shared_ptr<Vector> HF::occupation_a() const
{
  SharedVector occA = SharedVector(new Vector(nmopi_));
  for(int h=0; h < nirrep_;++h)
    for(int n=0; n < nalphapi()[h]; n++)
      occA->set(h, n, 1.0);

  return occA;
}

//  Returns a vector of the occupation of the b orbitals
boost::shared_ptr<Vector> HF::occupation_b() const
{
  SharedVector occB = SharedVector(new Vector(nmopi_));
  for(int h=0; h < nirrep_;++h)
    for(int n=0; n < nbetapi()[h]; n++)
      occB->set(h, n, 1.0);

  return occB;
}

void HF::diagonalize_F(const SharedMatrix& Fm, SharedMatrix& Cm, boost::shared_ptr<Vector>& epsm)
{
    //Form F' = X'FX for canonical orthogonalization
    diag_temp_->gemm(true, false, 1.0, X_, Fm, 0.0);
    diag_F_temp_->gemm(false, false, 1.0, diag_temp_, X_, 0.0);

    //Form C' = eig(F')
    diag_F_temp_->diagonalize(diag_C_temp_, epsm);

    //Form C = XC'
    Cm->gemm(false, false, 1.0, X_, diag_C_temp_, 0.0);
}

void HF::reset_SAD_occupation()
{
    // RHF style for now
    for (int h = 0; h < Da_->nirrep(); h++) {
        nalphapi_[h] = sad_nocc_[h];
        nbetapi_[h]  = sad_nocc_[h];
        doccpi_[h]   = sad_nocc_[h];
        soccpi_[h]   = 0;
    }
}
SharedMatrix HF::form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi)
{
    int* nsopi = Cso->rowspi();
    int* nmopi = Cso->colspi();
    int* nvirpi = new int[nirrep_];

    for (int h = 0; h < nirrep_; h++)
        nvirpi[h] = nmopi[h] - noccpi[h];

    SharedMatrix Fia(new Matrix("Fia (Some Basis)", nirrep_, noccpi, nvirpi));

    // Hack to get orbital e for this Fock
    SharedMatrix C2(new Matrix("C2", nirrep_, nsopi, nmopi));
    boost::shared_ptr<Vector> E2(new Vector("E2", nirrep_, nmopi));
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
    outfile->Printf( "\t");
    char** irrep_labels = molecule_->irrep_labels();
    int count = 0;
    for(; iter != vec.end(); ++iter){
        ++count;
        outfile->Printf( "%4s %-10.6f", irrep_labels[iter->second], iter->first);
        if(count == 4){
            outfile->Printf( "\n\t");
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
void HF::stability_analysis()
{
    throw PSIEXCEPTION("Stability analysis hasn't been implemented yet for this wfn type.");
}
}}
