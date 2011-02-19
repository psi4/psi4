/*
 *  hf.cpp
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */

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
#include <libparallel/parallel.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include <libmints/mints.h>

#include "hf.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

HF::HF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt),
      df_storage_(disk),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6),
      print_(3),
      add_external_potential_(false)
{
    common_init();
}

HF::HF(Options& options, shared_ptr<PSIO> psio)
    : Wavefunction(options, psio),
      df_storage_(disk),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6),
      print_(3),
      add_external_potential_(false)
{
    common_init();
}

HF::~HF()
{
}

void HF::common_init()
{
    scf_type_ = options_.get_str("SCF_TYPE");

    S_.reset(factory_.create_matrix("S"));
    Shalf_.reset(factory_.create_matrix("S^-1/2"));
    X_.reset(factory_.create_matrix("X"));
    Sphalf_.reset(factory_.create_matrix("S^+1/2"));
    H_.reset(factory_.create_matrix("One-electron Hamiltonion"));
    epsilon_a_.reset(factory_.create_vector());
    orbital_e_ = epsilon_a_;

    memset((void*) nsopi_, '\0', factory_.nirrep()*sizeof(int));
    memset((void*) nmopi_, '\0', factory_.nirrep()*sizeof(int));
    nmo_ = 0;
    nso_ = 0;
    nirrep_ = factory_.nirrep();
    int* dimpi = factory_.colspi();
    for (int h = 0; h< factory_.nirrep(); h++){
        nsopi_[h] = dimpi[h];
        nmopi_[h] = nsopi_[h]; //For now
        nso_ += nsopi_[h];
        nmo_ += nmopi_[h]; //For now (form_Shalf may change this, and will record things in the chkpt)
    }

    Eold_    = 0.0;
    E_       = 0.0;
    maxiter_ = 40;

    // Read information from input file
    maxiter_ = options_.get_int("MAXITER");

    // Read in DOCC and SOCC from memory
    int nirreps = factory_.nirrep();
    int ndocc = 0, nsocc = 0;
    input_docc_ = false;
    if (options_["DOCC"].has_changed()) {
        input_docc_ = true;
        if(options_["DOCC"].size() != nirreps)
            throw PSIEXCEPTION("The DOCC array has the wrong dimensions");
        for (int i=0; i<nirreps; ++i) {
            doccpi_[i] = options_["DOCC"][i].to_integer();
            ndocc += 2*doccpi_[i];
        }
    } else {
        for (int i=0; i<nirreps; ++i)
            doccpi_[i] = 0;
    }

    input_socc_ = false;
    if (options_["SOCC"].has_changed()) {
        input_socc_ = true;
        if(options_["SOCC"].size() != nirreps)
            throw PSIEXCEPTION("The SOCC array has the wrong dimensions");
        for (int i=0; i<nirreps; ++i) {
            soccpi_[i] = options_["SOCC"][i].to_integer();
            nsocc += soccpi_[i];
        }
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
            fprintf(outfile,"\tThere are an odd number of electrons - assuming doublet.\n"
                            "\tSpecify the multiplicity with the MULTP option in the\n"
                            "\tinput if this is incorrect\n\n");
        }else{
            multiplicity_ = 1;
            // There are an even number of electrons
            fprintf(outfile,"\tThere are an even number of electrons - assuming singlet.\n"
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

    // implicit case
    if (nirreps == 1) {
        doccpi_[0] = nbeta_;
        soccpi_[0] = nalpha_ - nbeta_;
    }


//  if (ndocc != 0 && nbeta_ != ndocc && nalpha_ != (ndocc + nsocc)) {
//      char *str = "Your DOCC, SOCC, charge, and multiplicity does not make sense.\n";
//      fprintf(outfile, str);
//      throw SanityCheckError(str, __FILE__, __LINE__);
//  }

    perturb_h_ = false;
    perturb_h_ = options_.get_bool("PERTURB_H");
    perturb_ = nothing;
    lambda_ = 0.0;
    if (perturb_h_) {
        string perturb_with;

        lambda_ = options_.get_double("LAMBDA");

        if (options_["PERTURB_WITH"].has_changed()) {
            perturb_with = options_.get_str("PERTURB_WITH");
            // Do checks to see what perturb_with is.
            if (perturb_with == "DIPOLE_X")
                perturb_ = dipole_x;
            else if (perturb_with == "DIPOLE_Y")
                perturb_ = dipole_y;
            else if (perturb_with == "DIPOLE_Z")
                perturb_ = dipole_z;
            else
                fprintf(outfile, "Unknown PERTURB_WITH. Applying no perturbation.\n");
        } else {
            fprintf(outfile, "PERTURB_H is true, but PERTURB_WITH not found, applying no perturbation.\n");
        }
    }

    // How much stuff shall we echo to the user?
    print_ = options_.get_int("PRINT");
    //fprintf(outfile,"  Print = %d\n",print_);

    //For HF algorithms, J and K are both required always.
    J_is_required_ = true;
    K_is_required_ = true;

    //Use schwarz sieve? default no
    schwarz_ = 0.0;
    if (options_["SCHWARZ_CUTOFF"].has_changed())
    {
        schwarz_ = options_.get_double("SCHWARZ_CUTOFF");
    }

    // Handle common diis info
    diis_enabled_ = true;
    min_diis_vectors_ = 4;

    // Allocate memory for DIISmin_diis_vectors_
    //  First, did the user request a different number of diis vectors?
    min_diis_vectors_ = options_.get_int("MIN_DIIS_VECTORS");
    max_diis_vectors_ = options_.get_int("MAX_DIIS_VECTORS");
    diis_start_ = options_.get_int("START_DIIS_ITER");
    diis_enabled_ = options_.get_bool("DIIS");

    // Don't perform DIIS if less than 2 vectors requested, or user requested a negative number
    if (min_diis_vectors_ < 2) {
        // disable diis
        diis_enabled_ = false;
    }

    initialized_diis_manager_ = false;

    // Save cartesian grid? Temporary until OEPROP is fully redone
    save_grid_ = false;
    if (options_.get_bool("SAVE_CARTESIAN_GRID")) {
        save_grid_ = true;
    }
    // Alloc memory for multipoles
    Dipole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Dipole X SO-basis")));
    Dipole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Dipole Y SO-basis")));
    Dipole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Dipole Z SO-basis")));
    Quadrupole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Quadrupole XX")));
    Quadrupole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Quadrupole XY")));
    Quadrupole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Quadrupole XZ")));
    Quadrupole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Quadrupole YY")));
    Quadrupole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Quadrupole YZ")));
    Quadrupole_.push_back(SharedSimpleMatrix(factory_.create_simple_matrix("Quadrupole ZZ")));

    if(Communicator::world->me() == 0)
        print_header();
    if (scf_type_ == "PK") {
        form_indexing();
        shared_ptr<MintsHelper> mints (new MintsHelper());
        mints->integrals();
    }
}
void HF::finalize()
{
    if (scf_type_ == "PK") {
        delete[] so2symblk_;
        delete[] so2index_;
        delete[] pk_symoffset_;
    }
    
    S_.reset();
    Shalf_.reset();
    Sphalf_.reset();
    X_.reset();
    H_.reset();
    Dipole_.clear();
    Quadrupole_.clear();

    // Clean up after DIIS
    if(initialized_diis_manager_)
        diis_manager_->delete_diis_file();
    diis_manager_.reset();
    initialized_diis_manager_ = false;

    // Close the chkpt
    psio_->close(PSIF_CHKPT, 1);
}
void HF::find_occupation(Vector & evals)
{
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<evals.nirrep(); ++h) {
        for (int i=0; i<evals.dimpi()[h]; ++i)
            pairs.push_back(make_pair(evals.get(h, i), h));
    }
    sort(pairs.begin(),pairs.end());

    if(!input_docc_){
        memset(doccpi_, 0, sizeof(int) * evals.nirrep());
        for (int i=0; i<nbeta_; ++i)
            doccpi_[pairs[i].second]++;
    }
    if(!input_socc_){
        memset(soccpi_, 0, sizeof(int) * evals.nirrep());
        for (int i=nbeta_; i<nalpha_; ++i)
            soccpi_[pairs[i].second]++;
    }

    if(print_>5 && Communicator::world->me() == 0){
        fprintf(outfile, "\tDOCC: [");
        for (int h=0; h<evals.nirrep(); ++h){
            fprintf(outfile, "%3d ", doccpi_[h]);
        }
        fprintf(outfile, "]\n");
        fprintf(outfile, "\tSOCC: [");
        for (int h=0; h<evals.nirrep(); ++h){
            fprintf(outfile, "%3d ", soccpi_[h]);
        }
        fprintf(outfile, "]\n");
    }

    for (int i=0; i<evals.nirrep(); ++i) {
        nalphapi_[i] = doccpi_[i] + soccpi_[i];
        nbetapi_[i]  = doccpi_[i];
    }
}

void HF::print_header()
{
    fprintf(outfile, " %s: by Justin Turney and Rob Parrish\n\n", options_.get_str("REFERENCE").c_str());

#ifdef _DEBUG
    fprintf(outfile, "  Debug version.\n\n");
#else
    fprintf(outfile, "  Release version.\n\n");
#endif

    molecule_->print();

    if(Communicator::world->me() == 0) {
        fprintf(outfile, "  Running in %s symmetry.\n\n", molecule_->point_group()->symbol());

        CharacterTable ct = molecule_->point_group()->char_table();

        fprintf(outfile, "  Input DOCC vector = (");
        for (int h=0; h<factory_.nirrep(); ++h) {
            fprintf(outfile, "%2d %3s ", doccpi_[h], ct.gamma(h).symbol());
        }
        fprintf(outfile, ")\n");
        fprintf(outfile, "  Input SOCC vector = (");
        for (int h=0; h<factory_.nirrep(); ++h) {
            fprintf(outfile, "%2d %3s ", soccpi_[h], ct.gamma(h).symbol());
        }

        fprintf(outfile, ")\n\n");
    }

    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", nuclearrep_);

    fprintf(outfile, "  Energy threshold  = %3.2e\n", energy_threshold_);
    fprintf(outfile, "  Density threshold = %3.2e\n\n", density_threshold_);
    fflush(outfile);
}

void HF::form_indexing()
{
    int h, i, ij, offset, pk_size;
    int nirreps = factory_.nirrep();
    int *opi = factory_.rowspi();

    so2symblk_ = new int[nso_];
    so2index_  = new int[nso_];

    ij = 0; offset = 0; pk_size = 0; pk_pairs_ = 0;
    for (h=0; h<nirreps; ++h) {
        for (i=0; i<opi[h]; ++i) {
            so2symblk_[ij] = h;
            so2index_[ij] = ij-offset;

            if (debug_ > 3 && Communicator::world->me() == 0)
                fprintf(outfile, "so2symblk_[%3d] = %3d, so2index_[%3d] = %3d\n", ij, so2symblk_[ij], ij, so2index_[ij]);

            ij++;
        }
        offset += opi[h];

        // Add up possible pair combinations that yield A1 symmetry
        pk_pairs_ += ioff[opi[h]];
    }

    // Compute the number of pairs in PK
    pk_size_ = INDEX2(pk_pairs_-1, pk_pairs_-1) + 1;

    // Compute PK symmetry mapping
    pk_symoffset_ = new int[nirreps];

    // Compute an offset in the PK matrix telling where a given symmetry block starts.
    pk_symoffset_[0] = 0;
    for (h=1; h<nirreps; ++h) {
        pk_symoffset_[h] = pk_symoffset_[h-1] + ioff[opi[h-1]];
    }
}

void HF::form_H()
{
    SharedMatrix kinetic(factory_.create_matrix("Kinetic Integrals"));
    SharedMatrix potential(factory_.create_matrix("Potential Integrals"));

    // Form the multipole integrals
    //form_multipole_integrals();

    // Load in kinetic and potential matrices
    double *integrals = init_array(ioff[nso_]);

    // Kinetic
    if (scf_type_ == "PK" || scf_type_ == "OUT_OF_CORE") {
        IWL::read_one(psio_.get(), PSIF_OEI, const_cast<char*>(PSIF_SO_T), integrals, ioff[nso_], 0, 0, outfile);
        kinetic->set(integrals);
        IWL::read_one(psio_.get(), PSIF_OEI, const_cast<char*>(PSIF_SO_V), integrals, ioff[nso_], 0, 0, outfile);
        potential->set(integrals);
    }
    else {
        // Integral factory
        shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
        shared_ptr<OneBodySOInt>    soT(integral->so_kinetic());
        shared_ptr<OneBodySOInt>    soV(integral->so_potential());

        soT->compute(kinetic);
        soV->compute(potential);
    }

    if (debug_ > 2 && Communicator::world->me() == 0)
        kinetic->print(outfile);

    if (debug_ > 2 && Communicator::world->me() == 0)
        potential->print(outfile);

    H_->copy(kinetic);
    H_->add(potential);

    if (debug_ > 2 && Communicator::world->me() == 0)
        H_->print(outfile);

    free(integrals);

    // if (perturb_h_) {
    //     if (perturb_ == dipole_x) {
    //         fprintf(outfile, "  Perturbing H by %f Dmx.\n", lambda_);
    //         H_.add(lambda_ * Dipole_[0]);
    //     } else if (perturb_ == dipole_y) {
    //         fprintf(outfile, "  Perturbing H by %f Dmy.\n", lambda_);
    //         H_.add(lambda_ * Dipole_[1]);
    //     } else if (perturb_ == dipole_z) {
    //         fprintf(outfile, "  Perturbing H by %f Dmz.\n", lambda_);
    //         H_.add(lambda_ * Dipole_[2]);
    //     }
    //     H_.print(outfile, "with perturbation");
    // }
}

void HF::form_Shalf()
{
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //          SYMMETRIC ORTHOGONALIZATION
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    // Overlap
    if (scf_type_ == "PK" || scf_type_ == "OUT_OF_CORE") {
        double *integrals = init_array(ioff[nso_]);
        IWL::read_one(psio_.get(), PSIF_OEI, const_cast<char*>(PSIF_SO_S), integrals, ioff[nso_], 0, 0, outfile);
        S_->set(integrals);
        free(integrals);
    }
    else {
        // Integral factory
        shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
        shared_ptr<OneBodySOInt>   so_overlap(integral->so_overlap());
        so_overlap->compute(S_);
    }
    // Form S^(-1/2) matrix
    Matrix eigvec;
    Matrix eigtemp;
    Matrix eigtemp2;
    Matrix eigvec_store;
    Vector eigval;
    Vector eigval_store;
    factory_.create_matrix(eigvec, "L");
    factory_.create_matrix(eigtemp, "Temp");
    factory_.create_matrix(eigtemp2);
    factory_.create_matrix(eigvec_store);
    factory_.create_vector(eigval);
    factory_.create_vector(eigval_store);

    //Used to do this 3 times, now only once
    S_->diagonalize(eigvec, eigval);
    eigvec_store.copy(eigvec);
    eigval_store.copy(eigval);

    // Convert the eigenvales to 1/sqrt(eigenvalues)
    int *dimpi = eigval.dimpi();
    double min_S = fabs(eigval.get(0,0));
    for (int h=0; h<eigval.nirrep(); ++h) {
        for (int i=0; i<dimpi[h]; ++i) {
            if (min_S > eigval.get(h,i))
                min_S = eigval.get(h,i);
            double scale = 1.0 / sqrt(eigval.get(h, i));
            eigval.set(h, i, scale);
        }
    }
    if (print_ && Communicator::world->me() == 0)
        fprintf(outfile,"\n  Minimum eigenvalue in the overlap matrix is %14.10E.\n",min_S);
    // Create a vector matrix from the converted eigenvalues
    eigtemp2.set(eigval);

    eigtemp.gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    Shalf_->gemm(false, false, 1.0, eigvec, eigtemp, 0.0);

    //Sphalf needs a fresh copy of the factorization of S
    eigvec.copy(eigvec_store);
    eigval.copy(eigval_store);
    // Convert the eigenvalues to sqrt(eigenvalues)
    for (int h=0; h<eigval.nirrep(); ++h) {
        for (int i=0; i<dimpi[h]; ++i) {
            double scale = sqrt(eigval.get(h, i));
            eigval.set(h, i, scale);
        }
    }
    // Create a vector matrix from the converted eigenvalues
    eigtemp2.set(eigval);

    // Works for diagonalize:
    eigtemp.gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    Sphalf_->gemm(false, false, 1.0, eigvec, eigtemp, 0.0);

    if (debug_ > 3 && Communicator::world->me() == 0) {
        Shalf_->print(outfile);
        Sphalf_->print(outfile);
    }
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    //
    //          CANONICAL ORTHOGONALIZATION
    //
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    //If symmetric orthogonalization will work, use it
    //Unless the user wants canonical
    double S_cutoff = options_.get_double("S_MIN_EIGENVALUE");
    if (min_S > S_cutoff && options_.get_str("S_ORTHOGONALIZATION") == "SYMMETRIC") {
        if (print_ && Communicator::world->me() == 0)
            fprintf(outfile,"  Using Symmetric Orthogonalization.\n");
        canonical_X_ = false;

        return;
    } else {
        if (print_ && Communicator::world->me() == 0)
            fprintf(outfile,"  Using Canonical Orthogonalization with cutoff of %14.10E.\n",S_cutoff);
        canonical_X_ = true;

        //Diagonalize S (or just get a fresh copy)
        eigvec.copy(eigvec_store);
        eigval.copy(eigval_store);
        int delta_mos = 0;
        for (int h=0; h<eigval.nirrep(); ++h) {
            //in each irrep, scale significant cols i  by 1.0/sqrt(s_i)
            int start_index = 0;
            for (int i=0; i<dimpi[h]; ++i) {
                if (S_cutoff  < eigval.get(h,i)) {
                    double scale = 1.0 / sqrt(eigval.get(h, i));
                    eigvec.scale_column(h, i, scale);
                } else {
                    start_index++;
                    nmopi_[h]--;
                    nmo_--;
                }
            }
            //eigvec.print(outfile);
            //Copy significant columns of eigvec into X in
            //descending order
            X_->zero();
            for (int i=0; i<dimpi[h]-start_index; ++i) {
                for (int m = 0; m < dimpi[h]; m++)
                    X_->set(h,m,i,eigvec.get(h,m,dimpi[h]-i-1));
            }
            //X_->print(outfile);
            if (print_>2 && Communicator::world->me() == 0)
                fprintf(outfile,"  Irrep %d, %d of %d possible MOs eliminated.\n",h,start_index,nsopi_[h]);

            delta_mos += start_index;
        }

        if (print_ && Communicator::world->me() == 0) {
            fprintf(outfile,"  Overall, %d of %d possible MOs eliminated.\n",delta_mos,nso_);

        }
        orbital_e_->init(eigvec.nirrep(), nmopi_);
        C_->init(eigvec.nirrep(),nsopi_,nmopi_,"MO coefficients");
    }
}

int *HF::compute_fcpi(int nfzc, SharedVector &eigvalues)
{
    int *frzcpi = new int[eigvalues->nirrep()];
    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<eigvalues->nirrep(); ++h) {
        for (int i=0; i<eigvalues->dimpi()[h]; ++i)
            pairs.push_back(make_pair(eigvalues->get(h, i), h));
        frzcpi[h] = 0;
    }
    sort(pairs.begin(),pairs.end());

    for (int i=0; i<nfzc; ++i)
        frzcpi[pairs[i].second]++;

    return frzcpi;
}

int *HF::compute_fvpi(int nfzv, SharedVector &eigvalues)
{
    int *frzvpi = new int[eigvalues->nirrep()];
    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<eigvalues->nirrep(); ++h) {
        for (int i=0; i<eigvalues->dimpi()[h]; ++i)
            pairs.push_back(make_pair(eigvalues->get(h, i), h));
        frzvpi[h] = 0;
    }
    sort(pairs.begin(),pairs.end(), greater<std::pair<double, int> >());

    for (int i=0; i<nfzv; ++i)
        frzvpi[pairs[i].second]++;

    return frzvpi;
}

void HF::form_multipole_integrals()
{
    // Initialize an integral object
    IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);

    // Get a dipole integral object
    OneBodyAOInt* dipole = integral.ao_dipole();
    OneBodyAOInt* quadrupole= integral.quadrupole();

    // Compute the dipole integrals
    dipole->compute(Dipole_);
    quadrupole->compute(Quadrupole_);

    delete quadrupole;
    delete dipole;

    // Get the nuclear contribution to the dipole
    nuclear_dipole_contribution_ = molecule_->nuclear_dipole_contribution();
    nuclear_quadrupole_contribution_ = molecule_->nuclear_quadrupole_contribution();

    // Save the dipole integrals
    Dipole_[0]->save(psio_, PSIF_OEI);
    Dipole_[1]->save(psio_, PSIF_OEI);
    Dipole_[2]->save(psio_, PSIF_OEI);

    Quadrupole_[0]->save(psio_, PSIF_OEI);
    Quadrupole_[1]->save(psio_, PSIF_OEI);
    Quadrupole_[2]->save(psio_, PSIF_OEI);
    Quadrupole_[3]->save(psio_, PSIF_OEI);
    Quadrupole_[4]->save(psio_, PSIF_OEI);
    Quadrupole_[5]->save(psio_, PSIF_OEI);
}

bool HF::load_or_compute_initial_C()
{
    bool ret = false;
    string prefix(chkpt_->build_keyword(const_cast<char*>("MO coefficients")));
    if (chkpt_->exist(const_cast<char*>(prefix.c_str()))) {
        // Read MOs from checkpoint and set C_ to them
        double **vectors;
        if(Communicator::world->me() == 0)
            vectors = chkpt_->rd_scf();
        else
            vectors = block_matrix(nso_, nmo_);
        Communicator::world->raw_bcast(&(vectors[0][0]), nso_*nmo_*sizeof(double));

        C_->set(const_cast<const double**>(vectors));
        free_block(vectors);

        form_D();

        // Read SCF energy from checkpoint file.
        if(Communicator::world->me() == 0)
            E_ = chkpt_->rd_escf();
        Communicator::world->bcast(E_);

        ret = true;
    } else {
        form_initial_C();
        form_D();
        // Compute an initial energy using H and D
        E_ = compute_initial_E();

        ret = false;
    }

    return ret;
}


}}
