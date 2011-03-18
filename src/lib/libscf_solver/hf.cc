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
#include "integralfunctors.h"
#include "pseudospectral.h"
#include "pkintegrals.h"
#include "df.h"

#include <libmints/mints.h>

#include "hf.h"

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

HF::HF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6),
      print_(1)
{
    common_init();
}

HF::HF(Options& options, shared_ptr<PSIO> psio)
    : Wavefunction(options, psio),
      nuclear_dipole_contribution_(3),
      nuclear_quadrupole_contribution_(6),
      print_(1)
{
    common_init();
}

HF::~HF()
{
}

void HF::common_init()
{
    integral_threshold_ = 1.0E-14;

    scf_type_ = options_.get_str("SCF_TYPE");

    S_.reset(factory_->create_matrix("S"));
    Shalf_.reset(factory_->create_matrix("S^-1/2"));
    X_.reset(factory_->create_matrix("X"));
    Sphalf_.reset(factory_->create_matrix("S^+1/2"));
    H_.reset(factory_->create_matrix("One-electron Hamiltonion"));

    memset((void*) nsopi_, '\0', factory_->nirrep()*sizeof(int));
    memset((void*) nmopi_, '\0', factory_->nirrep()*sizeof(int));
    nmo_ = 0;
    nso_ = 0;
    nirrep_ = factory_->nirrep();
    int* dimpi = factory_->colspi();
    for (int h = 0; h< factory_->nirrep(); h++){
        nsopi_[h] = dimpi[h];
        nmopi_[h] = nsopi_[h]; //For now
        nso_ += nsopi_[h];
        nmo_ += nmopi_[h]; //For now (form_Shalf may change this, and will record things in the chkpt)
    }

    // Form the SO lookup information
    so2symblk_ = new int[nso_];
    so2index_  = new int[nso_];
    size_t so_count = 0;
    size_t offset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < nsopi_[h]; ++i) {
            so2symblk_[so_count] = h;
            so2index_[so_count] = so_count-offset;
            ++so_count;
        }
        offset += nsopi_[h];
    }


    Eold_    = 0.0;
    E_       = 0.0;
    maxiter_ = 40;

    // Read information from input file
    maxiter_ = options_.get_int("MAXITER");

    // Read in DOCC and SOCC from memory
    int nirreps = factory_->nirrep();
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
        nalphapi_[0] = nalpha_;
        nbetapi_[0] = nbeta_;
    }

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
    if(options_["PRINT"].has_changed())
        print_ = options_.get_int("PRINT");

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

    print_header();

    // We need some integrals on disk for these cases
    if (scf_type_ == "PK" || scf_type_ == "OUT_OF_CORE"){
        shared_ptr<MintsHelper> mints (new MintsHelper());
        mints->integrals();
        if(scf_type_ == "PK") pk_integrals_ = shared_ptr<PKIntegrals>(new PKIntegrals(memory_, psio_, options_, nirrep_,
                                                                                     nsopi_, so2index_, so2symblk_));
    }else if (scf_type_ == "PSEUDOSPECTRAL"){
        if(nirrep_ > 1)
            throw PSIEXCEPTION("SCF TYPE " + scf_type_ + " cannot use symmetry yet. Add 'symmetry c1' to the molecule specification");
        df_ = shared_ptr<DFHF>(new DFHF(basisset_, psio_, options_));
        pseudospectral_ = shared_ptr<PseudospectralHF>(new PseudospectralHF(basisset_, psio_, options_));
    }else if (scf_type_ == "DF"){
        if(nirrep_ > 1)
            throw PSIEXCEPTION("SCF TYPE " + scf_type_ + " cannot use symmetry yet. Add 'symmetry c1' to the molecule specification");
        df_ = shared_ptr<DFHF>(new DFHF(basisset_, psio_, options_));
    }else if (scf_type_ == "DIRECT"){
        shared_ptr<IntegralFactory> integral = shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
        shared_ptr<TwoBodyAOInt> aoeri = shared_ptr<TwoBodyAOInt>(integral->eri());
        eri_ = shared_ptr<TwoBodySOInt>(new TwoBodySOInt(aoeri, integral));
    }
}

void HF::finalize()
{
    delete[] so2symblk_;
    delete[] so2index_;
    if (scf_type_ == "PK") pk_integrals_.reset();

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

    S_.reset();
    Shalf_.reset();
    Sphalf_.reset();
    X_.reset();
    H_.reset();

    // Close the chkpt
    if(psio_->open_check(PSIF_CHKPT))
        psio_->close(PSIF_CHKPT, 1);
}



void HF::find_occupation()
{
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<epsilon_a_->nirrep(); ++h) {
        for (int i=0; i<epsilon_a_->dimpi()[h]; ++i)
            pairs.push_back(make_pair(epsilon_a_->get(h, i), h));
    }
    sort(pairs.begin(),pairs.end());

    if(!input_docc_){
        memset(doccpi_, 0, sizeof(int) * epsilon_a_->nirrep());
        for (int i=0; i<nbeta_; ++i)
            doccpi_[pairs[i].second]++;
    }
    if(!input_socc_){
        memset(soccpi_, 0, sizeof(int) * epsilon_a_->nirrep());
        for (int i=nbeta_; i<nalpha_; ++i)
            soccpi_[pairs[i].second]++;
    }

    if(print_>5) {
        fprintf(outfile, "\tDOCC: [");
        for (int h=0; h<epsilon_a_->nirrep(); ++h){
            fprintf(outfile, "%3d ", doccpi_[h]);
        }
        fprintf(outfile, "]\n");
        fprintf(outfile, "\tSOCC: [");
        for (int h=0; h<epsilon_a_->nirrep(); ++h){
            fprintf(outfile, "%3d ", soccpi_[h]);
        }
        fprintf(outfile, "]\n");
    }

    for (int i=0; i<epsilon_a_->nirrep(); ++i) {
        nalphapi_[i] = doccpi_[i] + soccpi_[i];
        nbetapi_[i]  = doccpi_[i];
    }
}

void HF::print_header()
{
    fprintf(outfile, "\n %s: by Justin Turney and Rob Parrish\n\n", options_.get_str("REFERENCE").c_str());

#ifdef _DEBUG
    fprintf(outfile, "  Debug version.\n\n");
#else
    fprintf(outfile, "  Release version.\n\n");
#endif

    molecule_->print();

    fprintf(outfile, "  Running in %s symmetry.\n\n", molecule_->point_group()->symbol());

    CharacterTable ct = molecule_->point_group()->char_table();

    fprintf(outfile, "  Input DOCC vector = (");
    for (int h=0; h<factory_->nirrep(); ++h) {
        fprintf(outfile, "%2d %3s ", doccpi_[h], ct.gamma(h).symbol());
    }
    fprintf(outfile, ")\n");
    fprintf(outfile, "  Input SOCC vector = (");
    for (int h=0; h<factory_->nirrep(); ++h) {
        fprintf(outfile, "%2d %3s ", soccpi_[h], ct.gamma(h).symbol());
    }

    fprintf(outfile, ")\n\n");

    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", nuclearrep_);

    fprintf(outfile, "  Energy threshold  = %3.2e\n", energy_threshold_);
    fprintf(outfile, "  Density threshold = %3.2e\n\n", density_threshold_);
    fflush(outfile);
}



void HF::form_H()
{
    SharedMatrix kinetic(factory_->create_matrix("Kinetic Integrals"));
    SharedMatrix potential(factory_->create_matrix("Potential Integrals"));

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

    if (debug_ > 2)
        kinetic->print(outfile);

    if (debug_ > 2)
        potential->print(outfile);

    H_->copy(kinetic);
    H_->add(potential);

    if (debug_ > 2)
        H_->print(outfile);

    free(integrals);

    if (perturb_h_) {
        shared_ptr<IntegralFactory> ifact(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
        MultipoleSymmetry msymm(1, molecule_, ifact, factory_);
        vector<SharedMatrix> dipoles = msymm.create_matrices("Dipole");

        OneBodySOInt *so_dipole = ifact->so_dipole();
        so_dipole->compute(dipoles);

        if (perturb_ == dipole_x) {
            if (msymm.component_symmetry(0) != 0){
                fprintf(outfile, "  WARNING: You requested mu(x) perturbation, but mu(x) is not symmetric.\n");
            }
            else {
                fprintf(outfile, "  Perturbing H by %f mu(x).\n", lambda_);
                dipoles[0]->scale(lambda_);
                H_->add(dipoles[0]);
            }
        } else if (perturb_ == dipole_y) {
            if (msymm.component_symmetry(1) != 0){
                fprintf(outfile, "  WARNING: You requested mu(y) perturbation, but mu(y) is not symmetric.\n");
            }
            else {
                fprintf(outfile, "  Perturbing H by %f mu(y).\n", lambda_);
                dipoles[1]->scale(lambda_);
                H_->add(dipoles[1]);
            }
        } else if (perturb_ == dipole_z) {
            if (msymm.component_symmetry(2) != 0){
                fprintf(outfile, "  WARNING: You requested mu(z) perturbation, but mu(z) is not symmetric.\n");
            }
            else {
                fprintf(outfile, "  Perturbing H by %f mu(z).\n", lambda_);
                dipoles[2]->scale(lambda_);
                H_->add(dipoles[2]);
            }
        }
    }

    if (print_ > 3)
        H_->print(outfile);
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
    factory_->create_matrix(eigvec, "L");
    factory_->create_matrix(eigtemp, "Temp");
    factory_->create_matrix(eigtemp2);
    factory_->create_matrix(eigvec_store);
    factory_->create_vector(eigval);
    factory_->create_vector(eigval_store);

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
    if (print_)
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

    if (debug_ > 3) {
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
        if (print_)
            fprintf(outfile,"  Using Symmetric Orthogonalization.\n");
        canonical_X_ = false;

        return;
    } else {
        if (print_)
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
            if (print_>2)
                fprintf(outfile,"  Irrep %d, %d of %d possible MOs eliminated.\n",h,start_index,nsopi_[h]);

            delta_mos += start_index;
        }

        if (print_)
            fprintf(outfile,"  Overall, %d of %d possible MOs eliminated.\n",delta_mos,nso_);

        epsilon_a_->init(eigvec.nirrep(), nmopi_);
        Ca_->init(eigvec.nirrep(),nsopi_,nmopi_,"MO coefficients");
    }

    if (print_ > 3) {
        S_->print(outfile);
        Shalf_->print(outfile);
        if (canonical_X_)
            X_->print(outfile);
    }
}

void HF::compute_fcpi()
{
    int nfzc = molecule_->nfrozen_core();
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

void HF::compute_fvpi()
{
    int nfzv = options_.get_int("FREEZE_VIRT");
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

void HF::print_orbitals(const char* header, int *&irrep_count,
                        const std::vector< std::pair<double, int> >& evals,
                        int start, int end)
{
    char **labels = molecule_->irrep_labels();
    fprintf(outfile, "\t%-70s\n\n\t", header);
    int count = 0;
    for (int i = start; i < end; ++i) {
        int irrep = evals[i].second;
        fprintf(outfile, "%4d%-4s%10.6f  ", ++irrep_count[irrep], labels[irrep], evals[i].first);
        if (count++ % 3 == 2 && count != end)
            fprintf(outfile, "\n\t");
    }
    fprintf(outfile, "\n\n");
    for(int h = 0; h < nirrep_; ++h)
        delete [] labels[h];
}


bool HF::load_or_compute_initial_C()
{
    bool ret = false;
    bool loaded = false;

    //What does the user want?
    //Options will be:
    // ""-Either READ or CORE (try READ first)
    // "READ"-try to read MOs from checkpoint (restart style)
    // "BASIS2"-try to read MOs from checkpoint after basis cast up (in INPUT) NOT WORKING!
    // "DUAL_BASIS"-real the results of the DB computation from File 100, Temporary hack
    // "CORE"-CORE Hamiltonain
    // "GWH"-Generalized Wolfsberg-Helmholtz
    // "SAD"-Superposition of Atomic Denisties
    string guess_type = options_.get_str("GUESS");
    if (guess_type == "" || guess_type == "READ" || guess_type == "BASIS2") {

        // Read MOs from checkpoint and set C_ to them
        double **vectors;

        if (restricted()) {
            // Check to see if there are MOs already in the checkpoint file.
            // If so, read them in instead of forming them.
            string prefix(chkpt_->build_keyword(const_cast<char*>("MO coefficients")));

            if (chkpt_->exist(const_cast<char*>(prefix.c_str()))) {

                vectors = chkpt_->rd_scf();
                Ca_->set(const_cast<const double**>(vectors));
                free_block(vectors);

                double *orbitale;
                orbitale = chkpt_->rd_evals();
                epsilon_a_->set(orbitale);
                delete[] orbitale;

                loaded = true;
            }
        }
        else {
            string prefix(chkpt_->build_keyword(const_cast<char*>("Alpha MO coefficients")));
            if (chkpt_->exist(const_cast<char*>(prefix.c_str()))) {

                vectors = chkpt_->rd_alpha_scf();
                Ca_->set(const_cast<const double**>(vectors));
                free_block(vectors);

                double *orbitale;
                orbitale = chkpt_->rd_alpha_evals();
                epsilon_a_->set(orbitale);
                delete[] orbitale;
            }

            prefix = chkpt_->build_keyword(const_cast<char*>("Beta MO coefficients"));
            if (chkpt_->exist(const_cast<char*>(prefix.c_str()))) {
                vectors = chkpt_->rd_beta_scf();
                Ca_->set(const_cast<const double**>(vectors));
                free_block(vectors);

                double *orbitale;
                orbitale = chkpt_->rd_beta_evals();
                epsilon_b_->set(orbitale);
                delete[] orbitale;

                loaded = true;
            }
        }

        if (loaded) {
            //Try for existing MOs already. Deuces of loaded spades
           if (print_)
               fprintf(outfile, "  SCF Guess: Reading previous MOs.\n\n");

           // Read SCF energy from checkpoint file.
           E_ = chkpt_->rd_escf();

           ret = true;
        }
    }
    else if (guess_type == "DUAL_BASIS") {
         //Try for dual basis MOs,
        if (print_ && Communicator::world->me() == 0)
            fprintf(outfile, "  SCF Guess: Dual-Basis. Reading from File 100.\n\n");

        psio_->open(PSIF_SCF_DB_MOS,PSIO_OPEN_OLD);
        psio_->read_entry(PSIF_SCF_DB_MOS,"DB SCF Energy",(char *) &(E_),sizeof(double));
        psio_->read_entry(PSIF_SCF_DB_MOS,"DB NIRREPS",(char *) &(nirrep_),sizeof(int));
        psio_->read_entry(PSIF_SCF_DB_MOS,"DB DOCCPI",(char *) (doccpi_),8*sizeof(int));
        psio_->read_entry(PSIF_SCF_DB_MOS,"DB SOCCPI",(char *) (soccpi_),8*sizeof(int));
        psio_->read_entry(PSIF_SCF_DB_MOS,"DB NALPHAPI",(char *) (nalphapi_),8*sizeof(int));
        psio_->read_entry(PSIF_SCF_DB_MOS,"DB NBETAPI",(char *) (nbetapi_),8*sizeof(int));

        shared_ptr<Matrix> Ctemp(new Matrix("DUAL BASIS MOS", nirrep_, nsopi_, doccpi_));
        Ctemp->load(psio_, PSIF_SCF_DB_MOS, Matrix::SubBlocks);

        Ca_->zero();
        for (int h = 0; h < nirrep_; h++)
            for (int m = 0; m<nsopi_[h]; m++)
                for (int i = 0; i<doccpi_[h]; i++)
                    Ca_->set(h,m,i,Ctemp->get(h,m,i));

        psio_->close(PSIF_SCF_DB_MOS,1);
        ret = true;
    }
    else if (guess_type == "SAD") {
        if (print_)
            fprintf(outfile, "  SCF Guess: Superposition of Atomic Densities via on-the-fly atomic UHF.\n");

        //Superposition of Atomic Density (will be preferred when we can figure it out)
        compute_SAD_guess();
    }
    else if (guess_type == "GWH") {
        //Generalized Wolfsberg Helmholtz (Sounds cool, easy to code)
        if (print_)
            fprintf(outfile, "  SCF Guess: Generalized Wolfsberg-Helmholtz.\n\n");

        Fa_->zero(); //Try Fa_{mn} = S_{mn} (H_{mm} + H_{nn})/2
        int h, i, j;
        S_->print(outfile);
        int *opi = S_->rowspi();
        int nirreps = S_->nirrep();
        for (h=0; h<nirreps; ++h) {
            for (i=0; i<opi[h]; ++i) {
                for (j=0; j<opi[h]; ++j) {
                    Fa_->set(h,i,j,0.5*S_->get(h,i,j)*(H_->get(h,i,i)+H_->get(h,j,j)));
                }
            }
        }

        form_C();
    }
    else if (guess_type == "READ" || guess_type == "BASIS2") {
        throw std::invalid_argument("Checkpoint MOs requested, but do not exist!");
    }

    // If the user specified CORE or we tried to READ and we couldn't find the coefficients.
    if (guess_type == "CORE" || loaded == false) {
        //CORE is an old Psi standby, so we'll play this as spades
        if (print_)
            fprintf(outfile, "  SCF Guess: Core (One-Electron) Hamiltonian.\n\n");

        Fa_->copy(H_); //Try the core Hamiltonian as the Fock Matrix
        Fb_->copy(H_);

        if (options_.get_str("REFERENCE") == "ROHF")
            form_initial_C();
        else
            form_C();
    }

    return ret;
}


void HF::dump_to_checkpoint()
{
    if(!psio_->open_check(PSIF_CHKPT))
        psio_->open(PSIF_CHKPT, PSIO_OPEN_OLD);
    chkpt_->wt_nirreps(nirrep_);
    char **labels = molecule_->irrep_labels();
    chkpt_->wt_irr_labs(labels);
    for(int h = 0; h < nirrep_; ++h)
        delete [] labels[h];
    delete [] labels;
    chkpt_->wt_nmo(nmo_);
    chkpt_->wt_nso(nso_);
    chkpt_->wt_nao(basisset_->nao());
    chkpt_->wt_ref(0);
    chkpt_->wt_etot(E_);
    chkpt_->wt_escf(E_);
    chkpt_->wt_eref(E_);
    chkpt_->wt_enuc(molecule_->nuclear_repulsion_energy());
    chkpt_->wt_orbspi(nmopi_);
    chkpt_->wt_clsdpi(doccpi_);
    chkpt_->wt_openpi(soccpi_);
    chkpt_->wt_phase_check(0);
    chkpt_->wt_sopi(nsopi_);
    // Figure out frozen core orbitals
    int nfzc = molecule_->nfrozen_core();
    int nfzv = options_.get_int("FREEZE_VIRT");
    chkpt_->wt_nfzc(nfzc);
    chkpt_->wt_nfzv(nfzv);
    // These were computed by HF::finalize()
    chkpt_->wt_frzcpi(frzcpi_);
    chkpt_->wt_frzvpi(frzvpi_);

    int m = 0;
    for(int h = 0; h < nirrep_; ++h)
        if(soccpi_[h]) ++m;
    chkpt_->wt_iopen(m*(m+1)/2);

    if(options_.get_str("REFERENCE") == "UHF"){
        double* values = epsilon_a_->to_block_vector();
        chkpt_->wt_alpha_evals(values);
        free(values);
        values = epsilon_b_->to_block_vector();
        chkpt_->wt_beta_evals(values);
        free(values);
        double** vectors = Ca_->to_block_matrix();
        chkpt_->wt_alpha_scf(vectors);
        free_block(vectors);
        vectors = Cb_->to_block_matrix();
        chkpt_->wt_beta_scf(vectors);
        free_block(vectors);
    }else{
        // All other reference type yield restricted orbitals
        double* values = epsilon_a_->to_block_vector();
        chkpt_->wt_evals(values);
        free(values);
        double** vectors = Ca_->to_block_matrix();
        chkpt_->wt_scf(vectors);
        free_block(vectors);
        double *ftmp = Fa_->to_lower_triangle();
        chkpt_->wt_fock(ftmp);
        delete[] ftmp;
    }
    psio_->close(PSIF_CHKPT, 1);
}

double HF::compute_energy()
{
    std::string reference = options_.get_str("REFERENCE");

    bool converged = false, diis_iter = false;
    if (options_.get_str("GUESS") == "SAD")
        iteration_ = -1;
    else
        iteration_ = 0;

    form_H(); //Core Hamiltonian
    form_Shalf(); //Shalf Matrix

    // Form initial MO guess by user specified method
    // Check to see if there are MOs already in the checkpoint file.
    // If so, read them in instead of forming them, unless the user disagrees.
    load_or_compute_initial_C();

    // Maybe this goes here?
    if (options_.get_str("GUESS") == "SAD") {
        for (int h = 0 ; h < nirrep(); h++) {
            nalphapi_[h] = sad_nocc_[h];
        }
    }

    // Guess the occupation, if needed.
    find_occupation();

    form_D();

    // Compute an initial energy using H and D
    E_ = compute_initial_E();
    if (print_)
        fprintf(outfile, "  Initial %s energy: %20.14f\n\n", reference.c_str(), E_);

    fprintf(outfile, "                        Total Energy        Delta E      Density RMS\n\n");
    fflush(outfile);

    // SCF iterations
    do {
        iteration_++;

        save_density_and_energy();

        // Call any preiteration callbacks
        call_preiteration_callbacks();

        timer_on("Form G");
        form_G();
        timer_off("Form G");

//        if (print_>3) {
//            J_->print(outfile);
//            K_->print(outfile);
//            G_->print(outfile);
//        }

        form_F();
//        if (print_>3) {
//            Fa_->print(outfile);
//        }

        E_ = compute_E();

        timer_on("DIIS");
        if (diis_enabled_ && iteration_ > 0 && iteration_ >= diis_start_ )
            save_fock();
        if (diis_enabled_ == true && iteration_ >= diis_start_ + min_diis_vectors_ - 1) {
            diis_iter = diis();
        } else {
            diis_iter = false;
        }
        timer_off("DIIS");

//        if (print_>4 && diis_iter) {
//            fprintf(outfile,"  After DIIS:\n");
//            Fa_->print(outfile);
//        }
        fprintf(outfile, "   @%s iter %3d: %20.14f   % 10.5e   % 10.5e %s\n", reference.c_str(), iteration_, E_, E_ - Eold_, Drms_, diis_iter == false ? " " : "DIIS");
        fflush(outfile);

        form_C();
        form_D();

//        if (print_>2) {
//            Ca_->print(outfile);
//            D_->print(outfile);
//        }

        converged = test_convergency();

        // Call any postiteration callbacks
        call_postiteration_callbacks();

    } while (!converged && iteration_ < maxiter_ );

    if (converged) {
        // Need to recompute the Fock matrices, as they are modified during the SCF interation
        // and might need to be dumped to checkpoint later
        form_F();

        // Print the orbitals
        if(print_){
            fprintf(outfile, "\n\n\tOrbital Energies (a.u.)\n\t-----------------------\n\n");
            std::vector<std::pair<double, int> > aPairs;
            std::vector<std::pair<double, int> > bPairs;
            for (int h = 0; h < nirrep_; ++h) {
                for (int i=0; i < nmopi_[h]; ++i){
                    aPairs.push_back(make_pair(epsilon_a_->get(h, i), h));
                    bPairs.push_back(make_pair(epsilon_b_->get(h, i), h));
                }
            }
            sort(aPairs.begin(), aPairs.end());
            sort(bPairs.begin(), bPairs.end());
            int *irrep_count = new int[nirrep_];
            ::memset(irrep_count, 0, nirrep_ * sizeof(int));
            if(reference == "RHF"){
                print_orbitals("Doubly Occupied:", irrep_count, aPairs, 0, nalpha_);
                print_orbitals("Virtual:", irrep_count, aPairs, nalpha_, nmo_);
            }else if(reference == "UHF"){
                print_orbitals("Alpha Occupied:", irrep_count, aPairs, 0, nalpha_);
                print_orbitals("Alpha Virtual:", irrep_count, aPairs, nalpha_, nmo_);
                ::memset(irrep_count, 0, nirrep_ * sizeof(int));
                print_orbitals("Beta Occupied:", irrep_count, bPairs, 0, nbeta_);
                print_orbitals("Beta Virtual:", irrep_count, bPairs, nbeta_, nmo_);
            }else if(reference == "ROHF"){
                print_orbitals("Doubly Occupied:", irrep_count, aPairs, 0, nbeta_);
                print_orbitals("Singly Occupied:", irrep_count, aPairs, nbeta_, nalpha_);
                print_orbitals("Virtual:", irrep_count, aPairs, nalpha_, nmo_);
            }else{
                throw PSIEXCEPTION("Unknown reference in HF::print_orbitals");
            }

            char **labels = molecule_->irrep_labels();
            fprintf(outfile, "\tFinal Occupation by Irrep:\n");
            fprintf(outfile, "\t      ");
            for(int h = 0; h < nirrep_; ++h) fprintf(outfile, " %4s ", labels[h]); fprintf(outfile, "\n");
            fprintf(outfile, "\tDOCC [ ");
            for(int h = 0; h < nirrep_-1; ++h) fprintf(outfile, " %4d,", doccpi_[h]);
            fprintf(outfile, " %4d ]\n", doccpi_[nirrep_-1]);
            if(reference != "RHF"){
                fprintf(outfile, "\tSOCC [ ");
                for(int h = 0; h < nirrep_-1; ++h) fprintf(outfile, " %4d,", soccpi_[h]);
                fprintf(outfile, " %4d ]\n", soccpi_[nirrep_-1]);
            }
            for(int h = 0; h < nirrep_; ++h) delete[] labels[h]; delete[] labels;
            delete [] irrep_count;
        }

        fprintf(outfile, "\n  Energy converged.\n");
        fprintf(outfile, "\n  @%s Final Energy: %20.14f",reference.c_str(), E_);
        if (perturb_h_) {
            fprintf(outfile, " with %f perturbation", lambda_);
        }
        fprintf(outfile, "\n");
        save_information();
    } else {
        fprintf(outfile, "\n  Failed to converged.\n");
        E_ = 0.0;
        psio_->close(PSIF_CHKPT, 1);
    }

    //often, we're close!
    if (options_.get_bool("DUAL_BASIS"))
        save_dual_basis_projection();
    if (options_.get_str("SAPT") != "FALSE") //not a bool because it has types
        save_sapt_info();

    // Compute the final dipole.
//    compute_multipole();

    // Clean memory off, handle diis closeout, etc
    finalize();

    //fprintf(outfile,"\nComputation Completed\n");
    fflush(outfile);
    return E_;
}

}}
