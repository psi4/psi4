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
#include <libipv1/ip_lib.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>

#include "hf.h"

#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>
#include <libmints/molecule.h>

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
      addExternalPotential_(false)
{
    common_init();
}

HF::~HF()
{
    if (direct_integrals_ == false && ri_integrals_ == false) {
        delete[] so2symblk_;
        delete[] so2index_;
        delete[] pk_symoffset_;
    }
    free(zvals_);
}

void HF::common_init()
{
    S_.reset(factory_.create_matrix("S"));
    Shalf_.reset(factory_.create_matrix("S^-1/2"));
    Sphalf_.reset(factory_.create_matrix("S^+1/2"));
    H_.reset(factory_.create_matrix("One-electron Hamiltonion"));
    C_.reset(factory_.create_matrix("MO coefficients"));

    Eold_    = 0.0;
    E_       = 0.0;
    maxiter_ = 40;

    // Read information from input file
    maxiter_ = options_.get_int("MAXITER");

    // Read in DOCC and SOCC from memory
    int nirreps = factory_.nirreps();
    input_docc_ = false;
    if (options_["DOCC"].has_changed()) {
        input_docc_ = true;
        for (int i=0; i<nirreps; ++i)
            doccpi_[i] = options_["DOCC"][i].to_integer();
    } else {
        for (int i=0; i<nirreps; ++i)
            doccpi_[i] = 0;
    }
    input_socc_ = false;
    if (options_["SOCC"].has_changed()) {
        input_socc_ = true;
        for (int i=0; i<nirreps; ++i)
            soccpi_[i] = options_["SOCC"][i].to_integer();
    } else {
        for (int i=0; i<nirreps; ++i)
            soccpi_[i] = 0;
    }

    // Read information from checkpoint
    nuclearrep_ = chkpt_->rd_enuc();
    natom_ = chkpt_->rd_natom();
    zvals_ = chkpt_->rd_zvals();


    // Determine the number of electrons in the system
    charge_ = options_.get_int("CHARGE");
    nElec_  = 0;
    for (int i=0; i<natom_; ++i)
        nElec_ += (int)zvals_[i];
    nElec_ -= charge_;

    // If the user told us the multiplicity, read it from the input
    if(options_["MULTP"].has_changed()){
        multiplicity_ = options_.get_int("MULTP");
    }else{
        if(nElec_%2){
            multiplicity_ = 2;
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
    if(multiplicity_ - 1 > nElec_){
        char *str = new char[100];
        sprintf(str, "There are not enough electrons for multiplicity = %d, \n"
                     "please check your input and use the MULTP keyword", multiplicity_);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }
    if(multiplicity_%2 == nElec_%2){
        char *str = new char[100];
        sprintf(str, "A multiplicity of %d with %d electrons is impossible.\n"
                     "Please check your input and use the MULTP and/or CHARGE keywords",
                     multiplicity_, nElec_);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }

    nbeta_  = (nElec_ - multiplicity_ + 1)/2;
    nalpha_ = nbeta_ + multiplicity_ - 1;

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

    // Run integral direct? default no
    direct_integrals_ = false;
    direct_integrals_ = options_.get_bool("DIRECT");

    // How much stuff shall we echo to the user?
    print_ = options_.get_int("PRINT");

    //Run density fitting? default no
    ri_integrals_ = false;

    if (options_.get_bool("RI_HF"))
    {
        ri_integrals_ = true;
        direct_integrals_ = false;  // direct 4-index ints
    }
    
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
    num_diis_vectors_ = 4;

    // Allocate memory for DIISnum_diis_vectors_
    //  First, did the user request a different number of diis vectors?
    num_diis_vectors_ = options_.get_int("DIIS_VECTORS");
    diis_enabled_ = options_.get_bool("DIIS");

    // Don't perform DIIS if less than 2 vectors requested, or user requested a negative number
    if (num_diis_vectors_ < 2) {
        // disable diis
        diis_enabled_ = false;
    }

    // Initialize DIIS manager
    if (diis_enabled_)
        diis_manager_ = shared_ptr<DIISManager>(new DIISManager(num_diis_vectors_, "HF DIIS vector"));

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

    if(print_ > 1) print_header();
    if (direct_integrals_ == false && ri_integrals_ == false)
        form_indexing();
}

void HF::find_occupation(Vector & evals)
{
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<evals.nirreps(); ++h) {
        for (int i=0; i<evals.dimpi()[h]; ++i)
            pairs.push_back(make_pair(evals.get(h, i), h));
    }
    sort(pairs.begin(),pairs.end());

    if(!input_docc_){
        memset(doccpi_, 0, sizeof(int) * evals.nirreps());
        for (int i=0; i<nbeta_; ++i)
            doccpi_[pairs[i].second]++;
    }
    if(!input_socc_){
        memset(soccpi_, 0, sizeof(int) * evals.nirreps());
        for (int i=nbeta_; i<nalpha_; ++i)
            soccpi_[pairs[i].second]++;
    }

    if(print_>5){
        fprintf(outfile, "\tDOCC: [");
        for (int h=0; h<evals.nirreps(); ++h){
            fprintf(outfile, "%3d ", doccpi_[h]);
        }
        fprintf(outfile, "]\n");
        fprintf(outfile, "\tSOCC: [");
        for (int h=0; h<evals.nirreps(); ++h){
            fprintf(outfile, "%3d ", soccpi_[h]);
        }
        fprintf(outfile, "]\n");
    }

    for (int i=0; i<evals.nirreps(); ++i) {
        nalphapi_[i] = doccpi_[i] + soccpi_[i];
        nbetapi_[i]  = doccpi_[i];
    }
}

void HF::print_header()
{
    char *temp;
    char **temp2;
    char *reference;

    ip_string(const_cast<char*>("REFERENCE"), &reference, 0);

    fprintf(outfile, " %s: by Justin Turney and Rob Parrish\n\n", reference);
#ifdef _DEBUG
    fprintf(outfile, "  Debug version.\n");
#else
    fprintf(outfile, "  Release version.\n");
#endif
    temp = chkpt_->rd_sym_label();
    fprintf(outfile, "  Running in %s symmetry.\n", temp);
    free(temp);

    temp2 = chkpt_->rd_irr_labs();
    fprintf(outfile, "  Input DOCC vector = (");
    for (int h=0; h<factory_.nirreps(); ++h) {
        fprintf(outfile, "%2d %3s ", doccpi_[h], temp2[h]);
    }
    fprintf(outfile, ")\n");
    fprintf(outfile, "  Input SOCC vector = (");
    for (int h=0; h<factory_.nirreps(); ++h) {
        fprintf(outfile, "%2d %3s ", soccpi_[h], temp2[h]);
        free(temp2[h]);
    }
    free(temp2);

    fprintf(outfile, ")\n");
    fprintf(outfile, "  Nuclear repulsion = %20.15f\n", nuclearrep_);

    fprintf(outfile, "  Energy threshold  = %3.2e\n", energy_threshold_);
    fprintf(outfile, "  Density threshold = %3.2e\n\n", density_threshold_);
    fflush(outfile);
    free(reference);
}

void HF::form_indexing()
{
    int h, i, ij, offset, pk_size;
    int nirreps = factory_.nirreps();
    int *opi = factory_.rowspi();
    int nso;

    nso = chkpt_->rd_nso();
    so2symblk_ = new int[nso];
    so2index_  = new int[nso];

    ij = 0; offset = 0; pk_size = 0; pk_pairs_ = 0;
    for (h=0; h<nirreps; ++h) {
        for (i=0; i<opi[h]; ++i) {
            so2symblk_[ij] = h;
            so2index_[ij] = ij-offset;

            if (debug_ > 3)
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
    form_multipole_integrals();

    // Load in kinetic and potential matrices
    int nso = chkpt_->rd_nso();
    double *integrals = init_array(ioff[nso]);

    // Kinetic
    if (!direct_integrals_&&!ri_integrals_) {
        IWL::read_one(psio_.get(), PSIF_OEI, const_cast<char*>(PSIF_SO_T), integrals, ioff[nso], 0, 0, outfile);
        kinetic->set(integrals);
        IWL::read_one(psio_.get(), PSIF_OEI, const_cast<char*>(PSIF_SO_V), integrals, ioff[nso], 0, 0, outfile);
        potential->set(integrals);
    }
    else {
        IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
        shared_ptr<OneBodyInt> T(integral.kinetic());
        shared_ptr<OneBodyInt> V(integral.potential());

        T->compute(kinetic);
        V->compute(potential);
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
    int nso = chkpt_->rd_nso();

    // Overlap
    if (!direct_integrals_&&!ri_integrals_) {
        double *integrals = init_array(ioff[nso]);
        IWL::read_one(psio_.get(), PSIF_OEI, const_cast<char*>(PSIF_SO_S), integrals, ioff[nso], 0, 0, outfile);
        S_->set(integrals);
        free(integrals);
    }
    else {
        IntegralFactory integral(basisset_, basisset_, basisset_, basisset_);
        OneBodyInt *S = integral.overlap();
        S->compute(S_);
        delete S;
    }
    // Form S^(-1/2) matrix
    Matrix eigvec;
    Matrix eigtemp;
    Matrix eigtemp2;
    Vector eigval;
    factory_.create_matrix(eigvec, "L");
    factory_.create_matrix(eigtemp, "Temp");
    factory_.create_matrix(eigtemp2);
    factory_.create_vector(eigval);

    S_->diagonalize(eigvec, eigval);

    // Convert the eigenvales to 1/sqrt(eigenvalues)
    int *dimpi = eigval.dimpi();
    for (int h=0; h<eigval.nirreps(); ++h) {
        for (int i=0; i<dimpi[h]; ++i) {
            double scale = 1.0 / sqrt(eigval.get(h, i));
            eigval.set(h, i, scale);
        }
    }
    // Create a vector matrix from the converted eigenvalues
    eigtemp2.set(eigval);

    eigtemp.gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    Shalf_->gemm(false, false, 1.0, eigvec, eigtemp, 0.0);

    S_->diagonalize(eigvec, eigval);
    // Convert the eigenvalues to sqrt(eigenvalues)
    for (int h=0; h<eigval.nirreps(); ++h) {
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
}

int *HF::compute_fcpi(int nfzc, SharedVector &eigvalues)
{
    int *frzcpi = new int[eigvalues->nirreps()];
    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<eigvalues->nirreps(); ++h) {
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
    int *frzvpi = new int[eigvalues->nirreps()];
    // Print out orbital energies.
    std::vector<std::pair<double, int> > pairs;
    for (int h=0; h<eigvalues->nirreps(); ++h) {
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
    OneBodyInt* dipole = integral.dipole();
    OneBodyInt* quadrupole= integral.quadrupole();

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
        double **vectors = chkpt_->rd_scf();
        C_->set(const_cast<const double**>(vectors));
        free_block(vectors);

        form_D();

        // Read SCF energy from checkpoint file.
        E_ = chkpt_->rd_escf();

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
