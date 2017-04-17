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
#include "psi4/libmints/mintshelper.h"
#include "psi4/psifiles.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/sointegral_twobody.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/factory.h"
#include "psi4/libmints/3coverlap.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/sointegral_onebody.h"
#include "psi4/psi4-dec.h"
#include "x2cint.h"

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USING_dkh
#include <DKH/DKH_MANGLE.h>
#define F_DKH  DKH_MANGLE_MODULE(dkh_main, dkh, DKH_MAIN, DKH)

extern "C" {
    void F_DKH(double *S, double *V, double *T, double *pVp, int *nbf, int *dkh_order);
}
#endif

namespace psi {

/**
* IWLWriter functor for use with SO TEIs
**/
class IWLWriter
{
    IWL &writeto_;
    size_t count_;
    int &current_buffer_count_;

    Label *plabel_;
    Value *pvalue_;
public:

    IWLWriter(IWL &writeto) : writeto_(writeto), count_(0),
                              current_buffer_count_(writeto_.index())
    {
        plabel_ = writeto_.labels();
        pvalue_ = writeto_.values();
    }

    void operator()(int i, int j, int k, int l, int, int, int, int, int, int, int, int, double value)
    {
        int current_label_position = 4 * current_buffer_count_;

        // Save the labels
        plabel_[current_label_position++] = i;
        plabel_[current_label_position++] = j;
        plabel_[current_label_position++] = k;
        plabel_[current_label_position] = l;

        // Save the value
        pvalue_[current_buffer_count_++] = value;

        // Increment overall counter
        count_++;

        // If our IWL buffer is full dump to disk.
        if (current_buffer_count_ == writeto_.ints_per_buffer()) {
            writeto_.last_buffer() = 0;
            writeto_.buffer_count() = current_buffer_count_;
            writeto_.put();
            current_buffer_count_ = 0;
        }
    }

    size_t count() const
    { return count_; }
};

MintsHelper::MintsHelper(std::shared_ptr <BasisSet> basis, Options &options, int print, std::shared_ptr<BasisSet> ecpbasis)
        : options_(options), print_(print)
{
    init_helper(basis, ecpbasis);
}

MintsHelper::MintsHelper(std::shared_ptr <Wavefunction> wavefunction)
        : options_(wavefunction->options())
{
    init_helper(wavefunction);
}

MintsHelper::~MintsHelper()
{
}

void MintsHelper::init_helper(std::shared_ptr <Wavefunction> wavefunction)
{

    if (wavefunction->basisset().get() == 0) {
        outfile->Printf("  Wavefunction does not have a basisset!");
        throw PSIEXCEPTION("Wavefunction does not have a basisset, what did you do?!");
    }

    psio_ = wavefunction->psio();
    basisset_ = wavefunction->basisset();
    molecule_ = basisset_->molecule();

    // Make sure molecule is valid.
    molecule_->update_geometry();

    common_init();
}

void MintsHelper::init_helper(std::shared_ptr <BasisSet> basis,  std::shared_ptr<BasisSet> ecpbasis)
{
    basisset_ = basis;
	ecpbasis_ = ecpbasis; 
    molecule_ = basis->molecule();
    psio_ = _default_psio_lib_;

    // Make sure molecule is valid.
    molecule_->update_geometry();

    common_init();
}

void MintsHelper::common_init()
{
    // Print the molecule.
    if (print_)
        molecule_->print();

    // Print the basis set
    if (print_)
        basisset_->print_detail();

    // How many threads?
    nthread_ = 1;
    #ifdef _OPENMP
        nthread_ = Process::environment.get_n_threads();
    #endif

    // Create integral factory
    integral_ = std::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, ecpbasis_));

    // Get the SO basis object.
    sobasis_ = std::shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral_));

    // Obtain dimensions from the sobasis
    const Dimension dimension = sobasis_->dimension();

    // Create a matrix factory and initialize it
    factory_ = std::shared_ptr<MatrixFactory>(new MatrixFactory());
    factory_->init_with(dimension, dimension);

    // Integral cutoff
    cutoff_ = Process::environment.options.get_double("INTS_TOLERANCE");
}

std::shared_ptr <PetiteList> MintsHelper::petite_list() const
{
    std::shared_ptr <PetiteList> pt(new PetiteList(basisset_, integral_));
    return pt;
}

std::shared_ptr <PetiteList> MintsHelper::petite_list(bool val) const
{
    std::shared_ptr <PetiteList> pt(new PetiteList(basisset_, integral_, val));
    return pt;
}

std::shared_ptr <BasisSet> MintsHelper::basisset() const
{
    return basisset_;
}

std::shared_ptr <SOBasisSet> MintsHelper::sobasisset() const
{
    return sobasis_;
}

std::shared_ptr<BasisSet> MintsHelper::ecpbasisset() const
{
	return ecpbasis_; 
}

std::shared_ptr <MatrixFactory> MintsHelper::factory() const
{
    return factory_;
}

std::shared_ptr <IntegralFactory> MintsHelper::integral() const
{
    return integral_;
}

int MintsHelper::nbf() const
{
    return basisset_->nbf();
}

void MintsHelper::integrals()
{
    if (print_) { outfile->Printf(" MINTS: Wrapper to libmints.\n   by Justin Turney\n\n"); }

    // Get ERI object
    std::vector <std::shared_ptr<TwoBodyAOInt>> tb;
    for (int i = 0; i < nthread_; ++i)
        tb.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->eri()));
    std::shared_ptr <TwoBodySOInt> eri(new TwoBodySOInt(tb, integral_));

    //// Print out some useful information
    if (print_) {
        outfile->Printf("   Calculation information:\n");
        outfile->Printf("      Number of threads:              %4d\n", nthread_);
        outfile->Printf("      Number of atoms:                %4d\n", molecule_->natom());
        outfile->Printf("      Number of AO shells:            %4d\n", basisset_->nshell());
        outfile->Printf("      Number of SO shells:            %4d\n", sobasis_->nshell());
        outfile->Printf("      Number of primitives:           %4d\n", basisset_->nprimitive());
        outfile->Printf("      Number of atomic orbitals:      %4d\n", basisset_->nao());
        outfile->Printf("      Number of basis functions:      %4d\n\n", basisset_->nbf());
        outfile->Printf("      Number of irreps:               %4d\n", sobasis_->nirrep());
        outfile->Printf("      Integral cutoff                 %4.2e\n", cutoff_);
        outfile->Printf("      Number of functions per irrep: [");
        for (int i = 0; i < sobasis_->nirrep(); ++i) {
            outfile->Printf("%4d ", sobasis_->nfunction_in_irrep(i));
        }
        outfile->Printf("]\n\n");
    }

    // Compute one-electron integrals.
    one_electron_integrals();

    // Open the IWL buffer where we will store the integrals.
    IWL ERIOUT(psio_.get(), PSIF_SO_TEI, cutoff_, 0, 0);
    IWLWriter writer(ERIOUT);

    // Let the user know what we're doing.
    if (print_) { outfile->Printf("      Computing two-electron integrals..."); }

    SOShellCombinationsIterator shellIter(sobasis_, sobasis_, sobasis_, sobasis_);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        eri->compute_shell(shellIter, writer);
    }

    // Flush out buffers.
    ERIOUT.flush(1);

    // We just did all this work to create the file, let's keep it around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    if (print_) {
        outfile->Printf("done\n");
        outfile->Printf("      Computed %lu non-zero two-electron integrals.\n"
                                "        Stored in file %d.\n\n", writer.count(), PSIF_SO_TEI);
    }
}

void MintsHelper::integrals_erf(double w)
{
    double omega = (w == -1.0 ? options_.get_double("OMEGA_ERF") : w);

    IWL ERIOUT(psio_.get(), PSIF_SO_ERF_TEI, cutoff_, 0, 0);
    IWLWriter writer(ERIOUT);

    // Get ERI object
    std::vector <std::shared_ptr<TwoBodyAOInt>> tb;
    for (int i = 0; i < nthread_; ++i)
        tb.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->erf_eri(omega)));
    std::shared_ptr <TwoBodySOInt> erf(new TwoBodySOInt(tb, integral_));

    // Let the user know what we're doing.
    outfile->Printf("      Computing non-zero ERF integrals (omega = %.3f)...", omega);

    SOShellCombinationsIterator shellIter(sobasis_, sobasis_, sobasis_, sobasis_);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next())
        erf->compute_shell(shellIter, writer);

    // Flush the buffers
    ERIOUT.flush(1);

    // Keep the integrals around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    outfile->Printf("done\n");
    outfile->Printf("      Computed %lu non-zero ERF integrals.\n"
                            "        Stored in file %d.\n\n", writer.count(), PSIF_SO_ERF_TEI);
}

void MintsHelper::integrals_erfc(double w)
{
    double omega = (w == -1.0 ? options_.get_double("OMEGA_ERF") : w);

    IWL ERIOUT(psio_.get(), PSIF_SO_ERFC_TEI, cutoff_, 0, 0);
    IWLWriter writer(ERIOUT);

    // Get ERI object
    std::vector <std::shared_ptr<TwoBodyAOInt>> tb;
    for (int i = 0; i < nthread_; ++i)
        tb.push_back(std::shared_ptr<TwoBodyAOInt>(integral_->erf_complement_eri(omega)));
    std::shared_ptr <TwoBodySOInt> erf(new TwoBodySOInt(tb, integral_));

    // Let the user know what we're doing.
    outfile->Printf("      Computing non-zero ERFComplement integrals...");

    SOShellCombinationsIterator shellIter(sobasis_, sobasis_, sobasis_, sobasis_);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next())
        erf->compute_shell(shellIter, writer);

    // Flush the buffers
    ERIOUT.flush(1);

    // Keep the integrals around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    outfile->Printf("done\n");
    outfile->Printf("      Computed %lu non-zero ERFComplement integrals.\n"
                            "        Stored in file %d.\n\n", writer.count(), PSIF_SO_ERFC_TEI);
}


void MintsHelper::one_electron_integrals()
{
//    outfile->Printf( " OEINTS: Wrapper to libmints.\n   by Justin Turney\n\n");
//
//    // Print out some useful information
//    outfile->Printf( "   Calculation information:\n");
//    outfile->Printf( "      Number of atoms:                %4d\n", molecule_->natom());
//    outfile->Printf( "      Number of AO shells:            %4d\n", basisset_->nshell());
//    outfile->Printf( "      Number of SO shells:            %4d\n", sobasis_->nshell());
//    outfile->Printf( "      Number of primitives:           %4d\n", basisset_->nprimitive());
//    outfile->Printf( "      Number of atomic orbitals:      %4d\n", basisset_->nao());
//    outfile->Printf( "      Number of basis functions:      %4d\n\n", basisset_->nbf());
//    outfile->Printf( "      Number of irreps:               %4d\n", sobasis_->nirrep());
//    outfile->Printf( "      Number of functions per irrep: [");
//    for (int i=0; i<sobasis_->nirrep(); ++i) {
//        outfile->Printf( "%4d ", sobasis_->nfunction_in_irrep(i));
//    }
//    outfile->Printf( "]\n\n");

    // Compute and dump one-electron SO integrals.

    if (options_.get_str("RELATIVISTIC") == "NO" || options_.get_str("RELATIVISTIC") == "DKH") {
        // Overlap
        so_overlap()->save(psio_, PSIF_OEI);

        // Kinetic
        so_kinetic()->save(psio_, PSIF_OEI);

        // Potential -- DKH perturbation added to potential integrals if needed.
        so_potential()->save(psio_, PSIF_OEI);
    } else if (options_.get_str("RELATIVISTIC") == "X2C") {
        outfile->Printf(" OEINTS: Using relativistic (X2C) overlap, kinetic, and potential integrals.\n");

        if (!rel_basisset_){
            throw PSIEXCEPTION("OEINTS: X2C requested, but relativistic basis was not set.");
        }
        X2CInt x2cint;
        SharedMatrix so_overlap_x2c = so_overlap();
        SharedMatrix so_kinetic_x2c = so_kinetic();
        SharedMatrix so_potential_x2c = so_potential();
        x2cint.compute(basisset_, rel_basisset_, so_overlap_x2c, so_kinetic_x2c, so_potential_x2c);

        // Overlap
        so_overlap_x2c->save(psio_, PSIF_OEI);

        // Kinetic
        so_kinetic_x2c->save(psio_, PSIF_OEI);

        // Potential
        so_potential_x2c->save(psio_, PSIF_OEI);
    }

    // Dipoles
    std::vector <SharedMatrix> dipole_mats = so_dipole();
    for (SharedMatrix m : dipole_mats) {
        m->save(psio_, PSIF_OEI);
    }

    // Quadrupoles
    std::vector <SharedMatrix> quadrupole_mats = so_quadrupole();
    for (SharedMatrix m : quadrupole_mats) {
        m->save(psio_, PSIF_OEI);
    }

    if (print_) {
        outfile->Printf(" OEINTS: Overlap, kinetic, potential, dipole, and quadrupole integrals\n"
                                "         stored in file %d.\n\n", PSIF_OEI);
    }
}

void MintsHelper::integral_gradients()
{
    throw FeatureNotImplemented("libmints", "MintsHelper::integral_derivatives", __FILE__, __LINE__);
}

void MintsHelper::integral_hessians()
{
    throw FeatureNotImplemented("libmints", "MintsHelper::integral_hessians", __FILE__, __LINE__);
}

void MintsHelper::one_body_ao_computer(std::vector<std::shared_ptr<OneBodyAOInt>> ints,
                                       SharedMatrix out, bool symm) {
    // Grab basis info
    std::shared_ptr<BasisSet> bs1 = ints[0]->basis1();
    std::shared_ptr<BasisSet> bs2 = ints[0]->basis2();

    // Limit to the number of incoming onbody ints
    size_t nthread = nthread_;
    if (nthread > ints.size()) {
        nthread = ints.size();
    }

    // Grab the buffers
    std::vector<const double *> ints_buff(nthread);
    for (size_t thread = 0; thread < nthread; thread++) {
        ints_buff[thread] = ints[thread]->buffer();
    }

    double **outp = out->pointer();

    // Loop it
    #pragma omp parallel for schedule(guided) num_threads(nthread)
    for (size_t MU = 0; MU < bs1->nshell(); ++MU) {
        const size_t num_mu = bs1->shell(MU).nfunction();
        const size_t index_mu = bs1->shell(MU).function_index();

        size_t rank = 0;
#ifdef _OPENMP
        rank = omp_get_thread_num();
#endif

        if (symm) {
            // Triangular
            for (size_t NU = 0; NU <= MU; ++NU) {
                const size_t num_nu = bs2->shell(NU).nfunction();
                const size_t index_nu = bs2->shell(NU).function_index();

                ints[rank]->compute_shell(MU, NU);

                size_t index = 0;
                for (size_t mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                    for (size_t nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                        outp[nu][mu] = outp[mu][nu] = ints_buff[rank][index++];
                    }
                }
            }  // End NU
        }      // End Symm
        else {
            // Rectangular
            for (size_t NU = 0; NU < bs2->nshell(); ++NU) {
                const size_t num_nu = bs2->shell(NU).nfunction();
                const size_t index_nu = bs2->shell(NU).function_index();

                ints[rank]->compute_shell(MU, NU);

                size_t index = 0;
                for (size_t mu = index_mu; mu < (index_mu + num_mu); ++mu) {
                    for (size_t nu = index_nu; nu < (index_nu + num_nu); ++nu) {
                        // printf("%zu %zu | %zu %zu | %lf\n", MU, NU, mu, nu, ints_buff[rank][index]);
                        outp[mu][nu] = ints_buff[rank][index++];
                    }
                }
            }  // End NU
        }      // End Rectangular
    }          // End Mu
}
SharedMatrix MintsHelper::ao_overlap()
{
    // Overlap
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++){
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_overlap()));
    }
    SharedMatrix overlap_mat(new Matrix(PSIF_AO_S, basisset_->nbf(), basisset_->nbf()));
    one_body_ao_computer(ints_vec, overlap_mat, true);

    // Is this needed?
    overlap_mat->save(psio_, PSIF_OEI);
    return overlap_mat;
}

// JWM 4/3/2015
SharedMatrix MintsHelper::ao_overlap(std::shared_ptr <BasisSet> bs1, std::shared_ptr <BasisSet> bs2)
{
    // Overlap
    IntegralFactory factory(bs1, bs2);
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++){
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(factory.ao_overlap()));
    }
    SharedMatrix overlap_mat(new Matrix(PSIF_AO_S, bs1->nbf(), bs2->nbf()));
    one_body_ao_computer(ints_vec, overlap_mat, false);
    return overlap_mat;
}

SharedMatrix MintsHelper::ao_kinetic()
{
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++){
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_kinetic()));
    }
    SharedMatrix kinetic_mat(new Matrix("AO-basis Kinetic Ints", basisset_->nbf(), basisset_->nbf()));
    one_body_ao_computer(ints_vec, kinetic_mat, true);
    return kinetic_mat;
}

SharedMatrix MintsHelper::ao_kinetic(std::shared_ptr <BasisSet> bs1, std::shared_ptr <BasisSet> bs2)
{
    IntegralFactory factory(bs1, bs2);
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++){
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(factory.ao_kinetic()));
    }
    SharedMatrix kinetic_mat(new Matrix("AO-basis Kinetic Ints", bs1->nbf(), bs2->nbf()));
    one_body_ao_computer(ints_vec, kinetic_mat, false);
    return kinetic_mat;
}

SharedMatrix MintsHelper::ao_potential()
{
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++){
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_potential()));
    }
    SharedMatrix potential_mat(new Matrix("AO-basis Potential Ints", basisset_->nbf(), basisset_->nbf()));
    one_body_ao_computer(ints_vec, potential_mat, true);
    return potential_mat;
}

SharedMatrix MintsHelper::ao_potential(std::shared_ptr <BasisSet> bs1, std::shared_ptr <BasisSet> bs2)
{
    IntegralFactory factory(bs1, bs2);
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++){
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(factory.ao_potential()));
    }
    SharedMatrix potential_mat(new Matrix("AO-basis Potential Ints", bs1->nbf(), bs2->nbf()));
    one_body_ao_computer(ints_vec, potential_mat, false);
    return potential_mat;
}

SharedMatrix MintsHelper::ao_ecp()
{
	std::shared_ptr <OneBodyAOInt> ECP(integral_->ao_ecp());
	SharedMatrix ecp_mat(new Matrix("AO-basis ECP Ints", basisset_->nbf(), basisset_->nbf()));
	ECP->compute(ecp_mat);
	return ecp_mat;
}

SharedMatrix MintsHelper::ao_ecp(std::shared_ptr <BasisSet> bs1, std::shared_ptr <BasisSet> bs2, std::shared_ptr <BasisSet> bsecp)
{
    IntegralFactory factory(bs1, bs2, bs1, bs2, bsecp);
	std::shared_ptr <OneBodyAOInt> ECP(factory.ao_ecp());
	SharedMatrix ecp_mat(new Matrix("AO-basis ECP Ints", basisset_->nbf(), basisset_->nbf()));
	ECP->compute(ecp_mat);
	return ecp_mat; 
}

SharedMatrix MintsHelper::ao_pvp()
{
    std::vector<std::shared_ptr<OneBodyAOInt>> ints_vec;
    for (size_t i = 0; i < nthread_; i++){
        ints_vec.push_back(std::shared_ptr<OneBodyAOInt>(integral_->ao_rel_potential()));
    }
    SharedMatrix pVp_mat(new Matrix("AO-basis pVp Ints", basisset_->nbf(), basisset_->nbf()));
    one_body_ao_computer(ints_vec, pVp_mat, true);
    return pVp_mat;
}

SharedMatrix MintsHelper::ao_dkh(int dkh_order)
{
#ifdef USING_dkh
    SharedMatrix S = ao_overlap();
    SharedMatrix T = ao_kinetic();
    SharedMatrix Torig = T->clone();
    SharedMatrix V = ao_potential();
    SharedMatrix Vorig = V->clone();
    SharedMatrix pVp = ao_pvp();
    SharedMatrix H_dk = T->clone();
    H_dk->zero();

    double *Sp   = S->pointer()[0];
    double *Tp   = T->pointer()[0];
    double *Vp   = V->pointer()[0];
    double *pVpp = pVp->pointer()[0];

    if (dkh_order < 1)
        dkh_order = 2;
    if (dkh_order > 4)
        dkh_order = 4;

    outfile->Printf("    Computing %d-order Douglas-Kroll-Hess integrals.\n", dkh_order);

    int nbf = basisset_->nbf();

//    V->print();
//    T->print();

    // Call DKH code from Markus Reiher
    F_DKH(Sp, Vp, Tp, pVpp, &nbf, &dkh_order);

    H_dk->add(V);
    H_dk->add(T);
    H_dk->subtract(Vorig);
    H_dk->subtract(Torig);

    H_dk->set_name("AO-basis DKH Ints");
    //H_dk->print();

    return H_dk;
#else
    UNUSED(dkh_order);
    outfile->Printf("    Douglas-Kroll-Hess integrals requested but are not available.\n");
    throw PSIEXCEPTION("Douglas-Kroll-Hess integrals requested but were not compiled in.");
#endif
}

SharedMatrix MintsHelper::so_dkh(int dkh_order)
{
    SharedMatrix dkh = factory_->create_shared_matrix("SO Douglas-Kroll-Hess Integrals");
    dkh->apply_symmetry(ao_dkh(dkh_order), petite_list()->aotoso());
    return dkh;
}

SharedMatrix MintsHelper::ao_helper(const std::string &label, std::shared_ptr <TwoBodyAOInt> ints)
{
    std::shared_ptr <BasisSet> bs1 = ints->basis1();
    std::shared_ptr <BasisSet> bs2 = ints->basis2();
    std::shared_ptr <BasisSet> bs3 = ints->basis3();
    std::shared_ptr <BasisSet> bs4 = ints->basis4();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int nbf3 = bs3->nbf();
    int nbf4 = bs4->nbf();

    SharedMatrix I(new Matrix(label, nbf1 * nbf2, nbf3 * nbf4));
    double **Ip = I->pointer();
    const double *buffer = ints->buffer();


    for (int M = 0; M < bs1->nshell(); M++) {
        for (int N = 0; N < bs2->nshell(); N++) {
            for (int P = 0; P < bs3->nshell(); P++) {
                for (int Q = 0; Q < bs4->nshell(); Q++) {

                    ints->compute_shell(M, N, P, Q);

                    for (int m = 0, index = 0; m < bs1->shell(M).nfunction(); m++) {
                        for (int n = 0; n < bs2->shell(N).nfunction(); n++) {
                            for (int p = 0; p < bs3->shell(P).nfunction(); p++) {
                                for (int q = 0; q < bs4->shell(Q).nfunction(); q++, index++) {

                                    Ip[(bs1->shell(M).function_index() + m) * nbf2 + bs2->shell(N).function_index() + n]
                                    [(bs3->shell(P).function_index() + p) * nbf4 + bs4->shell(Q).function_index() + q]
                                            = buffer[index];

                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2, nbf3, nbf4};
    I->set_numpy_shape(nshape);

    return I;
}

SharedMatrix MintsHelper::ao_shell_getter(const std::string &label, std::shared_ptr <TwoBodyAOInt> ints, int M, int N, int P, int Q)
{
    int mfxn = basisset_->shell(M).nfunction();
    int nfxn = basisset_->shell(N).nfunction();
    int pfxn = basisset_->shell(P).nfunction();
    int qfxn = basisset_->shell(Q).nfunction();
    SharedMatrix I(new Matrix(label, mfxn * nfxn, pfxn * qfxn));
    double **Ip = I->pointer();
    const double *buffer = ints->buffer();

    ints->compute_shell(M, N, P, Q);

    for (int m = 0, index = 0; m < mfxn; m++) {
        for (int n = 0; n < nfxn; n++) {
            for (int p = 0; p < pfxn; p++) {
                for (int q = 0; q < qfxn; q++, index++) {
                    Ip[m * nfxn + n][p * qfxn + q] = buffer[index];
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{mfxn, nfxn, pfxn, qfxn};
    I->set_numpy_shape(nshape);

    return I;
}

SharedMatrix MintsHelper::ao_erf_eri(double omega)
{
    return ao_helper("AO ERF ERI Integrals", std::shared_ptr<TwoBodyAOInt>(integral_->erf_eri(omega)));
}

SharedMatrix MintsHelper::ao_eri()
{
    std::shared_ptr <TwoBodyAOInt> ints(integral_->eri());
    return ao_helper("AO ERI Tensor", ints);
}

SharedMatrix MintsHelper::ao_eri(std::shared_ptr <BasisSet> bs1,
                                 std::shared_ptr <BasisSet> bs2,
                                 std::shared_ptr <BasisSet> bs3,
                                 std::shared_ptr <BasisSet> bs4)
{
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr <TwoBodyAOInt> ints(intf.eri());
    return ao_helper("AO ERI Tensor", ints);
}

SharedMatrix MintsHelper::ao_eri_shell(int M, int N, int P, int Q)
{
    if (eriInts_ == 0) {
        eriInts_ = std::shared_ptr<TwoBodyAOInt>(integral_->eri());
    }
    return ao_shell_getter("AO ERI Tensor", eriInts_, M, N, P, Q);
}

SharedMatrix MintsHelper::ao_erfc_eri(double omega)
{
    std::shared_ptr <TwoBodyAOInt> ints(integral_->erf_complement_eri(omega));
    return ao_helper("AO ERFC ERI Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12(std::shared_ptr <CorrelationFactor> corr)
{
    std::shared_ptr <TwoBodyAOInt> ints(integral_->f12(corr));
    return ao_helper("AO F12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12(std::shared_ptr <CorrelationFactor> corr,
                                 std::shared_ptr <BasisSet> bs1,
                                 std::shared_ptr <BasisSet> bs2,
                                 std::shared_ptr <BasisSet> bs3,
                                 std::shared_ptr <BasisSet> bs4)
{
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr <TwoBodyAOInt> ints(intf.f12(corr));
    return ao_helper("AO F12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_scaled(std::shared_ptr <CorrelationFactor> corr)
{
    std::shared_ptr <TwoBodyAOInt> ints(integral_->f12_scaled(corr));
    return ao_helper("AO F12 Scaled Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_scaled(std::shared_ptr <CorrelationFactor> corr,
                                        std::shared_ptr <BasisSet> bs1,
                                        std::shared_ptr <BasisSet> bs2,
                                        std::shared_ptr <BasisSet> bs3,
                                        std::shared_ptr <BasisSet> bs4)
{
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr <TwoBodyAOInt> ints(intf.f12_scaled(corr));
    return ao_helper("AO F12 Scaled Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_squared(std::shared_ptr <CorrelationFactor> corr)
{
    std::shared_ptr <TwoBodyAOInt> ints(integral_->f12_squared(corr));
    return ao_helper("AO F12 Squared Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_squared(std::shared_ptr <CorrelationFactor> corr,
                                         std::shared_ptr <BasisSet> bs1,
                                         std::shared_ptr <BasisSet> bs2,
                                         std::shared_ptr <BasisSet> bs3,
                                         std::shared_ptr <BasisSet> bs4)
{
    IntegralFactory intf(bs1, bs2, bs3, bs4);
    std::shared_ptr <TwoBodyAOInt> ints(intf.f12_squared(corr));
    return ao_helper("AO F12 Squared Tensor", ints);
}

SharedMatrix MintsHelper::ao_3coverlap_helper(const std::string &label, std::shared_ptr<ThreeCenterOverlapInt> ints)
{
    std::shared_ptr <BasisSet> bs1 = ints->basis1();
    std::shared_ptr <BasisSet> bs2 = ints->basis2();
    std::shared_ptr <BasisSet> bs3 = ints->basis3();

    int nbf1 = bs1->nbf();
    int nbf2 = bs2->nbf();
    int nbf3 = bs3->nbf();

    SharedMatrix I(new Matrix(label, nbf1 * nbf2, nbf3));
    double **Ip = I->pointer();
    const double *buffer = ints->buffer();


    for (int M = 0; M < bs1->nshell(); M++) {
        for (int N = 0; N < bs2->nshell(); N++) {
            for (int P = 0; P < bs3->nshell(); P++) {

                ints->compute_shell(M, N, P);

                for (int m = 0, index = 0; m < bs1->shell(M).nfunction(); m++) {
                    for (int n = 0; n < bs2->shell(N).nfunction(); n++) {
                        for (int p = 0; p < bs3->shell(P).nfunction(); p++) {
                            Ip[(bs1->shell(M).function_index() + m) * nbf2 + bs2->shell(N).function_index() + n]
                            [(bs3->shell(P).function_index() + p)]
                                    = buffer[index];

                        }
                    }
                }
            }
        }
    }

    // Build numpy and final matrix shape
    std::vector<int> nshape{nbf1, nbf2, nbf3};
    I->set_numpy_shape(nshape);

    return I;
}
SharedMatrix MintsHelper::ao_3coverlap()
{
    std::vector<SphericalTransform> trans;
    for (int i=0; i<basisset_->max_am(); i++){
        trans.push_back(SphericalTransform(i));
    }
    std::shared_ptr<ThreeCenterOverlapInt> ints(new ThreeCenterOverlapInt(trans, basisset_, basisset_, basisset_));
    return ao_3coverlap_helper("AO 3-Center Overlap Tensor", ints);
}

SharedMatrix MintsHelper::ao_3coverlap(std::shared_ptr<BasisSet> bs1,
                                       std::shared_ptr<BasisSet> bs2,
                                       std::shared_ptr<BasisSet> bs3)
{

    int max_am = std::max(std::max(bs1->max_am(), bs2->max_am()), bs3->max_am());
    std::vector<SphericalTransform> trans;
    for (int i=0; i<max_am; i++){
        trans.push_back(SphericalTransform(i));
    }
    std::shared_ptr<ThreeCenterOverlapInt> ints(new ThreeCenterOverlapInt(trans, bs1, bs2, bs3));
    return ao_3coverlap_helper("AO 3-Center Overlap Tensor", ints);
}


SharedMatrix MintsHelper::ao_f12g12(std::shared_ptr <CorrelationFactor> corr)
{
    std::shared_ptr <TwoBodyAOInt> ints(integral_->f12g12(corr));
    return ao_helper("AO F12G12 Tensor", ints);
}

SharedMatrix MintsHelper::ao_f12_double_commutator(std::shared_ptr <CorrelationFactor> corr)
{
    std::shared_ptr <TwoBodyAOInt> ints(integral_->f12_double_commutator(corr));
    return ao_helper("AO F12 Double Commutator Tensor", ints);
}

SharedMatrix MintsHelper::mo_erf_eri(double omega, SharedMatrix C1, SharedMatrix C2,
                                     SharedMatrix C3, SharedMatrix C4)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_erf_eri(omega), C1, C2, C3, C4);
    mo_ints->set_name("MO ERF ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_erfc_eri(double omega, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_erfc_eri(omega), C1, C2, C3, C4);
    mo_ints->set_name("MO ERFC ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12(std::shared_ptr <CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_f12(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12_squared(std::shared_ptr <CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_f12_squared(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Squared Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12g12(std::shared_ptr <CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_f12g12(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12G12 Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_f12_double_commutator(std::shared_ptr <CorrelationFactor> corr, SharedMatrix C1, SharedMatrix C2, SharedMatrix C3, SharedMatrix C4)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_f12_double_commutator(corr), C1, C2, C3, C4);
    mo_ints->set_name("MO F12 Double Commutator Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_eri(SharedMatrix C1, SharedMatrix C2,
                                 SharedMatrix C3, SharedMatrix C4)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_eri(), C1, C2, C3, C4);
    mo_ints->set_name("MO ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_erf_eri(double omega, SharedMatrix Co, SharedMatrix Cv)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_erf_eri(omega), Co, Cv);
    mo_ints->set_name("MO ERF ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_eri(SharedMatrix Co, SharedMatrix Cv)
{
    SharedMatrix mo_ints = mo_eri_helper(ao_eri(), Co, Cv);
    mo_ints->set_name("MO ERI Tensor");
    return mo_ints;
}

SharedMatrix MintsHelper::mo_eri_helper(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2,
                                        SharedMatrix C3, SharedMatrix C4)
{
    int nso = basisset_->nbf();
    int n1 = C1->colspi()[0];
    int n2 = C2->colspi()[0];
    int n3 = C3->colspi()[0];
    int n4 = C4->colspi()[0];

    double **C1p = C1->pointer();
    double **C2p = C2->pointer();
    double **C3p = C3->pointer();
    double **C4p = C4->pointer();

    double **Isop = Iso->pointer();
    SharedMatrix I2(new Matrix("MO ERI Tensor", n1 * nso, nso * nso));
    double **I2p = I2->pointer();

    C_DGEMM('T', 'N', n1, nso * (ULI) nso * nso, nso, 1.0, C1p[0], n1, Isop[0], nso * (ULI) nso * nso, 0.0, I2p[0], nso * (ULI) nso * nso);

    Iso.reset();
    SharedMatrix I3(new Matrix("MO ERI Tensor", n1 * nso, nso * n3));
    double **I3p = I3->pointer();

    C_DGEMM('N', 'N', n1 * (ULI) nso * nso, n3, nso, 1.0, I2p[0], nso, C3p[0], n3, 0.0, I3p[0], n3);

    I2.reset();
    SharedMatrix I4(new Matrix("MO ERI Tensor", nso * n1, n3 * nso));
    double **I4p = I4->pointer();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int m = 0; m < nso; m++) {
                for (int n = 0; n < nso; n++) {
                    I4p[m * n1 + i][j * nso + n] = I3p[i * nso + m][n * n3 + j];
                }
            }
        }
    }

    I3.reset();
    SharedMatrix I5(new Matrix("MO ERI Tensor", n2 * n1, n3 * nso));
    double **I5p = I5->pointer();

    C_DGEMM('T', 'N', n2, n1 * (ULI) n3 * nso, nso, 1.0, C2p[0], n2, I4p[0], n1 * (ULI) n3 * nso, 0.0, I5p[0], n1 * (ULI) n3 * nso);

    I4.reset();
    SharedMatrix I6(new Matrix("MO ERI Tensor", n2 * n1, n3 * n4));
    double **I6p = I6->pointer();

    C_DGEMM('N', 'N', n2 * (ULI) n1 * n3, n4, nso, 1.0, I5p[0], nso, C4p[0], n4, 0.0, I6p[0], n4);

    I5.reset();
    SharedMatrix Imo(new Matrix("MO ERI Tensor", n1 * n2, n3 * n4));
    double **Imop = Imo->pointer();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int a = 0; a < n2; a++) {
                for (int b = 0; b < n4; b++) {
                    Imop[i * n2 + a][j * n4 + b] = I6p[a * n1 + i][j * n4 + b];
                }
            }
        }
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{n1, n2, n3, n4};
    Imo->set_numpy_shape(nshape);

    return Imo;
}

SharedMatrix MintsHelper::mo_eri_helper(SharedMatrix Iso, SharedMatrix Co, SharedMatrix Cv)
{
    int nso = basisset_->nbf();
    int nocc = Co->colspi()[0];
    int nvir = Cv->colspi()[0];

    double **Cop = Co->pointer();
    double **Cvp = Cv->pointer();

    double **Isop = Iso->pointer();
    SharedMatrix I2(new Matrix("MO ERI Tensor", nocc * nso, nso * nso));
    double **I2p = I2->pointer();

    C_DGEMM('T', 'N', nocc, nso * (ULI) nso * nso, nso, 1.0, Cop[0], nocc, Isop[0], nso * (ULI) nso * nso, 0.0, I2p[0], nso * (ULI) nso * nso);

    Iso.reset();
    SharedMatrix I3(new Matrix("MO ERI Tensor", nocc * nso, nso * nocc));
    double **I3p = I3->pointer();

    C_DGEMM('N', 'N', nocc * (ULI) nso * nso, nocc, nso, 1.0, I2p[0], nso, Cop[0], nocc, 0.0, I3p[0], nocc);

    I2.reset();
    SharedMatrix I4(new Matrix("MO ERI Tensor", nso * nocc, nocc * nso));
    double **I4p = I4->pointer();

    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int m = 0; m < nso; m++) {
                for (int n = 0; n < nso; n++) {
                    I4p[m * nocc + i][j * nso + n] = I3p[i * nso + m][n * nocc + j];
                }
            }
        }
    }

    I3.reset();
    SharedMatrix I5(new Matrix("MO ERI Tensor", nvir * nocc, nocc * nso));
    double **I5p = I5->pointer();

    C_DGEMM('T', 'N', nvir, nocc * (ULI) nocc * nso, nso, 1.0, Cvp[0], nvir, I4p[0], nocc * (ULI) nocc * nso, 0.0, I5p[0], nocc * (ULI) nocc * nso);

    I4.reset();
    SharedMatrix I6(new Matrix("MO ERI Tensor", nvir * nocc, nocc * nvir));
    double **I6p = I6->pointer();

    C_DGEMM('N', 'N', nvir * (ULI) nocc * nocc, nvir, nso, 1.0, I5p[0], nso, Cvp[0], nvir, 0.0, I6p[0], nvir);

    I5.reset();
    SharedMatrix Imo(new Matrix("MO ERI Tensor", nocc * nvir, nocc * nvir));
    double **Imop = Imo->pointer();

    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = 0; a < nvir; a++) {
                for (int b = 0; b < nvir; b++) {
                    Imop[i * nvir + a][j * nvir + b] = I6p[a * nocc + i][j * nvir + b];
                }
            }
        }
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{nocc, nvir, nocc, nvir};
    Imo->set_numpy_shape(nshape);

    return Imo;
}

SharedMatrix MintsHelper::mo_spin_eri(SharedMatrix Co, SharedMatrix Cv)
{
    int n1 = Co->colspi()[0];
    int n2 = Cv->colspi()[0];
    SharedMatrix mo_ints = mo_eri_helper(ao_eri(), Co, Cv);
    SharedMatrix mo_spin_ints = mo_spin_eri_helper(mo_ints, n1, n2);
    mo_ints.reset();
    mo_spin_ints->set_name("MO Spin ERI Tensor");
    return mo_spin_ints;
}

SharedMatrix MintsHelper::mo_spin_eri_helper(SharedMatrix Iso, int n1, int n2)
{
    int n12 = n1 * 2;
    int n22 = n2 * 2;

    double **Isop = Iso->pointer();
    SharedMatrix Ispin(new Matrix("MO ERI Tensor", 4 * n1 * n1, 4 * n2 * n2));
    double **Ispinp = Ispin->pointer();

    double first, second;
    int mask1, mask2;
    for (int i = 0; i < n12; i++) {
        for (int j = 0; j < n12; j++) {
            for (int k = 0; k < n22; k++) {
                for (int l = 0; l < n22; l++) {
                    mask1 = (i % 2 == k % 2) * (j % 2 == l % 2);
                    mask2 = (i % 2 == l % 2) * (j % 2 == k % 2);

                    first = Isop[i / 2 * n2 + k / 2][j / 2 * n2 + l / 2];
                    second = Isop[i / 2 * n2 + l / 2][j / 2 * n2 + k / 2];
                    Ispinp[i * n12 + j][k * n22 + l] = first * mask1 - second * mask2;
                }
            }
        }
    }
    // Build numpy and final matrix shape
    std::vector<int> nshape{n12, n12, n22, n22};
    Ispin->set_numpy_shape(nshape);

    return Ispin;
}

SharedMatrix MintsHelper::so_overlap()
{
    if (factory_->nirrep() == 1) {
        SharedMatrix ret = ao_overlap();
        ret->set_name(PSIF_SO_S);
        return ret;
    } else {
        SharedMatrix overlap_mat(factory_->create_matrix(PSIF_SO_S));
        overlap_mat->apply_symmetry(ao_overlap(), petite_list()->aotoso());
        return overlap_mat;
    }
}

SharedMatrix MintsHelper::so_kinetic()
{
    if (factory_->nirrep() == 1) {
        SharedMatrix ret = ao_kinetic();
        ret->set_name(PSIF_SO_T);
        return ret;
    } else {
        SharedMatrix kinetic_mat(factory_->create_matrix(PSIF_SO_T));
        kinetic_mat->apply_symmetry(ao_kinetic(), petite_list()->aotoso());
        return kinetic_mat;
    }
}

SharedMatrix MintsHelper::so_potential(bool include_perturbations)
{
    // No symmetry
    SharedMatrix potential_mat;
    if (factory_->nirrep() == 1) {
        potential_mat = ao_potential();
        potential_mat->set_name(PSIF_SO_V);
    } else {
        potential_mat = factory_->create_shared_matrix(PSIF_SO_V);
        potential_mat->apply_symmetry(ao_potential(), petite_list()->aotoso());
    }

    if (integral_->hasECP()) {
        outfile->Printf("Adding ECP terms to potential..\n");
        std::shared_ptr <OneBodySOInt> ECP(integral_->so_ecp());
        SharedMatrix ecp_mat(new Matrix("AO-basis ECP Ints", potential_mat->nrow(), potential_mat->ncol()));
        ECP->compute(ecp_mat);
        potential_mat->print();
        ecp_mat->print();
        potential_mat->add(ecp_mat);
    }

    // Handle addition of any perturbations here and not in SCF code.
    if (include_perturbations) {
        if (options_.get_bool("PERTURB_H")) {
            std::string perturb_with = options_.get_str("PERTURB_WITH");
            Vector3 lambda(0.0, 0.0, 0.0);

            if (perturb_with == "DIPOLE_X")
                lambda[0] = options_.get_double("PERTURB_MAGNITUDE");
            else if (perturb_with == "DIPOLE_Y")
                lambda[1] = options_.get_double("PERTURB_MAGNITUDE");
            else if (perturb_with == "DIPOLE_Z")
                lambda[2] = options_.get_double("PERTURB_MAGNITUDE");
            else if (perturb_with == "DIPOLE") {
                if(options_["PERTURB_DIPOLE"].size() !=3)
                    throw PSIEXCEPTION("The PERTURB dipole should have exactly three floating point numbers.");
                for(int n = 0; n < 3; ++n)
                    lambda[n] = options_["PERTURB_DIPOLE"][n].to_double();
            } else {
                outfile->Printf("  MintsHelper doesn't understand the requested perturbation, might be done in SCF.");
            }

            OperatorSymmetry msymm(1, molecule_, integral_, factory_);
            std::vector <SharedMatrix> dipoles = msymm.create_matrices("Dipole");
            OneBodySOInt *so_dipole = integral_->so_dipole();
            so_dipole->compute(dipoles);

            if (lambda[0] != 0.0) {
                if (msymm.component_symmetry(0) != 0) {
                    outfile->Printf("  WARNING: Requested mu(x) perturbation, but mu(x) is not symmetric.\n");
                } else {
                    outfile->Printf("  Perturbing V by %f mu(x).\n", lambda[0]);
                    dipoles[0]->scale(lambda[0]);
                    potential_mat->add(dipoles[0]);
                }
            }
            if (lambda[1] != 0.0) {
                if (msymm.component_symmetry(1) != 0) {
                    outfile->Printf("  WARNING: Requested mu(y) perturbation, but mu(y) is not symmetric.\n");
                } else {
                    outfile->Printf("  Perturbing V by %f mu(y).\n", lambda[1]);
                    dipoles[1]->scale(lambda[1]);
                    potential_mat->add(dipoles[1]);
                }
            }
            if (lambda[2] != 0.0) {
                if (msymm.component_symmetry(2) != 0) {
                    outfile->Printf("  WARNING: Requested mu(z) perturbation, but mu(z) is not symmetric.\n");
                } else {
                    outfile->Printf("  Perturbing V by %f mu(z).\n", lambda[2]);
                    dipoles[2]->scale(lambda[2]);
                    potential_mat->add(dipoles[2]);
                }
            }
        }

        if (options_.get_str("RELATIVISTIC") == "DKH") {
            int dkh_order = options_.get_int("DKH_ORDER");
            SharedMatrix dkh = so_dkh(dkh_order);

            outfile->Printf("    Adding Douglas-Kroll-Hess corrections to the potential integrals.\n");

            potential_mat->add(dkh);
        }
    }

    return potential_mat;
}

std::vector <SharedMatrix> MintsHelper::so_dipole()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(1, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> dipole = msymm.create_matrices("SO Dipole");

    std::shared_ptr <OneBodySOInt> ints(integral_->so_dipole());
    ints->compute(dipole);

    return dipole;
}

std::vector <SharedMatrix> MintsHelper::so_quadrupole()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> quadrupole = msymm.create_matrices("SO Quadrupole");

    std::shared_ptr <OneBodySOInt> ints(integral_->so_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}

std::vector <SharedMatrix> MintsHelper::so_traceless_quadrupole()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> quadrupole = msymm.create_matrices("SO Traceless Quadrupole");

    std::shared_ptr <OneBodySOInt> ints(integral_->so_traceless_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}

std::vector <SharedMatrix> MintsHelper::so_nabla()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(OperatorSymmetry::P, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> nabla = msymm.create_matrices("SO Nabla");

    std::shared_ptr <OneBodySOInt> ints(integral_->so_nabla());
    ints->compute(nabla);

    return nabla;
}

std::vector <SharedMatrix> MintsHelper::so_angular_momentum()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(OperatorSymmetry::L, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> am = msymm.create_matrices("SO Angular Momentum");

    std::shared_ptr <OneBodySOInt> ints(integral_->so_angular_momentum());
    ints->compute(am);

    return am;
}

std::vector <SharedMatrix> MintsHelper::ao_angular_momentum()
{
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> angmom;

    angmom.push_back(SharedMatrix(new Matrix("AO Lx", basisset_->nbf(), basisset_->nbf())));
    angmom.push_back(SharedMatrix(new Matrix("AO Ly", basisset_->nbf(), basisset_->nbf())));
    angmom.push_back(SharedMatrix(new Matrix("AO Lz", basisset_->nbf(), basisset_->nbf())));

    std::shared_ptr <OneBodyAOInt> ints(integral_->ao_angular_momentum());
    ints->compute(angmom);

    return angmom;
}

std::vector <SharedMatrix> MintsHelper::ao_dipole()
{
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> dipole;

    dipole.push_back(SharedMatrix(new Matrix("AO Mux", basisset_->nbf(), basisset_->nbf())));
    dipole.push_back(SharedMatrix(new Matrix("AO Muy", basisset_->nbf(), basisset_->nbf())));
    dipole.push_back(SharedMatrix(new Matrix("AO Muz", basisset_->nbf(), basisset_->nbf())));

    std::shared_ptr <OneBodyAOInt> ints(integral_->ao_dipole());
    ints->compute(dipole);

    return dipole;
}

std::vector <SharedMatrix> MintsHelper::ao_nabla()
{
    // Create a vector of matrices with the proper symmetry
    std::vector <SharedMatrix> nabla;

    nabla.push_back(SharedMatrix(new Matrix("AO Px", basisset_->nbf(), basisset_->nbf())));
    nabla.push_back(SharedMatrix(new Matrix("AO Py", basisset_->nbf(), basisset_->nbf())));
    nabla.push_back(SharedMatrix(new Matrix("AO Pz", basisset_->nbf(), basisset_->nbf())));

    std::shared_ptr <OneBodyAOInt> ints(integral_->ao_nabla());
    ints->compute(nabla);

    return nabla;
}

std::shared_ptr <CdSalcList> MintsHelper::cdsalcs(int needed_irreps,
                                                  bool project_out_translations,
                                                  bool project_out_rotations)
{
    return std::shared_ptr<CdSalcList>(new CdSalcList(molecule_, factory_,
                                                      needed_irreps,
                                                      project_out_translations,
                                                      project_out_rotations));
}

SharedMatrix MintsHelper::mo_transform(SharedMatrix Iso, SharedMatrix C1, SharedMatrix C2,
                                       SharedMatrix C3, SharedMatrix C4)
{
    // Attempts to transform integrals in the most efficient manner. Will transpose left, right
    // and left_right where left and right are (left|right) indices. Does not consider the fimal
    // perturbation eg (12|34) -> (13|24), therefore integrals of type (oo|vv) will not be computed
    // in the optimal order. However, the first transformed index is guaranteed to be the smallest.

    int nso = C1->rowspi()[0];

    // Check C dimensions
    int dim_check = 0;
    dim_check += (C2->rowspi()[0] != nso);
    dim_check += (C3->rowspi()[0] != nso);
    dim_check += (C4->rowspi()[0] != nso);

    if (dim_check) {
        throw PSIEXCEPTION("MO Transform: Eigenvector lengths of the C matrices are not identical.");
    }

    if ((Iso->nirrep()) > 1) {
        throw PSIEXCEPTION("MO Transform: The ERI has more than one irrep.");
    }

    // Make sure I is square and of correction dimesion
    int Irows = Iso->rowspi()[0];
    int Icols = Iso->colspi()[0];

    if (((nso * nso) != Irows) || ((nso * nso) != Icols)) {
        throw PSIEXCEPTION("MO Transform: ERI shape does match that of the C matrices.");
    }

    int n1 = C1->colspi()[0];
    int n2 = C2->colspi()[0];
    int n3 = C3->colspi()[0];
    int n4 = C4->colspi()[0];
    std::vector<int> nshape{n1, n2, n3, n4};

    double **C1p = C1->pointer();
    double **C2p = C2->pointer();
    double **C3p = C3->pointer();
    double **C4p = C4->pointer();

    int shape_left = n1 * n2;
    int shape_right = n3 * n4;

    // Transform smallest indices first
    bool transpose_left = (n1 > n2);
    bool transpose_right = (n3 > n4);
    bool transpose_left_right = ((transpose_left ? n2 : n1) > (transpose_right ? n4 : n3));

    int tmp_n;
    double **tmp_p;
    if (transpose_left) {
        tmp_n = n1;
        n1 = n2;
        n2 = tmp_n;
        tmp_p = C1p;
        C1p = C2p;
        C2p = tmp_p;
    }
    if (transpose_right) {
        tmp_n = n3;
        n3 = n4;
        n4 = tmp_n;
        tmp_p = C3p;
        C3p = C4p;
        C4p = tmp_p;
    }
    if (transpose_left_right) {
        tmp_n = n1;
        n1 = n3;
        n3 = tmp_n;
        tmp_n = n2;
        n2 = n4;
        n4 = tmp_n;
        tmp_p = C1p;
        C1p = C3p;
        C3p = tmp_p;
        tmp_p = C2p;
        C2p = C4p;
        C4p = tmp_p;
    }

    double **Isop = Iso->pointer();
    SharedMatrix I2(new Matrix("MO ERI Tensor", n1 * nso, nso * nso));
    double **I2p = I2->pointer();

    C_DGEMM('T', 'N', n1, nso * (ULI) nso * nso, nso, 1.0, C1p[0], n1, Isop[0], nso * (ULI) nso * nso, 0.0, I2p[0], nso * (ULI) nso * nso);

    Iso.reset();
    SharedMatrix I3(new Matrix("MO ERI Tensor", n1 * nso, nso * n3));
    double **I3p = I3->pointer();

    C_DGEMM('N', 'N', n1 * (ULI) nso * nso, n3, nso, 1.0, I2p[0], nso, C3p[0], n3, 0.0, I3p[0], n3);

    I2.reset();
    SharedMatrix I4(new Matrix("MO ERI Tensor", nso * n1, n3 * nso));
    double **I4p = I4->pointer();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int m = 0; m < nso; m++) {
                for (int n = 0; n < nso; n++) {
                    I4p[m * n1 + i][j * nso + n] = I3p[i * nso + m][n * n3 + j];
                }
            }
        }
    }

    I3.reset();
    SharedMatrix I5(new Matrix("MO ERI Tensor", n2 * n1, n3 * nso));
    double **I5p = I5->pointer();

    C_DGEMM('T', 'N', n2, n1 * (ULI) n3 * nso, nso, 1.0, C2p[0], n2, I4p[0], n1 * (ULI) n3 * nso, 0.0, I5p[0], n1 * (ULI) n3 * nso);

    I4.reset();
    SharedMatrix I6(new Matrix("MO ERI Tensor", n2 * n1, n3 * n4));
    double **I6p = I6->pointer();

    C_DGEMM('N', 'N', n2 * (ULI) n1 * n3, n4, nso, 1.0, I5p[0], nso, C4p[0], n4, 0.0, I6p[0], n4);

    I5.reset();
    SharedMatrix Imo(new Matrix("MO ERI Tensor", shape_left, shape_right));
    double **Imop = Imo->pointer();

    // Currently 2143, need to transform back
    int left, right;
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int a = 0; a < n2; a++) {

                // Tranpose left
                if (transpose_left) {
                    left = a * n1 + i;
                } else {
                    left = i * n2 + a;
                }

                for (int b = 0; b < n4; b++) {
                    right = j * n4 + b;

                    // Transpose right
                    if (transpose_right) {
                        right = b * n3 + j;
                    } else {
                        right = j * n4 + b;
                    }

                    // Transpose left_right
                    if (transpose_left_right) {
                        Imop[right][left] = I6p[a * n1 + i][j * n4 + b];
                    } else {
                        Imop[left][right] = I6p[a * n1 + i][j * n4 + b];
                    }
                }
            }
        }
    }
    // Build numpy and final matrix shape
    Imo->set_numpy_shape(nshape);

    return Imo;
}

void MintsHelper::play()
{
}

} // namespace psi
