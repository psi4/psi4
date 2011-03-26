#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>

#include <psifiles.h>
#include <libpsio/psio.hpp>
#include <libiwl/iwl.hpp>
#include <libciomr/libciomr.h>
#include "mints.h"

#include <libqt/qt.h>

#include <psi4-dec.h>

using namespace boost;

namespace psi {

void MintsHelper::init_helper()
{
    psio_ = shared_ptr<PSIO>(new PSIO());
    molecule_ = shared_ptr<Molecule>(Process::environment.molecule());

    if (molecule_.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Make sure molecule is valid.
    molecule_->update_geometry();

    // Print the molecule.
    if (print_)
        molecule_->print();

    // Read in the basis set
    shared_ptr<BasisSetParser> parser (new Gaussian94BasisSetParser());
    basisset_ = shared_ptr<BasisSet>(BasisSet::construct(parser, molecule_, "BASIS"));

    // Print the basis set
    if (print_)
        basisset_->print_detail();

    // Create integral factory
    integral_ = shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));

    // Get the SO basis object.
    sobasis_ = shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral_));

    // Obtain dimensions from the sobasis
    const Dimension dimension = sobasis_->dimension();

    // Create a matrix factory and initialize it
    factory_ = shared_ptr<MatrixFactory>(new MatrixFactory());
    factory_->init_with(dimension, dimension);
}

MintsHelper::MintsHelper(Options & options, int print): options_(options), print_(print)
{
    init_helper();
}

MintsHelper::MintsHelper() : options_(Process::environment.options), print_(1)
{
    init_helper();
}

MintsHelper::~MintsHelper()
{
}

void MintsHelper::integrals()
{
    timer_init();

    fprintf(outfile, " MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // Get ERI object
    shared_ptr<TwoBodyAOInt> tb(integral_->eri());
    shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, integral_));

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:                %4d\n", molecule_->natom());
    fprintf(outfile, "      Number of AO shells:            %4d\n", basisset_->nshell());
    fprintf(outfile, "      Number of SO shells:            %4d\n", sobasis_->nshell());
    fprintf(outfile, "      Number of primitives:           %4d\n", basisset_->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals:      %4d\n", basisset_->nao());
    fprintf(outfile, "      Number of basis functions:      %4d\n\n", basisset_->nbf());
    fprintf(outfile, "      Number of irreps:               %4d\n", sobasis_->nirrep());
    fprintf(outfile, "      Number of functions per irrep: [");
    for (int i=0; i<sobasis_->nirrep(); ++i) {
        fprintf(outfile, "%4d ", sobasis_->nfunction_in_irrep(i));
    }
    fprintf(outfile, "]\n\n");

    // Compute and dump one-electron SO integrals.

    // Overlap
    shared_ptr<OneBodySOInt> overlap(integral_->so_overlap());
    shared_ptr<Matrix>       overlap_mat(factory_->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);

    // Kinetic
    shared_ptr<OneBodySOInt> kinetic(integral_->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(factory_->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(psio_, PSIF_OEI);

    // Potential
    shared_ptr<OneBodySOInt> potential(integral_->so_potential());
    shared_ptr<Matrix>       potential_mat(factory_->create_matrix(PSIF_SO_V));
    potential->compute(potential_mat);
    potential_mat->save(psio_, PSIF_OEI);

    if (Process::environment.options.get_int("PRINT") > 3) {
        overlap_mat->print();
        kinetic_mat->print();
        potential_mat->print();
    }

    // Open the IWL buffer where we will store the integrals.
    IWL ERIOUT(psio_.get(), PSIF_SO_TEI, 0.0, 0, 0);
    IWLWriter writer(ERIOUT);

    // Let the user know what we're doing.
    fprintf(outfile, "      Computing integrals..."); fflush(outfile);

    SOShellCombinationsIterator shellIter(sobasis_, sobasis_, sobasis_, sobasis_);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        eri->compute_shell(shellIter, writer);
    }

    // Flush out buffers.
    ERIOUT.flush(1);

    // We just did all this work to create the file, let's keep it around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    fprintf(outfile, "done\n\n");
    fprintf(outfile, "      Computed %lu non-zero integrals.\n\n", writer.count());
    timer_done();
}

void MintsHelper::one_electron_integrals()
{
    fprintf(outfile, " OEINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:                %4d\n", molecule_->natom());
    fprintf(outfile, "      Number of AO shells:            %4d\n", basisset_->nshell());
    fprintf(outfile, "      Number of SO shells:            %4d\n", sobasis_->nshell());
    fprintf(outfile, "      Number of primitives:           %4d\n", basisset_->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals:      %4d\n", basisset_->nao());
    fprintf(outfile, "      Number of basis functions:      %4d\n\n", basisset_->nbf());
    fprintf(outfile, "      Number of irreps:               %4d\n", sobasis_->nirrep());
    fprintf(outfile, "      Number of functions per irrep: [");
    for (int i=0; i<sobasis_->nirrep(); ++i) {
        fprintf(outfile, "%4d ", sobasis_->nfunction_in_irrep(i));
    }
    fprintf(outfile, "]\n\n");

    // Compute and dump one-electron SO integrals.

    // Overlap
    shared_ptr<OneBodySOInt> overlap(integral_->so_overlap());
    shared_ptr<Matrix>       overlap_mat(factory_->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);

    // Kinetic
    shared_ptr<OneBodySOInt> kinetic(integral_->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(factory_->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(psio_, PSIF_OEI);

    // Potential
    shared_ptr<OneBodySOInt> potential(integral_->so_potential());
    shared_ptr<Matrix>       potential_mat(factory_->create_matrix(PSIF_SO_V));
    potential->compute(potential_mat);
    potential_mat->save(psio_, PSIF_OEI);
}

void MintsHelper::integral_gradients()
{
    throw FeatureNotImplemented("libmints", "MintsHelper::integral_derivatives", __FILE__, __LINE__);
}

void MintsHelper::integral_hessians()
{
    throw FeatureNotImplemented("libmints", "MintsHelper::integral_hessians", __FILE__, __LINE__);
}

shared_ptr<Matrix> MintsHelper::ao_overlap()
{
    // Overlap
    shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    shared_ptr<Matrix>       overlap_mat(factory_->create_matrix(PSIF_AO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);
    return overlap_mat;
}

shared_ptr<Matrix> MintsHelper::ao_kinetic()
{
    shared_ptr<OneBodyAOInt> T(integral_->ao_kinetic());
    shared_ptr<Matrix>       kinetic_mat(factory_->create_matrix());
    T->compute(kinetic_mat);
    kinetic_mat->save(psio_, PSIF_OEI);
    return kinetic_mat;
}

shared_ptr<Matrix> MintsHelper::ao_potential()
{
    shared_ptr<OneBodyAOInt> V(integral_->ao_potential());
    shared_ptr<Matrix>       potential_mat(factory_->create_matrix());
    V->compute(potential_mat);
    potential_mat->save(psio_, PSIF_OEI);
    return potential_mat;
}

shared_ptr<Matrix> MintsHelper::ao_erf_eri(double omega, double alpha, double beta)
{
    int nbf = basisset_->nbf();
    shared_ptr<Matrix> I(new Matrix("AO ERF ERI Integrals", nbf*nbf, nbf*nbf));

    shared_ptr<TwoBodyAOInt> ints(integral_->erf_eri(omega, alpha, beta));
    double** Ip = I->pointer();
    const double* buffer = ints->buffer();

    for (int M = 0; M < basisset_->nshell(); M++) {
        for (int N = 0; N < basisset_->nshell(); N++) {
            for (int P = 0; P < basisset_->nshell(); P++) {
                for (int Q = 0; Q < basisset_->nshell(); Q++) {

                    ints->compute_shell(M,N,P,Q);

                    for (int m = 0, index = 0; m < basisset_->shell(M)->nfunction(); m++) {
                        for (int n = 0; n < basisset_->shell(N)->nfunction(); n++) {
                            for (int p = 0; p < basisset_->shell(P)->nfunction(); p++) {
                                for (int q = 0; q < basisset_->shell(Q)->nfunction(); q++, index++) {

                                    Ip[(basisset_->shell(M)->function_index() + m)*nbf + basisset_->shell(N)->function_index() + n]
                                            [(basisset_->shell(P)->function_index() + p)*nbf + basisset_->shell(Q)->function_index() + q]
                                            = buffer[index];

                                } } } }

                } } } }


    return I;
}

shared_ptr<Matrix> MintsHelper::ao_eri()
{
    int nbf = basisset_->nbf();
    shared_ptr<Matrix> I(new Matrix("AO ERI Integrals", nbf*nbf, nbf*nbf));

    shared_ptr<TwoBodyAOInt> ints(integral_->eri());
    double** Ip = I->pointer();
    const double* buffer = ints->buffer();

    for (int M = 0; M < basisset_->nshell(); M++) {
        for (int N = 0; N < basisset_->nshell(); N++) {
            for (int P = 0; P < basisset_->nshell(); P++) {
                for (int Q = 0; Q < basisset_->nshell(); Q++) {

                    ints->compute_shell(M,N,P,Q);

                    for (int m = 0, index = 0; m < basisset_->shell(M)->nfunction(); m++) {
                        for (int n = 0; n < basisset_->shell(N)->nfunction(); n++) {
                            for (int p = 0; p < basisset_->shell(P)->nfunction(); p++) {
                                for (int q = 0; q < basisset_->shell(Q)->nfunction(); q++, index++) {

                                    Ip[(basisset_->shell(M)->function_index() + m)*nbf + basisset_->shell(N)->function_index() + n]
                                            [(basisset_->shell(P)->function_index() + p)*nbf + basisset_->shell(Q)->function_index() + q]
                                            = buffer[index];

                                } } } }

                } } } }


    return I;
}

shared_ptr<Matrix> MintsHelper::so_overlap()
{
    shared_ptr<OneBodySOInt> S(integral_->so_overlap());
    shared_ptr<Matrix>       overlap_mat(factory_->create_matrix(PSIF_SO_S));
    S->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);
    return overlap_mat;
}

shared_ptr<Matrix> MintsHelper::so_kinetic()
{
    shared_ptr<OneBodySOInt> T(integral_->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(factory_->create_matrix(PSIF_SO_T));
    T->compute(kinetic_mat);
    kinetic_mat->save(psio_, PSIF_OEI);
    return kinetic_mat;
}

shared_ptr<Matrix> MintsHelper::so_potential()
{
    shared_ptr<OneBodySOInt> V(integral_->so_potential());
    shared_ptr<Matrix>       potential_mat(factory_->create_matrix(PSIF_SO_V));
    V->compute(potential_mat);
    potential_mat->save(psio_, PSIF_OEI);
    return potential_mat;
}

std::vector<shared_ptr<Matrix> > MintsHelper::so_dipole()
{
    // The matrix factory can create matrices of the correct dimensions...
    MultipoleSymmetry msymm(1, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> dipole = msymm.create_matrices("SO Dipole");

    shared_ptr<OneBodySOInt> ints(integral_->so_dipole());
    ints->compute(dipole);

    return dipole;
}
std::vector<shared_ptr<Matrix> > MintsHelper::so_quadrupole()
{
    // The matrix factory can create matrices of the correct dimensions...
    MultipoleSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole = msymm.create_matrices("SO Quadrupole");

    shared_ptr<OneBodySOInt> ints(integral_->so_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}
std::vector<shared_ptr<Matrix> > MintsHelper::so_traceless_quadrupole()
{
    // The matrix factory can create matrices of the correct dimensions...
    MultipoleSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole = msymm.create_matrices("SO Traceless Quadrupole");

    shared_ptr<OneBodySOInt> ints(integral_->so_traceless_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}
std::vector<shared_ptr<Matrix> > MintsHelper::so_nabla()
{
    // The matrix factory can create matrices of the correct dimensions...
    MultipoleSymmetry msymm(1, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> nabla = msymm.create_matrices("SO Nabla");

    // Jet, here's a stub

    return nabla;
}
std::vector<shared_ptr<Matrix> > MintsHelper::so_angular_momentum()
{
    // The matrix factory can create matrices of the correct dimensions...
    MultipoleSymmetry msymm(1, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> angmom = msymm.create_matrices("SO Angular Momentum");

    // Jet, here's a stub

    return angmom;
}

}
