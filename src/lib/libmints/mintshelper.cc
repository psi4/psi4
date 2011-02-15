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

MintsHelper::MintsHelper(Options & options): options_(options)
{
}

MintsHelper::MintsHelper() : options_(Process::environment.options)
{
}

MintsHelper::~MintsHelper()
{
}

void MintsHelper::integrals()
{
    timer_init();
    shared_ptr<PSIO> psio(new PSIO());

    fprintf(outfile, " MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // We'll only be working with the active molecule.
    shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Make sure molecule is valid.
    molecule->update_geometry();

    // Print the molecule.
    molecule->print();

    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basisset = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));

    // Print the basis set
    basisset->print_detail();

    // Create integral factory
    shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset, basisset, basisset, basisset));

    // Get the SO basis object.
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basisset, integral));

    // Get ERI object
    shared_ptr<TwoBodyAOInt> tb(integral->eri());
    shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, integral));

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:                %4d\n", molecule->natom());
    fprintf(outfile, "      Number of AO shells:            %4d\n", basisset->nshell());
    fprintf(outfile, "      Number of SO shells:            %4d\n", sobasis->nshell());
    fprintf(outfile, "      Number of primitives:           %4d\n", basisset->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals:      %4d\n", basisset->nao());
    fprintf(outfile, "      Number of basis functions:      %4d\n\n", basisset->nbf());
    fprintf(outfile, "      Number of irreps:               %4d\n", sobasis->nirrep());
    fprintf(outfile, "      Number of functions per irrep: [");
    for (int i=0; i<sobasis->nirrep(); ++i) {
        fprintf(outfile, "%4d ", sobasis->nfunction_in_irrep(i));
    }
    fprintf(outfile, "]\n\n");

    // Compute and dump one-electron SO integrals.
    // Obtain dimensions from the sobasis
    const Dimension dimension = sobasis->dimension();

    // Create a matrix factory and initialize it
    shared_ptr<MatrixFactory> factory(new MatrixFactory());
    factory->init_with(dimension, dimension);

    // Overlap
    shared_ptr<OneBodySOInt> overlap(integral->so_overlap());
    shared_ptr<Matrix>       overlap_mat(factory->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio, PSIF_OEI);

    // Kinetic
    shared_ptr<OneBodySOInt> kinetic(integral->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(factory->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(psio, PSIF_OEI);

    // Potential
    shared_ptr<OneBodySOInt> potential(integral->so_potential());
    shared_ptr<Matrix>       potential_mat(factory->create_matrix(PSIF_SO_V));
    potential->compute(potential_mat);
    potential_mat->save(psio, PSIF_OEI);

    // Open the IWL buffer where we will store the integrals.
    IWL ERIOUT(psio.get(), PSIF_SO_TEI, 0.0, 0, 0);
    IWLWriter writer(ERIOUT);

    // Let the user know what we're doing.
    fprintf(outfile, "      Computing integrals..."); fflush(outfile);

    SOShellCombinationsIterator shellIter(sobasis, sobasis, sobasis, sobasis);
    for (shellIter.first(); shellIter.is_done() == false; shellIter.next()) {
        eri->compute_shell(shellIter, writer);
    }

    // Flush out buffers.
    ERIOUT.flush(1);

    // We just did all this work to create the file, let's keep it around
    ERIOUT.set_keep_flag(true);
    ERIOUT.close();

    fprintf(outfile, "done\n\n"); fflush(outfile);

    fprintf(outfile, "      Computed %lu integrals.\n\n", writer.count());
    timer_done();
}

void MintsHelper::one_electron_integrals()
{
    shared_ptr<PSIO> psio(new PSIO());

    fprintf(outfile, " OEINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // We'll only be working with the active molecule.
    shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Make sure molecule is valid.
    molecule->update_geometry();

    // Print the molecule.
    molecule->print();

    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basisset = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));

    // Print the basis set
    basisset->print_detail();

    // Create integral factory
    shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset, basisset, basisset, basisset));

    // Get the SO basis object.
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basisset, integral));

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:                %4d\n", molecule->natom());
    fprintf(outfile, "      Number of AO shells:            %4d\n", basisset->nshell());
    fprintf(outfile, "      Number of SO shells:            %4d\n", sobasis->nshell());
    fprintf(outfile, "      Number of primitives:           %4d\n", basisset->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals:      %4d\n", basisset->nao());
    fprintf(outfile, "      Number of basis functions:      %4d\n\n", basisset->nbf());
    fprintf(outfile, "      Number of irreps:               %4d\n", sobasis->nirrep());
    fprintf(outfile, "      Number of functions per irrep: [");
    for (int i=0; i<sobasis->nirrep(); ++i) {
        fprintf(outfile, "%4d ", sobasis->nfunction_in_irrep(i));
    }
    fprintf(outfile, "]\n\n");

    // Compute and dump one-electron SO integrals.
    // Obtain dimensions from the sobasis
    const Dimension dimension = sobasis->dimension();

    // Create a matrix factory and initialize it
    shared_ptr<MatrixFactory> factory(new MatrixFactory());
    factory->init_with(dimension, dimension);

    // Overlap
    shared_ptr<OneBodySOInt> overlap(integral->so_overlap());
    shared_ptr<Matrix>       overlap_mat(factory->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio, PSIF_OEI);

    // Kinetic
    shared_ptr<OneBodySOInt> kinetic(integral->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(factory->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(psio, PSIF_OEI);

    // Potential
    shared_ptr<OneBodySOInt> potential(integral->so_potential());
    shared_ptr<Matrix>       potential_mat(factory->create_matrix(PSIF_SO_V));
    potential->compute(potential_mat);
    potential_mat->save(psio, PSIF_OEI);

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
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));

    int nbf = basis->nbf();
    shared_ptr<Matrix> S(new Matrix("AO Overlap Integrals", nbf, nbf));

    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<OneBodyAOInt> ints(fact->ao_overlap());
    ints->compute(S);

    return S;

}

shared_ptr<Matrix> MintsHelper::ao_kinetic()
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));

    int nbf = basis->nbf();
    shared_ptr<Matrix> S(new Matrix("AO Kinetic Integrals", nbf, nbf));

    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<OneBodyAOInt> ints(fact->ao_kinetic());
    ints->compute(S);

    return S;

}

shared_ptr<Matrix> MintsHelper::ao_potential()
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));

    int nbf = basis->nbf();
    shared_ptr<Matrix> S(new Matrix("AO Potential Integrals", nbf, nbf));

    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<OneBodyAOInt> ints(fact->ao_potential());
    ints->compute(S);

    return S;

}

shared_ptr<Matrix> MintsHelper::ao_erf_eri(double omega, double alpha, double beta)
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));

    int nbf = basis->nbf();
    shared_ptr<Matrix> I(new Matrix("AO ERF ERI Integrals", nbf*nbf, nbf*nbf));

    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<TwoBodyAOInt> ints(fact->erf_eri(omega, alpha, beta));
    double** Ip = I->pointer();
    const double* buffer = ints->buffer();

    for (int M = 0; M < basis->nshell(); M++) {
    for (int N = 0; N < basis->nshell(); N++) {
    for (int P = 0; P < basis->nshell(); P++) {
    for (int Q = 0; Q < basis->nshell(); Q++) {

    ints->compute_shell(M,N,P,Q);

    for (int m = 0, index = 0; m < basis->shell(M)->nfunction(); m++) {
    for (int n = 0; n < basis->shell(N)->nfunction(); n++) {
    for (int p = 0; p < basis->shell(P)->nfunction(); p++) {
    for (int q = 0; q < basis->shell(Q)->nfunction(); q++, index++) {

    Ip[(basis->shell(M)->function_index() + m)*nbf + basis->shell(N)->function_index() + n]
      [(basis->shell(P)->function_index() + p)*nbf + basis->shell(Q)->function_index() + q]
        = buffer[index];

    } } } }

    } } } }


    return I;

}

shared_ptr<Matrix> MintsHelper::ao_eri()
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));

    int nbf = basis->nbf();
    shared_ptr<Matrix> I(new Matrix("AO ERI Integrals", nbf*nbf, nbf*nbf));

    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<TwoBodyAOInt> ints(fact->eri());
    double** Ip = I->pointer();
    const double* buffer = ints->buffer();

    for (int M = 0; M < basis->nshell(); M++) {
    for (int N = 0; N < basis->nshell(); N++) {
    for (int P = 0; P < basis->nshell(); P++) {
    for (int Q = 0; Q < basis->nshell(); Q++) {

    ints->compute_shell(M,N,P,Q);

    for (int m = 0, index = 0; m < basis->shell(M)->nfunction(); m++) {
    for (int n = 0; n < basis->shell(N)->nfunction(); n++) {
    for (int p = 0; p < basis->shell(P)->nfunction(); p++) {
    for (int q = 0; q < basis->shell(Q)->nfunction(); q++, index++) {

    Ip[(basis->shell(M)->function_index() + m)*nbf + basis->shell(N)->function_index() + n]
      [(basis->shell(P)->function_index() + p)*nbf + basis->shell(Q)->function_index() + q]
        = buffer[index];

    } } } }

    } } } }


    return I;

}

shared_ptr<Matrix> MintsHelper::so_overlap()
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));
    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basis, fact));
    const Dimension& dim = sobasis->dimension();
    shared_ptr<Matrix> S(new Matrix("SO Overlap Integrals", dim, dim));

    shared_ptr<OneBodySOInt> ints(fact->so_overlap());
    ints->compute(S);

    return S;

}

shared_ptr<Matrix> MintsHelper::so_kinetic()
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));
    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basis, fact));
    const Dimension& dim = sobasis->dimension();
    shared_ptr<Matrix> S(new Matrix("SO Kinetic Integrals", dim, dim));

    shared_ptr<OneBodySOInt> ints(fact->so_kinetic());
    ints->compute(S);

    return S;

}

shared_ptr<Matrix> MintsHelper::so_potential()
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));
    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basis, fact));
    const Dimension& dim = sobasis->dimension();
    shared_ptr<Matrix> S(new Matrix("SO Potential Integrals", dim, dim));

    shared_ptr<OneBodySOInt> ints(fact->so_potential());
    ints->compute(S);

    return S;
}

std::vector<shared_ptr<Matrix> > MintsHelper::so_dipole()
{
    shared_ptr<Molecule> molecule = Process::environment.molecule();
    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    molecule->update_geometry();
    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options_.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basis = BasisSet::construct(parser, molecule, options_.get_str("BASIS"));
    shared_ptr<IntegralFactory> fact(new IntegralFactory(basis, basis, basis, basis));
    shared_ptr<SOBasisSet> sobasis(new SOBasisSet(basis, fact));
    const Dimension& dim = sobasis->dimension();

    // The matrix factory can create matrices of the correct dimensions...
    shared_ptr<MatrixFactory> factory(new MatrixFactory);
    factory->init_with(dim, dim);

    MultipoleSymmetry msymm(1, molecule, fact, factory);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> dipole = msymm.create_matrices("SO Dipole");

    shared_ptr<OneBodySOInt> ints(fact->so_dipole());
    ints->compute(dipole);

    return dipole;
}

}
