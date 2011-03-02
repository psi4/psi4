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
        _psio = shared_ptr<PSIO>(new PSIO());
        _molecule = shared_ptr<Molecule>(Process::environment.molecule());

        if (_molecule.get() == 0) {
            fprintf(outfile, "  Active molecule not set!");
            throw PSIEXCEPTION("Active molecule not set!");
        }

        // Make sure molecule is valid.
        _molecule->update_geometry();

        // Print the molecule.
        _molecule->print();

        // Read in the basis set
        shared_ptr<BasisSetParser> parser (new Gaussian94BasisSetParser());
        _basisset = shared_ptr<BasisSet>(BasisSet::construct(parser, _molecule, "BASIS"));

        // Print the basis set
        _basisset->print_detail();

        // Create integral factory
        _integral = shared_ptr<IntegralFactory>(new IntegralFactory(_basisset, _basisset, _basisset, _basisset));

        // Get the SO basis object.
        _sobasis = shared_ptr<SOBasisSet>(new SOBasisSet(_basisset, _integral));

        // Obtain dimensions from the sobasis
        const Dimension dimension = _sobasis->dimension();

        // Create a matrix factory and initialize it
        _factory = shared_ptr<MatrixFactory>(new MatrixFactory());
        _factory->init_with(dimension, dimension);

    }


MintsHelper::MintsHelper(Options & options): options_(options)
{
    MintsHelper::init_helper();
}

MintsHelper::MintsHelper() : options_(Process::environment.options)
{
    MintsHelper::init_helper();
}


MintsHelper::~MintsHelper()
{
}

void MintsHelper::integrals()
{
    timer_init();

    fprintf(outfile, " MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // Get ERI object
    shared_ptr<TwoBodyAOInt> tb(_integral->eri());
    shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, _integral));

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:                %4d\n", _molecule->natom());
    fprintf(outfile, "      Number of AO shells:            %4d\n", _basisset->nshell());
    fprintf(outfile, "      Number of SO shells:            %4d\n", _sobasis->nshell());
    fprintf(outfile, "      Number of primitives:           %4d\n", _basisset->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals:      %4d\n", _basisset->nao());
    fprintf(outfile, "      Number of basis functions:      %4d\n\n", _basisset->nbf());
    fprintf(outfile, "      Number of irreps:               %4d\n", _sobasis->nirrep());
    fprintf(outfile, "      Number of functions per irrep: [");
    for (int i=0; i<_sobasis->nirrep(); ++i) {
        fprintf(outfile, "%4d ", _sobasis->nfunction_in_irrep(i));
    }
    fprintf(outfile, "]\n\n");

    // Compute and dump one-electron SO integrals.

    // Overlap
    shared_ptr<OneBodySOInt> overlap(_integral->so_overlap());
    shared_ptr<Matrix>       overlap_mat(_factory->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(_psio, PSIF_OEI);

    // Kinetic
    shared_ptr<OneBodySOInt> kinetic(_integral->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(_factory->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(_psio, PSIF_OEI);

    // Potential
    shared_ptr<OneBodySOInt> potential(_integral->so_potential());
    shared_ptr<Matrix>       potential_mat(_factory->create_matrix(PSIF_SO_V));
    potential->compute(potential_mat);
    potential_mat->save(_psio, PSIF_OEI);

    if (Process::environment.options.get_int("PRINT") > 3) {
        overlap_mat->print();
        kinetic_mat->print();
        potential_mat->print();
    }

    // Open the IWL buffer where we will store the integrals.
    IWL ERIOUT(_psio.get(), PSIF_SO_TEI, 0.0, 0, 0);
    IWLWriter writer(ERIOUT);

    // Let the user know what we're doing.
    fprintf(outfile, "      Computing integrals..."); fflush(outfile);

    SOShellCombinationsIterator shellIter(_sobasis, _sobasis, _sobasis, _sobasis);
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

    fprintf(outfile, " OEINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:                %4d\n", _molecule->natom());
    fprintf(outfile, "      Number of AO shells:            %4d\n", _basisset->nshell());
    fprintf(outfile, "      Number of SO shells:            %4d\n", _sobasis->nshell());
    fprintf(outfile, "      Number of primitives:           %4d\n", _basisset->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals:      %4d\n", _basisset->nao());
    fprintf(outfile, "      Number of basis functions:      %4d\n\n", _basisset->nbf());
    fprintf(outfile, "      Number of irreps:               %4d\n", _sobasis->nirrep());
    fprintf(outfile, "      Number of functions per irrep: [");
    for (int i=0; i<_sobasis->nirrep(); ++i) {
        fprintf(outfile, "%4d ", _sobasis->nfunction_in_irrep(i));
    }
    fprintf(outfile, "]\n\n");

    // Compute and dump one-electron SO integrals.

     // Overlap
    shared_ptr<OneBodySOInt> overlap(_integral->so_overlap());
    shared_ptr<Matrix>       overlap_mat(_factory->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(_psio, PSIF_OEI);

    // Kinetic
    shared_ptr<OneBodySOInt> kinetic(_integral->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(_factory->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(_psio, PSIF_OEI);

    // Potential
    shared_ptr<OneBodySOInt> potential(_integral->so_potential());
    shared_ptr<Matrix>       potential_mat(_factory->create_matrix(PSIF_SO_V));
    potential->compute(potential_mat);
    potential_mat->save(_psio, PSIF_OEI);
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
   shared_ptr<OneBodyAOInt> overlap(_integral->ao_overlap());
   shared_ptr<Matrix>       overlap_mat(_factory->create_matrix(PSIF_AO_S));
   overlap->compute(overlap_mat);
   overlap_mat->save(_psio, PSIF_OEI);

}

shared_ptr<Matrix> MintsHelper::ao_kinetic()
{

    shared_ptr<OneBodyAOInt> T(_integral->ao_kinetic());
    shared_ptr<Matrix>       kinetic_mat(_factory->create_matrix());
    T->compute(kinetic_mat);
    kinetic_mat->save(_psio, PSIF_OEI);

}

shared_ptr<Matrix> MintsHelper::ao_potential()
{

    shared_ptr<OneBodyAOInt> V(_integral->ao_potential());
    shared_ptr<Matrix>       potential_mat(_factory->create_matrix());
    V->compute(potential_mat);
    potential_mat->save(_psio, PSIF_OEI);

}

shared_ptr<Matrix> MintsHelper::ao_erf_eri(double omega, double alpha, double beta)
{

    int nbf = _basisset->nbf();
    shared_ptr<Matrix> I(new Matrix("AO ERF ERI Integrals", nbf*nbf, nbf*nbf));

    shared_ptr<TwoBodyAOInt> ints(_integral->erf_eri(omega, alpha, beta));
    double** Ip = I->pointer();
    const double* buffer = ints->buffer();

    for (int M = 0; M < _basisset->nshell(); M++) {
    for (int N = 0; N < _basisset->nshell(); N++) {
    for (int P = 0; P < _basisset->nshell(); P++) {
    for (int Q = 0; Q < _basisset->nshell(); Q++) {

    ints->compute_shell(M,N,P,Q);

    for (int m = 0, index = 0; m < _basisset->shell(M)->nfunction(); m++) {
    for (int n = 0; n < _basisset->shell(N)->nfunction(); n++) {
    for (int p = 0; p < _basisset->shell(P)->nfunction(); p++) {
    for (int q = 0; q < _basisset->shell(Q)->nfunction(); q++, index++) {

    Ip[(_basisset->shell(M)->function_index() + m)*nbf + _basisset->shell(N)->function_index() + n]
      [(_basisset->shell(P)->function_index() + p)*nbf + _basisset->shell(Q)->function_index() + q]
        = buffer[index];

    } } } }

    } } } }


    return I;

}

shared_ptr<Matrix> MintsHelper::ao_eri()
{

    int nbf = _basisset->nbf();
    shared_ptr<Matrix> I(new Matrix("AO ERI Integrals", nbf*nbf, nbf*nbf));

    shared_ptr<TwoBodyAOInt> ints(_integral->eri());
    double** Ip = I->pointer();
    const double* buffer = ints->buffer();

    for (int M = 0; M < _basisset->nshell(); M++) {
    for (int N = 0; N < _basisset->nshell(); N++) {
    for (int P = 0; P < _basisset->nshell(); P++) {
    for (int Q = 0; Q < _basisset->nshell(); Q++) {

    ints->compute_shell(M,N,P,Q);

    for (int m = 0, index = 0; m < _basisset->shell(M)->nfunction(); m++) {
    for (int n = 0; n < _basisset->shell(N)->nfunction(); n++) {
    for (int p = 0; p < _basisset->shell(P)->nfunction(); p++) {
    for (int q = 0; q < _basisset->shell(Q)->nfunction(); q++, index++) {

    Ip[(_basisset->shell(M)->function_index() + m)*nbf + _basisset->shell(N)->function_index() + n]
      [(_basisset->shell(P)->function_index() + p)*nbf + _basisset->shell(Q)->function_index() + q]
        = buffer[index];

    } } } }

    } } } }


    return I;

}

shared_ptr<Matrix> MintsHelper::so_overlap()
{

    shared_ptr<OneBodySOInt> S(_integral->so_overlap());
    shared_ptr<Matrix>       overlap_mat(_factory->create_matrix(PSIF_SO_S));
    S->compute(overlap_mat);
    overlap_mat->save(_psio, PSIF_OEI);

}

shared_ptr<Matrix> MintsHelper::so_kinetic()
{

    shared_ptr<OneBodySOInt> T(_integral->so_kinetic());
    shared_ptr<Matrix>       kinetic_mat(_factory->create_matrix(PSIF_SO_T));
    T->compute(kinetic_mat);
    kinetic_mat->save(_psio, PSIF_OEI);

}

shared_ptr<Matrix> MintsHelper::so_potential()
{

    shared_ptr<OneBodySOInt> V(_integral->so_potential());
    shared_ptr<Matrix>       potential_mat(_factory->create_matrix(PSIF_SO_V));
    V->compute(potential_mat);
    potential_mat->save(_psio, PSIF_OEI);

}

std::vector<shared_ptr<Matrix> > MintsHelper::so_dipole()
{

    // The matrix factory can create matrices of the correct dimensions...

    MultipoleSymmetry msymm(1, _molecule, _integral, _factory);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> dipole = msymm.create_matrices("SO Dipole");

    shared_ptr<OneBodySOInt> ints(_integral->so_dipole());
    ints->compute(dipole);

    return dipole;
}

}
