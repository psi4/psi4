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
#include <psiconfig.h>

#include "matrix_distributed.h"

using namespace boost;

namespace psi {

void MintsHelper::init_helper(boost::shared_ptr<Wavefunction> wavefunction)
{
//    if (Communicator::world->nproc() > 1)
//        throw PSIEXCEPTION("MintsHelper: You can only use one process with this.");

    if (wavefunction) {
        psio_ = wavefunction->psio();
        molecule_ = wavefunction->molecule();
    }
    else {
        psio_ = boost::shared_ptr<PSIO>(new PSIO());
        molecule_ = boost::shared_ptr<Molecule>(Process::environment.molecule());
    }

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
    if(wavefunction)
        basisset_ = wavefunction->basisset();
    else {
        boost::shared_ptr<BasisSetParser> parser (new Gaussian94BasisSetParser());
        basisset_ = boost::shared_ptr<BasisSet>(BasisSet::construct(parser, molecule_, "BASIS"));
    }

    // Print the basis set
    if (print_)
        basisset_->print_detail();

    // Create integral factory
    integral_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));

    // Get the SO basis object.
    sobasis_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(basisset_, integral_));

    // Obtain dimensions from the sobasis
    const Dimension dimension = sobasis_->dimension();

    // Create a matrix factory and initialize it
    factory_ = boost::shared_ptr<MatrixFactory>(new MatrixFactory());
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

MintsHelper::MintsHelper(boost::shared_ptr<Wavefunction> wavefunction) :options_(wavefunction->options())
{
    init_helper(wavefunction);
}

MintsHelper::~MintsHelper()
{
}

boost::shared_ptr<BasisSet> MintsHelper::basisset() const
{
    return basisset_;
}

boost::shared_ptr<SOBasisSet> MintsHelper::sobasisset() const
{
    return sobasis_;
}

boost::shared_ptr<MatrixFactory> MintsHelper::factory() const
{
    return factory_;
}

void MintsHelper::integrals()
{
    fprintf(outfile, " MINTS: Wrapper to libmints.\n   by Justin Turney\n\n");

    // Get ERI object
    std::vector<boost::shared_ptr<TwoBodyAOInt> > tb;
    for (int i=0; i<Communicator::world->nthread(); ++i)
        tb.push_back(boost::shared_ptr<TwoBodyAOInt>(integral_->eri()));
    boost::shared_ptr<TwoBodySOInt> eri(new TwoBodySOInt(tb, integral_));

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
    boost::shared_ptr<OneBodySOInt> overlap(integral_->so_overlap());
    boost::shared_ptr<Matrix>       overlap_mat(factory_->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);

    // Kinetic
    boost::shared_ptr<OneBodySOInt> kinetic(integral_->so_kinetic());
    boost::shared_ptr<Matrix>       kinetic_mat(factory_->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(psio_, PSIF_OEI);

    // Potential
    boost::shared_ptr<OneBodySOInt> potential(integral_->so_potential());
    boost::shared_ptr<Matrix>       potential_mat(factory_->create_matrix(PSIF_SO_V));
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
    boost::shared_ptr<OneBodySOInt> overlap(integral_->so_overlap());
    boost::shared_ptr<Matrix>       overlap_mat(factory_->create_matrix(PSIF_SO_S));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);

    // Kinetic
    boost::shared_ptr<OneBodySOInt> kinetic(integral_->so_kinetic());
    boost::shared_ptr<Matrix>       kinetic_mat(factory_->create_matrix(PSIF_SO_T));
    kinetic->compute(kinetic_mat);
    kinetic_mat->save(psio_, PSIF_OEI);

    // Potential
    boost::shared_ptr<OneBodySOInt> potential(integral_->so_potential());
    boost::shared_ptr<Matrix>       potential_mat(factory_->create_matrix(PSIF_SO_V));
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

boost::shared_ptr<Matrix> MintsHelper::ao_overlap()
{
    // Overlap
    boost::shared_ptr<OneBodyAOInt> overlap(integral_->ao_overlap());
    boost::shared_ptr<Matrix>       overlap_mat(new Matrix(PSIF_AO_S, basisset_->nbf (), basisset_->nbf ()));
    overlap->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);
    return overlap_mat;
}

boost::shared_ptr<Matrix> MintsHelper::ao_kinetic()
{
    boost::shared_ptr<OneBodyAOInt> T(integral_->ao_kinetic());
    boost::shared_ptr<Matrix>       kinetic_mat(new Matrix( basisset_->nbf (), basisset_->nbf ()));
    T->compute(kinetic_mat);
    return kinetic_mat;
}

boost::shared_ptr<Matrix> MintsHelper::ao_potential()
{
    boost::shared_ptr<OneBodyAOInt> V(integral_->ao_potential());
    boost::shared_ptr<Matrix>       potential_mat(new Matrix(basisset_->nbf (), basisset_->nbf ()));
    V->compute(potential_mat);
    return potential_mat;
}

boost::shared_ptr<Matrix> MintsHelper::ao_erf_eri(double omega)
{
    int nbf = basisset_->nbf();
    boost::shared_ptr<Matrix> I(new Matrix("AO ERF ERI Integrals", nbf*nbf, nbf*nbf));

    boost::shared_ptr<TwoBodyAOInt> ints(integral_->erf_eri(omega));
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

boost::shared_ptr<Matrix> MintsHelper::ao_eri()
{
    int nbf = basisset_->nbf();
    boost::shared_ptr<Matrix> I(new Matrix("AO ERI Tensor", nbf*nbf, nbf*nbf));

    boost::shared_ptr<TwoBodyAOInt> ints(integral_->eri());
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
boost::shared_ptr<Matrix> MintsHelper::mo_erf_eri(double omega, boost::shared_ptr<Matrix> C1, boost::shared_ptr<Matrix> C2,
                                                                boost::shared_ptr<Matrix> C3, boost::shared_ptr<Matrix> C4)
{
    boost::shared_ptr<Matrix> mo_ints = mo_eri_helper(ao_erf_eri(omega), C1, C2, C3, C4);
    mo_ints->set_name("MO ERF ERI Tensor");
    return mo_ints;
}
boost::shared_ptr<Matrix> MintsHelper::mo_eri(boost::shared_ptr<Matrix> C1, boost::shared_ptr<Matrix> C2,
                                              boost::shared_ptr<Matrix> C3, boost::shared_ptr<Matrix> C4)
{
    boost::shared_ptr<Matrix> mo_ints = mo_eri_helper(ao_eri(), C1, C2, C3, C4);
    mo_ints->set_name("MO ERI Tensor");
    return mo_ints;
}
boost::shared_ptr<Matrix> MintsHelper::mo_erf_eri(double omega, boost::shared_ptr<Matrix> Co, boost::shared_ptr<Matrix> Cv)
{
    boost::shared_ptr<Matrix> mo_ints = mo_eri_helper(ao_erf_eri(omega), Co, Cv);
    mo_ints->set_name("MO ERF ERI Tensor");
    return mo_ints;
}
boost::shared_ptr<Matrix> MintsHelper::mo_eri(boost::shared_ptr<Matrix> Co, boost::shared_ptr<Matrix> Cv)
{
    boost::shared_ptr<Matrix> mo_ints = mo_eri_helper(ao_eri(), Co, Cv);
    mo_ints->set_name("MO ERI Tensor");
    return mo_ints;
}
boost::shared_ptr<Matrix> MintsHelper::mo_eri_helper(boost::shared_ptr<Matrix> Iso, boost::shared_ptr<Matrix> C1, boost::shared_ptr<Matrix> C2,
                                                                                    boost::shared_ptr<Matrix> C3, boost::shared_ptr<Matrix> C4)
{
    int nso = basisset_->nbf();
    int n1 = C1->colspi()[0];
    int n2 = C2->colspi()[0];
    int n3 = C3->colspi()[0];
    int n4 = C4->colspi()[0];

    double** C1p = C1->pointer();
    double** C2p = C2->pointer();
    double** C3p = C3->pointer();
    double** C4p = C4->pointer();

    double** Isop = Iso->pointer();
    boost::shared_ptr<Matrix> I2(new Matrix("MO ERI Tensor", n1 * nso, nso * nso));
    double** I2p = I2->pointer();

    C_DGEMM('T','N',n1,nso * (ULI) nso * nso,nso,1.0,C1p[0],n1,Isop[0],nso * (ULI) nso * nso,0.0,I2p[0],nso * (ULI) nso * nso);

    Iso.reset();
    boost::shared_ptr<Matrix> I3(new Matrix("MO ERI Tensor", n1 * nso, nso * n3));
    double** I3p = I3->pointer();

    C_DGEMM('N','N',n1 * (ULI) nso * nso,n3,nso,1.0,I2p[0],nso,C3p[0],n3,0.0,I3p[0],n3);

    I2.reset();
    boost::shared_ptr<Matrix> I4(new Matrix("MO ERI Tensor", nso * n1, n3 * nso));
    double** I4p = I4->pointer();

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
    boost::shared_ptr<Matrix> I5(new Matrix("MO ERI Tensor", n2 * n1, n3 * nso));
    double** I5p = I5->pointer();

    C_DGEMM('T','N',n2,n1 * (ULI) n3 * nso, nso,1.0,C2p[0],n2,I4p[0],n1*(ULI)n3*nso,0.0,I5p[0],n1*(ULI)n3*nso);

    I4.reset();
    boost::shared_ptr<Matrix> I6(new Matrix("MO ERI Tensor", n2 * n1, n3 * n4));
    double** I6p = I6->pointer();

    C_DGEMM('N','N',n2 * (ULI) n1 * n3, n4, nso,1.0,I5p[0],nso,C4p[0],n4,0.0,I6p[0],n4);

    I5.reset();
    boost::shared_ptr<Matrix> Imo(new Matrix("MO ERI Tensor", n1 * n2, n3 * n4));
    double** Imop = Imo->pointer();

    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n3; j++) {
            for (int a = 0; a < n2; a++) {
                for (int b = 0; b < n4; b++) {
                    Imop[i * n2 + a][j * n4 + b] = I6p[a * n1 + i][j * n4 + b];
                }
            }
        }
    }

    return Imo;
}
boost::shared_ptr<Matrix> MintsHelper::mo_eri_helper(boost::shared_ptr<Matrix> Iso, boost::shared_ptr<Matrix> Co, boost::shared_ptr<Matrix> Cv)
{
    int nso = basisset_->nbf();
    int nocc = Co->colspi()[0];
    int nvir = Cv->colspi()[0];

    double** Cop = Co->pointer();
    double** Cvp = Cv->pointer();

    double** Isop = Iso->pointer();
    boost::shared_ptr<Matrix> I2(new Matrix("MO ERI Tensor", nocc * nso, nso * nso));
    double** I2p = I2->pointer();

    C_DGEMM('T','N',nocc,nso * (ULI) nso * nso,nso,1.0,Cop[0],nocc,Isop[0],nso * (ULI) nso * nso,0.0,I2p[0],nso * (ULI) nso * nso);

    Iso.reset();
    boost::shared_ptr<Matrix> I3(new Matrix("MO ERI Tensor", nocc * nso, nso * nocc));
    double** I3p = I3->pointer();

    C_DGEMM('N','N',nocc * (ULI) nso * nso,nocc,nso,1.0,I2p[0],nso,Cop[0],nocc,0.0,I3p[0],nocc);

    I2.reset();
    boost::shared_ptr<Matrix> I4(new Matrix("MO ERI Tensor", nso * nocc, nocc * nso));
    double** I4p = I4->pointer();

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
    boost::shared_ptr<Matrix> I5(new Matrix("MO ERI Tensor", nvir * nocc, nocc * nso));
    double** I5p = I5->pointer();

    C_DGEMM('T','N',nvir,nocc * (ULI) nocc * nso, nso,1.0,Cvp[0],nvir,I4p[0],nocc*(ULI)nocc*nso,0.0,I5p[0],nocc*(ULI)nocc*nso);

    I4.reset();
    boost::shared_ptr<Matrix> I6(new Matrix("MO ERI Tensor", nvir * nocc, nocc * nvir));
    double** I6p = I6->pointer();

    C_DGEMM('N','N',nvir * (ULI) nocc * nocc, nvir, nso,1.0,I5p[0],nso,Cvp[0],nvir,0.0,I6p[0],nvir);

    I5.reset();
    boost::shared_ptr<Matrix> Imo(new Matrix("MO ERI Tensor", nocc * nvir, nocc * nvir));
    double** Imop = Imo->pointer();

    for (int i = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++) {
            for (int a = 0; a < nvir; a++) {
                for (int b = 0; b < nvir; b++) {
                    Imop[i * nvir + a][j * nvir + b] = I6p[a * nocc + i][j * nvir + b];
                }
            }
        }
    }

    return Imo;
}
boost::shared_ptr<Matrix> MintsHelper::so_overlap()
{
    boost::shared_ptr<OneBodySOInt> S(integral_->so_overlap());
    boost::shared_ptr<Matrix>       overlap_mat(factory_->create_matrix(PSIF_SO_S));
    S->compute(overlap_mat);
    overlap_mat->save(psio_, PSIF_OEI);
    return overlap_mat;
}

boost::shared_ptr<Matrix> MintsHelper::so_kinetic()
{
    boost::shared_ptr<OneBodySOInt> T(integral_->so_kinetic());
    boost::shared_ptr<Matrix>       kinetic_mat(factory_->create_matrix(PSIF_SO_T));
    T->compute(kinetic_mat);
    kinetic_mat->save(psio_, PSIF_OEI);
    return kinetic_mat;
}

boost::shared_ptr<Matrix> MintsHelper::so_potential()
{
    boost::shared_ptr<OneBodySOInt> V(integral_->so_potential());
    boost::shared_ptr<Matrix>       potential_mat(factory_->create_matrix(PSIF_SO_V));
    V->compute(potential_mat);
    potential_mat->save(psio_, PSIF_OEI);
    return potential_mat;
}

std::vector<boost::shared_ptr<Matrix> > MintsHelper::so_dipole()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(1, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> dipole = msymm.create_matrices("SO Dipole");

    boost::shared_ptr<OneBodySOInt> ints(integral_->so_dipole());
    ints->compute(dipole);

    return dipole;
}
std::vector<boost::shared_ptr<Matrix> > MintsHelper::so_quadrupole()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole = msymm.create_matrices("SO Quadrupole");

    boost::shared_ptr<OneBodySOInt> ints(integral_->so_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}
std::vector<boost::shared_ptr<Matrix> > MintsHelper::so_traceless_quadrupole()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(2, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> quadrupole = msymm.create_matrices("SO Traceless Quadrupole");

    boost::shared_ptr<OneBodySOInt> ints(integral_->so_traceless_quadrupole());
    ints->compute(quadrupole);

    return quadrupole;
}

std::vector<boost::shared_ptr<Matrix> > MintsHelper::so_nabla()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(OperatorSymmetry::P, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> nabla = msymm.create_matrices("SO Nabla");

    boost::shared_ptr<OneBodySOInt> ints(integral_->so_nabla());
    ints->compute(nabla);

    return nabla;
}

std::vector<boost::shared_ptr<Matrix> > MintsHelper::so_angular_momentum()
{
    // The matrix factory can create matrices of the correct dimensions...
    OperatorSymmetry msymm(OperatorSymmetry::L, molecule_, integral_, factory_);
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> am = msymm.create_matrices("SO Angular Momentum");

    boost::shared_ptr<OneBodySOInt> ints(integral_->so_angular_momentum());
    ints->compute(am);

    return am;
}

std::vector<boost::shared_ptr<Matrix> > MintsHelper::ao_angular_momentum()
{
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> angmom;

    angmom.push_back(SharedMatrix(new Matrix("AO Lx", basisset_->nbf(), basisset_->nbf())));
    angmom.push_back(SharedMatrix(new Matrix("AO Ly", basisset_->nbf(), basisset_->nbf())));
    angmom.push_back(SharedMatrix(new Matrix("AO Lz", basisset_->nbf(), basisset_->nbf())));

    boost::shared_ptr<OneBodyAOInt> ints(integral_->ao_angular_momentum());
    ints->compute(angmom);

    return angmom;
}

std::vector<boost::shared_ptr<Matrix> > MintsHelper::ao_nabla()
{
    // Create a vector of matrices with the proper symmetry
    std::vector<SharedMatrix> nabla;

    nabla.push_back(SharedMatrix(new Matrix("AO Px", basisset_->nbf(), basisset_->nbf())));
    nabla.push_back(SharedMatrix(new Matrix("AO Py", basisset_->nbf(), basisset_->nbf())));
    nabla.push_back(SharedMatrix(new Matrix("AO Pz", basisset_->nbf(), basisset_->nbf())));

    boost::shared_ptr<OneBodyAOInt> ints(integral_->ao_nabla());
    ints->compute(nabla);

    return nabla;
}

void MintsHelper::play()
{
#ifdef HAVE_MADNESS
    int M = 10;
    int N = 10;
    int tile_size = 3; // This makes the tiles 3 x 3
    Distributed_Matrix t(M, N, tile_size, "T");
    t.print_all_tiles();
#endif
}

} // namespace psi
