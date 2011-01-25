#include <stdio.h>
#include <stdlib.h>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.h>
#include <libqt/qt.h>

#include <libmints/mints.h>
#include <boost/shared_ptr.hpp>
#include <libparallel/parallel.h>

#include "psi4-dec.h"

#include "deriv.h"

using namespace boost;

namespace psi { namespace deriv {

PsiReturnType deriv(Options & options)
{
    tstart();

    shared_ptr<PSIO> psio(new PSIO);
//    psiopp_ipv1_config(psio);
    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    fprintf(outfile, " DERIV: Wrapper to libmints.\n   by Justin Turney\n\n");

    // We'll only be working with the active molecule.
    shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!\n   Mints wrapper is not meant to be run with IPV1 inputs.");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Create a new matrix factory
    shared_ptr<MatrixFactory> factory(new MatrixFactory);

    // Initialize the factory with data from checkpoint
    factory->init_with_chkpt(chkpt);

    // Read in the basis set
    shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser(options.get_str("BASIS_PATH")));
    shared_ptr<BasisSet> basisset = BasisSet::construct(parser, molecule, options.get_str("BASIS"));

    // Print the molecule.
    basisset->molecule()->print();

    // Print out some useful information
    fprintf(outfile, "   Calculation information:\n");
    fprintf(outfile, "      Number of atoms:           %4d\n", molecule->natom());
    fprintf(outfile, "      Number of shells:          %4d\n", basisset->nshell());
    fprintf(outfile, "      Number of primitives:      %4d\n", basisset->nprimitive());
    fprintf(outfile, "      Number of atomic orbitals: %4d\n", basisset->nao());
    fprintf(outfile, "      Number of basis functions: %4d\n\n", basisset->nbf());

    // Form Q for RHF
    SharedSimpleMatrix Q(factory->create_simple_matrix("Q"));
    int *clsdpi = chkpt->rd_clsdpi();

    // Read in C coefficients
    SharedSimpleMatrix C(factory->create_simple_matrix("MO coefficients"));
    double **vectors = chkpt->rd_scf();
    if (vectors == NULL) {
        fprintf(stderr, "Could not find MO coefficients. Run scf first.\n");
        return Failure;
    }
    C->set(vectors);
    free_block(vectors);

    // Load in orbital energies
    double *etmp = chkpt->rd_evals();
    shared_ptr<SimpleMatrix> W(factory->create_simple_matrix("W"));

    int nso = chkpt->rd_nso();
    for (int m=0; m<nso; ++m) {
        for (int n=0; n<nso; ++n) {
            double sum=0.0;
            double qsum =0.0;
            for (int i=0; i<clsdpi[0]; ++i) {
                sum += C->get(m, i) * C->get(n, i) * etmp[i];
                qsum += C->get(m, i) * C->get(n, i);
            }
            W->set(m, n, sum);
            Q->set(m, n, qsum);
        }
    }
    Chkpt::free(etmp);

//    Q->print();
//    fprintf(outfile, "AO-basis\n");
//    W->print();

    SharedSimpleMatrix G;
    Deriv deriv(ref_rhf, factory, basisset);
    SharedSimpleMatrix WdS(deriv.overlap());
    SharedSimpleMatrix QdH(deriv.one_electron());
    SharedSimpleMatrix tb(deriv.two_body());
    deriv.compute(C, Q, G, W);

    SimpleMatrix enuc = basisset->molecule()->nuclear_repulsion_energy_deriv1();

    enuc.print_atom_vector();
    QdH->print_atom_vector();
    WdS->print_atom_vector();
    tb->print_atom_vector();

    SimpleMatrix scf_grad("SCF gradient", basisset->molecule()->natom(), 3);
    scf_grad.add(&enuc);
    scf_grad.add(QdH);
    scf_grad.add(WdS);
    scf_grad.add(tb);

    scf_grad.print_atom_vector();

    GradientWriter grad(basisset->molecule(), scf_grad);
    grad.write("psi.file11.dat");

//    SimpleMatrix enuc2 = basis->molecule()->nuclear_repulsion_energy_deriv2();
//    enuc2.print();

    // Shut down psi
    tstop();

    Chkpt::free(clsdpi);

    return Success;
}

}}
