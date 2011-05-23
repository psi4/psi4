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

    boost::shared_ptr<PSIO> psio(new PSIO);

    fprintf(outfile, " DERIV: Wrapper to libmints.\n   by Justin Turney\n\n");

    // We'll only be working with the active molecule.
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (!molecule) {
        throw PSIEXCEPTION("deriv: Active molecule not set!");
    }

    // Create a new matrix factory
    boost::shared_ptr<MatrixFactory> factory(new MatrixFactory);

    // Read in the basis set
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> basisset = BasisSet::construct(parser, molecule, "BASIS");

    // Initialize an integral object.
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset, basisset, basisset, basisset));

    // Create an SOBasisSet
    boost::shared_ptr<SOBasisSet> sobasisset(new SOBasisSet(basisset, integral));
    // TODO change this back to use sobasisset's petitelist when we get the SO derivs working.  We make a new one
    // here because we don't want the pure transform folded in when we use AO integrals.
    shared_ptr<PetiteList> petite(new PetiteList(basisset, integral, false));
    SharedMatrix usotoao(petite->sotoao());
//    SharedMatrix usotoao(sobasisset->petitelist()->sotoao());

    // From the SOBasisSet obtain the block dimensions
    const Dimension dimension = sobasisset->dimension();

    // Initialize the factory
    factory->init_with(dimension, dimension);

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
    int *clsdpi = Process::environment.reference_wavefunction()->doccpi();

    // Read in C coefficients
    SharedMatrix Cso = Process::environment.reference_wavefunction()->Ca();
    SharedSimpleMatrix simple_Cso(new SimpleMatrix("Cso",
                                                   Process::environment.reference_wavefunction()->nso(),
                                                   Process::environment.reference_wavefunction()->nmo()));

    SharedSimpleMatrix simple_usotoao(new SimpleMatrix("USO -> AO",
                                                       Process::environment.reference_wavefunction()->nso(),
                                                       basisset->nbf()));

    int sooffset = 0, mooffset = 0;
    for (int h=0; h<Cso->nirrep(); ++h) {
        for (int m=0; m<Cso->rowspi()[h]; ++m) {
            for (int n=0; n<basisset->nbf(); ++n) {
                simple_usotoao->set(sooffset+m, n, usotoao->get(h, m, n));
            }
            for (int n=0; n<Process::environment.reference_wavefunction()->nmopi()[h]; ++n) {
                simple_Cso->set(sooffset+m, n+mooffset, Cso->get(h, m, n));
            }
        }
        sooffset += Process::environment.reference_wavefunction()->nsopi()[h];
        mooffset += Cso->rowspi()[h];
    }

    int nso = Process::environment.reference_wavefunction()->nso();

    // Load in orbital energies
    SharedVector etmp = Process::environment.reference_wavefunction()->epsilon_a();
    boost::shared_ptr<SimpleMatrix> W(factory->create_simple_matrix("W"));

    int nbf = basisset->nbf();
    for (int m=0; m<nbf; ++m) {
        for (int n=0; n<nbf; ++n) {
            double sum=0.0;
            double qsum=0.0;

            int mooffset = 0;
            for (int h=0; h<sobasisset->nirrep(); ++h) {
                for (int i=0; i<clsdpi[h]; ++i) {
                    sum  += simple_Cso->get(m, i+mooffset) * simple_Cso->get(n, i+mooffset) * etmp->get(h, i);
                    qsum += simple_Cso->get(m, i+mooffset) * simple_Cso->get(n, i+mooffset);
                }
                mooffset += Process::environment.reference_wavefunction()->nmopi()[h];
            }
            W->set(m, n, sum);
            Q->set(m, n, qsum);
        }
    }

    Q->transform(simple_usotoao);
    W->transform(simple_usotoao);

    SharedSimpleMatrix G;

    Deriv deriv(ref_rhf, factory, basisset);
    SharedSimpleMatrix WdS(deriv.overlap());
    SharedSimpleMatrix QdH(deriv.one_electron());
    SharedSimpleMatrix tb(deriv.two_body());
    deriv.compute(Q, G, W);

    SimpleMatrix enuc = basisset->molecule()->nuclear_repulsion_energy_deriv1();

    enuc.print_atom_vector();
    QdH->print_atom_vector();
    WdS->print_atom_vector();
    tb->print_atom_vector();

    SimpleMatrix scf_grad("SCF gradient", basisset->molecule()->natom(), 3);
    SimpleMatrix temp("Temp SCF gradient", basisset->molecule()->natom(), 3);
    temp.add(&enuc);
    temp.add(QdH);
    temp.add(WdS);
    temp.add(tb);

    CharacterTable ct = basisset->molecule()->point_group()->char_table();

    // Obtain atom mapping of atom * symm op to atom
    int **atom_map = compute_atom_map(basisset->molecule());

    // Symmetrize the gradients to remove any noise
    for (int atom=0; atom<basisset->molecule()->natom(); ++atom) {
        for (int g=0; g<ct.order(); ++g) {

            int Gatom = atom_map[atom][g];

            SymmetryOperation so = ct.symm_operation(g);

            scf_grad.add(atom, 0, so(0, 0) * temp.get(Gatom, 0) / ct.order());
            scf_grad.add(atom, 1, so(1, 1) * temp.get(Gatom, 1) / ct.order());
            scf_grad.add(atom, 2, so(2, 2) * temp.get(Gatom, 2) / ct.order());
        }
    }

    // Delete the atom map.
    delete_atom_map(atom_map, basisset->molecule());

    // Print the atom vector
    scf_grad.print_atom_vector();

    GradientWriter grad(basisset->molecule(), scf_grad);
    grad.write("psi.file11.dat");

//    SimpleMatrix enuc2 = basis->molecule()->nuclear_repulsion_energy_deriv2();
//    enuc2.print();

    // Shut down psi
    tstop();

    return Success;
}

}}
