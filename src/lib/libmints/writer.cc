#include "writer.h"
#include <libmints/mints.h>
#include <psi4-dec.h>

#include <cstdio>

using namespace psi;
using namespace boost;

GradientWriter::GradientWriter(shared_ptr<Molecule> mol, const SimpleMatrix& grad)
    : molecule_(mol), gradient_(grad)
{
}

void GradientWriter::write(const std::string &filename)
{
    FILE *file11 = fopen(filename.c_str(), "a");
    int i;

    if (file11 == NULL)
        throw PSIEXCEPTION("GradientWriter::write: Unable to open file11.dat");

    fprintf(file11, "%-59.59s %-10.10s%-8.8s\n",
            molecule_->name().c_str(),
            "(wfn)",
            "(dertype)");

    fprintf(file11, "%5d%20.10lf\n", molecule_->natom(), Process::environment.globals["CURRENT ENERGY"]);

    for (i=0; i<molecule_->natom(); ++i) {
        fprintf(file11, "%20.10lf%20.10lf%20.10lf%20.10lf\n",
                double(molecule_->Z(i)), molecule_->x(i), molecule_->y(i), molecule_->z(i));
    }

    for (i=0; i<molecule_->natom(); ++i) {
        fprintf(file11, "                    %20.10lf%20.10lf%20.10lf\n",
                gradient_(i, 0), gradient_(i, 1), gradient_(i, 2));
    }

    fclose(file11);
}

MoldenWriter::MoldenWriter(boost::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction)
{

}

void MoldenWriter::write(const std::string &filename)
{
    FILE *molden = fopen(filename.c_str(), "a");
    int atom;

    fprintf(molden, "[Molden Format]\n");

    // Get the molecule for ease
    BasisSet& basisset = *wavefunction_->basisset().get();
    SOBasisSet& sobasisset = *wavefunction_->sobasisset().get();
    Molecule& mol = *basisset.molecule().get();

//    basisset.print_detail();

    // Print the molecule to molden
    fprintf(molden, "[Atoms] (AU)\n");
    for (atom=0; atom<mol.natom(); ++atom) {
        Vector3 coord = mol.xyz(atom);
        fprintf(molden, "%-2s  %2d  %3d   %20.12f %20.12f %20.12f\n",
                mol.symbol(atom).c_str(), atom+1, mol.Z(atom), coord[0], coord[1], coord[2]);
    }

    // Dump the basis set using code adapted from psi2molden
    fprintf(molden, "[GTO]\n");

    // For each atom
    for (atom=0; atom<mol.natom(); ++atom) {
        fprintf(molden, "  %3d 0\n", atom+1);

        // Go through all the shells on this center
        for (int shell=0; shell < basisset.nshell_on_center(atom); ++shell) {
            int overall_shell = basisset.shell_on_center(atom, shell);

            GaussianShell& gs = *basisset.shell(overall_shell).get();

            fprintf(molden, "%c %3d 1.00\n", gs.amchar(), gs.nprimitive());

            for (int prim=0; prim<gs.nprimitive(); ++prim) {
                fprintf(molden, "%20.10f %20.10f\n", gs.exp(prim), gs.coef(prim));
            }
        }

        // An empty line separates atoms
        fprintf(molden, "\n");
    }

    // For each atom...since we're dumping alpha and beta we need this again
    for (atom=0; atom<mol.natom(); ++atom) {
        fprintf(molden, "  %3d 0\n", atom+1);

        // Go through all the shells on this center
        for (int shell=0; shell < basisset.nshell_on_center(atom); ++shell) {
            int overall_shell = basisset.shell_on_center(atom, shell);

            GaussianShell& gs = *basisset.shell(overall_shell).get();

            fprintf(molden, "%c %3d 1.00\n", gs.amchar(), gs.nprimitive());

            for (int prim=0; prim<gs.nprimitive(); ++prim) {
                fprintf(molden, "%20.10f %20.10f\n", gs.exp(prim), gs.coef(prim));
            }
        }

        // An empty line separates atoms
        fprintf(molden, "\n");
    }

    // Convert Ca & Cb
    // make copies
    Matrix Ca(wavefunction_->Ca());
    Matrix Cb(wavefunction_->Cb());
    Vector& Ea = *wavefunction_->epsilon_a().get();
    Vector& Eb = *wavefunction_->epsilon_b().get();

    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = sobasisset.petitelist()->aotoso();
    // need dimensions
    const Dimension aos = sobasisset.petitelist()->AO_basisdim();
    const Dimension sos = sobasisset.petitelist()->SO_basisdim();

//    Ca.print();
//    Cb.print();
//    aotoso->print();

    Matrix Ca_ao_mo("Ca AO x MO", aos, sos);
    Matrix Cb_ao_mo("Cb AO x MO", aos, sos);

    // do the half transform
    Ca_ao_mo.gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo.gemm(false, false, 1.0, aotoso, Cb, 0.0);

//    Ca_ao_mo.print();
//    Cb_ao_mo.print();

    CharacterTable ct = mol.point_group()->char_table();

    // Dump MO's to the molden file
    fprintf(molden, "[MO]\n");
    // do alpha's
    for (int h=0; h<wavefunction_->nirrep(); ++h) {
        for (int n=0; n<wavefunction_->nmopi()[h]; ++n) {
            fprintf(molden, " Sym= %s\n", ct.gamma(h).symbol());
            fprintf(molden, " Ene= %20.10f\n", Ea.get(h, n));
            fprintf(molden, " Spin= Alpha\n");
            int occ = n < (wavefunction_->doccpi()[h] + wavefunction_->soccpi()[h]) ? 1.0 : 0.0;
            fprintf(molden, " Occup= %3.1d\n", occ);
            for (int so=0; so<wavefunction_->nso(); ++so)
                fprintf(molden, "%3d %20.12f\n", so+1, Ca_ao_mo.get(h, so, n));
        }
    }

    // do beta's
    for (int h=0; h<wavefunction_->nirrep(); ++h) {
        for (int n=0; n<wavefunction_->nmopi()[h]; ++n) {
            fprintf(molden, " Sym= %s\n", ct.gamma(h).symbol());
            fprintf(molden, " Ene= %20.10f\n", Eb.get(h, n));
            fprintf(molden, " Spin= Beta\n");
            int occ = n < (wavefunction_->doccpi()[h]) ? 1.0 : 0.0;
            fprintf(molden, " Occup= %3.1d\n", occ);
            for (int so=0; so<wavefunction_->nso(); ++so)
                fprintf(molden, "%3d %20.12f\n", so+1, Cb_ao_mo.get(h, so, n));
        }
    }

    fclose(molden);
}
