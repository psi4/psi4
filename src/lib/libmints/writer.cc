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

NBOWriter::NBOWriter(boost::shared_ptr<Wavefunction> wavefunction)
    : wavefunction_(wavefunction)
{


}

void NBOWriter::write(const std::string &filename)
{
    FILE *file47 = fopen(filename.c_str(), "a");

    //Get the basis set and molecule from the wavefuntion
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule& mol = *basisset.molecule().get();

    //print $GENNBO section of file
    //BOHR indicates atomic units for the coordinates
    //OPEN indicates that we'll provide separate alpha and beta matrices
    fprintf(file47, "$GENNBO NATOMS = %d NBAS = %d BOHR OPEN $END\n", mol.natom(), basisset.nbf());

    //print NBO section of file47
    fprintf(file47, "$NBO       $END\n");

    //Now print out the molecule
    fprintf(file47, "$COORD\n");
    fprintf(file47, "GENNBO expects one comment line here.\n");
    for( int i =0; i< mol.natom(); i++)
        {
            //the second mol.Z() should be modified when pseudopotentials are implemented
            fprintf(file47, "%2d  %2d  %20.12f %20.12f %20.12f\n", mol.Z(i), mol.Z(i), mol.x(i), mol.y(i), mol.z(i));
        }
    fprintf(file47, "$END\n");


    //To form the BASIS and CONTRACT sections, we need some memory
    int nshells = basisset.nshell(); //Total number of shells
    int nprim = basisset.nprimitive(); //Total number of primitives
    Vector centers(basisset.nbf());
    Vector labels(basisset.nbf());
    Vector components(nshells); //Functions per shell
    Vector angmom(nshells); //Angular momentum of shell
    Vector nprimitives(nshells); //Primitives per shell
    Vector exponents(nprim); //Exponents of primitives
    //Coefficient matrix; first row is S, second P, third D, fourth F
    Matrix coefficient(4, nprim);
    coefficient.zero();
    int fnindex = 0;
    int primindex = 0;

    //Loop over all the shells
    for( int i =0; i < nshells; i++)
        {
            shared_ptr<GaussianShell> gshell(basisset.shell(i));
            int nfns = gshell->nfunction(); //get number of functions in shell
            components.set(0, i, (double)nfns);
            int angm = gshell->am(); //get angular momentum of shell
            angmom.set(0, i, (double)angm);
            for( int j = 0; j< nfns; j++)
                {
                    centers.set (0, fnindex, (double)gshell->ncenter());
                    if(gshell->is_pure())
                        labels.set (0, fnindex, angm*100+51+j);
                    else
                        labels.set (0, fnindex, angm*100+j+1);
                    fnindex++;
                }
            int nshellprim = gshell->nprimitive();
            nprimitives.set (0, i, (double)nshellprim);
            for( int k =0; k < nshellprim; k++)
                {
                    exponents.set(0, primindex, gshell->exp(k));
                    coefficient.set (0, angm, primindex, gshell->coef(k));
                    primindex++;
                }
        }

    //Now, we print out the basis section
    fprintf(file47, "$BASIS\n");
    //The CENTER section
    fprintf(file47, "CENTER = ");
    for( int i =0; i < basisset.nbf(); i++)
        {
            fprintf(file47, "%d, ", (int)centers.get(0, i)+1);
            if((i+1)%10 == 0)
                fprintf(file47, "\n");
        }

    //The LABEL section
    fprintf(file47, "LABEL = ");
    for( int i =0; i < basisset.nbf(); i++)
        {
            fprintf(file47, "%d, ", (int)labels.get(0, i));
            if((i+1)%10 == 0)
                fprintf(file47, "\n");
        }
    fprintf(file47, "\n$END\n");

    //The CONTRACT heading
    fprintf(file47, "$CONTRACT \n");
    fprintf(file47, "NSHELL = %d\n", nshells);
    fprintf(file47, "NEXP = %d\n", nprim);

    //List of the number of functions per shell
    fprintf(file47, "NCOMP = ");
    for(int i = 0; i < nshells; i++)
        {
            fprintf(file47, "%d, ", (int)components.get(0, i));
            if((i+1)%10 == 0)
                fprintf(file47, "\n");
        }
    fprintf(file47, "\nNPRIM = ");
    for(int i =0; i < nshells; i++)
        {
            fprintf(file47, "%d, ", (int)nprimitives.get(0, i));
            if((i+1)%10 == 0)
                fprintf(file47, "\n");
        }

    fprintf(file47, "\nNPTR = ");
    int ptr = 1;
    for( int i =0; i < nshells; i++)
        {
            fprintf(file47, "%d, ", ptr);
            ptr += nprimitives.get(0, i);
            if((i+1)%10 == 0)
                fprintf(file47, "\n");

        }
    fprintf(file47, "\nEXP = ");
    for( int i =0; i < nprim; i++)
        {
            fprintf(file47, "%20.10f, ", exponents.get(0, i));
            if((i+1)%4 == 0)
                fprintf(file47, "\n");
        }

    fprintf(file47, "\nCS = ");
    for( int i =0; i < nprim; i++)
        {
            fprintf(file47, "%20.10f, ", coefficient.get (0, 0, i));
            if((i+1)%4 == 0)
                fprintf(file47, "\n");
        }

    fprintf(file47, "\nCP = ");
    for( int i =0; i < nprim; i++)
        {
            fprintf(file47, "%20.10f, ", coefficient.get (0, 1, i));
            if((i+1)%4 == 0)
                fprintf(file47, "\n");
        }

    fprintf(file47, "\nCD = ");
    for( int i =0; i < nprim; i++)
        {
            fprintf(file47, "%20.10f, ", coefficient.get (0, 2, i));
            if((i+1)%4 == 0)
                fprintf(file47, "\n");
        }

    fprintf(file47, "\nCF = ");
    for( int i =0; i < nprim; i++)
        {
            fprintf(file47, "%20.10f, ", coefficient.get (0, 3, i));
            if((i+1)%4 == 0)
                fprintf(file47, "\n");
        }
    fprintf(file47, "\n$END");

    //Matrix transformation information we'll need
    SharedMatrix sotoao = wavefunction_->sobasisset ()->petitelist()->sotoao();
    int nbf = basisset.nbf ();

    //Now we need the overlap matrix in the AO basis
    MintsHelper helper;
    SharedMatrix overlap = helper.ao_overlap();
    //Print overlap matrix
    fprintf(file47, "\n$OVERLAP \n");
    for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
                {
                    fprintf(file47, "%20.10f, ", overlap->get (0, i, j));
                    if(((nbf*i+j+1)%4)==0)
                        fprintf(file47, "\n");
                }
        }
    fprintf(file47, "\n$END");


    //Alpha Density Matrix
    SharedMatrix soalphadens = wavefunction_->Da();
    SharedMatrix alphadens(new Matrix(nbf, nbf));
    alphadens->remove_symmetry (soalphadens, sotoao);
    //Beta density
    SharedMatrix betadens(new Matrix(nbf, nbf));
    if(wavefunction_->restricted ())
        betadens->copy(alphadens);
    else
        {
        SharedMatrix sobetadens = wavefunction_->Db();
        betadens->remove_symmetry (sobetadens, sotoao);
        }
    //Now print the density matrix
    fprintf(file47, "\n$DENSITY\n ");
    for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
                {
                    fprintf(file47, "%20.10f, ", alphadens->get (0, i, j));
                    if(((nbf*i+j+1)%4)==0)
                        fprintf(file47, "\n");
                }
        }
    for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
                {
                    fprintf(file47, "%20.10f, ", betadens->get (0, i, j));
                    if(((nbf*i+j+1)%4)==0)
                        fprintf(file47, "\n");
                }
        }
    fprintf(file47, "\n$END");


    //Alpha Fock Matrix
    SharedMatrix alphasofock = wavefunction_->Fa();
    SharedMatrix alphafock(new Matrix(nbf, nbf));
    alphafock->remove_symmetry (alphasofock, sotoao);
    //Beta Fock
    SharedMatrix betafock(new Matrix(nbf, nbf));
    if(wavefunction_->restricted ())
        betafock->copy(alphafock);
    else
        {
            SharedMatrix betasofock = wavefunction_->Fb();
            betafock->remove_symmetry(betasofock, sotoao);
        }
    //Print the Fock matrix
    fprintf(file47, "\n$FOCK\n ");
    for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
                {
                    fprintf(file47, "%20.10f, ", alphafock->get (0, i, j));
                    if(((nbf*i+j+1)%4)==0)
                        fprintf(file47, "\n");
                }
        }
    for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
                {
                    fprintf(file47, "%20.10f, ", betafock->get (0, i, j));
                    if(((nbf*i+j+1)%4)==0)
                        fprintf(file47, "\n");
                }
        }
    fprintf(file47, "\n$END");

    //Alpha AO->MO transformation
    SharedMatrix soalphac = wavefunction_->Ca();
    SharedMatrix alphac(new Matrix(nbf, nbf));
    alphac->gemm(true, false, 1.00, sotoao, soalphac, 0.00);
    //Beta AO->MO transformation
    SharedMatrix betac(new Matrix(nbf, nbf));
    if(wavefunction_->restricted ())
        betac->copy (alphac);
    else
        {
            SharedMatrix sobetac = wavefunction_->Cb();
            betac->gemm(true, false, 1.00, sotoao, sobetac, 0.00);
        }
    //Print the AO->MO coefficients
    fprintf(file47, "\n$LCAOMO\n ");
    for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
                {
                    fprintf(file47, "%20.10f, ", alphac->get (0, i, j));
                    if(((nbf*i+j+1)%4)==0)
                        fprintf(file47, "\n");
                }
        }
    for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
                {
                    fprintf(file47, "%20.10f, ", betac->get (0,   i, j));
                    if(((nbf*i+j+1)%4)==0)
                        fprintf(file47, "\n");
                }
        }
    fprintf(file47, "\n$END\n");







    fclose(file47);


}






