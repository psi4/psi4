#include "writer.h"
#include "view.h"
#include <libmints/mints.h>
#include <psi4-dec.h>

#include <cstdio>

using namespace psi;
using namespace boost;

GradientWriter::GradientWriter(boost::shared_ptr<Molecule> mol, const Matrix& grad)
    : molecule_(mol), gradient_(grad)
{
}

void GradientWriter::write(const std::string &filename)
{
    FILE *file11 = fopen(filename.c_str(), "a");
    int i;

    if (file11 == NULL)
        throw PSIEXCEPTION("GradientWriter::write: Unable to open file11.dat");

    fprintf(file11, "%-59.59s %-10.10s%-9.9s\n",
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
                mol.symbol(atom).c_str(), atom+1, static_cast<int>(mol.Z(atom)), coord[0], coord[1], coord[2]);
    }

    // Dump the basis set using code adapted from psi2molden
    fprintf(molden, "[GTO]\n");

    // For each atom
    for (atom=0; atom<mol.natom(); ++atom) {
        fprintf(molden, "  %3d 0\n", atom+1);

        // Go through all the shells on this center
        for (int shell=0; shell < basisset.nshell_on_center(atom); ++shell) {
            int overall_shell = basisset.shell_on_center(atom, shell);

            const GaussianShell& gs = basisset.shell(overall_shell);

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

            const GaussianShell& gs = basisset.shell(overall_shell);

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

    boost::shared_ptr<PetiteList> pl(new PetiteList(wavefunction_->basisset(), wavefunction_->integral()));
    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();

    SharedMatrix Ca_ao_mo(new Matrix("Ca AO x MO", aos, sos));
    SharedMatrix Cb_ao_mo(new Matrix("Cb AO x MO", aos, sos));

    // do the half transform
    Ca_ao_mo->gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo->gemm(false, false, 1.0, aotoso, Cb, 0.0);

    aotoso->print();
    Ca_ao_mo->print();
    Cb_ao_mo->print();

    // The order Molden expects
    //     P: x, y, z
    //    5D: D 0, D+1, D-1, D+2, D-2
    //    6D: xx, yy, zz, xy, xz, yz
    //
    //    7F: F 0, F+1, F-1, F+2, F-2, F+3, F-3
    //   10F: xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
    //
    //    9G: G 0, G+1, G-1, G+2, G-2, G+3, G-3, G+4, G-4
    //   15G: xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy,
    //        xxyy xxzz yyzz xxyz yyxz zzxy
    // Since Molden doesn't handle higher than g we'll just leave them be.
    int molden_cartesian_order[][15] = {
        { 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // p
        { 0, 3, 4, 1, 5, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0 },    // d
        { 0, 4, 5, 3, 9, 6, 1, 8, 7, 2, 0, 0, 0, 0, 0 },    // f
        { 0, 3, 4, 9, 12, 10, 5, 13, 14, 7, 1, 6, 11, 8, 2} // g
    };

    int nirrep = Ca_ao_mo->nirrep();
    Dimension countpi(nirrep);
    Dimension zeropi(nirrep);
    Dimension ncartpi(nirrep);

    for(int i = 0; i < basisset.nshell(); i++) {
        int am = basisset.shell(i).am();

        int ncart = basisset.shell(i).nfunction();
        if(am == 1 || (am > 0 && am < 5 && basisset.shell(i).is_cartesian())) {
            for (int h=0; h<nirrep; ++h)
                ncartpi[h] = ncart;

            View block_a(Ca_ao_mo, ncartpi, Ca_ao_mo->colspi(), countpi, zeropi);
            View block_b(Cb_ao_mo, ncartpi, Cb_ao_mo->colspi(), countpi, zeropi);

            SharedMatrix temp_a = block_a();
            SharedMatrix temp_b = block_b();

            for( int j =0; j < ncart; j++) {
                for (int h=0; h < Ca_ao_mo->nirrep(); ++h) {
                    for (int k=0; k<Ca_ao_mo->coldim(h); ++k) {
                        fprintf(outfile, "am %d\n, from %d to %d\n", am, j, countpi[h] + molden_cartesian_order[am-1][j]);
                        Ca_ao_mo->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_a->get(h, j, k));
                        Cb_ao_mo->set(h, countpi[h] + molden_cartesian_order[am-1][j], k, temp_b->get(h, j, k));
                    }
                }
            }
        }

        for (int h=0; h<nirrep; ++h)
            countpi[h] += ncart;
    }

    if (basisset.has_puream()) {
        // Tell Molden to use spherical.  5d implies 5d and 7f.
        fprintf(molden, "[5d]\n[9g]\n\n");
    }
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
                fprintf(molden, "%3d %20.12f\n", so+1, Ca_ao_mo->get(h, so, n));
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
                fprintf(molden, "%3d %20.12f\n", so+1, Cb_ao_mo->get(h, so, n));
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
    int pure_order[][7] = {
        { 1, 0, 0, 0, 0, 0, 0},      // s
        { 101, 102, 103, 0, 0, 0, 0}, // p
        // z2  xz   yz  x2-y2 xy
        { 255, 252, 253, 254, 251, 0, 0}, // d
        //z(z2-r2), x(z2-r2), y(z2-r2) z(x2-y2), xyz, x(x2-y2), y(x2-y2)
        { 351, 352, 353, 354, 355, 356, 357 } //f
    };

    MintsHelper helper;
    SharedMatrix sotoao = helper.petite_list()->sotoao();

    FILE *file47 = fopen(filename.c_str(), "a");

    //Get the basis set and molecule from the wavefuntion
    BasisSet& basisset = *wavefunction_->basisset().get();
    Molecule& mol = *basisset.molecule().get();

    //NBO can only handle up to f functions
    if( basisset.max_am () > 3)
    {
        throw PSIEXCEPTION("NBO cannot handle angular momentum above f functions. \n");
    }
    //print $GENNBO section of file
    //BOHR indicates atomic units for the coordinates
    //OPEN indicates that we'll provide separate alpha and beta matrices
    fprintf(file47, "$GENNBO NATOMS = %d NBAS = %d BOHR BODM ", mol.natom(), basisset.nbf());

    //To make this more user-friendly in the case of RHF wavefunctions...
    bool open_shell = (wavefunction_->nalpha() != wavefunction_->nbeta());
    if(open_shell)
        fprintf(file47, "OPEN $END\n");
    else
        fprintf(file47, "$END\n");

    //print NBO section of file47; user can modify this to suit their needs
    fprintf(file47, "$NBO       $END\n");

    //Now print out the molecule
    fprintf(file47, "$COORD\n");
    fprintf(file47, "GENNBO expects one comment line here. So, here's a comment line.\n");
    for( int i =0; i< mol.natom(); i++)
    {
        //the second mol.Z() should be modified when pseudopotentials are implemented
        fprintf(file47, "%2d  %2d  %20.12f %20.12f %20.12f\n",
                static_cast<int>(mol.Z(i)), static_cast<int>(mol.Z(i)),
                mol.x(i), mol.y(i), mol.z(i));
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
        const GaussianShell& gshell = basisset.shell(i);
        int nfns = gshell.nfunction(); //get number of functions in shell
        components.set(0, i, nfns);
        int angm = gshell.am(); //get angular momentum of shell
        angmom.set(0, i, angm);
        for( int j = 0; j< nfns; j++)
        {
            centers.set (0, fnindex, gshell.ncenter());
            if(gshell.is_pure()) {
                //fprintf(outfile, "fnindex %d pure_order[%d][%d] %d\n", fnindex, angm, j, pure_order[angm][j]);
                labels.set (0, fnindex, pure_order[angm][j]);
            }
            else
                labels.set (0, fnindex, angm*100+j+1);
            fnindex++;
        }
        int nshellprim = gshell.nprimitive();
        nprimitives.set (0, i, nshellprim);
        for( int k =0; k < nshellprim; k++)
        {
            exponents.set(0, primindex, gshell.exp(k));
            coefficient.set (0, angm, primindex, gshell.coef(k));
            primindex++;
        }
    }

    //Now, we print out the basis section
    fprintf(file47, "$BASIS\n");
    //The CENTER section
    fprintf(file47, "CENTER = ");
    for( int i =0; i < basisset.nbf(); i++)
    {
        fprintf(file47, "%5d      ", (int)centers.get(0, i)+1);
        if((i+1)%10 == 0)
            fprintf(file47, "\n");
    }

    //The LABEL section
    fprintf(file47, "\nLABEL = ");
    for( int i =0; i < basisset.nbf(); i++)
    {
        fprintf(file47, "%5d      ", (int)labels.get(0, i));
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
        fprintf(file47, "%5d      ", (int)components.get(0, i));
        if((i+1)%10 == 0)
            fprintf(file47, "\n");
    }
    //List the number of primitives per shell
    fprintf(file47, "\nNPRIM = ");
    for(int i =0; i < nshells; i++)
    {
        fprintf(file47, "%d      ", (int)nprimitives.get(0, i));
        if((i+1)%10 == 0)
            fprintf(file47, "\n");
    }
    //location of the first exponent for each shell
    fprintf(file47, "\nNPTR = ");
    int ptr = 1;
    for( int i =0; i < nshells; i++)
    {
        fprintf(file47, "%d      ", ptr);
        ptr += nprimitives.get(0, i);
        if((i+1)%10 == 0)
            fprintf(file47, "\n");

    }
    //The exponents
    fprintf(file47, "\nEXP = \n");
    for( int i =0; i < nprim; i++)
    {
        fprintf(file47, "%20.10f ", exponents.get(0, i));
        if((i+1)%4 == 0)
            fprintf(file47, "\n");
    }
    //Coefficients for s functions
    fprintf(file47, "\nCS = \n");
    for( int i =0; i < nprim; i++)
    {
        fprintf(file47, "%20.10f ", coefficient.get (0, 0, i));
        if((i+1)%4 == 0)
            fprintf(file47, "\n");
    }
    //Coefficients for p functions
    fprintf(file47, "\nCP = \n");
    for( int i =0; i < nprim; i++)
    {
        fprintf(file47, "%20.10f ", coefficient.get (0, 1, i));
        if((i+1)%4 == 0)
            fprintf(file47, "\n");
    }
    //coefficients for d functions
    fprintf(file47, "\nCD = \n");
    for( int i =0; i < nprim; i++)
    {
        fprintf(file47, "%20.10f ", coefficient.get (0, 2, i));
        if((i+1)%4 == 0)
            fprintf(file47, "\n");
    }
    //coefficients for f functions
    fprintf(file47, "\nCF = \n");
    for( int i =0; i < nprim; i++)
    {
        fprintf(file47, "%20.10f ", coefficient.get (0, 3, i));
        if((i+1)%4 == 0)
            fprintf(file47, "\n");
    }
    fprintf(file47, "\n$END");

    //Matrix transformation information we'll need
    int nbf = basisset.nbf ();

    //Now we need the overlap matrix in the AO basis
    SharedMatrix overlap = helper.ao_overlap();
    //Print overlap matrix
    fprintf(file47, "\n$OVERLAP \n");
    for(int i =0; i < nbf; i++)
    {
        for(int j =0; j < nbf; j++)
        {
            fprintf(file47, "%20.10f ", overlap->get (0, i, j));
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
    SharedMatrix sobetadens = wavefunction_->Db();
    betadens->remove_symmetry (sobetadens, sotoao);
    //Now print the density matrix
    fprintf(file47, "\n$DENSITY\n ");
    if(wavefunction_->same_a_b_dens ())
    {
        SharedMatrix density(new Matrix(nbf, nbf));
        density->copy (alphadens);
        density->add (betadens);
        for( int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", density->get(0, i, j));
                if(((nbf*i+j+1)%4)==0)
                    fprintf(file47, "\n");
            }
        }
    }
    else
    {
        int count = 0;
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", alphadens->get (0, i, j));
                count++;
                if(count%4 ==0)
                    fprintf(file47, "\n");
            }
        }
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", betadens->get (0, i, j));
                count++;
                if(count%4 ==0)
                    fprintf(file47, "\n");
            }
        }
    }
    fprintf(file47, "\n$END");


    //Alpha Fock Matrix
    SharedMatrix alphasofock = wavefunction_->Fa();
    SharedMatrix alphafock(new Matrix(nbf, nbf));
    alphafock->remove_symmetry (alphasofock, sotoao);
    //Print the Fock matrix
    fprintf(file47, "\n$FOCK\n ");
    if(wavefunction_->same_a_b_dens ())
    {
        for(int i = 0; i < nbf; i++)
        {
            for(int j = 0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", alphafock->get (0, i, j));
                if(((nbf*i+j+1)%4)==0)
                    fprintf(file47, "\n");
            }
        }
    }

    else
    {
        //Beta Fock
        SharedMatrix betafock(new Matrix(nbf, nbf));
        SharedMatrix betasofock = wavefunction_->Fb();
        betafock->remove_symmetry(betasofock, sotoao);
        int count=0;
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", alphafock->get (0, i, j));
                count++;
                if(count%4 ==0)
                    fprintf(file47, "\n");
            }
        }
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", betafock->get (0, i, j));
                count++;
                if(count%4 ==0)
                    fprintf(file47, "\n");
            }
        }
    }
    fprintf(file47, "\n$END");

    //Alpha AO->MO transformation
    SharedMatrix soalphac = wavefunction_->Ca();
    SharedMatrix alphac(new Matrix(nbf, nbf));
    alphac->gemm(true, false, 1.00, sotoao, soalphac, 0.00);

    fprintf(file47, "\n$LCAOMO\n ");
    if(wavefunction_->same_a_b_orbs ())
    {
        for(int i = 0; i < nbf; i++)
        {
            for(int j = 0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", alphac->get(0, i, j));
                if(((nbf*i+j+1)%4)==0)
                    fprintf(file47, "\n");
            }
        }
    }

    else
    {
        //Beta AO->MO transformation
        SharedMatrix betac(new Matrix(nbf, nbf));
        SharedMatrix sobetac = wavefunction_->Cb();
        betac->gemm(true, false, 1.00, sotoao, sobetac, 0.00);

        //Print the AO->MO coefficients
        int count = 0;
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", alphac->get (0, i, j));
                count++;
                if(count%4==0)
                    fprintf(file47, "\n");
            }
        }
        for(int i =0; i < nbf; i++)
        {
            for(int j =0; j < nbf; j++)
            {
                fprintf(file47, "%20.10f ", betac->get (0,   i, j));
                count++;
                if(count%4 ==0)
                    fprintf(file47, "\n");
            }
        }
    }
    fprintf(file47, "\n$END\n");
    fclose(file47);
}


MOWriter::MOWriter(boost::shared_ptr<Wavefunction> wavefunction,Options&options)
    : wavefunction_(wavefunction), options_(options)
{
}

void MOWriter::write()
{
    // what kind of reference?
    bool isrestricted = true;
    if (options_.get_str("REFERENCE") == "UHF")  isrestricted = false;
    if (options_.get_str("REFERENCE") == "UKS")  isrestricted = false;
    if (options_.get_str("REFERENCE") == "CUHF") isrestricted = false;

    // Get the molecule for ease
    BasisSet& basisset = *wavefunction_->basisset().get();
    SOBasisSet& sobasisset = *wavefunction_->sobasisset().get();
    Molecule & mol = *basisset.molecule().get();

    // Convert Ca & Cb
    // make copies
    Matrix Ca(wavefunction_->Ca());
    Matrix Cb(wavefunction_->Cb());
    Vector& Ea = *wavefunction_->epsilon_a().get();
    Vector& Eb = *wavefunction_->epsilon_b().get();

    boost::shared_ptr<PetiteList> pl(new PetiteList(wavefunction_->basisset(), wavefunction_->integral()));

    // get the "aotoso" transformation matrix, ao by so
    SharedMatrix aotoso = pl->aotoso();
    // need dimensions
    const Dimension aos = pl->AO_basisdim();
    const Dimension sos = pl->SO_basisdim();

    SharedMatrix Ca_ao_mo(new Matrix("Ca AO x MO", aos, sos));
    SharedMatrix Cb_ao_mo(new Matrix("Cb AO x MO", aos, sos));

    // do the half transform
    Ca_ao_mo->gemm(false, false, 1.0, aotoso, Ca, 0.0);
    Cb_ao_mo->gemm(false, false, 1.0, aotoso, Cb, 0.0);

    int nirrep = Ca_ao_mo->nirrep();

    // order orbitals in terms of energy:
    int minorb;
    nmo = 0;
    for (int h = 0; h < nirrep; h++) nmo += wavefunction_->nmopi()[h];

    map  = new int[nmo];
    bool * skip = new bool[nmo];
    for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
    for (int orb = 0; orb < nmo; orb++) {

        int count = 0;
        double minen = 1.0e9;
        for (int h = 0; h < nirrep; h++) {
            for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {

                if ( skip[count] ) {
                   count++;
                   continue;
                }
                if ( Ea.get(h,n) <= minen ) {
                   minen = Ea.get(h,n);
                   minorb = count;
                }

                count++;
                
            }
        }
        map[ orb ] = minorb;
        skip[ minorb ] = true;
    }

    // reorder orbitals:
    nso = wavefunction_->nso();
    eps = new double[nmo];
    sym = new int[nmo];
    occ = new int[nmo];
    Ca_pointer = new double[nmo * nso];
    for (int i = 0; i < nmo * nso; i++) Ca_pointer[i] = 0.0;

    int count = 0;
    int extra = isrestricted ? 1 : 0;
    for (int h = 0; h < nirrep; h++) {
        double ** Ca_old = Ca_ao_mo->pointer(h);
        for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {
            occ[ count ] = n < ( wavefunction_->doccpi()[h] + wavefunction_->soccpi()[h] ) ? 1 : 0;
            occ[ count ] += n <  wavefunction_->doccpi()[h] ? extra : 0;
            eps[ count ] = Ea.get(h,n);
            sym[ count ] = h;
            for (int mu = 0; mu < nso; mu++) {
                Ca_pointer[mu*nmo + count] = Ca_old[mu][n];
            }
            count++;

        }
    }

    // dump to output file
    fprintf(outfile,"\n");
    if ( isrestricted ) 
       fprintf(outfile,"  ==> Molecular Orbitals <==\n");
    else
       fprintf(outfile,"  ==> Alpha-Spin Molecular Orbitals <==\n");
    fprintf(outfile,"\n");

    write_mos(mol);

    // now for beta spin
    if ( !isrestricted ) {
   
   
       // order orbitals in terms of energy
       for (int orb = 0; orb < nmo; orb++) skip[orb] = false;
       for (int orb = 0; orb < nmo; orb++) {
   
           int count = 0;
           double minen = 1.0e9;
           for (int h = 0; h < nirrep; h++) {
               for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {
   
                   if ( skip[count] ) {
                      count++;
                      continue;
                   }
   
                   if ( Eb.get(h,n) <= minen ) {
                      minen = Eb.get(h,n);
                      minorb = count;
                   }
                   count++;
                   
               }
           }
           map[ orb ] = minorb;
           skip[ minorb ] = true;
       }
   
   
       // reorder orbitals:
       for (int i = 0; i < nmo * nso; i++) Ca_pointer[i] = 0.0;
       count = 0;
       for (int h = 0; h < nirrep; h++) {
           double ** Ca_old = Cb_ao_mo->pointer(h);
           for (int n = 0; n<wavefunction_->nmopi()[h]; n++) {
               occ[ count ] = n < wavefunction_->doccpi()[h] ? 1 : 0;
               eps[ count ] = Eb.get(h,n);
               sym[ count ] = h;
               for (int mu = 0; mu < nso; mu++) {
                   Ca_pointer[mu*nmo + count] = Ca_old[mu][n];
               }
               count++;
   
           }
       }
   
       // dump to output file
       fprintf(outfile,"\n");
       fprintf(outfile,"  ==> Beta-Spin Molecular Orbitals <==\n");
       fprintf(outfile,"\n");
       write_mos(mol);
    }

    delete skip;
    delete occ;
    delete sym;
    delete eps;
    delete Ca_pointer;
}

void MOWriter::write_mos(Molecule & mol){

    CharacterTable ct = mol.point_group()->char_table();

    // print mos (5 columns)
    int ncols = 5;
    int ncolsleft = nmo % ncols;
    int nrows = (nmo - ncolsleft ) / ncols;

    // print the full rows:
    int count = 0;
    for (int i = 0; i < nrows; i++) {

        // print blank space
        fprintf(outfile,"     ");
        // print mo number
        for (int j = 0; j < ncols; j++){
            fprintf(outfile,"%13d",count+j+1);
        }
        fprintf(outfile,"\n");
        fprintf(outfile,"\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao number
            fprintf(outfile,"%5i",mu+1);
            for (int j = 0; j < ncols; j++){

                fprintf(outfile,"%13.7lf",Ca_pointer[ mu*nmo + map[count + j] ]);
            }
            fprintf(outfile,"\n");
        }
        fprintf(outfile,"\n");
        // print energy
        fprintf(outfile," Ene ");
        for (int j = 0; j < ncols; j++){
            fprintf(outfile,"%13.7lf",eps[ map[count + j] ]);
        }
        fprintf(outfile,"\n");
        // print symmetry
        fprintf(outfile," Sym ");
        for (int j = 0; j < ncols; j++){
            fprintf(outfile,"%13s",ct.gamma(sym[map[count+j]]).symbol());
        }
        fprintf(outfile,"\n");
        // print occupancy
        fprintf(outfile," Occ ");
        for (int j = 0; j < ncols; j++){
            fprintf(outfile,"%13d",occ[map[count+j]]);
        }
        fprintf(outfile,"\n");
        fprintf(outfile,"\n");
        fprintf(outfile,"\n");
        count+=ncols;
    }

    // print the partial rows:
    if ( ncolsleft > 0 ) {
    
        // print blank space
        fprintf(outfile,"     ");
        // print mo number
        for (int j = 0; j < ncolsleft; j++){
            fprintf(outfile,"%13d",count+j+1);
        }
        fprintf(outfile,"\n");
        fprintf(outfile,"\n");
        // print orbitals
        for (int mu = 0; mu < nso; mu++) {
            // print ao number
            fprintf(outfile,"%5i",mu+1);
            for (int j = 0; j < ncolsleft; j++){
                fprintf(outfile,"%13.7lf",Ca_pointer[ mu*nmo + map[count + j] ]);
            }
            fprintf(outfile,"\n");
        }
        fprintf(outfile,"\n");
        // print energy
        fprintf(outfile," Ene ");
        for (int j = 0; j < ncolsleft; j++){
            fprintf(outfile,"%13.7lf",eps[ map[count + j] ]);
        }
        fprintf(outfile,"\n");
        fprintf(outfile," Sym ");
        for (int j = 0; j < ncolsleft; j++){
            fprintf(outfile,"%13s",ct.gamma(sym[map[count+j]]).symbol());
        }
        fprintf(outfile,"\n");
        // print occupancy
        fprintf(outfile," Occ ");
        for (int j = 0; j < ncolsleft; j++){
            fprintf(outfile,"%13d",occ[map[count+j]]);
        }
        fprintf(outfile,"\n");
        fprintf(outfile,"\n");

    }
}

