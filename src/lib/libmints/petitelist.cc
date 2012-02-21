//
// petite.cc --- implementation of the PetiteList class
//
// Modified for PSI4.
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include "mints.h"

using namespace boost;

#if FC_SYMBOL==2
    #define F_DGESVD dgesvd_
#elif FC_SYMBOL==1
    #define F_DGESVD dgesvd
#elif FC_SYMBOL==3
    #define F_DGESVD DGESVD
#elif FC_SYMBOL==4
    #define F_DGESVD DGESVD_
#endif

extern int sing_(double *q, int *lq, int *iq, double *s, double *p,
     int *lp, int *ip, double *a, int *la, int *m, int *n, double *w);

extern "C" {
extern int F_DGESVD(const char *, const char *, int *, int *, double *, int *,
                    double *, double *, int *, double *, int *, double *, int *, int *);
}

namespace psi {

///////////////////////////////////////////////////////////////////////////////

contribution::contribution()
{

}

contribution::~contribution()
{

}

contribution::contribution(int b, double c) : bfn(b), coef(c)
{

}

///////////////////////////////////////////////////////////////////////////////

SO::SO() : len(0), length(0), cont(0)
{

}

SO::SO(int l) : len(0), length(0), cont(0)
{
    set_length(l);
}

SO::~SO()
{
    set_length(0);
}

SO& SO::operator =(const SO& so)
{
    set_length(so.length);
    length = so.length;
    for (int i=0; i<length; ++i)
        cont[i] = so.cont[i];
    return *this;
}

void SO::set_length(int l)
{
    len=l;
    length=l;
    if (cont) {
        delete[] cont;
        cont=0;
    }

    if (l)
        cont = new contribution[l];
}

void SO::reset_length(int l)
{
    length=l;

    if (l<=len)
        return;

    l=l+10;

    contribution *newcont = new contribution[l];

    if (cont) {
        for (int i=0; i<len; ++i)
            newcont[i] = cont[i];

        delete[] cont;
    }

    cont = newcont;
    len=l;
}

int SO::equiv(const SO& so)
{
    int i;

    if (so.length != length)
        return 0;

    double c=0;
    for (i=0; i < length; i++) {
        if (cont[i].bfn != so.cont[i].bfn)
            return 0;
        c += cont[i].coef*so.cont[i].coef;
    }

    // if the overlap == 1.0, they're equal (SO's should have been
    // normalized by now)
    if (fabs(fabs(c)-1.0) < 1.0e-3)
        return 1;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

SO_block::SO_block() : len(0), so(0)
{
}

SO_block::SO_block(int l) : len(0), so(0)
{
    set_length(l);
}

SO_block::~SO_block()
{
    set_length(0);
}

void
SO_block::set_length(int l)
{
    len=l;
    if (so) {
        delete[] so;
        so=0;
    }

    if (l)
        so = new SO[l];
}

void
SO_block::reset_length(int l)
{
    if (len == l) return;

    SO *newso = new SO[l];

    if (so) {
        for (int i=0; i < len; i++)
            newso[i] = so[i];

        delete[] so;
    }

    so=newso;
    len=l;
}

int
SO_block::add(SO& s, int i)
{
    // first check to see if s is already here
    for (int j=0; j < ((i < len) ? i : len); j++)
        if (so[j].equiv(s))
            return 0;

    if (i >= len)
        reset_length(i+1);
    so[i] = s;

    return 1;
}

void
SO_block::print(const char *title)
{
    int i,j;

    fprintf(outfile, "SO block %s\n", title);

    for (i=0; i < len; i++) {
        fprintf(outfile, "  SO %d\n", i+1);
        for (j=0; j < so[i].length; j++)
            fprintf(outfile, " %10d", so[i].cont[j].bfn);
        fprintf(outfile, "\n");

        for (j=0; j < so[i].length; j++)
            fprintf(outfile, " %10.7f", so[i].cont[j].coef);
        fprintf(outfile, "\n");
    }
}

///////////////////////////////////////////////////////////////////////////////

struct lin_comb {
    int ns;
    int f0;
    int mapf0;
    double **c;

    lin_comb(int ins, int if0, int imf0) : ns(ins), f0(if0), mapf0(imf0) {
        int i;

        c = new double*[ns];
        for (i=0; i < ns; i++) {
            c[i] = new double[ns];
            memset(c[i],0,sizeof(double)*ns);
        }
    }

    ~lin_comb() {
        if (c) {
            for (int i=0; i < ns; i++)
                if (c[i])
                    delete[] c[i];
                delete[] c;
            c=0;
        }
    }

    void print() const {
        int i;
        for (i=0; i < ns; i++)
            fprintf(outfile, " %10d", mapf0+i);
        fprintf(outfile, "\n");

        for (i=0; i < ns; i++) {
            fprintf(outfile, "%2d", f0+i);
            for (int j=0; j < ns; j++)
                fprintf(outfile, " %10.7f", c[i][j]);
            fprintf(outfile, "\n");
        }
    }
};

///////////////////////////////////////////////////////////////////////////////

int **compute_atom_map(const Molecule* molecule)
{
    // grab references to the Molecule
    const Molecule& mol = *molecule;

    // create the character table for the point group
    CharacterTable ct = mol.point_group()->char_table();

    int natom = mol.natom();
    int ng = ct.order();
    int **atom_map;
    atom_map = new int*[natom];
    for (int i=0; i < natom; i++)
        atom_map[i] = new int[ng];

    double np[3];
    SymmetryOperation so;

    // loop over all centers
    for (int i=0; i < natom; i++) {
        Vector3 ac(mol.xyz(i));

        // then for each symop in the pointgroup, transform the coordinates of
        // center "i" and see which atom it maps into
        for (int g=0; g < ng; g++) {
            so = ct.symm_operation(g);

            for (int ii=0; ii < 3; ii++) {
                np[ii] = 0;
                for (int jj=0; jj < 3; jj++)
                    np[ii] += so(ii,jj) * ac[jj];
            }

            atom_map[i][g] = mol.atom_at_position1(np, 0.05);
            if (atom_map[i][g] < 0) {
                fprintf(outfile, "ERROR: Symmetry operation %d did not map atom %d to another atom:\n", g, i+1);
                fprintf(outfile, "  Molecule:\n");
                mol.print();
                fprintf(outfile, "  attempted to find atom at\n");
                fprintf(outfile, "    %lf %lf %lf\n", np[0], np[1], np[2]);
                abort();
            }
        }
    }

    return atom_map;
}

int **compute_atom_map(const boost::shared_ptr<Molecule> &molecule)
{
    return compute_atom_map(molecule.get());
}

void delete_atom_map(int **atom_map, const Molecule* molecule)
{
    if (atom_map) {
        int natom = molecule->natom();
        for (int i=0; i < natom; i++)
            delete[] atom_map[i];
        delete[] atom_map;
    }
}

void delete_atom_map(int **atom_map, const boost::shared_ptr<Molecule> &molecule)
{
    delete_atom_map(atom_map, molecule.get());
}

int **compute_shell_map(int **atom_map, const boost::shared_ptr<BasisSet> &basis)
{
    int **shell_map;

    BasisSet& gbs = *basis.get();
    Molecule& mol = *gbs.molecule().get();

    // create the character table for the point group
    CharacterTable ct = mol.point_group()->char_table();

    int natom = mol.natom();
    int ng = ct.order();

    int nshell = basis->nshell();
    shell_map = new int*[nshell];
    for (int i=0; i < nshell; i++)
        shell_map[i] = new int[ng];

    for (int i=0; i<natom; i++) {
        // hopefully, shells on equivalent centers will be numbered in the same
        // order
        for (int s=0; s < gbs.nshell_on_center(i); s++) {
            int shellnum = gbs.shell_on_center(i,s);
            for (int g=0; g < ng; g++) {
                shell_map[shellnum][g] = gbs.shell_on_center(atom_map[i][g],s);
            }
        }
    }

    return shell_map;
}

void delete_shell_map(int **shell_map, const boost::shared_ptr<BasisSet> &basis)
{
    int nshell = basis->nshell();
    if (shell_map) {
        for (int i=0; i < nshell; i++)
            delete[] shell_map[i];
        delete[] shell_map;
    }
}

////////////////////////////////////////////////////////////////////////////

PetiteList::PetiteList(const boost::shared_ptr<BasisSet> &gbs, const boost::shared_ptr<IntegralFactory> &ints, bool include_pure_transform)
    : basis_(gbs), integral_(ints.get()), include_pure_transform_(include_pure_transform)
{
    init();
}

PetiteList::PetiteList(const boost::shared_ptr<BasisSet> &gbs, const IntegralFactory* ints, bool include_pure_transform)
    : basis_(gbs), integral_(ints), include_pure_transform_(include_pure_transform)
{
    init();
}

PetiteList::~PetiteList()
{
    if (p1_)
        delete[] p1_;

    if (lamij_)
        delete[] lamij_;

    if (nbf_in_ir_)
        delete[] nbf_in_ir_;

    if (atom_map_) {
        for (int i=0; i < natom_; i++)
            delete[] atom_map_[i];
        delete[] atom_map_;
    }

    if (shell_map_) {
        for (int i=0; i < nshell_; i++)
            delete[] shell_map_[i];
        delete[] shell_map_;
    }

    if (unique_shell_map_) {
        for (int i=0; i < nunique_shell_; i++)
            delete[] unique_shell_map_[i];
        delete[] unique_shell_map_;
    }

    if (stablizer_)
        delete[] stablizer_;

    natom_=0;
    nshell_=0;
    ng_=0;
    nblocks_=0;
    nirrep_=0;
    p1_=0;
    atom_map_=0;
    shell_map_=0;

    lamij_=0;
    nbf_in_ir_=0;
}

boost::shared_ptr<PetiteList> PetiteList::clone()
{
    return boost::shared_ptr<PetiteList>(new PetiteList(basis_, integral_));
}

int PetiteList::nfunction(int i) const
{
    return (c1_) ? basis_->nbf() : nbf_in_ir_[i];
}

void PetiteList::init()
{
    int i;

    // grab references to the Molecule and BasisSet for convenience
    BasisSet& gbs = *basis_.get();
    Molecule& mol = *gbs.molecule().get();

    // create the character table for the point group
    CharacterTable ct = mol.point_group()->char_table();

    group_ = ct.bits();

    // initialize private members
    c1_=0;
    ng_ = ct.order();
    natom_ = mol.natom();
    nshell_ = gbs.nshell();  // full number of shells
    nirrep_ = ct.nirrep();

    // if point group is C1, then zero everything
    if (ng_==1) {
        c1_=1;
        nblocks_=1;

        p1_=0;
        atom_map_=0;
        shell_map_=0;
        unique_shell_map_=0;
        lamij_=0;
        nbf_in_ir_=0;
        stablizer_=0;
        return;
    }

    // count the number of so shells
    nunique_shell_ = 0;
    for (i=0; i<mol.nunique(); i++) {
        nunique_shell_ += basis_->nshell_on_center(mol.unique(i));
    }

    // allocate storage for arrays
    p1_ = new char[nshell_];
    lamij_ = new char[i_offset64(nshell_)];

    atom_map_ = new int*[natom_];
    for (i=0; i < natom_; i++)
        atom_map_[i] = new int[ng_];

    shell_map_ = new int*[nshell_];
    for (i=0; i < nshell_; i++)
        shell_map_[i] = new int[ng_];

    unique_shell_map_ = new int*[nunique_shell_];
    for (i=0; i < nunique_shell_; i++)
        unique_shell_map_[i] = new int[ng_];

    stablizer_ = new unsigned short[natom_];

    // set up atom and shell mappings
    double np[3];
    SymmetryOperation so;

    max_stablizer_ = nirrep_ / mol.max_nequivalent();

    // loop over all centers
    for (i=0; i < natom_; i++) {
        Vector3 ac(mol.xyz(i));

        stablizer_[i] = 0;

        // then for each symop in the pointgroup, transform the coordinates of
        // center "i" and see which atom it maps into
        for (int g=0; g < ng_; g++) {
            so = ct.symm_operation(g);

            for (int ii=0; ii < 3; ii++) {
                np[ii] = 0;
                for (int jj=0; jj < 3; jj++)
                    np[ii] += so(ii,jj) * ac[jj];
            }

            atom_map_[i][g] = mol.atom_at_position1(np, 0.05);

            // We want the list of operations that keeps the atom the same that is not E.
            if (atom_map_[i][g] == i)
                stablizer_[i] |= so.bit();

            if (atom_map_[i][g] < 0) {
                fprintf(outfile, "ERROR: Symmetry operation %d did not map atom %d to another atom:\n", g, i+1);
                fprintf(outfile, "  Molecule:\n");
                mol.print();
                fprintf(outfile, "  attempted to find atom at\n");
                fprintf(outfile, "    %lf %lf %lf\n", np[0], np[1], np[2]);
                abort();
            }
        }

        // hopefully, shells on equivalent centers will be numbered in the same
        // order
        for (int s=0; s < gbs.nshell_on_center(i); s++) {
            int shellnum = gbs.shell_on_center(i,s);
            for (int g=0; g < ng_; g++) {
                shell_map_[shellnum][g] = gbs.shell_on_center(atom_map_[i][g],s);
            }
        }
    }

    int ushell=0;
    for (int i=0; i<mol.nunique(); ++i) {
        int atom = mol.unique(i);
        for (int s=0; s<gbs.nshell_on_center(atom); ++s) {
            for (int g=0; g<ng_; ++g) {
                unique_shell_map_[ushell][g] = gbs.shell_on_center(atom_map_[atom][g], s);
            }

            ++ushell;
        }
    }

    memset(p1_,0,nshell_);
    memset(lamij_,0,i_offset64(nshell_));

    // now we do p1_ and lamij_
    for (i=0; i < nshell_; i++) {
        int g;

        // we want the highest numbered shell in a group of equivalent shells
        for (g=0; g < ng_; g++)
            if (shell_map_[i][g] > i)
                break;

        if (g < ng_)
            continue;

        // i is in the group P1
        p1_[i] = 1;

        for (int j=0; j <= i; j++) {
            int ij = i_offset64(i)+j;
            int nij = 0;

            // test to see if IJ is in the group P2, if it is, then set lambda(ij)
            // equal to the number of equivalent shell pairs.  This number is
            // just the order of the group divided by the number of times ij is
            // mapped into itself
            int gg;
            for (gg=0; gg < ng_; gg++) {
                int gi = shell_map_[i][gg];
                int gj = shell_map_[j][gg];
                int gij = ij_offset64(gi,gj);
                if (gij > ij)
                    break;
                else if (gij == ij)
                    nij++;
            }

            if (gg < ng_)
                continue;

            lamij_[ij] = (char) (ng_/nij);
        }
    }

    // form reducible representation of the basis functions
    double *red_rep = new double[ng_];
    memset(red_rep,0,sizeof(double)*ng_);

    for (i=0; i < natom_; i++) {
        for (int g=0; g < ng_; g++) {
            so = ct.symm_operation(g);
            int j= atom_map_[i][g];

            if (i!=j)
                continue;

            for (int s=0; s < gbs.nshell_on_center(i); s++) {
                int am=gbs.shell(i,s).am();

                if (am==0)
                    red_rep[g] += 1.0;
                else {
                    ShellRotation r(am,so,integral_,gbs.shell(i,s).is_pure());
                    red_rep[g] += r.trace();
                }
            }
        }
    }

    // and then use projection operators to figure out how many SO's of each
    // symmetry type there will be
    nblocks_ = 0;
    nbf_in_ir_ = new int[nirrep_];
    for (i=0; i < nirrep_; i++) {
        double t=0;
        for (int g=0; g < ng_; g++)
            t += ct.gamma(i).character(g)*red_rep[g];

        nbf_in_ir_[i] = ((int) (t+0.5))/ng_;
        if (ct.gamma(i).complex()) {
            nblocks_++;
            nbf_in_ir_[i] *= 2;
        } else {
            nblocks_ += ct.gamma(i).degeneracy();
        }
    }

    delete[] red_rep;
}

Dimension PetiteList::AO_basisdim()
{
    int one = 1;
    int nbf = include_pure_transform_ ? basis_->nao() : basis_->nbf();
    Dimension ret(1, "AO Basis Dimension");
    ret[0] = nbf;
    return ret;
}

Dimension PetiteList::SO_basisdim()
{
    int i, j, ii;

    // grab reference to the basis set;
    BasisSet& gbs = *basis_.get();

    // ncomp is the number of symmetry blocks we have
    int ncomp = nblocks();

    Dimension ret(ncomp, "SO Basis Dimension");

    for (i=0; i<nirrep_; ++i) {
        int nbas = (c1_) ? gbs.nbf() : nbf_in_ir_[i];
        ret[i] = nbas;
    }

    return ret;
}

void PetiteList::print(FILE *out)
{
    int i;

    fprintf(out, "PetiteList:\n");

    if (c1_) {
        fprintf(out, "  is c1\n");
        return;
    }

    fprintf(out, "  natom_ = %d\n", natom_);
    fprintf(out, "  nshell_ = %d\n", nshell_);
    fprintf(out, "  ng_ = %d\n", ng_);
    fprintf(out, "  nirrep_ = %d\n", nirrep_);

    fprintf(out, "  atom_map_ = \n");
    for (i=0; i<natom_; ++i) {
        fprintf(out, "    ");
        for (int g=0; g<ng_; ++g)
            fprintf(out, "%5d ", atom_map_[i][g]);
        fprintf(outfile, "\n");
    }

    fprintf(out, "  stabilizer_ =\n");
    for (i=0; i<natom_; ++i)
        fprintf(out, "    %5d %5d\n", i, stablizer_[i]);

    fprintf(out, "  shell_map_ = \n");
    for (i=0; i<nshell_; ++i) {
        fprintf(out, "    ");
        for (int g=0; g<ng_; ++g)
            fprintf(out, "%5d ", shell_map_[i][g]);
        fprintf(outfile, "\n");
    }

    fprintf(out, "  p1_ =\n");
    for (i=0; i<nshell_; ++i)
        fprintf(out, "    %5d\n", p1_[i]);

    fprintf(out, "  lamij_ = \n");
    for (i=0; i<nshell_; ++i) {
        fprintf(out, "    ");
        for (int j=0; j<=i; ++j)
            fprintf(out, "%5d ", lamij_[i_offset64(i)+j]);
        fprintf(outfile, "\n");
    }

    fprintf(out, "\n");

    CharacterTable ct = basis_->molecule()->point_group()->char_table();
    for (i=0; i<nirrep_; ++i)
        fprintf(out, "%5d functions of %s symmetry\n", nbf_in_ir_[i], ct.gamma(i).symbol());
}

/**
 * This function forms the mapping info from Cartesian AOs, to symmetry adapted pure (or Cartesian if the
 * basis requires this) functions, storing the result in a sparse buffer.
 * @param include_cart_to_pure whether to fold the spherical transform coefficients in or not.  IF true, a
 *        pure/Cartesian AOs to pure/Cartesian (depending on the basis) SOs transformation is returned.
 *        If false, a Cartesian AOs to Cartesian SOs is returned.
 * @return A pointer to the newly-created sparse SO_Block (remember to delete it!).
 */
SO_block*
PetiteList::compute_aotoso_info()
{
    bool to_pure   = include_pure_transform_ && basis_->has_puream();
    bool from_cart = include_pure_transform_ || !basis_->has_puream();

    shared_ptr<Molecule> mol = basis_->molecule();
    CharacterTable ct = mol->point_group()->char_table();
    int nunique = mol->nunique();
    int maxam = basis_->max_am();
    int **atom_map = compute_atom_map(mol);
    unsigned int functions_per_irrep[8];
    SO_block *SOs = new SO_block[nirrep_];
    for(int h = 0; h < nirrep_; ++h){
        SOs[h].set_length(nfunction(h));
        functions_per_irrep[h] = 0;
    }
    double*** function_parities = new double**[nirrep_];
    for(int symop = 0; symop < nirrep_; ++symop){
        function_parities[symop] = new double*[maxam+1];
        SymmetryOperation so = ct.symm_operation(symop);
        // How does this symmetry operation affect the x, y and z coordinates?
        bool op_inverts_x = so(0, 0) < 0.0;
        bool op_inverts_y = so(1, 1) < 0.0;
        bool op_inverts_z = so(2, 2) < 0.0;
        for(int am = 0; am <= maxam; ++am){
            // This is always the number of Cartesian functions
            int nfunctions = from_cart ? (am + 1) * (am + 2) / 2 :  2 * am + 1;
            function_parities[symop][am] = new double[nfunctions];
            int bf = 0;
            if(from_cart){
                CartesianIter cart_it(am);
                for(cart_it.start(); cart_it; cart_it.next()){
                    // What's the parity of this basis function, under the symmetry operation?
                    bool x_is_odd = cart_it.a()%2;
                    bool y_is_odd = cart_it.b()%2;
                    bool z_is_odd = cart_it.c()%2;
                    double x = x_is_odd && op_inverts_x ? -1.0 : 1.0;
                    double y = y_is_odd && op_inverts_y ? -1.0 : 1.0;
                    double z = z_is_odd && op_inverts_z ? -1.0 : 1.0;
                    function_parities[symop][am][bf] = x * y * z;
                    ++bf;
                }
            }else{
                const SphericalTransform &trans = *integral_->spherical_transform(am);
                SphericalTransformIter iter(trans);
                for(iter.first(); !iter.is_done(); iter.next()){
                    int pure = iter.pureindex();
                    if(pure == bf){
                        // Only process this the first time we encounter each pure index
                        bool x_is_odd = iter.a()%2;
                        bool y_is_odd = iter.b()%2;
                        bool z_is_odd = iter.c()%2;
                        double x = x_is_odd && op_inverts_x ? -1.0 : 1.0;
                        double y = y_is_odd && op_inverts_y ? -1.0 : 1.0;
                        double z = z_is_odd && op_inverts_z ? -1.0 : 1.0;
//                        fprintf(outfile, "l = %d, Setting functionparities[%d][%d][%d] = %f\n", am, symop, am, bf, x*y*z);
                        function_parities[symop][am][bf] = x * y * z;
                        ++bf;
                    }
                }
            }
            if(bf != nfunctions){
                std::stringstream err;
                err << "form_ao_to_so_info(): BF count problem, expected " << nfunctions
                    << " symmetry adapted functions, but found " << bf;
                throw PSIEXCEPTION(err.str());
            }
        }
    }

    for(int uatom = 0; uatom < nunique; ++uatom){
        int atom = mol->unique(uatom);
        int nimages = mol->nequivalent(uatom);
        double norm = 1.0 / sqrt(nimages);
        int nshells = basis_->nshell_on_center(atom);
        for(int shell = 0; shell < nshells; ++shell){
//fprintf(outfile, "working on shell %d\n", shell);fflush(outfile);
            int abs_shell = basis_->shell_on_center(atom, shell);
            int ncart = basis_->shell(abs_shell).ncartesian();
            int npure = basis_->shell(abs_shell).nfunction();
            int src_nfunc = to_pure ? ncart : npure;
            int src_dimension = nimages * src_nfunc;
            int l = basis_->shell(abs_shell).am();

            // Store the coefficients for each symmetry
            SOCoefficients *coefficients_list;
            coefficients_list = new SOCoefficients[src_dimension];
            size_t salc_count = 0;
            for(int bf = 0; bf < src_nfunc; ++bf){
                size_t so_count = 0;
                for(int h = 0; h < nirrep_; ++h){
                    SOCoefficients coefficients;
                    IrreducibleRepresentation gamma = ct.gamma(h);
                    // The number of stabilizers, i.e. operations that don't move the atom
                    int nstab = 0;
                    // First, figure out the symmetrization, by applying symmetry operations.
                    for(int symop = 0; symop < nirrep_; ++symop){
                        int mapped_atom = atom_map[atom][symop];
                        int mapped_shell = basis_->shell_on_center(mapped_atom, shell);
                        int mapped_bf = (to_pure ? basis_->shell_to_ao_function(mapped_shell) + bf :
                                                        basis_->shell_to_basis_function(mapped_shell) + bf);
                        if (mapped_atom == atom) ++nstab;
                        coefficients.add_contribution(mapped_bf,
                                                      function_parities[symop][l][bf] * norm * gamma.character(symop),
                                                      h);
                    }
                    coefficients.delete_zeros();
                    if(coefficients.size()){
//coefficients.print();fflush(outfile);
                        // Normalize the SO
                        coefficients.scale_coefficients(1.0/(double)nstab);
                        // We've found a non-zero contribution
                        if(to_pure){
                            // Pure, contract with the Cart->Pure transform and add to the list
                            std::map<int, double>::const_iterator coef_iter;
                            std::map<int, double>::const_iterator stop = coefficients.coefficients.end();
                            int irrep = coefficients.irrep;
                            const SphericalTransform &trans = *integral_->spherical_transform(l);
                            SphericalTransformIter cart_iter(trans);
                            for(cart_iter.first(); !cart_iter.is_done(); cart_iter.next()){
                                int cart = cart_iter.cartindex();
                                if(cart == bf){
                                    int pure = cart_iter.pureindex();
                                    for(coef_iter = coefficients.coefficients.begin(); coef_iter!=stop; ++coef_iter){
                                        size_t address = nimages * pure + so_count;
                                        double val = coef_iter->second * cart_iter.coef();
//fprintf(outfile, "l %d C %d P %d V %f v %f\n", l, cart, pure, trans->coef(n), val);
//fprintf(outfile, "Adding %d, %f, %d to %d\n",coef_iter->first, val, irrep, address);fflush(outfile);
                                        coefficients_list[address].add_contribution(coef_iter->first, val, irrep);
                                    }
                                }
                            }
                            ++so_count;
                        }else{
                            // Cartesian, all we do is add it to the list
//                            coefficients.print();
//                            fprintf(outfile, "\n\n");
                            coefficients_list[salc_count] = coefficients;
                        }
                        ++salc_count;
                    }
                }
            }
            // Sanity check
            int expected = from_cart ? nimages * ncart : nimages * npure;
            if(salc_count != expected){
                std::stringstream err;
                err << "form_ao_to_so_info(): Expected " << expected
                    << " symmetry adapted functions, but found " << salc_count;
                throw PSIEXCEPTION(err.str());
            }
            for(int n = 0; n < nimages * npure; ++n){
                SO so;
                so.set_length(coefficients_list[n].size());
                int irrep = coefficients_list[n].irrep;
                std::map<int, double>::const_iterator iter;
                std::map<int, double>::const_iterator stop = coefficients_list[n].coefficients.end();
                salc_count = 0;
                for(iter = coefficients_list[n].coefficients.begin(); iter != stop; ++iter){
                    so.cont[salc_count].bfn  = iter->first;
                    so.cont[salc_count].coef = iter->second;
                    ++salc_count;
                }
                if (SOs[irrep].add(so, functions_per_irrep[irrep])) {
                    ++functions_per_irrep[irrep];
                }else{
                    throw PSIEXCEPTION("PetiteList::aotoso_info: internal error: impossible duplicate SO");
                }
            }

            delete[] coefficients_list;
        }
    }
    for(int h = 0; h < nirrep_; ++h){
//fprintf(outfile, "Coeffs for irrep %d\n",h);
//SOs[h].print("");
        if(!c1_ && functions_per_irrep[h] != nbf_in_ir_[h] && include_pure_transform_){
            std::stringstream err;
            err << "PetiteList::aotoso_info(): In irrep " << h << " found " <<
                functions_per_irrep[h] << " SOs, but expected " << nbf_in_ir_[h];
            throw PSIEXCEPTION(err.str());
        }
    }

    // Free the cached temporaries.
    delete_atom_map(atom_map, mol);
    for(int h = 0; h < nirrep_; ++h){
        for(int am = 0; am <= maxam; ++am){
            delete [] function_parities[h][am];
        }
        delete [] function_parities[h];
    }
    delete [] function_parities;

    return SOs;
}

SharedMatrix PetiteList::sotoao()
{
    return SharedMatrix(aotoso()->transpose());
}

SharedMatrix PetiteList::aotoso()
{
    Dimension aodim = AO_basisdim();
    Dimension sodim = SO_basisdim();

    SharedMatrix aoso(new Matrix("AO->SO matrix", aodim, sodim));

//    if (c1_) {
//        aoso->identity();
//        return aoso;
//    }

    SO_block* SOs = compute_aotoso_info();

    // There is an SO_block for each irrep
    for (int h=0; h < nblocks(); ++h) {
        // If the block is empty, don't do anything.
        if (sodim[h] == 0)
            continue;

        SO_block& sob = SOs[h];

        for (int j=0; j<sob.len; ++j) {
            SO& soj = sob.so[j];

            for (int i=0; i<soj.len; ++i) {
                int ii =  soj.cont[i].bfn;

                aoso->set(h, ii, j, soj.cont[i].coef);
            }
        }
    }

    delete[] SOs;
    return aoso;
}

SharedMatrix PetiteList::evecs_to_AO_basis(SharedMatrix soevecs)
{
    // if C1, then do nothing
    if (c1_)
        return SharedMatrix(new Matrix(soevecs));

    SharedMatrix result(new Matrix(soevecs->name(), AO_basisdim(), SO_basisdim()));

    // Currently expects the caller to transpose the soevecs.
    // This is reasonable because soevecs will be used in the transposed state elsewhere.
    result->gemm(false, false, 1.0, aotoso(), soevecs, 0.0);

    return result;
}

const char *labels[] = {
    " E ",
    "C2z",
    "C2y",
    "C2x",
    " i ",
    "Sxy",
    "Sxz",
    "Syz",
    " E "
};

void PetiteList::print_group(unsigned short group) const
{
    fprintf(outfile, "(group_ %d group %d) ", group_, group);
    fprintf(outfile, "%s ", labels[0]);
    for(int op = 1; op < 9; ++op){
        if (group & (1 << (op-1)))
            fprintf(outfile, "%s ", labels[op]);
    }
    fprintf(outfile, "\n");
}


} // end namespace psi
