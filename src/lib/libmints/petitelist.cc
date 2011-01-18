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
using namespace psi;

#if FC_SYMBOL==2
    #define F_DGESVD dgesvd_
#elif FC_SYMBOL==1
    #define F_DGESVD dgesvd
#elif FC_SYMBOL==3
    #define F_DGESVD DGESVD
#elif FC_SYMBOL==4
    #define F_DGESVD DGESVD_
#endif

extern "C" {
extern int F_DGESVD(const char *, const char *, int *, int *, double *, int *,
                    double *, double *, int *, double *, int *, double *, int *, int *);
}

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

int **compute_atom_map(const shared_ptr<BasisSet> &basis)
{
    // grab references to the Molecule and BasisSet for convenience
    BasisSet& gbs = *basis.get();
    Molecule& mol = *gbs.molecule().get();

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

void delete_atom_map(int **atom_map, const shared_ptr<BasisSet> &basis)
{
    if (atom_map) {
        int natom = basis->molecule()->natom();
        for (int i=0; i < natom; i++)
            delete[] atom_map[i];
        delete[] atom_map;
    }
}

int **compute_shell_map(int **atom_map, const shared_ptr<BasisSet> &basis)
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

void delete_shell_map(int **shell_map, const shared_ptr<BasisSet> &basis)
{
    int nshell = basis->nshell();
    if (shell_map) {
        for (int i=0; i < nshell; i++)
            delete[] shell_map[i];
        delete[] shell_map;
    }
}

////////////////////////////////////////////////////////////////////////////

PetiteList::PetiteList(const shared_ptr<BasisSet> &gbs, const shared_ptr<IntegralFactory> &ints)
    : basis_(gbs), integral_(ints)
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

    // initialize private members
    c1_=0;
    ng_ = ct.order();
    natom_ = mol.natom();
    nshell_ = gbs.nshell();
    nirrep_ = ct.nirrep();

    // if point group is C1, then zero everything
    if (ng_==1) {
        c1_=1;
        nblocks_=1;

        p1_=0;
        atom_map_=0;
        shell_map_=0;
        lamij_=0;
        nbf_in_ir_=0;
        return;
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

    // set up atom and shell mappings
    double np[3];
    SymmetryOperation so;

    // loop over all centers
    for (i=0; i < natom_; i++) {
        Vector3 ac(mol.xyz(i));

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
                int am=gbs.shell(i,s)->am();

                if (am==0)
                    red_rep[g] += 1.0;
                else {
                    ShellRotation r(am,so,integral_.get(),gbs.shell(i,s)->is_pure());
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
    int nbf = basis_->nbf();
    Dimension ret(1, "AO Basis Dimension");
    ret[0] = nbf;
    return ret;
}

Dimension PetiteList::SO_basisdim()
{
    int i, j, ii;

    // grab reference to the basis set;
    BasisSet& gbs = *basis_.get();

    // Create the character table for the point group
    CharacterTable ct = gbs.molecule()->point_group()->char_table();

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

SO_block* PetiteList::aotoso_info()
{
    int iuniq, i, j, d, ii, jj, g, s, c, ir;

    BasisSet& gbs = *basis_.get();
    Molecule& mol = *gbs.molecule().get();

    // create the character table for the point group
    CharacterTable ct = mol.point_group()->char_table();
    SymmetryOperation so;

    if (c1_) {
        SO_block *SOs = new SO_block[1];
        SOs[0].set_length(gbs.nbf());
        for (i=0; i < gbs.nbf(); i++) {
            SOs[0].so[i].set_length(1);
            SOs[0].so[i].cont[0].bfn=i;
            SOs[0].so[i].cont[0].coef=1.0;
        }
        return SOs;
    }

    // ncomp is the number of symmetry blocks we have. for point groups with
    // complex E representations, this will be cut in two eventually
    int ncomp=0;
    for (i=0; i < nirrep_; i++)
        ncomp += ct.gamma(i).degeneracy();

    // saoelem is the current SO in a block
    int *saoelem = new int[ncomp];
    memset(saoelem,0,sizeof(int)*ncomp);

    int *whichir = new int[ncomp];
    int *whichcmp = new int[ncomp];
    for (i=ii=0; i < nirrep_; i++) {
        for (int j=0; j < ct.gamma(i).degeneracy(); j++,ii++) {
            whichir[ii] = i;
            whichcmp[ii] = j;
        }
    }

    // SOs is an array of SO_blocks which holds the redundant SO's
    SO_block *SOs = new SO_block[ncomp];

    for (i=0; i < ncomp; i++) {
        ir = whichir[i];
        int len = (ct.gamma(ir).complex()) ? nbf_in_ir_[ir]/2 : nbf_in_ir_[ir];
        SOs[i].set_length(len);
    }

    // loop over all unique shells
    for (iuniq=0; iuniq < mol.nunique(); iuniq++) {
        int nequiv = mol.nequivalent(iuniq);
        i = mol.unique(iuniq);
        for (s=0; s < gbs.nshell_on_center(i); s++) {
            int shell_i = gbs.shell_on_center(i,s);

            // test to see if there are any high am cartesian functions in this
            // shell.  for now don't allow symmetry with cartesian functions...I
            // just can't seem to get them working.
            if (gbs.shell(i,s)->am() > 1 && gbs.shell(i,s)->is_cartesian()) {
                if (ng_ != nirrep_) {
                    throw PSIEXCEPTION("PetiteList::asotoso: cannot yet handle symmetry for angular momentum >= 2");
                }
            }

            // the functions do not mix between contractions
            // so the contraction loop can be done outside the symmetry
            // operation loop
            int bfn_offset_in_shell = 0;
            int nfuncuniq = gbs.shell(i,s)->nfunction();
            int nfuncall = nfuncuniq * nequiv;

            // allocate an array to store linear combinations of orbitals
            // this is destroyed by the SVD routine
            double **linorb = new double*[nfuncuniq];
            linorb[0] = new double[nfuncuniq*nfuncall];
            for (j=1; j<nfuncuniq; j++) {
                linorb[j] = &linorb[j-1][nfuncall];
            }

            // a copy of linorb
            double **linorbcop = new double*[nfuncuniq];
            linorbcop[0] = new double[nfuncuniq*nfuncall];
            for (j=1; j<nfuncuniq; j++) {
                linorbcop[j] = &linorbcop[j-1][nfuncall];
            }

            // allocate an array for the SVD routine
            double **u = new double*[nfuncuniq];
            u[0] = new double[nfuncuniq*nfuncuniq];
            for (j=1; j<nfuncuniq; j++) {
                u[j] = &u[j-1][nfuncuniq];
            }

            // loop through each irrep to form the linear combination
            // of orbitals of that symmetry
            int irnum = 0;
            for (ir=0; ir<ct.nirrep(); ir++) {
                int cmplx = (ct.complex() && ct.gamma(ir).complex());
                for (int comp=0; comp < ct.gamma(ir).degeneracy(); comp++) {
                    memset(linorb[0], 0, nfuncuniq*nfuncall*sizeof(double));
                    for (d=0; d < ct.gamma(ir).degeneracy(); d++) {
                        // if this is a point group with a complex E representation,
                        // then only do the first set of projections for E
                        if (d && cmplx) break;

                        // operate on each function in this contraction with each
                        // symmetry operation to form symmetrized linear combinations
                        // of orbitals

                        for (g=0; g<ng_; g++) {
                            // the character
                            double ccdg = ct.gamma(ir).p(comp,d,g);

                            so = ct.symm_operation(g);
                            int equivatom = atom_map_[i][g];
                            int atomoffset
                                    = gbs.molecule()->atom_to_unique_offset(equivatom);

                            ShellRotation rr
                                    = integral_->shell_rotation(gbs.shell(i,s)->am(),
                                                            so,gbs.shell(i,s)->is_pure());

                            for (ii=0; ii < rr.dim(); ii++) {
                                for (jj=0; jj < rr.dim(); jj++) {
                                    linorb[ii][atomoffset*nfuncuniq+jj] += ccdg * rr(ii,jj);
                                }
                            }
                        }
                    }
                    // find the linearly independent SO's for this irrep/component
                    memcpy(linorbcop[0], linorb[0], nfuncuniq*nfuncall*sizeof(double));
                    double *singval = new double[nfuncuniq];
                    double djunk=0; int ijunk=1;
                    int lwork = 5*nfuncall;
                    double *work = new double[lwork];
                    int info;
                    // solves At = V SIGMA Ut (since FORTRAN array layout is used)
                    F_DGESVD("N","A",&nfuncall,&nfuncuniq,linorb[0],&nfuncall,
                             singval, &djunk, &ijunk, u[0], &nfuncuniq,
                             work, &lwork, &info);
//                    C_DGESVD('A', 'N', nfuncuniq, nfuncall, linorb[0], nfuncall, singval, u[0], nfuncuniq, &djunk, ijunk, work, lwork);

                    // put the lin indep symm orbs into the so array
                    for (j=0; j<nfuncuniq; j++) {
                        if (singval[j] > 1.0e-6) {
                            SO tso;
                            tso.set_length(nfuncall);
                            int ll = 0, llnonzero = 0;
                            for (int k=0; k<nequiv; k++) {
                                for (int l=0; l<nfuncuniq; l++,ll++) {
                                    double tmp = 0.0;
                                    for (int m=0; m<nfuncuniq; m++) {
                                        tmp += u[m][j] * linorbcop[m][ll] / singval[j];
                                    }
                                    if (fabs(tmp) > 1.0e-6) {
                                        int equivatom = gbs.molecule()->equivalent(iuniq,k);
                                        tso.cont[llnonzero].bfn
                                                = l
                                                + bfn_offset_in_shell
                                                + gbs.shell_to_function(gbs.shell_on_center(equivatom,
                                                                                            s));
                                        tso.cont[llnonzero].coef = tmp;
                                        llnonzero++;
                                    }
                                }
                            }
                            tso.reset_length(llnonzero);
                            if (llnonzero == 0) {
                                throw PSIEXCEPTION("PetiteList::aotoso_info: internal error: no bfns in SO");
                            }
                            if (SOs[irnum+comp].add(tso,saoelem[irnum+comp])) {
                                saoelem[irnum+comp]++;
                            }
                            else {
                                throw PSIEXCEPTION("PetiteList::aotoso_info: internal error: impossible duplicate SO");
                            }
                        }
                    }
                    delete[] singval;
                    delete[] work;
                }
                irnum += ct.gamma(ir).degeneracy();
            }
            bfn_offset_in_shell += gbs.shell(i,s)->nfunction();

            delete[] linorb[0];
            delete[] linorb;
            delete[] linorbcop[0];
            delete[] linorbcop;
            delete[] u[0];
            delete[] u;
        }
    }

    for (i=0; i < ncomp; i++) {
        ir = whichir[i];
        int scal = ct.gamma(ir).complex() ? 2 : 1;

        if (saoelem[i] < nbf_in_ir_[ir]/scal) {
            // if we found too few, there are big problems

            fprintf(stderr, "trouble making SO's for irrep %s\n", ct.gamma(ir).symbol());
            fprintf(stderr, "  only found %d out of %d SO's\n", saoelem[i], nbf_in_ir_[ir]/scal);
            SOs[i].print("");
            throw PSIEXCEPTION("PetiteList::aotoso_info: trouble making SO's");

        } else if (saoelem[i] > nbf_in_ir_[ir]/scal) {
            // there are some redundant so's left...need to do something to get
            // the elements we want

            fprintf(stderr, "trouble making SO's for irrep %s\n", ct.gamma(ir).symbol());
            fprintf(stderr, "  found %d SO's, but there should only be %d\n", saoelem[i], nbf_in_ir_[ir]/scal);
            SOs[i].print("");
            throw PSIEXCEPTION("PetiteList::aotoso_info: trouble making SO's");
        }
    }

    if (ct.complex()) {
        SO_block *nSOs = new SO_block[nblocks_];

        int in=0;

        for (i=ir=0; ir < nirrep_; ir++) {
            if (ct.gamma(ir).complex()) {
                nSOs[in].set_length(nbf_in_ir_[ir]);
                int k;
                for (k=0; k < SOs[i].len; k++)
                    nSOs[in].add(SOs[i].so[k],k);
                i++;

                for (j=0; j < SOs[i].len; j++,k++)
                    nSOs[in].add(SOs[i].so[j],k);

                i++;
                in++;
            } else {
                for (j=0; j < ct.gamma(ir).degeneracy(); j++,i++,in++) {
                    nSOs[in].set_length(nbf_in_ir_[ir]);
                    for (int k=0; k < SOs[i].len; k++)
                        nSOs[in].add(SOs[i].so[k],k);
                }
            }
        }

        SO_block *tmp= SOs;
        SOs = nSOs;
        delete[] tmp;
    }

    delete[] saoelem;
    delete[] whichir;
    delete[] whichcmp;

    return SOs;
}

Matrix* PetiteList::aotoso()
{
    Dimension aodim = AO_basisdim();
    Dimension sodim = SO_basisdim();

    Matrix* aoso = new Matrix("AO->SO matrix", aodim, sodim);

    if (c1_) {
        aoso->identity();
        return aoso;
    }

    SO_block *sos = aotoso_info();

    // There is an SO_block for each irrep
    for (int h=0; h < nblocks(); ++h) {
        // If the block is empty, don't do anything.
        if (sodim[h] == 0)
            continue;

        SO_block& sob = sos[h];

        for (int j=0; j<sob.len; ++j) {
            SO& soj = sob.so[j];

            for (int i=0; i<soj.len; ++i) {
                int ii =  soj.cont[i].bfn;

                aoso->set(h, ii, j, soj.cont[i].coef);
            }
        }
    }

    delete[] sos;
    return aoso;
}
