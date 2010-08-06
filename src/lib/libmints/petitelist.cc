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
                    ShellRotation r(am,so,integral_,gbs.shell(i,s)->is_pure());
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
