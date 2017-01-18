/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

//
// pointgrp.cc
//
// Additional modifications made by Justin Turney <jturney@ccqc.uga.edu>
// for use in PSI4.
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <janssen@limitpt.com>
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
// the Free Software Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/corrtab.h"

using namespace std;
using namespace psi;

////////////////////////////////////////////////////////////////////////

CorrelationTable::CorrelationTable():
    n_(0),
    ngamma_(0),
    gamma_(0)
{
}

CorrelationTable::CorrelationTable(const std::shared_ptr<PointGroup>& group,
    const std::shared_ptr<PointGroup>& subgroup):
    n_(0),
    ngamma_(0),
    gamma_(0)
{
    int rc = initialize_table(group,subgroup);
    if (rc != 0) {
        // ExEnv::err0()
        //     << "ERROR: CorrelationTable: " << error(rc) << endl;
        // abort();
        outfile->Printf( "CorrelationTable error: %s\n", error(rc));
        outfile->Printf("group %s -> subgroup %s\n", group->symbol().c_str(), subgroup->symbol().c_str());
        throw PSIEXCEPTION("ERROR: CorrelationTable");
    }
}

CorrelationTable::~CorrelationTable()
{
    clear();
}

int
    CorrelationTable::initialize_table(const std::shared_ptr<PointGroup>& group,
    const std::shared_ptr<PointGroup>& subgroup)
{
    clear();

    group_ = group;
    subgroup_ = subgroup;

    int i, j, k, l;

    CharacterTable ct = group_->char_table();
    CharacterTable subct = subgroup_->char_table();

    n_ = ct.nirrep();
    subn_ = subct.nirrep();
    ngamma_ = new int[n_];
    gamma_ = new int*[n_];

    // CAN ONLY HANDLE NONDEGENERATE POINT GROUPS

    for (i=0; i<n_; i++) {
        ngamma_[i] = 0;
        gamma_[i] = 0;
    }

    // map the ops in the high order to low order groups
    int *so_to_subso = new int[ct.order()];
    int *subso_to_so = new int[subct.order()];
    for (i=0; i<subct.order(); i++) subso_to_so[i] = -1;
    for (i=0; i<ct.order(); i++) {
        SymmetryOperation so = ct.symm_operation(i);
        int found = 0;
        so_to_subso[i] = -1;
        for (j=0; j<subct.order(); j++) {
            SymmetryOperation subso = subct.symm_operation(j);
            double sumsquare = 0.0;
            for (k=0; k<3; k++) {
                for (l=0; l<3; l++) {
                    double diff = so(k,l)-subso(k,l);
                    sumsquare += diff*diff;
                }
            }
            if (sumsquare < 1.0e-12) {
                found++;
                //ExEnv::outn() << scprintf("symmop %d in %s is %d in %s",
                //                 i,ct.symbol(),j,subct.symbol()) << endl;
                so_to_subso[i] = j;
                subso_to_so[j] = i;
            }
        }
        if (found > 1) {
            delete[] so_to_subso;
            delete[] subso_to_so;
            return -1;
        }
    }
    for (i=0; i<subct.order(); i++) {
        if (subso_to_so[i] == -1) {
            delete[] so_to_subso;
            delete[] subso_to_so;
            return -2;
        }
    }

    // compute the correlation table
    for (i=0; i<ct.nirrep(); i++) {
        for (j=0; j<subct.nirrep(); j++) {
            double nmatch = 0.0;
            for (k=0; k<ct.order(); k++) {
                double chr = ct.gamma(i).character(k);
                if (so_to_subso[k] >= 0) {
                    double subchr = subct.gamma(j).character(so_to_subso[k]);
                    nmatch += subchr*chr;
                }
            }
            nmatch /= subct.order();
            int inmatch = (int)(nmatch+0.5);
            if (fabs(nmatch-inmatch)>1.0e-6) {
                delete[] so_to_subso;
                delete[] subso_to_so;
                return -4;
            }
            if (inmatch > 0) {
                int *newgamma = new int[ngamma_[i] + inmatch];
                memcpy(newgamma,gamma_[i],ngamma_[i]*sizeof(int));
                for (k=0; k<inmatch; k++) newgamma[ngamma_[i]+k] = j;
                ngamma_[i] += inmatch;
                delete[] gamma_[i];
                gamma_[i] = newgamma;
            }
        }
    }

    delete[] so_to_subso;
    delete[] subso_to_so;

    for (i=0; i<n(); i++) {
        int degen = ct.gamma(i).degeneracy();
        int subdegen = 0;
        for (j=0; j<ngamma(i); j++) {
            subdegen += subct.gamma(gamma(i,j)).degeneracy();
        }
        if (degen != subdegen) {
            return -3;
        }
    }

    return 0;
}

void
    CorrelationTable::clear()
{
    for (int i=0; i<n_; i++) {
        delete[] gamma_[i];
    }
    delete[] ngamma_;
    delete[] gamma_;
}

const char *
    CorrelationTable::error(int rc)
{
    if (rc == -1) {
        return "too many symop matches";
    }
    else if (rc == -2) {
        return "not a subgroup or wrong ref frame";
    }
    else if (rc == -3) {
        return "only groups with non-complex characters supported (degen sum)";
    }
    else if (rc == -4) {
        return "only groups with non-complex characters supported (nonint nirr)";
    }
    else if (rc != 0) {
        return "unknown error";
    }
    return "no error";
}

int
    CorrelationTable::degen(int i) const
{
    return group_->char_table().gamma(i).degeneracy();
}

int
    CorrelationTable::subdegen(int i) const
{
    return subgroup_->char_table().gamma(i).degeneracy();
}

// void
//     CorrelationTable::print(ostream &o) const
// {
//     o << indent
//         << "Correlation Table from "
//         << group_->symbol() << " to " << subgroup_->symbol() << ":" << endl;
//
//     CharacterTable ct = group_->char_table();
//     CharacterTable subct = subgroup_->char_table();
//
//     o << incindent;
//     for (int i=0; i<n(); i++) {
//         o << indent
//             << ct.gamma(i).symbol() << ":";
//         for (int j=0; j<ngamma(i); j++) {
//             o << indent
//                 << " " << subct.gamma(gamma(i,j)).symbol();
//         }
//         o << endl;
//     }
//     o << decindent;
// }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
