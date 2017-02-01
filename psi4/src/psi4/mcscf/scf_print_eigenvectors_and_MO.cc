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

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <cstdio>
#include "psi4/libmoinfo/libmoinfo.h"

#include "scf.h"

extern FILE* outfile;

using namespace std;

namespace psi{ namespace mcscf{

void SCF::print_eigenvectors_and_MO()
{
  typedef vector<std::pair<double, string> >::iterator vecstr_it;

  // Assumes the eigenvalues of some Fock operator
  // are in the SBlockVector epsilon
  vector<std::pair<double, string> > docc_evals;
  vector<std::pair<double, string> > actv_evals;
  vector<std::pair<double, string> > virt_evals;

  for(int h = 0; h < nirreps; ++h)
    for(int i = 0; i < docc[h]; ++i)
      docc_evals.push_back(
        make_pair(epsilon->get(h,i),moinfo_scf->get_irr_labs(h))
      );
  for(int h = 0; h < nirreps; ++h)
    for(int i = docc[h]; i < docc[h] + actv[h]; ++i)
      actv_evals.push_back(
        make_pair(epsilon->get(h,i),moinfo_scf->get_irr_labs(h))
      );
  for(int h = 0; h < nirreps; ++h)
    for(int i = docc[h] + actv[h]; i < sopi[h]; ++i)
      virt_evals.push_back(
        make_pair(epsilon->get(h,i),moinfo_scf->get_irr_labs(h))
      );

  sort(docc_evals.begin(),docc_evals.end());
  sort(actv_evals.begin(),actv_evals.end());
  sort(virt_evals.begin(),virt_evals.end());


  outfile->Printf("\n\n  =========================================================================");
  outfile->Printf("\n  Eigenvalues (Eh)");
  outfile->Printf("\n  -------------------------------------------------------------------------");

  int print_nrows = 3;
  outfile->Printf("\n  Doubly occupied orbitals");
  int printed     = 0;
  int mo          = 1;
  for(vecstr_it it = docc_evals.begin(); it!= docc_evals.end(); ++it){
    outfile->Printf("%s  %5d %13.6f %3s",printed++ % print_nrows == 0 ? "\n" : "",
                                       mo++,
                                       it->first,
                                       it->second.c_str());
  }

  if(!actv_evals.empty()){
    outfile->Printf("\n  Active orbitals");
    printed     = 0;
    for(vecstr_it it = actv_evals.begin(); it!= actv_evals.end(); ++it){
      outfile->Printf("%s  %5d %13.6f %3s",printed++ % print_nrows == 0 ? "\n" : "",
                                        mo++,
                                        it->first,
                                        it->second.c_str());
    }
  }

  outfile->Printf("\n  Unoccupied orbitals");
  printed     = 0;
  for(vecstr_it it = virt_evals.begin(); it!= virt_evals.end(); ++it){
    outfile->Printf("%s  %5d %13.6f %3s",printed++ % print_nrows == 0 ? "\n" : "",
                                       mo++,
                                       it->first,
                                       it->second.c_str());
  }
  outfile->Printf("\n  =========================================================================\n");
}

}} // End namespace
