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
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "psi4/psi4-dec.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "moinfo_base.h"



using namespace std;

namespace psi {

MOInfoBase::MOInfoBase(Wavefunction& ref_wfn_, Options& options_, bool silent_)
: options(options_), silent(silent_), ref_wfn(ref_wfn_)
{
    startup();
    charge       = ref_wfn.molecule()->molecular_charge();
    multiplicity = ref_wfn.molecule()->multiplicity();
}

MOInfoBase::~MOInfoBase()
{
  cleanup();
}

void MOInfoBase::startup()
{
  nso = 0;
  nmo = 0;
  ndocc = 0;
  nactv = 0;
  nael = 0;
  nbel = 0;
  nactive_ael = 0;
  nactive_bel = 0;
  wfn_sym = 0;

  guess_occupation = true;
  PSI_NULL(ioff);
  compute_ioff();
}

void MOInfoBase::cleanup()
{
    PSI_DELETE_ARRAY(ioff);
}

void MOInfoBase::read_data()
{

    nirreps        = ref_wfn.nirrep();
    nso            = ref_wfn.nso();
    // Read sopi and save as a STL vector
    sopi           = convert_int_array_to_vector(nirreps, ref_wfn.nsopi());
    irr_labs       = ref_wfn.molecule()->irrep_labels();
    nuclear_energy = ref_wfn.molecule()->nuclear_repulsion_energy();
}

void MOInfoBase::compute_number_of_electrons()
{
    int nel   = 0;
    int natom = ref_wfn.molecule()->natom();

    for(int i=0; i < natom;i++){
        nel += static_cast<int>(ref_wfn.molecule()->Z(i));
    }
    nel -= charge;

    // Check if the multiplicity makes sense
    if( ((nel+1-multiplicity) % 2) != 0)
        throw PSIEXCEPTION("\n\n  MOInfoBase: Wrong multiplicity.\n\n");
    nael = (nel + multiplicity -1)/2;
    nbel =  nel - nael;
}

void MOInfoBase::compute_ioff()
{
  ioff    = new size_t[IOFF];
  ioff[0] = 0;
  for(size_t i=1; i < IOFF; i++)
    ioff[i] = ioff[i-1] + i;
}

void MOInfoBase::read_mo_space(int nirreps_ref, int& n, intvec& mo, string labels)
{
    bool read = false;

    vector<string> label_vec = split(labels);
    for(unsigned int k = 0; k < label_vec.size(); ++k){
        // Does the array exist in the input?
        std::string &label = label_vec[k];
        if(!options[label].has_changed()) continue; // The user didn't specify this, it's just the default
        int size  = options[label].size();
        // Defaults is to set all to zero
        mo.assign(nirreps_ref,0);
        n = 0;
        if(read){
            outfile->Printf("\n\n  libmoinfo has found a redundancy in the input keywords %s , please fix it!",labels.c_str());

            exit(1);
        }else{
            read = true;
        }
        if(size==nirreps_ref){
            for(int i=0;i<size;i++){
                mo[i] = options[label][i].to_integer();
                n += mo[i];
            }
        }else{
            outfile->Printf("\n\n  The size of the %s array (%d) does not match the number of irreps (%d), please fix the input file",label_vec[k].c_str(),size,nirreps_ref);

            exit(1);
        }
    }
}

void MOInfoBase::print_mo_space(int& n, intvec& mo, std::string labels)
{
  outfile->Printf("\n  %s",labels.c_str());

  for(int i=nirreps;i<8;i++)
    outfile->Printf("     ");
  for(int i=0;i<nirreps;i++)
    outfile->Printf(" %3d ",mo[i]);
  outfile->Printf("  %3d",n);
}

void MOInfoBase::correlate(char *ptgrp, int irrep, int& nirreps_old, int& nirreps_new,int*& arr)
{ /* This is a hack from input!!! (ACS) */
  int  i;

  if (strcmp(ptgrp,"C1 ") == 0)
    nirreps_old = 1;
  else if (strcmp(ptgrp,"Cs ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"Ci ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"C2 ") == 0)
    nirreps_old = 2;
  else if (strcmp(ptgrp,"C2v") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"D2 ") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"C2h") == 0)
    nirreps_old = 4;
  else if (strcmp(ptgrp,"D2h") == 0)
    nirreps_old = 8;
  else {
    outfile->Printf("point group %s unknown.\n",ptgrp);
    exit(1);
  }

  arr = new int[nirreps_old];

  if (irrep == 0) { /* return identity */
    nirreps_new = nirreps_old;
    for (i=0; i<nirreps_old; ++i)
      arr[i] = i;
    return;
  }

  nirreps_new = nirreps_old / 2;
  if ((strcmp(ptgrp,"C1 ") == 0) || (strcmp(ptgrp,"Cs ") == 0) ||
      (strcmp(ptgrp,"Ci ") == 0) || (strcmp(ptgrp,"C2 ") == 0) ) {
        arr[0] = 0; arr[1] = 0;
  }
  else if ( (strcmp(ptgrp,"C2v") == 0) || (strcmp(ptgrp,"D2 ") == 0) ||
            (strcmp(ptgrp,"C2h") == 0) ) {
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
    }
  }
  else if (strcmp(ptgrp,"D2h") == 0) {
    /* 1,2,3 give C2h displaced geometries */
    if (irrep == 1) {
      arr[0] = 0;  arr[1] = 0; arr[2] = 1;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 2; arr[6] = 3;  arr[7] = 3;
    }
    else if (irrep == 2) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 0;  arr[3] = 1;
      arr[4] = 2;  arr[5] = 3; arr[6] = 2;  arr[7] = 3;
    }
    else if (irrep == 3) {
      arr[0] = 0;  arr[1] = 1; arr[2] = 1;  arr[3] = 0;
      arr[4] = 2;  arr[5] = 3; arr[6] = 3;  arr[7] = 2;
    }
    /* 4 gives D2 displaced geometries */
    else if (irrep == 4) { /* D2 */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 0;  arr[5] = 1; arr[6] = 2;  arr[7] = 3;
    }
    /* displacements along irreps 5,6,7 make C2v structures */
    /* care is taken to make sure definition of b1 and b2 will
       match those that input will generate - the following seems to work:
       b1u disp: has C2(z), b2 irrep symmetric wrt sigma(yz)
       b2u disp: has C2(y), b2 irrep symmetric wrt sigma(xy)
       b3u disp: has C2(x), b2 irrep symmetric wrt sigma(xz) */
    else if (irrep == 5) { /* b1u */
      arr[0] = 0;  arr[1] = 1; arr[2] = 2;  arr[3] = 3;
      arr[4] = 1;  arr[5] = 0; arr[6] = 3;  arr[7] = 2;
    }
    else if (irrep == 6) { /* b2u */
      arr[0] = 0;  arr[1] = 3; arr[2] = 1;  arr[3] = 2;
      arr[4] = 1;  arr[5] = 2; arr[6] = 0;  arr[7] = 3;
    }
    else if (irrep == 7) { /* b3u */
      arr[0] = 0;  arr[1] = 2; arr[2] = 3;  arr[3] = 1;
      arr[4] = 1;  arr[5] = 3; arr[6] = 2;  arr[7] = 0;
    }
  }
  else {
    outfile->Printf("Point group unknown for correlation table.\n");
  }
  return;
}

intvec MOInfoBase::convert_int_array_to_vector(int n, const int* array)
{
    // Read an integer array and save as a STL vector
    intvec stl_vector(&array[0],&array[n]);
    return(stl_vector);
}

}
