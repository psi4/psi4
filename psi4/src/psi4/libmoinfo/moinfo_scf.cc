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

#include "psi4/libmints/corrtab.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/liboptions/liboptions.h"

#include "moinfo_scf.h"



using namespace std;

namespace psi {

MOInfoSCF::MOInfoSCF(Wavefunction& ref_wfn_, Options& options_, bool silent_)
    : MOInfoBase(ref_wfn_, options_, silent_)
{
    read_data();
    // Determine the wave function irrep
    // The first irrep is 0
    bool wfn_sym_found = false;
    wfn_sym = 0;
    string wavefunction_sym_str = options.get_str("WFN_SYM");
    for(int h = 0; h < nirreps; ++h){
        string irr_label_str = irr_labs[h];
        to_upper(irr_label_str);
        trim_spaces(irr_label_str);
        if(wavefunction_sym_str == irr_label_str){
            wfn_sym = h;
            wfn_sym_found = true;
            break;
        }
        if(wavefunction_sym_str == to_string(h+1)){
            wfn_sym = h;
            wfn_sym_found = true;
            break;
        }
    }
    if(!wfn_sym_found)
        throw PSIEXCEPTION("Wavefuntion symmetry " + wavefunction_sym_str +
                           " is not a valid choice for this point group");

    compute_number_of_electrons();
    read_mo_spaces();
    print_mo();
}

MOInfoSCF::~MOInfoSCF()
{
}

void MOInfoSCF::read_mo_spaces()
{
  /*****************************************************
     See if we're in a subgroup for finite difference
     calculations, by looking to see what OptKing has
     written to the checkpoint file.  Reassign the
     occupation vectors as appropriate.  N.B. the
     SOCC and DOCC are handled by Input (ACS)
  *****************************************************/

    docc.resize(nirreps,0);
    actv.resize(nirreps,0);

    // Map the symmetry of the input occupations, to account for displacements
    std::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
    if(old_pg){
        // This is one of a series of displacements;  check the dimension against the parent point group
        int nirreps_ref = old_pg->char_table().nirrep();

        intvec docc_ref;
        intvec actv_ref;

        read_mo_space(nirreps_ref,ndocc,docc_ref,"DOCC");
        read_mo_space(nirreps_ref,nactv,actv_ref,"SOCC");

        // Build the correlation table between full, and subgroup
        std::shared_ptr<PointGroup> full = Process::environment.parent_symmetry();
        std::shared_ptr<PointGroup> sub =  ref_wfn.molecule()->point_group();
        CorrelationTable corrtab(full, sub);

        // Find the occupation in the subgroup
        for(int h = 0; h < nirreps_ref; ++h){
            int target = corrtab.gamma(h, 0);
            docc[target] += docc_ref[h];
            actv[target] += actv_ref[h];
        }
    }else{
        // For a single-point only
        read_mo_space(nirreps,ndocc,docc,"DOCC");
        read_mo_space(nirreps,nactv,actv,"SOCC");
//        read_mo_space(nirreps,nactv,actv,"ACTV ACTIVE SOCC");
    }

    nactive_ael = nael  - ndocc;
    nactive_bel = nbel  - ndocc;

    if((ndocc > 0) || (nactv > 0))
        guess_occupation = false;
}

void MOInfoSCF::print_mo()
{
    outfile->Printf("\n");
    outfile->Printf("\n  MOs per irrep:                ");

    for(int i=nirreps;i<8;i++)
        outfile->Printf("     ");
    for(int i=0;i<nirreps;i++)
        outfile->Printf("  %s",irr_labs[i]);
    outfile->Printf(" Total");
    outfile->Printf("\n  ----------------------------------------------------------------------------");
    print_mo_space(nso,sopi,"Total                         ");
    if(!guess_occupation){
        print_mo_space(ndocc,docc,"Doubly Occupied               ");
        print_mo_space(nactv,actv,"Active/Singly Occupied        ");
    }
    outfile->Printf("\n  ----------------------------------------------------------------------------");
    if(guess_occupation)
        outfile->Printf("\n\n  Guessing orbital occupation");

}

}
