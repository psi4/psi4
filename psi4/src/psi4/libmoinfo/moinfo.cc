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

// Standard Libraries
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

// STL
#include <numeric>

// PSI Libraries
#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/corrtab.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

#include "moinfo.h"



using namespace std;

namespace psi {

MOInfo::MOInfo(Wavefunction& ref_wfn_, Options& options_, bool silent_)
    : MOInfoBase(ref_wfn_, options_, silent_)
{
    /***************
    Set defaults
  ***************/
    //  OrbitalSpace focc_os("focc","f");
    //  OrbitalSpace docc_os("docc","d");
    //  OrbitalSpace actv_os("actv","a");
    //  OrbitalSpace extr_os("extr","e");
    //  OrbitalSpace fvir_os("fvir","x");
    //
    //  OrbitalSpace occ_os("occ","o");
    //  occ_os.add_subspace(docc);
    //  occ_os.add_subspace(actv);
    //
    //  OrbitalSpace vir_os("vir","v");
    //  vir_os.add_subspace(actv);
    //  vir_os.add_subspace(extr);
    //
    //  mo_spaces.add_subspace(focc);
    //  mo_spaces.add_subspace(occ);
    //  mo_spaces.add_subspace(vir);
    //  mo_spaces.add_subspace(fvir);
    //
    //  mo_spaces.print();


    if(options_["CORR_MULTP"].has_changed())
        multiplicity = options_.get_int("CORR_MULTP");

    no_damp_convergence = 1.0e-9;
    dgemm_timing        = 0.0;
    scf                 = NULL;

    nfocc = 0;
    nfvir = 0;
    nactv_docc = 0;
    nocc  = 0;
    nvir  = 0;
    nall  = 0;
    nextr = 0;
    root  = 0;

    read_info();
    read_mo_spaces();
    compute_mo_mappings();

    if(!silent){
        print_info();
        print_mo();
    }
}

MOInfo::~MOInfo()
{
    free_memory();
}

void MOInfo::read_info()
{
    /*
     * Read Nuclear, SCF and other stuff
     */
    read_data();
    nmo            = ref_wfn.nmo();
    compute_number_of_electrons();
    scf_energy     = ref_wfn.reference_energy();
    mopi           = convert_int_array_to_vector(nirreps, ref_wfn.nmopi());
    SharedMatrix matCa = ref_wfn.Ca();
    scf            = block_matrix(nso, nmo);
    unsigned int soOffset = 0;
    unsigned int moOffset = 0;
    for(int h = 0; h < nirreps; ++h){
        for(int so = 0; so < sopi[h]; ++so){
            for(int mo = 0; mo < mopi[h]; ++mo){
                scf[so+soOffset][mo+moOffset] = matCa->get(h, so, mo);
            }
        }
        soOffset += sopi[h];
        moOffset += mopi[h];
    }
    scf_irrep      = new double**[nirreps];
    for(int i=0;i<nirreps;i++){
        scf_irrep[i] = 0;
        if(sopi[i] && mopi[i]){
            scf_irrep[i] = block_matrix(sopi[i], mopi[i]);
            ::memcpy(scf_irrep[i][0], matCa->pointer(i)[0], mopi[i]*sopi[i]*sizeof(double));
        }
    }

    // Determine the wave function irrep
    // The defalut irrep is 0 (A)
    wfn_sym = 0;
    string wavefunction_sym_str = options.get_str("WFN_SYM");
    bool wfn_sym_found = false;

    std::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
    if(old_pg){
        for(int h = 0; h < nirreps; ++h){
            string irr_label_str = old_pg->char_table().gamma(h).symbol_ns();
            trim_spaces(irr_label_str);
            to_upper(irr_label_str);
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
    }else{
        for(int h = 0; h < nirreps; ++h){
            string irr_label_str = irr_labs[h];
            trim_spaces(irr_label_str);
            to_upper(irr_label_str);
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
    }
    if(!wfn_sym_found)
        throw PSIEXCEPTION("Wavefuntion symmetry " + wavefunction_sym_str +
                           " is not a valid choice for this point group");
    // The lowest root in the input is 1, here we subtract one
    root = options.get_int("FOLLOW_ROOT") - 1;
}

void MOInfo::setup_model_space()
{
    // NB This code could be places elsewhere
    references.clear();
    alpha_internal_excitations.clear();
    beta_internal_excitations.clear();
    sign_internal_excitations.clear();
    all_refs.clear();
    unique_refs.clear();
    closed_shell_refs.clear();
    unique_open_shell_refs.clear();

    build_model_space();
    print_model_space();
    make_internal_excitations();
}





/*!
    \fn MOInfo::print_info()
 */
void MOInfo::print_info()
{
    outfile->Printf("\n");
    outfile->Printf("\n  ==============================================================================");
    outfile->Printf("\n  System Info:");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
    outfile->Printf("\n  Nuclear Energy   = %-15.9f  SCF Energy       = %-15.9f",nuclear_energy,scf_energy);
    outfile->Printf("\n");
    outfile->Printf("\n  MOs and Symmetry:");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
    outfile->Printf("\n  nirreps          = %-10d       root             = %-10d",nirreps,root);
    outfile->Printf("\n  nso              = %-10d       nmo              = %-10d",nso,nmo);
    outfile->Printf("\n  nael             = %-10d       nbel             = %-10d",nael,nbel);
    outfile->Printf("\n  nactive_ael      = %-10d       nactive_bel      = %-10d",nactive_ael,nactive_bel);
    outfile->Printf("\n");
    outfile->Printf("\n  Details of the Computation:");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
}

/*!
    \fn MOInfo::read_mo_spaces()
 */
void MOInfo::read_mo_spaces()
{
    /*****************************************************
     See if we're in a subgroup for finite difference
     calculations, by looking to see what OptKing has
     written to the checkpoint file.  Reassign the
     occupation vectors as appropriate.  N.B. the
     SOCC and DOCC are handled by Input (ACS)
  *****************************************************/

    focc.assign(nirreps,0);
    docc.assign(nirreps,0);
    actv.assign(nirreps,0);
    fvir.assign(nirreps,0);
    extr.assign(nirreps,0);
    occ.assign(nirreps,0);
    vir.assign(nirreps,0);
    all.assign(nirreps,0);
    actv_docc.assign(nirreps,0);

    // Map the symmetry of the input occupations, to account for displacements
    std::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
    if(old_pg){
        // This is one of a series of displacements;  check the dimension against the parent point group
        int nirreps_ref = old_pg->char_table().nirrep();
        intvec focc_ref;
        intvec docc_ref;
        intvec actv_ref;
        intvec fvir_ref;
//        intvec actv_docc_ref;

        focc_ref = convert_int_array_to_vector(nirreps, ref_wfn.frzcpi());
        docc_ref = convert_int_array_to_vector(nirreps, ref_wfn.doccpi());
        actv_ref = convert_int_array_to_vector(nirreps, ref_wfn.soccpi());
        fvir_ref.assign(nirreps_ref,0);
//        actv_docc_ref.assign(nirreps_ref,0);

        for (int h = 0; h < nirreps_ref; h++)
            docc_ref[h] -= focc_ref[h];

        nfocc = std::accumulate( focc_ref.begin(), focc_ref.end(), 0 );
        ndocc = std::accumulate( docc_ref.begin(), docc_ref.end(), 0 );
        nactv = std::accumulate( actv_ref.begin(), actv_ref.end(), 0 );

        read_mo_space(nirreps_ref,nfocc,focc_ref,"FROZEN_DOCC");
        read_mo_space(nirreps_ref,ndocc,docc_ref,"RESTRICTED_DOCC");
        read_mo_space(nirreps_ref,nactv,actv_ref,"ACTIVE");
        read_mo_space(nirreps_ref,nfvir,fvir_ref,"FROZEN_UOCC");



        std::shared_ptr<PointGroup> full = Process::environment.parent_symmetry();
        std::shared_ptr<PointGroup> sub =  ref_wfn.molecule()->point_group();
        // Build the correlation table between full, and subgroup
        CorrelationTable corrtab(full, sub);

        // Find the wave function symmetry
        wfn_sym = corrtab.gamma(wfn_sym, 0);

        // Find the occupation in the subgroup
        for(int h = 0; h < nirreps_ref; ++h){
            int target = corrtab.gamma(h, 0);
            focc[target] += focc_ref[h];
            docc[target] += docc_ref[h];
            actv[target] += actv_ref[h];
            fvir[target] += fvir_ref[h];
        }
        //    read_mo_space(nirreps_ref,nactv_docc,actv_docc_ref,"ACTIVE_DOCC");
    }else{
        // For a single-point only
        outfile->Printf("\n  For a single-point only");

        focc = convert_int_array_to_vector(nirreps, ref_wfn.frzcpi());
        docc = convert_int_array_to_vector(nirreps, ref_wfn.doccpi());
        actv = convert_int_array_to_vector(nirreps, ref_wfn.soccpi());

        for (int h = 0; h < nirreps; h++)
            docc[h] -= focc[h];

        nfocc = std::accumulate( focc.begin(), focc.end(), 0 );
        ndocc = std::accumulate( docc.begin(), docc.end(), 0 );
        nactv = std::accumulate( actv.begin(), actv.end(), 0 );

        read_mo_space(nirreps,nfocc,focc,"FROZEN_DOCC");
        read_mo_space(nirreps,ndocc,docc,"RESTRICTED_DOCC");
        read_mo_space(nirreps,nactv,actv,"ACTIVE");
        read_mo_space(nirreps,nfvir,fvir,"FROZEN_UOCC");
//        read_mo_space(nirreps,nfocc,focc,"CORR_FOCC FROZEN_DOCC");
//        read_mo_space(nirreps,ndocc,docc,"CORR_DOCC RESTRICTED_DOCC");
//        read_mo_space(nirreps,nactv,actv,"CORR_ACTV ACTV ACTIVE");
//        read_mo_space(nirreps,nfvir,fvir,"CORR_FVIR FROZEN_UOCC");
        //    read_mo_space(nirreps,nactv_docc,actv_docc,"ACTIVE_DOCC");
    }

    // Compute the number of external orbitals per irrep
    nextr = 0;
    for(int h = 0; h < nirreps; ++h){
        extr[h]= mopi[h] - focc[h] - docc[h] - actv[h] - fvir[h];
        occ[h] = docc[h] + actv[h];
        vir[h] = actv[h] + extr[h];
        all[h] = mopi[h] - focc[h] - fvir[h];
        nextr += extr[h];
    }
    nall        = nmo  - nfocc - nfvir;
    nactive_ael = nael - ndocc - nfocc;
    nactive_bel = nbel - ndocc - nfocc;
    nocc        = ndocc + nactv;
    nvir        = nactv + nextr;

    bool active_space_problem = false;
    string error_msg;
    if(nactv < nactive_ael){
        error_msg += "\n  - the number of active orbitals (nactv = " + to_string(nactv) + ")";
        error_msg += " is smaller than the number of active alpha electrons (nactive_ael =" + to_string(nactive_ael) +")",
                active_space_problem = true;
    }
    if(nactv < nactive_bel){
        error_msg += "\n  - the number of active orbitals (nactv = " + to_string(nactv) + ")";
        error_msg += " is smaller than the number of active beta electrons (nactive_bel =" + to_string(nactive_bel) +")",
                active_space_problem = true;
    }
    if(active_space_problem){
        error_msg = "MOInfo found a problem with the definition of the active space:" + error_msg;
        throw PSIEXCEPTION(error_msg);
    }


    /*********************************************
    Define the symmetry of each  non-frozen MO
  **********************************************/
    all_sym.resize(nall);
    int index_mo  = 0;
    for(int h = 0; h < nirreps; ++h){
        for(int i  = 0; i < all[h]; ++i){
            all_sym[index_mo] = h;
            index_mo++;
        }
    }

    /***************************************************************
    Build the array that connects the non-frozen MOs (all) to the
    the complete list of MOs (mo). Used when frozen MOs are used.
  ****************************************************************/
    all_to_mo.resize(nall);
    int index_all = 0;
    index_mo  = 0;
    for(int h = 0; h < nirreps; ++h){
        index_mo += focc[h];
        for(int i = 0; i < all[h]; ++i){
            all_to_mo[index_all]=index_mo;
            index_all++;
            index_mo++;
        }
        index_mo += fvir[h];
    }

    // The mapping of the MOs (mo) to the non-frozen MOs (all)
    // Set size to nmo and all elements to -1
    mo_to_all.assign(nmo,-1);
    // Set the mappings
    for(int i = 0; i < nall; ++i)
        mo_to_all[all_to_mo[i]]=i;

    //  /***************************************************************
    //    Build the array that connects the generalized occupied MOs
    //    to the non-frozen MOs (all).
    //  ****************************************************************/
    //  occ_to_all.resize(nocc);
    //  int index_occ = 0;
    //  index_all  = 0;
    //  for(int h = 0; h < nirreps; ++h){
    //    for(int i = 0; i < occ[h]; ++i){
    //      occ_to_all[index_occ] = index_all;
    //      index_all++;
    //      index_occ++;
    //    }
    //    index_all += extr[h];
    //  }
}

/**
    MOInfo::print_mo_spaces()
 */
void MOInfo::print_mo()
{
    /// @todo implement me
    outfile->Printf("\n");
    outfile->Printf("\n  MOs per irrep:                  ");

    for(int i=nirreps;i<8;i++)
        outfile->Printf("     ");
    for(int i=0;i<nirreps;i++)
        outfile->Printf("  %s",irr_labs[i]);
    outfile->Printf(" Total");
    outfile->Printf("\n  ------------------------------------------------------------------------------");
    print_mo_space(nmo,mopi,"Total                           ");
    print_mo_space(nfocc,focc,"Frozen Occupied                 ");
    print_mo_space(ndocc,docc,"Doubly Occupied                 ");
    print_mo_space(nactv,actv,"Active                          ");
    if(nactv_docc > 0){
        print_mo_space(nactv_docc,actv_docc,"Active Doubly Occupied          ");
    }
    print_mo_space(nextr,extr,"External                        ");
    print_mo_space(nfvir,fvir,"Frozen Virtual                  ");

}

/**
 *   MOInfo::free_memory()
 */
void MOInfo::free_memory()
{
    if(scf != NULL)
        free_block(scf);
    for(int i=0;i<nirreps;i++)
        free_block(scf_irrep[i]);
    delete[] scf_irrep;
}

}
