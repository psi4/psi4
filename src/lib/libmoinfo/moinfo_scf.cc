#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <libciomr/libciomr.h>
#include <liboptions/liboptions.h>

#include "moinfo_scf.h"

extern FILE *outfile;

using namespace std;

namespace psi {

MOInfoSCF::MOInfoSCF(Options& options_, bool silent_)
    : MOInfoBase(options_, silent_)
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

    // For single-point geometry optimizations and frequencies
    // TODO restore this functionality, avoiding chkpt calls.
    const char *current_displacement_label = "";
    if(0/*chkpt->exist(current_displacement_label)*/){
        //    int disp_irrep = 0;
        // int   disp_irrep  = chkpt->rd_disp_irrep();
        //    char *save_prefix = chkpt->rd_prefix();
        //    int nirreps_ref;

        // read symmetry info and MOs for undisplaced geometry from
        // root section of checkpoint file
        //    chkpt->reset_prefix();
        //    chkpt->commit_prefix();

        //    char *ptgrp_ref;
        //char *ptgrp_ref = chkpt->rd_sym_label();

        // Lookup irrep correlation table
        //    int* correlation;
        //    correlate(ptgrp_ref, disp_irrep, nirreps_ref, nirreps,correlation);

        //    intvec docc_ref;
        //    intvec actv_ref;

        //    // build orbital information for current point group
        //    read_mo_space(nirreps_ref,ndocc,docc_ref,"DOCC");
        //    read_mo_space(nirreps_ref,nactv,actv_ref,"ACTV ACTIVE SOCC");

        //    for (int h=0; h < nirreps_ref; h++) {
        //      docc[ correlation[h] ] += docc_ref[h];
        //      actv[ correlation[h] ] += actv_ref[h];
        //    }

        //    wfn_sym = correlation[wfn_sym];
        //    chkpt->set_prefix(save_prefix);
        //    chkpt->commit_prefix();
        //    free(save_prefix);
        //    free(ptgrp_ref);
        //    delete [] correlation;
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

    //  free(current_displacement_label);
}

void MOInfoSCF::print_mo()
{
    fprintf(outfile,"\n");
    fprintf(outfile,"\n  MOs per irrep:                ");

    for(int i=nirreps;i<8;i++)
        fprintf(outfile,"     ");
    for(int i=0;i<nirreps;i++)
        fprintf(outfile,"  %s",irr_labs[i]);
    fprintf(outfile," Total");
    fprintf(outfile,"\n  ----------------------------------------------------------------------------");
    print_mo_space(nso,sopi,"Total                         ");
    if(!guess_occupation){
        print_mo_space(ndocc,docc,"Doubly Occupied               ");
        print_mo_space(nactv,actv,"Active/Singly Occupied        ");
    }
    fprintf(outfile,"\n  ----------------------------------------------------------------------------");
    if(guess_occupation)
        fprintf(outfile,"\n\n  Guessing orbital occupation");
    fflush(outfile);
}

}
