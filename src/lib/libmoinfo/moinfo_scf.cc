#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>

#include <libmints/corrtab.h>
#include <libmints/molecule.h>
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

    // Map the symmetry of the input occupations, to account for displacements
    boost::shared_ptr<PointGroup> old_pg = Process::environment.parent_symmetry();
    if(old_pg){
        // This is one of a series of displacements;  check the dimension against the parent point group
        int nirreps_ref = old_pg->char_table().nirrep();

        intvec docc_ref;
        intvec actv_ref;

        read_mo_space(nirreps_ref,ndocc,docc_ref,"DOCC");
        read_mo_space(nirreps_ref,nactv,actv_ref,"SOCC");

        // Build the correlation table between full, and subgroup
        boost::shared_ptr<PointGroup> full = Process::environment.parent_symmetry();
        boost::shared_ptr<PointGroup> sub =  Process::environment.molecule()->point_group();
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
