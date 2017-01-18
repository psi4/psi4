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

#include "psi4-dec.h"
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libplugin/plugin.h>
#include "psifiles.h"
#include "mp2.h"

INIT_PLUGIN

namespace psi{ namespace mollerplesset2{

extern "C" int
read_options(std::string name, Options &options){
    if(name == "MOLLERPLESSET2") {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- How to cache quantities within the DPD library -*/
        options.add_int("CACHELEV", 2);
        /*- The amount of memory available (in Mb) -*/
        options.add_int("MEMORY", 2000);
        /*- The Reference -*/
        options.add_str("REFERENCE", "");
    }
}


extern "C"
SharedWavefunction mollerplesset2(SharedWavefunction ref_wfn, Options &options)
{
    int print = options.get_int("PRINT");
    // This will print out all of the user-provided options for this module
    if (print > 2) options.print();
   
    psio = _default_psio_lib_;

    eSCF = ref_wfn->reference_energy();

    nirreps             = ref_wfn->nirrep();
    SharedVector aEvals = ref_wfn->epsilon_a();
    SharedVector bEvals = ref_wfn->epsilon_b();
    char **labels       = ref_wfn->molecule()->irrep_labels();
    aOccOrbsPI          = new int[nirreps];
    bOccOrbsPI          = new int[nirreps];
    aVirOrbsPI          = new int[nirreps];
    bVirOrbsPI          = new int[nirreps];
    mopi                = ref_wfn->nmopi();
    nmo                 = ref_wfn->nmo();
    clsdpi              = ref_wfn->doccpi();
    openpi              = ref_wfn->soccpi();
    frzcpi              = ref_wfn->frzcpi();
    frzvpi              = ref_wfn->frzvpi();

    numAOcc = 0; numBOcc = 0; numAVir = 0; numBVir = 0;
    int aOccCount = 0, bOccCount = 0, aVirCount = 0, bVirCount = 0;
    for(int h = 0; h < nirreps; ++h){
        aOccOrbsPI[h] = clsdpi[h] + openpi[h] - frzcpi[h];
        bOccOrbsPI[h] = clsdpi[h] - frzcpi[h];
        aVirOrbsPI[h] = mopi[h]   - clsdpi[h] - openpi[h] - frzvpi[h];
        bVirOrbsPI[h] = mopi[h]   - clsdpi[h] - frzvpi[h];
        numAOcc += aOccOrbsPI[h];
        numBOcc += bOccOrbsPI[h];
        numAVir += aVirOrbsPI[h];
        numBVir += bVirOrbsPI[h];
    }
    aOccEvals = new double[numAOcc];
    bOccEvals = new double[numBOcc];
    aVirEvals = new double[numAVir];
    bVirEvals = new double[numBVir];

    outfile->Printf("\n\n\tIrrep  Core  Docc  Socc  aOcc  aVir  bOcc  bVir\n");
    outfile->Printf("\t===============================================\n");
    for(int h = 0; h < nirreps; ++h){
       outfile->Printf("\t %3s   %3d   %3d   %3d   %3d   %3d   %3d   %3d\n",
                             labels[h], frzcpi[h], clsdpi[h], openpi[h], 
                             aOccOrbsPI[h], aVirOrbsPI[h], bOccOrbsPI[h], bVirOrbsPI[h]);
    }
    outfile->Printf("\t===============================================\n\n");

    aOccCount = 0; bOccCount = 0; aVirCount = 0; bVirCount = 0;
    for(int h = 0; h < nirreps; ++h){
        for(int a = frzcpi[h]; a < clsdpi[h] + openpi[h]; ++a) aOccEvals[aOccCount++] = aEvals->get(h, a);
        for(int b = frzcpi[h]; b < clsdpi[h]; ++b) bOccEvals[bOccCount++] = bEvals->get(h, b);
        for(int a = clsdpi[h] + openpi[h]; a < mopi[h]; ++a) aVirEvals[aVirCount++] = aEvals->get(h, a);
        for(int b = clsdpi[h]; b < mopi[h]; ++b) bVirEvals[bVirCount++] = bEvals->get(h, b);
    }

    if(print > 2){
        for(int i = 0; i < numAOcc; ++i)
            outfile->Printf("\taOccEvals[%2d] = %10.6f\n", i, aOccEvals[i]);
        outfile->Printf("\n");
        for(int i = 0; i < numBOcc; ++i)
            outfile->Printf("\tbOccEvals[%2d] = %10.6f\n", i, bOccEvals[i]);
        outfile->Printf("\n");
        for(int i = 0; i < numAVir; ++i)
            outfile->Printf("\taVirEvals[%2d] = %10.6f\n", i, aVirEvals[i]);
        outfile->Printf("\n");
        for(int i = 0; i < numBVir; ++i)
            outfile->Printf("\tbVirEvals[%2d] = %10.6f\n", i, bVirEvals[i]);
        outfile->Printf("\n");
    }

    double eMP2;
    if(options.get_str("REFERENCE") == "UHF"){
        eMP2 = plugin_mp2_unrestricted(ref_wfn, options);
    }else if(options.get_str("REFERENCE") == "ROHF"){
        eMP2 = plugin_mp2_unrestricted(ref_wfn, options);
    }else if(options.get_str("REFERENCE") == "RHF"){
        eMP2 = plugin_mp2_restricted(ref_wfn, options);
    }else{
         std::string str1 = options.get_str("REFERENCE");
         std::string str2 = "reference in mollerplesset2";
         throw FeatureNotImplemented(str1, str2, __FILE__, __LINE__);
    }

    Process::environment.globals["MP2 ENERGY"] = eMP2 + eSCF; 
    Process::environment.globals["CURRENT ENERGY"] = eMP2 + eSCF; 
    outfile->Printf("\n\t\tSCF Reference energy  = %20.16f\n", eSCF);
    outfile->Printf("\t\tMP2 Total Energy      = %20.16f\n\n", eMP2 + eSCF);

    return ref_wfn;
}

}} // End Namespaces
