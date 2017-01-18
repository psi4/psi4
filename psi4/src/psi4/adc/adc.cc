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

#include "psi4/psi4-dec.h"

#include "psi4/liboptions/liboptions.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/adc/adc.h"
#include "psi4/libmints/molecule.h"

namespace psi{ namespace adc {

ADCWfn::ADCWfn(SharedWavefunction ref_wfn, Options& options) :
    Wavefunction(options)
{

    shallow_copy(ref_wfn);
    reference_wavefunction_ = ref_wfn;

    char **irreps_      = molecule_->irrep_labels();
    aoccpi_             = new int[nirrep_];
    boccpi_             = new int[nirrep_];
    avirpi_             = new int[nirrep_];
    bvirpi_             = new int[nirrep_];

    int naocc = 0, nbocc = 0, navir = 0, nbvir = 0;
    int aoccount = 0, boccount = 0, avircount = 0, bvircount = 0;
    for(int h = 0; h < nirrep_; h++){
        aoccpi_[h] = doccpi_[h] + soccpi_[h] - frzcpi_[h];
        boccpi_[h] = doccpi_[h] - frzcpi_[h];
        avirpi_[h] = nmopi_[h]   - doccpi_[h] - soccpi_[h] - frzvpi_[h];
        bvirpi_[h] = nmopi_[h]   - doccpi_[h] - frzvpi_[h];

        naocc += aoccpi_[h];
        nbocc += boccpi_[h];
        navir += avirpi_[h];
        nbvir += bvirpi_[h];
    }
    aocce_ = new double[naocc];
    bocce_ = new double[nbocc];
    avire_ = new double[navir];
    bvire_ = new double[nbvir];

    nopen_ = 0;
    for(int h = 0;h < nirrep_;h++)
        nopen_ += soccpi_[h];
    if (nopen_)
        throw PSIEXCEPTION("Openshell calculation has not been implemented yet!");

    aoccount = 0, boccount = 0, avircount = 0, bvircount = 0;
    for(int h = 0; h < nirrep_; h++){
        for(int a = frzcpi_[h]; a < doccpi_[h]+soccpi_[h]; a++)
            aocce_[aoccount++] = epsilon_a_->get(h, a);

        for(int b = frzcpi_[h]; b < doccpi_[h]; b++)
            bocce_[boccount++] = epsilon_b_->get(h, b);

        for(int a = doccpi_[h]+soccpi_[h]; a < nmopi_[h]-frzvpi_[h]; a++)
            avire_[avircount++] = epsilon_a_->get(h, a);

        for(int b = doccpi_[h]; b < nmopi_[h]-frzvpi_[h]; b++)
            bvire_[bvircount++] = epsilon_b_->get(h, b);
    }

    outfile->Printf( "\n\n\tIrrep  Core  Docc  Socc  aOcc  aVir  bOcc  bVir  FVir\n");
    outfile->Printf(     "\t*****************************************************\n");
    for(int h = 0; h < nirrep_; h++){
        outfile->Printf( "\t %3s   %3d   %3d   %3d   %3d   %3d   %3d   %3d   %3d\n",
                irreps_[h], frzcpi_[h], doccpi_[h], soccpi_[h],
                aoccpi_[h], avirpi_[h], boccpi_[h], bvirpi_[h], frzvpi_[h]);
    }
    outfile->Printf(     "\t*****************************************************\n\n");


    conv_     = options_.get_double("NEWTON_CONVERGENCE");
    norm_tol_ = options_.get_double("NORM_TOLERANCE");
    pole_max_ = options_.get_int("POLE_MAXITER");
    sem_max_  = options_.get_int("SEM_MAXITER");
    num_amps_ = options_.get_int("NUM_AMPS_PRINT");

  if(options_["ROOTS_PER_IRREP"].size() > 0){
        int i = options_["ROOTS_PER_IRREP"].size();

        if(i != nirrep_){
            outfile->Printf( "dim of states_per_irrep vector must be %d\n", nirrep_);
            throw PsiException("adc input comparison error ROOTS_PER_IRREP and nirrep_", __FILE__, __LINE__);
        }
        rpi_ = options_.get_int_array("ROOTS_PER_IRREP");
    }
    else {
        rpi_ = new int [nirrep_];
        for(int h = 0;h < nirrep_;h++){
            rpi_[h] = 1;
        }
    }

  // Setting up dimensions for each irrep block and totoal dimension of S manifold.
    nxs_ = 0;
    nxspi_ = new int [nirrep_];
    poles_ = (struct pole**)malloc(nirrep_*sizeof(struct pole*));
    for(int h = 0;h < nirrep_;h++){
        nxspi_[h] = 0;
        poles_[h] = (struct pole*)malloc(rpi_[h]*sizeof(struct pole));
        for(int Go = 0;Go < nirrep_;Go++){
            nxspi_[h] += aoccpi_[Go] * avirpi_[Go^h];
        }
        nxs_ += nxspi_[h];
    }

    outfile->Printf( "\t==> Input Parameters <==\n");
    outfile->Printf( "\tNEWTON_CONV = %3g, NORM_TOL = %3g\n", conv_, norm_tol_);
    outfile->Printf( "\tPOLE_MAX    = %3d, SEM_MAX  = %3d\n\n", pole_max_, sem_max_);

    outfile->Printf( "\tNXS           = %d\n", nxs_);
//    outfile->Printf( "\tIRREP_XYZ     = [");
//    for(int i = 0;i < 3;i++)
//        outfile->Printf( " %3s ", irreps_[irrep_axis_[i]]);
//    outfile->Printf( "]\n");
    outfile->Printf( "\tNXS_PER_IRREP = [");
    for(int i = 0;i < nirrep_;i++){
        outfile->Printf( " %d ", nxspi_[i]);
    }
    outfile->Printf( "]\n");

    if(DEBUG_){
        outfile->Printf( "Debagging mode...\n");
        outfile->Printf( "\tNMO   = %3d, NXS = %3d\n", nmo_, nxs_);
        outfile->Printf( "\tNOPEN = %3d\n", nopen_);
        outfile->Printf( "\tROOTS_PER_IRREP = [");
        for(int i = 0;i < nirrep_;i++) outfile->Printf( "%3d", rpi_[i]);
        outfile->Printf( " ]\n");
    }

}

ADCWfn::~ADCWfn()
{
}

void ADCWfn::release_mem()
{
    free(poles_);
    delete _ints;
    delete aocce_;
    delete avire_;
    delete bocce_;
    delete bvire_;

    //omega_guess_.reset();
}

}} // End Namespaces
