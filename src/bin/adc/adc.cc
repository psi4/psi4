#include <psi4-dec.h>
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libtrans/integraltransform.h>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include "adc.h"

namespace psi{ namespace adc {

ADC::ADC(): Wavefunction(Process::environment.options, _default_psio_lib_)
{   
    
    copy(Process::environment.reference_wavefunction()); 
    
    char **irreps_      = Process::environment.molecule()->irrep_labels();
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
    
    fprintf(outfile, "\n\n\tIrrep  Core  Docc  Socc  aOcc  aVir  bOcc  bVir  FVir\n");
    fprintf(outfile,     "\t*****************************************************\n");
    for(int h = 0; h < nirrep_; h++){
        fprintf(outfile, "\t %3s   %3d   %3d   %3d   %3d   %3d   %3d   %3d   %3d\n",
                irreps_[h], frzcpi_[h], doccpi_[h], soccpi_[h], 
                aoccpi_[h], avirpi_[h], boccpi_[h], bvirpi_[h], frzvpi_[h]);
    }
    fprintf(outfile,     "\t*****************************************************\n\n");
    fflush(outfile);
  
    conv_     = options_.get_double("NEWTON_CONV");
    norm_tol_ = options_.get_double("NORM_TOLERANCE");
    pole_max_ = options_.get_int("POLE_MAX");
    sem_max_  = options_.get_int("SEM_MAX");
    num_amps_ = options_.get_int("NUM_AMPS");

  if(options_["STATES_PER_IRREP"].size() > 0){
        int i = options_["STATES_PER_IRREP"].size();
        
        if(i != nirrep_){
            fprintf(outfile, "dim of states_per_irrep vector must be %d\n", nirrep_);
            throw PsiException("adc input comparison error STATES_PER_IRREP and nirrep_", __FILE__, __LINE__);
        }
        rpi_ = options_.get_int_array("STATES_PER_IRREP");
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

    fprintf(outfile, "\t==> Input Parameters <==\n");
    fprintf(outfile, "\tNEWTON_CONV = %3g, NORM_TOL = %3g\n", conv_, norm_tol_);
    fprintf(outfile, "\tPOLE_MAX    = %3d, SEM_MAX  = %3d\n\n", pole_max_, sem_max_);
    
    fprintf(outfile, "\tNXS           = %d\n", nxs_);
//    fprintf(outfile, "\tIRREP_XYZ     = [");
//    for(int i = 0;i < 3;i++)
//        fprintf(outfile, " %3s ", irreps_[irrep_axis_[i]]);
//    fprintf(outfile, "]\n");
    fprintf(outfile, "\tNXS_PER_IRREP = [");
    for(int i = 0;i < nirrep_;i++){
        fprintf(outfile, " %d ", nxspi_[i]);
    }
    fprintf(outfile, "]\n");

    if(DEBUG_){
        fprintf(outfile, "Debagging mode...\n");
        fprintf(outfile, "\tNMO   = %3d, NXS = %3d\n", nmo_, nxs_);
        fprintf(outfile, "\tNOPEN = %3d\n", nopen_);
        fprintf(outfile, "\tSTATES_PER_IRREP = [");
        for(int i = 0;i < nirrep_;i++) fprintf(outfile, "%3d", rpi_[i]);
        fprintf(outfile, " ]\n");
    }
    
}

void
ADC::release_mem()
{
    free(poles_);
    delete _ints;
    delete aocce_;
    delete avire_;
    delete bocce_;
    delete bvire_;
    
    //omega_guess_.reset();
}

ADC::~ADC()
{
}
    
}} // End Namespaces
