#include "psi4-dec.h"
#include <libmints/mints.h>
#include <libciomr/libciomr.h>
#include <libplugin/plugin.h>
#include "dfadc.h"

namespace psi{ namespace plugin_dfadc {
    
void
DFADC::init()
{    
    copy(Process::environment.reference_wavefunction());
    density_fitted_ = true; // Yes, exactly!

    if(!(nirrep_ == 1))
        throw PSIEXCEPTION("DF-ADC can be applicable only for C1 case for now.\n");
    
    nopen_ = 0;
    for(int h = 0;h < nirrep_;h++)
        nopen_ += soccpi_[h];
    if (nopen_)
        throw PSIEXCEPTION("Openshell calculation has not been implemented yet.");
        
    naocc_ = doccpi_[0] - frzcpi_[0] + soccpi_[0];
    navir_ = nmopi_[0]  - doccpi_[0] - frzvpi_[0];

    occe_ = new double[naocc_];
    vire_ = new double[navir_];

    fprintf(outfile, "\n\n\tFzc   Docc  Socc  aOcc  aVir  bOcc  bVir  Fzv \n");
    fprintf(outfile,     "\t*****************************************************\n");
    fprintf(outfile, "\t %3d   %3d   %3d   %3d   %3d   %3d   %3d   %3d\n",
                frzcpi_[0], doccpi_[0], soccpi_[0], 
                naocc_, navir_, doccpi_[0], nmopi_[0]-doccpi_[0]-frzvpi_[0], frzvpi_[0]);
    fprintf(outfile,     "\t*****************************************************\n\n");
    fflush(outfile);
    
    int aoccount = 0, avircount = 0;
    for(int a = frzcpi_[0]; a < doccpi_[0]; a++)                      occe_[aoccount++]  = epsilon_a_->get(0, a);
    for(int a = doccpi_[0]+soccpi_[0]; a < nmopi_[0]-frzvpi_[0]; a++) vire_[avircount++] = epsilon_a_->get(0, a);

    occCa_ = SharedMatrix(new Matrix(nsopi_[0], naocc_));
    for(int i = 0;i < naocc_;i++)
        for(int MU = 0;MU < nsopi_[0];MU++)
            occCa_->set(MU, i, Ca_->get(0, MU, i+frzcpi_[0]));
    
    virCa_ = SharedMatrix(new Matrix(nsopi_[0], navir_));
    for(int a = 0;a < navir_;a++)
        for(int MU = 0;MU < nsopi_[0];MU++)
            virCa_->set(MU, a, Ca_->get(0, MU, a+frzcpi_[0]+naocc_));

    // Setting up RI basis
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    if(options_.get_str("RI_BASIS_ADC") == ""){
        basisset_->molecule()->set_basis_all_atoms(options_.get_str("BASIS")+"-RI", "RI-BASIS");
        fprintf(outfile, " No auxiliary basis selected, defaulting to %s-RI\n\n", options_.get_str("BASIS").c_str());
    }
    ribasis_ = BasisSet::construct(parser, molecule_, "RI_BASIS_ADC");
    boost::shared_ptr<BasisSet> _zero = BasisSet::zero_ao_basis_set();
    
    // Setting up the dipole integrals 
    boost::shared_ptr<IntegralFactory> ints(new IntegralFactory(basisset_, basisset_, basisset_, basisset_));
    boost::shared_ptr<OneBodyAOInt> dipole(ints->ao_dipole());
    
    std::vector<SharedMatrix> ao_dipole_ints;
    ao_dipole_ints.push_back(SharedMatrix(new Matrix("AO_dipole_X", basisset_->nbf(), basisset_->nbf())));
    ao_dipole_ints.push_back(SharedMatrix(new Matrix("AO_dipole_Y", basisset_->nbf(), basisset_->nbf())));
    ao_dipole_ints.push_back(SharedMatrix(new Matrix("AO dipole_Z", basisset_->nbf(), basisset_->nbf())));
    dipole->compute(ao_dipole_ints);
    
    dipole_ints_.push_back(SharedMatrix(new Matrix("MO_dipole_X", naocc_, navir_)));
    dipole_ints_.push_back(SharedMatrix(new Matrix("MO_dipole_Y", naocc_, navir_)));
    dipole_ints_.push_back(SharedMatrix(new Matrix("MO dipole_Z", naocc_, navir_)));    

    for(int I = 0;I < 3;I++)
        dipole_ints_[I]->transform(occCa_, ao_dipole_ints[I], virCa_);
    
} 

void
DFADC::release_mem()
{
    free(occe_);
    free(vire_);
    free(diag_);
    free(Ecis_);
    free(Bcis_);
    //free(poles_);
    free_block(Aij_);
    free_block(Aab_);
    free_block(Kiajb_);
}

}} // End Namespaces
