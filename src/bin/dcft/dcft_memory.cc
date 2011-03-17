#include "dcft.h"
#include "defines.h"
#include <vector>
#include <liboptions/liboptions.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libmints/wavefunction.h>
#include <libmints/molecule.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.hpp>

using namespace boost;

namespace psi{ namespace dcft{

/**
 * Sets aside memory for the various arrays and matrices needed, and obtains
 * important information from the checkpoint file.  Must be run after an SCF
 * computation for this reason.
 */
void
DCFTSolver::init_moinfo()
{
}


/**
 * Reads the orbital information that can be determined before the SCF procedure
 * and initializes SO matrices.
 */
void
DCFTSolver::init()
{
    nso_       = reference_wavefunction_->nso();
    nirrep_    = reference_wavefunction_->nirrep();
    nmo_       = reference_wavefunction_->nmo();
    _eNuc      = Process::environment.molecule()->nuclear_repulsion_energy();
    _scfEnergy = reference_wavefunction_->reference_energy();
    _nTriSo    = nso_ * (nso_ + 1) / 2;
    for(int h = 0; h < nirrep_; ++h){
        soccpi_[h] = reference_wavefunction_->soccpi()[h];
        doccpi_[h] = reference_wavefunction_->doccpi()[h];
        frzcpi_[h] = reference_wavefunction_->frzcpi()[h];
        frzvpi_[h] = reference_wavefunction_->frzvpi()[h];
        nmopi_[h]  = reference_wavefunction_->nmopi()[h];
        nsopi_[h]  = reference_wavefunction_->nsopi()[h];
    }

    _nAOccPI = new int[nirrep_];
    _nBOccPI = new int[nirrep_];
    _nAVirPI = new int[nirrep_];
    _nBVirPI = new int[nirrep_];
    nalpha_  = nbeta_ =  _nAVir = _nBVir = 0;
    for(int h = 0; h < nirrep_; ++h){
        _nAOccPI[h] = doccpi_[h] + soccpi_[h];
        _nBOccPI[h] = doccpi_[h];
        _nAVirPI[h] = nmopi_[h] - doccpi_[h] - soccpi_[h] - frzvpi_[h];
        _nBVirPI[h] = nmopi_[h] - doccpi_[h] - frzvpi_[h];
        for(int n = 0; n < _nAOccPI[h]; ++n) ++nalpha_;
        for(int n = 0; n < _nBOccPI[h]; ++n) ++nbeta_;
        for(int n = 0; n < _nAVirPI[h]; ++n) ++_nAVir;
        for(int n = 0; n < _nBVirPI[h]; ++n) ++_nBVir;
    }

    _aOccC     = shared_ptr<Matrix>(new Matrix("Alpha Occupied MO Coefficients", nirrep_, nsopi_, _nAOccPI));
    _bOccC     = shared_ptr<Matrix>(new Matrix("Beta Occupied MO Coefficients", nirrep_, nsopi_, _nBOccPI));
    _aVirC     = shared_ptr<Matrix>(new Matrix("Alpha Virtual MO Coefficients", nirrep_, nsopi_, _nAVirPI));
    _bVirC     = shared_ptr<Matrix>(new Matrix("Beta Virtual MO Coefficients", nirrep_, nsopi_, _nBVirPI));
    _aScfError = shared_ptr<Matrix>(new Matrix("Alpha SCF Error Vector", nirrep_, _nAOccPI, _nAVirPI));
    _bScfError = shared_ptr<Matrix>(new Matrix("Beta SCF Error Vector", nirrep_, _nBOccPI, _nBVirPI));
    _Fa        = shared_ptr<Matrix>(new Matrix("Alpha Fock Matrix", nirrep_, nsopi_, nsopi_));
    _Fb        = shared_ptr<Matrix>(new Matrix("Beta Fock Matrix", nirrep_, nsopi_, nsopi_));
    Ca_        = shared_ptr<Matrix>(new Matrix("Alpha MO Coefficients", nirrep_, nsopi_, nsopi_));
    Cb_        = shared_ptr<Matrix>(new Matrix("Beta MO Coefficients", nirrep_, nsopi_, nsopi_));
    _oldCa     = shared_ptr<Matrix>(new Matrix("Old Alpha MO Coefficients", nirrep_, nsopi_, nsopi_));
    _oldCb     = shared_ptr<Matrix>(new Matrix("Old Beta MO Coefficients", nirrep_, nsopi_, nsopi_));
    _aKappa    = shared_ptr<Matrix>(new Matrix("Alpha Kappa Matrix", nirrep_, nsopi_, nsopi_));
    _bKappa    = shared_ptr<Matrix>(new Matrix("Beta Kappa Matrix", nirrep_, nsopi_, nsopi_));
    _aGTau     = shared_ptr<Matrix>(new Matrix("Alpha External Potential Matrix", nirrep_, nsopi_, nsopi_));
    _bGTau     = shared_ptr<Matrix>(new Matrix("Beta External Potential Matrix", nirrep_, nsopi_, nsopi_));
    _aoS       = shared_ptr<Matrix>(new Matrix("SO Basis Overlap Integrals", nirrep_, nsopi_, nsopi_));
    _soH       = shared_ptr<Matrix>(new Matrix("SO basis one-electron integrals", nirrep_, nsopi_, nsopi_));
    _sHalfInv  = shared_ptr<Matrix>(new Matrix("SO Basis Inverse Square Root Overlap Matrix", nirrep_, nsopi_, nsopi_));
    epsilon_a_    = shared_ptr<Vector>(new Vector(nirrep_, nsopi_));
    epsilon_b_    = shared_ptr<Vector>(new Vector(nirrep_, nsopi_));

    _aTau    = new double**[nirrep_];
    _bTau    = new double**[nirrep_];
    for(int h = 0; h < nirrep_; ++h){
        _aTau[h] = block_matrix(nsopi_[h], nsopi_[h]);
        _bTau[h] = block_matrix(nsopi_[h], nsopi_[h]);
    }
    

    // Store the AO overlap matrix
    double *sArray = new double[_nTriSo];
    IWL::read_one(psio_.get(), PSIF_OEI, PSIF_SO_S, sArray, _nTriSo, 0, 0, outfile);
    _aoS->set(sArray);
    delete [] sArray;

    // Form S^(-1/2) matrix
    Matrix eigvec(nirrep_, nsopi_, nsopi_);
    Matrix eigtemp(nirrep_, nsopi_, nsopi_);
    Matrix eigtemp2(nirrep_, nsopi_, nsopi_);
    Vector eigval(nirrep_, nsopi_);
    _aoS->diagonalize(eigvec, eigval);
    // Convert the eigenvales to 1/sqrt(eigenvalues)
    for (int h=0; h < nirrep_; ++h) {
        for (int i=0; i < nsopi_[h]; ++i) {
            double scale = 1.0 / sqrt(eigval.get(h, i));
            eigval.set(h, i, scale);
        }
    }
    eigtemp2.set(eigval);
    eigtemp.gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    _sHalfInv->gemm(false, false, 1.0, eigvec, eigtemp, 0.0);
}


/**
 * Frees up the memory sequestered by the init_moinfo() and read_checkpoint() routines.
 */
void
DCFTSolver::free_moinfo()
{
    for(int h = 0; h < nirrep_; ++h){
        //free_block(_aOccC[h]);
        //free_block(_bOccC[h]);
        //free_block(_aVirC[h]);
        //free_block(_bVirC[h]);
        free_block(_aTau[h]);
        free_block(_bTau[h]);
    }
    //delete [] _aOccC;
    //delete [] _bOccC;
    //delete [] _aVirC;
    //delete [] _bVirC;
    delete [] _aTau;
    delete [] _bTau;
    delete [] _nAOccPI;
    delete [] _nBOccPI;
    delete [] _nAVirPI;
    delete [] _nBVirPI;
}

}}//Namespaces
