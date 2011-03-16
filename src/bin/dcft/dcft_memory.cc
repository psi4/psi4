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
DCFTSolver::read_checkpoint()
{
    _nSo       = Process::environment.reference_wavefunction()->nso();
    _soPI      = Process::environment.reference_wavefunction()->nsopi();
    _nIrreps   = Process::environment.reference_wavefunction()->nirrep();
    _nMo       = Process::environment.reference_wavefunction()->nmo();
    _moPI      = Process::environment.reference_wavefunction()->nmopi();
    _openPI    = Process::environment.reference_wavefunction()->soccpi();
    _clsdPI    = Process::environment.reference_wavefunction()->doccpi();
    _frzcPI    = Process::environment.reference_wavefunction()->frzcpi();
    _frzvPI    = Process::environment.reference_wavefunction()->frzvpi();
    _eNuc      = Process::environment.molecule()->nuclear_repulsion_energy();
    _scfEnergy = Process::environment.reference_wavefunction()->reference_energy();
    _nTriSo    = _nSo * (_nSo + 1) / 2;

    _nAOccPI = new int[_nIrreps];
    _nBOccPI = new int[_nIrreps];
    _nAVirPI = new int[_nIrreps];
    _nBVirPI = new int[_nIrreps];
    nalpha_  = nbeta_ =  _nAVir = _nBVir = 0;
    for(int h = 0; h < _nIrreps; ++h){
        _nAOccPI[h] = _clsdPI[h] + _openPI[h];
        _nBOccPI[h] = _clsdPI[h];
        _nAVirPI[h] = _moPI[h] - _clsdPI[h] - _openPI[h] - _frzvPI[h];
        _nBVirPI[h] = _moPI[h] - _clsdPI[h] - _frzvPI[h];
        for(int n = 0; n < _nAOccPI[h]; ++n) ++nalpha_;
        for(int n = 0; n < _nBOccPI[h]; ++n) ++nbeta_;
        for(int n = 0; n < _nAVirPI[h]; ++n) ++_nAVir;
        for(int n = 0; n < _nBVirPI[h]; ++n) ++_nBVir;
    }

    _aOccC     = shared_ptr<Matrix>(new Matrix("Alpha Occupied MO Coefficients", _nIrreps, _soPI, _nAOccPI));
    _bOccC     = shared_ptr<Matrix>(new Matrix("Beta Occupied MO Coefficients", _nIrreps, _soPI, _nBOccPI));
    _aVirC     = shared_ptr<Matrix>(new Matrix("Alpha Virtual MO Coefficients", _nIrreps, _soPI, _nAVirPI));
    _bVirC     = shared_ptr<Matrix>(new Matrix("Beta Virtual MO Coefficients", _nIrreps, _soPI, _nBVirPI));
    _aScfError = shared_ptr<Matrix>(new Matrix("Alpha SCF Error Vector", _nIrreps, _nAOccPI, _nAVirPI));
    _bScfError = shared_ptr<Matrix>(new Matrix("Beta SCF Error Vector", _nIrreps, _nBOccPI, _nBVirPI));
    _Fa        = shared_ptr<Matrix>(new Matrix("Alpha Fock Matrix", _nIrreps, _soPI, _soPI));
    _Fb        = shared_ptr<Matrix>(new Matrix("Beta Fock Matrix", _nIrreps, _soPI, _soPI));
//    _Ca        = shared_ptr<Matrix>(new Matrix("Alpha MO Coefficients", _nIrreps, _soPI, _soPI));
//    _Cb        = shared_ptr<Matrix>(new Matrix("Beta MO Coefficients", _nIrreps, _soPI, _soPI));
    _Ca        = Process::environment.reference_wavefunction()->Ca();
    _Cb        = Process::environment.reference_wavefunction()->Cb();
    _oldCa     = shared_ptr<Matrix>(new Matrix("Old Alpha MO Coefficients", _nIrreps, _soPI, _soPI));
    _oldCb     = shared_ptr<Matrix>(new Matrix("Old Beta MO Coefficients", _nIrreps, _soPI, _soPI));
    _aKappa    = shared_ptr<Matrix>(new Matrix("Alpha Kappa Matrix", _nIrreps, _soPI, _soPI));
    _bKappa    = shared_ptr<Matrix>(new Matrix("Beta Kappa Matrix", _nIrreps, _soPI, _soPI));
    _aGTau     = shared_ptr<Matrix>(new Matrix("Alpha External Potential Matrix", _nIrreps, _soPI, _soPI));
    _bGTau     = shared_ptr<Matrix>(new Matrix("Beta External Potential Matrix", _nIrreps, _soPI, _soPI));
    _aoS       = shared_ptr<Matrix>(new Matrix("SO Basis Overlap Integrals", _nIrreps, _soPI, _soPI));
    _soH       = shared_ptr<Matrix>(new Matrix("SO basis one-electron integrals", _nIrreps, _soPI, _soPI));
    _sHalfInv  = shared_ptr<Matrix>(new Matrix("SO Basis Inverse Square Root Overlap Matrix", _nIrreps, _soPI, _soPI));
    _aEvals    = shared_ptr<Vector>(new Vector(_nIrreps, _soPI));
    _bEvals    = shared_ptr<Vector>(new Vector(_nIrreps, _soPI));

    _aTau    = new double**[_nIrreps];
    _bTau    = new double**[_nIrreps];
    for(int h = 0; h < _nIrreps; ++h){
        _aTau[h] = block_matrix(_soPI[h], _soPI[h]);
        _bTau[h] = block_matrix(_soPI[h], _soPI[h]);
    }
    

    // Store the AO overlap matrix
    double *sArray = new double[_nTriSo];
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_S, sArray, _nTriSo, 0, 0, outfile);
    _aoS->set(sArray);
    delete [] sArray;

    // Form S^(-1/2) matrix
    Matrix eigvec(_nIrreps, _soPI, _soPI);
    Matrix eigtemp(_nIrreps, _soPI, _soPI);
    Matrix eigtemp2(_nIrreps, _soPI, _soPI);
    Vector eigval(_nIrreps, _soPI);
    _aoS->diagonalize(eigvec, eigval);
    // Convert the eigenvales to 1/sqrt(eigenvalues)
    for (int h=0; h < _nIrreps; ++h) {
        for (int i=0; i < _soPI[h]; ++i) {
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
    for(int h = 0; h < _nIrreps; ++h){
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
