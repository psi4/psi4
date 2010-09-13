#include "dcft.h"
#include "defines.h"
#include <vector>
#include <liboptions/liboptions.h>
#include <libpsio/psio.hpp>
#include <libtrans/integraltransform.h>
#include <libtrans/mospace.h>
#include <libdpd/dpd.h>
#include <libiwl/iwl.hpp>

namespace psi{ namespace dcft{

/**
 * Sets aside memory for the various arrays and matrices needed, and obtains
 * important information from the checkpoint file.  Must be run after an SCF
 * computation for this reason.
 */
void
DCFTSolver::init_moinfo()
{
    _nMo     = _chkpt->rd_nmo();
    _moPI    = _chkpt->rd_orbspi();
    _openPI  = _chkpt->rd_openpi();
    _clsdPI  = _chkpt->rd_clsdpi();
    _frzcPI  = _chkpt->rd_frzcpi();
    _frzvPI  = _chkpt->rd_frzvpi();
    _nAOcc   = _nBOcc =  _nAVir = _nBVir = 0;
    for(int h = 0; h < _nIrreps; ++h){
        _nAOccPI[h] = _clsdPI[h] + _openPI[h];
        _nBOccPI[h] = _clsdPI[h];
        _nAVirPI[h] = _moPI[h] - _clsdPI[h] - _openPI[h] - _frzvPI[h];
        _nBVirPI[h] = _moPI[h] - _clsdPI[h] - _frzvPI[h];
        for(int n = 0; n < _nAOccPI[h]; ++n) ++_nAOcc;
        for(int n = 0; n < _nBOccPI[h]; ++n) ++_nBOcc;
        for(int n = 0; n < _nAVirPI[h]; ++n) ++_nAVir;
        for(int n = 0; n < _nBVirPI[h]; ++n) ++_nBVir;
    }
    _aOccC.init(_nIrreps, _soPI, _nAOccPI, "Alpha Occupied MO Coefficients");
    _bOccC.init(_nIrreps, _soPI, _nBOccPI, "Beta Occupied MO Coefficients");
    _aVirC.init(_nIrreps, _soPI, _nAVirPI, "Alpha Virtual MO Coefficients");
    _bVirC.init(_nIrreps, _soPI, _nBVirPI, "Beta Virtual MO Coefficients");
    _aScfError.init(_nIrreps, _nAOccPI, _nAVirPI, "Alpha SCF Error Vector");
    _bScfError.init(_nIrreps, _nBOccPI, _nBVirPI, "Beta SCF Error Vector");

}


/**
 * Reads the orbital information that can be determined before the SCF procedure
 * and initializes SO matrices.
 */
void
DCFTSolver::read_checkpoint()
{
    _nSo     = _chkpt->rd_nso();
    _soPI    = _chkpt->rd_sopi();
    _nTriSo  = _nSo * (_nSo + 1) / 2;
    _nIrreps = _chkpt->rd_nirreps();

    _Fa.init(_nIrreps, _soPI, _soPI, "Alpha Fock Matrix");
    _Fb.init(_nIrreps, _soPI, _soPI, "Beta Fock Matrix");
    _Ca.init(_nIrreps, _soPI, _soPI, "Alpha MO Coefficients");
    _Cb.init(_nIrreps, _soPI, _soPI, "Beta MO Coefficients");
    _oldCa.init(_nIrreps, _soPI, _soPI, "Old Alpha MO Coefficients");
    _oldCb.init(_nIrreps, _soPI, _soPI, "Old Beta MO Coefficients");
    _aKappa.init(_nIrreps, _soPI, _soPI, "Alpha Kappa Matrix");
    _bKappa.init(_nIrreps, _soPI, _soPI, "Beta Kappa Matrix");
    _aGTau.init(_nIrreps, _soPI, _soPI, "Alpha External Potential Matrix");
    _bGTau.init(_nIrreps, _soPI, _soPI, "Beta External Potential Matrix");
    _aoS.init(_nIrreps, _soPI, _soPI, "SO Basis Overlap Integrals");
    _soH.init(_nIrreps, _soPI, _soPI, "SO basis one-electron integrals");
    _sHalfInv.init(_nIrreps, _soPI, _soPI, "SO Basis Inverse Square Root Overlap Matrix");
    _aEvals.init(_nIrreps, _soPI);
    _bEvals.init(_nIrreps, _soPI);

    _nAOccPI = new int[_nIrreps];
    _nBOccPI = new int[_nIrreps];
    _nAVirPI = new int[_nIrreps];
    _nBVirPI = new int[_nIrreps];

    _aTau    = new double**[_nIrreps];
    _bTau    = new double**[_nIrreps];
    for(int h = 0; h < _nIrreps; ++h){
        _aTau[h]    = block_matrix(_soPI[h], _soPI[h]);
        _bTau[h]    = block_matrix(_soPI[h], _soPI[h]);
    }
    
    // Read information from checkpoint
    _eNuc      = _chkpt->rd_enuc();
    double *zVals = _chkpt->rd_zvals();
    int    nAtoms = _chkpt->rd_natom();


    // Read in DOCC and SOCC from memory
    if(_options["DOCC"].has_changed()){
        _inputDocc = true;
        for (int h = 0; h < _nIrreps; ++h)
            _nAOccPI[h] = _nBOccPI[h] = _options["DOCC"][h].to_integer();
    }else{
        _inputDocc = false;
        for (int h = 0; h < _nIrreps; ++h)
            _nAOccPI[h] = _nBOccPI[h] = 0;
    }
    if(_options["SOCC"].has_changed()){
        _inputSocc = true;
        for (int h = 0; h < _nIrreps; ++h)
            _nAOccPI[h] += _options["SOCC"][h].to_integer();
    }else{
        _inputSocc = false;
    }

    // Determine the number of electrons in the system
    int multiplicity = 1;
    int charge = _options.get_int("CHARGE");
    int nElec  = 0;
    for (int i=0; i<nAtoms; ++i)
        nElec += (int)zVals[i];
    nElec -= charge;

    delete[] zVals;
    
    // If the user told us the multiplicity, read it from the input
    if(_options["MULTP"].has_changed()){
        multiplicity = _options.get_int("MULTP");
    }else{
        if(nElec % 2){
            multiplicity = 2;
            // There are an odd number of electrons
            fprintf(outfile,"\tThere are an odd number of electrons - assuming doublet.\n"
                            "\tSpecify the multiplicity with the MULTP option in the\n"
                            "\tinput if this is incorrect\n\n");
        }else{
            multiplicity = 1;
            // There are an even number of electrons
            fprintf(outfile,"\tThere are an even number of electrons - assuming singlet.\n"
                            "\tSpecify the multiplicity with the MULTP option in the\n"
                            "\tinput if this is incorrect\n\n");
        }
    }
    // Make sure that the multiplicity is reasonable
    if(multiplicity - 1 > nElec){
        char *str = new char[100];
        sprintf(str, "There are not enough electrons for multiplicity = %d, \n"
                     "please check your input and use the MULTP keyword", multiplicity);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }
    if(multiplicity % 2 == nElec % 2){
        char *str = new char[100];
        sprintf(str, "A multiplicity of %d with %d electrons is impossible.\n"
                     "Please check your input and use the MULTP and/or CHARGE keywords",
                     multiplicity, nElec);
        throw SanityCheckError(str, __FILE__, __LINE__);
        delete [] str;
    }
    _nBOcc  = (nElec - multiplicity + 1)/2;
    _nAOcc  = _nBOcc + multiplicity - 1;

    // Store the AO overlap matrix
    double *sArray = new double[_nTriSo];
    IWL::read_one(_psio.get(), PSIF_OEI, PSIF_SO_S, sArray, _nTriSo, 0, 0, outfile);
    _aoS.set(sArray);
    delete [] sArray;

    // Form S^(-1/2) matrix
    Matrix eigvec(_nIrreps, _soPI, _soPI);
    Matrix eigtemp(_nIrreps, _soPI, _soPI);
    Matrix eigtemp2(_nIrreps, _soPI, _soPI);
    Vector eigval(_nIrreps, _soPI);
    _aoS.diagonalize(eigvec, eigval);
    // Convert the eigenvales to 1/sqrt(eigenvalues)
    for (int h=0; h < _nIrreps; ++h) {
        for (int i=0; i < _soPI[h]; ++i) {
            double scale = 1.0 / sqrt(eigval.get(h, i));
            eigval.set(h, i, scale);
        }
    }
    eigtemp2.set(eigval);
    eigtemp.gemm(false, true, 1.0, eigtemp2, eigvec, 0.0);
    _sHalfInv.gemm(false, false, 1.0, eigvec, eigtemp, 0.0);
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
    delete [] _openPI;
    delete [] _clsdPI;
    delete [] _frzcPI;
    delete [] _frzvPI;
    delete [] _soPI;
    delete [] _moPI;
    delete [] _nAOccPI;
    delete [] _nBOccPI;
    delete [] _nAVirPI;
    delete [] _nBVirPI;
}

}}//Namespaces
