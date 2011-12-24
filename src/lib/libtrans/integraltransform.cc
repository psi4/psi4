#include "integraltransform.h"
#include "mospace.h"
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <psi4-dec.h>
#include <libciomr/libciomr.h>
#include <libmints/matrix.h>
#include <libmints/molecule.h>
#include <libmints/wavefunction.h>
#define EXTERN
#include <libdpd/dpd.gbl>

using namespace boost;
using namespace psi;

IntegralTransform::IntegralTransform(shared_ptr<Wavefunction> wfn,
                                     SpaceVec spaces,
                                     TransformationType transformationType,
                                     OutputType outputType,
                                     MOOrdering moOrdering,
                                     FrozenOrbitals frozenOrbitals,
                                     bool init):
            initialized_(false),
            _psio(_default_psio_lib_),
            wfn_(wfn),
            _transformationType(transformationType),
            _uniqueSpaces(spaces),
            _moOrdering(moOrdering),
            _outputType(outputType),
            _frozenOrbitals(frozenOrbitals),
            _alreadyPresorted(false),
            _dpdIntFile(PSIF_LIBTRANS_DPD),
            _aHtIntFile(PSIF_LIBTRANS_A_HT),
            _bHtIntFile(PSIF_LIBTRANS_B_HT),
            _nTriSo(0),
            _nTriMo(0),
            _nfzc(0),
            _nfzv(0),
            _spaces(0),
            _labels(0),
            _tolerance(1.0E-16),
            _moIntFileAA(0),
            _moIntFileAB(0),
            _moIntFileBB(0),
            _myDPDNum(1),
            _print(1),
            _zeros(0),
            _sosym(0),
            _cacheFiles(0),
            _cacheList(0),
            _Ca(wfn->Ca()),
            _Cb(wfn->Cb()),
            _keepIwlSoInts(false),
            _keepIwlMoTpdm(true),
            _keepDpdSoInts(false),
            _keepDpdMoTpdm(true),
            _keepHtInts(true),
            _keepHtTpdm(true),
            _tpdmAlreadyPresorted(false)
{
    // Implement set/get functions to customize any of this stuff.  Delayed initialization
    // is possible in case any of these variables need to be changed before setup.
    _printTei      = _print > 5;
    _useIWL        = _outputType == IWLAndDPD || _outputType == IWLOnly;
    _useDPD        = _outputType == IWLAndDPD || _outputType == DPDOnly;
    _iwlAAIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_AA_TEI;
    _iwlABIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_AB_TEI;
    _iwlBBIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_BB_TEI;

    _labels  = Process::environment.molecule()->irrep_labels();
    _nirreps = wfn->nirrep();
    _nmo     = wfn->nmo();
    _nso     = wfn->nso();
    _sopi    = wfn->nsopi();
    _mopi    = wfn->nmopi();
    _clsdpi  = wfn->doccpi();
    _openpi  = wfn->soccpi();
    _frzcpi  = wfn->frzcpi();
    _frzvpi  = wfn->frzvpi();

    common_moinfo_initialize();

    if(init) initialize();
}


IntegralTransform::IntegralTransform(SharedMatrix c,
                                     SharedMatrix i,
                                     SharedMatrix a,
                                     SharedMatrix v,
                                     SpaceVec spaces,
                                     TransformationType transformationType,
                                     OutputType outputType,
                                     MOOrdering moOrdering,
                                     FrozenOrbitals frozenOrbitals,
                                     bool init):
    initialized_(false),
    _psio(_default_psio_lib_),
    _transformationType(transformationType),
    _uniqueSpaces(spaces),
    _moOrdering(moOrdering),
    _outputType(outputType),
    _frozenOrbitals(frozenOrbitals),
    _alreadyPresorted(false),
    _dpdIntFile(PSIF_LIBTRANS_DPD),
    _aHtIntFile(PSIF_LIBTRANS_A_HT),
    _bHtIntFile(PSIF_LIBTRANS_B_HT),
    _nirreps(0),
    _nmo(0),
    _nso(0),
    _nTriSo(0),
    _nTriMo(0),
    _nfzc(0),
    _nfzv(0),
    _spaces(0),
    _labels(0),
    _tolerance(1.0E-16),
    _memory(250 * 1024 * 1024),
    _moIntFileAA(0),
    _moIntFileAB(0),
    _moIntFileBB(0),
    _myDPDNum(1),
    _print(1),
    _zeros(0),
    _sopi(0),
    _sosym(0),
    _mopi(0),
    _clsdpi(0),
    _openpi(0),
    _frzcpi(0),
    _frzvpi(0),
    _cacheFiles(0),
    _cacheList(0),
    _aCorrToPitzer(0),
    _bCorrToPitzer(0),
    _keepIwlSoInts(false),
    _keepIwlMoTpdm(true),
    _keepDpdSoInts(false),
    _keepDpdMoTpdm(true),
    _keepHtInts(true),
    _keepHtTpdm(true),
    _tpdmAlreadyPresorted(false)
{
//    wfn_           = NULL;
//    _Ca            = NULL;
//    _Cb            = NULL;
    _printTei      = _print > 5;
    _useIWL        = _outputType == IWLAndDPD || _outputType == IWLOnly;
    _useDPD        = _outputType == IWLAndDPD || _outputType == DPDOnly;
    _iwlAAIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_AA_TEI;
    _iwlABIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_AB_TEI;
    _iwlBBIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_BB_TEI;

    _nirreps = c->nirrep();
    _nmo     = c->ncol() + i->ncol() + a->ncol() + v->ncol();
    _nso     = i->nrow();
    _sopi    = i->rowspi();   // use i for this since there will always be occupied orbitals
    _mopi    = c->colspi() + i->colspi() + a->colspi() + v->colspi();
    _clsdpi  = i->colspi();
    _openpi  = Dimension(_nirreps); // This is the restricted constructor, there are no unpaired electrons
    _frzcpi  = c->colspi();
    _frzvpi  = v->colspi();

    // Need to smash together the C's only for them to be ripped apart elsewhere.
    std::vector<SharedMatrix > Cs;
    Cs.push_back(c); Cs.push_back(i); Cs.push_back(a); Cs.push_back(v);
    _Ca = Matrix::horzcat(Cs);

    common_moinfo_initialize();

    if(init) initialize();
}


/**
 * Sets up the DPD buffers and performs semicanonicalization, if necessary.
 */
void
IntegralTransform::initialize()
{
    process_spaces();

    // Set up the DPD library
    // TODO implement caching of files
    int numSpaces = _spacesUsed.size();
    int numIndexArrays = numSpaces * (numSpaces - 1) + 5 * numSpaces;
    _cacheFiles = init_int_array(PSIO_MAXUNIT);
    _cacheList  = init_int_matrix(numIndexArrays, numIndexArrays);
    int currentActiveDPD = psi::dpd_default;
    _memory = Process::environment.get_memory();
    dpd_init(_myDPDNum, _nirreps, _memory, 0, _cacheFiles,
            _cacheList, NULL, numSpaces, _spaceArrays);

    // We have to redefine the MO coefficients for a UHF-like treatment
    if(_transformationType == SemiCanonical){
        wfn_->semicanonicalize();
        _Cb = wfn_->Cb();
    }
    process_eigenvectors();

    // Return DPD control to the user
    dpd_set_default(currentActiveDPD);


    // Set up the correlated to Pitzer arrays.  These have to include the occupied core terms, because
    // the reference contributions are already folded into the TPDM.
    _aCorrToPitzer = new int[_nmo];
    if(_transformationType != Restricted){
        _bCorrToPitzer = new int[_nmo];
    }else{
        _bCorrToPitzer = _aCorrToPitzer;
    }
    size_t aCorrCount = 0;
    size_t bCorrCount = 0;
    size_t pitzerOffset = 0;

    // Frozen DOCC
    for(int h = 0; h < _nirreps; ++h){
        for(int n = 0; n < _frzcpi[h]; ++n){
            _aCorrToPitzer[aCorrCount++] = pitzerOffset + n;
            if(_transformationType != Restricted)
                _bCorrToPitzer[bCorrCount++] = pitzerOffset + n;
        }
        pitzerOffset += _mopi[h];
    }
    // Active OCC
    pitzerOffset = 0;
    for(int h = 0; h < _nirreps; ++h){
        for(int n = _frzcpi[h]; n < _clsdpi[h] + _openpi[h]; ++n){
            _aCorrToPitzer[aCorrCount++] = pitzerOffset + n;
        }
        if(_transformationType != Restricted)
            for(int n = _frzcpi[h]; n < _clsdpi[h]; ++n){
                _bCorrToPitzer[bCorrCount++] = pitzerOffset + n;
            }
        pitzerOffset += _mopi[h];
    }
    // Active VIR
    pitzerOffset = 0;
    for(int h = 0; h < _nirreps; ++h){
        for(int n = _clsdpi[h] + _openpi[h]; n < _mopi[h] - _frzvpi[h]; ++n){
            _aCorrToPitzer[aCorrCount++] = pitzerOffset + n;
        }
        if(_transformationType != Restricted)
            for(int n = _clsdpi[h]; n < _mopi[h] - _frzvpi[h]; ++n){
                _bCorrToPitzer[bCorrCount++] = pitzerOffset + n;
            }
        pitzerOffset += _mopi[h];
    }


    initialized_ = true;
}

boost::shared_ptr<PSIO>
IntegralTransform::get_psio() const
{
    return _psio;
}

void
IntegralTransform::set_psio(boost::shared_ptr<PSIO> psio)
{
    _psio = psio;
}

IntegralTransform::~IntegralTransform()
{
    if (initialized_) {
        dpd_close(_myDPDNum);
        free_int_matrix(_cacheList);
        free(_cacheFiles);
        free(_zeros);
    }
}

void IntegralTransform::check_initialized()
{
    if (initialized_ == false)
        throw PSIEXCEPTION("IntegralTransform::check_initialized: This instance is not initialized.");
}

