#include "integraltransform.h"
#include "mospace.h"
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#define EXTERN
#include <libdpd/dpd.gbl>

namespace psi{ namespace libtrans{

IntegralTransform::IntegralTransform(SpaceVec spaces,
                                     TransformationType transformationType,
                                     OutputType outputType,
                                     MOOrdering moOrdering,
                                     FrozenOrbitals frozenOrbitals,
                                     bool init):
            _transformationType(transformationType),
            _moOrdering(moOrdering),
            _outputType(outputType),
            _uniqueSpaces(spaces),
            _frozenOrbitals(frozenOrbitals),
            _chkpt(_default_chkpt_lib_),
            _psio(_default_psio_lib_),
            _Ca(NULL),
            _Cb(NULL)
{
    // Implement set/get functions to customize any of this stuff.  Delayed initialization
    // is possible in case any of these variables need to be changed before setup.
    _myDPDNum      = 1;
    _print         = 1;
    _memory        = 2000 * 1024 * 1024;
    _tolerance     = 1.0E-14;
    _keepDpdSoInts = false;
    _keepIwlSoInts = false;
    _keepHtInts    = true;
    _printTei      = _print > 5;
    _useIWL        = _outputType == IWLAndDPD || _outputType == IWLOnly;
    _useDPD        = _outputType == IWLAndDPD || _outputType == DPDOnly;
    _dpdIntFile    = PSIF_LIBTRANS_DPD;
    _aHtIntFile    = PSIF_LIBTRANS_A_HT;
    _bHtIntFile    = PSIF_LIBTRANS_B_HT;
    _iwlAAIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_AA_TEI;
    _iwlABIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_AB_TEI;
    _iwlBBIntFile  = _transformationType == Restricted ? PSIF_MO_TEI : PSIF_MO_BB_TEI;

    if(init) initialize();
}

/**
 * Sets up the DPD buffers and performs semicanonicalization, if necessary.
 */
void
IntegralTransform::initialize()
{
    timer_init();

    raid_checkpoint();
    process_spaces();
    
    // Set up the DPD library
    // TODO implement cacheing of files
    int numSpaces = _spacesUsed.size();
    int numIndexArrays = numSpaces * (numSpaces - 1) + 5 * numSpaces;
    _cacheFiles = init_int_array(PSIO_MAXUNIT);
    _cacheList  = init_int_matrix(numIndexArrays, numIndexArrays);
    int currentActiveDPD = psi::dpd_default;
    dpd_init(_myDPDNum, _nirreps, _memory, 0, _cacheFiles,
            _cacheList, NULL, numSpaces, _spaceArrays);

    // We have to redefine the MO coefficients for a UHF-like treatment
    if(_transformationType == SemiCanonical){
        _Ca = new double**[_nirreps];
        _Cb = _Ca;
        for(int h = 0; h < _nirreps; ++h){
            _Ca[h] = _chkpt->rd_scf_irrep(h);
        }
        // This will also build the UHF Fock matrix, which we need
        generate_oei();
        semicanonicalize();
        // This second call does everything that generate_oei() does, but also
        // sorts the two electron integrals for the transformation, which saves
        // us a pass through the SO integral file.
        presort_so_tei();
    }
    process_eigenvectors();

    // Return DPD control to the user
    dpd_set_default(currentActiveDPD);

}


IntegralTransform::~IntegralTransform()
{
    //TODO clean up everything (use valgrind)
    for(int h = 0; h < _nirreps; ++h) {
        if(_sopi[h] && _mopi[h]) {
            free_block(_Ca[h]);
        }
    }
    delete [] _Ca;
    // Restricted transformations never allocated the beta matrices
    if(_transformationType != Restricted){
        for(int h = 0 ; h < _nirreps; ++h) {
            if(_sopi[h] && _mopi[h]) {
                free_block(_Cb[h]);
            }
        }
        delete [] _Cb;
    }


    dpd_close(_myDPDNum);
    free_int_matrix(_cacheList);
    free(_cacheFiles);
    timer_done();
}

}} // End namespaces
