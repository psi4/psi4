#include "integraltransform.h"
#include "mospace.h"
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#define EXTERN
#include <libdpd/dpd.gbl>

namespace psi{ namespace libtrans{

IntegralTransform::IntegralTransform(Options &options,
                                     shared_ptr<MOSpace> s1,
                                     shared_ptr<MOSpace> s2,
                                     shared_ptr<MOSpace> s3,
                                     shared_ptr<MOSpace> s4,
                                     TransformationType transformationType,
                                     MOOrdering moOrdering,
                                     OutputType outputType,
                                     FrozenOrbitals frozenOrbitals):
            _transformationType(transformationType),
            _moinfo_initialized(false),
            _moOrdering(moOrdering),
            _outputType(outputType),
            _frozenOrbitals(frozenOrbitals),
            _myDPDNum(1)
{
    _options = options;
    // TODO make sure that these options can be parsed correctly
    _print   = options.get_int("PRINT");
    _memory  = options.get_int("MEMORY") * 1024 * 1024;
    _tolerance = 1.0E-14;
    // For now, just assume that the tei are to be kept
    _deleteIwlSoTei = false;
    _printTei = _print > 5;
    _useIWL    = _outputType == IWLAndDPD || _outputType == IWLOnly;
    _useDPD    = _outputType == IWLAndDPD || _outputType == DPDOnly;

    std::vector<shared_ptr<MOSpace> > spaces;
    spaces.push_back(s1);
    spaces.push_back(s2);
    spaces.push_back(s3);
    spaces.push_back(s4);
    shared_ptr<PSIO> psio(new PSIO); psiopp_ipv1_config(psio);
    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    process_spaces(spaces, psio, chkpt);

    timer_init();

    // Set up the DPD library
    // TODO implement cacheing of files
    int numSpaces = _spaceArrays.size()/2;
    int numIndexArrays = numSpaces * (numSpaces - 1) + 5 * numSpaces;
    _cacheFiles = init_int_array(PSIO_MAXUNIT);
    _cacheList  = init_int_matrix(numIndexArrays, numIndexArrays);
    int currentActiveDPD = psi::dpd_default;
    dpd_init(_myDPDNum, _nirreps, _memory, 0, _cacheFiles,
            _cacheList, NULL, numSpaces, _spaceArrays);
    // Return DPD control to the user
    dpd_set_default(currentActiveDPD);

    // We have to redefine the MO coefficients for a UHF-like treatment
    if(_transformationType == SemiCanonical){
        // This will also build the UHF Fock matrix, which we need
        presort_so_tei();
        semicanonicalize(psio, chkpt);
    }


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
