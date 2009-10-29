#include "integraltransform.h"
#include "mospace.h"
#include <libdpd/dpd.h>
#include <libqt/qt.h>

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
            _aFzcOp(NULL),
            _bFzcOp(NULL),
            _aFzcD(NULL),
            _bFzcD(NULL)
{
    _options = options;
    // TODO make sure that these options can be parsed correctly
    _print   = options.get_int("Print");
    _memory  = options.get_int("MEMORY") * 1024 * 1024;
    _tolerance = 1.0E-14;
    // For now, just assume that the tei are to be kept
    _deleteIwlSoTei = false;
    _printTei = _print > 5;
    
    // TODO implement semicanonicalization
    std::vector<shared_ptr<MOSpace> > spaces;
    spaces.push_back(s1);
    spaces.push_back(s2);
    spaces.push_back(s3);
    spaces.push_back(s4);
    process_spaces(spaces);

    timer_init();

    // Set up the DPD library
    // TODO implement cacheing of files
    int numSpaces = _spaceArrays.size()/2;
    int numIndexArrays = numSpaces * (numSpaces - 1) + 5 * numSpaces;
    _cacheFiles = init_int_array(PSIO_MAXUNIT);
    _cacheList  = init_int_matrix(numIndexArrays, numIndexArrays);
    dpd_init(1, _nirreps, _memory, 0, _cacheFiles,
            _cacheList, NULL, numSpaces, _spaceArrays);

}

IntegralTransform::~IntegralTransform()
{
    for(int h = 0; h < _nirreps; ++h) {
        if(_sopi[h] && _mopi[h]) {
            free_block(_Ca[h]);
        }
    }
    free_block(_fullCa);
    delete [] _Ca;
    // Restricted transformations never allocated the beta matrices
    if(_transformationType == Unrestricted){
        for(int h = 0 ; h < _nirreps; ++h) {
            if(_sopi[h] && _mopi[h]) {
                free_block(_Cb[h]);
            }
        }
        free_block(_fullCb);
        delete [] _Cb;
    }

    if(_aFzcD != NULL)  free(_aFzcD);
    if(_aFzcOp != NULL) free(_aFzcOp);
    if(_transformationType != Restricted){
        if(_bFzcD != NULL)  free(_bFzcD);
        if(_bFzcOp != NULL) free(_bFzcOp);
    }

    dpd_close(1);
    free_int_matrix(_cacheList);
    free(_cacheFiles);
    timer_done();
}

}} // End namespaces
