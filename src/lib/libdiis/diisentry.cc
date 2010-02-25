#include "diisentry.h"
#include <math.h>
#include <libqt/qt.h>

namespace psi{ namespace libdiis{

DIISEntry::DIISEntry(int ID, int orderAdded, int errorVectorSize, double *errorVector,
                     int vectorSize, double *vector):
        _vectorSize(vectorSize),
        _errorVectorSize(errorVectorSize),
        _vector(vector),
        _errorVector(errorVector),
        _ID(ID),
        _orderAdded(orderAdded)
{
    double sumSQ = C_DDOT(_errorVectorSize, _errorVector, 1, _errorVector, 1);
    _rmsError = sqrt(sumSQ / _errorVectorSize);
    _dotProducts[_ID] = sumSQ;
    _knownDotProducts[_ID] = true;
}


void
DIISEntry::dump_to_disk()
{

}


void
DIISEntry::read_from_disk()
{

}


DIISEntry::~DIISEntry()
{
    if(_vector != NULL)
        delete[] _vector;
    if(_errorVector != NULL)
        delete[] _errorVector;
}

}} // Namespaces