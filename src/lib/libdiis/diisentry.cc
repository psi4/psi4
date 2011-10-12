#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include "diisentry.h"
#include <math.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include <sstream>

using namespace std;
using namespace boost;

namespace psi{

DIISEntry::DIISEntry(std::string label, int ID, int orderAdded,
                     int errorVectorSize, double *errorVector,
                     int vectorSize, double *vector, boost::shared_ptr<PSIO> psio):
        _vectorSize(vectorSize),
        _errorVectorSize(errorVectorSize),
        _vector(vector),
        _errorVector(errorVector),
        _ID(ID),
        _orderAdded(orderAdded),
        _label(label),
        _psio(psio)
{
    double sumSQ = C_DDOT(_errorVectorSize, _errorVector, 1, _errorVector, 1);
    _rmsError = sqrt(sumSQ / _errorVectorSize);
    _dotProducts[_ID] = sumSQ;
    _knownDotProducts[_ID] = true;
    stringstream s;
    s << _label << ":entry " << ID;
    _label = s.str();
}

void
DIISEntry::open_psi_file()
{
    if (_psio->open_check(PSIF_LIBDIIS) == 0) {
        _psio->open(PSIF_LIBDIIS, PSIO_OPEN_OLD);
    }
}

void
DIISEntry::close_psi_file()
{
    if (_psio->open_check(PSIF_LIBDIIS) == 1) {
        _psio->close(PSIF_LIBDIIS, 1);
    }
}

void
DIISEntry::dump_vector_to_disk()
{
    string label = _label + " vector";
    open_psi_file();
    _psio->write_entry(PSIF_LIBDIIS, label.c_str(), (char*)_vector, _vectorSize*sizeof(double));
    free_vector_memory();
}


void
DIISEntry::read_vector_from_disk()
{
    if (_vector == NULL) {
        _vector = new double[_vectorSize];
        string label = _label + " vector";
        open_psi_file();
        _psio->read_entry(PSIF_LIBDIIS, label.c_str(), (char*)_vector, _vectorSize*sizeof(double));
    }
}


void
DIISEntry::dump_error_vector_to_disk()
{
    string label = _label + " error";
    open_psi_file();
    _psio->write_entry(PSIF_LIBDIIS, label.c_str(), (char*)_errorVector, _errorVectorSize*sizeof(double));
    free_error_vector_memory();
}


void
DIISEntry::read_error_vector_from_disk()
{
    if (_errorVector == NULL) {
        _errorVector = new double[_errorVectorSize];
        string label = _label + " error";
        open_psi_file();
        _psio->read_entry(PSIF_LIBDIIS, label.c_str(), (char*)_errorVector, _errorVectorSize*sizeof(double));
    }
}

void
DIISEntry::free_vector_memory()
{
    if (_vector)
        delete[] _vector;
    _vector = NULL;
}


void
DIISEntry::free_error_vector_memory()
{
    if (_errorVector)
        delete[] _errorVector;
    _errorVector = NULL;
}


DIISEntry::~DIISEntry()
{
    if(_vector != NULL)
        delete[] _vector;
    if(_errorVector != NULL)
        delete[] _errorVector;
}

} // Namespace
