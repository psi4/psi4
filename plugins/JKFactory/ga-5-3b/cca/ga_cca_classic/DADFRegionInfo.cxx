#include "DADFRegionInfo.h"

DADFRegionInfo::DADFRegionInfo(const int size) {
  _lowerBounds.resize( size );
  _upperBounds.resize( size );
  _strides.resize( size );
  _data = 0;
}

DADFRegionInfo::DADFRegionInfo(const DADFRegionInfo & original) 
  : _lowerBounds( original._lowerBounds ),
    _upperBounds( original._upperBounds ),
    _strides( original._strides ),
    _data( original._data )
{}

DADFRegionInfo::~DADFRegionInfo(){}

void DADFRegionInfo::setBounds(const int lower[], const int upper[]) {
  unsigned int i;

  for ( i=0; i < _lowerBounds.size(); ++i ) {
    _lowerBounds[i] = lower[i];
    _upperBounds[i] = upper[i];
  }
}

void DADFRegionInfo::getBounds(int lower[], int upper[]) {
  unsigned int i;

  for ( i=0; i < _lowerBounds.size(); ++i ) {
    lower[i] = _lowerBounds[i];
    upper[i] = _upperBounds[i];
  }
}

/** Compare region bounds. */
int DADFRegionInfo::compareBounds(const int lower[], const int upper[]) {
  unsigned int i; int minUpper, maxLower;
  bool identical = true;
  bool overlap = true;

  /** Ranges in an axis overlap if the minimum upper bound is greater
      than or equal to the maximum lower bound.  If overlap occurs in
      all axes, the regions overlap. Identity implies overlap. 
  */
  i = 0;
  while ( (identical || overlap) && i < _lowerBounds.size() ) {
    identical = identical && ( lower[i] == _lowerBounds[i]);
    identical = identical && ( upper[i] == _upperBounds[i]);

    minUpper = (upper[i] < _upperBounds[i] ) ? upper[i] : _upperBounds[i];
    maxLower = (lower[i] > _lowerBounds[i] ) ? lower[i] : _lowerBounds[i];
    overlap = overlap && ( minUpper >= maxLower);
    i++;
  }
  // If regions are identical they will also overlap.  Identity has priority.
  if ( identical ) { return 1; }
  if ( overlap ) { return -1; }

  // The regions are disjoint
  return 0;
}

void DADFRegionInfo::setDataLocation(void * dataPtr, 
				     const int strides[] ) {
  unsigned int i;

  _data = dataPtr;
  for ( i=0; i < _strides.size(); ++i ) { _strides[i] = strides[i]; }
}

void DADFRegionInfo::getDataLocation( void * & dataPtr, int strides[] ) {
  unsigned int i;

  dataPtr = _data;

  for ( i=0; i < _strides.size(); ++i ) { strides[i] = _strides[i]; }
}

void DADFRegionInfo::printRegionInfo() {
  unsigned int i;

  cerr << "(";

  for ( i=0; i < _lowerBounds.size()-1; ++i ) {
    cerr << _lowerBounds[i] << ", ";
  }
  cerr << _lowerBounds[i] << ")";

  cerr << " --> (";

  for ( i=0; i < _upperBounds.size()-1; ++i ) {
    cerr << _upperBounds[i] << ", ";
  }
  cerr << _upperBounds[i] << ")";

  if ( _data != 0 ) { 
    cerr << " at " << _data << " with strides (";

    for ( i=0; i < _strides.size()-1; ++i ) {
      cerr << _strides[i] << ", ";
    }
    cerr << _strides[i] << ")";

    // WARNING! Assumes double data type.  Not general!
    cerr << " first elem = " << *(static_cast<double *>(_data));
  }
}
