#include "DADFAxisInfo.h"

/******************************************************************************
 * CollapsedAxisInfo
 *
 * Very little to implement here, since a collapsed axis has no info.
 *****************************************************************************/

CollapsedAxisInfo::~CollapsedAxisInfo(){}

DistArrayTemplate::DistType CollapsedAxisInfo::getDistType() {
  return DistArrayTemplate::Collapsed;
}

void CollapsedAxisInfo::printAxisInfo() {
  cerr << "(collapsed axis has no info)";
}

/******************************************************************************
 * BlockAxisInfo
 *****************************************************************************/

BlockAxisInfo::BlockAxisInfo(){
  _blockSize = -1;
  _first = -1;
}

BlockAxisInfo::BlockAxisInfo(const BlockAxisInfo & original) 
  : _blockSize( original._blockSize ),
    _first( original._first)
{}

BlockAxisInfo::~BlockAxisInfo(){}

DistArrayTemplate::DistType BlockAxisInfo::getDistType() {
  return DistArrayTemplate::Block;
}

void BlockAxisInfo::setDistParameters(const int blockSize, const int first) {
  _blockSize = blockSize;
  _first     = first;
}

void BlockAxisInfo::getDistParameters(int & blockSize, int & first) {
  blockSize = _blockSize;
  first     = _first;
}

void BlockAxisInfo::printAxisInfo() {
  cerr << "block size = " << _blockSize << 
    ", first block on process = " << _first;
}

/******************************************************************************
 * GenBlockAxisInfo
 *****************************************************************************/

GenBlockAxisInfo::GenBlockAxisInfo(const int size) {
  _blockSizes.resize( size );
}

GenBlockAxisInfo::GenBlockAxisInfo(const GenBlockAxisInfo & original) 
  : _blockSizes( original._blockSizes )
{}

GenBlockAxisInfo::~GenBlockAxisInfo(){}

DistArrayTemplate::DistType GenBlockAxisInfo::getDistType() {
  return DistArrayTemplate::GenBlock;
}

void GenBlockAxisInfo::setDistParameters(const int blockSizes[]) {
  unsigned int i;

  for ( i=0; i < _blockSizes.size(); ++i ) {
    _blockSizes[i] = blockSizes[i];
  }
}

void GenBlockAxisInfo::getDistParameters(int blockSizes[]) {
  unsigned int i;

  for ( i=0; i < _blockSizes.size(); ++i ) {
    blockSizes[i] = _blockSizes[i];
  }
}

void GenBlockAxisInfo::printAxisInfo() {
  unsigned int i;

  cerr << "block sizes = (";

  for ( i=0; i < _blockSizes.size()-1; ++i ) {
    cerr << _blockSizes[i] << ", ";
  }
  cerr << _blockSizes[i] << ")";
}

/******************************************************************************
 * ImplicitAxisInfo
 *****************************************************************************/

ImplicitAxisInfo::ImplicitAxisInfo(const int size) {
  _map.resize( size );
}

ImplicitAxisInfo::ImplicitAxisInfo(const ImplicitAxisInfo & original) 
  : _map( original._map )
{}

ImplicitAxisInfo::~ImplicitAxisInfo(){}

DistArrayTemplate::DistType ImplicitAxisInfo::getDistType() {
  return DistArrayTemplate::Implicit;
}

void ImplicitAxisInfo::setDistParameters(const int map[]) {
  unsigned int i;

  for ( i=0; i < _map.size(); ++i ) {
    _map[i] = map[i];
  }
}

void ImplicitAxisInfo::getDistParameters(int map[]) {
  unsigned int i;

  for ( i=0; i < _map.size(); ++i ) {
    map[i] = _map[i];
  }
}

void ImplicitAxisInfo::printAxisInfo() {
  unsigned int i;

  cerr << "element to process map = (";

  for ( i=0; i < _map.size()-1; ++i ) {
    cerr << _map[i] << ", ";
  }
  cerr << _map[i] << ")";
}
