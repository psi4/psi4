#include <cstring>
#include <iostream>

#include <libqt/qt.h>
#include <libutil/libutil.h>

#include "vector_base.h"

#include <psi4-dec.h>

extern FILE* outfile;

namespace psi{ namespace mcscf{

extern MemoryManager* memory_manager;

VectorBase::VectorBase(int elements) : elements_(elements),vector_(NULL)
{
  allocate1(double,vector_,elements_);
}

VectorBase::~VectorBase()
{
  release1(vector_);
}

void VectorBase::print()
{
  fprintf(outfile,"\n  ");
  for(size_t i = 0 ; i < elements_; ++i){
    fprintf(outfile,"%10.6f",vector_[i]);
  }
}

void VectorBase::copy(VectorBase& source)
{
  for(size_t i = 0 ; i < elements_; ++i)
    vector_[i] = source.vector_[i];
}

}}
