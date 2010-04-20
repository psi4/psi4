#include <cstdlib>
#include <psifiles.h>
#include <cstdio>

#include "sblock_vector.h"

#include <psi4-dec.h>

namespace psi{ namespace mcscf{

SBlockVector::SBlockVector()
 : block_vector_(0)
{
}

SBlockVector::SBlockVector(std::string label, int nirreps, int*& rows_size)
 : block_vector_(0)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

SBlockVector::SBlockVector(std::string label, int nirreps, vecint& rows_size)
 : block_vector_(0)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

SBlockVector::SBlockVector(BlockVector* block_vector)
 : block_vector_(block_vector)
{
  block_vector_->add_reference();
}

SBlockVector::SBlockVector(const SBlockVector& src)
{
  block_vector_ = src.block_vector_;
  block_vector_->add_reference();
}

void SBlockVector::allocate(std::string label, int nirreps, int*& rows_size)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

void SBlockVector::allocate(std::string label, int nirreps, vecint& rows_size)
{
  block_vector_ = new BlockVector(label,nirreps,rows_size);
  block_vector_->add_reference();
}

SBlockVector& SBlockVector::operator= (const SBlockVector& src)
{
  // Make sure we don't copy ourself!
  if (block_vector_ == src.block_vector_) return *this;

  block_vector_->subtract_reference();  // Remove reference from existing object
  block_vector_ = src.block_vector_;
  block_vector_->add_reference();       // Add reference to our new object

  return *this;
}

void SBlockVector::check(const char* cstr)
{
  if(!is_allocated()){
    fprintf(outfile,"\n\n  Error: SBlockVector operation '%s' is using an uninitialized matrix",cstr);
    fflush(outfile);
    exit(PSI_RETURN_FAILURE);
  }
}

void SBlockVector::copy(SBlockVector& source)
{
  block_vector_->copy(*source.block_vector_);
}

}}
