#ifndef _psi_src_lib_libmemtrix_sblock_vector_h_
#define _psi_src_lib_libmemtrix_sblock_vector_h_

#include <string>

#include "block_vector.h"

namespace psi{ namespace mcscf{

// Smart version of BlockVector
class SBlockVector
{
public:
  SBlockVector();
  SBlockVector(std::string label, int nirreps, int*& rows_size);
  SBlockVector(std::string label, int nirreps, vecint& rows_size);
  ~SBlockVector() {if(block_vector_) block_vector_->subtract_reference();}

  // Manual allocation
  void allocate(std::string label, int nirreps, int*& rows_size);
  void allocate(std::string label, int nirreps, vecint& rows_size);

  // Copy constructor and assignment operator
  SBlockVector            (const   SBlockVector& src);
  SBlockVector& operator= (const   SBlockVector& src);

  // Allow access to the implementation object
  const BlockVector* operator-> () const {return block_vector_;}
  BlockVector*       operator-> ()       {return block_vector_;}

  // Access the implementation object
  BlockVector* getBlockVector() {return block_vector_;}

  // Checking functions
  bool is_allocated() {return (block_vector_);}
  void check(const char* cstr);

  void copy(SBlockVector& source);
private:
  SBlockVector(BlockVector* block_vector);

  BlockVector* block_vector_;
};

}}

#endif // _psi_src_lib_libmemtrix_sblock_vector_h_
