#ifndef _psi_src_lib_libmemtrix_vector_base_h_
#define _psi_src_lib_libmemtrix_vector_base_h_

#include <libutil/memory_manager.h>
#include <cstring> // for size_t

namespace psi{ namespace mcscf{

class VectorBase
{
public:
  VectorBase();
  VectorBase(int rows);
  ~VectorBase();

  //Inlines
  int     get_elements()    {return(elements_);}
  void    set(int i, double value) {vector_[i]  = value;}
  void    add(int i, double value) {vector_[i] += value;}
  double  get(int i)               {return(vector_[i]);}
  double* get_vector()  {return(vector_);}

  void print();
  void copy(VectorBase& source);
private:
  // Vector size
  size_t  elements_;

  // Vector data
  double* vector_;
};

}}

#endif // _psi_src_lib_libmemtrix_vector_base_h_
