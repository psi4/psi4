#ifndef _psi_src_lib_libmoinfo_model_space_h_
#define _psi_src_lib_libmoinfo_model_space_h_

/*! \file    model_space.h
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding the model space
*/

#include "slater_determinant.h"

namespace psi {

class MOInfo;

class ModelSpace{
public:
  ModelSpace(MOInfo* moinfo_obj_);
  ~ModelSpace();
  void print();
private:
  void startup();
  void cleanup();
  void build();
  void classify();

  int wfn_sym;
  std::vector<SlaterDeterminant> determinants;
  std::vector<int>  closed_to_all;  // closed-shell determinants
  std::vector<int>  opensh_to_all;  // open-shell   determinants
  std::vector<int>  unique_to_all;  // spin-unique  determinants
  MOInfo* moinfo_obj;
};

extern ModelSpace  *model_space;

}

#endif // _psi_src_lib_libmoinfo_model_space_h_
