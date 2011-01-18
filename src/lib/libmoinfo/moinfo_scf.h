#ifndef _psi_src_lib_libmoinfo_moinfo_scf_h_
#define _psi_src_lib_libmoinfo_moinfo_scf_h_

/*! \file
    \ingroup LIBMOINFO
    \brief   This class stores all the basic info regarding MOs
*/

#include <string>

#include "moinfo_base.h"

namespace psi {

class Chkpt;

class MOInfoSCF : public MOInfoBase {
public:
  MOInfoSCF(Options& options_, boost::shared_ptr<Chkpt> chkpt_,bool silent_ = false);
  ~MOInfoSCF();
private:
  void read_mo_spaces();
  void print_mo();
};

extern MOInfoSCF* moinfo_scf;

}

#endif // _psi_src_lib_libmoinfo_moinfo_scf_h_
