#ifndef _psi3_bin_optking_except_h_
#define _psi3_bin_optking_except_h_

namespace psi { namespace optking {

class bad_intco_io
{
  std::string keyword;
  int row;
public:
  bad_intco_io(std::string key_in=0, int row_in=0) : keyword(key_in), row(row_in) {}
  void mesg(void) {
    std::cout << "bad_intco_io exception: cannot read " << keyword <<
     " entry row " << row << std::endl;
  }
};

}}

#endif
