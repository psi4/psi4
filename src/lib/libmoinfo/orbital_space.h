#ifndef _psi_src_lib_libmoinfo_orbital_space_h_
#define _psi_src_lib_libmoinfo_orbital_space_h_

/*! \file
    \ingroup (PSIMRCC)
    \brief   This class is used to specify a set of orbitals
*/

#include <string>
#include <vector>

namespace psi{

class OrbitalSpace{
  typedef std::vector<OrbitalSpace>::iterator OSIt;
public:
  OrbitalSpace();
  explicit OrbitalSpace(std::string label_,std::string symbol_);
  void print(std::string str = "");
  void add_subspace(OrbitalSpace space);
  std::string get_label();
private:
  std::string   label;
  std::string   symbol;
  std::vector<std::string>   input_labels;
  int           n;
  std::vector<int>   orbpi;
  std::vector<OrbitalSpace> subspaces;
};

} /* End Namespace */

#endif // _psi_src_lib_libmoinfo_orbital_space_h_
