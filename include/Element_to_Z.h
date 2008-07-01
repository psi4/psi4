/*! 
  \file element_to_Z.h
  \ingroup
  \brief convert element name or symbol to atomic number
*/

#ifndef _psi_include_element_to_Z_h_
#define _psi_include_element_to_Z_h_

#include <string>
#include <map>

namespace psi {

class Element_to_Z {

  public:
    Element_to_Z() { loaded = false; };
    double operator[](const std::string & elem_sym);

  private:
    bool loaded;
    std::map<std::string,double> element_to_Z;

    void load_values(void);
};

}

#endif

