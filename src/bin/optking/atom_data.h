/*! 
  \file element_to_Z.h
  \ingroup
  \brief convert element name or symbol to atomic number
*/

#ifndef _element_to_Z_h_
#define _element_to_Z_h_

#include <iostream>
#include <string>
#include <map>

namespace opt {

// atomic number to the atomic symbol
extern const char *Z_to_symbol[];

// atomic number to default mass
extern const double Z_to_mass[];

// convert symbol or name to atomic number
extern double Element_to_Z(std::string lbl);

// convert isotope name to mass
extern double Isotope_to_mass(std::string lbl);

}

#endif

