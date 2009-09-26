#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>

#include "orbital_space.h"
#include <libutil/libutil.h>

extern FILE *outfile;

namespace psi {

using namespace std;

// The default space is all
OrbitalSpace::OrbitalSpace()
{
  label  = "all";
  symbol = "n";
}

OrbitalSpace::OrbitalSpace(std::string label_,std::string symbol_)
{
  label = label_;
  if(symbol_.size() == 1){
    symbol = symbol_;
  }else{
    std::string message = "Cannot create orbital space \"" + label_
                        + "\" with symbol \"" + symbol_
                        + "\". Please use only one character.";
    print_error(outfile,message,__FILE__,__LINE__);
  }
}

void OrbitalSpace::add_subspace(OrbitalSpace space)
{
  subspaces.push_back(space);
}

void OrbitalSpace::print(std::string str)
{
  fprintf(outfile,"\n%s[%s]=[%s]",str.c_str(),label.c_str(),symbol.c_str());
  if(subspaces.size() > 0){
    fprintf(outfile," = {");
    for(OSIt it = subspaces.begin(); it != subspaces.end(); ++it)
      it->print(str + "  ");
    fprintf(outfile,"\n%s}",str.c_str());
  }
  fflush(outfile);
}

}
