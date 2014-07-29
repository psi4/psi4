/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include "model_space.h"
#include "moinfo.h"
#include <cstdio>
#include "psi4-dec.h"
namespace psi {

extern FILE *outfile;

ModelSpace::ModelSpace(MOInfo* moinfo_obj_) : moinfo_obj(moinfo_obj_)
{
  startup();
  build();
  classify();
//  print();
}

ModelSpace::~ModelSpace()
{
  cleanup();
}

void ModelSpace::startup()
{
  wfn_sym = moinfo_obj->get_wfn_sym();
}

void ModelSpace::cleanup()
{
}

void ModelSpace::print()
{
  psi::fprintf(outfile,"\n\n  Model space:");
  psi::fprintf(outfile,"\n  ------------------------------------------------------------------------------");
  for(int mu = 0; mu < determinants.size(); ++mu){
    psi::fprintf(outfile,"\n  %2d %s",mu,determinants[mu].get_label().c_str());
  }
  psi::fprintf(outfile,"\n\n  Closed-shell to model space mapping");
  for(int mu = 0; mu < closed_to_all.size(); ++mu){
    psi::fprintf(outfile,"\n  %d -> %d",mu,closed_to_all[mu]);
  }
  psi::fprintf(outfile,"\n\n  Open-shell to model space mapping");
  for(int mu = 0; mu < opensh_to_all.size(); ++mu){
    psi::fprintf(outfile,"\n  %d -> %d",mu,opensh_to_all[mu]);
  }

}

}
