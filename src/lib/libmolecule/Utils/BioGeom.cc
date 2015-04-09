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

#include "LibFragMolecule.h"
#include "BioGeom.h"
#include "AutoFxnalGroup.h"
#include "AutoFxnalGroup/PrimRunner.h"
namespace psi {
namespace LibMolecule {
/*
typedef PrimRunner<Amide1C, Amide1N,Amide2, Amide2G, Amide3> FindAmide;
typedef PrimRunner<PNTerminus, NCTerminus, AABackBone,
                   PNTermGly,NCTermGly,Glycine> FindAAPieces;
typedef PrimRunner<ValineR, LeucineR, IsoleucineR, MethionineR, TyrosineR,
      TryptophanR,SerineR,AlanineR> FindAARs;
*/
BioGeom::BioGeom(const Molecule& Mol) :
      OrganicGeom(&Mol,true) {
      Graph Nodes=OrganicGeom::GetGroups();
      Graph::iterator It;
     /* SetRunner<FindAAPieces>::Run(Nodes);
      SetRunner<FindAmide>::Run(Nodes);
      SetRunner<FindAARs>::Run(Nodes);*/
      std::cout<<Nodes.PrintOut()<<std::endl;
      GraphItr It1=Nodes.PrimBegin(),It1End=Nodes.PrimEnd();
      for(;It1!=It1End;++It1)
         std::cout<<"Atom "<<(*(*It1))[0]<<" "<<
                 (*It1)->Type().GenMMType()<<std::endl;
      exit(1);
}


}} //End namespaces

