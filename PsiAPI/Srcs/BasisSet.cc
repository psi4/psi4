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
#include "BasisSet.h"
#include "GaussianBasis.h"
namespace PsiAPI{

void BasisSet::insert(const BasisSet& object){
   Fxns_.push_back(boost::shared_ptr<BasisSet> (object.Clone()));
}

size_t BasisSet::NPrims(size_t Rec)const{
   if(Rec==0)return NBasis();
   size_t Result=0;
   for(size_t i=0;i<Fxns_.size();i++)
      Result+=Fxns_[i]->NPrims(Rec-1);
   return Result;
}

BasisSet::operator std::string()const{
   std::stringstream Result;
   PsiMap<ParamType,double>::const_iterator ParamI=Params_.begin(),
                           ParamEnd=Params_.end();
   for(;ParamI!=ParamEnd;++ParamI)Result<<ParamI->second<<" ";
   const_iterator FxnI=begin(),FxnEnd=end();
   for(;FxnI!=FxnEnd;++FxnI)Result<<(std::string)(*(*FxnI))<<std::endl;
   return Result.str();
}

GaussianPrim::operator std::string()const{
         std::stringstream Result;
         Result<<"Exponent, Coef: "<<GetParam(EXPONENT)<<","
               <<GetParam(COEF);
         return Result.str();
      }

GaussianAO::operator std::string()const{
   std::stringstream Result;
   for(size_t i=0;i<NBasis();i++)
      Result<<(std::string)(*this)[i]<<std::endl;
   return Result.str();
}


GaussianShell::operator std::string()const{
   std::stringstream Result;
   Result<<"L="<<L()<<" "<<(IsPure()?"PURE":"CARTESIAN")<<std::endl;
   Result<<GaussianAO::operator std::string();
   return Result.str();
}

GaussianAtomicBasis::operator std::string()const{
   std::stringstream Result;
   Result<<"Center (a.u.): "<<GetParam(CENTERX)<<" "<<GetParam(CENTERY)
         <<" "<<GetParam(CENTERZ)<<std::endl;
   for(size_t i=0;i<NBasis();i++)
      Result<<(std::string)(*this)[i];
   return Result.str();
}

}


