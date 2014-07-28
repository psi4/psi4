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

#include "FragOptions.h"
#include "psi4-dec.h"
#include "BSSEer.h"
#include "Fragmenter.h"
#include "Capper.h"

namespace LibFrag{

std::string FragOptions::ToString(const FragMethods& F){
   std::string method;
   switch(F){
      case(USER_DEFINED):{method="user defined";break;}
      case(BOND_BASED):{method="bond based";break;}
      case(DISTANCE_BASED):{method="distance based";break;}
   }
   return method;
}

std::string FragOptions::ToString(const EmbedMethods& E){
   std::string method;
   switch(E){
      case(NO_EMBED):{method="no embedding";break;}
      case(POINT_CHARGE):{method="non-iterative point charges";break;}
      case(ITR_POINT_CHARGE):{method="iterative point charges";break;}
      case(DENSITY):{method="non-iterative density embedding";break;}
      case(ITR_DENSITY):{method="iterative density embedding";break;}
   }
   return method;
}

std::string FragOptions::ToString(const CapMethods& C){
   std::string method;
   switch(C){
      case(NO_CAPS):{method="no capping";break;}
      case(H_REPLACE):{method="strict hydrogen replacement";break;}
      case(H_SHIFTED):{method="shifted hydrogen replacement";break;}
   }
   return method;
}

std::string FragOptions::ToString(const BSSEMethods& B){
   std::string method;
   switch(B){
      case(NO_BSSE):{method="no BSSE correction";break;}
      case(FULL):{
         method="supersystem basis applied to all terms (FULL)";
         break;
      }
      case(MBCPN):{
         method="Many-body counterpoise correction trunctated at order n";
         break;
      }
      case(VMFCN):{
         method="Valiron-Mayer Functional counterpoise correction truncated at order n";
         break;
      }
   }
   return method;
}

void FragOptions::copy(const FragOptions& other){
   this->MBEOrder=other.MBEOrder;
      this->FMethod=other.FMethod;
      this->EMethod=other.EMethod;
      this->CMethod=other.CMethod;
      this->BMethod=other.BMethod;
}

void FragOptions::SetFMethod(const std::string& Frag){
   if(Frag=="user_defined"||Frag=="ud"||Frag=="user")
      FMethod=USER_DEFINED;
   else if(Frag=="bond_based"||Frag=="bond")
      FMethod=BOND_BASED;
   else if(Frag=="distance_based"||Frag=="distance"||Frag=="dist")
      FMethod=DISTANCE_BASED;
   else
      throw psi::PSIEXCEPTION("Unrecognized Fragmentation method");
}

void FragOptions::SetEMethod(const std::string& Embed){
   if(Embed=="none"||Embed=="no_embed"||Embed=="no"||Embed=="false")
      EMethod=NO_EMBED;
   else if(Embed=="point_charge"||Embed=="charges")EMethod=POINT_CHARGE;
   else if(Embed=="iterative"||Embed=="iterative_point_charge")
      EMethod=ITR_POINT_CHARGE;
   else if(Embed=="density")EMethod=DENSITY;
   else if(Embed=="itr_density"||Embed=="iterative_density")
      EMethod=ITR_DENSITY;
}

void FragOptions::SetCMethod(const std::string& Cap){
   if(Cap=="none"||Cap=="no_cap"||Cap=="no"||Cap=="false")CMethod=NO_CAPS;
   else if(Cap=="h_replace"||Cap=="replace")CMethod=H_REPLACE;
   else if(Cap=="h_shifted"||Cap=="shift")CMethod=H_SHIFTED;
   else
      throw psi::PSIEXCEPTION
      ("You requested an unrecognized capping method");
}

void FragOptions::SetBMethod(const std::string& BSSE){
   if(BSSE=="none"||BSSE=="false"||BSSE=="no_bsse"||BSSE=="no")
      BMethod=NO_BSSE;
   else if(BSSE=="full"||BSSE=="bettens")BMethod=FULL;
   else if(BSSE=="MBCPN")BMethod=MBCPN;
   else if(BSSE=="VMFCN")BMethod=VMFCN;
}


void FragOptions::PrintOptions(){
   std::string stars=
         "\n**************************************************************************";
   fprintf(psi::outfile,stars.c_str());
   fprintf(psi::outfile,
"\n******************** Many-Body Expansion (MBE) module ********************"
    );
   fprintf(psi::outfile,stars.c_str());
   fprintf(psi::outfile,
     "\n\nA %d-body expansion is being performed with the following options:\n"
         ,MBEOrder);
   fprintf(psi::outfile,"Fragmenting system via: %s\n",
         (ToString(FMethod)).c_str());
   if(EMethod!=NO_EMBED)fprintf(psi::outfile,"Embedding via: %s\n",
         (ToString(EMethod)).c_str());
   if(CMethod!=NO_CAPS)fprintf(psi::outfile,"Capping via: %s\n",
         (ToString(CMethod)).c_str());
   if(BMethod!=NO_BSSE)fprintf(psi::outfile,"BSSE Corrections via: %s\n",
         (ToString(BMethod)).c_str());
   fprintf(psi::outfile,
   "\n**************************************************************************\n"
    );
}

void FragOptions::DefaultOptions(){
   FMethod=USER_DEFINED;
   EMethod=NO_EMBED;
   CMethod=NO_CAPS;
   BMethod=NO_BSSE;
   MBEOrder=2;
}

boost::shared_ptr<Fragmenter> FragOptions::MakeFragFactory()const{
   boost::shared_ptr<Fragmenter> FragFactory;
   switch(FMethod){
      case(USER_DEFINED):{
          FragFactory=boost::shared_ptr<Fragmenter>(new UDFragmenter);
          break;
      }
      case(BOND_BASED):{
          FragFactory=boost::shared_ptr<Fragmenter>(new BondFragmenter);
          break;
      }
      case(DISTANCE_BASED):{
         FragFactory=boost::shared_ptr<Fragmenter>(new DistFragmenter);
         break;
      }
      default:{
         throw psi::PSIEXCEPTION("Unrecognized fragmentation method");
         break;
      }
   }
   return FragFactory;
}
boost::shared_ptr<BSSEer> FragOptions::MakeBSSEFactory(int natoms)const{
   boost::shared_ptr<BSSEer> BSSEFactory;
   switch(BMethod){
      case(NO_BSSE):{
         break;
      }
      case(FULL):{
       BSSEFactory=boost::shared_ptr<BSSEer>(new FullBSSE(natoms));
       break;
    }
    default:{
       throw psi::PSIEXCEPTION(
             "Unrecognized BSSE correction method or it's not coded yet");
       break;
    }
   }
   return BSSEFactory;
}
boost::shared_ptr<Capper> FragOptions::MakeCapFactory(SharedMol& AMol)const{
   boost::shared_ptr<Capper> Factory;
   switch(CMethod){
      case(NO_CAPS):{
         throw psi::PSIEXCEPTION(
               "You requested no caps be made (or didn't specify a method "
               "for caps), but bonds were broken.  I strongly advise you "
               "reconsider your position on this matter. In fact, I feel so"
               " strongly about this that I have crashed the program");
         break;
      }
      case(H_REPLACE):{
         Factory=boost::shared_ptr<Capper>(new ReplaceAndCap(AMol));
         break;
      }
      case(H_SHIFTED):{
         Factory=boost::shared_ptr<Capper>(new ShiftAndCap(AMol));
         break;
      }
      default:{
         throw psi::PSIEXCEPTION("Unrecognized or capping algorithm.");
         break;
      }
   }
   return Factory;
}

}//End namespace


