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
#include "Embedder.h"
#include "libparallel/TableSpecs.h"
#include <sstream>
#include "libmints/molecule.h"
namespace psi{
namespace LibFrag{

void FragOptions::SetDefaults(){
   Methods_[USER_DEFINED]="User Defined";
   Methods_[BOND_BASED]="Bond Based";
   Methods_[DISTANCE_BASED]="Distance Based";
   Converter_["none"]=USER_DEFINED;
   Converter_["user_defined"]=USER_DEFINED;
   Converter_["ud"]=USER_DEFINED;
   Converter_["user"]=USER_DEFINED;
   Converter_["bond_based"]=BOND_BASED;
   Converter_["bond"]=BOND_BASED;
   Converter_["distance_based"]=DISTANCE_BASED;
   Converter_["distance"]=DISTANCE_BASED;
   Converter_["dist"]=DISTANCE_BASED;
}
boost::shared_ptr<Fragmenter> FragOptions::MakeFactory(SharedMol& AMol)const{
   boost::shared_ptr<Fragmenter> FragFactory;
   switch(DaMethod_){
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

void EmbedOptions::SetDefaults(){
   Methods_[NO_EMBED]="None Specified";
   Methods_[POINT_CHARGE]="Atomic-Centered Point Charges";
   Methods_[ITR_POINT_CHARGE]="Iterative Atomic-Centered Poitn Charges";
   Methods_[DENSITY]="Frozen Density";
   Methods_[ITR_DENSITY]="Iterative Density";
   Converter_["none"]=NO_EMBED;
   Converter_["no_embed"]=NO_EMBED;
   Converter_["no"]=NO_EMBED;
   Converter_["false"]=NO_EMBED;
   Converter_["point_charge"]=POINT_CHARGE;
   Converter_["charges"]=POINT_CHARGE;
   Converter_["charge"]=POINT_CHARGE;
   Converter_["iterative"]=ITR_POINT_CHARGE;
   Converter_["itr_charges"]=ITR_POINT_CHARGE;
   Converter_["itr_charge"]=ITR_POINT_CHARGE;
   Converter_["iterative_point_charge"]=ITR_POINT_CHARGE;
   Converter_["density"]=DENSITY;
   Converter_["itr_density"]=ITR_DENSITY;
   Converter_["iterative_density"]=ITR_DENSITY;
}
boost::shared_ptr<Embedder> EmbedOptions::MakeFactory(SharedMol& AMol)const{
   boost::shared_ptr<Embedder> Factory;
   switch(DaMethod_){
      case(NO_EMBED):{
         Factory=boost::shared_ptr<Embedder>(new NullEmbedder());
         break;
      }
      case(POINT_CHARGE):{
         Factory=boost::shared_ptr<Embedder>(new APCEmbedder(AMol,false));
         break;
      }
      case(ITR_POINT_CHARGE):{
         Factory=boost::shared_ptr<Embedder>(new APCEmbedder(AMol,true));
         break;
      }
      default:{
         throw psi::PSIEXCEPTION("Unrecognized or un-coded Embedding method");
         break;
      }
   }
   return Factory;
}

void CapOptions::SetDefaults(){
   Methods_[NO_CAPS]="None Specified";
   Methods_[H_REPLACE]="Literal Replace Using H Atom";
   Methods_[H_SHIFTED]="Shifted H Atom Cap";
   Converter_["none"]=NO_CAPS;
   Converter_["no_cap"]=NO_CAPS;
   Converter_["no"]=NO_CAPS;
   Converter_["false"]=NO_CAPS;
   Converter_["h_replace"]=H_REPLACE;
   Converter_["replace"]=H_REPLACE;
   Converter_["shift"]=H_SHIFTED;
   Converter_["h_shifted"]=H_SHIFTED;
}
boost::shared_ptr<Capper> CapOptions::MakeFactory(SharedMol& AMol)const{
   boost::shared_ptr<Capper> Factory;
   switch(DaMethod_){
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

void BSSEOptions::SetDefaults(){
   Methods_[NO_BSSE]="None Specified";
   Methods_[FULL]="Full Supersystem Basis Set";
   Methods_[MBCPN]="Many-Body Counterpoise Correction";
   Methods_[VMFCN]="Truncated Valiron-Mayer Functional Counterpoise Correction";
   Converter_["none"]=NO_BSSE;
   Converter_["false"]=NO_BSSE;
   Converter_["no_bsse"]=NO_BSSE;
   Converter_["no"]=NO_BSSE;
   Converter_["full"]=FULL;
   Converter_["mbcpn"]=MBCPN;
   Converter_["vmfcn"]=VMFCN;
   DaMethod_=NO_BSSE;
}
boost::shared_ptr<BSSEer> BSSEOptions::MakeFactory(SharedMol& AMol)const{
  int natoms=AMol->natom();
   boost::shared_ptr<BSSEer> BSSEFactory;
   switch(DaMethod_){
      case(NO_BSSE):{
         BSSEFactory=boost::shared_ptr<NullBSSE>(new NullBSSE());
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

void LibFragOptions::PrintOptions(){
   outfile->MakeBanner("Fragmentation Method Module");
   outfile->Printf("\n");

   (*outfile)<<"Performing the following fragment method expansion:\n\n";
   int noptions=5;
   std::vector<std::string> Titles;
   Titles.push_back("Option Name");
   Titles.push_back("Choice");

   std::vector<std::string> OptionNames;
   std::vector<std::string> Choices;

   OptionNames.push_back("Expansion Order:");
   std::stringstream Order;
   Order<<MBEOrder_;
   Choices.push_back(Order.str());

   OptionNames.push_back("Fragmentation Method:");
   Choices.push_back(FOptions_.MethodName());

   OptionNames.push_back("Capping Method:");
   Choices.push_back(COptions_.MethodName());

   OptionNames.push_back("Embedding Method:");
   Choices.push_back(EOptions_.MethodName());

   OptionNames.push_back("BSSE Method:");
   Choices.push_back(BOptions_.MethodName());

   TableSpecs<std::string,std::string> Table(noptions);
   Table.SetTitles(Titles);
   Table.Init(&OptionNames[0],&Choices[0]);
   (*outfile)<<Table.Table();
   (*outfile)<<std::endl;
}




}}//End namespaces


