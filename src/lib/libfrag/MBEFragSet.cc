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


#include "MBEFragSet.h"
#include "MBEFrag.h"
#include "Embedder.h"
#include "Fragmenter.h"
#include "Capper.h"
#include "BSSEer.h"

namespace psi{
namespace LibFrag{

void MBEFragSet::Copy(const MBEFragSet& other){
   this->EmbedFactory_=other.EmbedFactory_;
   this->BSSEFactory_=other.BSSEFactory_;
   this->CapFactory_=other.CapFactory_;
   this->FragFactory_=other.FragFactory_;
   this->Frags_=other.Frags_;
   this->Properties_=other.Properties_;
}

boost::shared_ptr<Embedder> MBEFragSet::EmbedFactory(){
   return EmbedFactory_;
}

void MBEFragSet::Embed(){
   EmbedFactory_->Embed(Frags_);
}

MBEFragSet::MBEFragSet(SharedOptions& Options,
      SharedMol& AMol){
   FragFactory_=Options->FOptions().MakeFactory(AMol);
   Properties_=FragFactory_->Fragment(AMol,Frags_);
   if(this->Severed()){
      CapFactory_=Options->COptions().MakeFactory(AMol);
      CapFactory_->MakeCaps(Frags_);
   }
   EmbedFactory_=Options->EOptions().MakeFactory(AMol);
   BSSEFactory_=Options->BOptions().MakeFactory(AMol);
   BSSEFactory_->AddBSSEJobs(Frags_);

}

MBEFragSet::MBEFragSet(const MBEFragSet& Monomers,const int N){
   Properties_=Monomers.Properties_;
   Monomers.FragFactory_->MakeNMers(Monomers.Frags_,N,Frags_,
         this->Disjoint());
   if(Properties_.Severed_)
      Monomers.CapFactory_->MakeCaps(Frags_);
   Monomers.EmbedFactory_->Embed(Frags_);
   Monomers.BSSEFactory_->AddBSSEJobs(Frags_);
}

}}//End namespaces

