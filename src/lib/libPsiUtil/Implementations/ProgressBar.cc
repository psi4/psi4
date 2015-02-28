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

#include "../ProgressBar.h"
#include "psi4-dec.h"
namespace psi{

ProgressBar::ProgressBar(const NTask_t NTasks):
      LessThan50_(NTasks<50),NTasks_(NTasks),Current_(0),
       NChars_(1),Char_('*'),NStars_(0){
   if(!LessThan50_){
      Remainder_=(NTasks%50);
      Increment_=(NTasks-Remainder_)/50;
   }
   else{
      Remainder_=50%NTasks;
      NChars_=floor(50/NTasks);
      Increment_=1;
   }
   (*outfile)<<"0%  10%  20%  30%  40%  50%  60%  70%  80%  90%  100%\n";
   (*outfile)<<"|----|----|----|----|----|----|----|----|----|----|\n";
   (*outfile)<<"*";
}

ProgressBar& ProgressBar::operator++(){
   ++Current_;
   if(LessThan50_){
      (*outfile)<<std::string(
            (NStars_>NTasks_-Remainder_-1?NChars_+1:NChars_),Char_);
      NStars_++;
   }
   else if((Current_==Increment_+1&&Remainder_>0)||
      (Current_==Increment_&&Remainder_==0)){
      (*outfile)<<std::string(NChars_,Char_);
      if(Remainder_>0)Remainder_--;
      Current_=0;
      NStars_++;
   }
   return *this;
}

}//End namespace

