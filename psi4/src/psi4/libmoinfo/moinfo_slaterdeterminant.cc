/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <iostream>
#include "moinfo.h"



using namespace std;

namespace psi {

MOInfo::SlaterDeterminant::SlaterDeterminant(const MOInfo *_moinfo)
    : moinfo(_moinfo)
{
}


MOInfo::SlaterDeterminant::~SlaterDeterminant()
{
}


/**
 * @fn MOInfo::SlaterDeterminant::is_closed_shell()
 */
bool MOInfo::SlaterDeterminant::is_closed_shell()
{
  for(int i=0;i<moinfo->get_nall();i++)
    if(bits[i]!=bits[i+moinfo->get_nall()])
      return(false);
  return(true);
}


/**
 * @fn MOInfo::SlaterDeterminant::is_spin_flip(SlaterDeterminant& det)
 */
bool MOInfo::SlaterDeterminant::is_spin_flipped(SlaterDeterminant& det)
{
  for(int i=0;i<moinfo->get_nall();i++){
    if(bits[i]!=det.test(i+moinfo->get_nall()))
      return(false);
    if(bits[i+moinfo->get_nall()]!=det.test(i))
      return(false);
  }
  return(true);
}

/**
 * @fn MOInfo::SlaterDeterminant::print(int n)
 */
std::string MOInfo::SlaterDeterminant::get_label()
{
  std::string label;
  label = "|";
  int counter = 0;
  for(int h = 0; h < moinfo->get_nirreps(); ++h){
    label += "[";
    for(int i = 0; i < moinfo->get_docc(h); ++i){
      label += get_occupation_symbol(counter);
      counter++;
    }
    for(int i = 0; i < moinfo->get_actv(h); ++i){
      label += get_occupation_symbol(counter);
      counter++;
    }
    counter += moinfo->get_extr(h);
    label += "]";
  }
  label += ">";
  return label;
}

/**
 * @fn MOInfo::SlaterDeterminant::get_internal_excitations(...)
 */
void MOInfo::SlaterDeterminant::get_internal_excitations(SlaterDeterminant& det,double& sign,
                                                 vector<pair<int,int> >& alpha_operators,
                                                 vector<pair<int,int> >& beta_operators)
{
  int ann, cre;
  int nall = moinfo->get_nall();
  bitdet bits_exc = det.get_bits();
  bitdet bits_tmp = bits;
  sign = 1.0;
  // Find one set of excitations at a time
  ann = -1; cre = -1;
  while(cre<nall){
    while(++ann<nall)
      if(bits[ann] && !bits_exc[ann]) break;
    while(++cre<nall)
      if(!bits[cre] && bits_exc[cre]) break;
    if(cre<nall){
      alpha_operators.push_back(make_pair(moinfo->get_all_to_occ(ann),moinfo->get_all_to_vir(cre)));
      sign *= annihilate(bits_tmp,ann);
      sign *=     create(bits_tmp,cre);
    }
  }
  ann = -1; cre = -1;
  while(cre<nall){
    while(++ann<nall)
      if(bits[ann+nall] && !bits_exc[ann+nall]) break;
    while(++cre<nall)
      if(!bits[cre+nall] && bits_exc[cre+nall]) break;
    if(cre<nall){
      beta_operators.push_back(make_pair(moinfo->get_all_to_occ(ann),moinfo->get_all_to_vir(cre)));
      sign *= annihilate(bits_tmp,ann+nall);
      sign *=     create(bits_tmp,cre+nall);
    }
  }
}

double MOInfo::SlaterDeterminant::annihilate(bitdet& bits_det,int so)
{
  if(bits_det.test(so)){
    bits_det.flip(so);
    double sign = 1.0;
    for(int i=0;i<so;i++)
      if(bits_det.test(i)) sign *= -1.0;
    return(sign);
  }
  else
    return(0.0);
}

double MOInfo::SlaterDeterminant::create(bitdet& bits_det,int so)
{
  if(!bits_det.test(so)){
    bits_det.flip(so);
    double sign = 1.0;
    for(int i=0;i<so;i++)
      if(bits_det.test(i)) sign *= -1.0;
    return(sign);
  }
  else
    return(0.0);
}

vector<int> MOInfo::SlaterDeterminant::get_aocc()
{
  vector<int> aocc;
  for(int i=0;i<moinfo->get_nall();i++)
    if(bits[i])
      aocc.push_back(moinfo->get_all_to_occ(i));
  return(aocc);
}

vector<int> MOInfo::SlaterDeterminant::get_bocc()
{
  vector<int> bocc;
  for(int i=0;i<moinfo->get_nall();i++)
    if(bits[i+moinfo->get_nall()])
      bocc.push_back(moinfo->get_all_to_occ(i));
  return(bocc);
}

vector<int> MOInfo::SlaterDeterminant::get_avir()
{
  vector<int> avir;
  for(int i=0;i<moinfo->get_nall();i++)
    if(!bits[i])
      avir.push_back(moinfo->get_all_to_vir(i));
  return(avir);
}

vector<int> MOInfo::SlaterDeterminant::get_bvir()
{
  vector<int> bvir;
  for(int i=0;i<moinfo->get_nall();i++)
    if(!bits[i+moinfo->get_nall()])
      bvir.push_back(moinfo->get_all_to_vir(i));
  return(bvir);
}

vector<bool> MOInfo::SlaterDeterminant::get_is_aocc()
{
  std::vector<int>  aocc = get_aocc();
  std::vector<bool> is_aocc(moinfo->get_nocc(),false);
  for(size_t i = 0; i < aocc.size(); i++) is_aocc[aocc[i]]=true;
  return(is_aocc);
}

vector<bool> MOInfo::SlaterDeterminant::get_is_bocc()
{
  std::vector<int>  bocc = get_bocc();
  std::vector<bool> is_bocc(moinfo->get_nocc(),false);
  for(size_t i = 0; i < bocc.size(); i++) is_bocc[bocc[i]]=true;
  return(is_bocc);
}

vector<bool> MOInfo::SlaterDeterminant::get_is_avir()
{
  std::vector<int>  avir = get_avir();
  std::vector<bool> is_avir(moinfo->get_nvir(),false);
  for(size_t i = 0; i < avir.size(); i++) is_avir[avir[i]]=true;
  return(is_avir);
}

vector<bool> MOInfo::SlaterDeterminant::get_is_bvir()
{
  std::vector<int>  bvir = get_bvir();
  std::vector<bool> is_bvir(moinfo->get_nvir(),false);
  for(size_t i = 0; i < bvir.size(); i++) is_bvir[bvir[i]]=true;
  return(is_bvir);
}

char MOInfo::SlaterDeterminant::get_occupation_symbol(int i)
{
  int occupation=bits[i]+bits[i+moinfo->get_nall()];
  if(occupation==0)                     return('0');
  if(occupation==2)                     return('2');
  if((occupation==1) && bits.test(i))   return('+');
  if((occupation==1) && bits.test(i+moinfo->get_nall())) return('-');
  return(' ');
}

}
