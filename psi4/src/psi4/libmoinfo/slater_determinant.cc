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

#include "slater_determinant.h"



namespace psi {

//SlaterDeterminant::SlaterDeterminant()
//{
//  startup();
//}

//SlaterDeterminant::SlaterDeterminant(SlaterDeterminant& det)
//{
//  alfa_sym    = det.alfa_sym;
//  beta_sym    = det.beta_sym;
//  alfa_string = det.alfa_string;
//  beta_string = det.beta_string;
//  alfa_bits   = det.alfa_bits;
//  beta_bits   = det.beta_bits;
//  startup();
//}

SlaterDeterminant::SlaterDeterminant(int alfa_sym_,int beta_sym_,std::vector<bool> alfa_bits_,std::vector<bool> beta_bits_)
: alfa_sym(alfa_sym_),beta_sym(beta_sym_),
  alfa_string(-1),beta_string(-1),
  alfa_bits(alfa_bits_), beta_bits(beta_bits_)
{
  startup();
}

SlaterDeterminant::~SlaterDeterminant()
{
  cleanup();
}

void SlaterDeterminant::startup()
{
}

void SlaterDeterminant::cleanup()
{
}

bool SlaterDeterminant::is_closed_shell()
{
  return(alfa_bits == beta_bits);
}

std::string SlaterDeterminant::get_label()
{
  std::string label;
  label = "|";
  int max_i = alfa_bits.size();
  for(int i = 0; i < max_i; ++i)  label += get_occupation_symbol(i);
  label += ">";
  return label;
}

char SlaterDeterminant::get_occupation_symbol(int i)
{
  char symbol;
  if( alfa_bits[i] &&  beta_bits[i]) symbol = '2';
  if( alfa_bits[i] && !beta_bits[i]) symbol = '+';
  if(!alfa_bits[i] &&  beta_bits[i]) symbol = '-';
  if(!alfa_bits[i] && !beta_bits[i]) symbol = '0';
  return(symbol);
}

}
