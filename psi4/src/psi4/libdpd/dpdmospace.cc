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
#include <string>
#include <vector>

#include "dpdmospace.h"
using std::string;
using std::vector;
using std::cout;
namespace psi {

DPDMOSpace::DPDMOSpace(const char label, const string &indices, vector<int> orbspi)
{
  label_ = label;
  indices_ = dpd_split(indices);
  orbPI_ = orbspi;
  nIrrep_ = orbspi.size();
  nOrb_ = 0;
  for(int i=0; i < nIrrep_; i++)
    for(int j=0; j < orbPI_[i]; j++,nOrb_++)
      orbSym_.push_back(i);
}

DPDMOSpace::DPDMOSpace(const char label, const string &indices, Dimension orbspi)
{
  label_ = label;
  indices_ = dpd_split(indices);
  nIrrep_ = orbspi.n();
  for(int i=0; i < nIrrep_; i++) 
    orbPI_.push_back(orbspi[i]);
  nOrb_ = 0;
  for(int i=0; i < nIrrep_; i++)
    for(int j=0; j < orbPI_[i]; j++,nOrb_++)
      orbSym_.push_back(i);
}

DPDMOSpace::DPDMOSpace()
{

}

DPDMOSpace::~DPDMOSpace() { }

bool DPDMOSpace::operator==(const char *c)
{
  for(int iter=0; iter < this->indices_.size(); iter++)
    if(c == this->indices_[iter]) return true;

  return false;
}

bool DPDMOSpace::operator==(const string &c)
{
  for(int iter=0; iter < this->indices_.size(); iter++)
    if(c == this->indices_[iter]) return true;

  return false;
}

bool operator==(const string &c, const DPDMOSpace &rhs)
{
  for(int iter=0; iter < rhs.indices_.size(); iter++)
    if(c == rhs.indices_[iter]) return true;

  return false;
}

bool operator==(const char *c, const DPDMOSpace &rhs)
{
  for(int iter=0; iter < rhs.indices_.size(); iter++)
    if(c == rhs.indices_[iter]) return true;

  return false;
}

bool DPDMOSpace::operator==(DPDMOSpace &lhs)
{
  if(lhs.label() != label_) return false;
  if(lhs.nIrrep() != this->nIrrep_) return false;
  if(lhs.indices() != this->indices_) return false;
  if(lhs.nOrb() != this->nOrb_) return false;
  if(lhs.orbPI() != this->orbPI_) return false;
  if(lhs.orbSym() != this->orbSym_) return false;

  return true;
}

void DPDMOSpace::print()
{
  cout << "DPDMOSpace ID: " << label_ << "\n";
  cout << "DPDMOSpace nIrrep: " << nIrrep_ << "\n";
  cout << "DPDMOSpace nOrb: " << nOrb_ << "\n";
  cout << "DPDMOSpace indices: ";
  for(vector<string>::size_type iter = 0; iter != indices_.size(); iter++) cout << indices_[iter] << " ";
  cout << "\n";
  cout << "DPDMOSpace orbPI: ";
  for(vector<int>::size_type iter=0; iter != orbPI_.size(); iter++) cout << orbPI_[iter] << " ";
  cout << "\n";
  cout << "DPDMOSpace orbSym: ";
  for(vector<int>::size_type iter=0; iter != orbSym_.size(); iter++) cout << orbSym_[iter] << " ";
  cout << "\n";
}

// return vector of common indices between two DPDMOSpaces
vector<string> DPDMOSpace::overlap(DPDMOSpace &rhs)
{
  vector<string> v;
  for(int i=0; i < indices_.size(); i++)
    for(int j=0; j < rhs.indices_.size(); j++)
      if(indices_[i] == rhs.indices_[j]) v.push_back(indices_[i]);

  return v;
}

} // namespace psi
