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


#ifndef DPD_MOSPACE_H
#define DPD_MOSPACE_H

#include <vector>
#include <string>
#include <libmints/mints.h>

using namespace std;

namespace psi {

vector<string> dpd_split(const string &indices);

class DPDMOSpace {
  protected:
    // name of the space
    char label_;
    // list of allowed orbital-index labels
    vector<string> indices_;
    // number of irreps
    int nIrrep_;
    // number of orbitals
    int nOrb_;
    // number of orbitals per irrep
    vector<int> orbPI_;
    // irrep of each orbital
    vector<int> orbSym_;

  public:
    DPDMOSpace(const char label, const string &indices, vector<int> orbspi);
    DPDMOSpace(const char label, const string &indices, Dimension orbspi);
    DPDMOSpace();
    ~DPDMOSpace();

    char label() { return label_; }
    vector<string> indices() { return indices_; }
    int nIrrep() { return nIrrep_; }
    int nOrb() { return nOrb_; }
    vector<int> orbPI() { return orbPI_; }
    vector<int> orbSym() { return orbSym_; }

    void print();
    bool operator==(const char *c);
    bool operator==(const string &c);
    bool operator==(DPDMOSpace &lhs);
    friend bool operator==(const char *c, const DPDMOSpace &rhs);
    friend bool operator==(const string &c, const DPDMOSpace &rhs);
    vector<string> overlap(DPDMOSpace &rhs);
};

} // namespace psi

#endif // DPD_MOSPACE_H

