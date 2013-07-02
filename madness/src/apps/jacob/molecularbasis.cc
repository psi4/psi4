/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/

#include <moldft/molecularbasis.h>


std::ostream& operator<<(std::ostream& s, const ContractedGaussianShell& c) {
    static const char* tag[] = {"s","p","d","f","g"};
    char buf[32768];
    char* p = buf;
    const std::vector<double>& coeff = c.get_coeff();
    const std::vector<double>& expnt = c.get_expnt();

    p += sprintf(p,"%s [",tag[c.angular_momentum()]);
    for (int i=0; i<c.nprim(); ++i) {
        p += sprintf(p, "%.6f(%.6f)",coeff[i],expnt[i]);
        if (i != (c.nprim()-1)) p += sprintf(p, ", ");
    }
    p += sprintf(p, "]");
    s << buf;
    return s;
}

std::ostream& operator<<(std::ostream& s, const AtomicBasis& c) {
    const std::vector<ContractedGaussianShell>& shells = c.get_shells();
    for (int i=0; i<c.nshell(); ++i) {
        s << "     " << shells[i] << std::endl;
    }
    if (c.has_guess_info()) {
        s << "     " << "Guess density matrix" << std::endl;
        s << c.get_dmat();
    }

    return s;
}

std::ostream& operator<<(std::ostream& s, const AtomicBasisFunction& a) {
    a.print_me(s);
    return s;
}

void AtomicBasisFunction::print_me(std::ostream& s) const {
    s << "atomic basis function: center " << xx << " " << yy << " " << zz << " : ibf " << ibf << " nbf " << nbf << " : shell " << shell << std::endl;
}

/// Print basis info for atoms in the molecule (once for each unique atom type)
void AtomicBasisSet::print(const Molecule& molecule) const {
    molecule.print();
    std::cout << "\n " << name << " atomic basis set" << std::endl;
    for (int i=0; i<molecule.natom(); ++i) {
        const Atom& atom = molecule.get_atom(i);
        const unsigned int atn = atom.atomic_number;
        for (int j=0; j<i; ++j) {
            if (molecule.get_atom(j).atomic_number == atn)
                goto doneitalready;
        }
        std::cout << std::endl;
        std::cout << "   " <<  get_atomic_data(atn).symbol << std::endl;
        std::cout << ag[atn];
doneitalready:
        ;
    }
}

/// Print basis info for all supported atoms
void AtomicBasisSet::print_all() const {
    std::cout << "\n " << name << " atomic basis set" << std::endl;
    for (unsigned int i=0; i<ag.size(); ++i) {
        if (ag[i].nbf() > 0) {
            std::cout << "   " <<  get_atomic_data(i).symbol << std::endl;
            std::cout << ag[i];
        }
    }
}

