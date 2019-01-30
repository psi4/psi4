/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef PROCESS_H_
#define PROCESS_H_

#include <string>
#include <map>
#include "psi4/pragma.h"
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
class Molecule;
class Wavefunction;
class PointGroup;
class Matrix;
class Vector;

class PSI_API Process {
   public:
    class PSI_API Environment {
        std::map<std::string, std::string> environment_;
        size_t memory_;
        int nthread_;
        std::string datadir_;

        std::shared_ptr<Molecule> molecule_;
        SharedMatrix gradient_;
        std::shared_ptr<PointGroup> parent_symmetry_;

        std::shared_ptr<Molecule> legacy_molecule_;
        std::shared_ptr<Wavefunction> legacy_wavefunction_;

       public:
        void initialize();

        /// The symmetry of the molecule, before any displacements have been made
        std::shared_ptr<PointGroup> parent_symmetry() { return parent_symmetry_; }
        /// Set the "parent" symmetry
        void set_parent_symmetry(std::shared_ptr<PointGroup> pg) { parent_symmetry_ = pg; }

        /// Set active molecule
        void set_molecule(const std::shared_ptr<Molecule>& molecule);
        /// Return active molecule
        std::shared_ptr<Molecule> molecule() const;

        /// Temporary slots for legacy as a stop-gap
        /// Set active molecule
        void set_legacy_molecule(const std::shared_ptr<Molecule>& molecule);
        /// Return active molecule
        std::shared_ptr<Molecule> legacy_molecule() const;

        /// Set wavefunction
        void set_legacy_wavefunction(const std::shared_ptr<Wavefunction>& wavefunction);
        /// Get wavefunction
        std::shared_ptr<Wavefunction> legacy_wavefunction() const;

        /// Set gradient manually
        void set_gradient(const SharedMatrix g) { gradient_ = g; }
        /// Get gradient manually
        SharedMatrix gradient() const { return gradient_; }

        /// Map containing current energies
        std::map<std::string, double> globals;

        /// Map containing current arrays
        std::map<std::string, std::shared_ptr<Matrix> > arrays;

        /// Number of threads per process
        int get_n_threads() const;
        void set_n_threads(int nthread);

        /// Memory in bytes
        size_t get_memory() const;
        void set_memory(size_t m);

        /// PSIDATADIR
        std::string get_datadir() const { return datadir_; }
        void set_datadir(const std::string pdd) { datadir_ = pdd; }

        /// "Global" liboptions object.
        Options options;
    };

    static Environment environment;

    static Environment get_environment();
};
}

#endif /* PROCESS_H_ */
