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

namespace psi{
class Molecule;
class Wavefunction;
class PointGroup;
class ExternalPotential;
class Matrix;
class Vector;

namespace efp {
    class EFP;
}

class Process
{
public:
    class Environment
    {
        std::map<std::string, std::string> environment_;
        unsigned long int memory_;
        int nthread_;

        std::shared_ptr<Molecule> molecule_;
        SharedMatrix gradient_;
        std::shared_ptr<efp::EFP> efp_;
        SharedMatrix efp_torque_;
        std::shared_ptr<Vector> frequencies_;
        std::shared_ptr<PointGroup> parent_symmetry_;

        std::shared_ptr<Molecule> legacy_molecule_;
        std::shared_ptr<Wavefunction> legacy_wavefunction_;
    public:
        void initialize();

        /// The symmetry of the molecule, before any displacements have been made
        std::shared_ptr<PointGroup> parent_symmetry() { return parent_symmetry_; }
        /// Set the "parent" symmetry
        void set_parent_symmetry(std::shared_ptr<PointGroup> pg) { parent_symmetry_ = pg; }
        const std::string& operator()(const std::string& key) const;
        std::string operator()(const std::string& key);
        const std::string set(const std::string& key, const std::string& value);

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

        /// Set frequencies manually
        void set_frequencies(const std::shared_ptr<Vector> f) { frequencies_ = f; }
        /// Get frequencies manually
        std::shared_ptr<Vector> frequencies() const { return frequencies_; }

        /// Set EFP
        void set_efp(const std::shared_ptr<psi::efp::EFP>& efp) { efp_ = efp; }
        /// Get EFP
        std::shared_ptr<psi::efp::EFP> get_efp() const { return efp_; }

        /// Set EFP gradient manually
        void set_efp_torque(const SharedMatrix g) { efp_torque_ = g; }
        /// Get EFP gradient manually
        SharedMatrix efp_torque() const { return efp_torque_; }

        /// Map containing current energies
        std::map<std::string, double> globals;

        /// Map containing current arrays
        std::map<std::string, std::shared_ptr<Matrix> > arrays;

        /// Number of threads per process
        int get_n_threads() const;
        void set_n_threads(int nthread);

        /// Memory in bytes
        unsigned long int get_memory() const;
        void set_memory(unsigned long int m);

        /// "Global" liboptions object.
        Options options;
    };

    static Environment environment;

    static Environment get_environment();
};
}



#endif /* PROCESS_H_ */
