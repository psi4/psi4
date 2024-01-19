/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
#include <memory>
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
class Molecule;
class Wavefunction;
class Matrix;

class PSI_API Process {
   public:
    class PSI_API Environment {
        std::map<std::string, std::string> environment_;
        size_t memory_;
        int nthread_;
        std::string datadir_;

        /* Keep a shared_ptr to the default psio_manager
         *
         * If not, some destructors (like for Wavefunction) will call for the
         * default _psio_manager. However, the default psio manager is also global.
         * Destruction order is not specified for global variables, and therefore
         * the default psio manager can be destructed BEFORE destructors from within
         * wavefunction are called. By keeping a copy of the shared pointer here,
         * we can keep it alive as long as we need.
         *
         * HOWEVER, pay attention to the order. Order of destruction is reverse
         * of declaration within a class (guaranteed). So we put it early in the
         * class definition
         */
        std::shared_ptr<PSIOManager> _psio_manager_keepalive;

        std::shared_ptr<Molecule> molecule_;

       public:
        void initialize();

        /// Set active molecule
        void set_molecule(const std::shared_ptr<Molecule>& molecule);
        /// Return active molecule
        std::shared_ptr<Molecule> molecule() const;

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
}  // namespace psi

#endif /* PROCESS_H_ */
