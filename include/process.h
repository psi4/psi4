/*
 * JKFactory: Interface and code for highly parallel J and K
 *             builds.
 *
 *  Copyright (c) 2014 Ryan M. Richard
 *
 *  This file is part of JKFactory.
 *
 *  JKFactory is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef PROCESS_H_
#define PROCESS_H_

#include <string>
#include <map>
#include <liboptions/liboptions.h>
namespace psi{
class Molecule;
class Wavefunction;
class PointGroup;
class ExternalPotential;
class Matrix;
class Vector;
typedef boost::shared_ptr<Matrix> SharedMatrix;

class Process
{
public:
    class Environment
    {
        std::map<std::string, std::string> environment_;
        unsigned long int memory_;
        int nthread_;

        boost::shared_ptr<Molecule> molecule_;
        SharedMatrix gradient_;
        boost::shared_ptr<Vector> frequencies_;
        boost::shared_ptr<Wavefunction> wavefunction_;
        boost::shared_ptr<PointGroup> parent_symmetry_;
    public:
        void initialize();

        /// The symmetry of the molecule, before any displacements have been made
        boost::shared_ptr<PointGroup> parent_symmetry() { return parent_symmetry_; }
        /// Set the "parent" symmetry
        void set_parent_symmetry(boost::shared_ptr<PointGroup> pg) { parent_symmetry_ = pg; }
        const std::string& operator()(const std::string& key) const;
        std::string operator()(const std::string& key);
        const std::string& set(const std::string& key, const std::string& value);

        /// Set active molecule
        void set_molecule(const boost::shared_ptr<Molecule>& molecule);
        /// Return active molecule
        boost::shared_ptr<Molecule> molecule() const;

        /// Set wavefunction
        void set_wavefunction(const boost::shared_ptr<Wavefunction>& wavefunction);
        /// Get wavefunction
        boost::shared_ptr<Wavefunction> wavefunction() const;

        /// Set gradient manually
        void set_gradient(const SharedMatrix g) { gradient_ = g; }
        /// Get gradient manually
        SharedMatrix gradient() const { return gradient_; }

        /// Set frequencies manually
        void set_frequencies(const boost::shared_ptr<Vector> f) { frequencies_ = f; }
        /// Get frequencies manually
        boost::shared_ptr<Vector> frequencies() const { return frequencies_; }

        /// Map containing current energies
        std::map<std::string, double> globals;

        /// Number of threads per process
        int get_n_threads() const;
        void set_n_threads(int nthread);

        /// Memory in bytes
        unsigned long int get_memory() const;
        void set_memory(unsigned long int m);

        /// "Global" liboptions object.
        Options options;
    };

    class Arguments
    {
        std::vector<std::string> arguments_;
        int argc_;
        char **argv_;
    public:
        void initialize(int argc, char **argv);

        int argc() const;

        const std::string& operator()(int argc) const;
        std::string operator()(int argc);

        char** argv() const { return argv_; }
    };

    static Environment environment;
    static Arguments arguments;

    static Environment get_environment();
};
}



#endif /* PROCESS_H_ */
