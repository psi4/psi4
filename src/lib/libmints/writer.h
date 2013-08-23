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

#ifndef WRITER_H
#define WRITER_H

#include <boost/shared_ptr.hpp>
#include <libmints/vector.h>
#include <string>

namespace psi {

class Molecule;
class Matrix;
class BasisSet;
class Wavefunction;
class Options;

class GradientWriter
{
    boost::shared_ptr<Molecule> molecule_;
    const Matrix& gradient_;

public:
    GradientWriter(boost::shared_ptr<Molecule> mol, const Matrix& grad);

    void write(const std::string& filename);
};

class MoldenWriter
{
    boost::shared_ptr<Wavefunction> wavefunction_;

public:
    void write(const std::string &filename, boost::shared_ptr<Matrix> Ca, boost::shared_ptr<Matrix> Cb, boost::shared_ptr<Vector> Ea, boost::shared_ptr<Vector> Eb, boost::shared_ptr<Vector> OccA, boost::shared_ptr<Vector> OccB);
    MoldenWriter(boost::shared_ptr<Wavefunction> wavefunction);

};

class MOWriter
{
    boost::shared_ptr<Wavefunction> wavefunction_;
    Options & options_;
private:
    double * Ca_pointer, * eps;
    int * map, * sym, * occ, nmo, nso;
    void write_mos(Molecule & mol);
public:
    MOWriter(boost::shared_ptr<Wavefunction> wavefunction,Options&options);
    void write();
};

class NBOWriter
{
    boost::shared_ptr<Wavefunction> wavefunction_;

public:
    NBOWriter(boost::shared_ptr<Wavefunction> wavefunction);

    void write(const std::string &filename);
};

}

#endif // WRITER_H
