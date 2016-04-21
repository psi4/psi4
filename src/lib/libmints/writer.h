/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#ifndef WRITER_H
#define WRITER_H

#include <boost/shared_ptr.hpp>
#include <libmints/vector.h>
#include <string>
#include "typedefs.h"

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

class FCHKWriter
{
private:
    /*! \brief Extracts information from a wavefunction object, and writes it into a formatted FCHK file.  */
    boost::shared_ptr<Wavefunction> wavefunction_;
    FILE *chk_;
    void write_number(const char *label, int value);
    void write_number(const char *label, double value);
    void write_sym_matrix(const char *label, const SharedMatrix &mat);
    void write_matrix(const char *label, const SharedVector &mat);
    void write_matrix(const char *label, const SharedMatrix &mat);
    void write_matrix(const char *label, const std::vector<double> &mat);
    void write_matrix(const char *label, const std::vector<int> &mat);

public:
    FCHKWriter(boost::shared_ptr<Wavefunction> wavefunction);
    void write(const std::string &filename);
};

class MoldenWriter
{
    boost::shared_ptr<Wavefunction> wavefunction_;

public:
    void write(const std::string &filename, boost::shared_ptr<Matrix> Ca, boost::shared_ptr<Matrix> Cb, boost::shared_ptr<Vector> Ea, boost::shared_ptr<Vector> Eb, boost::shared_ptr<Vector> OccA, boost::shared_ptr<Vector> OccB);
    void writeNO(const std::string &filename, boost::shared_ptr<Matrix> Na, boost::shared_ptr<Matrix> Nb, boost::shared_ptr<Vector> Oa, boost::shared_ptr<Vector> Ob);
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