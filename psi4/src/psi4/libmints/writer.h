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

#ifndef WRITER_H
#define WRITER_H

#include "psi4/pragma.h"
#include <memory>
#include "psi4/libmints/vector.h"
#include <string>
#include "typedefs.h"

namespace psi {

class Molecule;
class Matrix;
class Wavefunction;

class PSI_API FCHKWriter {
   private:
    /*! \brief Extracts information from a wavefunction object, and writes it into a formatted FCHK
     * file.  */
    std::shared_ptr<Wavefunction> wavefunction_;
    FILE *chk_;
    void write_number(const char *label, int value);
    void write_number(const char *label, double value);
    void write_sym_matrix(const char *label, const SharedMatrix &mat);
    void write_matrix(const char *label, const SharedVector &mat);
    void write_matrix(const char *label, const SharedMatrix &mat);
    void write_matrix(const char *label, const std::vector<double> &mat);
    void write_matrix(const char *label, const std::vector<int> &mat);
    std::string postscf_density_label_;
    std::string spin_postscf_density_label_;

   public:
    SharedMatrix Dtot_ao;
    SharedMatrix Ca_ao;
    SharedMatrix Cb_ao;
    FCHKWriter(std::shared_ptr<Wavefunction> wavefunction);
    void write(const std::string &filename);
    void set_postscf_density_label(const std::string &label);
    const SharedMatrix SCF_Dtot() const { return Dtot_ao; }
};


class PSI_API MoldenWriter {
    std::shared_ptr<Wavefunction> wavefunction_;

   public:
    PSI_DEPRECATED(
        "Constructing an MoldenWriter and then calling write instead of using `wfn.write_molden(name)` "
        "is both buggy and deprecated, and as soon as 1.5 it will stop working")
    MoldenWriter(std::shared_ptr<Wavefunction> wavefunction);

    PSI_DEPRECATED(
        "Constructing an MoldenWriter and then calling write instead of using `wfn.write_molden(name)` "
        "is both buggy and deprecated, and as soon as 1.5 it will stop working")
    void write(const std::string &filename, std::shared_ptr<Matrix> Ca, std::shared_ptr<Matrix> Cb,
               std::shared_ptr<Vector> Ea, std::shared_ptr<Vector> Eb, std::shared_ptr<Vector> OccA,
               std::shared_ptr<Vector> OccB, bool dovirtual);
};


class PSI_API MOWriter {
    std::shared_ptr<Wavefunction> wavefunction_;
    bool restricted_;

   private:
    double *Ca_pointer, *eps;
    int *map, *sym, *occ, nmo, nso;
    void write_mos(Molecule &mol);

   public:
    MOWriter(std::shared_ptr<Wavefunction> wavefunction);
    void write();
};
}  // namespace psi

#endif  // WRITER_H
