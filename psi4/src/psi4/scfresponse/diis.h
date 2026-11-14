/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef DIIS_H
#define DIIS_H

namespace psi{ namespace scf_response{

class DIIS {

public:

    DIIS(int n);
    ~DIIS();

    /// write Fock matrix disk
    void WriteVector(double * vector);
    void WriteVector(double * vector1, double * vector2);

    /// write erro vector to disk
    void WriteErrorVector(double * vector);
    void WriteErrorVector(double * vector1, double * vector2);

    /// perform diis extrapolation
    void Extrapolate(double * vector);
    void Extrapolate(double * vector1, double * vector2);

    /// restart diis
    void restart();

private:

    /// temporary storage
    double * tmp1_;
    double * tmp2_;

    /// determine diis expansion coefficients
    void DIISCoefficients(int nvec);

    /// maximum number of diis vectors
    int maxdiis_;          

    /// dimension of each diis vector
    double * diisvec_;      

    /// dimension of each diis vector
    int dimdiis_;           

    /// current number of diis vectors
    int diis_iter_;         

    /// diis vector to be replaced
    int replace_diis_iter_;
};

}} // end of namespace

#endif
