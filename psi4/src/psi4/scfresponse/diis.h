/* 
 *  @BEGIN LICENSE
 * 
 *  Hilbert: a space for quantum chemistry plugins to Psi4 
 * 
 *  Copyright (c) 2020 by its authors (LICENSE).
 * 
 *  The copyrights for code used from other parties are included in
 *  the corresponding files.
 * 
 *  This program is free software: you can redistribute it and/or modify
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
 *  along with this program.  If not, see http://www.gnu.org/licenses/.
 * 
 *  @END LICENSE
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
