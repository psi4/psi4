/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#ifndef _PSI4_SRC_DCT_DFTENSOR_H_
#define _PSI4_SRC_DCT_DFTENSOR_H_

#include "psi4/libmints/matrix.h"

namespace psi {
namespace dct {

class DFTensor : public Matrix {
  protected:
      int nQ_;
      Dimension dim2_;
      Dimension dim3_;
  
  public:
      DFTensor() : Matrix(){}
      DFTensor(const std::string& name, const int Q, const Dimension& idx2, const Dimension& idx3);
      DFTensor(const std::string& name, const DFTensor& tensor);
      const int nQ() const { return nQ_; }
      const Dimension& idx2pi() const { return dim2_; }
      const Dimension& idx3pi() const { return dim3_; }
      /// A helper function needed by the constructor.
      static const Dimension setup(const Dimension& idx2, const Dimension& idx3);

      // Non-static functions operate on the DFTensor they're called on.
      void add_3idx_transpose_inplace();
      /// r(Q|pq) = \sum_Q B(P|rs) (rs|pq) without transpose flag.
      /// r(Q|pq) = \sum_Q B(P|rs) (pq|rs) with.
      void contract343(const DFTensor& b, dpdbuf4& G, bool transpose, double alpha, double beta);
      /// Transform b(Q|mu,nu) from SO basis to another basis with symmetry
      void three_idx_primary_transform_gemm(const DFTensor& three_idx, const Matrix& left, const Matrix& right,
                                            double alpha, double beta);

      // Static functions create a new DFTensor and operate on it.
      static DFTensor three_idx_primary_transform(const DFTensor& three_idx, const Matrix& left, const Matrix& right);
      /// r(Q|pq) = \sum_Q J(PQ) B(P|pq)
      static DFTensor contract233(const Matrix& J, const DFTensor& B);
      /// (Q) (p|q) -> (Q|pq)
      static DFTensor contract123(const Matrix& Q, const Matrix& G);

};

}  // namespace dct
}  // namespace psi

#endif  // Header guard

