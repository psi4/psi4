/*
 * deriv.h
 *
 *  Created on: Feb 24, 2009
 *      Author: jturney
 */

#ifndef _psi_src_lib_libmints_deriv_h_
#define _psi_src_lib_libmints_deriv_h_

#include <vector>
#include "matrix.h"

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi {

class Wavefunction;
class IntegralFactory;
class BasisSet;
class SOBasisSet;
class Molecule;
class CdSalcList;

class Deriv
{
    const boost::shared_ptr<Wavefunction> wfn_;
    boost::shared_ptr<IntegralFactory> integral_;
    boost::shared_ptr<BasisSet> basis_;
    boost::shared_ptr<SOBasisSet> sobasis_;
    boost::shared_ptr<MatrixFactory> factory_;
    boost::shared_ptr<Molecule> molecule_;

    CdSalcList cdsalcs_;

    SharedMatrix P_2_;
    SharedMatrix W_2_;
    SharedMatrix SCF_D_;

    int natom_;

    std::vector<SharedMatrix> dH_;
    std::vector<SharedMatrix> dS_;

    // Results go here.
    /// One-electron contribution to the gradient
    SharedMatrix opdm_contr_;
    /// Reference one-electron contribution to the gradient
    SharedMatrix opdm_ref_contr_;
    /// Overlap contribution to the gradient
    SharedMatrix x_contr_;
    /// Reference overlap contribution to the gradient
    SharedMatrix x_ref_contr_;
    /// Two-electron contribution to the gradient
    SharedMatrix tpdm_contr_;
    /// Reference two-electron contribution to the gradient
    SharedMatrix tpdm_ref_contr_;
    /// Final gradient
    SharedMatrix gradient_;

    // Applies all symmetry operations to symmetrize the gradient and returns
    // the shared pointer to the input matrix, for convenient printing
    SharedMatrix symmetrize_gradient(SharedMatrix grad);
public:
    /*!
     * Constructor for derivative object.
     *
     * \param wave Wavefunction object to compute derivative for
     * \param needed_irreps by default do A1 derivatives
     * \param project_out_translations remove translations from the CdSalcs
     * \param project_out_rotations remove rotations from the CdSalcs
     */
    Deriv(const boost::shared_ptr<Wavefunction>& wave,
          char needed_irreps=0x1,
          bool project_out_translations=true,
          bool project_out_rotations=true);

    SharedMatrix compute();

    const SharedMatrix& one_electron() {
        return opdm_contr_;
    }

    const SharedMatrix& lagrangian() {
        return x_contr_;
    }

    const SharedMatrix& two_body() {
        return tpdm_contr_;
    }

    const SharedMatrix& gradient() {
        return gradient_;
    }
};

}

#endif /* DERIV_H_ */
