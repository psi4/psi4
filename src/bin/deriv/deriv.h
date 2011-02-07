/*
 * deriv.h
 *
 *  Created on: Feb 24, 2009
 *      Author: jturney
 */

#ifndef DERIV_H_
#define DERIV_H_

#include <vector>
#include <libmints/mints.h>

namespace psi { namespace deriv {

class Deriv
{
    const boost::shared_ptr<BasisSet>& basis_;
    int natom_;
    const boost::shared_ptr<MatrixFactory>& factory_;

    CdSalcList cdsalcs_;

    // Reference type
    reftype ref_;

    std::vector<SharedMatrix> dH_;
    std::vector<SharedMatrix> dS_;

#if 0
    // AO version of the code.
    std::vector<SharedSimpleMatrix> dH_;
    std::vector<SharedSimpleMatrix> dS_;
#endif

    // Results go here.
    SharedSimpleMatrix QdH_;
    SharedSimpleMatrix WdS_;
    SharedSimpleMatrix tb_;

public:
    Deriv(reftype ref, const boost::shared_ptr<MatrixFactory>& factory, const boost::shared_ptr<BasisSet>& basis);

    void compute(const SharedMatrix& Q, const SharedMatrix& G, const SharedMatrix& W);

#if 0
    // AO version of the code
    void compute(SharedSimpleMatrix& Q, SharedSimpleMatrix& G, SharedSimpleMatrix& W);
#endif

    const SharedSimpleMatrix& one_electron() {
        return QdH_;
    }

    const SharedSimpleMatrix& overlap() {
        return WdS_;
    }

    const SharedSimpleMatrix& two_body() {
        return tb_;
    }
};

}}

#endif /* DERIV_H_ */
