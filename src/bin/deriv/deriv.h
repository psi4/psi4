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
    std::vector<SharedSimpleMatrix> dH_;
    std::vector<SharedSimpleMatrix> dS_;

    boost::shared_ptr<BasisSet> basis_;

    int natom_;

    boost::shared_ptr<MatrixFactory> factory_;

    // Results go here.
    SharedSimpleMatrix QdH_;
    SharedSimpleMatrix WdS_;
    SharedSimpleMatrix tb_;

    // Sharederence
    reftype ref_;

public:
    Deriv(reftype ref, boost::shared_ptr<MatrixFactory>& factory, boost::shared_ptr<BasisSet>& basis);

    void compute(SharedSimpleMatrix& Q, SharedSimpleMatrix& G, SharedSimpleMatrix& W);

    SharedSimpleMatrix& one_electron() {
        return QdH_;
    }

    SharedSimpleMatrix& overlap() {
        return WdS_;
    }

    SharedSimpleMatrix& two_body() {
        return tb_;
    }
};

}}

#endif /* DERIV_H_ */
