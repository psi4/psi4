#include "sointegral.h"
#include "twobody.h"
#include "basisset.h"
#include "gshell.h"
#include "integral.h"

#include <boost/shared_ptr.hpp>

namespace psi {

TwoBodySOInt::TwoBodySOInt(const boost::shared_ptr<TwoBodyInt> & tb)
    : tb_(tb)
{
    // This code will be able to handle transforming into symmetry orbitals
    // with different basis sets. Currently other code does not know how to
    // generate the needed AO->SO transformation structures for this. But once
    // that is completed this should work as is.

    // Since I'm working on an SOBasis class and if I decide that is how I want
    // to do things these might need to change. The SOBasis logic might be moved
    // into BasisSet.
    b1_ = tb->basis1();
    b2_ = tb->basis2();
    b3_ = tb->basis3();
    b4_ = tb->basis4();

    // Allocate accumulation buffer
    buffer_ = new double[INT_NFUNC(b1_->has_puream(), b1_->max_am())
                        *INT_NFUNC(b2_->has_puream(), b2_->max_am())
                        *INT_NFUNC(b3_->has_puream(), b3_->max_am())
                        *INT_NFUNC(b4_->has_puream(), b4_->max_am())];
}

TwoBodySOInt::~TwoBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<BasisSet> TwoBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<BasisSet> TwoBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<BasisSet> TwoBodySOInt::basis2() const
{
    return b2_;
}

boost::shared_ptr<BasisSet> TwoBodySOInt::basis3() const
{
    return b3_;
}

boost::shared_ptr<BasisSet> TwoBodySOInt::basis4() const
{
    return b4_;
}

void TwoBodySOInt::compute_shell(int ish, int jsh, int ksh, int lsh)
{
    const double *aobuf = tb_->buffer();

//    SOTransformShell *t1 = b1_->so_transform(ish);
//    SOTransformShell *t2 = b2_->so_transform(jsh);
//    SOTransformShell *t3 = b3_->so_transform(ksh);
//    SOTransformShell *t4 = b4_->so_transform(lsh);

    int nso1 = b1_->shell(ish)->nfunction();
    int nso2 = b2_->shell(jsh)->nfunction();
    int nso3 = b3_->shell(ksh)->nfunction();
    int nso4 = b4_->shell(lsh)->nfunction();

    memset(buffer_, 0, nso1*nso2*nso3*nso4*sizeof(double));

    int nao2 = b2_->shell(jsh)->ncartesian();
    int nao3 = b3_->shell(ksh)->ncartesian();
    int nao4 = b4_->shell(lsh)->ncartesian();

//    for (int
}

}
