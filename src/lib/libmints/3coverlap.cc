#include <stdexcept>
#include <libqt/qt.h>
#include "mints.h"

using namespace psi;

ThreeCenterOverlapInt::ThreeCenterOverlapInt(boost::shared_ptr<BasisSet> bs1,
                                             boost::shared_ptr<BasisSet> bs2,
                                             boost::shared_ptr<BasisSet> bs3)
    : bs1_(bs1), bs2_(bs2), bs3_(bs3)
{
    buffer_ = 0;

    // Allocate memory for buffer_ storage
    try {
        size_t size = INT_NCART(bs1->max_am()) * INT_NCART(bs2->max_am()) * INT_NCART(bs3->max_am());
        buffer_ = new double[size];
    }
    catch (std::bad_alloc& e) {
        fprintf(stderr, "Error allocating memory for buffer_\n");
        exit(EXIT_FAILURE);
    }
    memset(buffer_, 0, sizeof(double)*size);
}

ThreeCenterOverlapInt::~ThreeCenterOverlapInt()
{
    delete[] buffer_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis()
{
    return bs1_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis1()
{
    return bs1_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis2()
{
    return bs2_;
}

boost::shared_ptr<BasisSet> ThreeCenterOverlapInt::basis3()
{
    return bs3_;
}
