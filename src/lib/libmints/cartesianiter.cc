#include "cartesianiter.h"

using namespace psi;

CartesianIter::CartesianIter(int l) :
    a_(0), b_(0), c_(0), l_(l), bfn_(0)
{

}

CartesianIter::~CartesianIter()
{

}

void CartesianIter::start()
{
    bfn_ = b_ = c_ = 0;
    a_ = l_;
}

void CartesianIter::next()
{
    if (c_ < l_ - a_) {
        b_--;
        c_++;
    }
    else {
        a_--;
        c_ = 0;
        b_ = l_ - a_;
    }
    bfn_++;
}

CartesianIter::operator int() 
{
    return (a_ >= 0);
}

