#ifndef _psi_src_lib_libmints_cartesianiter_h_
#define _psi_src_lib_libmints_cartesianiter_h_

namespace psi {

/** CartesianIter gives the ordering of the Cartesian functions
    that is used in PSI4. */
class CartesianIter
{
protected:
    int a_;
    int b_;
    int c_;
    int l_;
    int bfn_;

public:
    /// Initialize the iterator for the given angular momentum.
    CartesianIter(int l);
    ~CartesianIter();

    /// Start the iteration.
    virtual void start();
    /// Move to the next Catesian function.
    virtual void next();
    /// Returns nonzero if the iterator currently holds valid data.
    virtual operator int();

    /// Returns the number of Cartesian functions.
    int n() const { return ((l_>=0)?((((l_)+2)*((l_)+1))>>1):0); }
    /// Returns the x exponent
    int a() const { return a_; }
    /// Returns the y exponent
    int b() const { return b_; }
    /// Returns the z exponent
    int c() const { return c_; }
    /// Return the angular momentum
    int l() const { return l_; }
    /// Returns a() if i==0, b() if i==1, and c() if i==2.
    int l(int i) const { return i ? (i==1 ? b_ : c_) : a_; }
    /** Returns the number of the current basis function within the shell.
        This starts at 0 and sequentially increases as next() is called. */
    int bfn() { return bfn_; }
};

}

#endif // _psi_src_lib_libmints_cartesianiter_h_

