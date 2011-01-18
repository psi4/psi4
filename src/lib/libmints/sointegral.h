#ifndef _psi_src_lib_libmints_sointegral_h_
#define _psi_src_lib_libmints_sointegral_h_

namespace boost {
template <class T>
class shared_ptr;
}

namespace psi {

class TwoBodyInt;
class BasisSet;

class TwoBodySOInt
{
protected:
    boost::shared_ptr<TwoBodyInt> tb_;

    boost::shared_ptr<BasisSet> b1_;
    boost::shared_ptr<BasisSet> b2_;
    boost::shared_ptr<BasisSet> b3_;
    boost::shared_ptr<BasisSet> b4_;

    double *buffer_;

public:
    TwoBodySOInt(const boost::shared_ptr<TwoBodyInt>& );
    virtual ~TwoBodySOInt();

    boost::shared_ptr<BasisSet> basis() const;
    boost::shared_ptr<BasisSet> basis1() const;
    boost::shared_ptr<BasisSet> basis2() const;
    boost::shared_ptr<BasisSet> basis3() const;
    boost::shared_ptr<BasisSet> basis4() const;

    const double *buffer() const { return buffer_; }

    virtual void compute_shell(int, int, int, int);
};

}

#endif // _psi_src_lib_libmints_sointegral_h_
