#ifndef _psi_src_lib_libmints_sointegral_h_
#define _psi_src_lib_libmints_sointegral_h_

#include "onebody.h"
#include "twobody.h"
#include "basisset.h"
#include "integral.h"
#include "sobasis.h"

namespace boost {
template <class T>
class shared_ptr;
}

namespace psi {

class Matrix;

class OneBodySOInt
{
protected:
    boost::shared_ptr<OneBodyInt> ob_;

    boost::shared_ptr<SOBasis> b1_;
    boost::shared_ptr<SOBasis> b2_;

    double *buffer_;

    int only_totally_symmetric_;
public:
    OneBodySOInt(const boost::shared_ptr<OneBodyInt>& , const boost::shared_ptr<IntegralFactory> &);
    virtual ~OneBodySOInt();

    boost::shared_ptr<SOBasis> basis() const;
    boost::shared_ptr<SOBasis> basis1() const;
    boost::shared_ptr<SOBasis> basis2() const;

    /**
     * Resulting integral buffer.
     */
    const double* buffer() const { return buffer_; }

    /**
     * Computes a one-electron integral matrix. Only works for symmetric operators
     * (multipole operators will not work).
     */
    void compute(boost::shared_ptr<Matrix> result);

    /**
     * Compute a given SO shell doublet placing the result into buffer_
     */
    virtual void compute_shell(int, int);

    /**
     * Currently these settings are ignored, placed in for future work.
     */
    int only_totally_symmetric() const { return only_totally_symmetric_; }
    void set_only_totally_symmetric(int i) { only_totally_symmetric_ = i; }
};

class TwoBodySOInt
{
protected:
    boost::shared_ptr<TwoBodyInt> tb_;

    boost::shared_ptr<SOBasis> b1_;
    boost::shared_ptr<SOBasis> b2_;
    boost::shared_ptr<SOBasis> b3_;
    boost::shared_ptr<SOBasis> b4_;

    double *buffer_;

public:
    TwoBodySOInt(const boost::shared_ptr<TwoBodyInt>&,
                 const boost::shared_ptr<IntegralFactory>&);
    virtual ~TwoBodySOInt();

    boost::shared_ptr<SOBasis> basis() const;
    boost::shared_ptr<SOBasis> basis1() const;
    boost::shared_ptr<SOBasis> basis2() const;
    boost::shared_ptr<SOBasis> basis3() const;
    boost::shared_ptr<SOBasis> basis4() const;

    const double *buffer() const { return buffer_; }

    virtual void compute_shell(int, int, int, int);
};

}

#endif // _psi_src_lib_libmints_sointegral_h_
