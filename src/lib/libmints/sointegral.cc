#include "sointegral.h"
#include "twobody.h"
#include "basisset.h"
#include "gshell.h"
#include "integral.h"
#include "sobasis.h"
#include "matrix.h"

#include <boost/shared_ptr.hpp>

#define DEBUG

namespace psi {

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyInt> & ob,
                           const boost::shared_ptr<IntegralFactory>& integral)
    : ob_(ob)
{
    b1_ = boost::shared_ptr<SOBasis>(new SOBasis(ob->basis1(), integral));

    b1_->print();

    if (ob->basis2() == ob->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasis>(new SOBasis(ob->basis2(), integral));

    only_totally_symmetric_ = 0;

    buffer_ = new double[INT_NFUNC(ob->basis1()->has_puream(), ob->basis1()->max_am())
                        *INT_NFUNC(ob->basis2()->has_puream(), ob->basis2()->max_am())];
}

OneBodySOInt::~OneBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<SOBasis> OneBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> OneBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> OneBodySOInt::basis2() const
{
    return b2_;
}

void OneBodySOInt::compute_shell(int ish, int jsh)
{
    const double *aobuf = ob_->buffer();

    const SOTransform &t1 = b1_->trans(ish);
    const SOTransform &t2 = b2_->trans(jsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);

    memset(buffer_, 0, nso1*nso2*sizeof(double));

    int nao2 = b2_->naofunction(jsh);

    // I want to test only calling compute_shell for the first t1 and t2 aoshell pair
    // and then using the transformation coefficients to obtain everything else.
    // Otherwise using the petite list doesn't save us any computational time
    // in computing the integrals, but does save us time when we use the integrals.

    // loop through the AO shells that make up this SO shell
    for (int i=0; i<t1.naoshell; ++i) {
        const SOTransformShell &s1 = t1.aoshell[i];
        for (int j=0; j<t2.naoshell; ++j) {
            const SOTransformShell &s2 = t2.aoshell[j];

            fprintf(outfile, "aoshells: 1 = %d   2 = %d\n", s1.aoshell, s2.aoshell);
            ob_->compute_shell(s1.aoshell, s2.aoshell);

            for (int itr=0; itr<s1.nfunc; ++itr) {
                const SOTransformFunction &ifunc = s1.func[itr];
                double icoef = ifunc.coef;
                int iaofunc = ifunc.aofunc;
                int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                int iaooff = iaofunc;
                int isooff = isofunc;

                for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                    const SOTransformFunction &jfunc = s2.func[jtr];
                    double jcoef = jfunc.coef * icoef;
                    int jaofunc = jfunc.aofunc;
                    int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                    int jaooff = iaooff*nao2 + jaofunc;
                    int jsooff = isooff*nso2 + jsofunc;

                    buffer_[jsooff] += jcoef * aobuf[jaooff];

#ifdef DEBUG
//                    if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
                        fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
                                isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
                                ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
                                jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
//                    }
#endif
                }
            }
        }
    }
}

void OneBodySOInt::compute(boost::shared_ptr<Matrix> result)
{
    // Do not worry about zeroing out result
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over the unique AO shells.
    for (int ish=0; ish<ns1; ++ish) {
        for (int jsh=0; jsh<ns2; ++jsh) {

            fprintf(outfile, "computing ish = %d jsh = %d\n", ish, jsh);

            const SOTransform &t1 = b1_->trans(ish);
            const SOTransform &t2 = b2_->trans(jsh);

            int nso1 = b1_->nfunction(ish);
            int nso2 = b2_->nfunction(jsh);

            memset(buffer_, 0, nso1*nso2*sizeof(double));

            int nao2 = b2_->naofunction(jsh);

            // I want to test only calling compute_shell for the first t1 and t2 aoshell pair
            // and then using the transformation coefficients to obtain everything else.
            // Otherwise using the petite list doesn't save us any computational time
            // in computing the integrals, but does save us time when we use the integrals.

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    fprintf(outfile, "aoshells: 1 = %d   2 = %d\n", s1.aoshell, s2.aoshell);
                    ob_->compute_shell(s1.aoshell, s2.aoshell);

                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        int iirrep = b1_->irrep(ish, isooff);

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            int jirrep = b2_->irrep(jsh, jsooff);

                            buffer_[jsooff] += jcoef * aobuf[jaooff];

                            // Check the irreps to ensure symmetric quantities.
//                            if (fabs(buffer_[jsooff]) > 1.0e-10 && iirrep == jirrep) {
                                fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
                                        ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
                                        jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
//                                result->add(ifunc.irrep, b1_->function_within_irrep(ish, isofunc), b2_->function_within_irrep(jsh, jsofunc), jcoef * aobuf[jaooff]);
//                            }
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

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
