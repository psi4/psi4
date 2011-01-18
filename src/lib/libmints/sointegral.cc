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

    //    b1_->print();

    if (ob->basis2() == ob->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasis>(new SOBasis(ob->basis2(), integral));

    only_totally_symmetric_ = 0;

    buffer_ = new double[INT_NCART(ob->basis1()->max_am())
            *INT_NCART(ob->basis2()->max_am())];
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

            //            fprintf(outfile, "aoshells: 1 = %d   2 = %d\n", s1.aoshell, s2.aoshell);
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
                    //                        fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
                    //                                isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
                    //                                ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
                    //                                jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
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

            //            fprintf(outfile, "computing ish = %d jsh = %d\n", ish, jsh);

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

                    //                    fprintf(outfile, "aoshells: 1 = %d   2 = %d\n", s1.aoshell, s2.aoshell);
                    ob_->compute_shell(s1.aoshell, s2.aoshell);
                    //                    for (int z=0; z < INT_NPURE(ob_->basis1()->shell(s1.aoshell)->am()) *
                    //                         INT_NPURE(ob_->basis2()->shell(s2.aoshell)->am()); ++z) {
                    //                        fprintf(outfile, "raw: %d -> %8.5f\n", z, aobuf[z]);
                    //                    }

                    for (int itr=0; itr<s1.nfunc; ++itr) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish, ifunc.irrep) + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        int iirrep = ifunc.irrep;

                        for (int jtr=0; jtr<s2.nfunc; ++jtr) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh, jfunc.irrep) + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            int jirrep = jfunc.irrep;

                            buffer_[jsooff] += jcoef * aobuf[jaooff];

                            //                            if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
                            //                                fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d\n",
                            //                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
                            //                                        ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
                            //                                        jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc));
                            //                                fprintf(outfile, "(%d|%d) += %8.5f * (%d|%d): %8.5f -> %8.5f\n",
                            //                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff]);
                            //                            }

                            // Check the irreps to ensure symmetric quantities.
                            if (ifunc.irrep == jfunc.irrep)
                                result->add(ifunc.irrep, b1_->function_within_irrep(ish, isofunc), b2_->function_within_irrep(jsh, jsofunc), jcoef * aobuf[jaooff]);
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

TwoBodySOInt::TwoBodySOInt(const boost::shared_ptr<TwoBodyInt> &tb,
                           const boost::shared_ptr<IntegralFactory>& integral)
    : tb_(tb)
{
    b1_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis1(), integral));
    b2_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis2(), integral));
    b3_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis3(), integral));
    b4_ = boost::shared_ptr<SOBasis>(new SOBasis(tb->basis4(), integral));

    // Allocate accumulation buffer
    buffer_ = new double[INT_NCART(tb->basis1()->max_am())
            *INT_NCART(tb->basis2()->max_am())
            *INT_NCART(tb->basis3()->max_am())
            *INT_NCART(tb->basis4()->max_am())];
}

TwoBodySOInt::~TwoBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis2() const
{
    return b2_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis3() const
{
    return b3_;
}

boost::shared_ptr<SOBasis> TwoBodySOInt::basis4() const
{
    return b4_;
}

void TwoBodySOInt::compute_shell(int ish, int jsh, int ksh, int lsh)
{
    fprintf(outfile, "computing shell (%d %d %d %d)\n", ish, jsh, ksh, lsh);

    const double *aobuff = tb_->buffer();

    const SOTransform &t1 = b1_->trans(ish);
    const SOTransform &t2 = b2_->trans(jsh);
    const SOTransform &t3 = b3_->trans(ksh);
    const SOTransform &t4 = b4_->trans(lsh);

    int nso1 = b1_->nfunction(ish);
    int nso2 = b2_->nfunction(jsh);
    int nso3 = b3_->nfunction(ksh);
    int nso4 = b4_->nfunction(lsh);

    int nao1 = b1_->naofunction(ish);
    int nao2 = b2_->naofunction(jsh);
    int nao3 = b3_->naofunction(ksh);
    int nao4 = b4_->naofunction(lsh);

    memset(buffer_, 0, nao1*nao2*nao3*nao4*sizeof(double));

    // loop through the ao shells that make up this so shell
    for (int i=0; i<t1.naoshell; i++) {
        const SOTransformShell &s1 = t1.aoshell[i];
        for (int j=0; j<t2.naoshell; j++) {
            const SOTransformShell &s2 = t2.aoshell[j];
            for (int k=0; k<t3.naoshell; k++) {
                const SOTransformShell &s3 = t3.aoshell[k];
                for (int l=0; l<t4.naoshell; l++) {
                    const SOTransformShell &s4 = t4.aoshell[l];
                    tb_->compute_shell(s1.aoshell, s2.aoshell, s3.aoshell, s4.aoshell);
                    for (int itr=0; itr<s1.nfunc; itr++) {
                        const SOTransformFunction &ifunc = s1.func[itr];
                        double icoef = ifunc.coef;
                        int iaofunc = ifunc.aofunc;
                        int isofunc = b1_->function_offset_within_shell(ish,
                                                                        ifunc.irrep)
                                + ifunc.sofunc;
                        int iaooff = iaofunc;
                        int isooff = isofunc;
                        for (int jtr=0; jtr<s2.nfunc; jtr++) {
                            const SOTransformFunction &jfunc = s2.func[jtr];
                            double jcoef = jfunc.coef * icoef;
                            int jaofunc = jfunc.aofunc;
                            int jsofunc = b2_->function_offset_within_shell(jsh,
                                                                            jfunc.irrep)
                                    + jfunc.sofunc;
                            int jaooff = iaooff*nao2 + jaofunc;
                            int jsooff = isooff*nso2 + jsofunc;
                            for (int ktr=0; ktr<s3.nfunc; ktr++) {
                                const SOTransformFunction &kfunc = s3.func[ktr];
                                double kcoef = kfunc.coef * jcoef;
                                int kaofunc = kfunc.aofunc;
                                int ksofunc = b3_->function_offset_within_shell(ksh,
                                                                                kfunc.irrep)
                                        + kfunc.sofunc;
                                int kaooff = jaooff*nao3 + kaofunc;
                                int ksooff = jsooff*nso3 + ksofunc;
                                for (int ltr=0; ltr<s4.nfunc; ltr++) {
                                    const SOTransformFunction &lfunc = s4.func[ltr];
                                    double lcoef = lfunc.coef * kcoef;
                                    int laofunc = lfunc.aofunc;
                                    int lsofunc = b4_->function_offset_within_shell(lsh,
                                                                                    lfunc.irrep)
                                            + lfunc.sofunc;
                                    int laooff = kaooff*nao4 + laofunc;
                                    int lsooff = ksooff*nso4 + lsofunc;
                                    buffer_[lsooff] += lcoef * aobuff[laooff];
                                    if (fabs(aobuff[laooff]*lcoef) > 1.0e-10) {
                                        fprintf(outfile, "(%d %d|%d %d) += %8.5lf * (%d %d|%d %d): %8.5lf -> %8.5lf\n",
                                                isofunc, jsofunc, ksofunc, lsofunc, lcoef,
                                                iaofunc, jaofunc, kaofunc, laofunc,
                                                aobuff[laooff], buffer_[lsooff]);

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

}
