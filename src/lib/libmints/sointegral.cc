#include <psi4-dec.h>

#include "sointegral.h"
#include "twobody.h"
#include "basisset.h"
#include "gshell.h"
#include "integral.h"
#include "sobasis.h"
#include "matrix.h"
#include "molecule.h"

#include <boost/shared_ptr.hpp>

#define DEBUG

namespace psi {

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyAOInt> & ob,
                           const boost::shared_ptr<IntegralFactory>& integral,
                           int deriv)
    : ob_(ob), integral_(integral.get()), deriv_(ob->deriv())
{
    b1_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob->basis1(), integral));

    if (ob->basis2() == ob->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob->basis2(), integral));

    size_t size = INT_NCART(ob->basis1()->max_am())
                 *INT_NCART(ob->basis2()->max_am());
    if (deriv_ == 1) {
        // Make sure the AOInt knows how to compute derivatives
        if (!ob->has_deriv1())
            throw PSIEXCEPTION("OneBodySOInt::OneBodySOInt: The AO integral object doesn't provide first derivatives.");
        size *= 3 * ob->basis1()->molecule()->natom();
    }
    if (deriv_ > 1)
        throw FeatureNotImplemented("libmints", "Symmetrized integral derivatives greater than first order not implemented.",
                                    __FILE__, __LINE__);

    buffer_ = new double[size];
}

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyAOInt> & ob,
                           const IntegralFactory* integral,
                           int deriv)
    : ob_(ob), integral_(integral), deriv_(deriv)
{
    b1_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob->basis1(), integral));

    if (ob->basis2() == ob->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob->basis2(), integral));

    size_t size = INT_NCART(ob->basis1()->max_am())
                 *INT_NCART(ob->basis2()->max_am());
    if (deriv_ == 1)
        size *= 3 * ob->basis1()->molecule()->natom();
    if (deriv_ > 1)
        throw FeatureNotImplemented("libmints", "Symmetrized integral derivatives greater than first order not implemented.",
                                    __FILE__, __LINE__);

    buffer_ = new double[size];
}

OneBodySOInt::~OneBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<SOBasisSet> OneBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> OneBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> OneBodySOInt::basis2() const
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

            const SOTransform &t1 = b1_->trans(ish);
            const SOTransform &t2 = b2_->trans(jsh);

            int nso1 = b1_->nfunction(ish);
            int nso2 = b2_->nfunction(jsh);

            memset(buffer_, 0, nso1*nso2*sizeof(double));

            int nao2 = b2_->naofunction(jsh);

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    ob_->compute_shell(s1.aoshell, s2.aoshell);

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

void OneBodySOInt::compute_deriv1(std::vector<boost::shared_ptr<Matrix> > result,
                                  const CdSalcList &cdsalcs)
{
    // Do not worry about zeroing out result.

    // Do some checks:
    if (deriv_ < 1)
        throw SanityCheckError("OneBodySOInt::compute_deriv1: integral object not created to handle derivatives.", __FILE__, __LINE__);

    if (result.size() != cdsalcs.ncd())
        throw SanityCheckError("OneBodySOInt::compute_deriv1: result vector size does not match SALC size.", __FILE__, __LINE__);

    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();
    int i_offset = 0;
    int natom = ob_->basis()->molecule()->natom();
    int threenatom = 3*natom;

    // Loop over unique AO shells.
    for (int ish=0; ish<ns1; ++ish) {
        const SOTransform& t1 = b1_->trans(ish);
        int nso1 = b1_->nfunction(ish);
        int nao1 = b1_->naofunction(ish);

        for (int jsh=0; jsh<ns1; ++jsh) {
            const SOTransform& t2= b2_->trans(jsh);
            int nso2 = b2_->nfunction(jsh);
            int nao2 = b2_->naofunction(jsh);

            int nao12 = nao1 * nao2;
            int nso12 = nso1 * nso2;

            // Clear out the memory we need.
            memset(buffer_, 0, threenatom*nso1*nso2*sizeof(double));

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];

                int atom1 = ob_->basis1()->shell(s1.aoshell)->ncenter();

                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];

                    int atom2 = ob_->basis2()->shell(s2.aoshell)->ncenter();

                    ob_->compute_shell(s1.aoshell, s2.aoshell);

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

                            // atom 1
                            //    x
                            buffer_[jsooff + (atom1+0)*nso12] += jcoef * aobuf[jaooff + (atom1+0)*nao12];
                            //    y
                            buffer_[jsooff + (atom1+1)*nso12] += jcoef * aobuf[jaooff + (atom1+1)*nao12];
                            //    z
                            buffer_[jsooff + (atom1+2)*nso12] += jcoef * aobuf[jaooff + (atom1+2)*nao12];

                            // atom 2
                            //    x
                            buffer_[jsooff + (atom2+0)*nso12] += jcoef * aobuf[jaooff + (atom2+0)*nao12];
                            //    y
                            buffer_[jsooff + (atom2+1)*nso12] += jcoef * aobuf[jaooff + (atom2+1)*nao12];
                            //    z
                            buffer_[jsooff + (atom2+2)*nso12] += jcoef * aobuf[jaooff + (atom2+2)*nao12];
                        }
                    }
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

TwoBodySOInt::TwoBodySOInt(const boost::shared_ptr<TwoBodyAOInt> &tb,
                           const boost::shared_ptr<IntegralFactory>& integral)
    : tb_(tb), integral_(integral)
{
    // Try to reduce some work:
    b1_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis1(), integral));

    if (tb->basis1() == tb->basis2())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis2(), integral));

    if (tb->basis1() == tb->basis3())
        b3_ = b1_;
    else
        b3_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis3(), integral));

    if (tb->basis3() == tb->basis4())
        b4_ = b3_;
    else
        b4_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(tb->basis4(), integral));

    // Allocate accumulation buffer
    buffer_ = new double[16*INT_NCART(tb->basis1()->max_am())
            *INT_NCART(tb->basis2()->max_am())
            *INT_NCART(tb->basis3()->max_am())
            *INT_NCART(tb->basis4()->max_am())];

    ::memset(iirrepoff_, 0, sizeof(int) * 8);
    ::memset(jirrepoff_, 0, sizeof(int) * 8);
    ::memset(kirrepoff_, 0, sizeof(int) * 8);
    ::memset(lirrepoff_, 0, sizeof(int) * 8);

    for (int h=1; h<b1_->nirrep(); ++h) {
        iirrepoff_[h] = iirrepoff_[h-1] + b1_->nfunction_in_irrep(h-1);
        jirrepoff_[h] = jirrepoff_[h-1] + b2_->nfunction_in_irrep(h-1);
        kirrepoff_[h] = kirrepoff_[h-1] + b3_->nfunction_in_irrep(h-1);
        lirrepoff_[h] = lirrepoff_[h-1] + b4_->nfunction_in_irrep(h-1);
    }
}

TwoBodySOInt::~TwoBodySOInt()
{
    delete[] buffer_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis1() const
{
    return b1_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis2() const
{
    return b2_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis3() const
{
    return b3_;
}

boost::shared_ptr<SOBasisSet> TwoBodySOInt::basis4() const
{
    return b4_;
}

}
