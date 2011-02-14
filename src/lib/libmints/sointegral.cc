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
                           const boost::shared_ptr<IntegralFactory>& integral)
    : ob_(ob), integral_(integral.get()), deriv_(ob->deriv())
{
    common_init();

}

OneBodySOInt::OneBodySOInt(const boost::shared_ptr<OneBodyAOInt> & ob,
                           const IntegralFactory* integral)
    : ob_(ob), integral_(integral), deriv_(ob->deriv())
{
    common_init();
}

OneBodySOInt::~OneBodySOInt()
{
    delete[] buffer_;
}

void OneBodySOInt::common_init()
{
    b1_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob_->basis1(), integral_));

    if (ob_->basis2() == ob_->basis1())
        b2_ = b1_;
    else
        b2_ = boost::shared_ptr<SOBasisSet>(new SOBasisSet(ob_->basis2(), integral_));

    int max1_nequivalent_atoms = ob_->basis1()->molecule()->max_nequivalent();
    int max2_nequivalent_atoms = ob_->basis2()->molecule()->max_nequivalent();

    size_ = max1_nequivalent_atoms*max2_nequivalent_atoms*
            INT_NCART(ob_->basis1()->max_am())*INT_NCART(ob_->basis2()->max_am());

    if (deriv_ == 1) {
        // Make sure the AOInt knows how to compute derivatives
        if (!ob_->has_deriv1())
            throw PSIEXCEPTION("OneBodySOInt::OneBodySOInt: The AO integral object doesn't provide first derivatives.");
        size_ *= 3 * ob_->basis1()->molecule()->natom();
    }
    if (deriv_ > 1)
        throw FeatureNotImplemented("libmints", "Symmetrized integral derivatives greater than first order not implemented.",
                                    __FILE__, __LINE__);

    // Grab the number of chunks from the ao object and enlarge the buffer to handle
    size_ *= ob_->nchunk();

    buffer_ = new double[size_];
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

            memset(buffer_, 0, size_*sizeof(double));

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
//                                fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d jaoff = %d jsooff = %d\n",
//                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff],
//                                        ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
//                                        jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc),
//                                        jaooff, jsooff);
//                                fprintf(outfile, "(%d|%d) += %8.5f * (%d|%d): %8.5f -> %8.5f\n",
//                                        isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff], buffer_[jsooff]);
//                                fflush(outfile);
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

void OneBodySOInt::compute(std::vector<boost::shared_ptr<Matrix> > results)
{
    // Do not worry about zeroing out result
    int nchunk = ob_->nchunk();
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
            int nso = nso1*nso2;

            // size_ includes nchunk
            memset(buffer_, 0, nchunk*nso*sizeof(double));

            int nao1 = b1_->naofunction(ish);
            int nao2 = b2_->naofunction(jsh);
            int nao = nao1*nao2;

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

                            // Handle chunks
                            for (int i=0; i<nchunk; ++i) {
                                double temp = jcoef * aobuf[jaooff + (i*nao)];

                                // If I get the if statements working below the next line is not needed at all.
                                buffer_[jsooff + (i*nso)] += temp;

                                int ijirrep = ifunc.irrep ^ jfunc.irrep;
                                if (ijirrep == results[i]->symmetry() && fabs(temp) > 1.0e-14) {

//                                    if (fabs(aobuf[jaooff]*jcoef) > 1.0e-10) {
//                                        fprintf(outfile, "(%2d|%2d) += %+6f * (%2d|%2d): %+6f -> %+6f iirrep = %d ifunc = %d, jirrep = %d jfunc = %d jaoff = %d jsooff = %d\n",
//                                                isofunc, jsofunc, jcoef, iaofunc, jaofunc, aobuf[jaooff + (i*nao)], buffer_[jsooff + (i*nso)],
//                                                ifunc.irrep, b1_->function_within_irrep(ish, isofunc),
//                                                jfunc.irrep, b2_->function_within_irrep(jsh, jsofunc),
//                                                jaooff + (i*nao), jsooff + (i*nso));
//                                        fflush(outfile);
//                                    }

                                    // Add the contribution to the matrix
                                    results[i]->add(ifunc.irrep,
                                                    b1_->function_within_irrep(ish, isofunc),
                                                    b2_->function_within_irrep(jsh, jsofunc),
                                                    temp);
                                }
                            }
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

    Molecule& mol = *ob_->basis1()->molecule().get();
    int ns1 = b1_->nshell();
    int ns2 = b2_->nshell();
    const double *aobuf = ob_->buffer();

    // Loop over unique SO shells.
    for (int ish=0; ish<ns1; ++ish) {
        const SOTransform& t1 = b1_->trans(ish);
        int nso1 = b1_->nfunction(ish);
        int nao1 = b1_->naofunction(ish);

//        fprintf(outfile, "nso1 = %d, nao1 = %d\n", nso1, nao1);

        for (int jsh=0; jsh<ns2; ++jsh) {
            const SOTransform& t2= b2_->trans(jsh);
            int nso2 = b2_->nfunction(jsh);
            int nao2 = b2_->naofunction(jsh);

            fprintf(outfile, "nso2 = %d, nao2 = %d\n", nso2, nao2);

            int nao12 = nao1 * nao2;
            int nso12 = nso1 * nso2;

            // Clear out the memory we need.
            memset(buffer_, 0, size_*sizeof(double));

            // loop through the AO shells that make up this SO shell
            // by the end of these 4 for loops we will have our final integral in buffer_
            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                int atom1 = 3*ob_->basis1()->shell(s1.aoshell)->ncenter();

                fprintf(outfile, "center1 = %d\n", atom1);
                for (int j=0; j<t2.naoshell; ++j) {
                    const SOTransformShell &s2 = t2.aoshell[j];
                    int atom2 = 3*ob_->basis2()->shell(s2.aoshell)->ncenter();

                    fprintf(outfile, "center2 = %d\n", atom2);
                    ob_->compute_shell_deriv1(s1.aoshell, s2.aoshell);

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

                            fprintf(outfile, "ish %d jsh %d\n\tx1aobuf %lf aopos %d sopos %d\n", ish, jsh,
                                    aobuf[jaooff + (atom1+0)*nao12], jaooff + (atom1+0)*nao12,
                                    jsooff + (atom1+0)*nso12);
                            fprintf(outfile, "\ty1aobuf %lf aopos %d sopos %d\n",
                                    aobuf[jaooff + (atom1+1)*nao12], jaooff + (atom1+1)*nao12,
                                    jsooff + (atom1+1)*nso12);
                            fprintf(outfile, "\tz1aobuf %lf aopos %d sopos %d\n",
                                    aobuf[jaooff + (atom1+2)*nao12], jaooff + (atom1+2)*nao12,
                                    jsooff + (atom1+2)*nso12);
                            fprintf(outfile, "\tx2aobuf %lf aopos %d sopos %d\n",
                                    aobuf[jaooff + (atom2+0)*nao12], jaooff + (atom2+0)*nao12,
                                    jsooff + (atom2+0)*nso12);
                            fprintf(outfile, "\ty2aobuf %lf aopos %d sopos %d\n",
                                    aobuf[jaooff + (atom2+1)*nao12], jaooff + (atom2+1)*nao12,
                                    jsooff + (atom2+1)*nso12);
                            fprintf(outfile, "\tz2aobuf %lf aopos %d sopos %d\n",
                                    aobuf[jaooff + (atom2+2)*nao12], jaooff + (atom2+2)*nao12,
                                    jsooff + (atom2+2)*nso12);
                            fprintf(outfile, "isooff %d jsooff %d iaooff %d jaooff %d\n", isooff, jsooff,
                                    iaooff, jaooff);
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


            // Ok, symmetrize the derivatives and add their contribution to the result matrix.
            // we'll go by components:

            fprintf(outfile, "ish %d jsh %d aoshell1 %d aoshell2 %d\n",
                    ish, jsh, t1.aoshell[0].aoshell, t2.aoshell[0].aoshell);
//            cdsalcs.print();
            fflush(outfile);

            for (int i=0; i<t1.naoshell; ++i) {
                const SOTransformShell &s1 = t1.aoshell[i];
                int atom1 = ob_->basis1()->shell(s1.aoshell)->ncenter();
                const CdSalcWRTAtom& cdsalc1 = cdsalcs.atom_salc(atom1);
                atom1 *= 3;

                for (int itr=0; itr<nso1; ++itr) {

                    int ifunc  = b1_->function(ish) + itr;
                    int iirrep = b1_->irrep(ifunc);
                    int irel   = b1_->function_within_irrep(ifunc);
                    int isooff = itr;

                    for (int jtr=0; jtr<nso2; ++jtr) {

                        int jfunc  = b2_->function(jsh) + jtr;
                        int jirrep = b2_->irrep(jfunc);
                        int jrel   = b2_->function_within_irrep(jfunc);
                        int jsooff = isooff*nso2 + jtr;

                        // atom 1
                        for (int nx=0; nx<cdsalc1.nx(); ++nx) {
                            const CdSalcWRTAtom::Component element = cdsalc1.x(nx);
                            double temp = element.coef * buffer_[jsooff + (atom1+0)*nso12];

                            fprintf(outfile, "x1 component temp %20.16f, coef %f buffer %f iirrep %d, jirrep %d, salc symmetry %d\n",
                                    temp, element.coef, buffer_[jsooff + (atom1+0)*nso12], iirrep, jirrep, element.irrep);
                            if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                fprintf(outfile, "\tadding iirrep %d, irel %d, jrel %d salc %d\n", iirrep, irel, jrel, element.salc);
                                result[element.salc]->add(iirrep, irel, jrel, temp);
                            }
                            fflush(outfile);
                        }
                        for (int ny=0; ny<cdsalc1.ny(); ++ny) {
                            const CdSalcWRTAtom::Component element = cdsalc1.y(ny);
                            double temp = element.coef * buffer_[jsooff + (atom1+1)*nso12];
                            fprintf(outfile, "y1 component temp %20.16f, coef %f buffer %f iirrep %d, jirrep %d, salc symmetry %d\n",
                                    temp, element.coef, buffer_[jsooff + (atom1+1)*nso12], iirrep, jirrep, element.irrep);
                            if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-14) {
                                fprintf(outfile, "\tadding iirrep %d, irel %d, jrel %d salc %d\n", iirrep, irel, jrel, element.salc);
                                result[element.salc]->add(iirrep, irel, jrel, temp);
                            }
                        }
                        for (int nz=0; nz<cdsalc1.nz(); ++nz) {
                            const CdSalcWRTAtom::Component element = cdsalc1.z(nz);
                            double temp = element.coef * buffer_[jsooff + (atom1+2)*nso12];
                            fprintf(outfile, "z1 component temp %20.16f, coef %f buffer %f iirrep %d, jirrep %d, salc symmetry %d\n",
                                    temp, element.coef, buffer_[jsooff + (atom1+2)*nso12], iirrep, jirrep, element.irrep);
                            if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-14) {
                                fprintf(outfile, "\tadding iirrep %d, irel %d, jrel %d salc %d\n", iirrep, irel, jrel, element.salc);
                                result[element.salc]->add(iirrep, irel, jrel, temp);
                            }
                        }
                    }
                }
            }

            // don't double add diagonals
            if (ish == jsh) continue;

            for (int j=0; j<t2.naoshell; ++j) {
                const SOTransformShell &s2 = t2.aoshell[j];
                int atom2 = ob_->basis2()->shell(s2.aoshell)->ncenter();
                const CdSalcWRTAtom& cdsalc2 = cdsalcs.atom_salc(atom2);
                atom2 *= 3;

                for (int itr=0; itr<nso1; ++itr) {

                    int ifunc  = b1_->function(ish) + itr;
                    int iirrep = b1_->irrep(ifunc);
                    int irel   = b1_->function_within_irrep(ifunc);
                    int isooff = itr;

                    for (int jtr=0; jtr<nso2; ++jtr) {

                        int jfunc  = b2_->function(jsh) + jtr;
                        int jirrep = b2_->irrep(jfunc);
                        int jrel   = b2_->function_within_irrep(jfunc);
                        int jsooff = isooff*nso2 + jtr;

                        // atom 2
                        for (int nx=0; nx<cdsalc2.nx(); ++nx) {
                            const CdSalcWRTAtom::Component element = cdsalc2.x(nx);
                            double temp = element.coef * buffer_[jsooff + (atom2+0)*nso12];

                            fprintf(outfile, "x2 component temp %20.16f, coef %f buffer %f iirrep %d, jirrep %d, salc symmetry %d\n",
                                    temp, element.coef, buffer_[jsooff + (atom2+0)*nso12], iirrep, jirrep, element.irrep);
                            if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-10) {
                                fprintf(outfile, "\tadding iirrep %d, irel %d, jrel %d salc %d\n", iirrep, irel, jrel, element.salc);
                                result[element.salc]->add(iirrep, irel, jrel, temp);
                            }
                            fflush(outfile);
                        }
                        for (int ny=0; ny<cdsalc2.ny(); ++ny) {
                            const CdSalcWRTAtom::Component element = cdsalc2.y(ny);
                            double temp = element.coef * buffer_[jsooff + (atom2+1)*nso12];
                            fprintf(outfile, "y2 component temp %20.16f, coef %f buffer %f iirrep %d, jirrep %d, salc symmetry %d\n",
                                    temp, element.coef, buffer_[jsooff + (atom2+1)*nso12], iirrep, jirrep, element.irrep);
                            if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-14) {
                                fprintf(outfile, "\tadding iirrep %d, irel %d, jrel %d salc %d\n", iirrep, irel, jrel, element.salc);
                                result[element.salc]->add(iirrep, irel, jrel, temp);
                            }
                        }
                        for (int nz=0; nz<cdsalc2.nz(); ++nz) {
                            const CdSalcWRTAtom::Component element = cdsalc2.z(nz);
                            double temp = element.coef * buffer_[jsooff + (atom2+2)*nso12];
                            fprintf(outfile, "z2 component temp %20.16f, coef %f buffer %f iirrep %d, jirrep %d, salc symmetry %d\n",
                                    temp, element.coef, buffer_[jsooff + (atom2+2)*nso12], iirrep, jirrep, element.irrep);
                            if ((iirrep ^ jirrep) == element.irrep && fabs(temp) > 1.0e-14) {
                                fprintf(outfile, "\tadding iirrep %d, irel %d, jrel %d salc %d\n", iirrep, irel, jrel, element.salc);
                                result[element.salc]->add(iirrep, irel, jrel, temp);
                            }
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

    int max1_nequivalent_atoms = tb->basis1()->molecule()->max_nequivalent();
    int max2_nequivalent_atoms = tb->basis2()->molecule()->max_nequivalent();
    int max3_nequivalent_atoms = tb->basis3()->molecule()->max_nequivalent();
    int max4_nequivalent_atoms = tb->basis4()->molecule()->max_nequivalent();

    int max_nequivalent_atom = max1_nequivalent_atoms * max2_nequivalent_atoms * max3_nequivalent_atoms * max4_nequivalent_atoms;

    // Allocate accumulation buffer
    size_ = max_nequivalent_atom*INT_NCART(tb->basis1()->max_am())
            *INT_NCART(tb->basis2()->max_am())
            *INT_NCART(tb->basis3()->max_am())
            *INT_NCART(tb->basis4()->max_am());
    buffer_ = new double[size_];

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
