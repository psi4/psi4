#include<libmints/mints.h>
#include<libmints/sieve.h>
#include <libmints/sointegral_twobody.h>
#include<libpsio/psio.hpp>
#include"coupledpair.h"
#include"blas.h"
#ifdef _OPENMP
    #include<omp.h>
#endif

using namespace psi;
using namespace boost;

namespace psi{ namespace cepa{
  long int Position(long int i,long int j);
}}

namespace psi{namespace cepa{

void Vabcd_direct(boost::shared_ptr<BasisSet>primary,int nso,double*c2,double*r2,double*temp,int o){

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif
    const double **buffer1 = new const double*[nthreads];
    const double **buffer2 = new const double*[nthreads];

    boost::shared_ptr<IntegralFactory> integral (new IntegralFactory(primary,primary,primary,primary));
    boost::shared_ptr<TwoBodyAOInt> * eri1 = new boost::shared_ptr<TwoBodyAOInt>[nthreads];
    boost::shared_ptr<TwoBodyAOInt> * eri2 = new boost::shared_ptr<TwoBodyAOInt>[nthreads];
    for (int i = 0; i < nthreads; i++) {
        eri1[i] = boost::shared_ptr<TwoBodyAOInt>(integral->eri());
        eri2[i] = boost::shared_ptr<TwoBodyAOInt>(integral->eri());
        buffer1[i] = eri1[i]->buffer();
        buffer2[i] = eri2[i]->buffer();
    }
    int nshell = primary->nshell();

    long int otri = o*(o+1)/2;

    // pack c2+ into temp:
    for (int i = 0; i < o; i++) {
        for (int j = i; j < o; j++) {
            long int ij = Position(i,j);
            for (int a = 0; a < nso; a++) {
                for (int b = a+1; b < nso; b++) {
                    temp[Position(a,b)*otri+ij] =
                       c2[a*o*o*nso+b*o*o+i*o+j]+c2[b*o*o*nso+a*o*o+i*o+j];
                }
                temp[Position(a,a)*otri+ij] =
                   c2[a*o*o*nso+a*o*o+i*o+j];
            }
        }
    }
    // pack c2- into temp:
    long int shift = nso*(nso+1)/2 * otri;
    for (int i = 0; i < o; i++) {
        for (int j = i; j < o; j++) {
            long int ij = Position(i,j) + shift;
            for (int a = 0; a < nso; a++) {
                for (int b = a; b < nso; b++) {
                    temp[Position(a,b)*otri+ij] =
                       c2[a*o*o*nso+b*o*o+i*o+j]-c2[b*o*o*nso+a*o*o+i*o+j];
                }
            }
        }
    }

    // zero c2 - this is where the result will go
    memset((void*)c2,'\0',o*o*nso*nso*sizeof(double));

    for (int C = 0; C < primary->nshell(); C++) {
        int nc = primary->shell(C).nfunction();
        int cstart = primary->shell(C).function_index();

        for (int D = C; D < primary->nshell(); D++) {

            int nd = primary->shell(D).nfunction();
            int dstart = primary->shell(D).function_index();

            // thread over A and B shells
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int AB = 0; AB < nshell*nshell; AB++) {

                int B = AB % nshell;
                int A = (AB - B)/nshell;

                if ( A > B ) continue;

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                int na = primary->shell(A).nfunction();
                int astart = primary->shell(A).function_index();
                int nb = primary->shell(B).nfunction();
                int bstart = primary->shell(B).function_index();

                // (AC|BD)
                eri1[thread]->compute_shell(A,C,B,D);

                // (AD|BC)
                eri2[thread]->compute_shell(A,D,B,C);

                for (long int a = 0; a < na; a++) {

                    int bbegin = ( A == B ) ? a : 0;

                    for (long int b = bbegin; b < nb; b++) {

                        double*myr2_plus  = c2 + Position(a+astart,b+bstart)*otri;
                        double*myr2_minus = myr2_plus + shift;

                        for (long int c = 0; c < nc; c++) {

                            int dbegin = ( C == D ) ? c : 0;

                            for (long int d = dbegin; d < nd; d++) {

                                double val1 = buffer1[thread][a*nc*nb*nd+c*nb*nd+b*nd+d];
                                double val2 = buffer2[thread][a*nc*nb*nd+d*nb*nc+b*nc+c];

                                double valplus  = val1 + val2;
                                double valminus = val1 - val2;

                                double * myc2_plus = temp + Position(c+cstart,d+dstart)*otri;

                                if (fabs(valplus) > 1e-12){
                                   F_DAXPY(otri,valplus,myc2_plus,1,myr2_plus,1);
                                }
                                if (fabs(valminus) > 1e-12){
                                   double * myc2_minus = myc2_plus + shift;
                                   F_DAXPY(otri,valminus,myc2_minus,1,myr2_minus,1);
                                }
                            }
                        }
                    }
                }

            }
        }
    }

    // unpack r2+
    long int index;
    for (int a = 0, index = 0; a < nso; a++) {
        for (int b = 0; b < nso; b++) {
            long int ab = Position(a,b)*otri;
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    r2[index++] = 0.5*c2[ab+Position(i,j)];
                }
            }
        }
    }
    // unpack r2-
    for (int a = 0, index = 0; a < nso; a++) {
        for (int b = 0; b < nso; b++) {
            int sg2 = a < b ? 1 : -1;
            long int ab = Position(a,b)*otri + shift;
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    int sg = i < j ? sg2 : -sg2;
                    r2[index++] += 0.5*sg*c2[ab+Position(i,j)];
                }
            }
        }
    }
}

}}
