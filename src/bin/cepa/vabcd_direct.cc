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

    // schwarz sieve
    //double cutoff = 1e-4;
    //double cutoff2 = cutoff*cutoff;
    //boost::shared_ptr<ERISieve> sieve(new ERISieve(primary,cutoff));
    //const std::vector<std::pair<int,int> >& shell_pairs = sieve->shell_pairs();
    //int npairs = shell_pairs.size();

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

    long int nsotri = nso*(nso+1)/2;
    long int otri = o*(o+1)/2;

    // pack c2+ into r2:
    for (int i = 0; i < o; i++) {
        for (int j = i; j < o; j++) {
            long int ij = Position(i,j);
            for (int a = 0; a < nso; a++) {
                for (int b = a+1; b < nso; b++) {
                    r2[Position(a,b)*otri+ij] =
                       c2[a*o*o*nso+b*o*o+i*o+j]+c2[b*o*o*nso+a*o*o+i*o+j];
                }
                r2[Position(a,a)*otri+ij] =
                   c2[a*o*o*nso+a*o*o+i*o+j];
            }
        }
    }
    // pack c2- into r2:
    long int shift1 = nsotri * otri;
    for (int i = 0; i < o; i++) {
        for (int j = i; j < o; j++) {
            long int ij = Position(i,j) + shift1;
            for (int a = 0; a < nso; a++) {
                for (int b = a; b < nso; b++) {
                    r2[Position(a,b)*otri+ij] =
                       c2[a*o*o*nso+b*o*o+i*o+j]-c2[b*o*o*nso+a*o*o+i*o+j];
                }
            }
        }
    }

    // zero c2 - this is where the result will go
    memset((void*)c2,'\0',o*o*nso*nso*sizeof(double));

    long int shift2 = 0;
    double * c2_plus  = r2;
    double * c2_minus = r2 + shift1;

    for (int A = 0; A < primary->nshell(); A++) {

        int na = primary->shell(A).nfunction();
        int astart = primary->shell(A).function_index();

        for (int B = A; B < primary->nshell(); B++) {

            int nb = primary->shell(B).nfunction();
            int bstart = primary->shell(B).function_index();

            long int shift3 = na*nb*nsotri;
            double *temp_plus  = temp;
            double *temp_minus = temp + shift3;
            memset((void*)temp_plus,'\0',na*nb*nsotri*sizeof(double));
            memset((void*)temp_minus,'\0',na*nb*nsotri*sizeof(double));

            // thread over C and D shells, fill v(ab,cd)
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (int CD = 0; CD < nshell*nshell; CD++) {

                //s1 = omp_get_wtime();

                int D = CD % nshell;
                int C = (CD - D)/nshell;

                if ( C > D ) continue;

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                int nc = primary->shell(C).nfunction();
                int cstart = primary->shell(C).function_index();

                int nd = primary->shell(D).nfunction();
                int dstart = primary->shell(D).function_index();

                // (AC|BD)
                eri1[thread]->compute_shell(A,C,B,D);

                // (AD|BC)
                eri2[thread]->compute_shell(A,D,B,C);

                //ints[thread] += omp_get_wtime() - s1;
                //s1 = omp_get_wtime();

                for (long int a = 0; a < na; a++) {

                    int bbegin = ( A == B ) ? a : 0;

                    for (long int b = bbegin; b < nb; b++) {

                        for (long int c = 0; c < nc; c++) {

                            int dbegin = ( C == D ) ? c : 0;

                            for (long int d = dbegin; d < nd; d++) {

                                double val1 = buffer1[thread][a*nc*nb*nd+c*nb*nd+b*nd+d];
                                double val2 = buffer2[thread][a*nc*nb*nd+d*nb*nc+b*nc+c];

                                long int cd = Position(c+cstart,d+dstart);
                                temp_plus[(a*nb+b)*nsotri+cd]  = val1+val2;
                                temp_minus[(a*nb+b)*nsotri+cd] = val1-val2;
                            }
                        }
                    }
                }


            } // end of CD ... v should be filled

            // r(ab,ij) += t(cd,ij) v2(ab,cd)

            double * r2_plus  = c2 + shift2;
            double * r2_minus = c2 + nso*nso*otri+shift2;

            F_DGEMM('n','n',otri,na*nb,nsotri,1.0,c2_plus,otri,temp_plus,nsotri,0.0,r2_plus,otri);
            F_DGEMM('n','n',otri,na*nb,nsotri,1.0,c2_minus,otri,temp_minus,nsotri,0.0,r2_minus,otri);

            shift2 += na*nb*otri;
        }
    }

    // unpack into temp
    shift2 = 0;
    double *temp_plus  = temp;
    double *temp_minus = temp + nso*nso*otri;

    memset((void*)r2,'\0',o*o*nso*nso*sizeof(double));
    memset((void*)temp_plus,'\0',nso*nso*otri*sizeof(double));
    memset((void*)temp_minus,'\0',nso*nso*otri*sizeof(double));

    for (int A = 0; A < primary->nshell(); A++) {

        int na = primary->shell(A).nfunction();
        int astart = primary->shell(A).function_index();

        for (int B = A; B < primary->nshell(); B++) {

            int nb = primary->shell(B).nfunction();
            int bstart = primary->shell(B).function_index();

            double * r2_plus  = c2 + shift2;
            double * r2_minus = c2 + shift2 + nso*nso*otri;

            for (int a = 0; a < na; a++) {

                int absa = a + astart;

                for (int b = 0; b < nb; b++) {

                    int absb = b + bstart;

                    if (absa > absb) continue;

                    F_DCOPY(otri,r2_plus+ (a*nb+b)*otri,1,temp_plus+ (absa*nso+absb)*otri,1);
                    F_DCOPY(otri,r2_minus+(a*nb+b)*otri,1,temp_minus+(absa*nso+absb)*otri,1);
                }
            }
            shift2 += na*nb*otri;
        }
    }

    long int index;
    // unpack r2+
    for (int a = 0, index = 0; a < nso; a++) {
        for (int b = 0; b < nso; b++) {
            long int ab = a < b ? (a*nso+b)*otri : (b*nso+a)*otri;
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    r2[index++] = 0.5*temp_plus[ab+Position(i,j)];
                }
            }
        }
    }
    // unpack r2-
    for (int a = 0, index = 0; a < nso; a++) {
        for (int b = 0; b < nso; b++) {
            int sg2 = a < b ? 1 : -1;
            long int ab = a < b ? (a*nso+b)*otri : (b*nso+a)*otri;
            for (int i = 0; i < o; i++) {
                for (int j = 0; j < o; j++) {
                    int sg = i < j ? sg2 : -sg2;
                    r2[index++] += 0.5*sg*temp_minus[ab+Position(i,j)];
               }
            }
        }
    }

}

}}
