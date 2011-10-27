#include <libmints/mints.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include "hamiltonian.h"
#include "jk.h"

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {

Hamiltonian::Hamiltonian(boost::shared_ptr<JK> jk) :
    jk_(jk)
{
}
Hamiltonian::~Hamiltonian()
{
}

RHamiltonian::RHamiltonian(boost::shared_ptr<JK> jk) :
    Hamiltonian(jk)
{
}
RHamiltonian::~RHamiltonian()
{
}
SharedMatrix RHamiltonian::explicit_hamiltonian()
{
    boost::shared_ptr<Vector> diag = diagonal();

    SharedMatrix H(new Matrix("Explicit Hamiltonian", diag->nirrep(), diag->dimpi(), diag->dimpi()));

    boost::shared_ptr<Vector> b(diag->clone());
    boost::shared_ptr<Vector> s(diag->clone());
    std::vector<boost::shared_ptr<Vector> > bb;
    std::vector<boost::shared_ptr<Vector> > ss;
    bb.push_back(b);
    ss.push_back(s);
    for (int h = 0; h < diag->nirrep(); h++) {
        for (int n = 0; n < diag->dimpi()[h]; n++) {
            b->zero();
            s->zero();
            b->set(h,n,1.0);
            product(bb,ss); 
            C_DCOPY(diag->dimpi()[h],s->pointer(h),1,H->pointer(h)[n],1);
        }
    }

    return H;
}

UHamiltonian::UHamiltonian(boost::shared_ptr<JK> jk) :
    Hamiltonian(jk)
{
}
UHamiltonian::~UHamiltonian()
{
}

MatrixRHamiltonian::MatrixRHamiltonian(SharedMatrix M) :
    RHamiltonian(boost::shared_ptr<JK>()), M_(M)
{
}
MatrixRHamiltonian::~MatrixRHamiltonian()
{
}
void MatrixRHamiltonian::print_header() const 
{
    if (print_) {
        fprintf(outfile, "  ==> MatrixRHamiltonian (by Rob Parrish) <== \n\n");
    }
}
boost::shared_ptr<Vector> MatrixRHamiltonian::diagonal()
{
    boost::shared_ptr<Vector> diag(new Vector("Matrix Diagonal", M_->nirrep(), M_->rowspi()));
    for (int h = 0; h < M_->nirrep(); ++h) {
        int n = M_->rowspi()[h];
        if (!n) continue;
        double** Mp = M_->pointer(h);
        double*  Dp = diag->pointer(h);
        for (int i = 0; i < n; ++i) {
            Dp[i] = Mp[i][i];  
        }
    }
    return diag;
}
void MatrixRHamiltonian::product(const std::vector<boost::shared_ptr<Vector> >& x, 
                              std::vector<boost::shared_ptr<Vector> >& b)
{
    for (int N = 0; N < x.size(); ++N) {
        for (int h = 0; h < M_->nirrep(); ++h) {
            int n = M_->rowspi()[h];
            if (!n) continue;
            double** Mp = M_->pointer(h);
            double*  xp = x[N]->pointer(h);    
            double*  bp = b[N]->pointer(h);    
            C_DGEMV('N',n,n,1.0,Mp[0],n,xp,1,0.0,bp,1);
        }
    }
}

MatrixUHamiltonian::MatrixUHamiltonian(std::pair<SharedMatrix, SharedMatrix > M) :
    UHamiltonian(boost::shared_ptr<JK>()), M_(M)
{
}
MatrixUHamiltonian::~MatrixUHamiltonian()
{
}
void MatrixUHamiltonian::print_header() const 
{
    if (print_) {
        fprintf(outfile, "  ==> MatrixUHamiltonian (by Rob Parrish) <== \n\n");
    }
}
std::pair<boost::shared_ptr<Vector>, boost::shared_ptr<Vector> > MatrixUHamiltonian::diagonal()
{
    boost::shared_ptr<Vector> diaga(new Vector("Alpha Matrix Diagonal", M_.first->nirrep(), M_.first->rowspi()));
    boost::shared_ptr<Vector> diagb(new Vector("Beta Matrix Diagonal", M_.first->nirrep(), M_.first->rowspi()));
    for (int h = 0; h < M_.first->nirrep(); ++h) {
        int n = M_.first->rowspi()[h];
        if (!n) continue;
        double** Map = M_.first->pointer(h);
        double*  Dap = diaga->pointer(h);
        double** Mbp = M_.second->pointer(h);
        double*  Dbp = diagb->pointer(h);
        for (int i = 0; i < n; ++i) {
            Dap[i] = Map[i][i];  
            Dbp[i] = Mbp[i][i];  
        }
    }
    return make_pair(diaga,diagb);
}
void MatrixUHamiltonian::product(const std::vector<std::pair<boost::shared_ptr<Vector>,boost::shared_ptr<Vector> > >& x, 
                              std::vector<std::pair<boost::shared_ptr<Vector>,boost::shared_ptr<Vector> > >& b)
{
    for (int N = 0; N < x.size(); ++N) {
        for (int h = 0; h < M_.first->nirrep(); ++h) {
            int n = M_.first->rowspi()[h];
            if (!n) continue;
            double** Map = M_.first->pointer(h);
            double*  xap = x[N].first->pointer(h);    
            double*  bap = b[N].first->pointer(h);    
            C_DGEMV('N',n,n,1.0,Map[0],n,xap,1,0.0,bap,1);
            double** Mbp = M_.second->pointer(h);
            double*  xbp = x[N].second->pointer(h);    
            double*  bbp = b[N].second->pointer(h);    
            C_DGEMV('N',n,n,1.0,Map[0],n,xap,1,0.0,bap,1);
        }
    }
}

CISRHamiltonian::CISRHamiltonian(boost::shared_ptr<JK> jk, 
                                 SharedMatrix Caocc, 
                                 SharedMatrix Cavir, 
                                 boost::shared_ptr<Vector> eps_aocc,
                                 boost::shared_ptr<Vector> eps_avir) : 
    RHamiltonian(jk), Caocc_(Caocc), Cavir_(Cavir), eps_aocc_(eps_aocc), eps_avir_(eps_avir), singlet_(true) 
{
}
CISRHamiltonian::~CISRHamiltonian()
{
}
void CISRHamiltonian::print_header() const 
{
    if (print_) {
        fprintf(outfile, "  ==> CISRHamiltonian (by Rob Parrish) <== \n\n");
    }
}
boost::shared_ptr<Vector> CISRHamiltonian::diagonal()
{
    int nirrep = eps_aocc_->nirrep();
    int* nov = new int[nirrep];
    ::memset((void*)nov,'\0',nirrep*sizeof(int));
    for (int symm = 0; symm < nirrep; ++symm) {
        for (int h = 0; h < nirrep; ++h) {
            nov[symm] += eps_aocc_->dimpi()[h] * eps_avir_->dimpi()[symm^h];
        }
    }

    boost::shared_ptr<Vector> diag(new Vector("CIS Diagonal", nirrep, nov));

    for (int symm = 0; symm < nirrep; ++symm) {
        long int offset = 0L;
        for (int h = 0; h < nirrep; ++h) {
            int nocc = eps_aocc_->dimpi()[h];
            int nvir = eps_avir_->dimpi()[symm^h];

            if (!nocc || !nvir) continue;
            
            double* eop = eps_aocc_->pointer(h); 
            double* evp = eps_avir_->pointer(symm^h); 
            double* dp  = diag->pointer(symm); 

            for (int i = 0; i < nocc; ++i) {
                for (int a = 0; a < nvir; ++a) {
                    dp[i * nvir + a + offset] = evp[a] - eop[i];
                }
            }
            offset += nocc * nvir;
        }
    }

    delete[] nov;
    return diag;
}
void CISRHamiltonian::product(const std::vector<boost::shared_ptr<Vector> >& x,
                                    std::vector<boost::shared_ptr<Vector> >& b)
{
    std::vector<SharedMatrix >& C_left = jk_->C_left();
    std::vector<SharedMatrix >& C_right = jk_->C_right();

    C_left.clear();
    C_right.clear();

    int nirrep = (x.size() ? x[0]->nirrep() : 0);

    for (int symm = 0; symm < nirrep; ++symm) {
        for (int N = 0; N < x.size(); ++N) {
            C_left.push_back(Caocc_);
            double* xp = x[N]->pointer(symm);

            std::stringstream ss;
            ss << "C_right, h = " << symm << ", N = " << N;
            SharedMatrix Cr(new Matrix(ss.str(), Caocc_->nirrep(), Caocc_->rowspi(), Caocc_->colspi(), symm));
            
            long int offset = 0L;
            for (int h = 0; h < Caocc_->nirrep(); ++h) {

                int nocc = Caocc_->colspi()[h];
                int nvir = Cavir_->colspi()[h^symm];
                int nso  = Cavir_->rowspi()[h^symm];

                if (!nso || !nocc || !nvir) continue;

                double** Cvp = Cavir_->pointer(h^symm);
                double** Crp = Cr->pointer(h^symm);

                C_DGEMM('N','T',nso,nocc,nvir,1.0,Cvp[0],nvir,&xp[offset],nvir,0.0,Crp[0],nocc);

                offset += nocc * nvir;
            }

            C_right.push_back(Cr); 
        }
    }

    jk_->compute();

    const std::vector<SharedMatrix >& J = jk_->J();
    const std::vector<SharedMatrix >& K = jk_->K();

    double* Tp = new double[Caocc_->max_nrow() * Caocc_->max_ncol()];

    for (int symm = 0; symm < nirrep; ++symm) {
        for (int N = 0; N < x.size(); ++N) {
    
            double* bp = b[N]->pointer(symm);
            double* xp = x[N]->pointer(symm);
            long int offset = 0L;

            for (int h = 0; h < Caocc_->nirrep(); ++h) {

                int nsoocc = Caocc_->rowspi()[h];
                int nocc = Caocc_->colspi()[h];
                int nsovir = Cavir_->rowspi()[h^symm];
                int nvir = Cavir_->colspi()[h^symm];
    
                if (!nsoocc || !nsovir || !nocc || !nvir) continue;
    
                double** Cop = Caocc_->pointer(h);
                double** Cvp = Cavir_->pointer(h^symm);
                double*  eop  = eps_aocc_->pointer(h);
                double*  evp  = eps_avir_->pointer(h^symm);

                // TODO
                double** Jp  = J[symm * x.size() + N]->pointer(h);
                double** Kp  = K[symm * x.size() + N]->pointer(h);
    
                // -(ij|ab)P_jb = C_im K_mn C_na
                C_DGEMM('T','N',nocc,nsovir,nsoocc,1.0,Cop[0],nocc,Kp[0],nsovir,0.0,Tp,nsovir);
                C_DGEMM('N','N',nocc,nvir,nsovir,-1.0,Tp,nsovir,Cvp[0],nvir,0.0,&bp[offset],nvir);
    
                //b[N]->print();

                if (singlet_) {
                    // 2(ia|jb)P_jb = C_im J_mn C_na
                    C_DGEMM('T','N',nocc,nsovir,nsoocc,1.0,Cop[0],nocc,Jp[0],nsovir,0.0,Tp,nsovir);
                    C_DGEMM('N','N',nocc,nvir,nsovir,2.0,Tp,nsovir,Cvp[0],nvir,1.0,&bp[offset],nvir);
                }
            
                //b[N]->print();
                
                for (int i = 0; i < nocc; ++i) {
                    for (int a = 0; a < nvir; ++a) {
                        bp[i * nvir + a + offset] += (evp[a] - eop[i]) * xp[i * nvir + a + offset];        
                    }
                }

                offset += nocc * nvir;
            }        
        }
    }

    delete[] Tp;

    if (debug_) {
        for (int N = 0; N < x.size(); N++) {
            x[N]->print();
            b[N]->print();
        }
    }
    
}
std::vector<SharedMatrix > CISRHamiltonian::unpack(boost::shared_ptr<Vector> eig)
{
    int nirrep = eig->nirrep();
    std::vector<SharedMatrix > t1;
    for (int symm = 0; symm < nirrep; ++symm) {
        SharedMatrix t(new Matrix("T",Caocc_->nirrep(), Caocc_->colspi(), Cavir_->colspi(), symm)); 
        long int offset = 0L;
        for (int h = 0; h < nirrep; ++h) {
            int nocc = Caocc_->colspi()[h];
            int nvir = Cavir_->colspi()[h^symm];
            
            if (!nocc || !nvir) continue;

            ::memcpy((void*)(t->pointer(h)[0]), (void*)(&eig->pointer(h)[offset]), sizeof(double) * nocc * nvir);    
        
            offset += nocc * nvir;
        }    

        t1.push_back(t);
    }

    return t1;
}
CPHFRHamiltonian::CPHFRHamiltonian(boost::shared_ptr<JK> jk, 
                                   SharedMatrix Caocc, 
                                   SharedMatrix Cavir, 
                                   boost::shared_ptr<Vector> eps_aocc,
                                   boost::shared_ptr<Vector> eps_avir) : 
    RHamiltonian(jk), Caocc_(Caocc), Cavir_(Cavir), eps_aocc_(eps_aocc), eps_avir_(eps_avir)
{
}
CPHFRHamiltonian::~CPHFRHamiltonian()
{
}
void CPHFRHamiltonian::print_header() const 
{
    if (print_) {
        fprintf(outfile, "  ==> CPHFRHamiltonian (by Rob Parrish) <== \n\n");
    }
}
boost::shared_ptr<Vector> CPHFRHamiltonian::diagonal()
{
    int nirrep = eps_aocc_->nirrep();
    int* nov = new int[nirrep];
    ::memset((void*)nov,'\0',nirrep*sizeof(int));
    for (int symm = 0; symm < nirrep; ++symm) {
        for (int h = 0; h < nirrep; ++h) {
            nov[symm] += eps_aocc_->dimpi()[h] * eps_avir_->dimpi()[symm^h];
        }
    }

    boost::shared_ptr<Vector> diag(new Vector("CPHF Diagonal", nirrep, nov));

    for (int symm = 0; symm < nirrep; ++symm) {
        long int offset = 0L;
        for (int h = 0; h < nirrep; ++h) {
            int nocc = eps_aocc_->dimpi()[h];
            int nvir = eps_avir_->dimpi()[symm^h];

            if (!nocc || !nvir) continue;
            
            double* eop = eps_aocc_->pointer(h); 
            double* evp = eps_avir_->pointer(symm^h); 
            double* dp  = diag->pointer(symm); 

            for (int i = 0; i < nocc; ++i) {
                for (int a = 0; a < nvir; ++a) {
                    dp[i * nvir + a + offset] = evp[a] - eop[i];
                }
            }
            offset += nocc * nvir;
        }
    }

    delete[] nov;
    return diag;
}
void CPHFRHamiltonian::product(const std::vector<boost::shared_ptr<Vector> >& x,
                                     std::vector<boost::shared_ptr<Vector> >& b)
{
    std::vector<SharedMatrix >& C_left = jk_->C_left();
    std::vector<SharedMatrix >& C_right = jk_->C_right();

    C_left.clear();
    C_right.clear();

    int nirrep = (x.size() ? x[0]->nirrep() : 0);

    for (int symm = 0; symm < nirrep; ++symm) {
        for (int N = 0; N < x.size(); ++N) {
            C_left.push_back(Caocc_);
            double* xp = x[N]->pointer(symm);

            std::stringstream ss;
            ss << "C_right, h = " << symm << ", N = " << N;
            SharedMatrix Cr(new Matrix(ss.str(), Caocc_->nirrep(), Caocc_->rowspi(), Caocc_->colspi(), symm));
            
            long int offset = 0L;
            for (int h = 0; h < Caocc_->nirrep(); ++h) {

                int nocc = Caocc_->colspi()[h];
                int nvir = Cavir_->colspi()[h^symm];
                int nso  = Cavir_->rowspi()[h^symm];

                if (!nso || !nocc || !nvir) continue;

                double** Cvp = Cavir_->pointer(h^symm);
                double** Crp = Cr->pointer(h^symm);

                C_DGEMM('N','T',nso,nocc,nvir,1.0,Cvp[0],nvir,&xp[offset],nvir,0.0,Crp[0],nocc);

                offset += nocc * nvir;
            }

            C_right.push_back(Cr); 
        }
    }
    
    jk_->compute();

    const std::vector<SharedMatrix >& J = jk_->J();
    const std::vector<SharedMatrix >& K = jk_->K();

    double* Tp = new double[Caocc_->max_nrow() * Caocc_->max_ncol()];

    for (int symm = 0; symm < nirrep; ++symm) {
        for (int N = 0; N < x.size(); ++N) {
    
            double* bp = b[N]->pointer(symm);
            double* xp = x[N]->pointer(symm);
            long int offset = 0L;

            for (int h = 0; h < Caocc_->nirrep(); ++h) {

                int nsoocc = Caocc_->rowspi()[h];
                int nocc = Caocc_->colspi()[h];
                int nsovir = Cavir_->rowspi()[h^symm];
                int nvir = Cavir_->colspi()[h^symm];
    
                if (!nsoocc || !nsovir || !nocc || !nvir) continue;
    
                double** Cop = Caocc_->pointer(h);
                double** Cvp = Cavir_->pointer(h^symm);
                double*  eop  = eps_aocc_->pointer(h);
                double*  evp  = eps_avir_->pointer(h^symm);
                double** Jp  = J[N]->pointer(h);
                double** Kp  = K[N]->pointer(h);
                double** K2p = K[N]->pointer(h^symm);
    
                // 4(ia|jb)P_jb = C_im J_mn C_na
                C_DGEMM('T','N',nocc,nsovir,nsoocc,1.0,Cop[0],nocc,Jp[0],nsovir,0.0,Tp,nsovir);
                C_DGEMM('N','N',nocc,nvir,nsovir,4.0,Tp,nsovir,Cvp[0],nvir,0.0,&bp[offset],nvir);
    
                // -(ib|ja)P_jb = C_in K_nm C_ma
                C_DGEMM('T','T',nocc,nsovir,nsoocc,1.0,Cop[0],nocc,K2p[0],nsoocc,0.0,Tp,nsovir);
                C_DGEMM('N','N',nocc,nvir,nsovir,-1.0,Tp,nsovir,Cvp[0],nvir,1.0,&bp[offset],nvir);
    
                // -(ij|ab)P_jb = C_im K_mn C_ra
                C_DGEMM('T','N',nocc,nsovir,nsoocc,1.0,Cop[0],nocc,Kp[0],nsovir,0.0,Tp,nsovir);
                C_DGEMM('N','N',nocc,nvir,nsovir,-1.0,Tp,nsovir,Cvp[0],nvir,1.0,&bp[offset],nvir);
    
                for (int i = 0; i < nocc; ++i) {
                    for (int a = 0; a < nvir; ++a) {
                        bp[i * nvir + a + offset] += (evp[a] - eop[i]) * xp[i * nvir + a + offset];        
                    }
                }

                offset += nocc * nvir;
            }        
        }
    }

    delete[] Tp;
}
std::vector<SharedMatrix > CPHFRHamiltonian::unpack(boost::shared_ptr<Vector> eig)
{
    int nirrep = eig->nirrep();
    std::vector<SharedMatrix > t1;
    for (int symm = 0; symm < nirrep; ++symm) {
        SharedMatrix t(new Matrix("T",Caocc_->nirrep(), Caocc_->colspi(), Cavir_->colspi(), symm)); 
        long int offset = 0L;
        for (int h = 0; h < nirrep; ++h) {
            int nocc = Caocc_->colspi()[h];
            int nvir = Cavir_->colspi()[h^symm];
            
            if (!nocc || !nvir) continue;

            ::memcpy((void*)(t->pointer(h)[0]), (void*)(&eig->pointer(h)[offset]), sizeof(double) * nocc * nvir);    
        
            offset += nocc * nvir;
        }    

        t1.push_back(t);
    }

    return t1;
}
 
}
