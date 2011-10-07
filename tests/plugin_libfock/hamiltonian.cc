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

UHamiltonian::UHamiltonian(boost::shared_ptr<JK> jk) :
    Hamiltonian(jk)
{
}
UHamiltonian::~UHamiltonian()
{
}

MatrixRHamiltonian::MatrixRHamiltonian(boost::shared_ptr<Matrix> M) :
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

MatrixUHamiltonian::MatrixUHamiltonian(std::pair<boost::shared_ptr<Matrix>, boost::shared_ptr<Matrix> > M) :
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
                                 boost::shared_ptr<Matrix> Caocc, 
                                 boost::shared_ptr<Matrix> Cavir, 
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
    for (int h = 0; h < nirrep; ++h) {
        nov[h] = eps_aocc_->dimpi()[h] * eps_avir_->dimpi()[h];
    }

    boost::shared_ptr<Vector> diag(new Vector("CIS Diagonal", nirrep, nov));

    for (int h = 0; h < nirrep; ++h) {
        int nocc = eps_aocc_->dimpi()[h];
        int nvir = eps_avir_->dimpi()[h];

        if (!nocc || !nvir) continue;
        
        double* eop = eps_aocc_->pointer(h); 
        double* evp = eps_avir_->pointer(h); 
        double* dp  = diag->pointer(h); 

        for (int i = 0; i < nocc; ++i) {
            for (int a = 0; a < nvir; ++a) {
                dp[i * nvir + a] = evp[a] - eop[i];
            }
        }
    }

    delete[] nov;
    return diag;
}
void CISRHamiltonian::product(const std::vector<boost::shared_ptr<Vector> >& x,
                                    std::vector<boost::shared_ptr<Vector> >& b)
{
    if (debug_) {
        fprintf(outfile, "   > CISRHamiltonian::product <\n\n");
        Caocc_->print();
        Cavir_->print();
        for (int N = 0; N < x.size(); N++) {
            x[N]->print();
            b[N]->print();
        }
        fflush(outfile);
    }

    std::vector<boost::shared_ptr<Matrix> >& C_left = jk_->C_left();
    std::vector<boost::shared_ptr<Matrix> >& C_right = jk_->C_right();

    C_left.clear();
    C_right.clear();

    for (int N = 0; N < x.size(); ++N) {
        C_left.push_back(Caocc_);
        boost::shared_ptr<Matrix> Cr(Caocc_->clone());
        std::stringstream ss;
        ss << "C_right " << N;
        Cr->set_name(ss.str());
        
        for (int h = 0; h < Caocc_->nirrep(); ++h) {
            int nso = Caocc_->rowspi()[h];
            int nocc = Caocc_->colspi()[h];
            int nvir = Cavir_->colspi()[h];

            if (!nso || !nocc || !nvir) continue;

            double* xp = x[N]->pointer(h);
            double** Cvp = Cavir_->pointer(h);
            double** Crp = Cr->pointer(h);

            C_DGEMM('N','T',nso,nocc,nvir,1.0,Cvp[0],nvir,xp,nvir,0.0,Crp[0],nocc);
        }

        C_right.push_back(Cr); 
    }
    
    jk_->compute();

    const std::vector<boost::shared_ptr<Matrix> >& J = jk_->J();
    const std::vector<boost::shared_ptr<Matrix> >& K = jk_->K();

    double* Tp = new double[Caocc_->max_nrow() * Caocc_->max_ncol()];

    for (int N = 0; N < x.size(); ++N) {

        for (int h = 0; h < Caocc_->nirrep(); ++h) {
            int nso = Caocc_->rowspi()[h];
            int nocc = Caocc_->colspi()[h];
            int nvir = Cavir_->colspi()[h];

            if (!nso || !nocc || !nvir) continue;

            double** Cop = Caocc_->pointer(h);
            double** Cvp = Cavir_->pointer(h);
            double*  eop  = eps_aocc_->pointer(h);
            double*  evp  = eps_avir_->pointer(h);
            double** Jp  = J[N]->pointer(h);
            double** Kp  = K[N]->pointer(h);
            double*  bp  = b[N]->pointer(h);
            double*  xp  = x[N]->pointer(h);

            C_DGEMM('T','N',nocc,nso,nso,1.0,Cop[0],nocc,Kp[0],nso,0.0,Tp,nso);
            C_DGEMM('N','N',nocc,nvir,nso,-1.0,Tp,nso,Cvp[0],nvir,0.0,bp,nvir);
            if (singlet_) {
                C_DGEMM('T','N',nocc,nso,nso,1.0,Cop[0],nocc,Jp[0],nso,0.0,Tp,nso);
                C_DGEMM('N','N',nocc,nvir,nso,2.0,Tp,nso,Cvp[0],nvir,1.0,bp,nvir);
            }

            for (int i = 0; i < nocc; ++i) {
                for (int a = 0; a < nvir; ++a) {
                    bp[i * nvir + a] += (evp[a] - eop[i]) * xp[i * nvir + a];        
                }
            }
        }        
    }

    delete[] Tp;
}
 
CPHFRHamiltonian::CPHFRHamiltonian(boost::shared_ptr<JK> jk, 
                                   boost::shared_ptr<Matrix> Caocc, 
                                   boost::shared_ptr<Matrix> Cavir, 
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
    for (int h = 0; h < nirrep; ++h) {
        nov[h] = eps_aocc_->dimpi()[h] * eps_avir_->dimpi()[h];
    }

    boost::shared_ptr<Vector> diag(new Vector("CPHF Diagonal", nirrep, nov));

    for (int h = 0; h < nirrep; ++h) {
        int nocc = eps_aocc_->dimpi()[h];
        int nvir = eps_avir_->dimpi()[h];

        if (!nocc || !nvir) continue;
        
        double* eop = eps_aocc_->pointer(h); 
        double* evp = eps_avir_->pointer(h); 
        double* dp  = diag->pointer(h); 

        for (int i = 0; i < nocc; ++i) {
            for (int a = 0; a < nvir; ++a) {
                dp[i * nvir + a] = evp[a] - eop[i];
            }
        }
    }

    delete[] nov;
    return diag;
}
void CPHFRHamiltonian::product(const std::vector<boost::shared_ptr<Vector> >& x,
                                     std::vector<boost::shared_ptr<Vector> >& b)
{
    std::vector<boost::shared_ptr<Matrix> >& C_left = jk_->C_left();
    std::vector<boost::shared_ptr<Matrix> >& C_right = jk_->C_right();

    C_left.clear();
    C_right.clear();

    for (int N = 0; N < x.size(); ++N) {
        C_left.push_back(Caocc_);
        boost::shared_ptr<Matrix> Cr(Caocc_->clone());
        std::stringstream ss;
        ss << "C_right " << N;
        Cr->set_name(ss.str());
        
        for (int h = 0; h < Caocc_->nirrep(); ++h) {
            int nso = Caocc_->rowspi()[h];
            int nocc = Caocc_->colspi()[h];
            int nvir = Cavir_->colspi()[h];

            if (!nso || !nocc || !nvir) continue;

            double* xp = x[N]->pointer(h);
            double** Cvp = Cavir_->pointer(h);
            double** Crp = Cr->pointer(h);

            C_DGEMM('N','T',nso,nocc,nvir,1.0,Cvp[0],nvir,xp,nvir,0.0,Crp[0],nocc);
        }

        C_right.push_back(Cr); 
    }
    
    jk_->compute();

    const std::vector<boost::shared_ptr<Matrix> >& J = jk_->J();
    const std::vector<boost::shared_ptr<Matrix> >& K = jk_->K();

    double* Tp = new double[Caocc_->max_nrow() * Caocc_->max_ncol()];

    for (int N = 0; N < x.size(); ++N) {

        for (int h = 0; h < Caocc_->nirrep(); ++h) {
            int nso = Caocc_->rowspi()[h];
            int nocc = Caocc_->colspi()[h];
            int nvir = Cavir_->colspi()[h];

            if (!nso || !nocc || !nvir) continue;

            double** Cop = Caocc_->pointer(h);
            double** Cvp = Cavir_->pointer(h);
            double*  eop  = eps_aocc_->pointer(h);
            double*  evp  = eps_avir_->pointer(h);
            double** Jp  = J[N]->pointer(h);
            double** Kp  = K[N]->pointer(h);
            double*  bp  = b[N]->pointer(h);
            double*  xp  = x[N]->pointer(h);

            C_DGEMM('T','N',nocc,nso,nso,1.0,Cop[0],nocc,Jp[0],nso,0.0,Tp,nso);
            C_DGEMM('N','N',nocc,nvir,nso,4.0,Tp,nso,Cvp[0],nvir,0.0,bp,nvir);

            C_DGEMM('T','T',nocc,nso,nso,1.0,Cop[0],nocc,Kp[0],nso,0.0,Tp,nso);
            C_DGEMM('N','N',nocc,nvir,nso,-1.0,Tp,nso,Cvp[0],nvir,1.0,bp,nvir);

            C_DGEMM('T','N',nocc,nso,nso,1.0,Cop[0],nocc,Kp[0],nso,0.0,Tp,nso);
            C_DGEMM('N','N',nocc,nvir,nso,-1.0,Tp,nso,Cvp[0],nvir,1.0,bp,nvir);

            for (int i = 0; i < nocc; ++i) {
                for (int a = 0; a < nvir; ++a) {
                    bp[i * nvir + a] += (evp[a] - eop[i]) * xp[i * nvir + a];        
                }
            }
        }        
    }

    delete[] Tp;
}
 
}
