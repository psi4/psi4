#include "mp2.h"
#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libqt/qt.h>
#include <psiconfig.h>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef HAVE_MKL
#include <mkl.h>
#endif

using namespace boost;
using namespace psi;

namespace psi { namespace dfcc {

MP2::MP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
  print_header();
  CC::print_header();

  common_init();
}
MP2::~MP2()
{
}
void MP2::common_init()
{
    mp2_algorithm_ = options_.get_str("MP2_ALGORITHM");
}
double MP2::compute_energy()
{
  energies_["MP2J Energy"] = 0.0;
  energies_["MP2K Energy"] = 0.0;

  if (mp2_algorithm_ == "MP2") {
    compute_MP2();
  } else if (mp2_algorithm_ == "DF") {
    compute_DF_MP2();
  } else if (mp2_algorithm_ == "PS") {
    compute_PS_MP2();
  } else if (mp2_algorithm_ == "PS1") {
    compute_DF_MP2J();
    compute_PS_DF_MP2K();
  } else if (mp2_algorithm_ == "PS2") {
    compute_DF_MP2J();
    compute_PS_MP2K();
  } else if (mp2_algorithm_ == "PS3") {
    compute_PS_MP2J();
    compute_PS_DF_MP2K();
  } else if (mp2_algorithm_ == "PS4") {
    compute_PS_MP2J();
    compute_PS_MP2K();
  } else if (mp2_algorithm_ == "TEST_DENOM") {
    test_denominators();
  } else if (mp2_algorithm_ == "TEST_DF") {
    test_df();
  } else if (mp2_algorithm_ == "TEST_PS") {
    test_ps();
  } else if (mp2_algorithm_ == "TEST_PS_OMEGA") {
    test_ps_omega();
  } else if (mp2_algorithm_ == "TEST_DPS_OMEGA") {
    test_dps_omega();
  }

  energies_["Opposite-Spin Energy"] = 0.5*energies_["MP2J Energy"];
  energies_["Same-Spin Energy"] = 0.5*energies_["MP2J Energy"] +  energies_["MP2K Energy"];
  energies_["Correlation Energy"] = energies_["MP2J Energy"] + energies_["MP2K Energy"];
  energies_["Total Energy"] = energies_["Reference Energy"] + energies_["Correlation Energy"];

  energies_["SCS Opposite-Spin Energy"] = 0.5*oss_*energies_["MP2J Energy"];
  energies_["SCS Same-Spin Energy"] = 0.5*sss_*energies_["MP2J Energy"] +  sss_*energies_["MP2K Energy"];
  energies_["SCS Correlation Energy"] = energies_["SCS Opposite-Spin Energy"] + energies_["SCS Same-Spin Energy"];
  energies_["SCS Total Energy"] = energies_["Reference Energy"] + energies_["SCS Correlation Energy"];

  print_energies();
  return energies_["Total Energy"];
}
void MP2::check_integrals(SharedMatrix E)
{
    fprintf(outfile, "  => Complete MO ERI Analysis <=\n\n");
    fprintf(outfile, "   -o : Active Occupied\n");   
    fprintf(outfile, "   -v : Active Virtual\n\n");   

    char type[4];
 
    double current; 
    int index;
    int bits[4];
    double sum[16];
    double max[16];
    ULI    n[16];

    ::memset(static_cast<void*>(sum), '\0', 16 * sizeof(double));
    ::memset(static_cast<void*>(max), '\0', 16 * sizeof(double));
    ::memset(static_cast<void*>(n), '\0', 16 * sizeof(ULI));

    double** Ip = E->pointer();    
    
    bits[0] = 0; 
    for (int a = nfocc_; a < nocc_; a++) {
        bits[1] = 0; 
        for (int b = nfocc_; b < nocc_; b++) {
            bits[2] = 0; 
            for (int c = nfocc_; c < nocc_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
            bits[2] = 1; 
            for (int c = nocc_; c < nocc_ + navir_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
        } 
        bits[1] = 1; 
        for (int b = nocc_; b < nocc_ + navir_; b++) {
            bits[2] = 0; 
            for (int c = nfocc_; c < nocc_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
            bits[2] = 1; 
            for (int c = nocc_; c < nocc_ + navir_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
        } 
    } 
    bits[0] = 1; 
    for (int a = nocc_; a < nocc_ + navir_; a++) {
        bits[1] = 0; 
        for (int b = nfocc_; b < nocc_; b++) {
            bits[2] = 0; 
            for (int c = nfocc_; c < nocc_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
            bits[2] = 1; 
            for (int c = nocc_; c < nocc_ + navir_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
        } 
        bits[1] = 1; 
        for (int b = nocc_; b < nocc_ + navir_; b++) {
            bits[2] = 0; 
            for (int c = nfocc_; c < nocc_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
            bits[2] = 1; 
            for (int c = nocc_; c < nocc_ + navir_; c++) {
                bits[3] = 0;
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nfocc_; d < nocc_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
                bits[3] = 1; 
                index = 8 * bits[0] + 4 * bits[1] + 2 * bits[2] + 1 * bits[3]; 
                for (int d = nocc_; d < nocc_ + navir_; d++) {
                    current = Ip[a * nmo_ + b][c * nmo_ + d];
                    if (max[index] < fabs(current))
                        max[index] = fabs(current);
                    sum[index] += current * current;
                    n[index]++; 
                } 
            } 
        } 
    } 

    for (int i = 0; i < 16; i++) {
        sum[i] = sqrt(sum[i] / (double) n[i]);
    } 

    fprintf(outfile, "     Type %14s %14s %14s\n", "Number", "Max Error" , "RMS Error"); 
    for (int i = 0; i < 16; i++) {
        type[0] = (i & 8 ? 'v' : 'o');
        type[1] = (i & 4 ? 'v' : 'o');
        type[2] = (i & 2 ? 'v' : 'o');
        type[3] = (i & 1 ? 'v' : 'o');
        fprintf(outfile, "  (%c%c|%c%c) %14ld %14.7E %14.7E\n", type[0], type[1], type[2], type[3], n[i], max[i], sum[i]); 
    }

}
void MP2::test_denominators()
{
    fprintf(outfile, "  ==> Test Denominator <==\n\n");

    boost::shared_ptr<Denominator> cholesky(Denominator::buildDenominator("CHOLESKY", evals_aocc_, evals_avir_,
        options_.get_double("DENOMINATOR_DELTA")));
    boost::shared_ptr<Denominator> laplace(Denominator::buildDenominator("LAPLACE", evals_aocc_, evals_avir_,
        options_.get_double("DENOMINATOR_DELTA")));
    if (debug_) {
        cholesky->debug(); 
        laplace->debug();
    } 
}
void MP2::test_ps_omega()
{
    fprintf(outfile, "  ==> Test PS Omega <==\n\n");

    double omega = options_.get_double("PS_OMEGA");
    if (!options_.get_bool("PS_USE_OMEGA")) {
        omega = -1.0;
    }

    boost::shared_ptr<PSTensorII> ps(new PSTensorII(basisset_, C_, nocc_, nvir_, naocc_, navir_, 0L, options_));
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    SharedMatrix I;
    if (omega >= 0.0) {
        I = mints->mo_erf_eri(omega, C_,C_);   
    } else {
        I = mints->mo_eri(C_,C_);   
    }
    SharedMatrix Ips = ps->Ipsmo();
    SharedMatrix E(new Matrix("Error in PS MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));

    E->copy(Ips);
    E->subtract(I);

    if (debug_) {
        I->print();
        Ips->print();
        E->print();
    }

    check_integrals(E);
}
void MP2::test_dps_omega()
{
    fprintf(outfile, "  ==> Test Double PS Omega <==\n\n");

    double omega = options_.get_double("PS_OMEGA");
    if (!options_.get_bool("PS_USE_OMEGA")) {
        omega = -1.0;
    }

    boost::shared_ptr<PSTensorII> ps(new PSTensorII(basisset_, C_, nocc_, nvir_, naocc_, navir_, 0L, options_));
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    SharedMatrix I;
    if (omega >= 0.0) {
        I = mints->mo_erf_eri(omega, C_,C_);   
    } else {
        I = mints->mo_eri(C_,C_);   
    }
    SharedMatrix Ips = ps->Idpsmo();
    SharedMatrix E(new Matrix("Error in DPS MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));

    E->copy(Ips);
    E->subtract(I);

    if (debug_) {
        I->print();
        Ips->print();
        E->print();
    }

    check_integrals(E);
}
void MP2::test_ps()
{
    fprintf(outfile, "  ==> Test PS <==\n\n");

    boost::shared_ptr<PSTensor> ps(new PSTensor(basisset_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix I = ps->Imo();
    SharedMatrix Ips = ps->Ipsmo();
    SharedMatrix E(new Matrix("Error in PS MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));

    E->copy(Ips);
    E->subtract(I);

    if (debug_) {
        I->print();
        Ips->print();
        E->print();
    }

    check_integrals(E);
}
void MP2::test_df()
{
    fprintf(outfile, "  ==> Test DF <==\n\n");

    boost::shared_ptr<DFTensor> df(new DFTensor(basisset_, ribasis_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix I = df->Imo();
    SharedMatrix Idf = df->Idfmo();
    SharedMatrix E(new Matrix("Error in DF MO ERI Tensor", nmo_ * nmo_, nmo_ * nmo_));

    E->copy(Idf);
    E->subtract(I);

    if (debug_) {
        I->print();
        Idf->print();
        E->print();
    }

    check_integrals(E);
}
void MP2::compute_MP2()
{
    fprintf(outfile, "  ==> Conventional MP2 <==\n\n");
    
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    SharedMatrix I = mints->mo_eri(C_aocc_, C_avir_);

    double** Ip = I->pointer();

    double E_MP2J = 0.0;
    double E_MP2K = 0.0;

    for (int i = 0; i < naocc_; i++) {
        for (int a = 0; a < nvir_; a++) {
            for (int j = 0; j < naocc_; j++) {
                for (int b = 0; b < navir_; b++) {
                    double denom = 1.0 / (evals_avirp_[a] + evals_avirp_[b] -
                        evals_aoccp_[i] - evals_aoccp_[j]);
                    double iajb = Ip[i * navir_ + a][j * navir_ + b];
                    double ibja = Ip[i * navir_ + b][j * navir_ + a];
                    E_MP2J -= 2.0 * iajb * iajb * denom;    
                    E_MP2K += 1.0 * iajb * ibja * denom;    
                }
            }
        }
    }

    energies_["MP2J Energy"] = E_MP2J;
    energies_["MP2K Energy"] = E_MP2K;
}
void MP2::compute_DF_MP2()
{
    fprintf(outfile, "  ==> DF-MP2 <==\n\n");

    boost::shared_ptr<DFTensor> df(new DFTensor(basisset_, ribasis_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix Qia = df->Qov();
    double** Qiap = Qia->pointer();
    int nQ = ribasis_->nbf();

    SharedMatrix I(new Matrix("DF (ia|jb)", naocc_ * navir_ , naocc_ * navir_));
    double** Ip = I->pointer();

    C_DGEMM('T','N', naocc_ * navir_, naocc_ * navir_, nQ, 1.0, Qiap[0], naocc_ * navir_,
        Qiap[0], naocc_ * navir_, 0.0, Ip[0], naocc_ * navir_);

    double E_MP2J = 0.0;
    double E_MP2K = 0.0;

    for (int i = 0; i < naocc_; i++) {
        for (int a = 0; a < nvir_; a++) {
            for (int j = 0; j < naocc_; j++) {
                for (int b = 0; b < navir_; b++) {
                    double denom = 1.0 / (evals_avirp_[a] + evals_avirp_[b] -
                        evals_aoccp_[i] - evals_aoccp_[j]);
                    double iajb = Ip[i * navir_ + a][j * navir_ + b];
                    double ibja = Ip[i * navir_ + b][j * navir_ + a];
                    E_MP2J -= 2.0 * iajb * iajb * denom;    
                    E_MP2K += 1.0 * iajb * ibja * denom;    
                }
            }
        }
    }

    energies_["MP2J Energy"] = E_MP2J;
    energies_["MP2K Energy"] = E_MP2K;
}
void MP2::compute_PS_MP2()
{
    fprintf(outfile, "  ==> PS-MP2 <==\n\n");

    boost::shared_ptr<PSTensor> ps(new PSTensor(basisset_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix Pia = ps->Aov();
    SharedMatrix Qi = ps->Qaocc();
    SharedMatrix Ra = ps->Ravir();
    double** Piap = Pia->pointer();
    double** Rap = Ra->pointer();
    double** Qip = Qi->pointer();
    int nP = Pia->rowspi()[0];

    SharedMatrix QRia(new Matrix("QR", nP, naocc_ * navir_));
    double** QRiap = QRia->pointer();

    for (int i = 0; i < naocc_; i++) {
        for (int a = 0; a < navir_; a++) {
            for (int P = 0; P < nP; P++) {
                QRiap[P][i * navir_ + a] = Qip[i][P] * Rap[a][P];
            }
        }
    }

    SharedMatrix I(new Matrix("PS (ia|jb)", naocc_ * navir_ , naocc_ * navir_));
    double** Ip = I->pointer();

    C_DGEMM('T','N', naocc_ * navir_, naocc_ * navir_, nP, 1.0, Piap[0], naocc_ * navir_,
        QRiap[0], naocc_ * navir_, 0.0, Ip[0], naocc_ * navir_);

    double E_MP2J = 0.0;
    double E_MP2K = 0.0;

    for (int i = 0; i < naocc_; i++) {
        for (int a = 0; a < nvir_; a++) {
            for (int j = 0; j < naocc_; j++) {
                for (int b = 0; b < navir_; b++) {
                    double denom = 1.0 / (evals_avirp_[a] + evals_avirp_[b] -
                        evals_aoccp_[i] - evals_aoccp_[j]);
                    double iajb = Ip[i * navir_ + a][j * navir_ + b];
                    double ibja = Ip[i * navir_ + b][j * navir_ + a];
                    E_MP2J -= 2.0 * iajb * iajb * denom;    
                    E_MP2K += 1.0 * iajb * ibja * denom;    
                }
            }
        }
    }

    energies_["MP2J Energy"] = E_MP2J;
    energies_["MP2K Energy"] = E_MP2K;
}
void MP2::compute_DF_MP2J()
{
    fprintf(outfile, "  ==> DF-MP2J <==\n\n");

    boost::shared_ptr<Denominator> denom(Denominator::buildDenominator(
        options_.get_str("DENOMINATOR_ALGORITHM"), evals_aocc_, evals_avir_,
        options_.get_double("DENOMINATOR_DELTA")));
    SharedMatrix tau = denom->denominator();
    double** taup = tau->pointer();
    int nW = denom->nvector();    

    boost::shared_ptr<DFTensor> df(new DFTensor(basisset_, ribasis_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix Qia = df->Qov();
    double** Qiap = Qia->pointer();
    int nQ = ribasis_->nbf();

    SharedMatrix Qiaw(new Matrix("(Q|ia)^w", nQ, naocc_*navir_));  
    double** Qiawp = Qiaw->pointer();

    SharedMatrix Z(new Matrix("Z^QQw", nQ, nQ));  
    double** Zp = Z->pointer();

    double E_MP2J = 0.0;

    for (int w = 0; w < nW; w++) {
        C_DCOPY(nQ * (ULI) naocc_ * navir_, Qiap[0], 1, Qiawp[0], 1);
        for (int ia = 0; ia < naocc_ * navir_; ia++) {
            C_DSCAL(nQ, taup[w][ia], &Qiawp[0][ia], naocc_ * navir_);    
        }
        C_DGEMM('N','T', nQ, nQ, naocc_ * navir_, 1.0, Qiap[0], naocc_ * navir_, Qiawp[0], naocc_ * navir_, 0.0, Zp[0], nQ);
        E_MP2J -= 2.0 * C_DDOT(nQ * (ULI) nQ, Zp[0], 1, Zp[0], 1);
    }
    
    energies_["MP2J Energy"] = E_MP2J;
}
void MP2::compute_PS_MP2J()
{
    fprintf(outfile, "  ==> PS-MP2J <==\n\n");

    boost::shared_ptr<Denominator> denom(Denominator::buildDenominator(
        options_.get_str("DENOMINATOR_ALGORITHM"), evals_aocc_, evals_avir_,
        options_.get_double("DENOMINATOR_DELTA")));
    SharedMatrix tau = denom->denominator();
    double** taup = tau->pointer();
    int nW = denom->nvector();    

    boost::shared_ptr<PSTensor> ps(new PSTensor(basisset_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix Pia = ps->Aov();
    SharedMatrix Qi = ps->Qaocc();
    SharedMatrix Ra = ps->Ravir();
    double** Piap = Pia->pointer();
    double** Rap = Ra->pointer();
    double** Qip = Qi->pointer();
    int nP = Pia->rowspi()[0];

    SharedMatrix X(new Matrix("X", nP, naocc_ * navir_));
    SharedMatrix Z(new Matrix("Z", nP, nP));
    
    double** Xp = X->pointer();
    double** Zp = Z->pointer();

    double E_MP2J = 0.0;

    for (int w = 0; w < nW; w++) {

        for (int i = 0; i < naocc_; i++) {
            for (int a = 0; a < navir_; a++) {
                for (int P = 0; P < nP; P++) {
                    Xp[P][i * navir_ + a] = Qip[i][P] * Rap[a][P] * taup[w][i * navir_ + a]; 
                }    
            }    
        }

        C_DGEMM('N','T', nP, nP, naocc_ * navir_, 1.0, Xp[0], naocc_ * navir_, 
            Piap[0], naocc_ * navir_, 0.0,  Zp[0], nP);

        // I think it's transposed
        for (int P = 0; P < nP; P++) {
            for (int Q = 0; Q < nP; Q++) {
                E_MP2J -= 2.0 * Zp[P][Q] * Zp[Q][P];
            }
        }
    }     

    energies_["MP2J Energy"] = E_MP2J;
}
void MP2::compute_PS_DF_MP2K()
{
    fprintf(outfile, "  ==> PS-DF-MP2K <==\n\n");

    boost::shared_ptr<Denominator> denom(Denominator::buildDenominator(
        options_.get_str("DENOMINATOR_ALGORITHM"), evals_aocc_, evals_avir_,
        options_.get_double("DENOMINATOR_DELTA")));
    SharedMatrix tau = denom->denominator();
    double** taup = tau->pointer();
    int nW = denom->nvector();    

    if (debug_) {
        tau->print();
    }

    boost::shared_ptr<DFTensor> df(new DFTensor(basisset_, ribasis_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix Qia = df->Qov();
    double** Qiap = Qia->pointer();
    int nQ = ribasis_->nbf();

    boost::shared_ptr<PSTensor> ps(new PSTensor(basisset_, C_, nocc_, nvir_, naocc_, navir_, options_));
    SharedMatrix Pia = ps->Aov();
    SharedMatrix Qi = ps->Qaocc();
    SharedMatrix Ra = ps->Ravir();
    double** Piap = Pia->pointer();
    double** Rap = Ra->pointer();
    double** Qip = Qi->pointer();
    int nP = Pia->rowspi()[0];

    SharedMatrix Qiaw(new Matrix("(Q|ia)^w", nQ, naocc_*navir_));  
    double** Qiawp = Qiaw->pointer();

    SharedMatrix X(new Matrix("X", nP, navir_ * nQ));
    SharedMatrix Y(new Matrix("Y", nP, naocc_ * nQ));
    SharedMatrix Z(new Matrix("Z", nP, naocc_ * navir_));
    
    double** Xp = X->pointer();
    double** Yp = Y->pointer();
    double** Zp = Z->pointer();

    SharedMatrix Temp(new Matrix("Transpose Temp", nQ, naocc_));
    double** Tempp = Temp->pointer();
    SharedMatrix Temp2(new Matrix("Transpose Temp", nP, navir_));
    double** Temp2p = Temp2->pointer();

    double E_MP2K = 0.0;

    for (int w = 0; w < nW; w++) {
        C_DCOPY(nQ * (ULI) naocc_ * navir_, Qiap[0], 1, Qiawp[0], 1);
        for (int ia = 0; ia < naocc_ * navir_; ia++) {
            C_DSCAL(nQ, taup[w][ia], &Qiawp[0][ia], naocc_ * navir_);    
        }
      
        for (int Q = 0; Q < nQ; Q++) {
            C_DGEMM('T','N', nP, navir_, naocc_, 1.0, Qip[0], nP, Qiawp[Q], navir_, 0.0, Temp2p[0], navir_);
            C_DCOPY(nP * (ULI) navir_, Temp2p[0], 1, &Xp[0][Q], nQ);      
        }
 
        C_DGEMM('T','T', nP, nQ * naocc_, navir_, 1.0, Rap[0], nP, Qiawp[0], navir_, 0.0, Yp[0], nQ * naocc_);  

        for (int P = 0; P < nP; P++) {
            C_DCOPY(nQ * (ULI) naocc_, Yp[P], 1, Tempp[0], 1);
            for (int Q = 0; Q < nQ; Q++) {
                C_DCOPY(naocc_, Tempp[Q], 1, &Yp[P][Q], nQ);
            }
        }

        for (int P = 0; P < nP; P++) {
            C_DGEMM('N','T', naocc_, navir_, nQ, 1.0, Yp[P], nQ, Xp[P], nQ, 0.0, Zp[P], navir_);
        }

        E_MP2K += C_DDOT(nP * (ULI) naocc_ * navir_, Zp[0], 1, Piap[0], 1);
    }     

    energies_["MP2K Energy"] = E_MP2K;
}
void MP2::compute_PS_MP2K()
{
    throw FeatureNotImplemented("dfcc::MP2", "compute_PS_MP2K", __FILE__, __LINE__);    

    fprintf(outfile, "  ==> PS-PS-MP2K <==\n\n");
}
void MP2::print_header()
{
    fprintf(outfile, "\t ********************************************************\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *                        DF-MP2                        *\n");
    fprintf(outfile, "\t *    2nd-Order Density-Fitted Moller-Plesset Theory    *\n");
    fprintf(outfile, "\t *        with Laplace and Pseudospectral Grids         *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t *            Rob Parrish and Ed Hohenstein             *\n");
    fprintf(outfile, "\t *                                                      *\n");
    fprintf(outfile, "\t ********************************************************\n");
    fprintf(outfile, "\n");
}
void MP2::print_energies()
{
    fprintf(outfile, "\n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ====================> MP2 Energies <==================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Reference Energy",         energies_["Reference Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2J Energy",              energies_["MP2J Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "MP2K Energy",              energies_["MP2K Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Same-Spin Energy",         energies_["Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Opposite-Spin Energy",     energies_["Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Correlation Energy",       energies_["Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "Total Energy",             energies_["Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t ==================> SCS-MP2 Energies <================== \n");
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Same-Spin Scale",      sss_);
    fprintf(outfile, "\t %-25s = %24.16f [-]\n", "SCS Opposite-Spin Scale",  oss_);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Same-Spin Energy",     energies_["SCS Same-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Opposite-Spin Energy", energies_["SCS Opposite-Spin Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Correlation Energy",   energies_["SCS Correlation Energy"]);
    fprintf(outfile, "\t %-25s = %24.16f [H]\n", "SCS Total Energy",         energies_["SCS Total Energy"]);
    fprintf(outfile, "\t----------------------------------------------------------\n");
    fprintf(outfile, "\n");
    fflush(outfile);
}

}}
