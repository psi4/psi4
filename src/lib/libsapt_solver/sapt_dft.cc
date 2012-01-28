//#include "sapt_dft.h"
//#include <libutil/quad.h>
//
//namespace psi { namespace sapt {
//
//MP2C::MP2C(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt)
//    : SAPT0(options, psio, chkpt)
//{
//}
//
//MP2C::~MP2C()
//{
//}
//
//double MP2C::compute_energy()
//{
//    print_header();
//
//    form_quadrature();
//    form_J();
//    form_W();
//    form_X0();   
//    std::vector<double> Edisp20U = casimirPolder();
//    form_XC();   
//    std::vector<double> Edisp20C = casimirPolder();
//
//    fprintf(outfile, "\n");
//    fprintf(outfile, "  ------------------------------------------------------------------------------\n");       
//    fprintf(outfile, "   =======================> CASIMIR-POLDER INTEGRATION <=======================\n");       
//    fprintf(outfile, "  ------------------------------------------------------------------------------\n");       
//    fprintf(outfile, "   Point     Omega           Weight          E_UCHF [mH]         E_TDDFT [mH]\n");       
//    fprintf(outfile, "  ------------------------------------------------------------------------------\n");       
//    fflush(outfile);
//   
//    double E_UCHF = 0.0;
//    double E_TDDFT = 0.0;
//    double E_DeltaMP2C = 0.0;
//    double E_MP2_int = Process::environment.globals["MP2C DIMER MP2 ENERGY"] -
//                       Process::environment.globals["MP2C MONOMER A MP2 ENERGY"] -
//                       Process::environment.globals["MP2C MONOMER B MP2 ENERGY"];
//    double E_MP2C_int = 0.0;
//
//    int n = 0;
//    for (quad_->reset(); !quad_->isDone(); quad_->nextPoint() ) {
//    
//        double omega = quad_->getPoint();
//        double weight = quad_->getWeight();
//   
//        double UCHF = Edisp20U[n];
//        double TDDFT = Edisp20C[n];
//
//        E_UCHF += Edisp20U[n] * weight;
//        E_TDDFT += Edisp20C[n] * weight;
// 
//        n++;
//        fprintf(outfile, "  %3d  %12.8E  %12.8E  %18.12f  %18.12f\n", \
//           n, omega, weight, UCHF*1000.0,TDDFT*1000.0);       
//
//    }
//
//    E_DeltaMP2C = E_TDDFT - E_UCHF; 
//    E_MP2C_int = E_MP2_int + E_DeltaMP2C; 
//
//    fprintf(outfile, "  ------------------------------------------------------------------------------\n");       
//    fprintf(outfile, "   @ UCHF Dispersion Energy:  %18.12f [mH] %18.12f [kcal]\n", E_UCHF*1000.0, \
//        E_UCHF*627.509);   
//    fprintf(outfile, "   @ TDDFT Dispersion Energy: %18.12f [mH] %18.12f [kcal]\n", E_TDDFT*1000.0, \
//        E_TDDFT*627.509);   
//    fprintf(outfile, "   @ Delta MP2C Energy:       %18.12f [mH] %18.12f [kcal]\n", E_DeltaMP2C*1000.0, \
//        E_DeltaMP2C*627.509);   
//    fprintf(outfile, "   @ MP20 Interaction Energy: %18.12f [mH] %18.12f [kcal]\n", E_MP2_int*1000.0, \
//        E_MP2_int*627.509);   
//    fprintf(outfile, "   @ MP2C Interaction Energy: %18.12f [mH] %18.12f [kcal]\n", E_MP2C_int*1000.0, \
//        E_MP2C_int*627.509);   
//    fprintf(outfile, "  ------------------------------------------------------------------------------\n\n");       
//
//    fprintf(outfile, "  --\"We're looking for a bug no one's ever seen before...some kinda smart bug.\"\n");
//    
//    fflush(outfile); 
//
//    Process::environment.globals["MP2C UCHF ENERGY"] = E_UCHF; 
//    Process::environment.globals["MP2C TDDFT ENERGY"] = E_TDDFT; 
//    Process::environment.globals["MP2C DELTA MP2C ENERGY"] = E_DeltaMP2C; 
//    Process::environment.globals["MP2C MP2 ENERGY"] = E_MP2_int; 
//    Process::environment.globals["MP2C MP2C ENERGY"] = E_MP2C_int; 
//    Process::environment.globals["CURRENT ENERGY"] = Process::environment.globals["MP2C MP2C ENERGY"];
// 
//    return E_MP2C_int; 
//}
//
//void MP2C::form_quadrature()
//{
//    quad_ = boost::shared_ptr<Quadrature>(new ChebyshevIIQuadrature(options_.get_int("OMEGA_POINTS"),
//        options_.get_double("OMEGA_CENTER"))); 
//
//    if (debug_) {
//        quad_->print();
//    }
//}
//
//void MP2C::form_J()
//{
//    int naux = ribasis_->nbf();
//
//    // ==> Allocation <== // 
//    Jinv_ = SharedMatrix(new Matrix("J^-1",naux,naux));
//    J_ = SharedMatrix(new Matrix("J",naux,naux));
//    double** Jp = Jinv_->pointer();
//
//    // ==> J <== //
//    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(ribasis_,zero_,ribasis_,zero_));
//    boost::shared_ptr<TwoBodyAOInt> Qint(fact->eri());
//    const double* buffer = Qint->buffer();
//    
//    for (int A = 0; A < ribasis_->nshell(); A++) {
//        for (int B = 0; B <= A; B++) {
//            int nA = ribasis_->shell(A)->nfunction();
//            int nB = ribasis_->shell(B)->nfunction();
//            int sA = ribasis_->shell(A)->function_index();
//            int sB = ribasis_->shell(B)->function_index();
//            Qint->compute_shell(A,0,B,0);
//            for (int a = 0, index = 0; a < nA; a++) {
//                int oa = sA + a;
//                for (int b = 0; b < nB; b++, index++) {
//                    int ob = sB + b;
//                    Jp[oa][ob] = buffer[index];
//                    Jp[ob][oa] = buffer[index];
//                }
//            }
//        }
//    }
//    Qint.reset();   
//    fact.reset();
// 
//    if (debug_) {
//        Jinv_->print(outfile, "Before Inversion");
//    }
//
//    // ==> J^-1 <== //
//    J_->copy(Jinv_);
//    Jinv_->power(-1.0,1.0E-12);
//    
//    if (debug_) {
//        Jinv_->print(outfile, "After Inversion");
//    }
//}
//
//void MP2C::form_W()
//{
//    int naux = ribasis_->nbf();
//
//    // ==> Allocation <== //
//    W_A_ = SharedMatrix (new Matrix("W_A",naux,naux));  
//    W_B_ = SharedMatrix (new Matrix("W_B",naux,naux)); 
//    double** W_Ap = W_A_->pointer();
//    double** W_Bp = W_B_->pointer();
//
//    // ==> D_A/D_B <== // 
//    SharedMatrix D_A(new Matrix("D_A",nso_,nso_));
//    SharedMatrix D_B(new Matrix("D_B",nso_,nso_));
//
//    double** D_Ap = D_A->pointer();
//    double** D_Bp = D_B->pointer();
//
//    C_DGEMM('N','T',nso_,nso_,noccA_, 1.0, CA_[0], nmo_, CA_[0], nmo_, 0.0, D_Ap[0], nso_);
//    C_DGEMM('N','T',nso_,nso_,noccB_, 1.0, CB_[0], nmo_, CB_[0], nmo_, 0.0, D_Bp[0], nso_);
//
//    if (debug_) {
//        D_A->print();
//        D_B->print();
//    }
//    
//    // ==> S <== //
//    SharedMatrix S_A = SharedMatrix (new Matrix("S",naux,naux));
//    SharedMatrix S = SharedMatrix (new Matrix("S",naux,naux));
//
//    boost::shared_ptr<IntegralFactory> Afact(new IntegralFactory(ribasis_,ribasis_,zero_,zero_));
//    boost::shared_ptr<OneBodyAOInt> Aint(Afact->ao_overlap());
//
//    Aint->compute(S_A);
//    S->copy(S_A);   
// 
//    Afact.reset();
//    Aint.reset();    
//
//    if (debug_) {
//        S_A->print();
//    }  
//
//    // ==> V <== //
//    SharedMatrix V (new Matrix("V",naux,naux));
//    boost::shared_ptr<Vector> Veig (new Vector("Veig",naux));
//    double** Vp = V->pointer();
//    double* Veigp = Veig->pointer();
//    double maxS = Veigp[naux - 1];
//
//    S_A->diagonalize(V,Veig);
//
//    if (debug_) {
//        V->eivprint(Veig);
//    }
//    
//    for (int k = 0; k < naux; k++) {
//        if (Veigp[k] < 1.0E-9 * maxS) {
//            Veigp[k] = 0.0;
//        } else {
//            Veigp[k] = pow(Veigp[k],-1.0 / 2.0);
//        } 
//        C_DSCAL(naux,Veigp[k],&Vp[0][k],naux);
//    }
//
//    if (debug_) {
//        V->print();
//
//        S_A->gemm(true,false,1.0,V,S,0.0);
//        SharedMatrix one (new Matrix("Assert I",naux,naux));
//        one->gemm(false,false,1.0,S_A,V,0.0);
//        one->print();
//    } 
//
//    Veig.reset();
//    S_A.reset(); 
//
//    // ==> c <== //
//    boost::shared_ptr<Vector> c_A(new Vector("c_A",naux));
//    boost::shared_ptr<Vector> c_B(new Vector("c_B",naux));
//    double* c_Ap = c_A->pointer();
//    double* c_Bp = c_B->pointer();
//
//    int maxA = ribasis_->max_function_per_shell();
//    SharedMatrix QmnA (new Matrix("(Q|mn) A", maxA, nso_ * (ULI) nso_));
//    double** QmnAp = QmnA->pointer();
//    
//    boost::shared_ptr<IntegralFactory> Afact2(new IntegralFactory(ribasis_, zero_, basisset_, basisset_));
//    boost::shared_ptr<TwoBodyAOInt> Aeri (Afact2->eri());
//    const double* Abuffer = Aeri->buffer();
//
//    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
//        int nQ = ribasis_->shell(Q)->nfunction();
//        int sQ = ribasis_->shell(Q)->function_index();
//
//        // ==> Integrals <== //
//        for (int M = 0; M < basisset_->nshell(); M++) {
//            int nM = basisset_->shell(M)->nfunction();
//            int sM = basisset_->shell(M)->function_index();
//            for (int N = 0; N <= M; N++) {
//                int nN = basisset_->shell(N)->nfunction();
//                int sN = basisset_->shell(N)->function_index();
//                Aeri->compute_shell(Q,0,M,N);
//                for (int dQ = 0, index = 0; dQ < nQ; dQ++) {
//                    int oQ = sQ + dQ;
//                    for (int dM = 0; dM < nM; dM++) {
//                        int oM = sM + dM;
//                        for (int dN = 0; dN < nN; dN++, index++) {
//                            int oN = sN + dN;
//                            QmnAp[dQ][oM*nso_ + oN] = Abuffer[index];
//                            QmnAp[dQ][oN*nso_ + oM] = Abuffer[index];
//                        }
//                    }
//                }
//            } 
//        } 
//
//        //if (debug_) {
//        //    QmnA->print();
//        //}
//   
//        // ==> (A|mn)D_mn <== //
//        C_DGEMV('N',nQ,nso_*(ULI)nso_,1.0,QmnAp[0],nso_*(ULI)nso_,D_Ap[0],1,0.0,&c_Ap[sQ],1);
//        C_DGEMV('N',nQ,nso_*(ULI)nso_,1.0,QmnAp[0],nso_*(ULI)nso_,D_Bp[0],1,0.0,&c_Bp[sQ],1);
//    }
//
//    Afact2.reset(); 
//    Aeri.reset();
//    QmnA.reset();
//
//    if (debug_) {
//        c_A->print();
//        c_B->print();
//    }
//
//    // ==> d <== //
//    boost::shared_ptr<Vector> d_A(new Vector("d_A",naux));
//    boost::shared_ptr<Vector> d_B(new Vector("d_B",naux));
//    double* d_Ap = d_A->pointer();
//    double* d_Bp = d_B->pointer();
//   
//    double** Jinvp = Jinv_->pointer();  
//    
//    C_DGEMV('N',naux,naux,1.0,Jinvp[0],naux,c_Ap,1,0.0,d_Ap,1);
//    C_DGEMV('N',naux,naux,1.0,Jinvp[0],naux,c_Bp,1,0.0,d_Bp,1);
//
//    if (debug_) {
//        d_A->print();
//        d_B->print();
//    }
//
//    c_A.reset();
//    c_B.reset();
//    
//    // ==> (M|\rho|Q) <== //
//    boost::shared_ptr<IntegralFactory> PQRfactory(new IntegralFactory(ribasis_, ribasis_, ribasis_, ribasis_));
//
//    SharedMatrix PQR (new Matrix("PQR",maxA*maxA,naux));
//    boost::shared_ptr<Vector> TempA (new Vector("TempA",maxA*maxA)); 
//    boost::shared_ptr<Vector> TempB (new Vector("TempB",maxA*maxA)); 
//    double** PQRp = PQR->pointer(); 
//    double* TempAp = TempA->pointer();
//    double* TempBp = TempB->pointer();
// 
//    boost::shared_ptr<ThreeCenterOverlapInt> o3(PQRfactory->overlap_3c());
//    const double* buffer = o3->buffer();
// 
//    // A bit naive at the moment (no sieves or threading)
//    for (int P=0; P < ribasis_->nshell(); ++P) {
//        int numP = ribasis_->shell(P)->nfunction();
//        for (int Q=0; Q< ribasis_->nshell(); ++Q) {
//            int numQ = ribasis_->shell(Q)->nfunction();
//            
//            // ==> (PQR) Integrals <== //
//            for (int R=0; R < ribasis_->nshell();  ++R) {
//                int numR = ribasis_->shell(R)->nfunction();
//
//                o3->compute_shell(P, Q, R);
//                
//                for (int p=0 ; p < numP; ++p) {
//                    int op = ribasis_->shell(P)->function_index() + p;
//                    for (int q=0; q < numQ; ++q) {
//                        int oq = ribasis_->shell(Q)->function_index() + q;
//                        for (int r=0; r < numR; ++r) {
//                            int oR = ribasis_->shell(R)->function_index() + r;
//                            PQRp[p*numQ + q][r + ribasis_->shell(R)->function_index()] = \
//                                buffer[p*numQ*numR+q*numR+r];
//                        }
//                    }
//                }
//            }
//
//            // ==> (PQR) d_R <== //
//            C_DGEMV('n', numP*numQ, naux, 1.0, PQRp[0], naux, d_Ap, 1, 0.0, TempAp, 1);
//            C_DGEMV('n', numP*numQ, naux, 1.0, PQRp[0], naux, d_Bp, 1, 0.0, TempBp, 1);
//
//            for (int p=0 ; p < numP; ++p) {
//                int op = ribasis_->shell(P)->function_index() + p;
//                for (int q=0; q < numQ; ++q) {
//                    int oq = ribasis_->shell(Q)->function_index() + q;
//                    W_Ap[op][oq] = TempAp[p*numQ + q];
//                    W_Bp[op][oq] = TempBp[p*numQ + q];
//                }
//            }
//        }
//    }
//
//    PQRfactory.reset();
//    o3.reset();
//    PQR.reset();
//    TempA.reset();
//    TempB.reset();
//    d_A.reset();
//    d_B.reset();
//    
//    if (debug_) {
//        W_A_->print();
//        W_B_->print();
//    }
//
//    // ==> M-tilde <== //
//    SharedMatrix T(new Matrix("T",naux,naux));
//    SharedMatrix T2(new Matrix("T2",naux,naux));
//    double** Tp = T->pointer();
//    double** T2p = T2->pointer();
//
//    T->gemm(true,false,1.0,V,W_A_,0.0);
//    W_A_->gemm(false,false,1.0,T,V,0.0);
//    T->gemm(true,false,1.0,V,W_B_,0.0);
//    W_B_->gemm(false,false,1.0,T,V,0.0);
//    
//    if (debug_) {
//        W_A_->print();
//        W_B_->print();
//    }
//
//    // ==> ALDAx Kernel <== //
//    double C_x = 3.0/8.0*pow(3.0,1.0/3.0)*pow(4.0,2.0/3.0)*pow(M_PI,-1.0/3.0); 
//    boost::shared_ptr<Vector> lambda(new Vector("Lambda",naux));
//    double* lambdap = lambda->pointer();
//
//    // Diagonalize
//    W_A_->diagonalize(T,lambda);
//
//    if (debug_) {
//        T->eivprint(lambda);
//    }
//
//    // Apply kernel
//    bool warn = false;
//    for (int k = 0; k < naux; k++) {
//        if (lambdap[k] < 1.0E-10) {
//            warn = true;
//            lambdap[k] = 0.0;
//        } else {
//            lambdap[k] = -8.0/9.0 * C_x * pow(lambdap[k], -2.0/3.0);
//        }
//    }
//    if (warn)       
//        fprintf(outfile, "  WARNING: Small/negative eigenvalue detected in (P|\\rho|Q)\n");
//   
//    T2->gemm(false,false,1.0,V,T,0.0);
//    T->gemm(false,false,1.0,S,T2,0.0);
//    T2->copy(T);
//    
//    for (int k = 0; k < naux; k++)
//        C_DSCAL(naux,lambdap[k],&Tp[0][k],naux);     
//
//    W_A_->gemm(false,true,1.0,T,T2,0.0);
//
//    // Diagonalize
//    W_B_->diagonalize(T,lambda);
//
//    if (debug_) {
//        T->eivprint(lambda);
//    }
//
//    // Apply kernel
//    warn = false;
//    for (int k = 0; k < naux; k++) {
//        if (lambdap[k] < 1.0E-10) {
//            warn = true;
//            lambdap[k] = 0.0;
//        } else {
//            lambdap[k] = -8.0/9.0 * C_x * pow(lambdap[k], -2.0/3.0);
//        }
//    }
//    if (warn)       
//        fprintf(outfile, "  WARNING: Small/negative eigenvalue detected in (P|\\rho|Q)\n");
//   
//    T2->gemm(false,false,1.0,V,T,0.0);
//    T->gemm(false,false,1.0,S,T2,0.0);
//    T2->copy(T);
//    
//    for (int k = 0; k < naux; k++)
//        C_DSCAL(naux,lambdap[k],&Tp[0][k],naux);     
//
//    W_B_->gemm(false,true,1.0,T,T2,0.0);
//
//    if (debug_) {
//        W_A_->print();
//        W_B_->print();
//    }
//
//    W_A_->add(J_); 
//    W_B_->add(J_); 
//
//    J_.reset();
//
//    if (debug_) {
//        W_A_->print();
//        W_B_->print();
//    }
//}
//
//void MP2C::form_X0()
//{
//    int NA = ribasis_->nbf(); 
//    int NB = ribasis_->nbf(); 
//
//    std::vector<double> omega;
//    std::vector<double> weight;
//
//    for (quad_->reset(); !quad_->isDone(); quad_->nextPoint()) {
//        omega.push_back(quad_->getPoint());
//        weight.push_back(quad_->getWeight());
//    }
//
//    // ==> X0_A <== //
//    int maxA = ribasis_->max_function_per_shell();
//    SharedMatrix QiaA (new Matrix("(Q|ia) A", NA, noccA_ * (ULI) nvirA_)); 
//    SharedMatrix QmnA (new Matrix("(Q|mn) A", maxA, nso_ * (ULI) nso_));
//    SharedMatrix QmiA (new Matrix("(Q|mi) A", maxA, noccA_ * (ULI) nso_));
//    double** QiaAp = QiaA->pointer();
//    double** QmnAp = QmnA->pointer();
//    double** QmiAp = QmiA->pointer();
//    
//    boost::shared_ptr<IntegralFactory> Afact(new IntegralFactory(ribasis_, zero_, basisset_, basisset_));
//    boost::shared_ptr<TwoBodyAOInt> Aeri (Afact->eri());
//    const double* Abuffer = Aeri->buffer();
//
//    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
//        int nQ = ribasis_->shell(Q)->nfunction();
//        int sQ = ribasis_->shell(Q)->function_index();
//
//        // ==> Integrals <== //
//        for (int M = 0; M < basisset_->nshell(); M++) {
//            int nM = basisset_->shell(M)->nfunction();
//            int sM = basisset_->shell(M)->function_index();
//            for (int N = 0; N <= M; N++) {
//                int nN = basisset_->shell(N)->nfunction();
//                int sN = basisset_->shell(N)->function_index();
//                Aeri->compute_shell(Q,0,M,N);
//                for (int dQ = 0, index = 0; dQ < nQ; dQ++) {
//                    int oQ = sQ + dQ;
//                    for (int dM = 0; dM < nM; dM++) {
//                        int oM = sM + dM;
//                        for (int dN = 0; dN < nN; dN++, index++) {
//                            int oN = sN + dN;
//                            QmnAp[dQ][oM*nso_ + oN] = Abuffer[index];
//                            QmnAp[dQ][oN*nso_ + oM] = Abuffer[index];
//                        }
//                    }
//                }
//            } 
//        } 
//
//        //if (debug_) {
//        //    QmnA->print();
//        //}
//    
//        // First Half-transform
//        C_DGEMM('N','N',nQ*(ULI)nso_,noccA_,nso_,1.0,QmnAp[0],nso_,CA_[0],nmo_,0.0,QmiAp[0],noccA_); 
//
//        // Second Half-transform
//        for (int dQ = 0; dQ < nQ; dQ++) {
//            int oQ = sQ + dQ;
//            C_DGEMM('T','N',noccA_,nvirA_,nso_,1.0,QmiAp[dQ],noccA_,&CA_[0][noccA_],nmo_,0.0,QiaAp[oQ],nvirA_); 
//        }        
//    }
//
//    Afact.reset(); 
//    Aeri.reset();
//
//    if (debug_) {
//        QiaA->print();
//    }
//
//    // ==> Lambda <== //
//    SharedMatrix LiaA (new Matrix("Lia A", omega.size(), noccA_*(ULI)nvirA_));
//    double** LiaAp = LiaA->pointer();
//    for (int om = 0; om < omega.size(); om++) {
//        double OM = omega[om];
//        for (int i = 0; i < noccA_; i++) {
//            for (int a = 0; a < nvirA_; a++) {
//                double eps_ia = evalsA_[a + noccA_] - evalsA_[i];
//                LiaAp[om][i * nvirA_ + a] = sqrt(4.0 * eps_ia / (eps_ia * eps_ia + OM * OM)); 
//            }
//        }
//    }
//
//    if (debug_) {
//        LiaA->print();
//    }
//
//    // ==> X0 <== //
//    for (int om = 0; om < omega.size(); om++) {
//        SharedMatrix XA (new Matrix("XA", NA, NA));
//
//        if (om == 0) {
//            for (ULI ia = 0; ia < noccA_*(ULI)nvirA_; ia++) {
//                C_DSCAL(NA,LiaAp[om][ia],&QiaAp[0][ia],noccA_*(ULI)nvirA_);
//            }
//        } else {
//            for (ULI ia = 0; ia < noccA_*(ULI)nvirA_; ia++) {
//                C_DSCAL(NA,LiaAp[om][ia]/LiaAp[om - 1][ia],&QiaAp[0][ia],noccA_*(ULI)nvirA_);
//            }
//        }
//
//        XA->gemm(false,true,1.0,QiaA,QiaA,0.0);
//
//        X_A_.push_back(XA);
//    }
//
//    QiaA.reset();
//    LiaA.reset();
//
//    // ==> X0_B <== //
//    int maxB = ribasis_->max_function_per_shell();
//    SharedMatrix QiaB (new Matrix("(Q|ia) B", NB, noccB_ * (ULI) nvirB_)); 
//    SharedMatrix QmnB (new Matrix("(Q|mn) B", maxB, nso_ * (ULI) nso_));
//    SharedMatrix QmiB (new Matrix("(Q|mi) B", maxB, noccB_ * (ULI) nso_));
//    double** QiaBp = QiaB->pointer();
//    double** QmnBp = QmnB->pointer();
//    double** QmiBp = QmiB->pointer();
//    
//    boost::shared_ptr<IntegralFactory> Bfact(new IntegralFactory(ribasis_, zero_, basisset_, basisset_));
//    boost::shared_ptr<TwoBodyAOInt> Beri (Bfact->eri());
//    const double* Bbuffer = Beri->buffer();
//
//    for (int Q = 0; Q < ribasis_->nshell(); Q++) {
//        int nQ = ribasis_->shell(Q)->nfunction();
//        int sQ = ribasis_->shell(Q)->function_index();
//
//        // ==> Integrals <== //
//        for (int M = 0; M < basisset_->nshell(); M++) {
//            int nM = basisset_->shell(M)->nfunction();
//            int sM = basisset_->shell(M)->function_index();
//            for (int N = 0; N <= M; N++) {
//                int nN = basisset_->shell(N)->nfunction();
//                int sN = basisset_->shell(N)->function_index();
//                Beri->compute_shell(Q,0,M,N);
//                for (int dQ = 0, index = 0; dQ < nQ; dQ++) {
//                    int oQ = sQ + dQ;
//                    for (int dM = 0; dM < nM; dM++) {
//                        int oM = sM + dM;
//                        for (int dN = 0; dN < nN; dN++, index++) {
//                            int oN = sN + dN;
//                            QmnBp[dQ][oM*nso_ + oN] = Bbuffer[index];
//                            QmnBp[dQ][oN*nso_ + oM] = Bbuffer[index];
//                        }
//                    }
//                }
//            } 
//        } 
//   
//        //if (debug_) {
//        //    QmnB->print();
//        //}
//
//        // First Half-transform
//        C_DGEMM('N','N',nQ*(ULI)nso_,noccB_,nso_,1.0,QmnBp[0],nso_,CB_[0],nmo_,0.0,QmiBp[0],noccB_); 
//
//        // Second Half-transform
//        for (int dQ = 0; dQ < nQ; dQ++) {
//            int oQ = sQ + dQ;
//            C_DGEMM('T','N',noccB_,nvirB_,nso_,1.0,QmiBp[dQ],noccB_,&CB_[0][noccB_],nmo_,0.0,QiaBp[oQ],nvirB_); 
//        }        
//    }
//
//    Bfact.reset(); 
//    Beri.reset();
//
//    if (debug_) {
//        QiaB->print();
//    }
//
//    // ==> Lambda <== //
//    SharedMatrix LiaB (new Matrix("Lia B", omega.size(), noccB_*(ULI)nvirB_));
//    double** LiaBp = LiaB->pointer();
//    for (int om = 0; om < omega.size(); om++) {
//        double OM = omega[om];
//        for (int i = 0; i < noccB_; i++) {
//            for (int a = 0; a < nvirB_; a++) {
//                double eps_ia = evalsB_[a + noccB_] - evalsB_[i];
//                LiaBp[om][i * nvirB_ + a] = sqrt(4.0 * eps_ia / (eps_ia * eps_ia + OM * OM)); 
//            }
//        }
//    }
//
//    if (debug_) {
//        LiaB->print();
//    }
//
//    // ==> X0 <== //
//    for (int om = 0; om < omega.size(); om++) {
//        SharedMatrix XB (new Matrix("XB", NB, NB));
//
//        if (om == 0) {
//            for (ULI ia = 0; ia < noccB_*(ULI)nvirB_; ia++) {
//                C_DSCAL(NB,LiaBp[om][ia],&QiaBp[0][ia],noccB_*(ULI)nvirB_);
//            }
//        } else {
//            for (ULI ia = 0; ia < noccB_*(ULI)nvirB_; ia++) {
//                C_DSCAL(NB,LiaBp[om][ia]/LiaBp[om - 1][ia],&QiaBp[0][ia],noccB_*(ULI)nvirB_);
//            }
//        }
//
//        XB->gemm(false,true,1.0,QiaB,QiaB,0.0);
//
//        X_B_.push_back(XB);
//    }
//
//    QiaB.reset();
//    LiaB.reset();
//}
//
//void MP2C::form_XC()
//{
//    int naux = ribasis_->nbf();
//    int* piv = new int[naux];
//    int lwork = 3*naux;
//    double* work = new double[lwork]; 
//
//    SharedMatrix T1(new Matrix("T1",naux,naux));
//    SharedMatrix T2(new Matrix("T2",naux,naux));
//
//    // ==> S^-1 W S^-1 => W <== //
//    T1->gemm(false,false,1.0,Jinv_,W_A_,0.0);
//    W_A_->gemm(false,false,1.0,T1,Jinv_,0.0);
//    T1->gemm(false,false,1.0,Jinv_,W_B_,0.0);
//    W_B_->gemm(false,false,1.0,T1,Jinv_,0.0);
//
//    if (debug_) {
//        Jinv_->print();
//        W_A_->print();
//        W_B_->print();
//    }
//
//    // ==> XC_A <== //
//    for (int om = 0; om < X_A_.size(); om++) {
//
//        if (debug_) 
//            X_A_[om]->print();
//
//        T1->gemm(false,false,-1.0,X_A_[om],W_A_,0.0);
//        if (debug_) 
//            T1->print();
//
//        double** T1p = T1->pointer();
//        for (int i = 0; i < naux; i++)
//            T1p[i][i] += 1.0;
//
//        if (debug_) 
//            T1->print();
//        
//        int error = C_DGETRF(naux,naux,T1p[0],naux,piv);
//        error |= C_DGETRI(naux,T1p[0],naux,piv,work,lwork);
//
//        if (error != 0)
//            throw PSIEXCEPTION("MP2C::form_XC: LU inverse failed in Dyson Equation");
//
//        if (debug_) 
//            T1->print();
//        
//        T2->copy(X_A_[om]);
//
//        X_A_[om]->gemm(false,false,1.0,T1,T2,0.0);
//        if (debug_) 
//            X_A_[om]->print();
//    }
//
//    // ==> XC_B <== //
//    for (int om = 0; om < X_A_.size(); om++) {
//
//        if (debug_) 
//            X_B_[om]->print();
//
//        T1->gemm(false,false,-1.0,X_B_[om],W_B_,0.0);
//        if (debug_) 
//            T1->print();
//
//        double** T1p = T1->pointer();
//        for (int i = 0; i < naux; i++)
//            T1p[i][i] += 1.0;
//
//        if (debug_) 
//            T1->print();
//        
//        int error = C_DGETRF(naux,naux,T1p[0],naux,piv);
//        error |= C_DGETRI(naux,T1p[0],naux,piv,work,lwork);
//
//        if (error != 0)
//            throw PSIEXCEPTION("MP2C::form_XC: LU inverse failed in Dyson Equation");
//
//        if (debug_) 
//            T1->print();
//        
//        T2->copy(X_B_[om]);
//
//        X_B_[om]->gemm(false,false,1.0,T1,T2,0.0);
//        if (debug_) 
//            X_B_[om]->print();
//    }
//
//    delete[] piv;
//    delete[] work;
//}
//
//std::vector<double> MP2C::casimirPolder()
//{
//    int naux = ribasis_->nbf(); 
//
//    SharedMatrix T1(new Matrix("T1",naux,naux));
//    SharedMatrix T2(new Matrix("T2",naux,naux));
//
//    std::vector<double> E;
//
//    for (int om = 0; om < X_A_.size(); om++) {
//
//        T1->gemm(false, false, 1.0, X_A_[om], Jinv_, 0.0);
//        T2->gemm(false, false, 1.0, Jinv_, X_B_[om], 0.0);
//
//        // Vector dot, and -1/(2\pi). 
//        // The quadrature weight goes in later
//        E.push_back(-1.0 / (2.0 * M_PI) * T1->vector_dot(T2));    
//    }
//
//    return E;
//}
//
//void MP2C::print_header()
//{
//  int nthread = 1;
//  #ifdef _OPENMP
//      nthread = omp_get_max_threads();
//  #endif
//
//  fprintf(outfile, "\n");
//  fprintf(outfile, "         ---------------------------------------------------------\n");
//  fprintf(outfile, "                                   MP2C\n");
//  fprintf(outfile, "                     by Rob Parrish, and Ed Hohenstein\n");
//  fprintf(outfile, "                      %3d Threads, %6ld MiB Core\n", nthread, memory_ / 1000000L);
//  fprintf(outfile, "         ---------------------------------------------------------\n");
//  fprintf(outfile, "\n");
//
//  fprintf(outfile,"\n");
//  fprintf(outfile,"    Orbital Information\n");
//  fprintf(outfile,"  -----------------------\n");
//  fprintf(outfile,"    NSO     = %9d\n",nso_);
//  fprintf(outfile,"    NMO     = %9d\n",nmo_);
//  fprintf(outfile,"    NRI     = %9d\n",ndf_);
//  fprintf(outfile,"    NOCC A  = %9d\n",noccA_);
//  fprintf(outfile,"    NOCC B  = %9d\n",noccB_);
//  fprintf(outfile,"    FOCC A  = %9d\n",foccA_);
//  fprintf(outfile,"    FOCC B  = %9d\n",foccB_);
//  fprintf(outfile,"    NVIR A  = %9d\n",nvirA_);
//  fprintf(outfile,"    NVIR B  = %9d\n",nvirB_);
//  fprintf(outfile,"\n");
//  fflush(outfile);
//}
//
//void MP2C::print_results()
//{
//}
//
//
//}}
