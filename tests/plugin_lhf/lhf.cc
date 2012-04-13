#include "lhf.h" 
#include <psifiles.h>
#include <libmints/view.h>
#include <libmints/mints.h>
#include <libfock/v.h>
#include <lib3index/3index.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <libfock/jk.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sstream>

using namespace psi;

namespace psi{ namespace scf{

LHF::LHF(Options& options) 
    : RKS(options, _default_psio_lib_)
{
    common_init();
}
LHF::~LHF()
{
}
void LHF::common_init()
{
    if (Ca_->nirrep() != 1){
        throw PSIEXCEPTION("LHF: Must Run in C1");    
    }

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    auxiliary_ = BasisSet::construct(parser, Process::environment.molecule(), "DF_BASIS_SCF");

    EXX_ = 0.0;
}
SharedMatrix LHF::build_Dij()
{
    int nso = KS::basisset_->nbf();
    int nmo = nmopi_[0];
    int ni  = nalphapi_[0];
    int naux = auxiliary_->nbf();    

    SharedMatrix Dij(new Matrix("D_A^ij", naux, ni * ni));
    SharedMatrix Dmj(new Matrix("D_A^mj", naux, nso * ni));
    SharedMatrix Dmn(new Matrix("D_A^mn", naux, nso * nso));
    
    boost::shared_ptr<FittingMetric> metric(new FittingMetric(auxiliary_)); 
    metric->form_full_eig_inverse();
    SharedMatrix J = metric->get_metric();
    metric.reset();

    SharedMatrix T(new Matrix("T", naux, naux));
   
    double** Dijp = Dij->pointer();
    double** Dmjp = Dmj->pointer();  
    double** Dmnp = Dmn->pointer();  
    double** Jp = J->pointer();
    double** Tp = T->pointer();
    double** Cp = Ca_->pointer();

    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), KS::basisset_, KS::basisset_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > ints;
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif 
    for (int i = 0; i < nthread; i++) {
        ints.push_back(boost::shared_ptr<TwoBodyAOInt>(factory->eri()));
    }

    boost::shared_ptr<BasisSet> primary = KS::basisset_;

    // Integrals
    #pragma omp parallel for num_threads(nthread)
    for (int P = 0; P < auxiliary_->nshell(); P++) {
        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
        boost::shared_ptr<TwoBodyAOInt> eri = ints[thread];
        const double* buffer = eri->buffer();
        
        int nP = auxiliary_->shell(P).nfunction();
        int oP = auxiliary_->shell(P).function_index();

        for (int R = 0; R < primary->nshell(); R++) {
            for (int S = 0; S <= R; S++) {
                int nR = primary->shell(R).nfunction();
                int oR = primary->shell(R).function_index();
                int nS = primary->shell(S).nfunction();
                int oS = primary->shell(S).function_index();
               
                eri->compute_shell(P,0,R,S);
            
                for (int p = 0; p < nP; p++) {
                    for (int r = 0; r < nR; r++) {
                        for (int s = 0; s < nS; s++) {
                            Dmnp[p + oP][(r + oR) * nso + (s + oS)] = 
                            Dmnp[p + oP][(s + oS) * nso + (r + oR)] = 
                            buffer[p*nR*nS + r*nS + s];
                        }            
                    }            
                }            
            }
        } 
    }  

    Dmn->print();

    // First-half transform
    C_DGEMM('N','N',naux*nso,ni,nso,1.0,Dmnp[0],nso,Cp[0],nmo,0.0,Dmjp[0],ni);
    
    Dmj->print();

    // Second-half transform
    #pragma omp parallel for num_threads(nthread)
    for (int P = 0; P < naux; P++) {
        C_DGEMM('T','N',ni,ni,nso,1.0,Dmjp[P],ni,Cp[0],nmo,0.0,Dijp[P],ni);
    }

    Dij->print();
    J->print();    

    // Fitting
    for (int ij = 0; ij < ni * ni; ij += naux) {
        int ncol = (ij + naux > (ni * ni) ? (ni * ni - ij) : naux);
        for (int P = 0; P < naux; P++) {
            ::memcpy((void*) Tp[P], (void*) &Dijp[P][ij], sizeof(double) * ncol);
        }
        C_DGEMM('N','N',naux,ncol,naux,1.0,Tp[0],naux,Jp[0],naux,0.0,&Dijp[0][ij],ni*ni);
    }

    Dij->print();

    return Dij;
}
void LHF::setup_V()
{
    EXX_ = 0.0;
    boost::shared_ptr<DFTGrid> grid = potential_->grid(); 
    V_X_ = SharedVector(new Vector("V_X", grid->npoints()));  
    V_S_ = SharedVector(new Vector("V_S", grid->npoints()));  
    V_C_ = SharedVector(new Vector("V_C", grid->npoints()));  
    V_ = SharedMatrix(new Matrix("V", KS::basisset_->nbf(), KS::basisset_->nbf()));

    // => Initial LDA Potential <= //

    // Globals/Grid
    double** Vp = V_->pointer();
    double* V_X = V_X_->pointer();
    boost::shared_ptr<RKSFunctions> properties = (static_cast<RV*>(potential_.get()))->properties(); 
    properties->reset_pointers(Da_, Ca_subset("SO","OCC"));
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    
    // Scratch
    int max_functions = grid->max_functions(); 
    int max_points = grid->max_points();
    SharedMatrix V_local(new Matrix("V Temp", max_functions, max_functions));
    double** V2p = V_local->pointer();
    SharedMatrix T_local(new Matrix("T Temp", max_points, max_functions));
    double** Tp = T_local->pointer();

    // Slater constant
    double K0 = 3.0 * pow(3.0 / (4.0 * M_PI), 1.0/3.0);

    unsigned long int offset = 0L;
    for (int Q = 0; Q < blocks.size(); Q++) {

        // Current block 
        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double *x = block->x();
        double *y = block->y();
        double *z = block->z();
        double *w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Basis points
        properties->computeProperties(block);
        double** phi = properties->basis_value("PHI")->pointer();
        double *rho_a = properties->property_value("RHO_A")->pointer();

        // Potential
        double* v_x = &V_X[offset];

        // Scalings 
        for (int P = 0; P < npoints; P++) {
            // Slater potential is 4/3 K_0 rho^1/3 
            v_x[P] -= 4.0 / 3.0 * K0 * pow(rho_a[P],1.0/3.0);
            // Form T_n^P = w_P v_P phi_n^P
            ::memset(static_cast<void*>(Tp[P]),'\0',nlocal*sizeof(double));
            C_DAXPY(nlocal,0.5 * v_x[P] * w[P], phi[P], 1, Tp[P], 1); 
            // Form E_XX
            EXX_ -= w[P] * K0 * pow(rho_a[P],4.0/3.0);
        }
         
        // Form V += phi_m^P T_n^P
        C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phi[0],max_functions,Tp[0],max_functions,0.0,V2p[0],max_functions);

        // Symmetrization (V is Hermitian)
        for (int m = 0; m < nlocal; m++) {
            for (int n = 0; n <= m; n++) {
                V2p[m][n] = V2p[n][m] = V2p[m][n] + V2p[n][m]; 
            }
        } 

        // Unpacking
        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
                Vp[mg][ng] += V2p[ml][nl];
                Vp[ng][mg] += V2p[ml][nl];
            }
            Vp[mg][mg] += V2p[ml][ml];
        }
         
        offset += npoints;
    }
}
void LHF::build_V_S()
{
    // This is gonna suck
    V_S_->zero();

    // Global sizing
    int ni = nalphapi_[0];      
    int nso = nsopi_[0];
    int nmo = nmopi_[0];
    int naux = auxiliary_->nbf();

    // Grid info
    boost::shared_ptr<DFTGrid> grid = potential_->grid(); 
    boost::shared_ptr<RKSFunctions> properties = (static_cast<RV*>(potential_.get()))->properties(); 
    properties->reset_pointers(Da_, Ca_subset("SO","OCC"));
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    int max_functions = grid->max_functions(); 
    int max_points = grid->max_points();

    // AO -> MO transform registers
    SharedMatrix P(new Matrix("Phi_i",max_points,ni));
    SharedMatrix R(new Matrix("V Phi_i",max_points,ni));
    SharedMatrix C2(new Matrix("C2",max_functions,ni));
    double** Pp = P->pointer();
    double** Rp = R->pointer();
    double** Cp = Ca_->pointer();
    double** C2p = C2->pointer();

    // Electrostatic Potential Registers
    SharedMatrix Vij(new Matrix("V_ij^P", max_points, ni * ni));
    SharedMatrix QAP(new Matrix("Q", max_points, naux));
    double** Vijp = Vij->pointer();
    double** Qp = QAP->pointer();

    // D_ij^P = (ij|A)(A|B)^-1
    SharedMatrix Dij = build_Dij();

    Dij->print();

    double** Dijp = Dij->pointer();

    // Electrostatic integral objects (threaded)
    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(auxiliary_, BasisSet::zero_ao_basis_set(), BasisSet::zero_ao_basis_set(), BasisSet::zero_ao_basis_set()));
    std::vector<boost::shared_ptr<PseudospectralInt> > ints;
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_get_max_threads();
    #endif 
    for (int i = 0; i < nthread; i++) {
        ints.push_back(boost::shared_ptr<PseudospectralInt>(static_cast<PseudospectralInt*>(factory->ao_pseudospectral())));
    }

    // Target
    double* V_S = V_S_->pointer();

    unsigned long int offset = 0L;
    for (int Q = 0; Q < blocks.size(); Q++) {

        // Current block 
        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double *x = block->x();
        double *y = block->y();
        double *z = block->z();
        double *w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Basis points
        properties->computeProperties(block);
        double** phi = properties->basis_value("PHI")->pointer();

        // AO -> MO transform
        for (int ml = 0; ml < function_map.size(); ml++) {
            int mg = function_map[ml];
            ::memcpy((void*) C2p[ml], (void*) Cp[mg], sizeof(double) * ni); 
        } 
        C_DGEMM('N','N',npoints,ni,nlocal,1.0,phi[0],max_functions,C2p[0],ni,0.0,Pp[0],ni);

        // Form Q integrals (N^2, high prefactor)
        #pragma omp parallel for num_threads(nthread)
        for (int A = 0; A < npoints; A++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            // Set to do the current point
            boost::shared_ptr<PseudospectralInt> eri = ints[thread];
            const double* buffer = eri->buffer();
            eri->set_point(x[A],y[A],z[A]);
            double* reg = Qp[A];

            for (int B = 0; B < auxiliary_->nshell(); B++) {
                int nB = auxiliary_->shell(B).nfunction();
                int oB = auxiliary_->shell(B).function_index();
                eri->compute_shell(B,0);
                ::memcpy((void*) &reg[oB], (void*) buffer, sizeof(double) * nB);        
            }
        }        

        // Form Vijp (N^4 Rate-Limiter)
        C_DGEMM('N','N',npoints,ni*ni,naux,1.0,Qp[0],naux,Dijp[0],ni*ni,0.0,Vijp[0],ni*ni);

        // Potential
        C_DGEMV('N',npoints * ni,ni,1.0,Vijp[0],ni,Pp[0],1,0.0,Rp[0],1);

        // Summation
        double* v_s = &V_S[offset];
        for (int A = 0; A < npoints; A++) {
            v_s[A] = C_DDOT(ni,Pp[A],1,Rp[A],1);
        }
         
        offset += npoints;
    }
}
void LHF::build_V_C()
{
    // This only partially sucks
    V_C_->zero();

    // Global sizing
    int ni = nalphapi_[0];      
    int nso = nsopi_[0];
    int nmo = nmopi_[0];

    // Build the matrix Delta = V_LHF - V_NL from last iteration's values 
    SharedMatrix Delta3(V_->clone());
    Delta3->add(K_);
    SharedMatrix Delta2(new Matrix("Temp", ni, nso));
    SharedMatrix Delta(new Matrix("V_X - V_NL", ni, ni));
    C_DGEMM('T','N',ni,nso,nso,1.0,Ca_->pointer()[0],nmo,Delta3->pointer()[0],nso,0.0,Delta2->pointer()[0],nso);
    C_DGEMM('N','N',ni,ni,nso,1.0,Delta2->pointer()[0],nso,Ca_->pointer()[0],nmo,0.0,Delta->pointer()[0],ni);
    double** Dp = Delta->pointer();
    double   Df = Dp[ni - 1][ni - 1];
    
    // Grid info
    boost::shared_ptr<DFTGrid> grid = potential_->grid(); 
    boost::shared_ptr<RKSFunctions> properties = (static_cast<RV*>(potential_.get()))->properties(); 
    properties->reset_pointers(Da_, Ca_subset("SO","OCC"));
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    int max_functions = grid->max_functions(); 
    int max_points = grid->max_points();

    // AO -> MO transform registers
    SharedMatrix P(new Matrix("Phi_i",max_points,ni));
    SharedMatrix R(new Matrix("V Phi_i",max_points,ni));
    SharedMatrix C2(new Matrix("C2",max_functions,ni));
    double** Pp = P->pointer();
    double** Rp = R->pointer();
    double** Cp = Ca_->pointer();
    double** C2p = C2->pointer();

    // Target
    double* V_C = V_C_->pointer();

    unsigned long int offset = 0L;
    for (int Q = 0; Q < blocks.size(); Q++) {

        // Current block 
        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double *x = block->x();
        double *y = block->y();
        double *z = block->z();
        double *w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Basis points
        properties->computeProperties(block);
        double** phi = properties->basis_value("PHI")->pointer();

        // AO -> MO transform
        for (int ml = 0; ml < function_map.size(); ml++) {
            int mg = function_map[ml];
            ::memcpy((void*) C2p[ml], (void*) Cp[mg], sizeof(double) * ni); 
        } 
        C_DGEMM('N','N',npoints,ni,nlocal,1.0,phi[0],max_functions,C2p[0],ni,0.0,Pp[0],ni);

        // Potential
        C_DGEMM('N','N',npoints,ni,ni,1.0,Pp[0],ni,Dp[0],ni,0.0,Rp[0],ni);

        // Summation, including HOMO restriction
        double* v_c = &V_C[offset];
        for (int A = 0; A < npoints; A++) {
            v_c[A] = C_DDOT(ni,Pp[A],1,Rp[A],1) - Df * Pp[A][ni-1] * Pp[A][ni - 1];   
        }
         
        offset += npoints;
    }
}
void LHF::build_V_X()
{
    // Globals/Grid
    EXX_ = 0.0;
    V_->zero();
    double** Vp = V_->pointer();
    double* V_X = V_X_->pointer();
    boost::shared_ptr<DFTGrid> grid = potential_->grid(); 
    boost::shared_ptr<RKSFunctions> properties = (static_cast<RV*>(potential_.get()))->properties(); 
    properties->reset_pointers(Da_, Ca_subset("SO","OCC"));
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    
    // Scratch
    int max_functions = grid->max_functions(); 
    int max_points = grid->max_points();
    SharedMatrix V_local(new Matrix("V Temp", max_functions, max_functions));
    double** V2p = V_local->pointer();
    SharedMatrix T_local(new Matrix("T Temp", max_points, max_functions));
    double** Tp = T_local->pointer();

    // Add V_S + V_C to make V_S (sans 2 / rho)
    ::memcpy((void*) V_X, (void*) V_S_->pointer(), sizeof(double) * grid->npoints());
    C_DAXPY(grid->npoints(), 1.0, V_C_->pointer(), 1, V_X,1);

    unsigned long int offset = 0L;
    for (int Q = 0; Q < blocks.size(); Q++) {

        // Current block 
        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double *x = block->x();
        double *y = block->y();
        double *z = block->z();
        double *w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        // Basis points
        properties->computeProperties(block);
        double** phi = properties->basis_value("PHI")->pointer();
        double *rho_a = properties->property_value("RHO_A")->pointer();

        // Potential
        double* v_x = &V_X[offset];

        // Scalings 
        for (int P = 0; P < npoints; P++) {
            // Apply 2 / rho = 1 / rho_s
            v_x[P] *= 1.0 / rho_a[P];
            // Form T_n^P = w_P v_P phi_n^P
            ::memset(static_cast<void*>(Tp[P]),'\0',nlocal*sizeof(double));
            C_DAXPY(nlocal,0.5 * v_x[P] * w[P], phi[P], 1, Tp[P], 1); 
            // Form E_XX
            EXX_ += 0.5 * w[P] * v_x[P] * rho_a[P];
        }
         
        // Form V += phi_m^P T_n^P
        C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phi[0],max_functions,Tp[0],max_functions,0.0,V2p[0],max_functions);

        // Symmetrization (V is Hermitian)
        for (int m = 0; m < nlocal; m++) {
            for (int n = 0; n <= m; n++) {
                V2p[m][n] = V2p[n][m] = V2p[m][n] + V2p[n][m]; 
            }
        } 

        // Unpacking
        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
                Vp[mg][ng] += V2p[ml][nl];
                Vp[ng][mg] += V2p[ml][nl];
            }
            Vp[mg][mg] += V2p[ml][ml];
        }
         
        offset += npoints;
    }
     
}
void LHF::form_V()
{
    if (!V_X_) {
        // Build V from LDA kernel
        setup_V();
    } else {
        // Form the Slater Potential on the grid
        build_V_S(); 
        // Form the local correction on the grid
        build_V_C();
        // Build V_X = <m|v_x|n>
        build_V_X();
    }
}
void LHF::form_G()
{
    // From the V matrix
    form_V();

    // Push the C matrix on
    std::vector<SharedMatrix> & C = jk_->C_left();
    C.clear();
    C.push_back(Ca_subset("SO", "OCC"));
    
    // Run the JK object
    jk_->compute();

    // Pull the J and K matrices off
    const std::vector<SharedMatrix> & J = jk_->J();
    const std::vector<SharedMatrix> & K = jk_->K();

    J_ = J[0];
    K_ = K[0];

    // F = 2 J + V
    G_->zero();
    G_->add(J_);
    G_->scale(2.0);
    G_->add(V_);
}
double LHF::compute_E()
{
    // E_LHF = 2.0 D*H + 2.0 D*J - E_X
    double one_electron_E = 2.0*D_->vector_dot(H_);
    double coulomb_E = 2.0*D_->vector_dot(J_);
    double X_E = EXX_; 

    double Etotal = 0.0;
    Etotal += nuclearrep_;
    Etotal += one_electron_E;
    Etotal += coulomb_E;
    Etotal += X_E; 

    if (debug_) {
        fprintf(outfile, "   => Energetics <=\n\n");
        fprintf(outfile, "    Nuclear Repulsion Energy = %24.14f\n", nuclearrep_);
        fprintf(outfile, "    One-Electron Energy =      %24.14f\n", one_electron_E);
        fprintf(outfile, "    Coulomb Energy =           %24.14f\n", coulomb_E);
        fprintf(outfile, "    X Functional Energy  =     %24.14f\n", X_E); 
    }

    return Etotal;
}

}} // Namespaces
