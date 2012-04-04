#include <libmints/mints.h>
#include <libfunctional/superfunctional.h>
#include <libqt/qt.h>
#include <psi4-dec.h>

#include "cubature.h"
#include "points.h"
#include "v.h"

#include <sstream>

using namespace psi;

namespace psi {

VBase::VBase(boost::shared_ptr<SuperFunctional> functional,
    boost::shared_ptr<BasisSet> primary,
    Options& options):
    functional_(functional), primary_(primary), options_(options)
{
    common_init();
}
VBase::~VBase()
{
}
void VBase::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");
}
boost::shared_ptr<VBase> VBase::build_V(Options& options, const std::string& type)
{
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> primary = BasisSet::construct(parser, Process::environment.molecule(), "BASIS");

    int depth = 1; // By default, do first partials of the kernel
    if (type == "RK" || type == "UK")
        depth = 2;

    int block_size = options.get_int("DFT_BLOCK_MAX_POINTS");
    boost::shared_ptr<SuperFunctional> functional = SuperFunctional::current(options,block_size,depth);

    boost::shared_ptr<VBase> v;
    if (type == "RV") {
        v = boost::shared_ptr<VBase>(new RV(functional,primary,options));
    } else if (type == "UV") {
        v = boost::shared_ptr<VBase>(new UV(functional,primary,options));
    } else if (type == "RK") {
        v = boost::shared_ptr<VBase>(new RK(functional,primary,options));
    } else if (type == "UK") { 
        v = boost::shared_ptr<VBase>(new UK(functional,primary,options));
    } else {
        throw PSIEXCEPTION("V: V type is not recognized");    
    }

    return v;
}
void VBase::compute_D()
{
    // Allocate D if needed
    if (D_.size() != C_.size()) {
        D_.clear();
        for (int A = 0; A < C_.size(); A++) {
            std::stringstream ss;
            ss << "D (SO) " << A;
            D_.push_back(SharedMatrix(new Matrix(ss.str(), C_[A]->rowspi(), C_[A]->rowspi())));
        }
    }
    
    for (int A = 0; A < C_.size(); A++) {
        SharedMatrix C = C_[A];    
        SharedMatrix D = D_[A];    
        D->zero();
        for (int h = 0; h < C->nirrep(); h++) {
            int nso = C->rowspi()[h];
            int nmo = C->colspi()[h];
            if (!nso || !nmo) continue;
            double** Cp = C->pointer(h);
            double** Dp = D->pointer(h);
            C_DGEMM('N','T',nso,nso,nmo,1.0,Cp[0],nmo,Cp[0],nmo,0.0,Dp[0],nso);
        }
    }

    // Allocate P_SO_ if needed
    if (P_.size()) {
        bool same = true;
        if (P_SO_.size() != P_.size()) {
            same = false;   
        } else {
            for (int A = 0; A < P_.size(); A++) {
                if (P_[A]->symmetry() != P_SO_[A]->symmetry())
                    same = false;
            }
        }

        if (!same) {
            P_SO_.clear();
            for (int A = 0; A < P_.size(); A++) {
                std::stringstream ss;
                ss << "P (SO) " << A;
                P_SO_.push_back(SharedMatrix(new Matrix(ss.str(),C_[0]->rowspi(),C_[0]->rowspi(),P_[A]->symmetry())));
            }
        } 

        int maxi = Caocc_[0]->max_ncol();
        if (Caocc_.size() > 1) 
            maxi = (Caocc_[1]->max_ncol() > maxi ? Caocc_[1]->max_ncol() : maxi);
        double* temp = new double[maxi * (ULI) Caocc_[0]->max_nrow()];
        for (int A = 0; A < P_.size(); A++) {
            SharedMatrix Cl = Caocc_[A % Caocc_.size()];
            SharedMatrix Cr = Cavir_[A % Cavir_.size()];
            SharedMatrix P = P_[A];
            SharedMatrix P_SO = P_SO_[A];
            int symm = P->symmetry();
            for (int h = 0; h < P->nirrep(); h++) {
                int nl = Cl->rowspi()[h];
                int nr = Cr->rowspi()[h^symm];
                int ni = Cl->colspi()[h];
                int na = Cr->colspi()[h^symm];
                if (!nl || !nr || !ni || !na) continue;
                double** Clp   = Cl->pointer(h);
                double** Crp   = Cr->pointer(h^symm);
                double** Pp    = P->pointer(h);
                double** P_SOp = P_SO->pointer(h);

                C_DGEMM('N','T',ni,nr,na,1.0,Pp[0],na,Crp[0],na,0.0,temp,nr);
                C_DGEMM('N','N',nl,nr,ni,1.0,Clp[0],ni,temp,nr,0.0,P_SOp[0],nr);
            }
        }    
        delete[] temp;
    }    
}   
void VBase::USO2AO()
{
    // Build AO2USO matrix, if needed
    if (!AO2USO_ && (C_[0]->nirrep() != 1)) {
        boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(primary_,primary_,primary_,primary_));
        boost::shared_ptr<PetiteList> pet(new PetiteList(primary_, integral));
        AO2USO_ = SharedMatrix(pet->aotoso());
    }

    // Allocate V if needed
    if (P_.size()) {
        bool same = true;
        if (V_.size() != P_.size()) {
            same = false;   
        } else {
            for (int A = 0; A < P_.size(); A++) {
                if (P_[A]->symmetry() != V_[A]->symmetry())
                    same = false;
            }
        } 
        if (!same) {
            V_.clear();
            for (int A = 0; A < P_.size(); A++) {
                std::stringstream ss1;
                ss1 << "V (SO) " << A;
                V_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss1.str(), C_[0]->rowspi(), C_[0]->rowspi(),P_[A]->symmetry())));
            } 
        }
    } else {
        if (V_.size() != C_.size()) {
            V_.clear();
            for (int A = 0; A < C_.size(); A++) {
                std::stringstream ss1;
                ss1 << "V (SO) " << A;
                V_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss1.str(), C_[0]->rowspi(), C_[0]->rowspi())));
            } 
            
        }
    }

    // If C1, just assign pointers
    if (C_[0]->nirrep() == 1) {
        C_AO_ = C_;
        D_AO_ = D_;
        V_AO_ = V_;
        // Clear V_AO_ out
        for (int A = 0; A < V_AO_.size(); A++) {
            V_AO_[A]->zero();
        }
        P_AO_ = P_SO_;
        for (int A = 0; A < P_AO_.size(); A++) {
            P_AO_[A]->hermitivitize();
        }
        return;    
    }

    if (P_.size()) {
        if (V_.size() != P_.size()) {
            V_AO_.clear();
            P_AO_.clear();
            for (int A = 0; A < P_.size(); A++) {
                std::stringstream ss2;
                ss2 << "V (AO) " << A;
                V_AO_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss2.str(), C_[0]->nrow(), C_[0]->nrow())));
                std::stringstream ss3;
                ss3 << "P (AO) " << A;
                P_AO_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss3.str(), C_[0]->nrow(), C_[0]->nrow())));
            } 
        }
    } else {
        if (V_AO_.size() != C_.size()) {
            V_AO_.clear();
            P_AO_.clear();
            for (int A = 0; A < C_.size(); A++) {
                std::stringstream ss2;
                ss2 << "V (AO) " << A;
                V_AO_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss2.str(), C_[0]->nrow(), C_[0]->nrow())));
                std::stringstream ss3;
                ss3 << "P (AO) " << A;
                P_AO_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss3.str(), C_[0]->nrow(), C_[0]->nrow())));
            } 
            
        }
    }

    // Clear V_AO_ out
    for (int A = 0; A < V_AO_.size(); A++) {
        V_AO_[A]->zero();
    }

    // Allocate C_AO/D_AO if needed
    if (C_AO_.size() != C_.size()) {
        C_AO_.clear();
        D_AO_.clear();
        for (int A = 0; A < C_.size(); A++) {
            std::stringstream ss1;
            ss1 << "C (AO) " << A;
            C_AO_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss1.str(), C_[0]->nrow(), C_[0]->ncol())));
            std::stringstream ss2;
            ss2 << "D (AO) " << A;
            D_AO_.push_back(boost::shared_ptr<Matrix>(new Matrix(ss2.str(), C_[0]->nrow(), C_[0]->nrow())));
        }    
    }

    // C_AO (Order is not important, just KE Density)
    for (int A = 0; A < C_.size(); A++) {
        SharedMatrix C = C_[A];
        SharedMatrix C_AO = C_AO_[A];
        int offset = 0;
        for (int h = 0; h < C->nirrep(); h++) {
            int nocc = C->colspi()[h];
            int nso = C->rowspi()[h];
            int nao = AO2USO_->rowspi()[h];
            int nmo = C_AO->colspi()[0];
            if (!nocc || !nso || !nao || !nmo) continue;
            double** Cp = C->pointer(h);
            double** C_AOp = C_AO->pointer(0);
            double** Up = AO2USO_->pointer(h);
            C_DGEMM('N','N',nao,nocc,nso,1.0,Up[0],nso,Cp[0],nocc,0.0,&C_AOp[0][offset],nmo);
            offset += nocc;
        }
    }
    
    // D_AO
    double* temp = new double[AO2USO_->max_nrow() * (ULI) AO2USO_->max_ncol()];
    for (int A = 0; A < D_AO_.size(); A++) {
        SharedMatrix D = D_[A];
        SharedMatrix D_AO = D_AO_[A];
        D_AO->zero();
        for (int h = 0; h < D->nirrep(); h++) {
            int symm = D->symmetry();
            int nao = AO2USO_->rowspi()[h];
            int nsol = AO2USO_->colspi()[h];
            int nsor = AO2USO_->colspi()[h^symm];
            if (!nao || !nsol || !nsor) continue;
            double** Dp = D->pointer(h);
            double** D_AOp = D_AO->pointer();
            double** Ulp = AO2USO_->pointer(h);
            double** Urp = AO2USO_->pointer(h^symm);
            C_DGEMM('N','N',nao,nsor,nsol,1.0,Ulp[0],nsol,Dp[0],nsor,0.0,temp,nsor);
            C_DGEMM('N','T',nao,nao,nsor,1.0,temp,nsor,Urp[0],nsor,1.0,D_AOp[0],nao);
        }
    }    

    // P_AO
    for (int A = 0; A < P_SO_.size(); A++) {
        SharedMatrix D = P_SO_[A];
        SharedMatrix D_AO = P_AO_[A];
        D_AO->zero();
        for (int h = 0; h < D->nirrep(); h++) {
            int symm = D->symmetry();
            int nao = AO2USO_->rowspi()[h];
            int nsol = AO2USO_->colspi()[h];
            int nsor = AO2USO_->colspi()[h^symm];
            if (!nao || !nsol || !nsor) continue;
            double** Dp = D->pointer(h);
            double** D_AOp = D_AO->pointer();
            double** Ulp = AO2USO_->pointer(h);
            double** Urp = AO2USO_->pointer(h^symm);
            C_DGEMM('N','N',nao,nsor,nsol,1.0,Ulp[0],nsol,Dp[0],nsor,0.0,temp,nsor);
            C_DGEMM('N','T',nao,nao,nsor,1.0,temp,nsor,Urp[0],nsor,1.0,D_AOp[0],nao);
        }
        D_AO->hermitivitize();
    }
    delete[] temp;
}
void VBase::AO2USO()
{
    if (C_[0]->nirrep() == 1) {
        // V_ is already assigned
        return;    
    }

    double* temp = new double[AO2USO_->max_nrow() * (ULI) AO2USO_->max_ncol()];
    for (int A = 0; A < V_AO_.size(); A++) {
        SharedMatrix V = V_[A];
        SharedMatrix V_AO = V_AO_[A];
        for (int h = 0; h < V->nirrep(); h++) {
            int symm = V->symmetry();
            int nao = AO2USO_->rowspi()[h];
            int nsol = AO2USO_->colspi()[h];
            int nsor = AO2USO_->colspi()[h^symm];
            if (!nao || !nsol || !nsor) continue;
            double** Vp = V->pointer(h);
            double** V_AOp = V_AO->pointer();
            double** Ulp = AO2USO_->pointer(h);
            double** Urp = AO2USO_->pointer(h^symm);
            C_DGEMM('N','N',nao,nsor,nao,1.0,V_AOp[0],nao,Urp[0],nsor,0.0,temp,nsor);
            C_DGEMM('T','N',nsol,nsor,nao,1.0,Ulp[0],nsol,temp,nsor,0.0,Vp[0],nsor);
        }
    }    
    delete[] temp;
}
void VBase::initialize()
{
    timer_on("V: Grid");
    grid_ = boost::shared_ptr<DFTGrid>(new DFTGrid(primary_->molecule(),primary_,options_));
    timer_off("V: Grid");
}
void VBase::compute()
{
    timer_on("V: D");
    compute_D();
    timer_off("V: D");

    timer_on("V: USO2AO");
    USO2AO();
    timer_off("V: USO2AO");

    timer_on("V: V");
    compute_V();
    timer_off("V: V");

    timer_on("V: AO2USO");
    AO2USO();
    timer_off("V: AO2USO");
}
SharedMatrix VBase::compute_gradient()
{
    throw PSIEXCEPTION("VBase: gradient not implemented for this V instance.");
}
void VBase::finalize()
{
    grid_.reset();
}
void VBase::print_header() const
{
    fprintf(outfile, "  ==> DFT Potential <==\n\n");
    functional_->print(outfile, print_);  
    grid_->print(outfile,print_);
}

RV::RV(boost::shared_ptr<SuperFunctional> functional,
    boost::shared_ptr<BasisSet> primary,
    Options& options) : VBase(functional,primary,options)
{
}
RV::~RV()
{
}
void RV::initialize()
{
    VBase::initialize();
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions(); 
    properties_ = boost::shared_ptr<PointFunctions>(new RKSFunctions(primary_,max_points,max_functions));
    properties_->set_ansatz(functional_->ansatz());
}
void RV::finalize()
{
    properties_.reset();
    VBase::finalize();
}
void RV::print_header() const
{
    VBase::print_header();
}
void RV::compute_V()
{
    if ((D_AO_.size() != 1) || (V_AO_.size() != 1) || (C_AO_.size() != 1))
        throw PSIEXCEPTION("V: RKS should have only one D/C/V Matrix"); 

    // Setup the pointers
    SharedMatrix D_AO = D_AO_[0];
    SharedMatrix C_AO = C_AO_[0];
    SharedMatrix V_AO = V_AO_[0];
    properties_->set_pointers(D_AO,C_AO);

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions(); 
    int max_points = grid_->max_points();

    // Local/global V matrices
    SharedMatrix V_local(new Matrix("V Temp", max_functions, max_functions));
    double** V2p = V_local->pointer();
    double** Vp = V_AO->pointer();

    // Scratch
    std::vector<SharedMatrix> scratch = properties_->scratch();
    SharedMatrix T_local = scratch[0];
    double** Tp = T_local->pointer();

    // Traverse the blocks of points
    double functionalq = 0.0;
    double rhoaq       = 0.0;
    double rhoaxq      = 0.0;
    double rhoayq      = 0.0;
    double rhoazq      = 0.0;

    boost::shared_ptr<Vector> QT(new Vector("Quadrature Temp", max_points));
    double *restrict QTp = QT->pointer();
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid_->blocks();

    for (int Q = 0; Q < blocks.size(); Q++) {

        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double *restrict x = block->x();
        double *restrict y = block->y();
        double *restrict z = block->z();
        double *restrict w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        timer_on("Properties");
        properties_->compute_points(block);
        timer_off("Properties");
        timer_on("Functional");
        std::map<std::string, SharedVector>& vals = functional_->compute_functional(properties_->point_values(), npoints); 
        timer_off("Functional");

        if (debug_ > 4) {
            block->print(outfile, debug_);
            properties_->print(outfile, debug_);
        }

        timer_on("V_XC");
        double** phi = properties_->basis_value("PHI")->pointer();
        double *restrict rho_a = properties_->point_value("RHO_A")->pointer();
        double *restrict zk = vals["V"]->pointer(); 
        double *restrict v_rho_a = vals["V_RHO_A"]->pointer();

        // => Quadrature values <= //
        functionalq += C_DDOT(npoints,w,1,zk,1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_a[P];
        }
        rhoaq       += C_DDOT(npoints,w,1,rho_a,1);
        rhoaxq      += C_DDOT(npoints,QTp,1,x,1);
        rhoayq      += C_DDOT(npoints,QTp,1,y,1);
        rhoazq      += C_DDOT(npoints,QTp,1,z,1);

        // => LSDA contribution (symmetrized) <= //
        timer_on("LSDA");
        for (int P = 0; P < npoints; P++) {
            ::memset(static_cast<void*>(Tp[P]),'\0',nlocal*sizeof(double));
            C_DAXPY(nlocal,0.5 * v_rho_a[P] * w[P], phi[P], 1, Tp[P], 1); 
        }
        timer_off("LSDA");
        
        // => GGA contribution (symmetrized) <= // 
        if (ansatz >= 1) {
            timer_on("GGA");
            double** phix = properties_->basis_value("PHI_X")->pointer();
            double** phiy = properties_->basis_value("PHI_Y")->pointer();
            double** phiz = properties_->basis_value("PHI_Z")->pointer();
            double *restrict rho_ax = properties_->point_value("RHO_AX")->pointer();
            double *restrict rho_ay = properties_->point_value("RHO_AY")->pointer();
            double *restrict rho_az = properties_->point_value("RHO_AZ")->pointer();
            double *restrict v_sigma_aa = vals["V_GAMMA_AA"]->pointer(); 
            double *restrict v_sigma_ab = vals["V_GAMMA_AB"]->pointer(); 

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_ax[P] + v_sigma_ab[P] * rho_ax[P]), phix[P], 1, Tp[P], 1); 
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_ay[P] + v_sigma_ab[P] * rho_ay[P]), phiy[P], 1, Tp[P], 1); 
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_az[P] + v_sigma_ab[P] * rho_az[P]), phiz[P], 1, Tp[P], 1); 
            }        
            timer_off("GGA");
        }

        // Single GEMM slams GGA+LSDA together (man but GEM's hot!)
        timer_on("LSDA");
        C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phi[0],max_functions,Tp[0],max_functions,0.0,V2p[0],max_functions);

        // Symmetrization (V is Hermitian)
        for (int m = 0; m < nlocal; m++) {
            for (int n = 0; n <= m; n++) {
                V2p[m][n] = V2p[n][m] = V2p[m][n] + V2p[n][m]; 
            }
        } 
        timer_off("LSDA");

        // => Meta contribution <= //
        if (ansatz >= 2) {
            timer_on("Meta");
            double** phix = properties_->basis_value("PHI_X")->pointer();
            double** phiy = properties_->basis_value("PHI_Y")->pointer();
            double** phiz = properties_->basis_value("PHI_Z")->pointer();
            double *restrict v_tau_a = vals["V_TAU_A"]->pointer(); 
            
            // \nabla x
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tp[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_a[P] * w[P], phix[P], 1, Tp[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phix[0],max_functions,Tp[0],max_functions,1.0,V2p[0],max_functions);

            // \nabla y
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tp[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_a[P] * w[P], phiy[P], 1, Tp[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phiy[0],max_functions,Tp[0],max_functions,1.0,V2p[0],max_functions);

            // \nabla z
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tp[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_a[P] * w[P], phiz[P], 1, Tp[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phiz[0],max_functions,Tp[0],max_functions,1.0,V2p[0],max_functions);
            timer_off("Meta");
        }       
 
        // => Unpacking <= //
        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
                Vp[mg][ng] += V2p[ml][nl];
                Vp[ng][mg] += V2p[ml][nl];
            }
            Vp[mg][mg] += V2p[ml][ml];
        }
        timer_off("V_XC");
    } 
   
    quad_values_["FUNCTIONAL"] = functionalq;
    quad_values_["RHO_A"]      = rhoaq; 
    quad_values_["RHO_AX"]     = rhoaxq; 
    quad_values_["RHO_AY"]     = rhoayq; 
    quad_values_["RHO_AZ"]     = rhoazq; 
    quad_values_["RHO_B"]      = rhoaq; 
    quad_values_["RHO_BX"]     = rhoaxq; 
    quad_values_["RHO_BY"]     = rhoayq; 
    quad_values_["RHO_BZ"]     = rhoazq; 
 
    if (debug_) {
        fprintf(outfile, "   => Numerical Integrals <=\n\n");
        fprintf(outfile, "    Functional Value:  %24.16E\n",quad_values_["FUNCTIONAL"]);
        fprintf(outfile, "    <\\rho_a>        :  %24.16E\n",quad_values_["RHO_A"]);
        fprintf(outfile, "    <\\rho_b>        :  %24.16E\n",quad_values_["RHO_B"]);
        fprintf(outfile, "    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",quad_values_["RHO_AX"],quad_values_["RHO_AY"],quad_values_["RHO_AZ"]);
        fprintf(outfile, "    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",quad_values_["RHO_BX"],quad_values_["RHO_BY"],quad_values_["RHO_BZ"]);
    }
}
SharedMatrix RV::compute_gradient()
{
    compute_D();
    USO2AO();

    if ((D_AO_.size() != 1) && (C_AO_.size() != 1))
        throw PSIEXCEPTION("V: RKS should have only one D/C/V Matrix"); 

    // Build the target gradient Matrix
    int natom = primary_->molecule()->natom();
    SharedMatrix G(new Matrix("XC Gradient", natom,3));
    double** Gp = G->pointer();

    // Set Hessian derivative level in properties
    int old_deriv = properties_->deriv(); 
    properties_->set_deriv((functional_->is_gga() ? 2 : 1));

    // Setup the pointers
    SharedMatrix D_AO = D_AO_[0];
    SharedMatrix C_AO = C_AO_[0];
    properties_->set_pointers(D_AO,C_AO);

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions(); 
    int max_points = grid_->max_points();

    // Scratch
    std::vector<SharedMatrix> scratch = properties_->scratch();
    SharedMatrix T_local = scratch[0];
    SharedMatrix U_local(T_local->clone());
    double** Tp = T_local->pointer();
    double** Up = U_local->pointer();
    std::vector<SharedMatrix> Dscratch = properties_->D_scratch();
    SharedMatrix D_local = Dscratch[0];
    double** Dp = D_local->pointer();

    // Traverse the blocks of points
    double functionalq = 0.0;
    double rhoaq       = 0.0;
    double rhoaxq      = 0.0;
    double rhoayq      = 0.0;
    double rhoazq      = 0.0;

    boost::shared_ptr<Vector> QT(new Vector("Quadrature Temp", max_points));
    double* QTp = QT->pointer();
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid_->blocks();

    for (int Q = 0; Q < blocks.size(); Q++) {

        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        timer_on("Properties");
        properties_->compute_points(block);
        timer_off("Properties");
        timer_on("Functional");
        std::map<std::string, SharedVector>& vals = functional_->compute_functional(properties_->point_values(), npoints); 
        timer_off("Functional");

        double** phi = properties_->basis_value("PHI")->pointer();
        double** phi_x = properties_->basis_value("PHI_X")->pointer();
        double** phi_y = properties_->basis_value("PHI_Y")->pointer();
        double** phi_z = properties_->basis_value("PHI_Z")->pointer();
        double* rho_a = properties_->point_value("RHO_A")->pointer();
        double* zk = vals["V"]->pointer(); 
        double* v_rho_a = vals["V_RHO_A"]->pointer();

        // => Quadrature values <= //
        functionalq += C_DDOT(npoints,w,1,zk,1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_a[P];
        }
        rhoaq       += C_DDOT(npoints,w,1,rho_a,1);
        rhoaxq      += C_DDOT(npoints,QTp,1,x,1);
        rhoayq      += C_DDOT(npoints,QTp,1,y,1);
        rhoazq      += C_DDOT(npoints,QTp,1,z,1);

        // => LSDA Contribution <= //
        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_a[P], phi[P], 1, Tp[P], 1);
        }
    
        // => GGA Contribution (Term 1) <= //
        if (functional_->is_gga()) {
            double* rho_ax = properties_->point_value("RHO_AX")->pointer();
            double* rho_ay = properties_->point_value("RHO_AY")->pointer();
            double* rho_az = properties_->point_value("RHO_AZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_ax[P]), phi_x[P], 1, Tp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_ay[P]), phi_y[P], 1, Tp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_az[P]), phi_z[P], 1, Tp[P], 1);
            }
    
        }

        // => Synthesis <= //
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,Tp[0],max_functions,Dp[0],max_functions,0.0,Up[0],max_functions);
    
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints,&Up[0][ml],max_functions,&phi_x[0][ml],max_functions);
            Gp[A][1] += C_DDOT(npoints,&Up[0][ml],max_functions,&phi_y[0][ml],max_functions);
            Gp[A][2] += C_DDOT(npoints,&Up[0][ml],max_functions,&phi_z[0][ml],max_functions);
        }          
        
        // => GGA Contribution (Term 2) <= //
        if (functional_->is_gga()) {
            double** phi_xx = properties_->basis_value("PHI_XX")->pointer();
            double** phi_xy = properties_->basis_value("PHI_XY")->pointer();
            double** phi_xz = properties_->basis_value("PHI_XZ")->pointer();
            double** phi_yy = properties_->basis_value("PHI_YY")->pointer();
            double** phi_yz = properties_->basis_value("PHI_YZ")->pointer();
            double** phi_zz = properties_->basis_value("PHI_ZZ")->pointer();
            double* rho_ax = properties_->point_value("RHO_AX")->pointer();
            double* rho_ay = properties_->point_value("RHO_AY")->pointer();
            double* rho_az = properties_->point_value("RHO_AZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();

            C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi[0],max_functions,Dp[0],max_functions,0.0,Up[0],max_functions);
            
            // x
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_ax[P]), Up[P], 1, Tp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xx[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xz[0][ml],max_functions);
            }          
            
            // y
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_ay[P]), Up[P], 1, Tp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_yy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
            }          
            
            // z
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_az[P]), Up[P], 1, Tp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_xz[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tp[0][ml],max_functions,&phi_zz[0][ml],max_functions);
            }          
    
        }
        
        // => Meta Contribution <= //
        if (functional_->is_meta()) {
            throw PSIEXCEPTION("RV: Meta Gradients not implemented");
        }

    } 
   
    quad_values_["FUNCTIONAL"] = functionalq;
    quad_values_["RHO_A"]      = rhoaq; 
    quad_values_["RHO_AX"]     = rhoaxq; 
    quad_values_["RHO_AY"]     = rhoayq; 
    quad_values_["RHO_AZ"]     = rhoazq; 
    quad_values_["RHO_B"]      = rhoaq; 
    quad_values_["RHO_BX"]     = rhoaxq; 
    quad_values_["RHO_BY"]     = rhoayq; 
    quad_values_["RHO_BZ"]     = rhoazq; 
 
    if (debug_) {
        fprintf(outfile, "   => XC Gradient: Numerical Integrals <=\n\n");
        fprintf(outfile, "    Functional Value:  %24.16E\n",quad_values_["FUNCTIONAL"]);
        fprintf(outfile, "    <\\rho_a>        :  %24.16E\n",quad_values_["RHO_A"]);
        fprintf(outfile, "    <\\rho_b>        :  %24.16E\n",quad_values_["RHO_B"]);
        fprintf(outfile, "    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",quad_values_["RHO_AX"],quad_values_["RHO_AY"],quad_values_["RHO_AZ"]);
        fprintf(outfile, "    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",quad_values_["RHO_BX"],quad_values_["RHO_BY"],quad_values_["RHO_BZ"]);
    }

    properties_->set_deriv(old_deriv);

    // RKS
    G->scale(2.0);

    return G;
}

UV::UV(boost::shared_ptr<SuperFunctional> functional,
    boost::shared_ptr<BasisSet> primary,
    Options& options) : VBase(functional,primary,options)
{
}
UV::~UV()
{
}
void UV::initialize()
{
    VBase::initialize();
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions(); 
    properties_ = boost::shared_ptr<PointFunctions>(new UKSFunctions(primary_,max_points,max_functions));
    properties_->set_ansatz(functional_->ansatz());
}
void UV::finalize()
{
    properties_.reset();
    VBase::finalize();
}
void UV::print_header() const
{
    VBase::print_header();
}
void UV::compute_V()
{
    if ((D_AO_.size() != 2) || (V_AO_.size() != 2) || (C_AO_.size() != 2))
        throw PSIEXCEPTION("V: UKS should have two D/C/V Matrices"); 

    // Setup the pointers
    SharedMatrix Da_AO = D_AO_[0];
    SharedMatrix Ca_AO = C_AO_[0];
    SharedMatrix Va_AO = V_AO_[0];
    SharedMatrix Db_AO = D_AO_[1];
    SharedMatrix Cb_AO = C_AO_[1];
    SharedMatrix Vb_AO = V_AO_[1];
    properties_->set_pointers(Da_AO,Ca_AO,Db_AO,Cb_AO);

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Local/global V matrices
    SharedMatrix Va_local(new Matrix("Va Temp", max_functions, max_functions));
    double** Va2p = Va_local->pointer();
    double** Vap = Va_AO->pointer();
    SharedMatrix Vb_local(new Matrix("Vb Temp", max_functions, max_functions));
    double** Vb2p = Vb_local->pointer();
    double** Vbp = Vb_AO->pointer();

    // Scratch
    std::vector<SharedMatrix> scratch = properties_->scratch();
    SharedMatrix Ta_local = scratch[0];
    SharedMatrix Tb_local = scratch[1]; 
    double** Tap = Ta_local->pointer();
    double** Tbp = Tb_local->pointer();

    // Traverse the blocks of points
    double functionalq = 0.0;
    double rhoaq       = 0.0;
    double rhoaxq      = 0.0;
    double rhoayq      = 0.0;
    double rhoazq      = 0.0;
    double rhobq       = 0.0;
    double rhobxq      = 0.0;
    double rhobyq      = 0.0;
    double rhobzq      = 0.0;
    boost::shared_ptr<Vector> QTa(new Vector("Quadrature Temp", max_points));
    double* QTap = QTa->pointer();
    boost::shared_ptr<Vector> QTb(new Vector("Quadrature Temp", max_points));
    double* QTbp = QTb->pointer();
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid_->blocks();
    for (int Q = 0; Q < blocks.size(); Q++) {

        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        timer_on("Properties");
        properties_->compute_points(block);
        timer_off("Properties");
        timer_on("Functional");
        std::map<std::string, SharedVector>& vals = functional_->compute_functional(properties_->point_values(), npoints); 
        timer_off("Functional");

        if (debug_ > 3) {
            block->print(outfile, debug_);
            properties_->print(outfile, debug_);
        }

        timer_on("V_XC");
        double** phi = properties_->basis_value("PHI")->pointer();
        double *restrict rho_a = properties_->point_value("RHO_A")->pointer();
        double *restrict rho_b = properties_->point_value("RHO_B")->pointer();
        double *restrict zk = vals["V"]->pointer(); 
        double *restrict v_rho_a = vals["V_RHO_A"]->pointer(); 
        double *restrict v_rho_b = vals["V_RHO_B"]->pointer(); 

        // => Quadrature values <= //
        functionalq += C_DDOT(npoints,w,1,zk,1);
        for (int P = 0; P < npoints; P++) {
            QTap[P] = w[P] * rho_a[P];
            QTbp[P] = w[P] * rho_b[P];
        }
        rhoaq       += C_DDOT(npoints,w,1,rho_a,1);
        rhoaxq      += C_DDOT(npoints,QTap,1,x,1);
        rhoayq      += C_DDOT(npoints,QTap,1,y,1);
        rhoazq      += C_DDOT(npoints,QTap,1,z,1);
        rhobq       += C_DDOT(npoints,w,1,rho_b,1);
        rhobxq      += C_DDOT(npoints,QTbp,1,x,1);
        rhobyq      += C_DDOT(npoints,QTbp,1,y,1);
        rhobzq      += C_DDOT(npoints,QTbp,1,z,1);

        // => LSDA contribution (symmetrized) <= //
        timer_on("LSDA");
        for (int P = 0; P < npoints; P++) {
            ::memset(static_cast<void*>(Tap[P]),'\0',nlocal*sizeof(double));
            ::memset(static_cast<void*>(Tbp[P]),'\0',nlocal*sizeof(double));
            C_DAXPY(nlocal,0.5 * v_rho_a[P] * w[P], phi[P], 1, Tap[P], 1); 
            C_DAXPY(nlocal,0.5 * v_rho_b[P] * w[P], phi[P], 1, Tbp[P], 1); 
        }
        timer_off("LSDA");
        
        // => GGA contribution (symmetrized) <= // 
        if (ansatz >= 1) {
            timer_on("GGA");
            double** phix = properties_->basis_value("PHI_X")->pointer();
            double** phiy = properties_->basis_value("PHI_Y")->pointer();
            double** phiz = properties_->basis_value("PHI_Z")->pointer();
            double *restrict rho_ax = properties_->point_value("RHO_AX")->pointer();
            double *restrict rho_ay = properties_->point_value("RHO_AY")->pointer();
            double *restrict rho_az = properties_->point_value("RHO_AZ")->pointer();
            double *restrict rho_bx = properties_->point_value("RHO_BX")->pointer();
            double *restrict rho_by = properties_->point_value("RHO_BY")->pointer();
            double *restrict rho_bz = properties_->point_value("RHO_BZ")->pointer();
            double *restrict v_sigma_aa = vals["V_GAMMA_AA"]->pointer(); 
            double *restrict v_sigma_ab = vals["V_GAMMA_AB"]->pointer(); 
            double *restrict v_sigma_bb = vals["V_GAMMA_BB"]->pointer(); 

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_ax[P] + v_sigma_ab[P] * rho_bx[P]), phix[P], 1, Tap[P], 1); 
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_ay[P] + v_sigma_ab[P] * rho_by[P]), phiy[P], 1, Tap[P], 1); 
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_aa[P] * rho_az[P] + v_sigma_ab[P] * rho_bz[P]), phiz[P], 1, Tap[P], 1); 
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_bb[P] * rho_bx[P] + v_sigma_ab[P] * rho_ax[P]), phix[P], 1, Tbp[P], 1); 
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_bb[P] * rho_by[P] + v_sigma_ab[P] * rho_ay[P]), phiy[P], 1, Tbp[P], 1); 
                C_DAXPY(nlocal,w[P] * (2.0 * v_sigma_bb[P] * rho_bz[P] + v_sigma_ab[P] * rho_az[P]), phiz[P], 1, Tbp[P], 1); 
            }        
            timer_off("GGA");
        }

        timer_on("LSDA");
        // Single GEMM slams GGA+LSDA together (man but GEM's hot!)
        C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phi[0],max_functions,Tap[0],max_functions,0.0,Va2p[0],max_functions);
        C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phi[0],max_functions,Tbp[0],max_functions,0.0,Vb2p[0],max_functions);

        // Symmetrization (V is Hermitian) 
        for (int m = 0; m < nlocal; m++) {
            for (int n = 0; n <= m; n++) {
                Va2p[m][n] = Va2p[n][m] = Va2p[m][n] + Va2p[n][m]; 
                Vb2p[m][n] = Vb2p[n][m] = Vb2p[m][n] + Vb2p[n][m]; 
            }
        }
        timer_off("LSDA");
        
        // => Meta contribution <= //
        if (ansatz >= 2) {
            timer_on("Meta");
            double** phix = properties_->basis_value("PHI_X")->pointer();
            double** phiy = properties_->basis_value("PHI_Y")->pointer();
            double** phiz = properties_->basis_value("PHI_Z")->pointer();
            double *restrict v_tau_a = vals["V_TAU_A"]->pointer(); 
            double *restrict v_tau_b = vals["V_TAU_B"]->pointer(); 
           
            // Alpha 
            // \nabla x
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tap[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_a[P] * w[P], phix[P], 1, Tap[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phix[0],max_functions,Tap[0],max_functions,1.0,Va2p[0],max_functions);

            // \nabla y
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tap[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_a[P] * w[P], phiy[P], 1, Tap[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phiy[0],max_functions,Tap[0],max_functions,1.0,Va2p[0],max_functions);

            // \nabla z
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tap[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_a[P] * w[P], phiz[P], 1, Tap[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phiz[0],max_functions,Tap[0],max_functions,1.0,Va2p[0],max_functions);
           
            // Beta 
            // \nabla x
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tbp[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_b[P] * w[P], phix[P], 1, Tbp[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phix[0],max_functions,Tbp[0],max_functions,1.0,Vb2p[0],max_functions);

            // \nabla y
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tbp[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_b[P] * w[P], phiy[P], 1, Tbp[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phiy[0],max_functions,Tbp[0],max_functions,1.0,Vb2p[0],max_functions);

            // \nabla z
            for (int P = 0; P < npoints; P++) {
                ::memset(static_cast<void*>(Tbp[P]),'\0',nlocal*sizeof(double));
                C_DAXPY(nlocal,v_tau_b[P] * w[P], phiz[P], 1, Tbp[P], 1); 
            }        
            C_DGEMM('T','N',nlocal,nlocal,npoints,1.0,phiz[0],max_functions,Tbp[0],max_functions,1.0,Vb2p[0],max_functions);
            timer_off("Meta");
        }       
 
        // => Unpacking <= //
        for (int ml = 0; ml < nlocal; ml++) {
            int mg = function_map[ml];
            for (int nl = 0; nl < ml; nl++) {
                int ng = function_map[nl];
                Vap[mg][ng] += Va2p[ml][nl];
                Vap[ng][mg] += Va2p[ml][nl];
                Vbp[mg][ng] += Vb2p[ml][nl];
                Vbp[ng][mg] += Vb2p[ml][nl];
            }
            Vap[mg][mg] += Va2p[ml][ml];
            Vbp[mg][mg] += Vb2p[ml][ml];
        }
        timer_off("V_XC");
    } 
   
    quad_values_["FUNCTIONAL"] = functionalq;
    quad_values_["RHO_A"]      = rhoaq; 
    quad_values_["RHO_AX"]     = rhoaxq; 
    quad_values_["RHO_AY"]     = rhoayq; 
    quad_values_["RHO_AZ"]     = rhoazq; 
    quad_values_["RHO_B"]      = rhobq; 
    quad_values_["RHO_BX"]     = rhobxq; 
    quad_values_["RHO_BY"]     = rhobyq; 
    quad_values_["RHO_BZ"]     = rhobzq; 
 
    if (debug_) {
        fprintf(outfile, "   => Numerical Integrals <=\n\n");
        fprintf(outfile, "    Functional Value:  %24.16E\n",quad_values_["FUNCTIONAL"]);
        fprintf(outfile, "    <\\rho_a>        :  %24.16E\n",quad_values_["RHO_A"]);
        fprintf(outfile, "    <\\rho_b>        :  %24.16E\n",quad_values_["RHO_B"]);
        fprintf(outfile, "    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",quad_values_["RHO_AX"],quad_values_["RHO_AY"],quad_values_["RHO_AZ"]);
        fprintf(outfile, "    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",quad_values_["RHO_BX"],quad_values_["RHO_BY"],quad_values_["RHO_BZ"]);
    }
}
SharedMatrix UV::compute_gradient()
{
    compute_D();
    USO2AO();

    if ((D_AO_.size() != 2) && (C_AO_.size() != 2))
        throw PSIEXCEPTION("V: UKS should have two D/C/V Matrices"); 

    // Build the target gradient Matrix
    int natom = primary_->molecule()->natom();
    SharedMatrix G(new Matrix("XC Gradient", natom,3));
    double** Gp = G->pointer();

    // Set Hessian derivative level in properties
    int old_deriv = properties_->deriv(); 
    properties_->set_deriv((functional_->is_gga() ? 2 : 1));

    // Setup the pointers
    SharedMatrix Da_AO = D_AO_[0];
    SharedMatrix Ca_AO = C_AO_[0];
    SharedMatrix Db_AO = D_AO_[1];
    SharedMatrix Cb_AO = C_AO_[1];
    properties_->set_pointers(Da_AO,Ca_AO, Db_AO, Cb_AO);

    // What local XC ansatz are we in?
    int ansatz = functional_->ansatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions(); 
    int max_points = grid_->max_points();

    // Scratch
    std::vector<SharedMatrix> scratch = properties_->scratch();
    SharedMatrix Ta_local = scratch[0];
    SharedMatrix Ua_local(Ta_local->clone());
    SharedMatrix Tb_local = scratch[1];
    SharedMatrix Ub_local(Tb_local->clone());
    double** Tap = Ta_local->pointer();
    double** Uap = Ua_local->pointer();
    double** Tbp = Tb_local->pointer();
    double** Ubp = Ub_local->pointer();
    std::vector<SharedMatrix> Dscratch = properties_->D_scratch();
    SharedMatrix Da_local = Dscratch[0];
    SharedMatrix Db_local = Dscratch[1];
    double** Dap = Da_local->pointer();
    double** Dbp = Db_local->pointer();

    // Traverse the blocks of points
    boost::shared_ptr<Vector> QT(new Vector("Quadrature Temp", max_points));
    double* QTp = QT->pointer();
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid_->blocks();

    for (std::map<std::string, double>::const_iterator it = quad_values_.begin(); it != quad_values_.end(); ++it) {
        quad_values_[(*it).first] = 0.0;
    }

    for (int Q = 0; Q < blocks.size(); Q++) {

        boost::shared_ptr<BlockOPoints> block = blocks[Q];
        int npoints = block->npoints();
        double* x = block->x();
        double* y = block->y();
        double* z = block->z();
        double* w = block->w();
        const std::vector<int>& function_map = block->functions_local_to_global();
        int nlocal = function_map.size();

        timer_on("Properties");
        properties_->compute_points(block);
        timer_off("Properties");
        timer_on("Functional");
        std::map<std::string, SharedVector>& vals = functional_->compute_functional(properties_->point_values(), npoints); 
        timer_off("Functional");

        double** phi = properties_->basis_value("PHI")->pointer();
        double** phi_x = properties_->basis_value("PHI_X")->pointer();
        double** phi_y = properties_->basis_value("PHI_Y")->pointer();
        double** phi_z = properties_->basis_value("PHI_Z")->pointer();
        double* rho_a = properties_->point_value("RHO_A")->pointer();
        double* rho_b = properties_->point_value("RHO_B")->pointer();
        double* zk = vals["V"]->pointer(); 
        double* v_rho_a = vals["V_RHO_A"]->pointer();
        double* v_rho_b = vals["V_RHO_B"]->pointer();

        // => Quadrature values <= //
        quad_values_["FUNCTIONAL"] += C_DDOT(npoints,w,1,zk,1); 
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_a[P];
        }
        quad_values_["RHO_A"] += C_DDOT(npoints,w,1,rho_a,1);
        quad_values_["RHO_AX"] += C_DDOT(npoints,QTp,1,x,1);
        quad_values_["RHO_AY"] += C_DDOT(npoints,QTp,1,y,1);
        quad_values_["RHO_AZ"] += C_DDOT(npoints,QTp,1,z,1);
        for (int P = 0; P < npoints; P++) {
            QTp[P] = w[P] * rho_b[P];
        }
        quad_values_["RHO_B"] += C_DDOT(npoints,w,1,rho_b,1);
        quad_values_["RHO_BX"] += C_DDOT(npoints,QTp,1,x,1);
        quad_values_["RHO_BY"] += C_DDOT(npoints,QTp,1,y,1);
        quad_values_["RHO_BZ"] += C_DDOT(npoints,QTp,1,z,1);
    
        // => LSDA Contribution <= //
        for (int P = 0; P < npoints; P++) {
            ::memset((void*) Tap[P], '\0', sizeof(double) * nlocal);
            ::memset((void*) Tbp[P], '\0', sizeof(double) * nlocal);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_a[P], phi[P], 1, Tap[P], 1);
            C_DAXPY(nlocal, -2.0 * w[P] * v_rho_b[P], phi[P], 1, Tbp[P], 1);
        }
    
        // => GGA Contribution (Term 1) <= //
        if (functional_->is_gga()) {
            double* rho_ax = properties_->point_value("RHO_AX")->pointer();
            double* rho_ay = properties_->point_value("RHO_AY")->pointer();
            double* rho_az = properties_->point_value("RHO_AZ")->pointer();
            double* rho_bx = properties_->point_value("RHO_BX")->pointer();
            double* rho_by = properties_->point_value("RHO_BY")->pointer();
            double* rho_bz = properties_->point_value("RHO_BZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();
            double* v_gamma_bb = vals["V_GAMMA_BB"]->pointer();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_bx[P]), phi_x[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_by[P]), phi_y[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_bz[P]), phi_z[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bx[P] + v_gamma_ab[P] * rho_ax[P]), phi_x[P], 1, Tbp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_by[P] + v_gamma_ab[P] * rho_ay[P]), phi_y[P], 1, Tbp[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bz[P] + v_gamma_ab[P] * rho_az[P]), phi_z[P], 1, Tbp[P], 1);
            }
    
        }

        // => Synthesis <= //
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,Tap[0],max_functions,Dap[0],max_functions,0.0,Uap[0],max_functions);
        C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,Tbp[0],max_functions,Dbp[0],max_functions,0.0,Ubp[0],max_functions);
    
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary_->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints,&Uap[0][ml],max_functions,&phi_x[0][ml],max_functions);
            Gp[A][1] += C_DDOT(npoints,&Uap[0][ml],max_functions,&phi_y[0][ml],max_functions);
            Gp[A][2] += C_DDOT(npoints,&Uap[0][ml],max_functions,&phi_z[0][ml],max_functions);
            Gp[A][0] += C_DDOT(npoints,&Ubp[0][ml],max_functions,&phi_x[0][ml],max_functions);
            Gp[A][1] += C_DDOT(npoints,&Ubp[0][ml],max_functions,&phi_y[0][ml],max_functions);
            Gp[A][2] += C_DDOT(npoints,&Ubp[0][ml],max_functions,&phi_z[0][ml],max_functions);
        }          
        
        // => GGA Contribution (Term 2) <= //
        if (functional_->is_gga()) {
            double** phi_xx = properties_->basis_value("PHI_XX")->pointer();
            double** phi_xy = properties_->basis_value("PHI_XY")->pointer();
            double** phi_xz = properties_->basis_value("PHI_XZ")->pointer();
            double** phi_yy = properties_->basis_value("PHI_YY")->pointer();
            double** phi_yz = properties_->basis_value("PHI_YZ")->pointer();
            double** phi_zz = properties_->basis_value("PHI_ZZ")->pointer();
            double* rho_ax = properties_->point_value("RHO_AX")->pointer();
            double* rho_ay = properties_->point_value("RHO_AY")->pointer();
            double* rho_az = properties_->point_value("RHO_AZ")->pointer();
            double* rho_bx = properties_->point_value("RHO_BX")->pointer();
            double* rho_by = properties_->point_value("RHO_BY")->pointer();
            double* rho_bz = properties_->point_value("RHO_BZ")->pointer();
            double* v_gamma_aa = vals["V_GAMMA_AA"]->pointer();
            double* v_gamma_ab = vals["V_GAMMA_AB"]->pointer();
            double* v_gamma_bb = vals["V_GAMMA_BB"]->pointer();

            C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi[0],max_functions,Dap[0],max_functions,0.0,Uap[0],max_functions);
            C_DGEMM('N','N',npoints,nlocal,nlocal,1.0,phi[0],max_functions,Dbp[0],max_functions,0.0,Ubp[0],max_functions);
            
            // x
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tap[P], '\0', sizeof(double) * nlocal);
                ::memset((void*) Tbp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P] + v_gamma_ab[P] * rho_bx[P]), Uap[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bx[P] + v_gamma_ab[P] * rho_ax[P]), Ubp[P], 1, Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_xx[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_xz[0][ml],max_functions);
                Gp[A][0] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_xx[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_xz[0][ml],max_functions);
            }          
            
            // y
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tap[P], '\0', sizeof(double) * nlocal);
                ::memset((void*) Tbp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P] + v_gamma_ab[P] * rho_by[P]), Uap[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_by[P] + v_gamma_ab[P] * rho_ay[P]), Ubp[P], 1, Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_yy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][0] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_xy[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_yy[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
            }          
            
            // z
            for (int P = 0; P < npoints; P++) {
                ::memset((void*) Tap[P], '\0', sizeof(double) * nlocal);
                ::memset((void*) Tbp[P], '\0', sizeof(double) * nlocal);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P] + v_gamma_ab[P] * rho_bz[P]), Uap[P], 1, Tap[P], 1);
                C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_bb[P] * rho_bz[P] + v_gamma_ab[P] * rho_az[P]), Ubp[P], 1, Tbp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary_->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_xz[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tap[0][ml],max_functions,&phi_zz[0][ml],max_functions);
                Gp[A][0] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_xz[0][ml],max_functions);
                Gp[A][1] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_yz[0][ml],max_functions);
                Gp[A][2] += C_DDOT(npoints,&Tbp[0][ml],max_functions,&phi_zz[0][ml],max_functions);
            }          
    
        }
        
        // => Meta Contribution <= //
        if (functional_->is_meta()) {
            throw PSIEXCEPTION("RV: Meta Gradients not implemented");
        }

    } 
 
    if (debug_) {
        fprintf(outfile, "   => XC Gradient: Numerical Integrals <=\n\n");
        fprintf(outfile, "    Functional Value:  %24.16E\n",quad_values_["FUNCTIONAL"]);
        fprintf(outfile, "    <\\rho_a>        :  %24.16E\n",quad_values_["RHO_A"]);
        fprintf(outfile, "    <\\rho_b>        :  %24.16E\n",quad_values_["RHO_B"]);
        fprintf(outfile, "    <\\vec r\\rho_a>  : <%24.16E,%24.16E,%24.16E>\n",quad_values_["RHO_AX"],quad_values_["RHO_AY"],quad_values_["RHO_AZ"]);
        fprintf(outfile, "    <\\vec r\\rho_b>  : <%24.16E,%24.16E,%24.16E>\n\n",quad_values_["RHO_BX"],quad_values_["RHO_BY"],quad_values_["RHO_BZ"]);
    }

    properties_->set_deriv(old_deriv);

    return G;
}

RK::RK(boost::shared_ptr<SuperFunctional> functional,
    boost::shared_ptr<BasisSet> primary,
    Options& options) : RV(functional,primary,options)
{
}
RK::~RK()
{
}
void RK::print_header() const
{
    VBase::print_header();
}
void RK::compute_V()
{
    // TODO
}

UK::UK(boost::shared_ptr<SuperFunctional> functional,
    boost::shared_ptr<BasisSet> primary,
    Options& options) : UV(functional,primary,options)
{
}
UK::~UK()
{
}
void UK::print_header() const
{
    VBase::print_header();
}
void UK::compute_V()
{
    // TODO
}

}
