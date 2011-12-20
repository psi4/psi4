#include <libfunctional/superfunctional.h>
#include <libmints/mints.h>
#include <libmints/points.h>
#include <libmints/cubature.h>
#include <libqt/qt.h>
#include <math.h>
#include "dft.h"
#include <psiconfig.h>

namespace psi {
namespace scf {

KSPotential::KSPotential(boost::shared_ptr<functional::SuperFunctional> functional,
    boost::shared_ptr<Molecule> molecule,
    boost::shared_ptr<BasisSet> primary,
    Options& options) :
    functional_(functional), molecule_(molecule),
    primary_(primary), options_(options)
{
    common_init();
}
KSPotential::~KSPotential()
{
}
void KSPotential::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");

    buildGrid();
    buildAO2USO();
}
void KSPotential::buildGrid()
{
    timer_on("Grid");
    grid_ = boost::shared_ptr<DFTGrid>(new DFTGrid(molecule_,primary_,options_));
    grid_->print(outfile,options_.get_int("PRINT"));
    timer_off("Grid");
}
void KSPotential::buildAO2USO()
{
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PetiteList> pet(new PetiteList(primary_, integral));
    if (pet->nirrep() != 1)
        AO2USO_ = SharedMatrix(pet->aotoso());
}
double KSPotential::quadrature_value(const std::string& key)
{
    return quad_values_[key];
}
void KSPotential::print(FILE* out, int print)
{
    fprintf(outfile, "  ==> KS Potential <==\n\n");
    molecule_->print(); 
    primary_->print_by_level(out,print);
    functional_->print(out, print);
    grid_->print(out, print);
}

RKSPotential::RKSPotential(boost::shared_ptr<functional::SuperFunctional> functional,
    boost::shared_ptr<Molecule> molecule,
    boost::shared_ptr<BasisSet> primary,
    Options& options) :
    KSPotential(functional,molecule,primary,options)
{
    buildProperties();
}
RKSPotential::~RKSPotential()
{
}
void RKSPotential::buildProperties()
{
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions();
    properties_ = boost::shared_ptr<RKSFunctions>(new RKSFunctions(primary_,max_points,max_functions));
    properties_->set_ansatz(functional_->getLocalAnsatz());
}
void RKSPotential::buildPotential(SharedMatrix D_USO, SharedMatrix C_USO, boost::shared_ptr<Dimension> noccpi)
{
    // Build D_AO_, C_AO_, allocate V_AO_ if need be
    timer_on("USO2AO");
    USO2AO(D_USO,C_USO,noccpi);
    timer_off("USO2AO");

    // Setup the pointers
    properties_->reset_pointers(D_AO_,C_AO_);

    // What local XC ansatz are we in?
    int ansatz = functional_->getLocalAnsatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions(); 
    int max_points = grid_->max_points();

    // Local/global V matrices
    SharedMatrix V_local(new Matrix("V Temp", max_functions, max_functions));
    double** V2p = V_local->pointer();
    double** Vp = V_AO_->pointer();

    // Scratch
    SharedMatrix T_local = properties_->scratch();
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
        properties_->computeProperties(block);
        timer_off("Properties");
        timer_on("Functional");
        functional_->computeRKSFunctional(properties_); 
        timer_off("Functional");

        if (debug_ > 4) {
            block->print(outfile, debug_);
            properties_->print(outfile, debug_);
        }

        timer_on("V_XC");
        double** phi = properties_->basis_value("PHI")->pointer();
        double *restrict rho_a = properties_->property_value("RHO_A")->pointer();
        double *restrict zk = functional_->getFunctionalValue();
        double *restrict v_rho_a = functional_->getV_RhoA();

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
            double *restrict rho_ax = properties_->property_value("RHO_AX")->pointer();
            double *restrict rho_ay = properties_->property_value("RHO_AY")->pointer();
            double *restrict rho_az = properties_->property_value("RHO_AZ")->pointer();
            double *restrict v_sigma_aa = functional_->getV_GammaAA();

            for (int P = 0; P < npoints; P++) {
                C_DAXPY(nlocal,v_sigma_aa[P] * rho_ax[P] * w[P], phix[P], 1, Tp[P], 1); 
                C_DAXPY(nlocal,v_sigma_aa[P] * rho_ay[P] * w[P], phiy[P], 1, Tp[P], 1); 
                C_DAXPY(nlocal,v_sigma_aa[P] * rho_az[P] * w[P], phiz[P], 1, Tp[P], 1); 
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
            double *restrict v_tau_a = functional_->getV_TauA();
            
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
    
    // V_AO_->V_USO_
    timer_on("AO2USO");
    AO2USO();
    timer_off("AO2USO");
}
void RKSPotential::USO2AO(SharedMatrix D_USO, SharedMatrix C_USO, boost::shared_ptr<Dimension> noccpi)
{
    int ansatz = functional_->getLocalAnsatz();
    int nirrep = D_USO->nirrep();
    int nao = primary_->nbf();
   
    // Allocate AO V matrix, if needed 
    if (V_AO_.get() == NULL) {
        V_AO_ = SharedMatrix(new Matrix("V (AO)", nao, nao)); 
    }
    V_AO_->zero();  
 
    // Move D_USO -> D_AO_ 
    if (nirrep == 1) {
        D_AO_ = D_USO;
    } else {
        if (D_AO_.get() == NULL) {
            D_AO_ = SharedMatrix(new Matrix("D(AO)", nao, nao));
        }

        D_AO_->zero();
        double** D_AOp = D_AO_->pointer();
        for (int h = 0; h < nirrep; h++) {
            int nso = AO2USO_->colspi()[h];
            if (nso == 0) continue;
            double** D_USOp = D_USO->pointer(h);
            double** Up = AO2USO_->pointer(h);
        
            SharedMatrix T(new Matrix("Temp", nao,nso));
            double** Tp = T->pointer();

            C_DGEMM('N','N',nao,nso,nso,1.0,Up[0],nso,D_USOp[0],nso,0.0,Tp[0],nso);
            C_DGEMM('N','T',nao,nao,nso,1.0,Tp[0],nso,Up[0],nso,1.0,D_AOp[0],nao);
        }
    }
    
    // Move C_USO -> C_AO_ if meta ansatz 
    if (ansatz >= 2) {
        int nocc = noccpi->sum();
        C_AO_ = SharedMatrix(new Matrix("Cocc (AO)", nao, nocc));
        double** C_AOp = C_AO_->pointer();

        if (nirrep == 1) {
            double** C_USOp = C_USO->pointer();

            for (int m = 0; m < nao; m++) {
                ::memcpy(static_cast<void*>(C_AOp), static_cast<void*>(C_USOp), nocc*sizeof(double));
            }
        
        } else {
            int offset = 0;
            for (int h = 0; h < nirrep; h++) {
                int nso = AO2USO_->colspi()[h];
                int nmo = C_USO->colspi()[h];
                int na  = (*noccpi)[h];
    
                if (nso == 0 || nmo == 0 || na == 0) continue;

                double** C_USOp = C_USO->pointer(h);
                double** Up = AO2USO_->pointer(h);
        
                C_DGEMM('N','N',nao,na,nso,1.0,Up[0],nso,C_USOp[0],nmo,0.0,&C_AOp[0][offset],nocc);

                offset += na; 
            } 
        } 
    }
}
void RKSPotential::AO2USO()
{
    if (AO2USO_.get() == NULL) {
        V_USO_ = V_AO_;
        return;
    }

    int nirrep = AO2USO_->nirrep();
    int nao = V_AO_->colspi()[0];
    int* dimpi = AO2USO_->colspi();
    
    // Allocate USO V matrix, if needed
    if (V_USO_.get() == NULL) {
        V_USO_ = SharedMatrix(new Matrix("V (USO)", nirrep, dimpi, dimpi)); 
    }

    // Move V_AO_->V_USO_
    double** V_AOp = V_AO_->pointer();
    for (int h = 0; h < nirrep; h++) {
        int nso = dimpi[h];
        if (nso == 0) continue;
        double** V_USOp = V_USO_->pointer(h);
        double** Up = AO2USO_->pointer(h);
        
        SharedMatrix T(new Matrix("Temp", nao,nso));
        double** Tp = T->pointer();

        C_DGEMM('N','N',nao,nso,nao,1.0,V_AOp[0],nao,Up[0],nso,0.0,Tp[0],nso); 
        C_DGEMM('T','N',nso,nso,nao,1.0,Up[0],nso,Tp[0],nso,0.0,V_USOp[0],nso); 
        
    } 
}
void RKSPotential::print(FILE* out, int print)
{
    fprintf(outfile, "  ==> RKS Potential <==\n\n");
    properties_->print(out,print);

    KSPotential::print(out,print);
}

UKSPotential::UKSPotential(boost::shared_ptr<functional::SuperFunctional> functional,
    boost::shared_ptr<Molecule> molecule,
    boost::shared_ptr<BasisSet> primary,
    Options& options) :
    KSPotential(functional,molecule,primary,options)
{
    buildProperties();
}
UKSPotential::~UKSPotential()
{
}
void UKSPotential::buildProperties()
{
    int max_points = grid_->max_points();
    int max_functions = grid_->max_functions(); 
    properties_ = boost::shared_ptr<UKSFunctions>(new UKSFunctions(primary_,max_points,max_functions));
    properties_->set_ansatz(functional_->getLocalAnsatz());
}
void UKSPotential::buildPotential(SharedMatrix Da_USO, SharedMatrix Ca_USO, boost::shared_ptr<Dimension> napi,
                                  SharedMatrix Db_USO, SharedMatrix Cb_USO, boost::shared_ptr<Dimension> nbpi)
{
    // Build D_AO_, C_AO_, allocate V_AO_ if need be
    timer_on("USO2AO");
    USO2AO(Da_USO,Ca_USO,napi,Db_USO,Cb_USO,nbpi);
    timer_off("USO2AO");

    // Setup the pointers
    properties_->reset_pointers(Da_AO_,Ca_AO_,Db_AO_,Cb_AO_);

    // What local XC ansatz are we in?
    int ansatz = functional_->getLocalAnsatz();

    // How many functions are there (for lda in Vtemp, T)
    int max_functions = grid_->max_functions();
    int max_points = grid_->max_points();

    // Local/global V matrices
    SharedMatrix Va_local(new Matrix("Va Temp", max_functions, max_functions));
    double** Va2p = Va_local->pointer();
    double** Vap = Va_AO_->pointer();
    SharedMatrix Vb_local(new Matrix("Vb Temp", max_functions, max_functions));
    double** Vb2p = Vb_local->pointer();
    double** Vbp = Vb_AO_->pointer();

    // Scratch
    SharedMatrix Ta_local = properties_->scratchA();
    double** Tap = Ta_local->pointer();
    SharedMatrix Tb_local = properties_->scratchB();
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
        properties_->computeProperties(block);
        timer_off("Properties");
        timer_on("Functional");
        functional_->computeUKSFunctional(properties_); 
        timer_off("Functional");

        if (debug_ > 3) {
            block->print(outfile, debug_);
            properties_->print(outfile, debug_);
        }

        timer_on("V_XC");
        double** phi = properties_->basis_value("PHI")->pointer();
        double *restrict rho_a = properties_->property_value("RHO_A")->pointer();
        double *restrict rho_b = properties_->property_value("RHO_B")->pointer();
        double *restrict zk = functional_->getFunctionalValue();
        double *restrict v_rho_a = functional_->getV_RhoA();
        double *restrict v_rho_b = functional_->getV_RhoB();

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
            double *restrict rho_ax = properties_->property_value("RHO_AX")->pointer();
            double *restrict rho_ay = properties_->property_value("RHO_AY")->pointer();
            double *restrict rho_az = properties_->property_value("RHO_AZ")->pointer();
            double *restrict rho_bx = properties_->property_value("RHO_BX")->pointer();
            double *restrict rho_by = properties_->property_value("RHO_BY")->pointer();
            double *restrict rho_bz = properties_->property_value("RHO_BZ")->pointer();
            double *restrict v_sigma_aa = functional_->getV_GammaAA();
            double *restrict v_sigma_ab = functional_->getV_GammaAB();
            double *restrict v_sigma_bb = functional_->getV_GammaBB();

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
            double *restrict v_tau_a = functional_->getV_TauA();
            double *restrict v_tau_b = functional_->getV_TauB();
           
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
    // V_AO_->V_USO_
    timer_on("AO2USO");
    AO2USO();
    timer_off("AO2USO");
}
void UKSPotential::USO2AO(SharedMatrix Da_USO, SharedMatrix Ca_USO, boost::shared_ptr<Dimension> napi,
                          SharedMatrix Db_USO, SharedMatrix Cb_USO, boost::shared_ptr<Dimension> nbpi)
{
    int ansatz = functional_->getLocalAnsatz();
    int nirrep = Da_USO->nirrep();
    int nao = primary_->nbf();
   
    // Allocate AO V matrix, if needed 
    if (Va_AO_.get() == NULL) {
        Va_AO_ = SharedMatrix(new Matrix("Va (AO)", nao, nao)); 
        Vb_AO_ = SharedMatrix(new Matrix("Vb (AO)", nao, nao)); 
    }
    Va_AO_->zero();  
    Vb_AO_->zero();  
 
    // Move D_USO -> D_AO_ 
    if (nirrep == 1) {
        Da_AO_ = Da_USO;
        Db_AO_ = Db_USO;
    } else {
        if (Da_AO_.get() == NULL) {
            Da_AO_ = SharedMatrix(new Matrix("Da (AO)", nao, nao));
            Db_AO_ = SharedMatrix(new Matrix("Db (AO)", nao, nao));
        }

        Da_AO_->zero();
        double** Da_AOp = Da_AO_->pointer();
        for (int h = 0; h < nirrep; h++) {
            int nso = AO2USO_->colspi()[h];
            if (nso == 0) continue;
            double** Da_USOp = Da_USO->pointer(h);
            double** Up = AO2USO_->pointer(h);
        
            SharedMatrix T(new Matrix("Temp", nao,nso));
            double** Tp = T->pointer();

            C_DGEMM('N','N',nao,nso,nso,1.0,Up[0],nso,Da_USOp[0],nso,0.0,Tp[0],nso);
            C_DGEMM('N','T',nao,nao,nso,1.0,Tp[0],nso,Up[0],nso,1.0,Da_AOp[0],nao);
        }

        Db_AO_->zero();
        double** Db_AOp = Db_AO_->pointer();
        for (int h = 0; h < nirrep; h++) {
            int nso = AO2USO_->colspi()[h];
            if (nso == 0) continue;
            double** Db_USOp = Db_USO->pointer(h);
            double** Up = AO2USO_->pointer(h);
        
            SharedMatrix T(new Matrix("Temp", nao,nso));
            double** Tp = T->pointer();

            C_DGEMM('N','N',nao,nso,nso,1.0,Up[0],nso,Db_USOp[0],nso,0.0,Tp[0],nso);
            C_DGEMM('N','T',nao,nao,nso,1.0,Tp[0],nso,Up[0],nso,1.0,Db_AOp[0],nao);
        }
    }
    
    // Move C_USO -> C_AO_ if meta ansatz 
    if (ansatz >= 2) {
        int nalpha = napi->sum();
        int nbeta  = nbpi->sum();
        Ca_AO_ = SharedMatrix(new Matrix("Caocc (AO)", nao, nalpha));
        Cb_AO_ = SharedMatrix(new Matrix("Cbocc (AO)", nao, nbeta));
        double** Ca_AOp = Ca_AO_->pointer();
        double** Cb_AOp = Cb_AO_->pointer();

        if (nirrep == 1) {
            double** Ca_USOp = Ca_USO->pointer();
            double** Cb_USOp = Cb_USO->pointer();

            for (int m = 0; m < nao; m++) {
                ::memcpy(static_cast<void*>(Ca_AOp), static_cast<void*>(Ca_USOp), nalpha*sizeof(double));
                ::memcpy(static_cast<void*>(Cb_AOp), static_cast<void*>(Cb_USOp), nbeta *sizeof(double));
            }
        
        } else {
            int offset = 0;
            for (int h = 0; h < nirrep; h++) {
                int nso = AO2USO_->colspi()[h];
                int nmo = Ca_USO->colspi()[h];
                int na  = (*napi)[h];
    
                if (nso == 0 || nmo == 0 || na == 0) continue;

                double** Ca_USOp = Ca_USO->pointer(h);
                double** Up = AO2USO_->pointer(h);
        
                C_DGEMM('N','N',nao,na,nso,1.0,Up[0],nso,Ca_USOp[0],nmo,0.0,&Ca_AOp[0][offset],nalpha);

                offset += na; 
            } 

            offset = 0;
            for (int h = 0; h < nirrep; h++) {
                int nso = AO2USO_->colspi()[h];
                int nmo = Ca_USO->colspi()[h];
                int nb  = (*nbpi)[h];
    
                if (nso == 0 || nmo == 0 || nb == 0) continue;

                double** Cb_USOp = Cb_USO->pointer(h);
                double** Up = AO2USO_->pointer(h);
        
                C_DGEMM('N','N',nao,nb,nso,1.0,Up[0],nso,Cb_USOp[0],nmo,0.0,&Cb_AOp[0][offset],nbeta);

                offset += nb; 
            } 
        } 
    }
}
void UKSPotential::AO2USO()
{
    if (AO2USO_.get() == NULL) {
        Va_USO_ = Va_AO_;
        Vb_USO_ = Vb_AO_;
        return;
    }

    int nirrep = AO2USO_->nirrep();
    int nao = Va_AO_->colspi()[0];
    int* dimpi = AO2USO_->colspi();
    
    // Allocate USO V matrix, if needed
    if (Va_USO_.get() == NULL) {
        Va_USO_ = SharedMatrix(new Matrix("Va (USO)", nirrep, dimpi, dimpi)); 
        Vb_USO_ = SharedMatrix(new Matrix("Vb (USO)", nirrep, dimpi, dimpi)); 
    }

    // Move V_AO_->V_USO_
    double** Va_AOp = Va_AO_->pointer();
    for (int h = 0; h < nirrep; h++) {
        int nso = dimpi[h];
        if (nso == 0) continue;
        double** Va_USOp = Va_USO_->pointer(h);
        double** Up = AO2USO_->pointer(h);
        
        SharedMatrix T(new Matrix("Temp", nao,nso));
        double** Tp = T->pointer();

        C_DGEMM('N','N',nao,nso,nao,1.0,Va_AOp[0],nao,Up[0],nso,0.0,Tp[0],nso); 
        C_DGEMM('T','N',nso,nso,nao,1.0,Up[0],nso,Tp[0],nso,0.0,Va_USOp[0],nso); 
        
    } 

    // Move V_AO_->V_USO_
    double** Vb_AOp = Vb_AO_->pointer();
    for (int h = 0; h < nirrep; h++) {
        int nso = dimpi[h];
        if (nso == 0) continue;
        double** Vb_USOp = Vb_USO_->pointer(h);
        double** Up = AO2USO_->pointer(h);
        
        SharedMatrix T(new Matrix("Temp", nao,nso));
        double** Tp = T->pointer();

        C_DGEMM('N','N',nao,nso,nao,1.0,Vb_AOp[0],nao,Up[0],nso,0.0,Tp[0],nso); 
        C_DGEMM('T','N',nso,nso,nao,1.0,Up[0],nso,Tp[0],nso,0.0,Vb_USOp[0],nso); 
        
    } 
}
void UKSPotential::print(FILE* out, int print)
{
    fprintf(outfile, "  ==> UKS Potential <==\n\n");
    properties_->print(out,print);

    KSPotential::print(out,print);
}


}} // Namespace psi::scf
