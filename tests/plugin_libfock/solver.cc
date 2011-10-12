#include <libmints/mints.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include "solver.h"
#include "hamiltonian.h"
#include "jk.h"

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {

Solver::Solver()
{
    common_init();
}
Solver::~Solver()
{
}
void Solver::common_init()
{
    print_ = 1;
    debug_ = 0;
    // Unlimited default
    memory_ = 0L; 
    converged_ = false;
    iteration_ = 0;
    convergence_ = 0.0;
    criteria_ = 1.0E-6;
    maxiter_ = 100;
}

RSolver::RSolver(boost::shared_ptr<RHamiltonian> H) :
    Solver(), H_(H)
{
}
RSolver::~RSolver()
{
}

USolver::USolver(boost::shared_ptr<UHamiltonian> H) :
    Solver(), H_(H)
{
}
USolver::~USolver()
{
}

CGRSolver::CGRSolver(boost::shared_ptr<RHamiltonian> H,
          std::vector<boost::shared_ptr<Vector> >& b):
    RSolver(H), b_(b), precondition_(true)
{
}
CGRSolver::~CGRSolver()
{
}
boost::shared_ptr<CGRSolver> CGRSolver::build_solver(Options& options,
    boost::shared_ptr<RHamiltonian> H)
{
    std::vector<boost::shared_ptr<Vector> > b;
    boost::shared_ptr<CGRSolver> solver(new CGRSolver(H,b));

    if (options["PRINT"].has_changed()) {
        solver->set_print(options.get_int("PRINT"));
    }
    if (options["DEBUG"].has_changed()) {
        solver->set_debug(options.get_int("DEBUG"));
    }
    if (options["SOLVER_MAXITER"].has_changed()) {
        solver->set_maxiter(options.get_int("SOLVER_MAXITER"));
    }
    if (options["SOLVER_CONVERGENCE"].has_changed()) {
        solver->set_convergence(options.get_double("SOLVER_CONVERGENCE"));
    }

    return solver;
}
void CGRSolver::print_header() const
{
    if (print_) {
    	fprintf(outfile, "  ==> CGRSolver (by Rob Parrish) <==\n");
        fprintf(outfile, "   Number of roots    = %d\n", b_.size());
        fprintf(outfile, "   Preconditioning    = %s\n\n", (precondition_ ? "Yes" : "No"));
        fprintf(outfile, "   Convergence cutoff = %9.0E\n", criteria_);
        fprintf(outfile, "   Maximum iterations = %9d\n\n", maxiter_); 
    }
}
unsigned long int CGRSolver::memory_estimate()
{
    unsigned long int dimension = 0L;
    if (!diag_) diag_ = H_->diagonal();
    for (int h = 0; h < diag_->nirrep(); h++) {
        dimension += diag_->dimpi()[h];
    }
    return (6L * b_.size()) * dimension;
}
void CGRSolver::initialize()
{
    finalize();

    int nvec = b_.size();
    for (int N = 0; N < nvec; ++N) {
        std::stringstream xs;
        xs << "Solution Vector " << N+1;
        x_.push_back(boost::shared_ptr<Vector>(new Vector(xs.str(),b_[0]->nirrep(),b_[0]->dimpi())));
        std::stringstream Aps;
        Aps << "Product Vector " << N+1;
        Ap_.push_back(boost::shared_ptr<Vector>(new Vector(Aps.str(),b_[0]->nirrep(),b_[0]->dimpi())));
        std::stringstream zs;
        zs << "Z Vector " << N+1;
        z_.push_back(boost::shared_ptr<Vector>(new Vector(zs.str(),b_[0]->nirrep(),b_[0]->dimpi())));
        std::stringstream rs;
        rs << "Residual Vector " << N+1;
        r_.push_back(boost::shared_ptr<Vector>(new Vector(rs.str(),b_[0]->nirrep(),b_[0]->dimpi())));
        std::stringstream ps;
        ps << "Conjugate Vector " << N+1;
        p_.push_back(boost::shared_ptr<Vector>(new Vector(ps.str(),b_[0]->nirrep(),b_[0]->dimpi())));
        alpha_.push_back(0.0);
        beta_.push_back(0.0);
        r_nrm2_.push_back(0.0);
        z_r_.push_back(0.0);
        r_converged_.push_back(false);
    }

    diag_ = H_->diagonal();
}
void CGRSolver::solve()
{
    iteration_ = 0;
    converged_ = false;
    nconverged_ = 0;
    convergence_ = 0.0;

    if (print_) {
        fprintf(outfile, "  => Iterations <=\n\n");
        fprintf(outfile, "  %4s %11s %11s %11s\n", "Iter", "Max |R|_2", "Converged", "Remaining");
        fflush(outfile);
    }

    guess();
    products_x();
    residual();
    update_z();
    update_p();

    do {
        iteration_++;

        products_p();
        alpha();
        update_x();
        update_r();
        check_convergence();
        if (print_) {
            fprintf(outfile, "  %4d %11.3E %11d %11d\n", iteration_, convergence_, nconverged_, b_.size() - nconverged_);
            fflush(outfile);
        }
        update_z();
        beta();
        update_p();

    } while (iteration_ < maxiter_ && !converged_);

    if (print_) {
        fprintf(outfile, "\n");
    }

    if (!converged_) {
        throw PSIEXCEPTION("CGRSolver did not converge");
    } else if (print_) {
        fprintf(outfile, "    CGRSolver converged.\n\n");
    }

    fflush(outfile);
}
void CGRSolver::finalize()
{
    Ap_.clear();
    z_.clear();
    r_.clear();
    p_.clear();
    alpha_.clear();
    beta_.clear();
    r_nrm2_.clear();
    z_r_.clear();
    r_converged_.clear();
    diag_.reset();
}
void CGRSolver::guess()
{
    for (int N = 0; N < b_.size(); ++N) {
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* bp = b_[N]->pointer();
            double* xp = x_[N]->pointer();
            double* dp = diag_->pointer();
            for (int i = 0; i < n; ++i) {
                xp[i] = bp[i] / dp[i];
            }
        }
    }
}
void CGRSolver::residual()
{
    for (int N = 0; N < b_.size(); ++N) {
        x_[N]->copy(b_[N].get());
        x_[N]->scale(-1.0);
        x_[N]->add(Ap_[N]);
    }
}
void CGRSolver::products_x()
{
    H_->product(x_,Ap_);
}
void CGRSolver::products_p()
{
    std::vector<boost::shared_ptr<Vector> > p;
    std::vector<boost::shared_ptr<Vector> > Ap;

    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        p.push_back(p_[N]);
        Ap.push_back(Ap_[N]);
    }

    H_->product(p,Ap);
}
void CGRSolver::alpha()
{
    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        z_r_[N] = 0.0;
        double p_Ap = 0.0;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer();
            double* zp = z_[N]->pointer();
            double* pp = p_[N]->pointer();
            double* App = Ap_[N]->pointer();
            z_r_[N] += C_DDOT(n,rp,1,zp,1);
            p_Ap += C_DDOT(n,pp,1,App,1);
        }
        alpha_[N] = z_r_[N] / p_Ap;
    }
}
void CGRSolver::update_x()
{
    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* xp = x_[N]->pointer();
            double* pp = p_[N]->pointer();
            C_DAXPY(n,alpha_[N],pp,1,xp,1);
        }
    }
}
void CGRSolver::update_r()
{
    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer();
            double* App = Ap_[N]->pointer();
            C_DAXPY(n,-alpha_[N],App,1,rp,1);
        }
    }
}
void CGRSolver::check_convergence()
{
    convergence_ = 0.0;
    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        r_nrm2_[N] = 0.0;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer();
            r_nrm2_[N] += C_DDOT(n,rp,1,rp,1);
        }
        sqrt(r_nrm2_[N]);
        if (convergence_ < r_nrm2_[N]) {
            convergence_ = r_nrm2_[N];
        }
        if (r_nrm2_[N] < criteria_) {
            r_converged_[N] = true;
            nconverged_++;
        }
    }
    if (nconverged_ = b_.size()) {
        converged_ = true;
    }
}
void CGRSolver::update_z()
{
    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* zp = z_[N]->pointer();
            double* xp = x_[N]->pointer();
            double* dp = diag_->pointer();
            if (precondition_) {
                for (int i = 0; i < n; ++i) {
                    zp[i] = xp[i] / dp[i];
                }
            } else { 
                for (int i = 0; i < n; ++i) {
                    zp[i] = xp[i];
                }
            }
        }
    }
}
void CGRSolver::beta()
{
    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        double zr = 0.0;
        for (int h = 0; h < b_[N]->nirrep(); ++h) {
            int n = b_[N]->dimpi()[h];
            if (!n) continue;
            double* rp = r_[N]->pointer();
            double* zp = z_[N]->pointer();
            zr += C_DDOT(n,rp,1,zp,1);
        }
        beta_[N] = zr / z_r_[N];
    }
}
void CGRSolver::update_p()
{
    for (int N = 0; N < b_.size(); ++N) {
        if (r_converged_[N]) continue;
        p_[N]->scale(beta_[N]);
        p_[N]->add(z_[N]);
    }
}

DLRSolver::DLRSolver(boost::shared_ptr<RHamiltonian> H) :
    RSolver(H),
    nroot_(1),
    nguess_(1),
    norm_(1.E-6),
    max_subspace_(6),
    min_subspace_(2),
    nsubspace_(0),
    nconverged_(0)
{
}
DLRSolver::~DLRSolver()
{
}
boost::shared_ptr<DLRSolver> DLRSolver::build_solver(Options& options,
    boost::shared_ptr<RHamiltonian> H)
{
    boost::shared_ptr<DLRSolver> solver(new DLRSolver(H));

    if (options["PRINT"].has_changed()) {
        solver->set_print(options.get_int("PRINT"));
    }
    if (options["DEBUG"].has_changed()) {
        solver->set_debug(options.get_int("DEBUG"));
    }
    if (options["SOLVER_MAXITER"].has_changed()) {
        solver->set_maxiter(options.get_int("SOLVER_MAXITER"));
    }
    if (options["SOLVER_CONVERGENCE"].has_changed()) {
        solver->set_convergence(options.get_double("SOLVER_CONVERGENCE"));
    }
    if (options["SOLVER_N_ROOT"].has_changed()) {
        solver->set_nroot(options.get_int("SOLVER_N_ROOT"));
    }
    if (options["SOLVER_N_GUESS"].has_changed()) {
        solver->set_nguess(options.get_int("SOLVER_N_GUESS"));
    }
    if (options["SOLVER_MIN_SUBSPACE"].has_changed()) {
        solver->set_min_subspace(options.get_int("SOLVER_MIN_SUBSPACE"));
    }
    if (options["SOLVER_MAX_SUBSPACE"].has_changed()) {
        solver->set_max_subspace(options.get_int("SOLVER_MAX_SUBSPACE"));
    }
    if (options["SOLVER_NORM"].has_changed()) {
        solver->set_norm(options.get_double("SOLVER_NORM"));
    }

    return solver;
}
void DLRSolver::print_header() const 
{
    if (print_) {
        fprintf(outfile, "  ==> DLRSolver (by Rob Parrish) <== \n\n");
        fprintf(outfile, "   Number of roots         = %9d\n", nroot_);
        fprintf(outfile, "   Number of guess vectors = %9d\n", nguess_);
        fprintf(outfile, "   Maximum subspace size   = %9d\n", max_subspace_);
        fprintf(outfile, "   Minimum subspace size   = %9d\n", min_subspace_);
        fprintf(outfile, "   Subspace expansion norm = %9.0E\n", norm_);
        fprintf(outfile, "   Convergence cutoff      = %9.0E\n", criteria_);
        fprintf(outfile, "   Maximum iterations      = %9d\n\n", maxiter_); 
    }
}
unsigned long int DLRSolver::memory_estimate()
{
    unsigned long int dimension = 0L;
    if (!diag_) diag_ = H_->diagonal();
    for (int h = 0; h < diag_->nirrep(); h++) {
        dimension += diag_->dimpi()[h];
    }
    return (2L * max_subspace_ + 3L * nroot_ + 1L) * dimension;
}
void DLRSolver::initialize()
{
    finalize();

    c_.clear();
    E_.clear();

    diag_ = H_->diagonal();
} 
void DLRSolver::solve()
{
    iteration_ = 0;
    converged_ = false;
    nconverged_ = 0;
    convergence_ = 0.0;

    if (print_) {
        fprintf(outfile, "  => Iterations <=\n\n"); 
        fprintf(outfile, "   %4s  %10s  %10s  %11s\n", "Iter", "NConverged", "NSubspace", "Residual"); 
        fflush(outfile);
    }

    // Compute the first set of sigma vectors
    guess();
    sigma();
    
    do {
        iteration_++;
        
        // Compute subspace Hamiltonian
        subspaceHamiltonian();
        // Diagonalize subspace Hamiltonian
        subspaceDiagonalization();
        // Find eigenvectors 
        eigenvecs();
        // Find eigenvalues
        eigenvals();   
        // Find residuals, update convergence 
        residuals();

        fprintf(outfile, "   %4d  %10d  %10d  %11.3E\n", iteration_, nconverged_, 
            nsubspace_, convergence_);
        fflush(outfile);       
 
        // Check for convergence
        if (converged_) break;

        // Find correctors 
        correctors();
        // Collapse subspace if needed
        subspaceCollapse();
        // Orthogonalize/add significant correctors 
        subspaceExpansion();
        // Compute new sigma vectors 
        sigma();

    } while (iteration_ < maxiter_ && !converged_);

    if (print_) {
        fprintf(outfile, "\n");
    }

    if (!converged_) {
        throw PSIEXCEPTION("DLRSolver did not converge");
    } else if (print_) {
        fprintf(outfile, "    DLRSolver converged.\n\n");
    }

    fflush(outfile);
}
void DLRSolver::finalize()
{
    b_.clear();
    s_.clear();
    G_.reset();
    a_.reset();
    l_.reset();
    r_.clear();
    n_.clear();
    d_.clear();
    diag_.reset();
} 
void DLRSolver::guess()
{
    for (int i = 0; i < nguess_; ++i) {
        std::stringstream ss;
        ss << "Subspace Vector " << i;
        b_.push_back(boost::shared_ptr<Vector>(new Vector(ss.str(), diag_->nirrep(), diag_->dimpi())));
    } 

    for (int h = 0; h < diag_->nirrep(); ++h) {
        int n = diag_->dimpi()[h];
        if (!n) continue;

        std::vector<std::pair<double, int> > d;
        for (int i = 0; i < n; ++i) {
            d.push_back(make_pair(diag_->get(h,i),i));
        }
        std::sort(d.begin(), d.end());

        for (int i = 0; (i < nguess_) && (i < n); ++i) {
            b_[i]->set(h,d[i].second,1.0);    
        } 
    }

    nsubspace_ = nguess_;

    if (debug_) {
        fprintf(outfile, "   > Guess <\n\n");
        diag_->print();
        for (int i = 0; i < b_.size(); ++i) {
            b_[i]->print();
        }
        fflush(outfile);
    }
}
void DLRSolver::sigma()
{
    int n = b_.size() - s_.size();
    int offset = s_.size();
    for (int i = 0; i < n; i++) {   
        std::stringstream s;
        s << "Sigma Vector " << (i + offset);
        s_.push_back(boost::shared_ptr<Vector>(new Vector(s.str(), diag_->nirrep(), diag_->dimpi())));
    }
    
    std::vector<boost::shared_ptr<Vector> > x;
    std::vector<boost::shared_ptr<Vector> > b;

    for (int i = offset; i < offset + n; i++) {
        x.push_back(b_[i]);
        b.push_back(s_[i]);
    }

    H_->product(x,b);

    if (debug_) { 
        fprintf(outfile, "   > Sigma <\n\n");
        for (int i = 0; i < s_.size(); i++) {
            s_[i]->print();
        }
        fflush(outfile);
    }
}
void DLRSolver::subspaceHamiltonian()
{
    int n = s_.size();
    int nirrep = diag_->nirrep();
    int* npi = new int[nirrep];
    for (int h = 0; h < nirrep; ++h) {
        npi[h] = n;
    }

    G_ = boost::shared_ptr<Matrix> (new Matrix("Subspace Hamiltonian",nirrep,npi,npi));
    delete[] npi;

    for (int h = 0; h < nirrep; ++h) {

        int dimension = diag_->dimpi()[h];

        if (!dimension) continue;

        double** Gp = G_->pointer(h);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                Gp[i][j] = Gp[j][i] = C_DDOT(dimension,b_[i]->pointer(h),1,s_[j]->pointer(h),1);
            }
        }
    }

    if (debug_) { 
        fprintf(outfile, "   > SubspaceHamiltonian <\n\n");
        G_->print();
        fflush(outfile);
    }
}
void DLRSolver::subspaceDiagonalization()
{
    int n = s_.size();
    int nirrep = diag_->nirrep();
    int* npi = new int[nirrep];
    for (int h = 0; h < nirrep; ++h) {
        npi[h] = n;
    }

    boost::shared_ptr<Matrix> G2(G_->clone());
    a_ = boost::shared_ptr<Matrix> (new Matrix("Subspace Eigenvectors",nirrep,npi,npi));
    l_ = boost::shared_ptr<Vector> (new Vector("Subspace Eigenvalues",nirrep,npi));
    delete[] npi;
   
    G2->diagonalize(a_,l_); 

    if (debug_) { 
        fprintf(outfile, "   > SubspaceDiagonalize <\n\n");
        a_->print();
        l_->print();
        fflush(outfile);
    }
}
void DLRSolver::eigenvecs()
{
    if (c_.size() != nroot_) {
        c_.clear();
        for (int m = 0; m < nroot_; ++m) {
            std::stringstream s;
            s << "Eigenvector " << m;
            boost::shared_ptr<Vector> c(new Vector(s.str().c_str(),diag_->nirrep(), diag_->dimpi()));
            c_.push_back(c);
        }
    }

    for (int h = 0; h < diag_->nirrep(); ++h) {

        int dimension = diag_->dimpi()[h];
        
        if (!dimension) continue;

        double** ap = a_->pointer(h);
        for (int m = 0; m < nroot_; m++) {
            double* cp = c_[m]->pointer(h);
            ::memset((void*) cp, '\0', dimension*sizeof(double));
            for (int i = 0; i < b_.size(); i++) {
                double* bp = b_[i]->pointer(h);
                C_DAXPY(dimension,ap[i][m],bp,1,cp,1);
            }       
        } 
    }

    if (debug_) { 
        fprintf(outfile, "   > Eigenvectors <\n\n");
        for (int m = 0; m < c_.size(); m++) {
            c_[m]->print();
        }
        fflush(outfile);
    }
}
void DLRSolver::eigenvals()
{
    E_.clear();
    E_.resize(nroot_);

    for (int h = 0; h < diag_->nirrep(); ++h) {

        double** ap = a_->pointer(h);
        double** Gp = G_->pointer(h);

        for (int k = 0; k < nroot_; k++) {
            double E = 0.0;
            for (int i = 0; i < b_.size(); i++) {
                for (int j = 0; j < nsubspace_; j++) {
                    E += ap[i][k] * ap[j][k] * Gp[i][j];    
                }
            }
            E_[k].push_back(E);
        }
    }

    if (debug_) { 
        fprintf(outfile, "   > Eigenvalues <\n\n");
        for (int m = 0; m < E_.size(); m++) {
            for (int h = 0; h < E_[0].size(); ++h) {
                fprintf(outfile, "    Eigenvalue %d, Irrep %d = %24.16E\n", m, h, E_[m][h]);
            }
        }
        fprintf(outfile, "\n");
        fflush(outfile);
    }
}
void DLRSolver::residuals()
{
    n_.resize(nroot_);
    nconverged_ = 0;

    if (r_.size() != nroot_) {
        r_.clear();
        for (int k = 0; k < nroot_; ++k) {
            // Residual k
            std::stringstream s;
            s << "Residual Vector " << k;
            r_.push_back(boost::shared_ptr<Vector> (new Vector(s.str().c_str(),diag_->nirrep(),diag_->dimpi())));
        }
    }

    for (int k = 0; k < nroot_; k++) {

        double R2 = 0.0;

        for (int h = 0; h < diag_->nirrep(); ++h) {
        
        int dimension = diag_->dimpi()[h];
        if (!dimension) continue;
    
            double** ap = a_->pointer(h);
            double*  lp = l_->pointer(h);
            double*  rp = r_[k]->pointer(h);

            ::memset((void*)rp, '\0', dimension*sizeof(double));

            for (int i = 0; i < b_.size(); i++) {
                double* bp = b_[i]->pointer(h);
                double* sp = s_[i]->pointer(h);
                C_DAXPY(dimension,ap[i][k],sp,1,rp,1);
                C_DAXPY(dimension,-lp[k]*ap[i][k],bp,1,rp,1);
            }       

            R2 += C_DDOT(dimension,rp,1,rp,1);

        }

        // Residual norm k
        double rnorm = sqrt(R2);
        n_[k] = rnorm;
        if (rnorm < criteria_) {
            nconverged_++;
        }
    }

    // Global convergence check
    convergence_ = 0.0;
    for (int k = 0; k < nroot_; k++) {
        if (convergence_ < n_[k]) 
            convergence_ = n_[k]; 
    }

    if (nconverged_ == nroot_) converged_ = true; 
    if (debug_) { 
        fprintf(outfile, "   > Residuals <\n\n");
        for (int i = 0; i < r_.size(); i++) {
            r_[i]->print();
        }
        for (int i = 0; i < n_.size(); i++) {
            fprintf(outfile, "    Residual %d = %24.16E\n", i, n_[i]);
        }
        fprintf(outfile,"\n");
        fprintf(outfile,"    %d of %d roots converged, we are %s\n\n",
            nconverged_, nroot_, (converged_ ? "converged": "not converged"));
        fflush(outfile);
    }
}
void DLRSolver::correctors()
{
    // Only add correctors for roots that are not converged
    d_.clear();

    for (int k = 0; k < nroot_; k++) {

        // Do not attempt to add a corrector if root is already converged
        if (n_[k] < criteria_) continue;

        std::stringstream s;
        s << "Corrector Vector " << k;
        boost::shared_ptr<Vector> d(new Vector(s.str().c_str(),diag_->nirrep(),diag_->dimpi()));

        for (int h = 0; h < diag_->nirrep(); ++h) {

            int dimension = diag_->dimpi()[h];
            if (!dimension) continue;

            double* hp = diag_->pointer(h);
            double lambda = E_[k][h];
            double* dp = d->pointer(h);
            double* rp = r_[k]->pointer(h);

            for (int m = 0; m < dimension; m++) {
                dp[m] = rp[m] / (lambda - hp[m]);
            }

            // Substitute r for this vector, if norm is bad
            double norm = sqrt(C_DDOT(dimension, dp, 1, dp, 1));
            if (norm != norm || isinf(norm) || norm < criteria_) {
                C_DCOPY(dimension,rp,1,dp,1);
                norm = sqrt(C_DDOT(dimension, dp, 1, dp, 1));
            } 

            // Normalize the correctors
            C_DSCAL(dimension, 1.0 / norm, dp, 1); 
        }

        d_.push_back(d);
    }
    if (debug_) { 
        fprintf(outfile, "   > Correctors <\n\n");
        for (int i = 0; i < d_.size(); i++) {
            d_[i]->print();
        }
        fflush(outfile);
    }
}
void DLRSolver::subspaceExpansion()
{
    if (debug_) {
        fprintf(outfile, "   > SubspaceExpansion <\n\n");
    }

    // Which vectors are significant?
    std::vector<bool> sig(d_.size());
    for (int i = 0; i < d_.size(); ++i) {
        sig[i] = false;
    }

    // Orthonormalize d_ via Modified Gram-Schmidt
    for (int h = 0; h < diag_->nirrep(); ++h) {

        int dimension = diag_->dimpi()[h];
        if (!dimension) continue;

        // Remove the projection of d on b from b
        for (int i = 0; i < d_.size(); ++i) {
            for (int j = 0; j < b_.size(); ++j) {
                double* dp = d_[i]->pointer(h);
                double* bp = b_[j]->pointer(h);
                
                double r_ji = C_DDOT(dimension,dp,1,bp,1);
                C_DAXPY(dimension,-r_ji,bp,1,dp,1);
            } 
        }

        // Remove the self-projection of d on d from d
        for (int i = 0; i < d_.size(); ++i) {
            double* dip = d_[i]->pointer(h);
            double r_ii = sqrt(C_DDOT(dimension,dip,1,dip,1));
            C_DSCAL(dimension,(r_ii > norm_ ? 1.0 / r_ii : 0.0), dip,1);
            for (int j = i + 1; j < d_.size(); ++j) {
                double* djp = d_[j]->pointer(h);
                double r_ij = C_DDOT(dimension,djp,1,dip,1);
                C_DAXPY(dimension,-r_ij,dip,1,djp,1);
            }  
            if (r_ii > norm_) {
                sig[i] = sig[i] | true;    
            } 
        }
    }

    // Add signinficant vectors
    for (int i = 0; i < d_.size(); ++i) {
        if (sig[i]) {
            b_.push_back(d_[i]);
        }
    } 

    nsubspace_ = b_.size();

    if (debug_) {
        fprintf(outfile, "Final subspace after addition\n\n");
        for (int i = 0; i < b_.size(); i++) {
            b_[i]->print();
        }
        fflush(outfile);
    }
}
void DLRSolver::subspaceCollapse()
{
    if (nsubspace_ <= max_subspace_) return;

    std::vector<boost::shared_ptr<Vector> > s2;
    std::vector<boost::shared_ptr<Vector> > b2;

    for (int k = 0; k < min_subspace_; ++k) {
        std::stringstream bs;
        bs << "Subspace Vector " << k; 
        b2.push_back(boost::shared_ptr<Vector>(new Vector(bs.str(), diag_->nirrep(), diag_->dimpi())));
        std::stringstream ss;
        ss << "Sigma Vector " << k; 
        s2.push_back(boost::shared_ptr<Vector>(new Vector(ss.str(), diag_->nirrep(), diag_->dimpi())));
    }

    int n = a_->rowspi()[0];
    for (int k = 0; k < min_subspace_; ++k) {
        for (int h = 0; h < diag_->nirrep(); ++h) {
            int dimension = diag_->dimpi()[h];
            if (!dimension) continue;
            
            double** ap = a_->pointer(h);
            double*  bp = b_[k]->pointer(h);
            double*  sp = s_[k]->pointer(h);
            double*  b2p = b2[k]->pointer(h);
            double*  s2p = s2[k]->pointer(h);

            for (int i = 0; i < n; ++i) {
                C_DAXPY(dimension,ap[i][k],sp,1,s2p,1);
                C_DAXPY(dimension,ap[i][k],bp,1,b2p,1);
            }
        }        
    }

    s_ = s2;
    b_ = b2;
    nsubspace_ = b_.size();

    if (debug_) {
        fprintf(outfile, "   > SubspaceCollapse <\n\n");
        for (int i = 0; i < b_.size(); i++) {
            b_[i]->print();
        }
        for (int i = 0; i < s_.size(); i++) {
            s_[i]->print();
        }
    }
}

}
