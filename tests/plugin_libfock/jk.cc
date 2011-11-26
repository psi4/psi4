#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libpsio/psio.hpp>
#include <libpsio/aiohandler.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include "sieve.h"
#include "jk.h"
#include"gpudfjkhelper.h"

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;
using namespace psi;

namespace psi {

JK::JK(std::vector<SharedMatrix >& C_left, 
   std::vector<SharedMatrix >& C_right,
   boost::shared_ptr<BasisSet> primary) : 
   C_left_(C_left), C_right_(C_right), lr_symmetric_(false), primary_(primary)
{
    common_init();
}
JK::JK(std::vector<SharedMatrix >& C_symm,
   boost::shared_ptr<BasisSet> primary) :
   C_left_(C_symm), C_right_(C_symm), lr_symmetric_(true), primary_(primary)
{
    common_init();
}
JK::~JK()
{
}
boost::shared_ptr<JK> JK::build_JK(Options& options, bool lr_symmetric)
{
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> primary = BasisSet::construct(parser, Process::environment.molecule(), "BASIS");

    if (options.get_str("SCF_TYPE") == "GPUDF") {

        boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, Process::environment.molecule(), "RI_BASIS_SCF");

        GPUDFJK* jk;

        if (!lr_symmetric) {
            std::vector<SharedMatrix > C_left;
            std::vector<SharedMatrix > C_right;
            jk = new GPUDFJK(C_left, C_right, primary, auxiliary);
        } else {
            std::vector<SharedMatrix > C_left;
            jk = new GPUDFJK(C_left, primary, auxiliary);
        } 
             
        if (options["OMP_N_THREAD"].has_changed())
            jk->set_omp_nthread(options.get_int("OMP_N_THREAD"));
        if (options["SCHWARZ_CUTOFF"].has_changed())
            jk->set_cutoff(options.get_double("SCHWARZ_CUTOFF"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["FITTING_CONDITION"].has_changed())
            jk->set_condition(options.get_double("FITTING_CONDITION"));
        if (options["FITTING_ALGORITHM"].has_changed())
            jk->set_algorithm(options.get_int("FITTING_ALGORITHM"));

        return boost::shared_ptr<JK>(jk);

    } else if (options.get_str("SCF_TYPE") == "DF") {

        boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, Process::environment.molecule(), "RI_BASIS_SCF");

        DFJK* jk;

        if (!lr_symmetric) {
            std::vector<SharedMatrix > C_left;
            std::vector<SharedMatrix > C_right;
            jk = new DFJK(C_left, C_right, primary, auxiliary);
        } else {
            std::vector<SharedMatrix > C_left;
            jk = new DFJK(C_left, primary, auxiliary);
        } 
             
        if (options["OMP_N_THREAD"].has_changed())
            jk->set_omp_nthread(options.get_int("OMP_N_THREAD"));
        if (options["SCHWARZ_CUTOFF"].has_changed())
            jk->set_cutoff(options.get_double("SCHWARZ_CUTOFF"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["FITTING_CONDITION"].has_changed())
            jk->set_condition(options.get_double("FITTING_CONDITION"));
        if (options["FITTING_ALGORITHM"].has_changed())
            jk->set_algorithm(options.get_int("FITTING_ALGORITHM"));

        return boost::shared_ptr<JK>(jk);

    } else if (options.get_str("SCF_TYPE") == "DIRECT") {

        DirectJK* jk;

        if (!lr_symmetric) {
            std::vector<SharedMatrix > C_left;
            std::vector<SharedMatrix > C_right;
            jk = new DirectJK(C_left, C_right, primary);
        } else {
            std::vector<SharedMatrix > C_left;
            jk = new DirectJK(C_left, primary);
        } 
             
        if (options["OMP_N_THREAD"].has_changed())
            jk->set_omp_nthread(options.get_int("OMP_N_THREAD"));
        if (options["SCHWARZ_CUTOFF"].has_changed())
            jk->set_cutoff(options.get_double("SCHWARZ_CUTOFF"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));

        return boost::shared_ptr<JK>(jk);

    } else {
        throw PSIEXCEPTION("JK::build_JK: Unknown SCF Type");
    }
}
void JK::common_init()
{
    print_ = 1;
    debug_ = 0;
    // 256 MB default
    memory_ = 32000000L; 
    omp_nthread_ = 1;
    #ifdef _OPENMP
	omp_nthread_ = omp_get_max_threads();
    #endif
    cutoff_ = 0.0;

    do_J_ = true;
    do_K_ = true;
    do_wK_ = false;

    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(primary_,primary_,primary_,primary_));
    boost::shared_ptr<PetiteList> pet(new PetiteList(primary_, integral));
    AO2USO_ = SharedMatrix(pet->aotoso());
}
unsigned long int JK::memory_overhead()
{
    unsigned long int mem = 0L;

    int JKwKD_factor = 1;
    if (do_J_) JKwKD_factor++;
    if (do_K_) JKwKD_factor++;
    if (do_wK_) JKwKD_factor++;

    int C_factor = 1;
    if (!lr_symmetric_) C_factor++;

    // USO Quantities
    for (int N = 0; N < C_left_.size(); N++) {
        int symml = C_left_[N]->symmetry();
        int symmr = C_right_[N]->symmetry();
        for (int h = 0; h < C_left_[N]->nirrep(); h++) {
            int nbfl = C_left_[N]->rowspi()[h];
            int nbfr = C_right_[N]->rowspi()[h];
            int nocc = C_left_[N]->colspi()[symml^h];

            mem += C_factor * (unsigned long int) nocc * (nbfl + nbfr) / 2L + JKwKD_factor * (unsigned long int) nbfl * nbfr;
        }
    }

    // AO Copies
    if (C1() && C_left_.size() && C_left_[0]->nirrep() != 1) {
        int nbf = primary_->nbf();
        for (int N = 0; N < C_left_.size(); N++) {
            int nocc = 0;
            for (int h = 0; h < C_left_[N]->nirrep(); h++) {
                nocc += C_left_[N]->colspi()[h];
            }
            mem += C_factor * (unsigned long int) nocc * nbf + JKwKD_factor * (unsigned long int) nbf * nbf;
        }
    }

    return mem;
}
void JK::compute_D()
{
    /// Make sure the memory is there
    if (D_.size() != C_left_.size()) {
        D_.clear();
	for (int N = 0; N < C_left_.size(); ++N) {
            std::stringstream s;
	    s << "D " << N << " (SO)";
            D_.push_back(SharedMatrix(new Matrix(s.str(),C_left_[N]->nirrep(), C_left_[N]->rowspi(), C_right_[N]->rowspi(), C_left_[N]->symmetry() ^ C_right_[N]->symmetry())));
	}
    }

    for (int N = 0; N < D_.size(); ++N) {
        int symm = D_[N]->symmetry();
        for (int h = 0; h < D_[N]->nirrep(); ++h) {

            int nsol = C_left_[N]->rowspi()[h];
            int nocc = C_left_[N]->colspi()[h];
            int nsor = C_right_[N]->rowspi()[h^symm];

            if (!nsol || !nsor || !nocc) continue;

            double** Dp = D_[N]->pointer(h^symm);
            double** Clp = C_left_[N]->pointer(h);
            double** Crp = C_right_[N]->pointer(h^symm);

            C_DGEMM('N','T', nsol, nsor, nocc, 1.0, Clp[0], nocc, Crp[0], nocc, 0.0, Dp[0], nsor);
        }
    }
}
void JK::allocate_JK()
{
    // Allocate J/K in the case that the algorithm uses USOs, so AO2USO will not allocate. 
    if (J_.size() != D_.size()) {
	J_.clear();    
	K_.clear();    
	wK_.clear();    
    	for (int N = 0; N < D_.size() && do_J_; ++N) {
            std::stringstream s;
	    s << "J " << N << " (SO)";
            J_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(), D_[N]->symmetry())));
    	}
    	for (int N = 0; N < D_.size() && do_K_; ++N) {
            std::stringstream s;
	    s << "K " << N << " (SO)";
            K_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(), D_[N]->symmetry())));
    	}
    	for (int N = 0; N < D_.size() && do_wK_; ++N) {
            std::stringstream s;
	    s << "wK " << N << " (SO)";
            wK_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(), D_[N]->symmetry())));
    	}
    }
    
    // Zero out J/K for compute_JK()
    for (int N = 0; N < D_.size(); ++N) {
        if (do_J_) J_[N]->zero();
        if (do_K_) K_[N]->zero();
        if (do_wK_) wK_[N]->zero();
    }	
}
void JK::USO2AO()
{
    // If C1, C_ao and D_ao are equal to C and D
    if (AO2USO_->nirrep() == 1) {
        C_left_ao_ = C_left_;
        C_right_ao_ = C_right_;
        D_ao_ = D_;
    }

    // Whether C1 or not, allocate J_ao and K_ao
    if (D_.size() != J_ao_.size()) {
	J_ao_.clear();    
	K_ao_.clear();    
	wK_ao_.clear();    

    	for (int N = 0; N < D_.size() && do_J_; ++N) {
            std::stringstream s;
	    s << "J " << N << " (AO)";
            J_ao_.push_back(SharedMatrix(new Matrix(s.str(),AO2USO_->rowspi()[0], AO2USO_->rowspi()[0])));
    	}
    	for (int N = 0; N < D_.size() && do_K_; ++N) {
            std::stringstream s;
	    s << "K " << N << " (AO)";
            K_ao_.push_back(SharedMatrix(new Matrix(s.str(),AO2USO_->rowspi()[0], AO2USO_->rowspi()[0])));
    	}
    	for (int N = 0; N < D_.size() && do_wK_; ++N) {
            std::stringstream s;
	    s << "wK " << N << " (AO)";
            wK_ao_.push_back(SharedMatrix(new Matrix(s.str(),AO2USO_->rowspi()[0], AO2USO_->rowspi()[0])));
    	}
    }	    

    // Zero out J_ao/K_ao for compute_JK()
    for (int N = 0; N < D_.size(); ++N) {
        if (do_J_) J_ao_[N]->zero();
        if (do_K_) K_ao_[N]->zero();
        if (do_wK_) K_ao_[N]->zero();
    }	

    if (AO2USO_->nirrep() == 1) return;

    // If not C1, allocate C_ao and D_ao
    if (D_.size() != C_left_ao_.size()) {
	C_left_ao_.clear();    
	C_right_ao_.clear();    
	D_ao_.clear();    
    	for (int N = 0; N < D_.size(); ++N) {
            std::stringstream s;
	    s << "D " << N << " (AO)";
            D_ao_.push_back(SharedMatrix(new Matrix(s.str(),AO2USO_->rowspi()[0], AO2USO_->rowspi()[0])));
    	}
    	for (int N = 0; N < D_.size(); ++N) {
            std::stringstream s;
	    s << "C Left " << N << " (AO)";
	    int ncol = 0;
	    for (int h = 0; h < C_left_[N]->nirrep(); ++h) ncol += C_left_[N]->colspi()[h];
            C_left_ao_.push_back(SharedMatrix(new Matrix(s.str(),AO2USO_->rowspi()[0], ncol)));
    	}
    	for (int N = 0; (N < D_.size()) && (!lr_symmetric_); ++N) {
            std::stringstream s;
	    s << "C Right " << N << " (AO)";
	    int ncol = 0;
	    for (int h = 0; h < C_right_[N]->nirrep(); ++h) ncol += C_right_[N]->colspi()[h];
            C_right_ao_.push_back(SharedMatrix(new Matrix(s.str(),AO2USO_->rowspi()[0], ncol)));
    	}
    }	    

    // If not C1, check for changes in nocc
    for (int N = 0; N < D_.size(); ++N) {
	int ncol = 0;
	for (int h = 0; h < C_left_[N]->nirrep(); ++h) ncol += C_left_[N]->colspi()[h];
	int ncol_old = C_left_ao_[N]->colspi()[0];
	if (ncol != ncol_old) {
            std::stringstream s;
	    s << "C Left " << N << " (AO)";
            C_left_ao_[N] = SharedMatrix(new Matrix(s.str(),AO2USO_->rowspi()[0], ncol));
	    if (!lr_symmetric_) {
                std::stringstream s2;
	        s2 << "C Right " << N << " (AO)";
                C_right_ao_[N] = SharedMatrix(new Matrix(s2.str(),AO2USO_->rowspi()[0], ncol));
			
            }
	}
	if (lr_symmetric_) C_right_ao_ = C_left_ao_;
    }	    

    // Transform D
    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    for (int N = 0; N < D_.size(); ++N) {
	D_ao_[N]->zero();    
        int symm = D_[N]->symmetry();
	for (int h = 0; h < AO2USO_->nirrep(); ++h) {    
            int nao = AO2USO_->rowspi()[0];
            int nsol = AO2USO_->colspi()[h];
            int nsor = AO2USO_->colspi()[h^symm];
	    if (!nsol || !nsor) continue;
	    double** Ulp = AO2USO_->pointer(h);
	    double** Urp = AO2USO_->pointer(h^symm);
	    double** DSOp = D_[N]->pointer(h^symm);
	    double** DAOp = D_ao_[N]->pointer();
	    C_DGEMM('N','T',nsol,nao,nsor,1.0,DSOp[0],nsor,Urp[0],nsor,0.0,temp,nao);
	    C_DGEMM('N','N',nao,nao,nsol,1.0,Ulp[0],nsol,temp,nao,1.0,DAOp[0],nao);
	}
    }	    
    delete[] temp;	    

    // Transform C
    for (int N = 0; N < D_.size(); ++N) {
        int offset = 0;
	for (int h = 0; h < AO2USO_->nirrep(); ++h) {
            int nao = AO2USO_->rowspi()[0];
            int nso = AO2USO_->colspi()[h];
	    int ncol = C_left_ao_[N]->colspi()[0];
	    int ncolspi = C_left_[N]->colspi()[h];
	    if (nso == 0 || ncolspi == 0) continue;
	    double** Up = AO2USO_->pointer(h);
	    double** CAOp = C_left_ao_[N]->pointer();
	    double** CSOp = C_left_[N]->pointer(h);
	    C_DGEMM('N','N',nao,ncolspi,nso,1.0,Up[0],nso,CSOp[0],ncolspi,0.0,&CAOp[0][offset],ncol);
	    offset += ncolspi;
	}
    }	    
    for (int N = 0; (N < D_.size()) && (!lr_symmetric_); ++N) {
        int offset = 0;
        int symm = D_[N]->symmetry();
	for (int h = 0; h < AO2USO_->nirrep(); ++h) {
            int nao = AO2USO_->rowspi()[0];
            int nso = AO2USO_->colspi()[h];
	    int ncol = C_right_ao_[N]->colspi()[0];
	    int ncolspi = C_right_[N]->colspi()[h^symm];
	    if (nso == 0 || ncolspi == 0) continue;
	    double** Up = AO2USO_->pointer(h);
	    double** CAOp = C_right_ao_[N]->pointer();
	    double** CSOp = C_right_[N]->pointer(h);
	    C_DGEMM('N','N',nao,ncolspi,nso,1.0,Up[0],nso,CSOp[0],ncolspi,0.0,&CAOp[0][offset],ncol);
	    offset += ncolspi;
	}
    }	    
}
void JK::AO2USO()
{
    // If already C1, J/K are J_ao/K_ao
    if (AO2USO_->nirrep() == 1) {
	J_ = J_ao_;	    
	K_ = K_ao_;	    
	wK_ = wK_ao_;	    
	return;
    }

    // If not C1, J/K must be allocated
    if (J_.size() != J_ao_.size()) {
	J_.clear();    
	K_.clear();    
	wK_.clear();    
    	for (int N = 0; N < D_.size() && do_J_; ++N) {
            std::stringstream s;
	    s << "J " << N << " (SO)";
            J_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(), D_[N]->symmetry())));
    	}
    	for (int N = 0; N < D_.size() && do_K_; ++N) {
            std::stringstream s;
	    s << "K " << N << " (SO)";
            K_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(), D_[N]->symmetry())));
    	}
    	for (int N = 0; N < D_.size() && do_wK_; ++N) {
            std::stringstream s;
	    s << "wK " << N << " (SO)";
            wK_.push_back(SharedMatrix(new Matrix(s.str(),D_[N]->nirrep(), D_[N]->rowspi(), D_[N]->rowspi(), D_[N]->symmetry())));
    	}
    }

    // Transform
    double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
    for (int N = 0; N < D_.size(); ++N) {
        int symm = D_[N]->symmetry();
	for (int h = 0; h < AO2USO_->nirrep(); ++h) {    
            int nao = AO2USO_->rowspi()[0];
            int nsol = AO2USO_->colspi()[h];
            int nsor = AO2USO_->colspi()[h^symm];

	    if (!nsol || !nsor) continue;

	    double** Ulp = AO2USO_->pointer(h);
	    double** Urp = AO2USO_->pointer(h^symm);

            if (do_J_) {
	        double** JAOp = J_ao_[N]->pointer();
	        double** JSOp = J_[N]->pointer(h);
	        C_DGEMM('N','N',nao,nsor,nao,1.0,JAOp[0],nao,Urp[0],nsor,0.0,temp,nsor);
	        C_DGEMM('T','N',nsol,nsor,nao,1.0,Ulp[0],nsol,temp,nsor,0.0,JSOp[0],nsor);
            }
            if (do_K_) {
	        double** KAOp = K_ao_[N]->pointer();
	        double** KSOp = K_[N]->pointer(h);
	        C_DGEMM('N','N',nao,nsor,nao,1.0,KAOp[0],nao,Urp[0],nsor,0.0,temp,nsor);
	        C_DGEMM('T','N',nsol,nsor,nao,1.0,Ulp[0],nsol,temp,nsor,0.0,KSOp[0],nsor);
            }
            if (do_wK_) {
	        double** wKAOp = wK_ao_[N]->pointer();
	        double** wKSOp = wK_[N]->pointer(h);
	        C_DGEMM('N','N',nao,nsor,nao,1.0,wKAOp[0],nao,Urp[0],nsor,0.0,temp,nsor);
	        C_DGEMM('T','N',nsol,nsor,nao,1.0,Ulp[0],nsol,temp,nsor,0.0,wKSOp[0],nsor);
            }
	}
    }
    delete[] temp;
}
void JK::initialize()
{
    preiterations();
}
void JK::compute()
{
    compute_D();
    if (C1()) USO2AO();
    else allocate_JK();
    compute_JK();
    if (C1()) AO2USO();

    if (debug_ > 3) {
        fprintf(outfile, "   > JK <\n\n");
        for (int N = 0; N < C_left_.size(); N++) {
            C_left_ao_[N]->print(outfile);
            C_left_[N]->print(outfile);
            C_right_ao_[N]->print(outfile);
            C_right_[N]->print(outfile);
            D_ao_[N]->print(outfile);
            D_[N]->print(outfile);
            J_ao_[N]->print(outfile);
            J_[N]->print(outfile);
            K_ao_[N]->print(outfile);
            K_[N]->print(outfile);
        }
        fflush(outfile);
    }
}
void JK::finalize()
{
    postiterations();
}

DirectJK::DirectJK(std::vector<SharedMatrix >& C_left, 
   std::vector<SharedMatrix >& C_right,
   boost::shared_ptr<BasisSet> primary) : 
   JK(C_left,C_right,primary)	
{
    common_init();
}
DirectJK::DirectJK(std::vector<SharedMatrix >& C_symm,
   boost::shared_ptr<BasisSet> primary) :
   JK(C_symm,primary)
{
    common_init();
}
DirectJK::~DirectJK()
{
}
void DirectJK::common_init()
{
}
void DirectJK::print_header() const
{
    if (print_) {
    	fprintf(outfile, "  ==> DirectJK: Integral-Direct J/K Matrices <==\n");
	fprintf(outfile, "                  by Rob Parrish\n\n");

        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        fprintf(outfile, "    OpenMP threads:    %11d\n", omp_nthread_);
        fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n\n", cutoff_);
    }
}
void DirectJK::preiterations()
{
    sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));
    factory_= boost::shared_ptr<IntegralFactory>(new IntegralFactory(primary_,primary_,primary_,primary_));
    eri_.clear();
    for (int thread = 0; thread < omp_nthread_; thread++) {
        eri_.push_back(boost::shared_ptr<TwoBodyAOInt>(factory_->eri()));        
    } 
}
void DirectJK::compute_JK()
{
    // Correctness always counts
    const double* buffer = eri_[0]->buffer();
    for (int M = 0; M < primary_->nshell(); ++M) {
    for (int N = 0; N < primary_->nshell(); ++N) {
    for (int R = 0; R < primary_->nshell(); ++R) {
    for (int S = 0; S < primary_->nshell(); ++S) {
        
        eri_[0]->compute_shell(M,N,R,S);
        
        int nM = primary_->shell(M)->nfunction();
        int nN = primary_->shell(N)->nfunction();
        int nR = primary_->shell(R)->nfunction();
        int nS = primary_->shell(S)->nfunction();

        int sM = primary_->shell(M)->function_index();
        int sN = primary_->shell(N)->function_index();
        int sR = primary_->shell(R)->function_index();
        int sS = primary_->shell(S)->function_index();
     
        for (int oM = 0, index = 0; oM < nM; oM++) {
        for (int oN = 0; oN < nN; oN++) {
        for (int oR = 0; oR < nR; oR++) {
        for (int oS = 0; oS < nS; oS++, index++) {

            double val = buffer[index];
    
            int m = oM + sM;     
            int n = oN + sN;     
            int r = oR + sR;     
            int s = oS + sS;     
   
            if (do_J_) {
                for (int N = 0; N < J_ao_.size(); N++) {
                    J_ao_[N]->add(0,m,n, D_ao_[N]->get(0,r,s)*val);
                }
            }
   
            if (do_K_) {
                for (int N = 0; N < K_ao_.size(); N++) {
                    K_ao_[N]->add(0,m,s, D_ao_[N]->get(0,n,r)*val);
                }
            }

        }}}} 

    }}}}

    // Faster version, not finished
    /**
    sieve_->set_sieve(cutoff_);
    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    unsigned long int nMN = shell_pairs.size();
    unsigned long int nMNRS = nMN * nMN;
    int nthread = eri_.size();

    #pragma omp parallel for schedule(dynamic,30) num_threads(nthread)
    for (unsigned long int index = 0L; index < nMNRS; ++index) {

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        const double* buffer = eri_[thread]->buffer();

        unsigned long int MN = index / nMN;        
        unsigned long int RS = index % nMN;
        if (MN < RS) continue;       

        int M = shell_pairs[MN].first;
        int N = shell_pairs[MN].second;
        int R = shell_pairs[RS].first;
        int S = shell_pairs[RS].second;

        eri_[thread]->compute_shell(M,N,R,S);
        
        int nM = primary_->shell(M)->nfunction();
        int nN = primary_->shell(N)->nfunction();
        int nR = primary_->shell(R)->nfunction();
        int nS = primary_->shell(S)->nfunction();

        int sM = primary_->shell(M)->function_index();
        int sN = primary_->shell(N)->function_index();
        int sR = primary_->shell(R)->function_index();
        int sS = primary_->shell(S)->function_index();

        for (int oM = 0, index = 0; oM < nM; oM++) {
        for (int oN = 0; oN < nN; oN++) {
        for (int oR = 0; oR < nR; oR++) {
        for (int oS = 0; oS < nS; oS++, index++) {

            int m = oM + sM;     
            int n = oN + sN;     
            int r = oR + sR;     
            int s = oS + sS;     

            if ((n > m) || (s > r) || ((r*(r+1) >> 1) + s > (m*(m+1) >> 1) + n)) continue;            

            double val = buffer[index];

            if (do_J_) {
                for (int N = 0; N < J_ao_.size(); N++) {
                    double** Dp = D_ao_[N]->pointer();
                    double** Jp = J_ao_[N]->pointer();

                    // I've given you all the unique ones
                    // Make sure to use #pragma omp atomic 
                    // TODO
                }
            }

            if (do_K_) {
                for (int N = 0; N < K_ao_.size(); N++) {
                    double** Dp = D_ao_[N]->pointer();
                    double** Kp = J_ao_[N]->pointer();

                    // I've given you all the unique ones
                    // Make sure to use #pragma omp atomic 
                    // TODO
                }
            }

        }}}} 
         
    }
    **/
}
void DirectJK::postiterations()
{
    sieve_.reset();
    factory_.reset();
    eri_.clear();
}

DFJK::DFJK(std::vector<SharedMatrix >& C_left, 
   std::vector<SharedMatrix >& C_right,
   boost::shared_ptr<BasisSet> primary, 
   boost::shared_ptr<BasisSet> auxiliary) :
   JK(C_left,C_right,primary), auxiliary_(auxiliary)	
{
    common_init();
}
DFJK::DFJK(std::vector<SharedMatrix >& C_symm,
   boost::shared_ptr<BasisSet> primary,
   boost::shared_ptr<BasisSet> auxiliary) :
   JK(C_symm,primary), auxiliary_(auxiliary)	
{
    common_init();
}
DFJK::~DFJK()
{
}
void DFJK::common_init()
{
    algorithm_ = 0;
    condition_ = 1.0E-12;
    unit_ = PSIF_DFSCF_BJ;
    is_core_ = true;
    psio_ = PSIO::shared_object();
}
void DFJK::print_header() const
{
    if (print_) {
    	fprintf(outfile, "  ==> DFJK: Density-Fitted J/K Matrices <==\n");
	fprintf(outfile, "             by Rob Parrish\n\n");

        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        fprintf(outfile, "    OpenMP threads:    %11d\n", omp_nthread_);
        fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        fprintf(outfile, "    Algorithm:         %11s\n",  (is_core_ ? "Core" : "Disk")); 
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        fprintf(outfile, "    Fitting Condition: %11.0E\n\n", condition_);

        fprintf(outfile, "   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level(outfile, print_);
    }
}
bool DFJK::is_core()
{
    int ntri = sieve_->function_pairs().size();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    // Two is for buffer space in fitting
    return (three_memory + 2L*two_memory < memory_); 
}
unsigned long int DFJK::memory_temp()
{
    unsigned long int mem = 0L;

    // J Overhead (Jtri, Dtri, d)
    mem += 2L * sieve_->function_pairs().size() + auxiliary_->nbf();
    // K Overhead (C_temp, Q_temp)
    mem += omp_nthread_ * (unsigned long int) primary_->nbf() * (auxiliary_->nbf() + max_nocc());

    // Sort temp in K
    mem += (algorithm_ == 1 ? 1L : 0L) * primary_->nbf() * (ULI) primary_->nbf(); 

    return mem;  
} 
int DFJK::max_rows()
{
    // Start with all memory
    unsigned long int mem = memory_;
    // Subtract J/K/C/D overhead
    mem -= memory_overhead();
    // Subtract threading temp overhead
    mem -= memory_temp();
    
    // How much will each row cost?
    unsigned long int row_cost = 0L;
    // Copies of E tensor
    row_cost += (lr_symmetric_ ? 1L : 2L) * max_nocc() * primary_->nbf();
    // Slices of Qmn tensor, including AIO buffer (NOTE: AIO not implemented yet)
    row_cost += (is_core_ ? 1L : 1L) * sieve_->function_pairs().size(); 
    // Slices of unpacked Qmn tensor if algorithm2
    row_cost += (algorithm_ == 1 ? 1L : 0L) * primary_->nbf() * (ULI) primary_->nbf(); 
    
    unsigned long int max_rows = mem / row_cost;

    if (max_rows > (unsigned long int) auxiliary_->nbf())
        max_rows = (unsigned long int) auxiliary_->nbf();
    if (max_rows < 1L)
        max_rows = 1L;

    return (int) max_rows;
}
int DFJK::max_nocc()
{
    int max_nocc = 0;
    for (int N = 0; N < C_left_ao_.size(); N++) {
        max_nocc = (C_left_ao_[N]->colspi()[0] > max_nocc ? C_left_ao_[N]->colspi()[0] : max_nocc);
    }
    return max_nocc;
}
void DFJK::initialize_temps()
{
    J_temp_ = boost::shared_ptr<Vector>(new Vector("Jtemp", sieve_->function_pairs().size()));
    D_temp_ = boost::shared_ptr<Vector>(new Vector("Dtemp", sieve_->function_pairs().size()));
    d_temp_ = boost::shared_ptr<Vector>(new Vector("dtemp", max_rows_)); 
    
    for (int thread = 0; thread < omp_nthread_; thread++) {
        C_temp_.push_back(SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf())));
        Q_temp_.push_back(SharedMatrix(new Matrix("Qtemp", max_rows_, primary_->nbf())));
    } 

    E_left_ = SharedMatrix(new Matrix("E_left", primary_->nbf(), max_rows_ * max_nocc_));
    if (lr_symmetric_) 
        E_right_ = E_left_;
    else 
        E_right_ = boost::shared_ptr<Matrix>(new Matrix("E_right", primary_->nbf(), max_rows_ * max_nocc_));

    if (algorithm_ == 1) {
        Qmn2_ = boost::shared_ptr<Matrix>(new Matrix("Qmn2", max_rows_, primary_->nbf() * (ULI) primary_->nbf()));
        sort_ = boost::shared_ptr<Vector>(new Vector("Sort", primary_->nbf() * (ULI) primary_->nbf()));
    }
}
void DFJK::free_temps()
{
    J_temp_.reset(); 
    D_temp_.reset(); 
    d_temp_.reset(); 
    E_left_.reset(); 
    E_right_.reset(); 
    C_temp_.clear(); 
    Q_temp_.clear(); 
    Qmn2_.reset();
    sort_.reset();
}
void DFJK::preiterations()
{
    // DF requires constant sieve, must be static througout object life
    if (!sieve_) {
        sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));
    }

    // Core or disk?
    is_core_ =  is_core();

    if (is_core_)
        initialize_JK_core();
    else 
        initialize_JK_disk();
}
void DFJK::compute_JK()
{
    max_nocc_ = max_nocc();
    max_rows_ = max_rows();
    initialize_temps();
    if (is_core_) 
        manage_JK_core();
    else
        manage_JK_disk();
    free_temps();
}
void DFJK::postiterations()
{
    psio_->close(unit_, 1);
    Qmn_.reset();
}
void DFJK::initialize_JK_core()
{
    int ntri = sieve_->function_pairs().size();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_nthread_;
    #endif
    int rank = 0;

    Qmn_ = SharedMatrix(new Matrix("Qmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qmn_->pointer();

    // Try to load
    // TODO
    if (false) {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        return;
    } 
    
    //Get a TEI for each thread
    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }

    const std::vector<long int>& schwarz_shell_pairs = sieve_->shell_pairs_reverse();
    const std::vector<long int>& schwarz_fun_pairs = sieve_->function_pairs_reverse();

    int numP,Pshell,MU,NU,P,PHI,mu,nu,nummu,numnu,omu,onu;
    //The integrals (A|mn)
    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_->shell(MU)->nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU)->nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell)->nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU)->function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU)->function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell)->function_index() + P;
                                    Qmnp[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    delete []buffer;
    delete []eri;

    boost::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true)); 
    Jinv->form_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    ULI max_cols = (memory_-three_memory-two_memory) / auxiliary_->nbf();
    if (max_cols < 1)
        max_cols = 1;
    if (max_cols > ntri)
        max_cols = ntri;
    SharedMatrix temp(new Matrix("Qmn buffer", auxiliary_->nbf(), max_cols));
    double** tempp = temp->pointer();

    int nblocks = ntri / max_cols;
    if ((ULI)nblocks*max_cols != ntri) nblocks++;

    int ncol = 0;
    int col = 0;
    for (int block = 0; block < nblocks; block++) {

        ncol = max_cols;
        if (col + ncol > ntri)
            ncol = ntri - col;

        C_DGEMM('N','N',auxiliary_->nbf(), ncol, auxiliary_->nbf(), 1.0,
            Jinvp[0], auxiliary_->nbf(), &Qmnp[0][col], ntri, 0.0,
            tempp[0], max_cols);

        for (int Q = 0; Q < auxiliary_->nbf(); Q++) {
            C_DCOPY(ncol, tempp[Q], 1, &Qmnp[Q][col], 1);
        }

        col += ncol;
    }

    //Qmn_->print();

    psio_->open(unit_,PSIO_OPEN_NEW);
    psio_->write_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
}
void DFJK::initialize_JK_disk()
{
    // Try to load
    // TODO

    int nso = primary_->nbf();
    int nshell = primary_->nshell();
    int naux = auxiliary_->nbf();
    
    // ==> Schwarz Indexing <== //
    const std::vector<std::pair<int,int> >& schwarz_shell_pairs = sieve_->shell_pairs();
    const std::vector<std::pair<int,int> >& schwarz_fun_pairs = sieve_->function_pairs();
    int nshellpairs = schwarz_shell_pairs.size();
    int ntri = schwarz_fun_pairs.size(); 
    const std::vector<long int>&  schwarz_shell_pairs_r = sieve_->shell_pairs_reverse();
    const std::vector<long int>&  schwarz_fun_pairs_r = sieve_->function_pairs_reverse();

    // ==> Memory Sizing <== //
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI buffer_memory = memory_ - 2*two_memory; // Two is for buffer space in fitting

    //fprintf(outfile, "Buffer memory = %ld words\n", buffer_memory);

    //fprintf(outfile,"Schwarz Shell Pairs:\n");
    //for (int MN = 0; MN < nshellpairs; MN++) {
    //    fprintf(outfile,"  %3d: (%3d,%3d)\n", MN, schwarz_shell_pairs[2*MN], schwarz_shell_pairs[2*MN + 1]);   
    //}

    //fprintf(outfile,"Schwarz Function Pairs:\n");
    //for (int MN = 0; MN < ntri; MN++) {
    //    fprintf(outfile,"  %3d: (%3d,%3d)\n", MN, schwarz_fun_pairs[2*MN], schwarz_fun_pairs[2*MN + 1]);   
    //}

    //fprintf(outfile,"Schwarz Reverse Shell Pairs:\n");
    //for (int MN = 0; MN < primary_->nshell() * (primary_->nshell() + 1) / 2; MN++) {
    //    fprintf(outfile,"  %3d: %4ld\n", MN, schwarz_shell_pairs_r[MN]);
    //}

    //fprintf(outfile,"Schwarz Reverse Function Pairs:\n");
    //for (int MN = 0; MN < primary_->nbf() * (primary_->nbf() + 1) / 2; MN++) {
    //    fprintf(outfile,"  %3d: %4ld\n", MN, schwarz_fun_pairs_r[MN]);
    //}

    // Find out exactly how much memory per MN shell 
    boost::shared_ptr<IntVector> MN_mem(new IntVector("Memory per MN pair", nshell * (nshell + 1) / 2));
    int *MN_memp = MN_mem->pointer();

    for (int mn = 0; mn < ntri; mn++) { 
        int m = schwarz_fun_pairs[mn].first;
        int n = schwarz_fun_pairs[mn].second;

        int M = primary_->function_to_shell(m);
        int N = primary_->function_to_shell(n);
        
        MN_memp[M * (M + 1) / 2 + N] += naux; 
    } 
  
    //MN_mem->print(outfile);
 
    // Figure out exactly how much memory per M row  
    ULI* M_memp = new ULI[nshell];
    memset(static_cast<void*>(M_memp), '\0', nshell*sizeof(ULI));

    for (int M = 0; M < nshell; M++) {
        for (int N = 0; N <= M; N++) {
            M_memp[M] += MN_memp[M * (M + 1) / 2 + N];
        }
    } 

    //fprintf(outfile,"  # Memory per M row #\n\n");
    //for (int M = 0; M < nshell; M++) 
    //    fprintf(outfile,"   %3d: %10ld\n", M+1,M_memp[M]);
    //fprintf(outfile,"\n");

    // Find and check the minimum required memory for this problem
    ULI min_mem = naux*(ULI) ntri;
    for (int M = 0; M < nshell; M++) {
        if (min_mem > M_memp[M])
            min_mem = M_memp[M];
    }

    if (min_mem > buffer_memory) {
        std::stringstream message;
        message << "SCF::DF: Disk based algorithm requires 2 (A|B) fitting metrics and an (A|mn) chunk on core." << std::endl;
        message << "         This is 2Q^2 + QNP doubles, where Q is the auxiliary basis size, N is the" << std::endl;       
        message << "         primary basis size, and P is the maximum number of functions in a primary shell." << std::endl; 
        message << "         For this problem, that is " << ((8L*(min_mem + 2*two_memory))) << " bytes before taxes,";   
        message << ((80L*(min_mem + 2*two_memory) / 7L)) << " bytes after taxes. " << std::endl;                         
    
        throw PSIEXCEPTION(message.str());           
    } 
   
    // ==> Reduced indexing by M <== //
 
    // Figure out the MN start index per M row
    boost::shared_ptr<IntVector> MN_start(new IntVector("MUNU start per M row", nshell));
    int* MN_startp = MN_start->pointer();
    
    MN_startp[0] = schwarz_shell_pairs_r[0];
    int M_index = 1;
    for (int MN = 0; MN < nshellpairs; MN++) {
        if (schwarz_shell_pairs[MN].first == M_index) {
            MN_startp[M_index] = MN;
            M_index++;
        }
    }

    // Figure out the mn start index per M row
    boost::shared_ptr<IntVector> mn_start(new IntVector("munu start per M row", nshell));
    int* mn_startp = mn_start->pointer();
    
    mn_startp[0] = schwarz_fun_pairs[0].first;
    int m_index = 1;
    for (int mn = 0; mn < ntri; mn++) {
        if (primary_->function_to_shell(schwarz_fun_pairs[mn].first) == m_index) {
            mn_startp[m_index] = mn;
            m_index++;
        }
    }

    // Figure out the MN columns per M row 
    boost::shared_ptr<IntVector> MN_col(new IntVector("MUNU cols per M row", nshell));
    int* MN_colp = MN_col->pointer();
   
    for (int M = 1; M < nshell; M++) {
        MN_colp[M - 1] = MN_startp[M] - MN_startp[M - 1]; 
    }
    MN_colp[nshell - 1] = nshellpairs - MN_startp[nshell - 1];
    
    // Figure out the mn columns per M row 
    boost::shared_ptr<IntVector> mn_col(new IntVector("munu cols per M row", nshell));
    int* mn_colp = mn_col->pointer();
   
    for (int M = 1; M < nshell; M++) {
        mn_colp[M - 1] = mn_startp[M] - mn_startp[M - 1]; 
    }
    mn_colp[nshell - 1] = ntri - mn_startp[nshell - 1];

    //MN_start->print(outfile);
    //MN_col->print(outfile);
    //mn_start->print(outfile);
    //mn_col->print(outfile);
 
    // ==> Block indexing <== // 
    // Sizing by block 
    std::vector<int> MN_start_b;
    std::vector<int> MN_col_b;
    std::vector<int> mn_start_b;
    std::vector<int> mn_col_b; 

    // Determine MN and mn block starts
    // also MN and mn block cols
    int nblock = 1;
    ULI current_mem = 0L;
    MN_start_b.push_back(0);
    mn_start_b.push_back(0);
    MN_col_b.push_back(0);
    mn_col_b.push_back(0);
    for (int M = 0; M < nshell; M++) {
        if (current_mem + M_memp[M] > buffer_memory) {
            MN_start_b.push_back(MN_startp[M]);
            mn_start_b.push_back(mn_startp[M]);
            MN_col_b.push_back(0);
            mn_col_b.push_back(0);
            nblock++;
            current_mem = 0L;
        }
        MN_col_b[nblock - 1] += MN_colp[M];
        mn_col_b[nblock - 1] += mn_colp[M];
        current_mem += M_memp[M];
    } 

    //fprintf(outfile,"Block, MN start, MN cols, mn start, mn cols\n"); 
    //for (int block = 0; block < nblock; block++) {
    //    fprintf(outfile,"  %3d: %12d %12d %12d %12d\n", block, MN_start_b[block], MN_col_b[block], mn_start_b[block], mn_col_b[block]);  
    //}
    //fflush(outfile);
 
    // Full sizing not required any longer 
    MN_mem.reset();
    MN_start.reset();
    MN_col.reset(); 
    mn_start.reset();
    mn_col.reset(); 
    delete[] M_memp;     

    // ==> Buffer allocation <== //
    int max_cols = 0;
    for (int block = 0; block < nblock; block++) {
        if (max_cols < mn_col_b[block])
            max_cols = mn_col_b[block];
    }

    // Primary buffer
    Qmn_ = SharedMatrix(new Matrix("(Q|mn) (Disk Chunk)", naux, max_cols));
    // Fitting buffer
    SharedMatrix Amn (new Matrix("(Q|mn) (Buffer)",naux,naux));
    double** Qmnp = Qmn_->pointer();
    double** Amnp = Amn->pointer();

    // ==> Prestripe/Jinv <== // 
    psio_->open(unit_,PSIO_OPEN_NEW);
    boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

    // Dispatch the prestripe
    aio->zero_disk(unit_,"(Q|mn) Integrals",naux,ntri);
    
    // Form the J symmetric inverse
    boost::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true)); 
    Jinv->form_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    // Synch up
    aio->synchronize();

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        nthread = omp_nthread_;
    #endif

    // ==> ERI initialization <== //
    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->eri());
        buffer[Q] = eri[Q]->buffer();
    }
    
    // ==> Main loop <== //
    for (int block = 0; block < nblock; block++) {
        int MN_start_val = MN_start_b[block]; 
        int mn_start_val = mn_start_b[block]; 
        int MN_col_val = MN_col_b[block]; 
        int mn_col_val = mn_col_b[block]; 
    
        // ==> (A|mn) integrals <== //
        #pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {

            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int MU = schwarz_shell_pairs[MUNU + 0].first;
            int NU = schwarz_shell_pairs[MUNU + 0].second;
            int nummu = primary_->shell(MU)->nfunction();
            int numnu = primary_->shell(NU)->nfunction();
            int mu = primary_->shell(MU)->function_index();
            int nu = primary_->shell(NU)->function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P)->nfunction();
                int p = auxiliary_->shell(P)->function_index();
                eri[rank]->compute_shell(P,0,MU,NU);
                for (int dm = 0; dm < nummu; dm++) {
                    int omu = mu + dm;
                    for (int dn = 0; dn < numnu;  dn++) {
                        int onu = nu + dn;
                        if (omu >= onu && schwarz_fun_pairs_r[omu*(omu+1)/2 + onu] >= 0) {
                            int delta = schwarz_fun_pairs_r[omu*(omu+1)/2 + onu] - mn_start_val;
                            for (int dp = 0; dp < nump; dp ++) {
                                int op = p + dp;
                                Qmnp[op][delta] = buffer[rank][dp*nummu*numnu + dm*numnu + dn];
                            }
                        }
                    }
                }
            }
        }

        // ==> (Q|mn) fitting <== //
        for (int mn = 0; mn < mn_col_val; mn+=naux) {
            int cols = naux;
            if (mn + naux >= mn_col_val)
                cols = mn_col_val - mn;

            for (int Q = 0; Q < naux; Q++)
                C_DCOPY(cols,&Qmnp[Q][mn],1,Amnp[Q],1);

            C_DGEMM('N','N',naux,cols,naux,1.0,Jinvp[0],naux,Amnp[0],naux,0.0,&Qmnp[0][mn],max_cols);
        }

        // ==> Disk striping <== //
        psio_address addr;
        for (int Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri + mn_start_val)*sizeof(double));
            psio_->write(unit_,"(Q|mn) Integrals", (char*)Qmnp[Q],mn_col_val*sizeof(double),addr,&addr);
        }
    }

    // ==> Close out <== //
    Qmn_.reset();
    delete[] eri;
}
void DFJK::unpack_Qmn(double **Qmnp, int naux)
{
     double* temp = new double[naux];

     double** Qmn2p = Qmn2_->pointer();
     int nbf = primary_->nbf();
     const std::vector<std::pair<int,int> >& pairs = sieve_->function_pairs();
     ULI nbf2 = primary_->nbf() * (ULI) primary_->nbf();
     ULI ntri = pairs.size(); 
     ::memset((void*) Qmn2p[0], '\0', sizeof(double)*naux*nbf2);
     for (ULI mn = 0; mn < pairs.size(); mn++) {

         int m = pairs[mn].first;
         int n = pairs[mn].second;

         if (m != n) {
             C_DCOPY(naux, &Qmnp[0][mn], ntri,  &Qmn2p[0][m*(ULI)nbf + n], nbf2);
             C_DCOPY(naux, &Qmnp[0][mn], ntri,  &Qmn2p[0][n*(ULI)nbf + m], nbf2);
         } else {
             C_DCOPY(naux, &Qmnp[0][mn], ntri,  &Qmn2p[0][m*(ULI)nbf + n], nbf2);
         }
         
    }
    
    delete[] temp;
}
void DFJK::manage_JK_core()
{
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);

        if (algorithm_ == 1) unpack_Qmn(&Qmn_->pointer()[Q],naux);

        block_J(&Qmn_->pointer()[Q],naux);
        block_K(&Qmn_->pointer()[Q],naux);
    } 
}
void DFJK::manage_JK_disk()
{
    // TODO AIO 
    int ntri = sieve_->function_pairs().size();
    Qmn_ = SharedMatrix(new Matrix("(Q|mn) Block", max_rows_, ntri));
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);
        psio_address addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri) * sizeof(double));
        psio_->read(unit_,"(Q|mn) Integrals", (char*)(Qmn_->pointer()[0]),sizeof(double)*naux*ntri,addr,&addr);

        if (algorithm_ == 1) unpack_Qmn(&Qmn_->pointer()[0],naux);

        block_J(&Qmn_->pointer()[0],naux);
        block_K(&Qmn_->pointer()[0],naux);
    } 
}
void DFJK::block_J(double** Qmnp, int naux)
{
    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    unsigned long int num_nm = function_pairs.size();

    for (int N = 0; N < J_ao_.size(); N++) {
        
        double** Dp   = D_ao_[N]->pointer();
        double** Jp   = J_ao_[N]->pointer();
        double*  J2p  = J_temp_->pointer();
        double*  D2p  = D_temp_->pointer();
        double*  dp   = d_temp_->pointer();

        for (unsigned long int mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            D2p[mn] = (m == n ? Dp[m][n] : Dp[m][n] + Dp[n][m]);
        }

        C_DGEMV('N',naux,num_nm,1.0,Qmnp[0],num_nm,D2p,1,0.0,dp,1);
        C_DGEMV('T',naux,num_nm,1.0,Qmnp[0],num_nm,dp,1,0.0,J2p,1);
        
        for (unsigned long int mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            Jp[m][n] += J2p[mn];
            Jp[n][m] += (m == n ? 0.0 : J2p[mn]);
        }
    }
}
void DFJK::block_K(double** Qmnp, int naux)
{
    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    const std::vector<long int>& function_pairs_reverse = sieve_->function_pairs_reverse();
    unsigned long int num_nm = function_pairs.size();

    if (algorithm_ == 0) {

        for (int N = 0; N < K_ao_.size(); N++) {

            int nbf = C_left_ao_[N]->rowspi()[0];
            int nocc = C_left_ao_[N]->colspi()[0];
            
            double** Clp  = C_left_ao_[N]->pointer();
            double** Crp  = C_right_ao_[N]->pointer();
            double** Elp  = E_left_->pointer();
            double** Erp  = E_right_->pointer();
            double** Kp   = K_ao_[N]->pointer();

            if (N == 0 || C_left_[N].get() != C_left_[N-1].get()) {
                
                #pragma omp parallel for schedule (dynamic)
                for (int m = 0; m < nbf; m++) {

                    int thread = 0;
                    #ifdef _OPENMP
                        thread = omp_get_thread_num();
                    #endif

                    double** Ctp = C_temp_[thread]->pointer();
                    double** QSp = Q_temp_[thread]->pointer();

                    const std::vector<int>& pairs = sieve_->function_to_function()[m];
                    int rows = pairs.size();

                    for (int i = 0; i < rows; i++) {
                        int n = pairs[i];
                        long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                        C_DCOPY(naux,&Qmnp[0][ij],num_nm,&QSp[0][i],nbf);
                        C_DCOPY(nocc,Clp[n],1,&Ctp[0][i],nbf);
                    }
                    
                    C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Elp[0][m*(ULI)nocc*naux],naux);
                }

            }

            if (!lr_symmetric_ && (N == 0 || C_right_[N].get() != C_right_[N-1].get())) {
                
                #pragma omp parallel for schedule (dynamic)
                for (int m = 0; m < nbf; m++) {

                    int thread = 0;
                    #ifdef _OPENMP
                        thread = omp_get_thread_num();
                    #endif

                    double** Ctp = C_temp_[thread]->pointer();
                    double** QSp = Q_temp_[thread]->pointer();

                    const std::vector<int>& pairs = sieve_->function_to_function()[m];
                    int rows = pairs.size();

                    for (int i = 0; i < rows; i++) {
                        int n = pairs[i];
                        long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                        C_DCOPY(naux,&Qmnp[0][ij],num_nm,&QSp[0][i],nbf);
                        C_DCOPY(nocc,Crp[n],1,&Ctp[0][i],nbf);
                    }
                    
                    C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Erp[0][m*(ULI)nocc*naux],naux);
                }
            }

            C_DGEMM('N','T',nbf,nbf,naux*nocc,1.0,Elp[0],naux*nocc,Erp[0],naux*nocc,1.0,Kp[0],nbf);        
        } 

    } else {

        for (int N = 0; N < K_ao_.size(); N++) {

            int nbf = C_left_ao_[N]->rowspi()[0];
            int nbf2 = nbf * (ULI) nbf; 
            int nocc = C_left_ao_[N]->colspi()[0];
            
            double*  sortp  = sort_->pointer();
            double** Qmn2p  = Qmn2_->pointer();
            double** Clp  = C_left_ao_[N]->pointer();
            double** Crp  = C_right_ao_[N]->pointer();
            double** Elp  = E_left_->pointer();
            double** Erp  = E_right_->pointer();
            double** Kp   = K_ao_[N]->pointer();

            if (N == 0 || C_left_[N].get() != C_left_[N-1].get()) {

                C_DGEMM('N','N', naux*(ULI) nbf, nocc, nbf, 1.0, Qmn2p[0], nbf, Clp[0], nocc, 0.0, Elp[0], nocc);

                for (int Q = 0; Q < naux; Q++) {    
                    ::memcpy((void*) sortp, (void*) &Elp[0][Q*(ULI)nocc*nbf], sizeof(double) * nocc * nbf);                   
                    for (int m = 0; m < nbf; m++) {
                        C_DCOPY(nocc,&sortp[m*nocc],1,&Elp[0][Q*(ULI)nocc*nbf + m], nbf);
                    }                    
                } 

            }

            if (!lr_symmetric_ && (N == 0 || C_right_[N].get() != C_right_[N-1].get())) {

                C_DGEMM('N','N', naux*(ULI) nbf, nocc, nbf, 1.0, Qmn2p[0], nbf, Crp[0], nocc, 0.0, Erp[0], nocc);

                for (int Q = 0; Q < naux; Q++) {    
                    ::memcpy((void*) sortp, (void*) &Erp[0][Q*(ULI)nocc*nbf], sizeof(double) * nocc * nbf);                   
                    for (int m = 0; m < nbf; m++) {
                        C_DCOPY(nocc,&sortp[m*nocc],1,&Erp[0][Q*(ULI)nocc*nbf + m], nbf);
                    }                    
                } 
    
            }

            C_DGEMM('T','N',nbf,nbf,naux*nocc,1.0,Elp[0],nbf,Erp[0],nbf,1.0,Kp[0],nbf);        

        }
    
    }
}

GPUDFJK::GPUDFJK(std::vector<SharedMatrix >& C_left, 
   std::vector<SharedMatrix >& C_right,
   boost::shared_ptr<BasisSet> primary, 
   boost::shared_ptr<BasisSet> auxiliary) :
   DFJK(C_left, C_right, primary, auxiliary)
{
    common_init();
}
GPUDFJK::GPUDFJK(std::vector<SharedMatrix >& C_symm,
   boost::shared_ptr<BasisSet> primary,
   boost::shared_ptr<BasisSet> auxiliary) :
   DFJK(C_symm, primary, auxiliary)
{
    common_init();
}
GPUDFJK::~GPUDFJK()
{
}
void GPUDFJK::common_init()
{
}
void GPUDFJK::preiterations()
{
    DFJK::preiterations();
    helper_ = boost::shared_ptr<GPUDFJKHelper>(new GPUDFJKHelper);
    helper_->Initialize(max_rows(),max_nocc(),primary_->nbf());
}
void GPUDFJK::postiterations()
{
    DFJK::postiterations();
    helper_->Finalize();
}

void GPUDFJK::print_header() const
{
    if (print_) {
    	fprintf(outfile, "  ==> GPUDFJK: Density-Fitted J/K Matrices <==\n");
	fprintf(outfile, "       by Eugene DePrince and Rob Parrish\n\n");

        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        fprintf(outfile, "    OpenMP threads:    %11d\n", omp_nthread_);
        fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        fprintf(outfile, "    Algorithm:         %11s\n",  (is_core_ ? "Core" : "Disk")); 
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        fprintf(outfile, "    Fitting Condition: %11.0E\n\n", condition_);

        fprintf(outfile, "   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level(outfile, print_);
    }
}
void GPUDFJK::block_J(double** Qmnp, int naux)
{
    // TODO Replace all this

    // J_mn = (mn|Q) (Q|ls) D_ls
    // (Q|ls) D_ls -> d_Q 
    // (Q|mn) d_Q  -> J_mn
    //
    // mn is stored in sparse packed triangular form
    // This can be used throughout this term

    // (Q|mn) is naux x num_nm
    // sieve tells you where relevent stuff is in nm
    // J is symmetric, D need not be

    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    unsigned long int num_nm = function_pairs.size();

    for (int N = 0; N < J_ao_.size(); N++) {
        
        double** Dp   = D_ao_[N]->pointer();
        double** Jp   = J_ao_[N]->pointer();
        double*  J2p  = J_temp_->pointer();
        double*  D2p  = D_temp_->pointer();
        double*  dp   = d_temp_->pointer();

        for (unsigned long int mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            D2p[mn] = (m == n ? Dp[m][n] : Dp[m][n] + Dp[n][m]);
        }

        C_DGEMV('N',naux,num_nm,1.0,Qmnp[0],num_nm,D2p,1,0.0,dp,1);
        C_DGEMV('T',naux,num_nm,1.0,Qmnp[0],num_nm,dp,1,0.0,J2p,1);
        
        for (unsigned long int mn = 0; mn < num_nm; ++mn) {
            int m = function_pairs[mn].first;
            int n = function_pairs[mn].second;
            Jp[m][n] += J2p[mn];
            Jp[n][m] += (m == n ? 0.0 : J2p[mn]);
        }
    }
}
void GPUDFJK::block_K(double** Qmnp, int naux)
{
    // TODO Replace all this

    // K_mn = C_left_il (ml|Q) (Q|sn) C_right_is
    // (Q|lm) C_left_il -> E_left_miQ
    // (Q|lm) C_right_il -> E_right_miQ
    // E_left_miQ E_right_niQ -> K_mn
    //
    // K is not symmetric
    // The (Q|mn) tensor must be unpacked to do the
    // first GEMM
    // Sparsity can still be used in the first GEMM

    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    const std::vector<long int>& function_pairs_reverse = sieve_->function_pairs_reverse();
    unsigned long int num_nm = function_pairs.size();

    for (int N = 0; N < K_ao_.size(); N++) {

        int nbf = C_left_ao_[N]->rowspi()[0];
        int nocc = C_left_ao_[N]->colspi()[0];
        
        double** Clp  = C_left_ao_[N]->pointer();
        double** Crp  = C_right_ao_[N]->pointer();
        double** Elp  = E_left_->pointer();
        double** Erp  = E_right_->pointer();
        double** Kp   = K_ao_[N]->pointer();

        if (N == 0 || C_left_ao_[N].get() != C_left_ao_[N-1].get()) {
            
            #pragma omp parallel for schedule (dynamic)
            for (int m = 0; m < nbf; m++) {

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                double** Ctp = C_temp_[thread]->pointer();
                double** QSp = Q_temp_[thread]->pointer();

                const std::vector<int>& pairs = sieve_->function_to_function()[m];
                int rows = pairs.size();

                for (int i = 0; i < rows; i++) {
                    int n = pairs[i];
                    long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                    C_DCOPY(naux,&Qmnp[0][ij],num_nm,&QSp[0][i],nbf);
                    C_DCOPY(nocc,Clp[n],1,&Ctp[0][i],nbf);
                }
                
                helper_->GPU_DGEMM_2DTile('T','N',naux,nocc,rows,1.0,QSp[0],nbf,Ctp[0],nbf,0.0,&Elp[0][m*(ULI)nocc*naux],naux,thread);

            }

        }

        if (!lr_symmetric_ && (N == 0 || C_right_[N].get() != C_right_[N-1].get())) {
            
            #pragma omp parallel for schedule (dynamic)
            for (int m = 0; m < nbf; m++) {

                int thread = 0;
                #ifdef _OPENMP
                    thread = omp_get_thread_num();
                #endif

                double** Ctp = C_temp_[thread]->pointer();
                double** QSp = Q_temp_[thread]->pointer();

                const std::vector<int>& pairs = sieve_->function_to_function()[m];
                int rows = pairs.size();

                for (int i = 0; i < rows; i++) {
                    int n = pairs[i];
                    long int ij = function_pairs_reverse[(m >= n ? (m * (m + 1L) >> 1) + n : (n * (n + 1L) >> 1) + m)];
                    C_DCOPY(naux,&Qmnp[0][ij],num_nm,&QSp[0][i],nbf);
                    C_DCOPY(nocc,Crp[n],1,&Ctp[0][i],nbf);
                }
                
                helper_->GPU_DGEMM_2DTile('T','N',naux,nocc,rows,1.0,QSp[0],nbf,Ctp[0],nbf,0.0,&Erp[0][m*(ULI)nocc*naux],naux,thread);

            }
        }

        helper_->GPU_DGEMM_2DTile('T','N',nbf,nbf,naux*nocc,1.0,Erp[0],naux*nocc,Elp[0],naux*nocc,1.0,Kp[0],nbf,0);

    } 
}

}
