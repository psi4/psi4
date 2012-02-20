#include <libmints/mints.h>
#include <lib3index/3index.h>
#include <libpsio/psio.hpp>
#include <libpsio/aiohandler.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <psifiles.h>
#include <libmints/sieve.h>
#include <libiwl/iwl.hpp>
#include "jk.h"
#include "cubature.h"
#include "points.h"

#include <sstream>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PRINT_ME fprintf(stdout,"If in %d\n",__LINE__);
#define PRINT_ME

using namespace std;
using namespace psi;

namespace psi {

JK::JK( boost::shared_ptr<BasisSet> primary) :
    primary_(primary)
{
    common_init();
}
JK::~JK()
{
}
boost::shared_ptr<JK> JK::build_JK()
{
    Options& options = Process::environment.options;
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> primary = BasisSet::construct(parser, Process::environment.molecule(), "BASIS");

    if (options.get_str("SCF_TYPE") == "DF") {

        if (options.get_str("DF_BASIS_SCF") == "") {
            primary->molecule()->set_basis_all_atoms(options.get_str("BASIS") + "-JKFIT", "DF_BASIS_SCF");
        }
        boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, primary->molecule(), "DF_BASIS_SCF");

        DFJK* jk = new DFJK(primary,auxiliary);

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        if (options["DF_INTS_IO"].has_changed())
            jk->set_df_ints_io(options.get_str("DF_INTS_IO"));
        if (options["DF_FITTING_CONDITION"].has_changed())
            jk->set_condition(options.get_double("DF_FITTING_CONDITION"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));

        return boost::shared_ptr<JK>(jk);

    } else if (options.get_str("SCF_TYPE") == "PK") {

        DiskJK* jk = new DiskJK(primary);;

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));

        return boost::shared_ptr<JK>(jk);

    } else if (options.get_str("SCF_TYPE") == "DIRECT") {

        DirectJK* jk = new DirectJK(primary);

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));

        return boost::shared_ptr<JK>(jk);

    } else if (options.get_str("SCF_TYPE") == "PK") {
    
        PSJK* jk = new PSJK(primary,options);

        if (options["INTS_TOLERANCE"].has_changed())
            jk->set_cutoff(options.get_double("INTS_TOLERANCE"));
        if (options["PRINT"].has_changed())
            jk->set_print(options.get_int("PRINT"));
        if (options["DEBUG"].has_changed())
            jk->set_debug(options.get_int("DEBUG"));
        if (options["BENCH"].has_changed())
            jk->set_bench(options.get_int("BENCH"));
        if (options["DF_INTS_NUM_THREADS"].has_changed())
            jk->set_df_ints_num_threads(options.get_int("DF_INTS_NUM_THREADS"));
        if (options["PS_THETA"].has_changed())
            jk->set_theta(options.get_double("PS_THETA"));
        if (options["PS_DEALIASING"].has_changed())
            jk->set_dealiasing(options.get_str("PS_DEALIASING"));
        if (options["DEALIAS_BASIS_SCF"].has_changed()) {
            if (options.get_str("DEALIAS_BASIS_SCF") != "") {
                boost::shared_ptr<BasisSet> dealias = BasisSet::construct(parser, primary->molecule(), "DEALIAS_BASIS_SCF");
                jk->set_dealias_basis(dealias);
            }
        }

        return boost::shared_ptr<JK>(jk);

    } else {
        throw PSIEXCEPTION("JK::build_JK: Unknown SCF Type");
    }
}
void JK::common_init()
{
    print_ = 1;
    debug_ = 0;
    bench_ = 0;

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
    bool same = true;
    if (C_left_.size() != D_.size()) {
        same = false;
    } else {
        for (int N = 0; N < D_.size(); N++) {
            if (D_[N]->symmetry() != (C_left_[N]->symmetry() ^ C_right_[N]->symmetry()))
                same = false;
        }
    }

    if (!same) {
        D_.clear();
        for (int N = 0; N < C_left_.size(); ++N) {
            std::stringstream s;
            s << "D " << N << " (SO)";
            D_.push_back(SharedMatrix(new Matrix(s.str(),C_left_[N]->nirrep(), C_left_[N]->rowspi(), C_right_[N]->rowspi(), C_left_[N]->symmetry() ^ C_right_[N]->symmetry())));
        }
    }

    for (int N = 0; N < D_.size(); ++N) {
        int symm = D_[N]->symmetry();
        D_[N]->zero();
        for (int h = 0; h < D_[N]->nirrep(); ++h) {

            int nsol = C_left_[N]->rowspi()[h];
            int nocc = C_left_[N]->colspi()[h];
            int nsor = C_right_[N]->rowspi()[h^symm];

            if (!nsol || !nsor || !nocc) continue;

            double** Dp = D_[N]->pointer(h);
            double** Clp = C_left_[N]->pointer(h);
            double** Crp = C_right_[N]->pointer(h^symm);

            C_DGEMM('N','T', nsol, nsor, nocc, 1.0, Clp[0], nocc, Crp[0], nocc, 0.0, Dp[0], nsor);
        }
    }
}
void JK::allocate_JK()
{
    // Allocate J/K in the case that the algorithm uses USOs, so AO2USO will not allocate. 
    bool same = true;
    if (J_.size() != D_.size()) {
        same = false;
    } else {
        for (int N = 0; N < D_.size(); N++) {
            if (D_[N]->symmetry() != J_[N]->symmetry())
                same = false;
        }
    }

    if (!same) {
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
    allocate_JK();

    // If C1, C_ao and D_ao are equal to C and D
    if (AO2USO_->nirrep() == 1) {
        C_left_ao_ = C_left_;
        C_right_ao_ = C_right_;
        D_ao_ = D_;
        J_ao_ = J_;
        K_ao_ = K_;
        wK_ao_ = wK_; 
        return;
    }

    if (J_ao_.size() != D_.size()) {
        J_ao_.clear();    
        K_ao_.clear();    
        wK_ao_.clear();    
        D_ao_.clear();
        int nao = AO2USO_->rowspi()[0];
        for (int N = 0; N < D_.size() && do_J_; ++N) {
            std::stringstream s;
            s << "J " << N << " (AO)";
            J_ao_.push_back(SharedMatrix(new Matrix(s.str(),nao,nao)));
        }
        for (int N = 0; N < D_.size() && do_K_; ++N) {
            std::stringstream s;
            s << "K " << N << " (AO)";
            K_ao_.push_back(SharedMatrix(new Matrix(s.str(),nao,nao)));
        }
        for (int N = 0; N < D_.size() && do_wK_; ++N) {
            std::stringstream s;
            s << "wK " << N << " (AO)";
            wK_ao_.push_back(SharedMatrix(new Matrix(s.str(),nao,nao)));
        }
        for (int N = 0; N < D_.size(); ++N) {
            std::stringstream s;
            s << "D " << N << " (AO)";
            D_ao_.push_back(SharedMatrix(new Matrix(s.str(),nao,nao)));
        }
    }

    // Always reallocate C matrices, the occupations are tricky
    C_left_ao_.clear();
    C_right_ao_.clear();
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

    // Alias pointers if lr_symmetric_
    if (lr_symmetric_) {
        C_right_ao_ = C_left_ao_;
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

    for (int N = 0; N < D_.size(); ++N) {
        if (do_J_) J_ao_[N]->zero();
        if (do_K_) K_ao_[N]->zero();
        if (do_wK_) wK_ao_[N]->zero();
    }
}
void JK::AO2USO()
{
    // If already C1, J/K are J_ao/K_ao, pointers are already aliased
    if (AO2USO_->nirrep() == 1) {
        return;
    }

    // If not C1, J/K/wK are already allocated

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
    if (C_left_.size() && !C_right_.size()) {
        lr_symmetric_ = true;
        C_right_ = C_left_;
    } else {
        lr_symmetric_ = false;
    }

    timer_on("JK: D");
    compute_D();
    timer_off("JK: D");
    if (C1()) {
        timer_on("JK: USO2AO");
        USO2AO();
        timer_off("JK: USO2AO");
    } else {
        allocate_JK();
    }

    timer_on("JK: JK");
    compute_JK();
    timer_off("JK: JK");

    if (C1()) {
        timer_on("JK: AO2USO");
        AO2USO();
        timer_off("JK: AO2USO");
    }

    if (debug_ > 3) {
        fprintf(outfile, "   > JK <\n\n");
        for (int N = 0; N < C_left_.size(); N++) {
            if (C1()) {
                C_left_ao_[N]->print(outfile);
                C_right_ao_[N]->print(outfile);
                D_ao_[N]->print(outfile);
                J_ao_[N]->print(outfile);
                K_ao_[N]->print(outfile);
            }
            C_left_[N]->print(outfile);
            C_right_[N]->print(outfile);
            D_[N]->print(outfile);
            J_[N]->print(outfile);
            K_[N]->print(outfile);
        }
        fflush(outfile);
    }

    if (lr_symmetric_) {
        C_right_.clear();
    }
}
void JK::finalize()
{
    postiterations();
}

DiskJK::DiskJK(boost::shared_ptr<BasisSet> primary) :
   JK(primary)
{
    common_init();
}
DiskJK::~DiskJK()
{
}
void DiskJK::common_init()
{
}
void DiskJK::print_header() const
{
    if (print_) {
        fprintf(outfile, "  ==> DiskJK: Disk-Based J/K Matrices <==\n\n");

        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n\n", cutoff_);
    }
}
void DiskJK::preiterations()
{
    boost::shared_ptr<MintsHelper> mints(new MintsHelper());
    mints->integrals();

    boost::shared_ptr<SOBasisSet> bas = mints->sobasisset();

    so2symblk_ = new int[primary_->nbf()];
    so2index_  = new int[primary_->nbf()];
    size_t so_count = 0;
    size_t offset = 0;
    for (int h = 0; h < bas->nirrep(); ++h) {
        for (int i = 0; i < bas->dimension()[h]; ++i) {
            so2symblk_[so_count] = h;
            so2index_[so_count] = so_count-offset;
            ++so_count;
        }
        offset += bas->dimension()[h];
    }
}
void DiskJK::compute_JK()
{
    boost::shared_ptr<PSIO> psio(new PSIO());
    IWL *iwl = new IWL(psio.get(), PSIF_SO_TEI, cutoff_, 1, 1);
    Label *lblptr = iwl->labels();
    Value *valptr = iwl->values();
    int labelIndex, pabs, qabs, rabs, sabs, prel, qrel, rrel, srel, psym, qsym, rsym, ssym;
    double value;
    bool lastBuffer;
    do{
        lastBuffer = iwl->last_buffer();
        for(int index = 0; index < iwl->buffer_count(); ++index){
            labelIndex = 4*index;
            pabs  = abs((int) lblptr[labelIndex++]);
            qabs  = (int) lblptr[labelIndex++];
            rabs  = (int) lblptr[labelIndex++];
            sabs  = (int) lblptr[labelIndex++];
            prel  = so2index_[pabs];
            qrel  = so2index_[qabs];
            rrel  = so2index_[rabs];
            srel  = so2index_[sabs];
            psym  = so2symblk_[pabs];
            qsym  = so2symblk_[qabs];
            rsym  = so2symblk_[rabs];
            ssym  = so2symblk_[sabs];
            value = (double) valptr[index];

            int pqsym = psym ^ qsym;
            int rssym = rsym ^ ssym;
            int qrsym = qsym ^ rsym;
            int pssym = psym ^ ssym;
            int prsym = psym ^ rsym;
            int qssym = qsym ^ ssym;


            for (int N = 0; N < J_.size(); N++) {

                SharedMatrix J = J_[N];
                SharedMatrix K = K_[N];
                SharedMatrix D = D_[N];

                int sym = J->symmetry();

                /* (pq|rs) */
                if(pqsym == rssym && pqsym == sym){
                    J->add(rsym, rrel, srel, D->get(psym, prel, qrel) * value);
                    PRINT_ME
                }

                if(qrsym == pssym && qrsym == sym){
                    K->add(qsym, qrel, rrel, D->get(psym, prel, srel) * value);
                    PRINT_ME
                }

                if(pabs!=qabs && rabs!=sabs && (pabs!=rabs || qabs!=sabs)){
                    /* (pq|sr) */
                    if(qssym == prsym && qssym == sym){
                        K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        PRINT_ME
                    }

                    /* (qp|rs) */
                    if(rssym == pqsym && rssym == sym){
                        J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        PRINT_ME
                    }

                    if(prsym == qssym && prsym == sym){
                        K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        PRINT_ME
                    }

                    /* (qp|sr) */
                    if(pssym == qrsym && pssym == sym){
                        K->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        PRINT_ME
                    }

                    /* (rs|pq) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        PRINT_ME
                    }

                    if(pssym == qrsym && pssym == sym){
                        K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        PRINT_ME
                    }

                    /* (sr|pq) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(psym, prel, qrel, D->get(ssym, srel, rrel) * value);
                        PRINT_ME
                    }

                    if(prsym == qssym && prsym == sym){
                        K->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        PRINT_ME
                    }

                    /* (rs|qp) */
                    if(qssym == prsym && qssym == sym){
                        K->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        PRINT_ME
                    }

                    /* (sr|qp) */
                    if(qrsym == pssym && qrsym == sym){
                        K->add(rsym, rrel, qrel, D->get(ssym, srel, prel) * value);
                        PRINT_ME
                    }
                }else if(pabs!=qabs && rabs!=sabs && pabs==rabs && qabs==sabs){
                    /* (pq|sr) */
                    if(qssym == prsym && qssym == sym){
                        K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        PRINT_ME
                    }
                    /* (qp|rs) */
                    if(rssym == pqsym && rssym == sym){
                        J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        PRINT_ME
                    }
                    if(prsym == qssym && prsym == sym){
                        K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        PRINT_ME
                    }

                    /* (qp|sr) */
                    if(pssym == qrsym && pssym == sym){
                        K->add(psym, prel, srel, D->get(qsym, qrel, rrel) * value);
                        PRINT_ME
                    }
                }else if(pabs!=qabs && rabs==sabs){
                    /* (qp|rs) */
                    if(rssym == pqsym && rssym == sym){
                        J->add(rsym, rrel, srel, D->get(qsym, qrel, prel) * value);
                        PRINT_ME
                    }

                    if(prsym == qssym && prsym == sym){
                        K->add(psym, prel, rrel, D->get(qsym, qrel, srel) * value);
                        PRINT_ME
                    }

                    /* (rs|pq) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        PRINT_ME
                    }

                    if(pssym == qrsym && pssym == sym){
                        K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        PRINT_ME
                    }

                    /* (rs|qp) */
                    if(qssym == prsym && qssym == sym){
                        K->add(ssym, srel, qrel, D->get(rsym, rrel, prel) * value);
                        PRINT_ME
                    }
                }else if(pabs==qabs && rabs!=sabs){
                    /* (pq|sr) */
                    if(qssym == prsym && qssym == sym){
                        K->add(qsym, qrel, srel, D->get(psym, prel, rrel) * value);
                        PRINT_ME
                    }

                    /* (rs|pq) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        PRINT_ME
                    }

                    if(pssym == qrsym && pssym == sym){
                        K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        PRINT_ME
                    }

                    /* (sr|pq) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(psym, prel, qrel, D->get(ssym, srel, rrel) * value);
                        PRINT_ME
                    }

                    if(prsym == qssym && prsym == sym){
                        K->add(rsym, rrel, prel, D->get(ssym, srel, qrel) * value);
                        PRINT_ME
                    }
                }else if(pabs==qabs && rabs==sabs && (pabs!=rabs || qabs!=sabs)){
                    /* (rs|pq) */
                    if(pqsym == rssym && pqsym == sym){
                        J->add(psym, prel, qrel, D->get(rsym, rrel, srel) * value);
                        PRINT_ME
                    }

                    if(pssym == qrsym && pssym == sym){
                        K->add(ssym, srel, prel, D->get(rsym, rrel, qrel) * value);
                        PRINT_ME
                    }
                }
            }
        } /* end loop through current buffer */
        if(!lastBuffer) iwl->fetch();
    }while(!lastBuffer);
    iwl->set_keep_flag(1);

    for (int N = 0; N < J_.size(); N++) {
        J_[N]->copy_lower_to_upper();
        if (K_[N]->symmetry()) K_[N]->transpose_this();
    }

}
void DiskJK::postiterations()
{
    delete[] so2symblk_;
    delete[] so2index_;
}

DirectJK::DirectJK(boost::shared_ptr<BasisSet> primary) :
   JK(primary)
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
        fprintf(outfile, "  ==> DirectJK: Integral-Direct J/K Matrices <==\n\n");

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

        int nM = primary_->shell(M).nfunction();
        int nN = primary_->shell(N).nfunction();
        int nR = primary_->shell(R).nfunction();
        int nS = primary_->shell(S).nfunction();

        int sM = primary_->shell(M).function_index();
        int sN = primary_->shell(N).function_index();
        int sR = primary_->shell(R).function_index();
        int sS = primary_->shell(S).function_index();

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

DFJK::DFJK(boost::shared_ptr<BasisSet> primary,
   boost::shared_ptr<BasisSet> auxiliary) :
   JK(primary), auxiliary_(auxiliary)
{
    common_init();
}
DFJK::~DFJK()
{
}
void DFJK::common_init()
{
    df_ints_num_threads_ = 1;
    #ifdef _OPENMP
        df_ints_num_threads_ = omp_get_max_threads();
    #endif
    df_ints_io_ = "NONE";
    condition_ = 1.0E-12;
    unit_ = PSIF_DFSCF_BJ;
    is_core_ = true;
    psio_ = PSIO::shared_object();
}
void DFJK::print_header() const
{
    if (print_) {
        fprintf(outfile, "  ==> DFJK: Density-Fitted J/K Matrices <==\n\n");

        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        fprintf(outfile, "    OpenMP threads:    %11d\n", omp_nthread_);
        fprintf(outfile, "    Integrals threads: %11d\n", df_ints_num_threads_);
        fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        fprintf(outfile, "    Algorithm:         %11s\n",  (is_core_ ? "Core" : "Disk"));
        fprintf(outfile, "    Integral Cache:    %11s\n",  df_ints_io_.c_str());
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        fprintf(outfile, "    Fitting Condition: %11.0E\n\n", condition_);

        fprintf(outfile, "   => Auxiliary Basis Set <=\n\n");
        auxiliary_->print_by_level(outfile, print_);
    }
}
bool DFJK::is_core() const
{
    int ntri = sieve_->function_pairs().size();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    // Two is for buffer space in fitting
    if (do_wK_)
        return (3L*three_memory + 2L*two_memory < memory_);
    else
        return (three_memory + 2L*two_memory < memory_);
}
unsigned long int DFJK::memory_temp()
{
    unsigned long int mem = 0L;

    // J Overhead (Jtri, Dtri, d)
    mem += 2L * sieve_->function_pairs().size() + auxiliary_->nbf();
    // K Overhead (C_temp, Q_temp)
    mem += omp_nthread_ * (unsigned long int) primary_->nbf() * (auxiliary_->nbf() + max_nocc());

    return mem;
}
int DFJK::max_rows()
{
    // Start with all memory
    unsigned long int mem = memory_;
    // Subtract J/K/wK/C/D overhead
    mem -= memory_overhead();
    // Subtract threading temp overhead
    mem -= memory_temp();

    // How much will each row cost?
    unsigned long int row_cost = 0L;
    // Copies of E tensor
    row_cost += (lr_symmetric_ ? 1L : 2L) * max_nocc() * primary_->nbf();
    // Slices of Qmn tensor, including AIO buffer (NOTE: AIO not implemented yet)
    row_cost += (is_core_ ? 1L : 1L) * sieve_->function_pairs().size();

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


    #ifdef _OPENMP
    int temp_nthread = omp_get_max_threads();
    omp_set_num_threads(omp_nthread_);
        C_temp_.resize(omp_nthread_);
        Q_temp_.resize(omp_nthread_);
        #pragma omp parallel
        {
            C_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf()));
            Q_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Qtemp", max_rows_, primary_->nbf()));
        }
    omp_set_num_threads(temp_nthread);
    #else
        for (int thread = 0; thread < omp_nthread_; thread++) {
            C_temp_.push_back(SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf())));
            Q_temp_.push_back(SharedMatrix(new Matrix("Qtemp", max_rows_, primary_->nbf())));
        }
    #endif

    E_left_ = SharedMatrix(new Matrix("E_left", primary_->nbf(), max_rows_ * max_nocc_));
    if (lr_symmetric_)
        E_right_ = E_left_;
    else
        E_right_ = boost::shared_ptr<Matrix>(new Matrix("E_right", primary_->nbf(), max_rows_ * max_nocc_));

}
void DFJK::initialize_w_temps()
{
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);

    #ifdef _OPENMP
    int temp_nthread = omp_get_max_threads();
    omp_set_num_threads(omp_nthread_);
        C_temp_.resize(omp_nthread_);
        Q_temp_.resize(omp_nthread_);
        #pragma omp parallel
        {
            C_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf()));
            Q_temp_[omp_get_thread_num()] = SharedMatrix(new Matrix("Qtemp", max_rows_w, primary_->nbf()));
        }
    omp_set_num_threads(temp_nthread);
    #else
        for (int thread = 0; thread < omp_nthread_; thread++) {
            C_temp_.push_back(SharedMatrix(new Matrix("Ctemp", max_nocc_, primary_->nbf())));
            Q_temp_.push_back(SharedMatrix(new Matrix("Qtemp", max_rows_w, primary_->nbf())));
        }
    #endif

    E_left_  = SharedMatrix(new Matrix("E_left", primary_->nbf(), max_rows_w * max_nocc_));
    E_right_ = SharedMatrix(new Matrix("E_right", primary_->nbf(), max_rows_w * max_nocc_));

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
}
void DFJK::free_w_temps()
{
    E_left_.reset();
    E_right_.reset();
    C_temp_.clear();
    Q_temp_.clear();
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

    if (do_wK_) {
        if (is_core_)
            initialize_wK_core();
        else
            initialize_wK_disk();
    }
}
void DFJK::compute_JK()
{
    max_nocc_ = max_nocc();
    max_rows_ = max_rows();

    if (do_J_ || do_K_) {
        initialize_temps();
        if (is_core_)
            manage_JK_core();
        else
            manage_JK_disk();
        free_temps();
    }

    if (do_wK_) {
        initialize_w_temps();
        if (is_core_)
            manage_wK_core();
        else
            manage_wK_disk();
        free_w_temps();
        // Bring the wK matrices back to Hermitian
        if (lr_symmetric_) {
            for (int N = 0; N < wK_ao_.size(); N++) {
                wK_ao_[N]->hermitivitize();
            }
        }
    }
}
void DFJK::postiterations()
{
    Qmn_.reset();
    Qlmn_.reset();
    Qrmn_.reset();
}
void DFJK::initialize_JK_core()
{
    int ntri = sieve_->function_pairs().size();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif
    int rank = 0;

    Qmn_ = SharedMatrix(new Matrix("Qmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qmn_->pointer();

    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_,1);
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

    timer_on("JK: (A|mn)");

    //The integrals (A|mn)
    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell).nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU).function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU).function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell).function_index() + P;
                                    Qmnp[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("JK: (A|mn)");

    delete []buffer;
    delete []eri;

    timer_on("JK: (A|Q)^-1/2");

    boost::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true));
    Jinv->form_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    timer_off("JK: (A|Q)^-1/2");

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

    timer_on("JK: (Q|mn)");

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

    timer_off("JK: (Q|mn)");
    //Qmn_->print();

    if (df_ints_io_ == "SAVE") {
        psio_->open(unit_,PSIO_OPEN_NEW);
        psio_->write_entry(unit_, "(Q|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_,1);
    }
}
void DFJK::initialize_JK_disk()
{
    // Try to load
    if (df_ints_io_ == "LOAD") {
        return;
    }

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

    timer_on("JK: (A|Q)^-1");

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

    timer_off("JK: (A|Q)^-1");

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
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

        timer_on("JK: (A|mn)");

        #pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {

            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int MU = schwarz_shell_pairs[MUNU + 0].first;
            int NU = schwarz_shell_pairs[MUNU + 0].second;
            int nummu = primary_->shell(MU).nfunction();
            int numnu = primary_->shell(NU).nfunction();
            int mu = primary_->shell(MU).function_index();
            int nu = primary_->shell(NU).function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P).nfunction();
                int p = auxiliary_->shell(P).function_index();
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

        timer_off("JK: (A|mn)");

        // ==> (Q|mn) fitting <== //

        timer_on("JK: (Q|mn)");

        for (int mn = 0; mn < mn_col_val; mn+=naux) {
            int cols = naux;
            if (mn + naux >= mn_col_val)
                cols = mn_col_val - mn;

            for (int Q = 0; Q < naux; Q++)
                C_DCOPY(cols,&Qmnp[Q][mn],1,Amnp[Q],1);

            C_DGEMM('N','N',naux,cols,naux,1.0,Jinvp[0],naux,Amnp[0],naux,0.0,&Qmnp[0][mn],max_cols);
        }

        timer_off("JK: (Q|mn)");

        // ==> Disk striping <== //

        timer_on("JK: (Q|mn) Write");

        psio_address addr;
        for (int Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri + mn_start_val)*sizeof(double));
            psio_->write(unit_,"(Q|mn) Integrals", (char*)Qmnp[Q],mn_col_val*sizeof(double),addr,&addr);
        }

        timer_off("JK: (Q|mn) Write");
    }

    // ==> Close out <== //
    Qmn_.reset();
    delete[] eri;

    psio_->close(unit_,1);
}
void DFJK::initialize_wK_core()
{
    int naux = auxiliary_->nbf();
    int ntri = sieve_->function_pairs().size();
    ULI three_memory = ((ULI)auxiliary_->nbf())*ntri;
    ULI two_memory = ((ULI)auxiliary_->nbf())*auxiliary_->nbf();

    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif
    int rank = 0;

    // Check that the right integrals are using the correct omega
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        double check_omega;
        psio_->read_entry(unit_, "Omega", (char*)&check_omega, sizeof(double));
        if (check_omega != omega_) {
            rebuild_wK_disk();
        }
        psio_->close(unit_,1);
    }

    Qlmn_ = SharedMatrix(new Matrix("Qlmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmnp = Qlmn_->pointer();

    Qrmn_ = SharedMatrix(new Matrix("Qrmn (Fitted Integrals)",
        auxiliary_->nbf(), ntri));
    double** Qmn2p = Qrmn_->pointer();

    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->read_entry(unit_, "Left (Q|w|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->read_entry(unit_, "Right (Q|w|mn) Integrals", (char*) Qmn2p[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->close(unit_,1);
        return;
    }

    // => Left Integrals <= //

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

    timer_on("JK: (A|mn)^L");

    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell).nfunction();
                    eri[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU).function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU).function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell).function_index() + P;
                                    Qmn2p[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("JK: (A|mn)^L");

    delete []buffer;
    delete []eri;

    // => Fitting <= //

    timer_on("JK: (A|Q)^-1");

    // Fitting metric
    boost::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true));
    Jinv->form_full_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    timer_off("JK: (A|Q)^-1");

    timer_on("JK: (Q|mn)^L");

    // Fitting in one GEMM (being a clever bastard)
    C_DGEMM('N','N',naux,ntri,naux,1.0,Jinvp[0],naux,Qmn2p[0],ntri,0.0,Qmnp[0],ntri);

    timer_off("JK: (Q|mn)^L");

    // => Right Integrals <= //

    const double **buffer2 = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri2 = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri2[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        buffer2[Q] = eri2[Q]->buffer();
    }

    //The integrals (A|w|mn)

    timer_on("JK: (A|mn)^R");

    #pragma omp parallel for private (numP, Pshell, MU, NU, P, PHI, mu, nu, nummu, numnu, omu, onu, rank) schedule (dynamic) num_threads(nthread)
    for (MU=0; MU < primary_->nshell(); ++MU) {
        #ifdef _OPENMP
            rank = omp_get_thread_num();
        #endif
        nummu = primary_->shell(MU).nfunction();
        for (NU=0; NU <= MU; ++NU) {
            numnu = primary_->shell(NU).nfunction();
            if (schwarz_shell_pairs[MU*(MU+1)/2+NU] > -1) {
                for (Pshell=0; Pshell < auxiliary_->nshell(); ++Pshell) {
                    numP = auxiliary_->shell(Pshell).nfunction();
                    eri2[rank]->compute_shell(Pshell, 0, MU, NU);
                    for (mu=0 ; mu < nummu; ++mu) {
                        omu = primary_->shell(MU).function_index() + mu;
                        for (nu=0; nu < numnu; ++nu) {
                            onu = primary_->shell(NU).function_index() + nu;
                            if(omu>=onu && schwarz_fun_pairs[omu*(omu+1)/2+onu] > -1) {
                                for (P=0; P < numP; ++P) {
                                    PHI = auxiliary_->shell(Pshell).function_index() + P;
                                    Qmn2p[PHI][schwarz_fun_pairs[omu*(omu+1)/2+onu]] = buffer2[rank][P*nummu*numnu + mu*numnu + nu];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    timer_off("JK: (A|mn)^R");

    delete []buffer2;
    delete []eri2;

    // Try to save
    if (df_ints_io_ == "SAVE") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        psio_->write_entry(unit_, "Left (Q|w|mn) Integrals", (char*) Qmnp[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->write_entry(unit_, "Right (Q|w|mn) Integrals", (char*) Qmn2p[0], sizeof(double) * ntri * auxiliary_->nbf());
        psio_->write_entry(unit_, "Omega", (char*) &omega_, sizeof(double));
        psio_->close(unit_,1);
    }
}
void DFJK::initialize_wK_disk()
{
    // Try to load
    if (df_ints_io_ == "LOAD") {
        psio_->open(unit_,PSIO_OPEN_OLD);
        double check_omega;
        psio_->read_entry(unit_, "Omega", (char*)&check_omega, sizeof(double));
        if (check_omega != omega_) {
            rebuild_wK_disk();
        }
        psio_->close(unit_,1);
    }

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
    ULI buffer_memory = memory_ - 2L*two_memory; // Two is for buffer space in fitting

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
    psio_->open(unit_,PSIO_OPEN_OLD);
    boost::shared_ptr<AIOHandler> aio(new AIOHandler(psio_));

    // Dispatch the prestripe
    aio->zero_disk(unit_,"Left (Q|w|mn) Integrals",naux,ntri);

    // Form the J full inverse
    boost::shared_ptr<FittingMetric> Jinv(new FittingMetric(auxiliary_, true));
    Jinv->form_full_eig_inverse();
    double** Jinvp = Jinv->get_metric()->pointer();

    // Synch up
    aio->synchronize();

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
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

        timer_on("JK: (A|mn)^L");

        #pragma omp parallel for schedule(guided) num_threads(nthread)
        for (int MUNU = MN_start_val; MUNU < MN_start_val + MN_col_val; MUNU++) {

            int rank = 0;
            #ifdef _OPENMP
                rank = omp_get_thread_num();
            #endif

            int MU = schwarz_shell_pairs[MUNU + 0].first;
            int NU = schwarz_shell_pairs[MUNU + 0].second;
            int nummu = primary_->shell(MU).nfunction();
            int numnu = primary_->shell(NU).nfunction();
            int mu = primary_->shell(MU).function_index();
            int nu = primary_->shell(NU).function_index();
            for (int P = 0; P < auxiliary_->nshell(); P++) {
                int nump = auxiliary_->shell(P).nfunction();
                int p = auxiliary_->shell(P).function_index();
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

        timer_off("JK: (A|mn)^L");

        // ==> (Q|mn) fitting <== //

        timer_on("JK: (Q|mn)^L");

        for (int mn = 0; mn < mn_col_val; mn+=naux) {
            int cols = naux;
            if (mn + naux >= mn_col_val)
                cols = mn_col_val - mn;

            for (int Q = 0; Q < naux; Q++)
                C_DCOPY(cols,&Qmnp[Q][mn],1,Amnp[Q],1);

            C_DGEMM('N','N',naux,cols,naux,1.0,Jinvp[0],naux,Amnp[0],naux,0.0,&Qmnp[0][mn],max_cols);
        }

        timer_off("JK: (Q|mn)^L");

        // ==> Disk striping <== //

        timer_on("JK: (Q|mn)^L Write");

        psio_address addr;
        for (int Q = 0; Q < naux; Q++) {
            addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri + mn_start_val)*sizeof(double));
            psio_->write(unit_,"Left (Q|w|mn) Integrals", (char*)Qmnp[Q],mn_col_val*sizeof(double),addr,&addr);
        }

        timer_off("JK: (Q|mn)^L Write");

    }

    Qmn_.reset();
    delete[] eri;

    // => Right Integrals <= //

    const double **buffer2 = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri2 = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri2[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        buffer2[Q] = eri2[Q]->buffer();
    }

    ULI maxP = auxiliary_->max_function_per_shell();
    ULI max_rows = memory_ / ntri;
    max_rows = (max_rows > naux ? naux : max_rows);
    max_rows = (max_rows < maxP ? maxP : max_rows);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        int nQ = auxiliary_->shell(Q).nfunction();
        if (counter + nQ > max_rows) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(auxiliary_->nshell());

    SharedMatrix Amn2(new Matrix("(A|mn) Block", max_rows, ntri));
    double** Amn2p = Amn2->pointer();
    psio_address next_AIA = PSIO_ZERO;

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = auxiliary_->shell(Qstart).function_index();
        int nrows  = (Qstop == auxiliary_->nshell() ?
                     auxiliary_->nbf() -
                     auxiliary_->shell(Qstart).function_index() :
                     auxiliary_->shell(Qstop).function_index() -
                     auxiliary_->shell(Qstart).function_index());

        // Compute TEI tensor block (A|mn)

        timer_on("JK: (Q|mn)^R");

        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = auxiliary_->shell(Q).nfunction();
            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();

            int sq =  auxiliary_->shell(Q).function_index();
            int sm =  primary_->shell(M).function_index();
            int sn =  primary_->shell(N).function_index();

            eri2[thread]->compute_shell(Q,0,M,N);

            for (int om = 0; om < nm; om++) {
                for (int on = 0; on < nn; on++) {
                    long int m = sm + om;
                    long int n = sn + on;
                    if (m >= n && schwarz_fun_pairs_r[m*(m+1)/2 + n] >= 0) {
                        long int delta = schwarz_fun_pairs_r[m*(m+1)/2 + n];
                        for (int oq = 0; oq < nq; oq++) {
                            Amn2p[sq + oq - qoff][delta] =
                            buffer2[thread][oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }

        timer_off("JK: (Q|mn)^R");

        // Dump block to disk
        timer_on("JK: (Q|mn)^R Write");

        psio_->write(unit_,"Right (Q|w|mn) Integrals",(char*)Amn2p[0],sizeof(double)*nrows*ntri,next_AIA,&next_AIA);

        timer_off("JK: (Q|mn)^R Write");

    }
    Amn2.reset();
    delete[] eri2;

    psio_->write_entry(unit_, "Omega", (char*) &omega_, sizeof(double));
    psio_->close(unit_,1);
}
void DFJK::rebuild_wK_disk()
{
    // Already open
    fprintf(outfile, "    Rebuilding (Q|w|mn) Integrals (new omega)\n\n");

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

    // ==> Thread setup <== //
    int nthread = 1;
    #ifdef _OPENMP
        nthread = df_ints_num_threads_;
    #endif

    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    boost::shared_ptr<IntegralFactory> rifactory(new IntegralFactory(auxiliary_, zero, primary_, primary_));
    const double **buffer2 = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *eri2 = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        eri2[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory->erf_eri(omega_));
        buffer2[Q] = eri2[Q]->buffer();
    }

    ULI maxP = auxiliary_->max_function_per_shell();
    ULI max_rows = memory_ / ntri;
    max_rows = (max_rows > naux ? naux : max_rows);
    max_rows = (max_rows < maxP ? maxP : max_rows);

    // Block extents
    std::vector<int> block_Q_starts;
    int counter = 0;
    block_Q_starts.push_back(0);
    for (int Q = 0; Q < auxiliary_->nshell(); Q++) {
        int nQ = auxiliary_->shell(Q).nfunction();
        if (counter + nQ > max_rows) {
            counter = 0;
            block_Q_starts.push_back(Q);
        }
        counter += nQ;
    }
    block_Q_starts.push_back(auxiliary_->nshell());

    SharedMatrix Amn2(new Matrix("(A|mn) Block", max_rows, ntri));
    double** Amn2p = Amn2->pointer();
    psio_address next_AIA = PSIO_ZERO;

    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();
    const size_t npairs = shell_pairs.size();

    // Loop over blocks of Qshell
    for (int block = 0; block < block_Q_starts.size() - 1; block++) {

        // Block sizing/offsets
        int Qstart = block_Q_starts[block];
        int Qstop  = block_Q_starts[block+1];
        int qoff   = auxiliary_->shell(Qstart).function_index();
        int nrows  = (Qstop == auxiliary_->nshell() ?
                     auxiliary_->nbf() -
                     auxiliary_->shell(Qstart).function_index() :
                     auxiliary_->shell(Qstop).function_index() -
                     auxiliary_->shell(Qstart).function_index());

        // Compute TEI tensor block (A|mn)

        timer_on("JK: (Q|mn)^R");

        #pragma omp parallel for schedule(dynamic) num_threads(nthread)
        for (long int QMN = 0L; QMN < (Qstop - Qstart) * (ULI) npairs; QMN++) {

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int Q =  QMN / npairs + Qstart;
            int MN = QMN % npairs;

            std::pair<int,int> pair = shell_pairs[MN];
            int M = pair.first;
            int N = pair.second;

            int nq = auxiliary_->shell(Q).nfunction();
            int nm = primary_->shell(M).nfunction();
            int nn = primary_->shell(N).nfunction();

            int sq =  auxiliary_->shell(Q).function_index();
            int sm =  primary_->shell(M).function_index();
            int sn =  primary_->shell(N).function_index();

            eri2[thread]->compute_shell(Q,0,M,N);

            for (int om = 0; om < nm; om++) {
                for (int on = 0; on < nn; on++) {
                    long int m = sm + om;
                    long int n = sn + on;
                    if (m >= n && schwarz_fun_pairs_r[m*(m+1)/2 + n] >= 0) {
                        long int delta = schwarz_fun_pairs_r[m*(m+1)/2 + n];
                        for (int oq = 0; oq < nq; oq++) {
                            Amn2p[sq + oq - qoff][delta] =
                            buffer2[thread][oq * nm * nn + om * nn + on];
                        }
                    }
                }
            }
        }

        timer_off("JK: (Q|mn)^R");

        // Dump block to disk
        timer_on("JK: (Q|mn)^R Write");

        psio_->write(unit_,"Right (Q|w|mn) Integrals",(char*)Amn2p[0],sizeof(double)*nrows*ntri,next_AIA,&next_AIA);

        timer_off("JK: (Q|mn)^R Write");

    }
    Amn2.reset();
    delete[] eri2;

    psio_->write_entry(unit_, "Omega", (char*) &omega_, sizeof(double));
    // No need to close
}
void DFJK::manage_JK_core()
{
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);

        if (do_J_) {
            timer_on("JK: J");
            block_J(&Qmn_->pointer()[Q],naux);
            timer_off("JK: J");
        }
        if (do_K_) {
            timer_on("JK: K");
            block_K(&Qmn_->pointer()[Q],naux);
            timer_off("JK: K");
        }
    }
}
void DFJK::manage_JK_disk()
{
    int ntri = sieve_->function_pairs().size();
    Qmn_ = SharedMatrix(new Matrix("(Q|mn) Block", max_rows_, ntri));
    psio_->open(unit_,PSIO_OPEN_OLD);
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_ ? auxiliary_->nbf() - Q : max_rows_);
        psio_address addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri) * sizeof(double));

        timer_on("JK: (Q|mn) Read");
        psio_->read(unit_,"(Q|mn) Integrals", (char*)(Qmn_->pointer()[0]),sizeof(double)*naux*ntri,addr,&addr);
        timer_off("JK: (Q|mn) Read");

        if (do_J_) {
            timer_on("JK: J");
            block_J(&Qmn_->pointer()[0],naux);
            timer_off("JK: J");
        }
        if (do_K_) {
            timer_on("JK: K");
            block_K(&Qmn_->pointer()[0],naux);
            timer_off("JK: K");
        }
    }
    psio_->close(unit_,1);
    Qmn_.reset();
}
void DFJK::manage_wK_core()
{
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_w) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_w ? auxiliary_->nbf() - Q : max_rows_w);

        timer_on("JK: wK");
        block_wK(&Qlmn_->pointer()[Q],&Qrmn_->pointer()[Q],naux);
        timer_off("JK: wK");
    }
}
void DFJK::manage_wK_disk()
{
    int max_rows_w = max_rows_ / 2;
    max_rows_w = (max_rows_w < 1 ? 1 : max_rows_w);
    int ntri = sieve_->function_pairs().size();
    Qlmn_ = SharedMatrix(new Matrix("(Q|mn) Block", max_rows_w, ntri));
    Qrmn_ = SharedMatrix(new Matrix("(Q|mn) Block", max_rows_w, ntri));
    psio_->open(unit_,PSIO_OPEN_OLD);
    for (int Q = 0 ; Q < auxiliary_->nbf(); Q += max_rows_w) {
        int naux = (auxiliary_->nbf() - Q <= max_rows_w ? auxiliary_->nbf() - Q : max_rows_w);
        psio_address addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri) * sizeof(double));

        timer_on("JK: (Q|mn)^L Read");
        psio_->read(unit_,"Left (Q|w|mn) Integrals", (char*)(Qlmn_->pointer()[0]),sizeof(double)*naux*ntri,addr,&addr);
        timer_off("JK: (Q|mn)^L Read");

        addr = psio_get_address(PSIO_ZERO, (Q*(ULI) ntri) * sizeof(double));

        timer_on("JK: (Q|mn)^R Read");
        psio_->read(unit_,"Right (Q|w|mn) Integrals", (char*)(Qrmn_->pointer()[0]),sizeof(double)*naux*ntri,addr,&addr);
        timer_off("JK: (Q|mn)^R Read");

        timer_on("JK: wK");
        block_wK(&Qlmn_->pointer()[0],&Qrmn_->pointer()[0],naux);
        timer_off("JK: wK");
    }
    psio_->close(unit_,1);
    Qlmn_.reset();
    Qrmn_.reset();
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

        timer_on("JK: J1");
        C_DGEMV('N',naux,num_nm,1.0,Qmnp[0],num_nm,D2p,1,0.0,dp,1);
        timer_off("JK: J1");

        timer_on("JK: J2");
        C_DGEMV('T',naux,num_nm,1.0,Qmnp[0],num_nm,dp,1,0.0,J2p,1);
        timer_off("JK: J2");

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

    for (int N = 0; N < K_ao_.size(); N++) {

        int nbf = C_left_ao_[N]->rowspi()[0];
        int nocc = C_left_ao_[N]->colspi()[0];

        if (!nocc) continue;

        double** Clp  = C_left_ao_[N]->pointer();
        double** Crp  = C_right_ao_[N]->pointer();
        double** Elp  = E_left_->pointer();
        double** Erp  = E_right_->pointer();
        double** Kp   = K_ao_[N]->pointer();

        if (N == 0 || C_left_[N].get() != C_left_[N-1].get()) {

            timer_on("JK: K1");

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

                if (nocc > 1) {
                    C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Elp[0][m*(ULI)nocc*naux],naux);
                } else {
                    C_DGEMV('N',naux,rows,1.0,QSp[0],nbf,Clp[0],1,0.0,&Elp[0][m*(ULI)nocc*naux],1);
                }
            }

            timer_off("JK: K1");

        }

        if (!lr_symmetric_ && (N == 0 || C_right_[N].get() != C_right_[N-1].get())) {

            timer_on("JK: K1");

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

                if (nocc > 1) {
                    C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Erp[0][m*(ULI)nocc*naux],naux);
                } else {
                    C_DGEMV('N',naux,rows,1.0,QSp[0],nbf,Clp[0],1,0.0,&Erp[0][m*(ULI)nocc*naux],1);
                }
            }

            timer_off("JK: K1");

        }

        timer_on("JK: K2");
        C_DGEMM('N','T',nbf,nbf,naux*nocc,1.0,Elp[0],naux*nocc,Erp[0],naux*nocc,1.0,Kp[0],nbf);
        timer_off("JK: K2");
    }

}
void DFJK::block_wK(double** Qlmnp, double** Qrmnp, int naux)
{
    const std::vector<std::pair<int, int> >& function_pairs = sieve_->function_pairs();
    const std::vector<long int>& function_pairs_reverse = sieve_->function_pairs_reverse();
    unsigned long int num_nm = function_pairs.size();

    for (int N = 0; N < wK_ao_.size(); N++) {

        int nbf = C_left_ao_[N]->rowspi()[0];
        int nocc = C_left_ao_[N]->colspi()[0];

        if (!nocc) continue;

        double** Clp  = C_left_ao_[N]->pointer();
        double** Crp  = C_right_ao_[N]->pointer();
        double** Elp  = E_left_->pointer();
        double** Erp  = E_right_->pointer();
        double** wKp   = wK_ao_[N]->pointer();

        if (N == 0 || C_left_[N].get() != C_left_[N-1].get()) {

            timer_on("JK: wK1");

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
                    C_DCOPY(naux,&Qlmnp[0][ij],num_nm,&QSp[0][i],nbf);
                    C_DCOPY(nocc,Clp[n],1,&Ctp[0][i],nbf);
                }

                C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Elp[0][m*(ULI)nocc*naux],naux);
            }

            timer_off("JK: wK1");

        }

        timer_on("JK: wK1");

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
                C_DCOPY(naux,&Qrmnp[0][ij],num_nm,&QSp[0][i],nbf);
                C_DCOPY(nocc,Crp[n],1,&Ctp[0][i],nbf);
            }

            C_DGEMM('N','T',nocc,naux,rows,1.0,Ctp[0],nbf,QSp[0],nbf,0.0,&Erp[0][m*(ULI)nocc*naux],naux);
        }

        timer_off("JK: wK1");

        timer_on("JK: wK2");
        C_DGEMM('N','T',nbf,nbf,naux*nocc,1.0,Elp[0],naux*nocc,Erp[0],naux*nocc,1.0,wKp[0],nbf);
        timer_off("JK: wK2");
    }
}

PSJK::PSJK(boost::shared_ptr<BasisSet> primary, 
    Options& options) :
    JK(primary_), options_(options)
{
    common_init();
}
PSJK::~PSJK()
{
}
void PSJK::common_init()
{
    df_ints_num_threads_ = 1;
    #ifdef _OPENMP
        df_ints_num_threads_ = omp_get_max_threads();
    #endif
    unit_ = PSIF_DFSCF_BJ;
    theta_ = 0.3;
    psio_ = PSIO::shared_object();
    dealiasing_ = "QUADRATURE";
}
void PSJK::print_header() const
{
    if (print_) {
        fprintf(outfile, "  ==> PSJK: Pseudospectral J/K Matrices <==\n\n");

        fprintf(outfile, "    J tasked:          %11s\n", (do_J_ ? "Yes" : "No"));
        fprintf(outfile, "    K tasked:          %11s\n", (do_K_ ? "Yes" : "No"));
        fprintf(outfile, "    wK tasked:         %11s\n", (do_wK_ ? "Yes" : "No"));
        fprintf(outfile, "    OpenMP threads:    %11d\n", omp_nthread_);
        fprintf(outfile, "    Integrals threads: %11d\n", df_ints_num_threads_);
        fprintf(outfile, "    Memory (MB):       %11ld\n", (memory_ *8L) / (1024L * 1024L));
        fprintf(outfile, "    Schwarz Cutoff:    %11.0E\n", cutoff_);
        fprintf(outfile, "    Theta:             %11.3E\n", theta_);
        fprintf(outfile, "    Dealiasing:        %11s\n", dealiasing_.c_str());
        fprintf(outfile, "\n");

        if (dealiasing_ == "DEALIAS") {
            fprintf(outfile, "   => Dealias Basis Set <=\n\n");
            dealias_->print_by_level(outfile,print_);
        }
    }
}
void PSJK::preiterations()
{
    // PS requires constant sieve, must be static througout object life
    if (!sieve_) {
        sieve_ = boost::shared_ptr<ERISieve>(new ERISieve(primary_, cutoff_));
    }

    build_QR();

    if (do_J_ || do_K_) {
        build_Amn_disk(theta_,"(A|mn) JK");    
    }

    if (do_wK_) {
        throw FeatureNotImplemented("PSJK", "wK", __FILE__, __LINE__);
    }
}
void PSJK::postiterations()
{
    Q_.reset();
    R_.reset();
    grid_.reset();
}
void PSJK::compute_JK()
{
    int nbf  = primary_->nbf();
    int naux = grid_->rowspi()[0];
    int ntri = sieve_->function_pairs().size();
    int maxrows = max_rows();

    // Long Range
    SharedMatrix Amn(new Matrix("Amn",maxrows,ntri)); 
    double** Amnp = Amn->pointer();    

    V_ = SharedMatrix(new Matrix("V", maxrows,nbf));
    W_ = SharedMatrix(new Matrix("W", maxrows,nbf));
    d_ = SharedVector(new Vector("d", maxrows));;
    J_temp_ = SharedVector(new Vector("J temp", ntri));

    if (do_J_ || do_K_) {
        // Long Range
        psio_->open(unit_, PSIO_OPEN_OLD);
        psio_address addr = PSIO_ZERO;
        for (int Pstart = 0; Pstart < naux; Pstart += maxrows) {
            int nrows = (naux - Pstart >= maxrows ? naux - Pstart : maxrows);
            psio_->read(unit_,"(A|mn) JK", (char*) Amnp[0], sizeof(double)*nrows*ntri,addr,&addr); 
            if (do_J_) {
                block_J(Amnp,Pstart,nrows,J_ao_);
            }
            if (do_K_) {
                block_K(Amnp,Pstart,nrows,K_ao_);
            }
        }
        psio_->close(unit_, 1);
    }
    
    Amn.reset();
    V_.reset();
    W_.reset();
    d_.reset();
    J_temp_.reset();
    
    // Short Range 
    build_JK_SR(); 

    // TODO wK
}
void PSJK::build_QR()
{
    boost::shared_ptr<PseudospectralGrid> grid(new PseudospectralGrid(primary_->molecule(), 
        primary_, options_));
    int npoints = grid->npoints(); 
    int nbf = primary_->nbf();
    int max_points = grid->max_points();
    int max_functions = grid->max_functions();

    // Grid
    grid_ = SharedMatrix(new Matrix("xyzw", npoints, 4));
    double** gridp = grid_->pointer();
    double* x = grid->x();
    double* y = grid->y();
    double* z = grid->z();
    double* w = grid->w();
    for (int P = 0; P < npoints; P++) {
        gridp[P][0] = x[P]; 
        gridp[P][1] = y[P]; 
        gridp[P][2] = z[P]; 
        gridp[P][3] = w[P]; 
    } 

    // R (Collocation)
    R_ = SharedMatrix(new Matrix("R", nbf, npoints));
    double** Rp = R_->pointer();
    boost::shared_ptr<PointFunctions> points(new PointFunctions(primary_, max_points, max_functions));
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid->blocks();
    int offset = 0;
    for (int index = 0; index < blocks.size(); index++) {
        points->computePoints(blocks[index]);
        SharedMatrix phi = points->basis_value("PHI");
        double** phip = phi->pointer();
        const std::vector<int>& funmap = blocks[index]->functions_local_to_global();
        int nP = blocks[index]->npoints();

        for (int i = 0; i < funmap.size(); i++) {
            int iglobal = funmap[i];
            C_DCOPY(nP,&phip[0][i],max_functions,&Rp[iglobal][offset],1);
        }
        
        offset += nP;
    }

    points.reset();
    grid.reset();
    
    // Q (Quadrature Rule)
    if (dealiasing_ == "QUADRATURE") {
        for (int P = 0; P < npoints; P++) {
            C_DSCAL(nbf,sqrt(gridp[P][3]),&Rp[0][P],npoints);
        }
        Q_ = R_;
    } else if (dealiasing_ == "RENORMALIZATION") {

        for (int P = 0; P < npoints; P++) {
            C_DSCAL(nbf,sqrt(gridp[P][3]),&Rp[0][P],npoints);
        }

        bool warning;     
        SharedMatrix Rplus = R_->pseudoinverse(1.0E-10, &warning);
        if (warning) {
            fprintf(outfile, "    Warning, Renormalization had to be conditioned.\n\n");
        } 

        boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_));
        boost::shared_ptr<OneBodyAOInt> ints(factory->ao_overlap());
        SharedMatrix S(new Matrix("S", nbf, nbf));
        ints->compute(S);
        ints.reset();
        
        Q_ = SharedMatrix(R_->clone());
        Q_->set_name("Q"); 
        
        double** Qp = Q_->pointer();
        double** Rp = Rplus->pointer();
        double** Sp = S->pointer();
        
        C_DGEMM('N','N',nbf,npoints,nbf,1.0,Sp[0],nbf,Rp[0],npoints,0.0,Qp[0],npoints);

    } else if (dealiasing_ == "DEALIAS") {
        // TODO 
        throw FeatureNotImplemented("PSJK", "Dealiasing", __FILE__, __LINE__);
    }
}
void PSJK::build_Amn_disk(double theta, const std::string& entry)
{
    boost::shared_ptr<IntegralFactory> factory(new IntegralFactory(primary_));
    std::vector<boost::shared_ptr<PseudospectralInt> > ints;
    for (int i = 0; i < df_ints_num_threads_; i++) {
        ints.push_back(boost::shared_ptr<PseudospectralInt>(static_cast<PseudospectralInt*>(factory->ao_pseudospectral())));
        ints[i]->set_omega(theta);        
    }

    int maxrows = max_rows();
    int ntri  = sieve_->function_pairs().size();
    int naux  = grid_->rowspi()[0];

    const std::vector<long int>& function_pairs_r = sieve_->function_pairs_reverse();
    const std::vector<std::pair<int,int> >& shell_pairs = sieve_->shell_pairs();

    SharedMatrix Amn(new Matrix("Amn", maxrows, ntri));
    double** Amnp = Amn->pointer();
    double** Gp = grid_->pointer();

    psio_->open(unit_,PSIO_OPEN_OLD);
    psio_address addr = PSIO_ZERO;
    
    for (int Pstart = 0; Pstart < naux; Pstart += maxrows) {
        int nrows =  (Pstart + maxrows >= naux ? naux - Pstart : maxrows);

        #pragma omp parallel for num_threads(df_ints_num_threads_)
        for (int row = 0; row < (unsigned int)nrows; row++) {
            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif
            boost::shared_ptr<PseudospectralInt> eri = ints[thread];
            const double* buffer = eri->buffer();
            eri->set_point(Gp[row + Pstart][0], Gp[row + Pstart][1], Gp[row + Pstart][2]);
            double* Amp = Amnp[row]; 
            
            for (int ind = 0; ind < shell_pairs.size(); ind++) {
                int M = shell_pairs[ind].first;
                int N = shell_pairs[ind].second;
                
                eri->compute_shell(M,N);

                int nM = primary_->shell(M).nfunction();
                int nN = primary_->shell(N).nfunction();
                int oM = primary_->shell(M).function_index();
                int oN = primary_->shell(N).function_index();

                for (int m = 0 ; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        int am = m + oM;
                        int an = n + oN;
                        if (am >= an) { 
                            int mn = (am *(am + 1) >> 1) + an;
                            long int mn_local = function_pairs_r[mn];
                            if (mn >= 0) {
                                Amp[mn_local] = buffer[m * nN + n];
                            }
                        }
                    }
                }
            }
        }
        psio_->write(unit_, entry.c_str(), (char*) Amnp[0], sizeof(double) * ntri * nrows, addr, &addr);
    }

    psio_->close(unit_, 1);
}
int PSJK::max_rows() 
{
    // TODO
    return 1;
}
void PSJK::block_J(double** Amnp, int Pstart, int nP, const std::vector<SharedMatrix>& J)
{
    int nbf = primary_->nbf();
    int npoints = grid_->rowspi()[0];
    const std::vector<std::pair<int,int> >& funmap = sieve_->function_pairs();
    int ntri = funmap.size();

    double** Qp = Q_->pointer();
    double** Rp = R_->pointer();
    double** Vp = V_->pointer();
    double*  dp = d_->pointer();
    double* J2p = J_temp_->pointer();

    for (int N = 0; N < D_ao_.size(); N++) {
        
        double** Dp = D_ao_[N]->pointer();
        double** Jp = J[N]->pointer();

        C_DGEMM('T','N',nP,nbf,nbf,1.0,&Qp[0][Pstart],npoints,Dp[0],nbf,0.0,Vp[0],nbf);
        for (int P = 0; P < nP; P++) {
            dp[P] = C_DDOT(nbf,&Rp[0][Pstart + P],npoints,Vp[P],1);
        }

        C_DGEMV('T',nP,ntri,1.0,Amnp[0],ntri,dp,1,0.0,J2p,1);

        for (long int mn = 0; mn < ntri; mn++) {
            int m = funmap[mn].first;
            int n = funmap[mn].second;
            Jp[m][n] += J2p[mn];
            Jp[n][m] += (m == n ? 0.0 : J2p[mn]);
        }
    }
}
void PSJK::block_K(double** Amnp, int Pstart, int nP, const std::vector<SharedMatrix>& K)
{
    int nbf = primary_->nbf();
    int npoints = grid_->rowspi()[0];
    const std::vector<std::pair<int,int> >& funmap = sieve_->function_pairs();
    int ntri = funmap.size();

    double** Qp = Q_->pointer();
    double** Rp = R_->pointer();
    double** Vp = V_->pointer();
    double** Wp = W_->pointer();

    for (int N = 0; N < D_ao_.size(); N++) {

        double** Dp = D_ao_[N]->pointer();
        double** Kp = K[N]->pointer();

        C_DGEMM('T','N',nP,nbf,nbf,1.0,&Qp[0][Pstart],npoints,Dp[0],nbf,0.0,Vp[0],nbf);

        W_->zero();
        #pragma omp parallel for num_threads(omp_nthread_)
        for (int P = 0; P < nP; P++) {
            double* Arp = Amnp[P];
            double* Wrp = Wp[P];
            double* Vrp = Vp[P];
            for (long int mn = 0; mn < ntri; mn++) {
                int m = funmap[mn].first;
                int n = funmap[mn].second;
                double Aval = Arp[mn];
                Wrp[m] += Aval * Vrp[n]; 
                Wrp[n] += (m == n ? 0.0 : Aval * Vrp[m]); 
            }
        }

        C_DGEMM('T','N',nbf,nbf,nP,1.0,Vp[0],nbf,Wp[0],nbf,1.0,Kp[0],nbf);
    }
}
void PSJK::build_JK_SR()
{
}

}
