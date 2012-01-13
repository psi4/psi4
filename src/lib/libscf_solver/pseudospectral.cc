#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#include "pseudospectral.h"
#include <psiconfig.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

#include <libmints/mints.h>
#include <lib3index/3index.h>

using namespace std;
using namespace psi;
using namespace boost;

namespace psi { namespace scf {

PseudospectralHF::PseudospectralHF(boost::shared_ptr<BasisSet> basis, SharedMatrix Da,
SharedMatrix Ja, SharedMatrix Ka, boost::shared_ptr<PSIO> psio, Options& opt) :
    primary_(basis), Da_(Da), Ja_(Ja), Ka_(Ka), psio_(psio), options_(opt), restricted_(true)
{
    common_init();
}

PseudospectralHF::PseudospectralHF(boost::shared_ptr<BasisSet> basis, boost::shared_ptr<PSIO> psio, Options& opt) :
    primary_(basis), psio_(psio), options_(opt), restricted_(true)
{
    common_init();
}

PseudospectralHF::~PseudospectralHF()
{
}
void PseudospectralHF::common_init()
{
    print_ = options_.get_int("PRINT");
    fprintf(outfile, " PseudospectalHF: Pseudospectral SCF Algorithms (In Progress).\n");
    fprintf(outfile, "   by Rob Parrish\n\n");

    // How many doubles do we have?
    memory_ = Process::environment.get_memory() / 8L;
    memory_ = (unsigned long int) 0.7 * memory_;

    // Build auxiliary basis from options
    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    auxiliary_ = BasisSet::construct(parser, primary_->molecule(), "DF_BASIS_SCF");
    parser.reset();

    // Build a fitting metric object
    Jinv_ = boost::shared_ptr<FittingMetric>(new FittingMetric(auxiliary_));
    Jinv_->form_cholesky_factor();

    // Build a Schwarz sieve object
    schwarz_ = boost::shared_ptr<SchwarzSieve>(new SchwarzSieve(primary_, options_.get_double("INTS_TOLERANCE")));

    // Build a set of thread-local integrators
    int nthread = 1;
    #ifdef _OMP
        nthread = omp_get_max_threads();
    #endif
    eri_.resize(nthread);
    pot_.resize(nthread);

    #pragma omp parallel
    {
        int thread = 0;
        #ifdef _OMP
            thread = omp_get_thread_num();
        #endif

        boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(auxiliary_, zero, primary_, primary_));
        eri_[thread] = boost::shared_ptr<TwoBodyAOInt>(fact->eri());
        fact.reset();

        boost::shared_ptr<IntegralFactory> fact2(new IntegralFactory(primary_, primary_, zero, zero));
        pot_[thread] = boost::shared_ptr<PseudospectralInt>(static_cast<PseudospectralInt*>(fact2->ao_pseudospectral()));
    }

    #if 0
    // Build the peudospectral grid points and X matrix
    boost::shared_ptr<Integrator> quad = Integrator::build_integrator(primary_->molecule(), psio_, options_);
    quad->buildGrid(5000);
    P_ = quad->getNPoints();

    boost::shared_ptr<BasisPoints> points(new BasisPoints(primary_, 5000));
    points->setToComputePoints(true);
    double** bpoints = points->getPoints();

    points_ = SharedMatrix(new Matrix("Pseudospectral Points", P_, 3));
    double** rp = points_->pointer(0);

    X_ = SharedMatrix(new Matrix("Pseudospectral X", P_, primary_->nbf()));
    double** Xp = X_->pointer(0);

    unsigned long int counter = 0L;
    for (int P = 0; P < quad->getNBlocks(); P++) {
        boost::shared_ptr<GridBlock> block = quad->getBlock(P);
        double* x = block->getX();
        double* y = block->getY();
        double* z = block->getZ();
        double* w = block->getZ();
        int n = block->getTruePoints();

        // Compute the basis points
        points->computePoints(block);

        // Copy the points in
        for (int i = 0; i < n; i++) {
            rp[counter][0] = x[i];
            rp[counter][1] = y[i];
            rp[counter][2] = z[i];
            double weight = sqrt(w[i]);
            for (int Q = 0; Q < primary_->nbf(); Q++)
                Xp[counter][Q] = weight * bpoints[i][Q];
            counter++;
        }
    }
    #endif
}
void PseudospectralHF::form_G_RHF()
{
    form_J_DF_RHF();
    form_K_PS_RHF();
}
void PseudospectralHF::form_J_DF_RHF()
{
    // Form c_A = (A|mn) D_{mn}
    int max_p = 0;
    for (int P = 0; P < primary_->nshell(); P++)
        if (max_p < primary_->shell(P)->nfunction())
            max_p = primary_->shell(P)->nfunction();
    int nca = max_p*max_p;

    int nthread = 1;
    #ifdef _OMP
        nthread = omp_get_max_threads();
    #endif

    double** Dp = Da_->pointer();
    double** Jp = Ja_->pointer();

    double** Db = new double*[nthread];
    double** db = new double*[nthread];
    double** A = new double*[nthread];

    bool allocated = false;
    int thread = 0;

    int Mshell = primary_->nshell();
    int Ashell = auxiliary_->nshell();
    int naux = auxiliary_->nbf();

    unsigned long int npairs = schwarz_->get_nshell_pairs();
    int* pairs = schwarz_->get_schwarz_shells();

    #ifdef HAVE_MKL
        int mkl_n = mkl_get_max_threads();
        mkl_set_num_threads(1);
    #endif

    // TODO use larger blocks than this, DGEMM will be inefficient, even on 1 thread
    #pragma omp parallel for firstprivate(allocated) private(thread) schedule(guided) num_threads(nthread)
    for (unsigned long int PQ = 0; PQ < npairs; PQ++) {
        if (!allocated) {
            #ifdef _OMP
                thread = omp_get_thread_num();
            #endif
            Db[thread] = new double[nca];
            db[thread] = new double[naux];
            memset(static_cast<void*>(db[thread]), '\0', naux*sizeof(double));
            A[thread] = new double[naux*nca];
            allocated = true;
        }

        int P = pairs[2*PQ];
        int Q = pairs[2*PQ + 1];
        int np = primary_->shell(P)->nfunction();
        int nq = primary_->shell(Q)->nfunction();
        int op = primary_->shell(P)->function_index();
        int oq = primary_->shell(Q)->function_index();

        unsigned long int offset = 0L;
        const double* buffer = eri_[thread]->buffer();

        // Block over all A for this PQ shell
        for (int Ap = 0; Ap < Ashell; Ap++) {
            // Compute the integrals
            eri_[thread]->compute_shell(Ap, 0 , P, Q);
            int na = auxiliary_->shell(Ap)->nfunction();

            memcpy(static_cast<void*>(&A[thread][offset]), buffer, np*nq*na*sizeof(double));

        }

        // grab the relevant bit of the density matrix
        for (int p = 0; p < np; p++) {
            memcpy(static_cast<void*>(&Db[thread][p*nq]), static_cast<void*> (&Dp
                [op + p][oq]), nq*sizeof(double));
        }

        // Permutational symmetry
        if (P != Q) C_DSCAL(np*nq, 2.0, Db[thread], 1);

        // do the GEMV
        C_DGEMV('N', naux, np*nq, 1.0, A[thread], np*nq, Db[thread], 1, 1.0, db[thread], 1);

    }

    #ifdef HAVE_MKL
        mkl_set_num_threads(mkl_n);
    #endif

    // Add all the threads into the one element
    double* d = new double[naux];
    memset(static_cast<void*>(d), '\0', naux*sizeof(double));
    for (int t = 0; t < nthread; t++) {
        C_DAXPY(naux, 1.0, db[t], 1, d, 1);
        delete[] A[thread];
        delete[] db[thread];
        delete[] Db[thread];
    }
    delete[] Db;
    delete[] db;
    delete[] A;

    // Form d_B = J_AB^-1 c_A
    double** Wp = Jinv_->get_metric()->pointer();
    int info = C_DPOTRS('U',naux,1,Wp[0],naux,d,naux);

    // Form d_B = (A|mn) D_{mn}
    #ifdef HAVE_MKL
        mkl_set_num_threads(1);
    #endif

    allocated = false;

    // TODO use larger blocks than this, DGEMM will be inefficient, even on 1 thread
    #pragma omp parallel for firstprivate(allocated) private(thread) schedule(guided) num_threads(nthread)
    for (unsigned long int PQ = 0; PQ < npairs; PQ++) {
        if (!allocated) {
            #ifdef _OMP
                thread = omp_get_thread_num();
            #endif
            db[thread] = new double[naux];
            memcpy(static_cast<void*>(d), static_cast<void*>(db[thread]), naux*sizeof(double));
            Db[thread] = new double[nca];
            A[thread] = new double[naux*nca];
            allocated = true;
        }

        int P = pairs[2*PQ];
        int Q = pairs[2*PQ + 1];
        int np = primary_->shell(P)->nfunction();
        int nq = primary_->shell(Q)->nfunction();
        int op = primary_->shell(P)->function_index();
        int oq = primary_->shell(Q)->function_index();

        unsigned long int offset = 0L;
        const double* buffer = eri_[thread]->buffer();

        // Block over all A for this PQ shell
        for (int Ap = 0; Ap < Ashell; Ap++) {
            // Compute the integrals
            eri_[thread]->compute_shell(Ap, 0 , P, Q);
            int na = auxiliary_->shell(Ap)->nfunction();

            memcpy(static_cast<void*>(&A[thread][offset]), buffer, np*nq*na*sizeof(double));

        }

        // do the GEMV
        C_DGEMV('T', np*nq, naux, 1.0, A[thread], np*nq, db[thread], 1, 0.0, Db[thread], 1);

        // Fill the relevant bit of the J matrix
        for (int p = 0; p < np; p++) {
            memcpy(static_cast<void*>(&Jp[op + p][oq]), static_cast<void*> (&Db[thread][p*nq]), nq*sizeof(double));
        }

    }

    for (int m = 1; m < primary_->nbf(); m++)
        for (int n = m; n < primary_->nbf(); n++)
            Jp[m][n] = Jp[n][m];

    #ifdef HAVE_MKL
        mkl_set_num_threads(mkl_n);
    #endif

    delete[] d;
    delete[] db;
    delete[] A;

}
void PseudospectralHF::form_K_PS_RHF()
{
    // Cheesy N^3 approach for now
    // We can exploit sparsity like crazy
    SharedMatrix Q(new Matrix("Q", P_, primary_->nbf()));

    C_DGEMM('N', 'N', P_, primary_->nbf(), primary_->nbf(), 1.0, X_->pointer()[0], primary_->nbf(),
        Da_->pointer()[0], primary_->nbf(), 0.0, Q->pointer()[0], primary_->nbf());

    int nthread = 1;
    #ifdef _OMP
        nthread = omp_get_max_threads();
    #endif

    SharedMatrix Z(new Matrix("Z", P_, primary_->nbf()));
    std::vector<SharedMatrix > A;
    A.resize(nthread);

    #pragma omp parallel
    {
        int thread = 0;
        #ifdef _OMP
            thread = omp_thread_num();
        #endif
        A[thread] = SharedMatrix(new Matrix("A", primary_->nbf(), primary_->nbf()));
    }

    #ifdef HAVE_MKL
        int mkl_n = mkl_get_max_threads();
        mkl_set_num_threads(1);
    #endif

    double** Rp = points_->pointer();
    double** Zp = Z->pointer();
    double** Qp = Q->pointer();

    #pragma omp parallel for schedule(guided) num_threads(nthread)
    for (unsigned long int P = 0L; P < P_; P++) {

        int thread = 0;
        #ifdef _OMP
            thread = omp_thread_num();
        #endif

        pot_[thread]->set_point(Rp[P][0], Rp[P][1], Rp[P][2]);
        pot_[thread]->compute(A[thread]);

        double** Ap = A[thread]->pointer();
        C_DGEMV('N', primary_->nbf(), primary_->nbf(), 1.0, Ap[0], primary_->nbf(), Qp[P],
            1, 0.0, Zp[P], 1);
    }
    A.clear();

    C_DGEMM('T', 'N', primary_->nbf(), primary_->nbf(), P_, 1.0, Z->pointer()[0], primary_->nbf(),
        X_->pointer()[0], primary_->nbf(), 0.0, Ka_->pointer()[0], primary_->nbf());

    double** Kp = Ka_->pointer();
    for (int m = 1; m < primary_->nbf(); m++) {
        for (int n = m; n < primary_->nbf(); n++) {
            Kp[m][n] = 0.5 * Kp[n][m] + 0.5 * Kp[m][n];
            Kp[n][m] = Kp[m][n];
        }
    }

    #ifdef HAVE_MKL
        mkl_set_num_threads(mkl_n);
    #endif

}

void PseudospectralHF::form_J_DF_UHF()
{
    return; // TODO implement UHF!
}

void PseudospectralHF::form_K_PS_UHF()
{
    return; // TODO implement UHF!
}

}}
