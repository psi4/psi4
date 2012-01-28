#include "3index.h"

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include <libmints/mints.h>
#include <libmints/psimath.h>

//MKL Header
#ifdef HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

namespace psi {

FittingMetric::FittingMetric()
{
    Options& opt = Process::environment.options;
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    if (molecule.get() == 0) {
        fprintf(outfile, "  Active molecule not set!");
        throw PSIEXCEPTION("Active molecule not set!");
    }

    // Make sure molecule is valid.
    molecule->update_geometry();

    // Grab the auxiliary bases (use the SCF tags for now)
    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    aux_ = BasisSet::construct(parser, molecule, "DF_BASIS_SCF");
    if (opt.get_str("POISSON_BASIS_SCF") != "") {
        pois_ =  BasisSet::construct(parser, molecule, "POISSON_BASIS_SCF");
        is_poisson_ = true;
    } else {
        is_poisson_ = false;
    }
    is_inverted_ = false;
    force_C1_ = false;
}
FittingMetric::FittingMetric(boost::shared_ptr<BasisSet> aux, bool force_C1) :
    aux_(aux), is_poisson_(false), is_inverted_(false), force_C1_(force_C1), omega_(0.0)
{
}
FittingMetric::FittingMetric(boost::shared_ptr<BasisSet> aux, double omega, bool force_C1) :
    aux_(aux), is_poisson_(false), is_inverted_(false), force_C1_(force_C1), omega_(omega)
{
}
FittingMetric::FittingMetric(boost::shared_ptr<BasisSet> aux, boost::shared_ptr<BasisSet> pois, bool force_C1) :
    aux_(aux), pois_(pois), is_poisson_(true), is_inverted_(false), force_C1_(force_C1), omega_(0.0)
{
}

FittingMetric::~FittingMetric()
{
}

void FittingMetric::form_fitting_metric()
{
    is_inverted_ = false;
    algorithm_ = "NONE";

    // Sizing/symmetry indexing
    boost::shared_ptr<IntegralFactory> auxfact(new IntegralFactory(aux_, aux_, aux_, aux_));
    boost::shared_ptr<PetiteList> auxpet(new PetiteList(aux_, auxfact));
    boost::shared_ptr<IntegralFactory> poisfact;
    boost::shared_ptr<PetiteList> poispet;
    if (is_poisson_) {
        poisfact = boost::shared_ptr<IntegralFactory>(new IntegralFactory(pois_, pois_, pois_, pois_));
        poispet = boost::shared_ptr<PetiteList>(new PetiteList(pois_, poisfact));
    }

    int naux = 0;
    int ngaussian = 0;
    int npoisson = 0;
    Dimension nauxpi(auxpet->nirrep(), "Fitting Metric Dimensions");
    for (int h = 0; h < auxpet->nirrep(); h++) {
        naux += auxpet->SO_basisdim()[h];
        ngaussian += auxpet->SO_basisdim()[h];
        nauxpi[h] = auxpet->SO_basisdim()[h];
        if (is_poisson_) {
            naux += poispet->SO_basisdim()[h];
            npoisson += poispet->SO_basisdim()[h];
            nauxpi[h] += poispet->SO_basisdim()[h];
        }
    }
    Dimension ngauspi = auxpet->SO_basisdim();

    // Build the full DF/Poisson matrix in the AO basis first
    SharedMatrix AOmetric(new Matrix("AO Basis DF Metric", naux, naux));
    double** W = AOmetric->pointer(0);
    boost::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    // Only thread if not already in parallel (handy for local fitting)
    int nthread = 1;
    #ifdef _OPENMP
        if (!omp_in_parallel()) {
            nthread = omp_get_max_threads();
        }
    #endif

    // == (A|B) Block == //
    IntegralFactory rifactory_J(aux_, zero, aux_, zero);
    const double **Jbuffer = new const double*[nthread];
    boost::shared_ptr<TwoBodyAOInt> *Jint = new boost::shared_ptr<TwoBodyAOInt>[nthread];
    for (int Q = 0; Q<nthread; Q++) {
        if (omega_ > 0.0) {
            Jint[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory_J.erf_eri(omega_));
        } else {
            Jint[Q] = boost::shared_ptr<TwoBodyAOInt>(rifactory_J.eri());
        }
        Jbuffer[Q] = Jint[Q]->buffer();
    }

    #pragma omp parallel for schedule (dynamic) num_threads(nthread)
    for (int MU=0; MU < aux_->nshell(); ++MU) {
        int nummu = aux_->shell(MU)->nfunction();

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        for (int NU=0; NU <= MU; ++NU) {
            int numnu = aux_->shell(NU)->nfunction();

            Jint[thread]->compute_shell(MU, 0, NU, 0);

            int index = 0;
            for (int mu=0; mu < nummu; ++mu) {
                int omu = aux_->shell(MU)->function_index() + mu;

                for (int nu=0; nu < numnu; ++nu, ++index) {
                    int onu = aux_->shell(NU)->function_index() + nu;

                    W[omu][onu] = Jbuffer[thread][index];
                    W[onu][omu] = Jbuffer[thread][index];
                }
            }
        }
    }
    delete[] Jbuffer;
    delete[] Jint;

    if (is_poisson_) {
        // == (AB) Block == //
        IntegralFactory rifactory_RP(pois_, aux_,  zero, zero);
        const double **Obuffer = new const double*[nthread];
        boost::shared_ptr<OneBodyAOInt> *Oint = new boost::shared_ptr<OneBodyAOInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            Oint[Q] = boost::shared_ptr<OneBodyAOInt>(rifactory_RP.ao_overlap());
            Obuffer[Q] = Oint[Q]->buffer();
        }

        #pragma omp parallel for schedule (dynamic) num_threads(nthread)
        for (int NU=0; NU < pois_->nshell(); ++NU) {
            int numnu = pois_->shell(NU)->nfunction();

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            for (int MU=0; MU < aux_->nshell(); ++MU) {
                int nummu = aux_->shell(MU)->nfunction();

                Oint[thread]->compute_shell(NU, MU);

                int index = 0;
                for (int nu=0; nu < numnu; ++nu) {
                    int onu = pois_->shell(NU)->function_index() + nu;

                    for (int mu=0; mu < nummu; ++mu, ++index) {
                        int omu = aux_->shell(MU)->function_index() + mu;

                        W[omu][onu + ngaussian] = Obuffer[thread][index];
                        W[onu + ngaussian][omu] = Obuffer[thread][index];
                    }
                }
            }
        }

        // == (A|T|B) Block == //
        IntegralFactory rifactory_P(pois_, pois_,  zero, zero);
        const double **Tbuffer = new const double*[nthread];
        boost::shared_ptr<OneBodyAOInt> *Tint = new boost::shared_ptr<OneBodyAOInt>[nthread];
        for (int Q = 0; Q<nthread; Q++) {
            Tint[Q] = boost::shared_ptr<OneBodyAOInt>(rifactory_P.ao_kinetic());
            Tbuffer[Q] = Tint[Q]->buffer();
        }

        #pragma omp parallel for schedule (dynamic) num_threads(nthread)
        for (int MU=0; MU < pois_->nshell(); ++MU) {
            int nummu = pois_->shell(MU)->nfunction();

            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            for (int NU=0; NU <= MU ; ++NU) {
                int numnu = pois_->shell(NU)->nfunction();

                Tint[thread]->compute_shell(MU, NU);

                int index = 0;
                for (int mu=0; mu < nummu; ++mu) {
                    int omu = pois_->shell(MU)->function_index() + mu;

                    for (int nu=0; nu < numnu; ++nu, ++index) {
                        int onu = pois_->shell(NU)->function_index() + nu;

                        // These integrals are (A | -1/2 \nabla^2 | B), and should be (A | - 1 / (4 * pi) \nabla^2 |B)
                        // So a factor of 1 / (2 * PI) is the difference
                        W[omu + ngaussian][onu + ngaussian] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                        W[onu + ngaussian][omu + ngaussian] = 1.0 / (2.0 * M_PI) * Tbuffer[thread][index];
                    }
                }
            }
        }
        delete[] Tbuffer;
        delete[] Tint;
    }

    // If C1, form indexing and exit immediately (multiplying by 1 is not so gratifying)
    if (auxpet->nirrep() == 1 || force_C1_ == true) {
        metric_ = AOmetric;
        metric_->set_name("SO Basis Fitting Metric");
        pivots_ = boost::shared_ptr<IntVector>(new IntVector(naux));
        rev_pivots_ = boost::shared_ptr<IntVector>(new IntVector(naux));
        int* piv = pivots_->pointer();
        int* rpiv = pivots_->pointer();
        for (int Q = 0; Q < naux; Q++) {
            piv[Q] = Q;
            rpiv[Q] = Q;
        }
        return;
    }

    // Get the similarity transform objects
    SharedMatrix auxAO2USO(auxpet->sotoao());
    //auxAO2USO->print();
    SharedMatrix poisAO2USO;
    if (is_poisson_) {
        poisAO2USO = SharedMatrix(poispet->sotoao());
        //poisAO2USO->print();
    }

    // Allocate the fitting metric
    metric_ = SharedMatrix(new Matrix("SO Basis Fitting Metric", nauxpi, nauxpi));
    SharedMatrix Temp;
    double** Temp1;

    // Transform AO to SO
    for (int h = 0; h < auxpet->nirrep(); h++) {

        // Gaussian-Gaussian part
        double** J = metric_->pointer(h);
        double** auxU = auxAO2USO->pointer(h);

        if (ngauspi[h] != 0) {
            Temp = SharedMatrix(new Matrix("Temp", ngauspi[h], ngaussian));
            Temp1 = Temp->pointer();
            C_DGEMM('N', 'N', ngauspi[h], ngaussian, ngaussian, 1.0, auxU[0], ngaussian, W[0], naux, 0.0, Temp1[0], ngaussian);
            C_DGEMM('N', 'T', ngauspi[h], ngauspi[h], ngaussian, 1.0, Temp1[0], ngaussian, auxU[0], ngaussian, 0.0, J[0], nauxpi[h]);
            Temp.reset();
        }

        if (is_poisson_ && poispet->SO_basisdim()[h] != 0) {
           Dimension npoispi = poispet->SO_basisdim();
           double** poisU = poisAO2USO->pointer(h);

            // Gaussian-Poisson part
            if (ngauspi[h] != 0) {
                Temp = SharedMatrix(new Matrix("Temp", ngauspi[h], npoisson));
                Temp1 = Temp->pointer();
                C_DGEMM('N', 'N', ngauspi[h], npoisson, ngaussian, 1.0, auxU[0], ngaussian, &W[0][ngaussian], naux, 0.0, Temp1[0], npoisson);
                C_DGEMM('N', 'T', ngauspi[h], npoispi[h], npoisson, 1.0, Temp1[0], npoisson, poisU[0], npoisson, 0.0, &J[0][ngauspi[h]], nauxpi[h]);
                for (int Q = 0; Q < ngauspi[h]; Q++)
                    for (int P =0; P < npoispi[h]; P++)
                        J[P + ngauspi[h]][Q] = J[Q][P + ngauspi[h]];
                Temp.reset();
            }

            // Poisson-Poisson part
            unsigned long int AOoffset = ngaussian*(unsigned long int)naux + (unsigned long int) ngaussian;
            unsigned long int SOoffset = ngauspi[h]*(unsigned long int)nauxpi[h] + (unsigned long int) ngauspi[h];
            Temp = SharedMatrix(new Matrix("Temp", npoispi[h], npoisson));
            Temp1 = Temp->pointer();
            C_DGEMM('N', 'N', npoispi[h], npoisson, npoisson, 1.0, poisU[0], npoisson, &W[0][AOoffset], naux, 0.0, Temp1[0], npoisson);
            C_DGEMM('N', 'T', npoispi[h], npoispi[h], npoisson, 1.0, Temp1[0], npoisson, poisU[0], npoisson, 0.0, &J[0][SOoffset], nauxpi[h]);
            Temp.reset();

        }
    }

    // Form indexing
    pivots_ = boost::shared_ptr<IntVector>(new IntVector(nauxpi.n(), nauxpi.pointer()));
    rev_pivots_ = boost::shared_ptr<IntVector>(new IntVector(nauxpi.n(), nauxpi.pointer()));
    for (int h = 0; h < auxpet->nirrep(); h++) {
        int* piv = pivots_->pointer(h);
        int* rpiv = pivots_->pointer(h);
        for (int Q = 0; Q < nauxpi[h]; Q++) {
            piv[Q] = Q;
            rpiv[Q] = Q;
        }
    }


}
void FittingMetric::form_cholesky_inverse()
{
    is_inverted_ = true;
    algorithm_ = "CHOLESKY";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {

        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        int info = C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
        for (int A = 0; A < metric_->colspi()[h]; A++)
            for (int B = 0; B < A; B++)
                J[A][B] = 0.0;
    }
    metric_->set_name("SO Basis Fitting Inverse (Cholesky)");
}
void FittingMetric::form_QR_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "QR";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {

        if (metric_->colspi()[h] == 0) continue;

//        metric_->print();

        double** J = metric_->pointer(h);
        int n = metric_->colspi()[h];

        // Copy the J matrix to R (actually R')
        SharedMatrix R(new Matrix("R", n, n));
        double** Rp = R->pointer();
        C_DCOPY(n*(unsigned long int)n, J[0], 1, Rp[0], 1);

        // QR Decomposition
        double* tau = new double[n];

        // First, find out how much workspace to provide
        double work_size;
        C_DGEQRF(n,n,Rp[0],n,tau,&work_size, -1);

        // Now, do the QR decomposition
        int lwork = (int)work_size;
        double *work = new double[lwork];
        C_DGEQRF(n,n,Rp[0],n,tau,work,lwork);
        delete[] work;

        // Copy Jcopy to Q (actually Q')
        SharedMatrix Q(new Matrix("Q", n, n));
        double** Qp = Q->pointer();
        C_DCOPY(n*(unsigned long int)n, Rp[0], 1, Qp[0], 1);

        // Put R in the upper triangle where it belongs
        for (int i = 1; i < n; i++)
            for (int j = 0; j < i; j++) {
                Rp[j][i] = 0.0;
            }

        // First, find out how much workspace to provide
        C_DORGQR(n,n,n,Qp[0],n,tau,&work_size,-1);

        // Now, form Q
        lwork = (int)work_size;
        work = new double[lwork];
        C_DORGQR(n,n,n,Qp[0],n,tau,work,lwork);
        delete[] work;

        //Q->print();
        //R->print();

        // Find the number of significant basis functions
        int nsig = 0;
        double R_max = fabs(Rp[0][0]);
        for (int A = 0; A < n; A++) {
            if ((fabs(Rp[A][A]) / R_max) < tol)
                break;
            nsig++;
        }

        // Transform into the reduced basis
        // Just use R's memory, don't need it anymore
        C_DGEMM('N','N',nsig,n,n,1.0,Qp[0],n,J[0],n,0.0,Rp[0],n);
        C_DGEMM('N','T',nsig,nsig,n,1.0,Rp[0],n,Qp[0],n,0.0,J[0],nsig);

        // Find the Cholesky factor in the reduced basis
        C_DPOTRF('L',nsig,J[0],nsig);

        // Backsolve the triangular factor against the change of basis matrix
        C_DTRSM('L','U','N','N',nsig,n,1.0,J[0],nsig,Qp[0],n);

        // Zero out the metric
        memset(static_cast<void*>(J[0]), '\0', n*(unsigned long int)n);

        // Copy the top bit in
        C_DCOPY(n*(unsigned long int)nsig, Qp[0], 1, J[0], 1);

        delete[] tau;
    }
    metric_->set_name("SO Basis Fitting Inverse (QR)");
}
void FittingMetric::form_eig_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "EIG";

    form_fitting_metric();

    //metric_->print();

    for (int h = 0; h < metric_->nirrep(); h++) {

        if (metric_->colspi()[h] == 0) continue;

        double** J = metric_->pointer(h);
        int n = metric_->colspi()[h];

        // Copy J to W
        SharedMatrix W(new Matrix("W", n, n));
        double** Wp = W->pointer();
        C_DCOPY(n*(unsigned long int)n,J[0],1,Wp[0],1);

        double* eigval = new double[n];
        int lwork = n * 3;
        double* work = new double[lwork];
        int stat = C_DSYEV('v','u',n,Wp[0],n,eigval,work,lwork);
        delete[] work;

        SharedMatrix Jcopy(new Matrix("Jcopy", n, n));
        double** Jcopyp = Jcopy->pointer();

        C_DCOPY(n*(unsigned long int)n,Wp[0],1,Jcopyp[0],1);

        // Now form Jp^{-1/2} = U(T)*j'^{-1/2}*U,
        // where j'^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of J'
        double max_J = eigval[n-1];

        int nsig = 0;
        for (int ind=0; ind<n; ind++) {
            if (eigval[ind] / max_J < tol || eigval[ind] <= 0.0)
                eigval[ind] = 0.0;
            else {
                nsig++;
                eigval[ind] = 1.0 / sqrt(eigval[ind]);
            }
            // scale one set of eigenvectors by the diagonal elements j^{-1/2}
            C_DSCAL(n, eigval[ind], Wp[ind], 1);
        }
        delete[] eigval;

        C_DGEMM('T','N',n,n,n,1.0,Jcopyp[0],n,Wp[0],n,0.0,J[0],n);

    }
    metric_->set_name("SO Basis Fitting Inverse (Eig)");
}
void FittingMetric::form_full_eig_inverse(double tol)
{
    is_inverted_ = true;
    algorithm_ = "EIG";

    form_fitting_metric();

    //metric_->print();

    for (int h = 0; h < metric_->nirrep(); h++) {

        if (metric_->colspi()[h] == 0) continue;

        double** J = metric_->pointer(h);
        int n = metric_->colspi()[h];

        // Copy J to W
        SharedMatrix W(new Matrix("W", n, n));
        double** Wp = W->pointer();
        C_DCOPY(n*(unsigned long int)n,J[0],1,Wp[0],1);

        double* eigval = new double[n];
        int lwork = n * 3;
        double* work = new double[lwork];
        int stat = C_DSYEV('v','u',n,Wp[0],n,eigval,work,lwork);
        delete[] work;

        SharedMatrix Jcopy(new Matrix("Jcopy", n, n));
        double** Jcopyp = Jcopy->pointer();

        C_DCOPY(n*(unsigned long int)n,Wp[0],1,Jcopyp[0],1);

        // Now form Jp^{-1/2} = U(T)*j'^{-1/2}*U,
        // where j'^{-1/2} is the diagonal matrix of the inverse square roots
        // of the eigenvalues, and U is the matrix of eigenvectors of J'
        double max_J = eigval[n-1];

        int nsig = 0;
        for (int ind=0; ind<n; ind++) {
            if (eigval[ind] / max_J < tol || eigval[ind] <= 0.0)
                eigval[ind] = 0.0;
            else {
                nsig++;
                eigval[ind] = 1.0 / eigval[ind];
            }
            // scale one set of eigenvectors by the diagonal elements j^{-1/2}
            C_DSCAL(n, eigval[ind], Wp[ind], 1);
        }
        delete[] eigval;

        C_DGEMM('T','N',n,n,n,1.0,Jcopyp[0],n,Wp[0],n,0.0,J[0],n);

    }
    metric_->set_name("SO Basis Fitting Inverse (Eig)");
}
void FittingMetric::form_full_inverse()
{
    is_inverted_ = true;
    algorithm_ = "FULL";

    form_fitting_metric();

    pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {

        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        int info = C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
        // Inverse
        info = C_DPOTRI('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);

        for (int A = 0; A < metric_->colspi()[h]; A++)
            for (int B = 0; B < A; B++)
                J[A][B] = J[B][A];
    }
    metric_->set_name("SO Basis Fitting Inverse (Full)");
}
void FittingMetric::form_cholesky_factor()
{
    is_inverted_ = true;
    algorithm_ = "CHOLESKY";

    form_fitting_metric();

    //pivot();
    for (int h = 0; h < metric_->nirrep(); h++) {

        if (metric_->colspi()[h] == 0) continue;

        // Cholesky Decomposition
        double** J = metric_->pointer(h);
        int info = C_DPOTRF('L', metric_->colspi()[h], J[0], metric_->colspi()[h]);
    }
    metric_->set_name("SO Basis Cholesky Factor (Full)");
}
void FittingMetric::pivot()
{
    for (int h = 0; h < metric_->nirrep(); h++) {

        if (metric_->colspi()[h] == 0) continue;

        double** J = metric_->pointer(h);
        int* P = pivots_->pointer(h);
        int norbs = metric_->colspi()[h];
        double* Temp = new double[norbs];

        // Pivot
        double max;
        int Temp_p;
        int pivot;
        for (int i = 0; i<norbs-1; i++) {
            max = 0.0;
            //Where's the pivot diagonal?
            for (int j = i; j<norbs; j++)
                if (max <= fabs(J[j][j])) {
                    max = fabs(J[j][j]);
                    pivot = j;
                }

            //Rows
            C_DCOPY(norbs,&J[pivot][0],1,Temp,1);
            C_DCOPY(norbs,&J[i][0],1,&J[pivot][0],1);
            C_DCOPY(norbs,Temp,1,&J[i][0],1);

            //Columns
            C_DCOPY(norbs,&J[0][pivot],norbs,Temp,1);
            C_DCOPY(norbs,&J[0][i],norbs,&J[0][pivot],norbs);
            C_DCOPY(norbs,Temp,1,&J[0][i],norbs);

            Temp_p = P[i];
            P[i] = P[pivot];
            P[pivot] = Temp_p;
        }
        delete[] Temp;

        int* R = rev_pivots_->pointer(h);
        for (int i = 0; i < norbs; i++)
            R[P[i]] = i;
    }
}

}
