#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>

#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <psifiles.h>
#include <psiconfig.h>

#include "hf.h"
#include "rhf.h"
#include "uhf.h"
#include "rohf.h"
#include "df.h"
#include <lib3index/3index.h>

//MKL Header
#include <psiconfig.h>
#if HAVE_MKL
#include <mkl.h>
#endif

//OpenMP Header
//_OPENMP is defined by the compiler if it exists
#if _OPENMP
#include <omp.h>
#endif

#include <libmints/mints.h>

using namespace boost;
using namespace std;
using namespace psi;

namespace psi { namespace scf {

DFHFLocalizer::DFHFLocalizer(int nocc, boost::shared_ptr<BasisSet> primary, boost::shared_ptr<BasisSet> auxiliary) :
    nocc_(nocc), primary_(primary), auxiliary_(auxiliary)
{
    common_init();
}
DFHFLocalizer::~DFHFLocalizer()
{
    free_int_matrix(domains_);
    free_int_matrix(old_domains_);
}
void DFHFLocalizer::common_init()
{
    debug_ = false;

    // Defaults (tight)
    schwarz_ = boost::shared_ptr<SchwarzSieve>(new SchwarzSieve(primary_, 0.0));
    charge_cutoff_ = 0.05;
    charge_leak_cutoff_ = 0.01;
    R_ext_ = 5.0;
    l_drop_ = 0;
    alpha_ = 0.0;

    int n = primary_->nbf();

    zero_ = BasisSet::zero_ao_basis_set();

    factory_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(primary_, primary_, zero_, zero_));
    auxfactory_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(auxiliary_, zero_, primary_, primary_));
    Jfactory_ = boost::shared_ptr<IntegralFactory>(new IntegralFactory(auxiliary_, zero_, auxiliary_, zero_));

    C_ = SharedMatrix(new Matrix("C Occupied (Localized AO Basis)", n, nocc_));

    Sp12_ = SharedMatrix(new Matrix("S^+1/2 (AO Basis)", n, n));
    boost::shared_ptr<OneBodyAOInt> overlap(factory_->ao_overlap());
    overlap->compute(Sp12_);

    double** Sp = Sp12_->pointer();
    SharedMatrix W(new Matrix("W", n, n));
    double** Wp = W->pointer();
    C_DCOPY(n*(unsigned long int)n,Sp[0],1,Wp[0],1);

    double* eigval = new double[n];
    int lwork = n * 3;
    double* work = new double[lwork];
    int stat = C_DSYEV('v','u',n,Wp[0],n,eigval,work,lwork);
    delete[] work;

    SharedMatrix Jcopy(new Matrix("Jcopy", n, n));
    double** Jcopyp = Jcopy->pointer();

    C_DCOPY(n*(unsigned long int)n,Wp[0],1,Jcopyp[0],1);

    double max_J = eigval[n-1];
    int nsig = 0;
    for (int ind=0; ind<n; ind++) {
        eigval[ind] = sqrt(eigval[ind]);
        C_DSCAL(n, eigval[ind], Wp[ind], 1);
    }
    delete[] eigval;

    C_DGEMM('T','N',n,n,n,1.0,Jcopyp[0],n,Wp[0],n,0.0,Sp[0],n);

    Itemp_ = SharedMatrix(new Matrix("Lowdin Charges (basis functions x orbitals)", n, nocc_));
    I_ = SharedMatrix(new Matrix("Lowdin Charges (atoms x orbitals)", primary_->molecule()->natom(), nocc_));
    domains_ = init_int_matrix(primary_->molecule()->natom(), nocc_);
    old_domains_ = init_int_matrix(primary_->molecule()->natom(), nocc_);
    domains_changed_.resize(nocc_);
    globalChange();

    Jint_ = boost::shared_ptr<TwoBodyAOInt> (Jfactory_->eri());
    Amnint_ = boost::shared_ptr<TwoBodyAOInt> (auxfactory_->eri());

}
void DFHFLocalizer::set_schwarz_cutoff(double d)
{
    schwarz_->form_schwarz_sieve(d);
    globalChange();
}
void DFHFLocalizer::globalChange()
{
    for (int i = 0; i < nocc_; i++)
        domains_changed_[i] = true;
    for (int s = 0; s < superdomains_changed_.size(); s++)
        superdomains_changed_[s] = true;
}
void DFHFLocalizer::globalReset()
{
    for (int i = 0; i < nocc_; i++)
        domains_changed_[i] = false;
    for (int s = 0; s < superdomains_changed_.size(); s++)
        superdomains_changed_[s] = false;
}
void DFHFLocalizer::choleskyLocalize(SharedMatrix D)
{
    if (debug_)
        D->print();

    SharedMatrix Dc = D->partial_cholesky_factorize(0.0,false);
    double** Dcp = Dc->pointer();

    int rank = Dc->colspi()[0];
    if (rank > nocc_)
        rank = nocc_;

    if (debug_)
        Dc->print();

    // Copy in
    C_->zero();
    double** Cp = C_->pointer();
    for (int m = 0; m < D->colspi()[0]; m++) {
        C_DCOPY(rank, Dcp[m], 1, Cp[m], 1);
    }

    if (debug_)
        C_->print();

    // Domain identification
    lowdinDomains();
    delocalDomains();
    superDomains();
    basisDomains();
    checkForChanges();
}
void DFHFLocalizer::lowdinDomains()
{
    double cutoff = charge_cutoff_;
    double R = R_ext_;

    boost::shared_ptr<Molecule> mol = primary_->molecule();
    int n = primary_->nbf();
    int nA = mol->natom();
    int nocc = nocc_;

    double** Cp = C_->pointer();
    double** Sp = Sp12_->pointer();
    double** Itempp = Itemp_->pointer();
    double** Ip = I_->pointer();

    // ==> LOWDIN CHARGES <== //

    // S^+1/2*C
    C_DGEMM('N','N',n,nocc,n,1.0,Sp[0],n,Cp[0],nocc,0.0,Itempp[0],nocc);

    // [S^+1/2]^2
    for (ULI i = 0; i < n*(ULI)nocc; i++)
        Itempp[0][i] *= Itempp[0][i];

    if (debug_)
        Itemp_->print();

    // Assign basis functions to atoms
    I_->zero();
    for (int m = 0; m < n; m++) {
        C_DAXPY(nocc, 1.0, Itempp[m], 1, Ip[primary_->shell(primary_->function_to_shell(m)).ncenter()],1);
    }

    if (debug_)
        I_->print();

    // ==> ATOM DOMAINS <== //

    // Primary domain selection
    // Save the old one first
    memcpy(static_cast<void*>(old_domains_[0]), static_cast<void*>(domains_[0]), sizeof(int)*nocc*mol->natom());
    memset(static_cast<void*>(domains_[0]), '\0', sizeof(int)*nocc*mol->natom());
    for (int m = 0; m < nA; m++) {
        for (int i = 0; i < nocc; i++) {
            if (Ip[m][i] >= cutoff)
                domains_[m][i] = 1;
        }
    }

    // Secondary domain selection
    for (int m = 0; m < nA; m++) {
        for (int i = 0; i < nocc; i++) {
            if (domains_[m][i])
                continue;
            for (int p = 0; p < nA; p++) {
                if ((domains_[p][i] == 1) && (mol->xyz(m).distance(mol->xyz(p)) <= R)) {
                    domains_[m][i] = 2;
                    continue;
                }
            }
        }
    }

    // Leak check selection
    for (int i = 0; i < nocc; i++) {
        double charge = 0.0;
        for (int A = 0; A < nA; A++) {
            if (domains_[A][i])
                charge += Ip[A][i];
        }
        charge = 1.0 - charge;
        if (charge > charge_leak_cutoff_) {
            for (int A = 0; A < nA; A++) {
                if (!domains_[A][i])
                    domains_[A][i] = 5;
            }
        }
    }
}
void DFHFLocalizer::delocalDomains()
{
    int nA = primary_->molecule()->natom();
    int nmax = (int)round(alpha_*nA);

    local_set_.clear();
    delocal_set_.clear();

    for (int i = 0; i < nocc_; i++) {
        int count = 0;
        for (int A = 0; A < nA; A++) {
            if (domains_[A][i])
                count++;
        }
        if (count > nmax) {
            delocal_set_.insert(i);
            for (int A = 0; A < nA; A++) {
                if (!domains_[A][i])
                    domains_[A][i] = 4;
            }
        }
        else {
            local_set_.insert(i);
        }
    }

    if (debug_) {
        fprintf(outfile," Domains (atoms x orbitals):\n");
        print_int_mat(domains_, nA, nocc_, outfile);
        fprintf(outfile," Old Domains (atoms x orbitals):\n");
        print_int_mat(old_domains_, nA, nocc_, outfile);
        fprintf(outfile, "\n");
        fprintf(outfile, " Localized Orbitals:\n");
        for (std::set<int>::iterator it = local_set_.begin(); it != local_set_.end(); it++) {
            fprintf(outfile, "  %4d\n", (*it));
        }
        fprintf(outfile, " Delocalized Orbitals:\n");
        for (std::set<int>::iterator it = delocal_set_.begin(); it != delocal_set_.end(); it++) {
            fprintf(outfile, "  %4d\n", (*it));
        }
    }
}
void DFHFLocalizer::superDomains()
{
    // For now every domain is a superdomain
    // TODO: intelligently unify the domains
    superdomains_.clear();
    superdomains_.resize(local_set_.size());
    for (int i = 0; i < nocc_; i++) {
        if (local_set_.count(i) != 0)
            superdomains_[i].push_back(i);
    }
}
void DFHFLocalizer::basisDomains()
{
    int nocc = nocc_;

    // ==> BASIS DOMAINS <== //

    primary_shells_.clear();
    primary_funs_.clear();
    primary_starts_.clear();
    primary_shells_.resize(nocc);
    primary_funs_.resize(nocc);
    primary_starts_.resize(nocc);

    for (int i = 0; i < nocc; i++) {
        primary_shells_[i].clear();
        primary_funs_[i].clear();
        primary_starts_[i].clear();
        for (int P = 0; P < primary_->nshell(); P++) {
            int center = primary_->shell_to_center(P);
            if (domains_[center][i]) {
                primary_shells_[i].push_back(P);
                int start = primary_->shell(P).function_index();
                primary_starts_[i].push_back(primary_funs_[i].size());
                for (int p = 0; p < primary_->shell(P).nfunction(); p++)
                    primary_funs_[i].push_back(start + p);
            }
        }
    }

    if (debug_) {
        for (int i = 0; i < nocc; i++) {
            fprintf(outfile, " Primary shells/starts for orbital %d\n", i);
            for (int P = 0; P < primary_shells_[i].size(); P++) {
                fprintf(outfile, "   %4d: %4d / %4d\n", P, primary_shells_[i][P], primary_starts_[i][P]);
            }
            fprintf(outfile, " Primary functions for orbital %d\n", i);
            for (int P = 0; P < primary_funs_[i].size(); P++) {
                fprintf(outfile, "   %4d: %4d \n", P, primary_funs_[i][P]);
            }
        }
    }

    // Highest am to keep for this iteration
    int max_am = auxiliary_->max_am() - l_drop_;

    auxiliary_shells_.clear();
    auxiliary_funs_.clear();
    auxiliary_starts_.clear();
    auxiliary_shells_.resize(nocc);
    auxiliary_funs_.resize(nocc);
    auxiliary_starts_.resize(nocc);

    for (int i = 0; i < nocc; i++) {
        auxiliary_shells_[i].clear();
        auxiliary_funs_[i].clear();
        auxiliary_starts_[i].clear();
        for (int P = 0; P < auxiliary_->nshell(); P++) {

            if (auxiliary_->shell(P).am() > max_am) continue;

            int center = auxiliary_->shell_to_center(P);
            if (domains_[center][i]) {
                auxiliary_shells_[i].push_back(P);
                int start = auxiliary_->shell(P).function_index();
                auxiliary_starts_[i].push_back(auxiliary_funs_[i].size());
                for (int p = 0; p < auxiliary_->shell(P).nfunction(); p++)
                    auxiliary_funs_[i].push_back(start + p);
            }
        }
    }

    if (debug_) {
        for (int i = 0; i < nocc; i++) {
            fprintf(outfile, " Auxiliary shells/starts for orbital %d\n", i);
            for (int P = 0; P < auxiliary_shells_[i].size(); P++) {
                fprintf(outfile, "   %4d: %4d / %4d\n", P, auxiliary_shells_[i][P], auxiliary_starts_[i][P]);
            }
            fprintf(outfile, " Auxiliary functions for orbital %d\n", i);
            for (int P = 0; P < auxiliary_funs_[i].size(); P++) {
                fprintf(outfile, "   %4d: %4d \n", P, auxiliary_funs_[i][P]);
            }
        }
    }
}
void DFHFLocalizer::checkForChanges()
{
    for (int i = 0; i < nocc_; i++) {
        if (domains_changed_[i]) continue;
        for (int A = 0; A < primary_->molecule()->natom(); A++) {
            if ((domains_[A][i] >  0 && old_domains_[A][i] == 0) ||
                (domains_[A][i] == 0 && old_domains_[A][i] >  0)) {
                domains_changed_[i] = true;
                break;
            }
        }
    }
    superdomains_changed_.resize(superdomains_.size());
    for (int s = 0; s < superdomains_.size(); s++) {
        superdomains_changed_[s] = false;
        for (int q = 0; q < superdomains_[s].size(); q++) {
            if (domains_changed_[superdomains_[s][q]]) {
                superdomains_changed_[s] = true;
                break;
            }
        }
    }
    if (debug_) {
        fprintf(outfile, "Domain status:\n");
        for (int i = 0; i < nocc_; i++)
            fprintf(outfile, "  %4d: %s\n", i, domains_changed_[i] ? "CHANGED" : "SAME");
        fprintf(outfile, "Superdomain status:\n");
        for (int i = 0; i < superdomains_.size(); i++)
            fprintf(outfile, "  %4d: %s\n", i, superdomains_changed_[i] ? "CHANGED" : "SAME");
    }
}
SharedMatrix DFHFLocalizer::computeJCholesky(int i)
{
    // Sizing
    std::stringstream name;
    name << "J Cholesky Factor for Orbital ";
    name << i;
    int n = auxiliary_funs_[i].size();
    int nshell = auxiliary_shells_[i].size();
    SharedMatrix J(new Matrix(name.str().c_str(), n, n));
    double** Jp = J->pointer();
    const double* buffer = Jint_->buffer();

    // ==> J <== //

    for (int Pl = 0; Pl < nshell; Pl++) {
        int P = auxiliary_shells_[i][Pl];
        int pstart = auxiliary_starts_[i][Pl];
        int np = auxiliary_->shell(P).nfunction();
        for (int Ql = 0; Ql <= Pl; Ql++) {
            int Q = auxiliary_shells_[i][Ql];
            int qstart = auxiliary_starts_[i][Ql];
            int nq = auxiliary_->shell(Q).nfunction();

            Jint_->compute_shell(P,0,Q,0);

            for (int op = 0, index = 0; op < np; op++) {
                int p = op + pstart;
                for (int oq = 0; oq < nq; oq++, index++) {
                    int q = oq + qstart;
                    Jp[p][q] = buffer[index];
                    Jp[q][p] = buffer[index];
                }
            }
        }
    }

    if (debug_)
        J->print();

    // ==> J CHOLESKY <== //

    // Choleskify
    C_DPOTRF('U', n, Jp[0], n);

    // Zero the upper part out
    for (int mu = 0; mu < n -1; mu++)
        memset(static_cast<void*>(&Jp[mu][mu+1]), '\0', (n - mu- 1)*sizeof(double));

    if (debug_)
        J->print();

    return J;
}
SharedMatrix DFHFLocalizer::computeAmn(int i)
{
    // Sizing
    std::stringstream name;
    name << "(A|mn) for Orbital ";
    name << i;
    int naux = auxiliary_funs_[i].size();
    int nauxshell = auxiliary_shells_[i].size();
    int n = primary_funs_[i].size();
    int nshell = primary_shells_[i].size();
    SharedMatrix Amn(new Matrix(name.str().c_str(), naux, n*(ULI)n));
    double** Amnp = Amn->pointer();
    const double* buffer = Amnint_->buffer();

    // ==> (A|mn) <== //

    long int* schwarz_shells = schwarz_->get_schwarz_shells_reverse();
    for (int Pl =0; Pl < nauxshell; Pl++) {
        int P = auxiliary_shells_[i][Pl];
        int pstart = auxiliary_starts_[i][Pl];
        int np = auxiliary_->shell(P).nfunction();
        for (int MUl =0; MUl < nshell; MUl++) {
            int MU = primary_shells_[i][MUl];
            int mustart = primary_starts_[i][MUl];
            int nmu = primary_->shell(MU).nfunction();
            for (int NUl =0; NUl <= MUl; NUl++) {
                int NU = primary_shells_[i][NUl];
                int nustart = primary_starts_[i][NUl];
                int nnu = primary_->shell(NU).nfunction();

                if (schwarz_shells[MU*(MU+1)/2 + NU] < 0) continue;

                //printf("(%2d|%2d,%2d)\n", P, MU, NU);
                Amnint_->compute_shell(P,0,MU,NU);

                for (int op = 0, index = 0; op < np; op++) {
                    int p = op + pstart;
                    for (int omu = 0; omu < nmu; omu++) {
                        int mu = omu + mustart;
                        for (int onu = 0; onu < nnu; onu++, index++) {
                            int nu = onu + nustart;
                            if (p == 0) {
                                //printf(" Adding (%d|%d,%d), which has o values of (%d,%d,%d) to %4d and %4d\n", p, mu, nu, op, omu, onu, mu*n + nu, nu*n + mu);
                            }
                            Amnp[p][mu*n + nu] = buffer[index];
                            Amnp[p][nu*n + mu] = buffer[index];
                        }
                    }
                }
            }
        }
    }

    if (debug_)
        Amn->print();

    return Amn;
}
SharedMatrix DFHFLocalizer::computeK(int i, SharedMatrix J, SharedMatrix Amn)
{
    int n = primary_funs_[i].size();
    int naux = auxiliary_funs_[i].size();

    // Result (partial K)
    std::stringstream kname;
    kname << "K for Orbital ";
    kname << i;
    SharedMatrix K(new Matrix(kname.str().c_str(), n, n));
    double** Kp = K->pointer();

    // (A|mi) <- C_in (A|mn)
    // (Q|mi) will be stored here as well
    std::stringstream Aname;
    Aname << "(A|mi) for Orbital ";
    Aname << i;
    SharedMatrix Ami(new Matrix(Aname.str().c_str(), naux, n));
    double** Amip = Ami->pointer();

    // Local C matrix
    std::stringstream Cname;
    Cname << "C_mi for Orbital ";
    Cname << i;
    boost::shared_ptr<Vector> C(new Vector(Cname.str().c_str(), n));
    double* Cp = C->pointer();
    for (int m = 0; m < n; m++) {
        Cp[m] = C_->get(0, primary_funs_[i][m], i);
    }

    if (debug_)
        C->print();

    // (A|mn)
    double** Amnp = Amn->pointer();

    // Cholesky factor for J
    double** Jp = J->pointer();

    // ==> HALF TRANSFORM <== //

    C_DGEMV('N', naux*(ULI)n, n, 1.0, Amnp[0], n, Cp, 1, 0.0, Amip[0], 1);

    if (debug_)
        Ami->print();

    // ==> FITTING <== //

    C_DTRSM('L','L','N','N',naux,n,1.0,Jp[0],naux,Amip[0],n);

    if (debug_)
        Ami->print(outfile, "after Fitting");

    // ==> CONTRACTION <== //

    C_DGEMM('T','N',n,n,naux,1.0,Amip[0],n,Amip[0],n,0.0,Kp[0],n);

    if (debug_)
        K->print();

    return K;
}
void DFHFLocalizer::computeKLocal(SharedMatrix Kglobal)
{
    double** Kp = Kglobal->pointer();

    // Resize in case the number of superdomains has changed
    JCholesky_.resize(superdomains_.size());
    for (int s = 0; s < superdomains_.size(); s++) {
        // Representative orbital of the superdomain
        int index = superdomains_[s][0];

        // Recompute J cholesky for this superdomain, if needed
        if (superdomains_changed_[s])
           JCholesky_[s] = computeJCholesky(index);

        // Compute (A|mn) for this superdomain (recomputed in each iteration,
        // but only once per
        SharedMatrix Amn = computeAmn(index);

        for (int q = 0; q < superdomains_[s].size(); q++) {
            // Compute the actual local K
            SharedMatrix Klocal = computeK(superdomains_[s][q], JCholesky_[s], Amn);

            // Add the contribution in
            double** Klp = Klocal->pointer();
            int n = primary_funs_[index].size();
            for (int m = 0; m < n; m++) {
                for (int p = 0; p < n; p++) {
                    Kp[primary_funs_[index][m]][primary_funs_[index][p]] += Klp[m][p];
                }
            }
        }
    }

    if (debug_)
        Kglobal->print(outfile, "after local part");

    // All the relevant bits are computed
    globalReset();
}

}}
