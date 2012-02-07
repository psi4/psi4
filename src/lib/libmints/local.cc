#include <utility>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

#include "mints.h"
#include "local.h"
#include <physconst.h>

using namespace boost;
using namespace psi;

namespace psi {

Local::Local(boost::shared_ptr<BasisSet> basisset, SharedMatrix C_USO) :
    basisset_(basisset), C_USO_(C_USO), print_(0), debug_(0)
{
    common_init();
}
Local::Local(boost::shared_ptr<BasisSet> basisset, boost::shared_ptr<BasisSet> auxiliary, SharedMatrix C_USO) :
    basisset_(basisset), auxiliary_(auxiliary), C_USO_(C_USO), print_(0), debug_(0)
{
    common_init();
}
Local::~Local()
{
}
void Local::common_init()
{
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));

    if (C_USO_->nirrep() > 1) {
        boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral));
        AO2USO_ = pet->aotoso();
    }

    C_AO_ = C_AO();
    L_AO_ = C_AO_;

    int nso = basisset_->nbf();
    S_ = SharedMatrix(new Matrix("S",nso,nso));
    X_ = SharedMatrix(new Matrix("S^+1/2",nso,nso));
    boost::shared_ptr<OneBodyAOInt> Sint(integral->ao_overlap());
    Sint->compute(S_);
    X_->copy(S_);
    X_->power(+1.0/2.0);
}
void Local::print(FILE* out)
{
    fprintf(out, "  ==> Localization <==\n\n");

    basisset_->print_by_level(out,3);
    if (auxiliary_.get())
        auxiliary_->print_by_level(out,3);

    C_USO_->print(out);
    C_AO_->print(out);
    L_AO_->print(out);

    fprintf(out, "  => Metrics <=\n\n");
    fprintf(out, "  Boys:                %11.3E\n", boys_metric());
    fprintf(out, "  Edmiston-Ruedenberg: %11.3E\n", er_metric());
    fprintf(out, "  Pipek-Mezey:         %11.3E\n\n", pm_metric());

    if (Q_.get())
        Q_->print(out);

    if (domains_.size()) {
        fprintf(out, "  => Primary Domains <=\n\n");
        for (int i = 0; i < L_AO_->colspi()[0]; i++) {
            domains_[i]->print(out, i);
        }
    }

    if (auxiliary_domains_.size()) {
        fprintf(out, "  => Auxiliary Domains <=\n\n");
        for (int i = 0; i < L_AO_->colspi()[0]; i++) {
            auxiliary_domains_[i]->print(out, i);
        }
    }

    fflush(out);
}
double Local::er_metric()
{
    int nso = L_AO_->rowspi()[0];
    int nmo = L_AO_->colspi()[0];

    boost::shared_ptr<Vector> iiii(new Vector("(ii|ii)", nmo));
    double* ip = iiii->pointer();

    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<TwoBodyAOInt> eri(integral->eri());
    const double* buffer = eri->buffer();
    double** Lp = L_AO_->pointer();

    for (int M = 0; M < basisset_->nshell(); M++) {
        for (int N = 0; N < basisset_->nshell(); N++) {
            for (int R = 0; R < basisset_->nshell(); R++) {
                for (int S = 0; S < basisset_->nshell(); S++) {

                int nM = basisset_->shell(M).nfunction();
                int nN = basisset_->shell(N).nfunction();
                int nR = basisset_->shell(R).nfunction();
                int nS = basisset_->shell(S).nfunction();
                int mstart = basisset_->shell(M).function_index();
                int nstart = basisset_->shell(N).function_index();
                int rstart = basisset_->shell(R).function_index();
                int sstart = basisset_->shell(S).function_index();

                eri->compute_shell(M,N,R,S);

                for (int oM = 0, index = 0; oM < nM; oM++) {
                    for (int oN = 0; oN < nN; oN++) {
                        for (int oR = 0; oR < nR; oR++) {
                            for (int oS = 0; oS < nS; oS++, index++) {
                                for (int i = 0; i < nmo; i++) {
                                    ip[i] += Lp[oM + mstart][i] * Lp[oN + nstart][i] *
                                             buffer[index] *
                                             Lp[oR + rstart][i] * Lp[oS + sstart][i];
                                }
                }}}}
    }}}}

    double metric = 0.0;

    for (int i = 0; i < nmo; i++) {
        metric += ip[i];
    }

    return metric;
}
double Local::pm_metric()
{
    int nso = L_AO_->rowspi()[0];
    int nmo = L_AO_->colspi()[0];
    int natom = basisset_->molecule()->natom();

    SharedMatrix Q = mulliken_charges(L_AO_);

    double metric = C_DDOT(nmo * (ULI) natom, Q->pointer()[0], 1, Q->pointer()[0], 1);

    return metric;
}
double Local::boys_metric()
{
    double metric = 0.0;

    int nso = L_AO_->rowspi()[0];
    int nmo = L_AO_->colspi()[0];

    // Build dipole integrals
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));
    boost::shared_ptr<OneBodyAOInt> Dint(integral->ao_dipole());

    std::vector<SharedMatrix > dipole;
    dipole.push_back(SharedMatrix(new Matrix("Dipole X", nso, nso)));
    dipole.push_back(SharedMatrix(new Matrix("Dipole Y", nso, nso)));
    dipole.push_back(SharedMatrix(new Matrix("Dipole Z", nso, nso)));

    Dint->compute(dipole);

    SharedMatrix XC(new Matrix("XC", nso, nmo));

    // X contribution
    C_DGEMM('N','N',nso,nmo,nso,1.0,dipole[0]->pointer()[0],nso,L_AO_->pointer()[0],nmo,0.0,XC->pointer()[0],nmo);

    for (int i = 0; i < nmo; i++) {
        double cont = C_DDOT(nso,&L_AO_->pointer()[0][i],nmo,&dipole[0]->pointer()[0][i],nmo);
        metric += cont * cont;
    }

    // Y contribution
    C_DGEMM('N','N',nso,nmo,nso,1.0,dipole[1]->pointer()[0],nso,L_AO_->pointer()[0],nmo,0.0,XC->pointer()[0],nmo);

    for (int i = 0; i < nmo; i++) {
        double cont = C_DDOT(nso,&L_AO_->pointer()[0][i],nmo,&dipole[1]->pointer()[0][i],nmo);
        metric += cont * cont;
    }

    // Z contribution
    C_DGEMM('N','N',nso,nmo,nso,1.0,dipole[2]->pointer()[0],nso,L_AO_->pointer()[0],nmo,0.0,XC->pointer()[0],nmo);

    for (int i = 0; i < nmo; i++) {
        double cont = C_DDOT(nso,&L_AO_->pointer()[0][i],nmo,&dipole[2]->pointer()[0][i],nmo);
        metric += cont * cont;
    }

    return metric;
}
SharedMatrix Local::C_USO()
{
    return C_USO_;
}
SharedMatrix Local::C_AO()
{
    if (!AO2USO_.get() || AO2USO_->nirrep() == 1)
        return C_USO_;

    int nao = AO2USO_->rowspi()[0];
    int nmo = 0;
    for (int h = 0; h < AO2USO_->nirrep(); h++)
        nmo += C_USO_->colspi()[h];

    SharedMatrix C = SharedMatrix(new Matrix("C (C1 Symmetry)",nao,nmo));
    double** Cp = C->pointer();

    int counter = 0;
    for (int h = 0; h < AO2USO_->nirrep(); h++) {
        int nsopi = AO2USO_->colspi()[h];
        int nmopi = C->colspi()[h];
        if (nsopi == 0 || nmopi == 0) continue;
        double** Ca = C_USO_->pointer(h);
        double** X = AO2USO_->pointer(h);

        C_DGEMM('N','N',nao,nmopi,nsopi,1.0,X[0],nsopi,Ca[0],nmopi,0.0,&Cp[0][counter],nmo);

        counter += nmopi;
    }
    return C;
}
SharedMatrix Local::L_AO()
{
    return L_AO_;
}
SharedMatrix Local::AO2USO()
{
    return AO2USO_;
}
void Local::localize(const std::string& algorithm, double conv)
{
    if (algorithm == "CHOLESKY")
        localize_cholesky(conv);
    else if (algorithm == "PM")
        localize_pm(conv);
    else if (algorithm == "BOYS")
        localize_boys(conv);
    else if (algorithm == "ER")
        localize_er(conv);
    else
        throw PSIEXCEPTION("Localization algorithm not recognized");
}
void Local::localize_cholesky(double conv)
{
    if (print_) {
        fprintf(outfile, "  ==> Localization: Cholesky <==\n\n");
    }

    L_AO_ = SharedMatrix(C_AO()->clone());
    L_AO_->set_name("L Cholesky (C1 Symmetry)");

    int nso = L_AO_->rowspi()[0];
    int nmo = L_AO_->colspi()[0];

    SharedMatrix D(new Matrix("D",nso,nso));
    double** Lp = L_AO_->pointer();
    double** Dp = D->pointer();

    C_DGEMM('N','T',nso,nso,nmo,1.0,Lp[0],nmo,Lp[0],nmo,0.0,Dp[0],nso);

    if (debug_)
        D->print();

    SharedMatrix L = D->partial_cholesky_factorize();

    if (debug_)
        L->print();

    int nmo2 = L->colspi()[0];
    double** Lp2 = L->pointer();

    if (nmo2 < nmo)
        throw PSIEXCEPTION("Local: Cholesky factor has smaller numerical rank than nmo!");

    for (int i = 0; i < nmo; i++) {
        C_DCOPY(nso, &Lp2[0][i], nmo2, &Lp[0][i], nmo);
    }

    if (debug_) {
        C_DGEMM('N','T',nso,nso,nmo,-1.0,Lp[0],nmo,Lp[0],nmo,1.0,Dp[0],nso);
        D->print(outfile, "Residual");
    }
}
void Local::localize_pm(double conv)
{
    throw FeatureNotImplemented("psi::Local","localize_pm",__FILE__,__LINE__);
    L_AO_ = SharedMatrix(C_AO()->clone());
    L_AO_->set_name("L Pipek-Mezey (C1 Symmetry)");

}
void Local::localize_boys(double conv)
{
    throw FeatureNotImplemented("psi::Local","localize_boys",__FILE__,__LINE__);
    L_AO_ = SharedMatrix(C_AO()->clone());
    L_AO_->set_name("L Boys (C1 Symmetry)");

}
void Local::localize_er(double conv)
{
    throw FeatureNotImplemented("psi::Local","localize_er",__FILE__,__LINE__);
    L_AO_ = SharedMatrix(C_AO()->clone());
    L_AO_->set_name("L Edmiston-Ruedenberg (C1 Symmetry)");

}
SharedMatrix Local::lowdin_charges(SharedMatrix C)
{
    boost::shared_ptr<Molecule> molecule = basisset_->molecule();
    int nmo = C->colspi()[0];
    int nso = C->rowspi()[0];
    int natom = molecule->natom();
    SharedMatrix Q(new Matrix("Q: Gross Lowdin charges (nmo x natom)", nmo, natom));

    SharedMatrix XC(new Matrix("XC", nso, nmo));

    double** Cp  = C->pointer();
    double** Xp  = X_->pointer();
    double** XCp = XC->pointer();

    C_DGEMM('N','N',nso,nmo,nso,1.0,Xp[0],nso,Cp[0],nmo,0.0,XCp[0],nmo);

    double** Qp = Q->pointer();
    for (int i = 0; i < nmo; i++) {
        for (int m = 0; m < nso; m++) {
            int atom = basisset_->shell_to_center(basisset_->function_to_shell(m));
            Qp[i][atom] += XCp[m][i] * XCp[m][i];
        }
    }

    return Q;
}
SharedMatrix Local::mulliken_charges(SharedMatrix C)
{
    boost::shared_ptr<Molecule> molecule = basisset_->molecule();
    int nmo = C->colspi()[0];
    int nso = C->rowspi()[0];
    int natom = molecule->natom();
    SharedMatrix Q(new Matrix("Q: Gross Mulliken charges (nmo x natom)", nmo, natom));

    SharedMatrix XC(new Matrix("XC", nso, nmo));

    double** Cp  = C->pointer();
    double** Sp  = S_->pointer();
    double** XCp = XC->pointer();

    C_DGEMM('N','N',nso,nmo,nso,1.0,Sp[0],nso,Cp[0],nmo,0.0,XCp[0],nmo);

    double** Qp = Q->pointer();
    for (int i = 0; i < nmo; i++) {
        for (int m = 0; m < nso; m++) {
            int atom = basisset_->shell_to_center(basisset_->function_to_shell(m));
            Qp[i][atom] += XCp[m][i] * Cp[m][i];
        }
    }

    return Q;
}
void Local::compute_boughton_pulay_domains(double Qcutoff)
{
    boost::shared_ptr<Molecule> molecule = basisset_->molecule();
    int nmo = L_AO_->colspi()[0];
    int nso = L_AO_->rowspi()[0];
    int natom = molecule->natom();

    // Build Mulliken charges for current local coefficients
    Q_ = mulliken_charges(L_AO_);
    double** Qp = Q_->pointer();

    // Clear domains
    domains_.clear();
    auxiliary_domains_.clear();

    double** Sp = S_->pointer();
    double** Lp = L_AO_->pointer();

    // Premultiply SL
    SharedMatrix SL(new Matrix("SL", nso, nmo));
    double** SLp = SL->pointer();
    C_DGEMM('N','N',nso,nmo,nso,1.0,Sp[0],nso,Lp[0],nmo,0.0,SLp[0],nmo);

    // Build each domain
    for (int i = 0; i < nmo; i++) {
        std::set<int> total_atoms;

        // Rank atoms in this domain by Lowdin charges
        std::vector<std::pair<double, int> > charges;
        for (int A = 0; A < natom; A++) {
            charges.push_back(std::make_pair(Qp[i][A], A));
        }
        std::sort(charges.begin(), charges.end(), std::greater<std::pair<double, int> >());

        // Add atoms until domain becomes complete enough
        std::vector<int> atoms_in_domain;
        std::vector<int> funs_in_domain;
        for (int atom = 0; atom < natom; atom++) {
            // Add the current highest-charge atom to the domain
            int current_atom = charges[atom].second;
            atoms_in_domain.push_back(current_atom);

            // Add the functions from that atom to the funs_in_domain list
            for (int M = 0; M < basisset_->nshell(); M++) {
                if (basisset_->shell_to_center(M) != current_atom) continue;
                int nM = basisset_->shell(M).nfunction();
                int mstart = basisset_->shell(M).function_index();
                for (int om = 0; om < nM; om++) {
                    funs_in_domain.push_back(om + mstart);
                }
            }

            // Figure out how many functions there are
            int nfun = funs_in_domain.size();

            // Temps
            SharedMatrix Smn(new Matrix("Smn", nfun, nfun));
            boost::shared_ptr<Vector> A(new Vector("A", nfun));
            double** Smnp = Smn->pointer();
            double* Ap = A->pointer();

            // Place the proper overlap elements
            for (int m = 0; m < nfun; m++) {
                for (int n = 0; n < nfun; n++) {
                    Smnp[m][n] = Sp[funs_in_domain[m]][funs_in_domain[n]];
                }
            }

            // Place the proper SL elements
            for (int m = 0; m < nfun; m++) {
                Ap[m] = SLp[funs_in_domain[m]][i];
            }

            // Find the A vector
            int info1 = C_DPOTRF('L', nfun, Smnp[0], nfun);
            if (info1 != 0) {
                throw PSIEXCEPTION("Local: Boughton Pulay Domains: Cholesky Failed!");
            }
            int info2 = C_DPOTRS('L', nfun, 1, Smnp[0], nfun, Ap, nfun);
            if (info2 != 0) {
                throw PSIEXCEPTION("Local: Boughton Pulay Domains: Cholesky Solve Failed!");
            }

            // Compute the incompleteness metric
            double incompleteness = 1.0 - C_DDOT(nfun, Ap, 1, Ap, 1);
            if (incompleteness < 1.0 - Qcutoff)
                break;
        }

        // rebuild a set<int>, which OrbtialDomain expects
        for (int A = 0; A < atoms_in_domain.size(); A++) {
            total_atoms.insert(atoms_in_domain[A]);
        }

        domains_.push_back(OrbitalDomain::buildOrbitalDomain(basisset_, total_atoms));
        if (auxiliary_.get() != NULL)
            auxiliary_domains_.push_back(OrbitalDomain::buildOrbitalDomain(auxiliary_, total_atoms));
    }
}
void Local::compute_polly_domains(double Qcutoff, double Rext, double Qcheck)
{
    boost::shared_ptr<Molecule> molecule = basisset_->molecule();
    int nmo = L_AO_->colspi()[0];
    int nso = L_AO_->rowspi()[0];
    int natom = molecule->natom();

    // Build Lowdin charges for current local coefficients
    Q_ = lowdin_charges(L_AO_);
    double** Qp = Q_->pointer();

    // Clear domains
    domains_.clear();
    auxiliary_domains_.clear();

    // Build each domain
    for (int i = 0; i < nmo; i++) {
        std::set<int> charge_atoms;
        std::set<int> range_atoms;
        std::set<int> total_atoms;

        // Add atoms if they are greater than the charge cutoff
        for (int A = 0; A < natom; A++) {
            if (Qp[i][A] > Qcutoff) {
                charge_atoms.insert(A);
            }
        }

        // Add atoms if they are inside the extended range parameter
        for (std::set<int>::const_iterator it = charge_atoms.begin();
            it != charge_atoms.end(); it++) {
            Vector3 v = molecule->xyz(*it);
            for (int A = 0; A < natom; A++) {
                if (v.distance(molecule->xyz(A)) < Rext / _bohr2angstroms) {
                    range_atoms.insert(A);
                }
            }
        }

        // Merge the two lists
        for (std::set<int>::const_iterator it = charge_atoms.begin();
            it != charge_atoms.end(); it++) {
            total_atoms.insert((*it));
        }
        for (std::set<int>::const_iterator it = range_atoms.begin();
            it != range_atoms.end(); it++) {
            total_atoms.insert((*it));
        }

        // Check total charge to be sure the domain did not sneak by
        double charge_check = 0.0;
        for (std::set<int>::const_iterator it = total_atoms.begin();
            it != total_atoms.end(); it++) {
            charge_check += Qp[i][(*it)];
        }

        // If the domain sneaks by, add all atoms to the domain
        if (charge_check < Qcheck) {
            total_atoms.clear();
            for (int A = 0; A < natom; A++) {
                total_atoms.insert(A);
            }
        }

        domains_.push_back(OrbitalDomain::buildOrbitalDomain(basisset_, total_atoms));
        if (auxiliary_.get() != NULL)
            auxiliary_domains_.push_back(OrbitalDomain::buildOrbitalDomain(auxiliary_, total_atoms));
    }
}
OrbitalDomain::OrbitalDomain()
{
}
OrbitalDomain::~OrbitalDomain()
{
}
void OrbitalDomain::print(FILE* out, int label)
{
    fprintf(out, "  => Orbital Domain %d: %ld Atoms <=\n\n", label, atoms_.size());

    fprintf(out, "    Atoms Local to Global:\n");
    for (int A = 0; A < atoms_local_to_global_.size(); A++) {
        fprintf(out, "    %4d -> %4d\n", A, atoms_local_to_global_[A]);
    }
    fprintf(out, "    \n");

    fprintf(out, "    Atoms Global to Local:\n");
    for (int A = 0; A < atoms_global_to_local_.size(); A++) {
        int a = atoms_local_to_global_[A];
        fprintf(out, "    %4d -> %4d\n", a, atoms_global_to_local_[a]);
    }
    fprintf(out, "    \n");

    fprintf(out, "    Shells Local to Global:\n");
    for (int A = 0; A < shells_local_to_global_.size(); A++) {
        fprintf(out, "    %4d -> %4d\n", A, shells_local_to_global_[A]);
    }
    fprintf(out, "    \n");

    fprintf(out, "    Shells Global to Local:\n");
    for (int A = 0; A < shells_global_to_local_.size(); A++) {
        int a = shells_local_to_global_[A];
        fprintf(out, "    %4d -> %4d\n", a, shells_global_to_local_[a]);
    }
    fprintf(out, "    \n");

    fprintf(out, "    Functions Local to Global:\n");
    for (int A = 0; A < functions_local_to_global_.size(); A++) {
        fprintf(out, "    %4d -> %4d\n", A, functions_local_to_global_[A]);
    }
    fprintf(out, "    \n");

    fprintf(out, "    Functions Global to Local:\n");
    for (int A = 0; A < functions_global_to_local_.size(); A++) {
        int a = functions_local_to_global_[A];
        fprintf(out, "    %4d -> %4d\n", a, functions_global_to_local_[a]);
    }
    fprintf(out, "    \n");
}
boost::shared_ptr<OrbitalDomain> OrbitalDomain::buildOrbitalDomain(boost::shared_ptr<BasisSet> basis, const std::set<int>& atoms)
{
    boost::shared_ptr<OrbitalDomain> domain(new OrbitalDomain());
    domain->atoms_ = atoms;

    boost::shared_ptr<Molecule> molecule = basis->molecule();
    int natom = molecule->natom();
    int nso = basis->nbf();
    int nshell = basis->nshell();

    for (std::set<int>::const_iterator it = atoms.begin();
        it != atoms.end(); it++) {
        domain->atoms_local_to_global_.push_back((*it));
    }

    std::sort(domain->atoms_local_to_global_.begin(), domain->atoms_local_to_global_.end());

    for (int A = 0; A < atoms.size(); A++) {
        domain->atoms_global_to_local_[domain->atoms_local_to_global_[A]] = A;
    }

    int shell_counter = 0;
    int fun_counter = 0;
    for (int M = 0; M < nshell; M++) {
        int atom = basis->shell_to_center(M);
        if (!atoms.count(atom)) continue;

        domain->shells_local_to_global_.push_back(M);
        domain->shells_global_to_local_[M] = shell_counter++;

        int mstart = basis->shell(M).function_index();
        int nM = basis->shell(M).nfunction();

        for (int om = 0; om < nM; om++) {
            domain->functions_local_to_global_.push_back(om + mstart);
            domain->functions_global_to_local_[om + mstart] = fun_counter++;
        }
    }

    return domain;
}

} // Namespace psi

