#include <libciomr/libciomr.h>

#include "mints.h"
#include "local.h"

using namespace boost;
using namespace psi;

namespace psi {

Local::Local(boost::shared_ptr<BasisSet> basisset, boost::shared_ptr<Matrix> C_USO) :
    basisset_(basisset), C_USO_(C_USO), print_(0), debug_(0)
{
}
Local::Local(boost::shared_ptr<BasisSet> basisset, boost::shared_ptr<BasisSet> auxiliary, boost::shared_ptr<Matrix> C_USO) :
    basisset_(basisset), auxiliary_(auxiliary), C_USO_(C_USO), print_(0), debug_(0)
{
}
Local::~Local()
{
}
void Local::common_init()
{
    boost::shared_ptr<IntegralFactory> integral(new IntegralFactory(basisset_,basisset_,basisset_,basisset_));

    boost::shared_ptr<PetiteList> pet(new PetiteList(basisset_, integral));
    AO2USO_ = pet->aotoso();

    C_AO_ = C_AO();
    L_AO_ = C_AO_;

    int nso = basisset_->nbf();
    X_ = boost::shared_ptr<Matrix>(new Matrix("S^-1/2",nso,nso));
    boost::shared_ptr<OneBodyAOInt> Sint(integral->ao_overlap());
    Sint->compute(X_);
    X_->power(-1.0/2.0);
}
boost::shared_ptr<Matrix> Local::C_USO()
{
    return C_USO_;
}
boost::shared_ptr<Matrix> Local::C_AO()
{
    if (AO2USO_->nirrep() == 1)
        return C_USO_;

    int nao = AO2USO_->rowspi()[0];
    int nmo = 0;
    for (int h = 0; h < AO2USO_->nirrep(); h++) 
        nmo += C_USO_->colspi()[h];

    boost::shared_ptr<Matrix> C = boost::shared_ptr<Matrix>(new Matrix("C (C1 Symmetry)",nao,nmo));
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
boost::shared_ptr<Matrix> Local::L_AO()
{
    return L_AO_;
}
boost::shared_ptr<Matrix> Local::AO2USO()
{
    return AO2USO_;
}
void Local::localize(const std::string& algorithm, double conv)
{
    if (algorithm == "CHOLESKY")
        localize_cholesky(conv);    
    else if (algorithm == "CHOLESKY")
        localize_pipek_mezey(conv);    
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

    L_AO_ = boost::shared_ptr<Matrix>(C_AO()->clone());

    int nso = L_AO_->rowspi()[0];
    int nmo = L_AO_->colspi()[0];

    boost::shared_ptr<Matrix> D(new Matrix("D",nso,nso));
    double** Lp = L_AO_->pointer();
    double** Dp = D->pointer();
    
    C_DGEMM('N','T',nso,nso,nmo,1.0,Lp[0],nmo,Lp[0],nmo,0.0,Dp[0],nso);

    D->partial_cholesky_factorize();

    for (int i = 0; i < nmo; i++) {
        C_DCOPY(nso, &Dp[0][i], nso, &Lp[0][i], nmo);
    }
}
void Local::localize_pipek_mezey(double conv) 
{
    throw FeatureNotImplemented("psi::Local","localize_pipek_mezey",__FILE__,__LINE__);
    L_AO_ = boost::shared_ptr<Matrix>(C_AO()->clone());

}
void Local::localize_boys(double conv) 
{
    throw FeatureNotImplemented("psi::Local","localize_boys",__FILE__,__LINE__);
    L_AO_ = boost::shared_ptr<Matrix>(C_AO()->clone());

}
void Local::localize_er(double conv) 
{
    throw FeatureNotImplemented("psi::Local","localize_er",__FILE__,__LINE__);
    L_AO_ = boost::shared_ptr<Matrix>(C_AO()->clone());

}
boost::shared_ptr<Matrix> Local::lowdin_charges(boost::shared_ptr<Matrix> C)
{
    boost::shared_ptr<Molecule> molecule = basisset_->molecule(); 
    int nmo = C->colspi()[0];
    int nso = C->rowspi()[0];
    int natom = molecule->natom();
    boost::shared_ptr<Matrix> Q(new Matrix("Q: Gross Lowdin charges (nmo x natom)", nmo, natom));

    boost::shared_ptr<Matrix> XC(new Matrix("XC", nso, nmo));

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
void Local::compute_boughton_pulay_domains(double Qcutoff) 
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
        std::set<int> total_atoms;

        // TODO

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
                if (v.distance(molecule->xyz(A)) < Rext) {
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

        int mstart = basis->shell(M)->function_index();
        int nM = basis->shell(M)->nfunction();

        for (int om = 0; om < nM; om++) {
            domain->functions_local_to_global_.push_back(om + mstart);
            domain->functions_global_to_local_[om + mstart] = fun_counter++;
        } 
    }

    return domain;
}


}

