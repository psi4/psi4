#include <libciomr/libciomr.h>

#include "mints.h"
#include "extern.h"

using namespace boost;
using namespace psi;

ExternalPotential::ExternalPotential() :
    print_(1), debug_(0)
{
}
ExternalPotential::~ExternalPotential()
{
}
void ExternalPotential::clear()
{
    charges_.clear();
    bases_.clear();
}
void ExternalPotential::addCharge(double Z, double x, double y, double z)
{
    charges_.push_back(make_tuple(Z,x,y,z));
}
void ExternalPotential::addBasis(boost::shared_ptr<BasisSet> basis, SharedVector coefs)
{
    bases_.push_back(std::make_pair(basis, coefs));
}
void ExternalPotential::print(FILE* out) const
{
    fprintf(out, "  ==> External Potential Field: %s <==\n\n", name_.c_str());

    // Charges
    fprintf(out, "   => Charges [a.u.] <=\n\n");
    fprintf(out, "    %8s %8s %8s %8s\n","Z","x","y","z");
    fprintf(out, "    ---------------------------------------------\n");
    for (int i = 0 ; i < charges_.size(); i++) {
        fprintf(out, "   %8.5f %8.5f %8.5f %8.5f\n",
            get<0>(charges_[i]), get<1>(charges_[i]), get<2>(charges_[i]), get<3>(charges_[i]));
    }
    fprintf(out,"\n");

    // Bases 
    fprintf(out, "   => Diffuse Bases <=\n\n");
    for (int i = 0; i < bases_.size(); i++) {
        fprintf(out, "    Molecule %d\n\n", i+1);
        bases_[i].first->molecule()->print();
        fprintf(out, "    Basis %d\n\n", i+1);
        bases_[i].first->print_by_level(out, print_);
        fprintf(out, "    Density Coefficients %d\n\n", i+1);
        bases_[i].second->print();
    } 

    fflush(out);
}
SharedMatrix ExternalPotential::computePotentialMatrix(shared_ptr<BasisSet> basis)
{
    int n = basis->nbf();
    SharedMatrix V(new Matrix("External Potential",n,n));
    boost::shared_ptr<IntegralFactory> fact(new IntegralFactory(basis,basis,basis,basis));

    // Monopoles
    SharedMatrix V_charge(new Matrix("External Potential (Charges)", n, n));

    SharedMatrix Zxyz(new Matrix("Charges (Z,x,y,z)", charges_.size(), 4));
    double** Zxyzp = Zxyz->pointer();
    for (int i = 0; i < charges_.size(); i++) {
        Zxyzp[i][0] = get<0>(charges_[i]);
        Zxyzp[i][1] = get<1>(charges_[i]);
        Zxyzp[i][2] = get<2>(charges_[i]);
        Zxyzp[i][3] = get<3>(charges_[i]);
    }

    boost::shared_ptr<PotentialInt> pot(static_cast<PotentialInt*>(fact->ao_potential()));
    pot->set_charge_field(Zxyz);
    pot->compute(V_charge);

    V->add(V_charge);
    V_charge.reset();
    pot.reset();

    // Diffuse Bases 
    for (int ind = 0; ind < bases_.size(); ind++) {

        boost::shared_ptr<BasisSet> aux = bases_[ind].first;
        SharedVector d = bases_[ind].second;

        // TODO thread this
        boost::shared_ptr<IntegralFactory> fact2(new IntegralFactory(basis,basis,aux,BasisSet::zero_ao_basis_set()));
        boost::shared_ptr<TwoBodyAOInt> eri(fact2->eri());
        const double* buffer = eri->buffer();

        double** Vp = V->pointer();
        double*  dp = d->pointer();

        for (int Q = 0; Q < aux->nshell(); Q++) {
            for (int M = 0; M < basis->nshell(); M++) {
                for (int N = 0; N < basis->nshell(); N++) {
                    int numQ = basis->shell(Q)->nfunction();
                    int numM = basis->shell(M)->nfunction();
                    int numN = basis->shell(N)->nfunction();
                    int Qstart = basis->shell(Q)->function_index();
                    int Mstart = basis->shell(M)->function_index();
                    int Nstart = basis->shell(N)->function_index();
                    
                    eri->compute_shell(M,N,Q,0);
            
                    for (int om = 0, index = 0; om < numQ; om++) {
                        for (int on = 0; on < numN; on++) {
                            for (int oq = 0; oq < numQ; oq++, index++) {
                                Vp[om + Mstart][on + Nstart] += dp[oq + Qstart] * buffer[index];
                            }
                        }
                    }
                }
            }
        } 
    }

    return V;
}


