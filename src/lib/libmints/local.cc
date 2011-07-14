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
boost::shared_ptr<Matrix> Local::localize(const std::string& algorithm, double conv)
{
    if (algorithm == "CHOLESKY")
        return localize_cholesky(conv);    
    else if (algorithm == "CHOLESKY")
        return localize_pipek_mezey(conv);    
    else if (algorithm == "BOYS")
        return localize_boys(conv);    
    else if (algorithm == "ER")
        return localize_er(conv);    
    else 
        throw PSIEXCEPTION("Localization algorithm not recognized");
}
boost::shared_ptr<Matrix> Local::localize_cholesky(double conv) 
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

    return L_AO_;
}
boost::shared_ptr<Matrix> Local::localize_pipek_mezey(double conv) 
{
    throw FeatureNotImplemented("psi::Local","localize_pipek_mezey",__FILE__,__LINE__);
    L_AO_ = boost::shared_ptr<Matrix>(C_AO()->clone());

    return L_AO_;
}
boost::shared_ptr<Matrix> Local::localize_boys(double conv) 
{
    throw FeatureNotImplemented("psi::Local","localize_boys",__FILE__,__LINE__);
    L_AO_ = boost::shared_ptr<Matrix>(C_AO()->clone());

    return L_AO_;
}
boost::shared_ptr<Matrix> Local::localize_er(double conv) 
{
    throw FeatureNotImplemented("psi::Local","localize_er",__FILE__,__LINE__);
    L_AO_ = boost::shared_ptr<Matrix>(C_AO()->clone());

    return L_AO_;
}

}

