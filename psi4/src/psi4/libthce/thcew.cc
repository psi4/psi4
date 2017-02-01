/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>



#include "psi4/libmints/sieve.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "thce.h"
#include "thcew.h"
#include "laplace.h"
#include "lreri.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/matrix.h"

#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;


namespace psi {

THCEW::THCEW() :
    Wavefunction(Process::environment.options)
{
    common_init();
}
THCEW::~THCEW()
{
}
void THCEW::common_init()
{
    print_ = options_.get_int("PRINT");
    debug_ = options_.get_int("DEBUG");


    //reference_wavefunction_ = Process::environment.legacy_wavefunction();
    throw PSIEXCEPTION("Rob: I broke your code. Check your email.");
    if (!reference_wavefunction_) {
        throw PSIEXCEPTION("THCEW: Run SCF first");
    }

    if (options_.get_str("REFERENCE") == "ROHF" || options_.get_str("REFERENCE") == "CUHF")
        throw PSIEXCEPTION("Does not currently work for non RHF references. Blame DGAS.");
        // reference_wavefunction_->semicanonicalize();

    shallow_copy(reference_wavefunction_);

    thce_ = std::shared_ptr<THCE>(new THCE());
}
RTHCEW::RTHCEW() :
    THCEW()
{
    common_init();
}
RTHCEW::~RTHCEW()
{
}
void RTHCEW::common_init()
{
    int nso  = nso_;
    int nmo  = nmo_;
    int nocc = doccpi_.sum();
    int nvir = nmo - nocc;
    int nfocc = frzcpi_.sum();
    int nfvir = frzvpi_.sum();
    int naocc = nocc - nfocc;
    int navir = nvir - nfvir;
    int nact  = naocc + navir;

    thce_->new_dimension("nso",nso);
    thce_->new_dimension("nmo",nmo);
    thce_->new_dimension("nocc",nocc);
    thce_->new_dimension("nfocc",nfocc);
    thce_->new_dimension("naocc",naocc);
    thce_->new_dimension("navir",navir);
    thce_->new_dimension("nfvir",nfvir);
    thce_->new_dimension("nvir",nvir);
    thce_->new_dimension("nact",nact);

    std::shared_ptr<Matrix> C = Ca_subset("AO","ALL");
    std::shared_ptr<Tensor> Ca = CoreTensor::build("Cmo","nmo",nmo,"nso",nso);
    double* C1p = C->pointer()[0];
    double* Cap = Ca->pointer();

    for (int m = 0; m < nso; m++) {
        for (int n = 0; n < nmo; n++) {
            *(Cap + n * nso + m) = *(C1p + m * nmo + n);
        }
    }

    thce_->add_tensor("Cmo", Ca);
    thce_->add_tensor("Cocc", CoreTensor::build("Cocc" ,"nocc" ,nocc ,"nso",nso,Cap + 0              * nso,true));
    thce_->add_tensor("Cfocc",CoreTensor::build("Cfocc","nfocc",nfocc,"nso",nso,Cap + 0              * nso,true));
    thce_->add_tensor("Caocc",CoreTensor::build("Caocc","naocc",naocc,"nso",nso,Cap + nfocc          * nso,true));
    thce_->add_tensor("Cavir",CoreTensor::build("Cavir","navir",navir,"nso",nso,Cap + nocc           * nso,true));
    thce_->add_tensor("Cfvir",CoreTensor::build("Cfvir","nfvir",nfvir,"nso",nso,Cap + (nocc + navir) * nso,true));
    thce_->add_tensor("Cvir", CoreTensor::build("Cvir" ,"nvir" ,nvir ,"nso",nso,Cap + nocc           * nso,true));

    std::shared_ptr<Vector> eps = epsilon_a_subset("AO","ALL");
    std::shared_ptr<Tensor> eps_a = CoreTensor::build("eps_mo","nmo",nmo,eps->pointer(),false);
    double* eps_ap = eps_a->pointer();

    thce_->add_tensor("eps_mo", eps_a);
    thce_->add_tensor("eps_occ", CoreTensor::build("eps_occ" ,"nocc" ,nocc ,eps_ap + 0             ,true));
    thce_->add_tensor("eps_focc",CoreTensor::build("eps_focc","nfocc",nfocc,eps_ap + 0             ,true));
    thce_->add_tensor("eps_aocc",CoreTensor::build("eps_aocc","naocc",naocc,eps_ap + nfocc         ,true));
    thce_->add_tensor("eps_avir",CoreTensor::build("eps_avir","navir",navir,eps_ap + nocc          ,true));
    thce_->add_tensor("eps_fvir",CoreTensor::build("eps_fvir","nfvir",nfvir,eps_ap + (nocc + navir),true));
    thce_->add_tensor("eps_vir", CoreTensor::build("eps_vir" ,"nvir" ,nvir ,eps_ap + nocc          ,true));

}
void RTHCEW::build_laplace(double delta, double omega)
{
    std::shared_ptr<Vector> eps_aocc = epsilon_a_subset("AO","ACTIVE_OCC");
    std::shared_ptr<Vector> eps_avir = epsilon_a_subset("AO","ACTIVE_VIR");
    std::shared_ptr<LaplaceDenom> laplace(new LaplaceDenom(eps_aocc,eps_avir,delta,omega,2));
    laplace->compute("pi_i","pi_a");
    thce_->add_tensor("pi_i",laplace->tau_occ());
    thce_->add_tensor("pi_a",laplace->tau_vir());
    (*thce_)["pi_i"]->dimensions()[0] = "nw";
    (*thce_)["pi_i"]->dimensions()[1] = "naocc";
    (*thce_)["pi_a"]->dimensions()[0] = "nw";
    (*thce_)["pi_a"]->dimensions()[1] = "navir";
    thce_->new_dimension("nw",laplace->npoints());
}
void RTHCEW::build_df_ia(std::shared_ptr<BasisSet> auxiliary)
{
    std::shared_ptr<DFERI> dferi = DFERI::build(basisset_,auxiliary,options_,reference_wavefunction_);
    dferi->add_pair_space("Bia","ACTIVE_OCC","ACTIVE_VIR");
    dferi->compute();
    std::shared_ptr<Tensor> Bia = dferi->ints()["Bia"];
    dferi.reset();

    Bia->dimensions()[0] = "naocc";
    Bia->dimensions()[1] = "navir";
    Bia->dimensions()[2] = "naux";
    thce_->new_dimension("naux", Bia->sizes()[2]);
    thce_->add_tensor("Bia",Bia);
}
void RTHCEW::build_df_act(std::shared_ptr<BasisSet> auxiliary)
{
    std::shared_ptr<DFERI> dferi = DFERI::build(basisset_,auxiliary,options_,reference_wavefunction_);
    dferi->add_pair_space("Bii","ACTIVE_OCC","ACTIVE_OCC");
    dferi->add_pair_space("Bia","ACTIVE_OCC","ACTIVE_VIR");
    dferi->add_pair_space("Bai","ACTIVE_VIR","ACTIVE_OCC");
    dferi->add_pair_space("Baa","ACTIVE_VIR","ACTIVE_VIR");
    dferi->compute();
    std::shared_ptr<Tensor> Bii = dferi->ints()["Bii"];
    std::shared_ptr<Tensor> Bia = dferi->ints()["Bia"];
    std::shared_ptr<Tensor> Bai = dferi->ints()["Bai"];
    std::shared_ptr<Tensor> Baa = dferi->ints()["Baa"];
    dferi.reset();

    Bii->dimensions()[0] = "naocc";
    Bii->dimensions()[1] = "naocc";
    Bii->dimensions()[2] = "naux";
    Bia->dimensions()[0] = "naocc";
    Bia->dimensions()[1] = "navir";
    Bia->dimensions()[2] = "naux";
    Bai->dimensions()[0] = "navir";
    Bai->dimensions()[1] = "naocc";
    Bai->dimensions()[2] = "naux";
    Baa->dimensions()[0] = "navir";
    Baa->dimensions()[1] = "navir";
    Baa->dimensions()[2] = "naux";
    thce_->new_dimension("naux", Bia->sizes()[2]);
    thce_->add_tensor("Bii",Bii);
    thce_->add_tensor("Bia",Bia);
    thce_->add_tensor("Bai",Bai);
    thce_->add_tensor("Baa",Baa);
}
void RTHCEW::build_df_pp(std::shared_ptr<BasisSet> auxiliary)
{
    std::shared_ptr<DFERI> dferi = DFERI::build(basisset_,auxiliary,options_,reference_wavefunction_);
    dferi->add_pair_space("Bpp","ACTIVE_ALL","ACTIVE_ALL");
    dferi->compute();
    std::shared_ptr<Tensor> Bpp = dferi->ints()["Bpp"];
    dferi.reset();

    Bpp->dimensions()[0] = "nact";
    Bpp->dimensions()[1] = "nact";
    Bpp->dimensions()[2] = "naux";
    thce_->new_dimension("naux", Bpp->sizes()[2]);
    thce_->add_tensor("Bpp",Bpp);
}
void RTHCEW::build_lsthc_ia(std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> X)
{
    std::shared_ptr<LSTHCERI> lsthceri = LSTHCERI::build(basisset_,auxiliary,X,options_,reference_wavefunction_);
    lsthceri->add_eri_space("ovov","ACTIVE_OCC","ACTIVE_VIR","ACTIVE_OCC","ACTIVE_VIR");
    lsthceri->compute();
    std::vector<std::shared_ptr<Tensor> > ovov = lsthceri->ints()["ovov"];
    lsthceri.reset();

    std::shared_ptr<Tensor> Xi = ovov[0];
    std::shared_ptr<Tensor> Xa = ovov[1];
    std::shared_ptr<Tensor> Z  = ovov[2];
    std::shared_ptr<Tensor> L  = ovov[5];
    std::shared_ptr<Tensor> S  = ovov[7];

    Xi->dimensions()[0] = "naocc";
    Xi->dimensions()[1] = "ngrid";
    Xa->dimensions()[0] = "navir";
    Xa->dimensions()[1] = "ngrid";
    Z->dimensions()[0] = "ngrid";
    Z->dimensions()[1] = "ngrid";
    L->dimensions()[0] = "ngrid";
    L->dimensions()[1] = "naux";
    S->dimensions()[0] = "ngrid";
    S->dimensions()[1] = "ngrid";

    int ngrid = L->sizes()[0];
    thce_->new_dimension("ngrid", ngrid);
    int naux = L->sizes()[1];
    thce_->new_dimension("naux", naux);

    thce_->add_tensor("Xi",Xi);
    thce_->add_tensor("Xa",Xa);
    thce_->add_tensor("Ziaia",Z);
    thce_->add_tensor("Lia",L);
    thce_->add_tensor("Sia",S);
}
void RTHCEW::build_lsthc_act(std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> X)
{
    std::shared_ptr<LSTHCERI> lsthceri = LSTHCERI::build(basisset_,auxiliary,X,options_,reference_wavefunction_);
    lsthceri->add_eri_space("oooo","ACTIVE_OCC","ACTIVE_OCC","ACTIVE_OCC","ACTIVE_OCC");
    lsthceri->add_eri_space("ooov","ACTIVE_OCC","ACTIVE_OCC","ACTIVE_OCC","ACTIVE_VIR");
    lsthceri->add_eri_space("oovv","ACTIVE_OCC","ACTIVE_OCC","ACTIVE_VIR","ACTIVE_VIR");
    lsthceri->add_eri_space("ovov","ACTIVE_OCC","ACTIVE_VIR","ACTIVE_OCC","ACTIVE_VIR");
    lsthceri->add_eri_space("ovvv","ACTIVE_OCC","ACTIVE_VIR","ACTIVE_VIR","ACTIVE_VIR");
    lsthceri->add_eri_space("vvvv","ACTIVE_VIR","ACTIVE_VIR","ACTIVE_VIR","ACTIVE_VIR");
    lsthceri->compute();
    std::vector<std::shared_ptr<Tensor> > oooo = lsthceri->ints()["oooo"];
    std::vector<std::shared_ptr<Tensor> > ooov = lsthceri->ints()["ooov"];
    std::vector<std::shared_ptr<Tensor> > oovv = lsthceri->ints()["oovv"];
    std::vector<std::shared_ptr<Tensor> > ovov = lsthceri->ints()["ovov"];
    std::vector<std::shared_ptr<Tensor> > ovvv = lsthceri->ints()["ovvv"];
    std::vector<std::shared_ptr<Tensor> > vvvv = lsthceri->ints()["vvvv"];
    lsthceri.reset();

    std::shared_ptr<Tensor> Xi = ovov[0];
    std::shared_ptr<Tensor> Xa = ovov[1];

    std::shared_ptr<Tensor> Lii  = oooo[5];
    std::shared_ptr<Tensor> Sii  = oooo[7];
    std::shared_ptr<Tensor> Lia  = ovov[5];
    std::shared_ptr<Tensor> Sia  = ovov[7];
    std::shared_ptr<Tensor> Laa  = vvvv[5];
    std::shared_ptr<Tensor> Saa  = vvvv[7];

    std::shared_ptr<Tensor> Ziiii  = oooo[2];
    std::shared_ptr<Tensor> Ziiia  = ooov[2];
    std::shared_ptr<Tensor> Ziiaa  = oovv[2];
    std::shared_ptr<Tensor> Ziaia  = ovov[2];
    std::shared_ptr<Tensor> Ziaaa  = ovvv[2];
    std::shared_ptr<Tensor> Zaaaa  = vvvv[2];

    Xi->dimensions()[0] = "naocc";
    Xi->dimensions()[1] = "ngrid";
    Xa->dimensions()[0] = "navir";
    Xa->dimensions()[1] = "ngrid";
    Lii->dimensions()[0] = "ngrid";
    Lii->dimensions()[1] = "naux";
    Sii->dimensions()[0] = "ngrid";
    Sii->dimensions()[1] = "ngrid";
    Lia->dimensions()[0] = "ngrid";
    Lia->dimensions()[1] = "naux";
    Sia->dimensions()[0] = "ngrid";
    Sia->dimensions()[1] = "ngrid";
    Laa->dimensions()[0] = "ngrid";
    Laa->dimensions()[1] = "naux";
    Saa->dimensions()[0] = "ngrid";
    Saa->dimensions()[1] = "ngrid";
    Ziiii->dimensions()[0] = "ngrid";
    Ziiii->dimensions()[1] = "ngrid";
    Ziiia->dimensions()[0] = "ngrid";
    Ziiia->dimensions()[1] = "ngrid";
    Ziiaa->dimensions()[0] = "ngrid";
    Ziiaa->dimensions()[1] = "ngrid";
    Ziaia->dimensions()[0] = "ngrid";
    Ziaia->dimensions()[1] = "ngrid";
    Ziaaa->dimensions()[0] = "ngrid";
    Ziaaa->dimensions()[1] = "ngrid";
    Zaaaa->dimensions()[0] = "ngrid";
    Zaaaa->dimensions()[1] = "ngrid";

    int ngrid = Lia->sizes()[0];
    thce_->new_dimension("ngrid", ngrid);
    int naux = Lia->sizes()[1];
    thce_->new_dimension("naux", naux);

    thce_->add_tensor("Xi",Xi);
    thce_->add_tensor("Xa",Xa);
    thce_->add_tensor("Lia",Lia);
    thce_->add_tensor("Sia",Sia);
    thce_->add_tensor("Ziiii",Ziiii);
    thce_->add_tensor("Ziiia",Ziiia);
    thce_->add_tensor("Ziiaa",Ziiaa);
    thce_->add_tensor("Ziaia",Ziaia);
    thce_->add_tensor("Ziaaa",Ziaaa);
    thce_->add_tensor("Zaaaa",Zaaaa);
}
void RTHCEW::build_lsthc_pp(std::shared_ptr<BasisSet> auxiliary, std::shared_ptr<Matrix> X)
{
    std::shared_ptr<LSTHCERI> lsthceri = LSTHCERI::build(basisset_,auxiliary,X,options_,reference_wavefunction_);
    lsthceri->add_eri_space("pppp","ACTIVE_ALL","ACTIVE_ALL","ACTIVE_ALL","ACTIVE_ALL");
    lsthceri->compute();
    std::vector<std::shared_ptr<Tensor> > pppp = lsthceri->ints()["pppp"];
    lsthceri.reset();

    std::shared_ptr<Tensor> X2 = pppp[0];
    std::shared_ptr<Tensor> Z  = pppp[2];
    std::shared_ptr<Tensor> L  = pppp[5];
    std::shared_ptr<Tensor> S  = pppp[7];

    int naocc = thce_->dimensions()["naocc"];
    int navir = thce_->dimensions()["naocc"];
    int ngrid = L->sizes()[0];
    thce_->new_dimension("ngrid", ngrid);
    int naux = L->sizes()[1];
    thce_->new_dimension("naux", naux);

    std::shared_ptr<Tensor> Xi = CoreTensor::build("Xi","naocc",naocc,"ngrid",ngrid,X2->pointer() + 0L,false);
    std::shared_ptr<Tensor> Xa = CoreTensor::build("Xa","navir",navir,"ngrid",ngrid,X2->pointer() + naocc * (size_t) ngrid,false);

    Z->dimensions()[0] = "ngrid";
    Z->dimensions()[1] = "ngrid";
    L->dimensions()[0] = "ngrid";
    L->dimensions()[1] = "naux";
    S->dimensions()[0] = "ngrid";
    S->dimensions()[1] = "ngrid";

    thce_->add_tensor("Xi",Xi);
    thce_->add_tensor("Xa",Xa);
    thce_->add_tensor("Zpppp",Z);
    thce_->add_tensor("Lpp",L);
    thce_->add_tensor("Spp",S);
}
void RTHCEW::build_meth_ia(std::shared_ptr<Matrix> X)
{
    std::shared_ptr<LSTHCERI> lsthceri = LSTHCERI::build(basisset_,std::shared_ptr<BasisSet>(),X,options_,reference_wavefunction_);
    lsthceri->add_eri_space("ovov","ACTIVE_OCC","ACTIVE_VIR","ACTIVE_OCC","ACTIVE_VIR");
    lsthceri->compute_meth();
    std::vector<std::shared_ptr<Tensor> > ovov = lsthceri->meths()["ovov"];
    lsthceri.reset();

    std::shared_ptr<Tensor> Xi = ovov[0];
    std::shared_ptr<Tensor> Xa = ovov[1];
    std::shared_ptr<Tensor> S  = ovov[2];

    Xi->dimensions()[0] = "naocc";
    Xi->dimensions()[1] = "namp";
    Xa->dimensions()[0] = "navir";
    Xa->dimensions()[1] = "namp";
    S->dimensions()[0] = "namp";
    S->dimensions()[1] = "namp";

    int namp = S->sizes()[0];
    thce_->new_dimension("namp", namp);

    thce_->add_tensor("Ti",Xi);
    thce_->add_tensor("Ta",Xa);
    thce_->add_tensor("STia",S);
}

} // Namespace
