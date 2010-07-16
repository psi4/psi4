#include <libmints/basisset.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/overlap.h>
#include <libmints/kinetic.h>
#include <libmints/potential.h>
#include <libmints/electrostatic.h>
#include <libmints/integral.h>
#include <libmints/dipole.h>
#include <libmints/quadrupole.h>
#include <libmints/eri.h>
#include <libmints/electricfield.h>

using namespace psi;

namespace psi {
template <class T>
static void swap(T& x, T& y) {
    T tmp=x; x = y; y = tmp;
}
}

/** Initialize IntegralFactory object given a GaussianBasisSet for each center. */
IntegralFactory::IntegralFactory(shared_ptr<BasisSet> bs1, shared_ptr<BasisSet> bs2,
                shared_ptr<BasisSet> bs3, shared_ptr<BasisSet> bs4)
{
    set_basis(bs1, bs2, bs3, bs4);
}

IntegralFactory::~IntegralFactory()
{

}

void IntegralFactory::set_basis(shared_ptr<BasisSet> bs1, shared_ptr<BasisSet> bs2,
                shared_ptr<BasisSet> bs3, shared_ptr<BasisSet> bs4)
{
    bs1_ = bs1;
    bs2_ = bs2;
    bs3_ = bs3;
    bs4_ = bs4;

    // Find the max am
    shared_ptr<BasisSet> max12(bs1_->max_am() > bs2_->max_am() ? bs1_ : bs2_);
    shared_ptr<BasisSet> max34(bs3_->max_am() > bs4_->max_am() ? bs3_ : bs4_);
    shared_ptr<BasisSet> max1234(max12->max_am() > max34->max_am() ? max12 : max34);

    init_spherical_harmonics(max1234->max_am());
}

OneBodyInt* IntegralFactory::overlap(int deriv)
{
    return new OverlapInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodyInt* IntegralFactory::kinetic(int deriv)
{
    return new KineticInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodyInt* IntegralFactory::potential(int deriv)
{
    return new PotentialInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodyInt* IntegralFactory::electrostatic()
{
    return new ElectrostaticInt(spherical_transforms_, bs1_, bs2_, 0);
}

OneBodyInt* IntegralFactory::dipole(int deriv)
{
    return new DipoleInt(spherical_transforms_, bs1_, bs2_, deriv);
}

OneBodyInt* IntegralFactory::quadrupole()
{
    return new QuadrupoleInt(spherical_transforms_, bs1_, bs2_);
}

OneBodyInt* IntegralFactory::electric_field()
{
    return new ElectricFieldInt(spherical_transforms_, bs1_, bs2_);
}

TwoBodyInt* IntegralFactory::eri(int deriv, double schwarz)
{
    return new ERI(bs1_, bs2_, bs3_, bs4_, deriv, schwarz);
}

void IntegralFactory::init_spherical_harmonics(int max_am)
{
    for (int i=0; i<=max_am; ++i)
        spherical_transforms_.push_back(SphericalTransform(i));
}

ShellCombinationsIterator IntegralFactory::shells_iterator()
{
    return ShellCombinationsIterator(bs1_, bs2_, bs3_, bs4_);
}

IntegralsIterator ShellCombinationsIterator::integrals_iterator()
{
    return IntegralsIterator(bs1_->shell(p()), bs2_->shell(q()), bs3_->shell(r()), bs4_->shell(s()));
}
IntegralsIterator IntegralFactory::integrals_iterator(int p, int q, int r, int s)
{
    return IntegralsIterator(bs1_->shell(p), bs2_->shell(q), bs3_->shell(r), bs4_->shell(s));
}

/*
void ShellCombinationsIterator::generate_combinations(BasisSet*bs1, BasisSet*bs2, BasisSet*bs3, BasisSet*bs4)
{

    for (usii=0; usii<bs1->nshell(); usii++) {
        for (usjj=0; usjj<=usii; usjj++) {
            for (uskk=0; uskk<=usjj; uskk++) {
                for (usll=0; usll<=uskk; usll++) {
                    // Decide what shell quartets out of (ij|kl), (ik|jl), and (il|jk) are unique
                    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
                    if (usii == usjj && usii == uskk || usjj == uskk && usjj == usll)
                        num_unique_pk = 1;
                    else if (usii == uskk || usjj == usll) {
                        num_unique_pk = 2;
                        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
                    }
                    else if (usjj == uskk) {
                        num_unique_pk = 2;
                        usi_arr[1] = usii; usj_arr[1] = usll; usk_arr[1] = usjj; usl_arr[1] = uskk;
                    }
                    else if (usii == usjj || uskk == usll) {
                        num_unique_pk = 2;
                        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
                    }
                    else {
                        num_unique_pk = 3;
                        usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
                        usi_arr[2] = usii; usj_arr[2] = usll; usk_arr[2] = usjj; usl_arr[2] = uskk;
                    }

                    // For each num_unique_pk
                    for (int upk=0; upk < num_unique_pk; ++upk) {
                        usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];

                        // Sort shells based on AM, save ERI some work doing permutation resorting.
                        if (bs1->shell(usi)->am() < bs2->shell(usj)->am()) {
                            swap(usi, usj);
                        }
                        if (bs3->shell(usk)->am() < bs4->shell(usl)->am()) {
                            swap(usk, usl);
                        }
                        if (bs1->shell(usi)->am() + bs2->shell(usj)->am() >
                            bs3->shell(usk)->am() + bs4->shell(usl)->am()) {
                            swap(usi, usk);
                            swap(usj, usl);
                        }

                        ShellQuartet q;
                        q.P = usi; q.Q = usj; q.R = usk; q.S = usl; q.end_of_PK = false;

                        if (upk == num_unique_pk - 1) {
                            // If this is the last unique shell flag it as end of a pk block.
                            q.end_of_PK = true;
                        }
                        unique_quartets_.push_back(q);
                    }
                }
            }
        }
    }
}*/

void ShellCombinationsIterator::first(){
    usii = usjj = uskk = usll = upk = 0;
    done = false;

    num_unique_pk = 1;
    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;

    int usi, usj, usk, usl;
    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];

    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->shell(usi)->am() < bs2_->shell(usj)->am()) {
        swap(usi, usj);
    }
    if (bs3_->shell(usk)->am() < bs4_->shell(usl)->am()) {
        swap(usk, usl);
    }
    if (bs1_->shell(usi)->am() + bs2_->shell(usj)->am() >
        bs3_->shell(usk)->am() + bs4_->shell(usl)->am()) {
        swap(usi, usk);
        swap(usj, usl);
    }

    current.P = usi; current.Q = usj; current.R = usk; current.S = usl; current.end_of_PK = false;

    if (upk == num_unique_pk - 1) {
        // If this is the last unique shell flag it as end of a pk block.
        current.end_of_PK = true;
    }
    else{
        current.end_of_PK = false;
    }

}


void ShellCombinationsIterator::next(){
    ++upk;
    if(upk >= num_unique_pk){
        upk = 0;
        ++usll;
        if (usll > uskk){
            ++uskk;
            usll = 0;
            if(uskk > usjj){
                ++usjj;
                uskk = 0;
                if(usjj > usii){
                    ++usii;
                    usjj = 0;
                    if(usii >= bs1_->nshell()){
                        done = true;
                        return;
                    }
                }
            }
        }
        usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
        if (usii == usjj && usii == uskk || usjj == uskk && usjj == usll)
            num_unique_pk = 1;
        else if (usii == uskk || usjj == usll) {
            num_unique_pk = 2;
            usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
        }
        else if (usjj == uskk) {
            num_unique_pk = 2;
            usi_arr[1] = usii; usj_arr[1] = usll; usk_arr[1] = usjj; usl_arr[1] = uskk;
        }
        else if (usii == usjj || uskk == usll) {
            num_unique_pk = 2;
            usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
        }
        else {
            num_unique_pk = 3;
            usi_arr[1] = usii; usj_arr[1] = uskk; usk_arr[1] = usjj; usl_arr[1] = usll;
            usi_arr[2] = usii; usj_arr[2] = usll; usk_arr[2] = usjj; usl_arr[2] = uskk;
        }
    }



    int usi, usj, usk, usl;
    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];


    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->shell(usi)->am() < bs2_->shell(usj)->am()) {
        swap(usi, usj);
    }
    if (bs3_->shell(usk)->am() < bs4_->shell(usl)->am()) {
        swap(usk, usl);
    }
    if (bs1_->shell(usi)->am() + bs2_->shell(usj)->am() >
        bs3_->shell(usk)->am() + bs4_->shell(usl)->am()) {
        swap(usi, usk);
        swap(usj, usl);
    }

    current.P = usi; current.Q = usj; current.R = usk; current.S = usl; current.end_of_PK = false;

    if (upk == num_unique_pk - 1) {
        // If this is the last unique shell flag it as end of a pk block.
        current.end_of_PK = true;
    }
    else{
        current.end_of_PK = false;
    }

}




void IntegralsIterator::first(){
    current.i = 0 + fii;
    current.j = 0 + fij;
    current.k = 0 + fik;
    current.l = 0 + fil;
    current.index = 0;
    if (usi == usj && usk == usl && usi == usk) {     // (aa|aa) case
    }
    else if(usi== usk && usj == usl){
        if (current.i < current.j) {
            swap(current.i, current.j);
            swap(current.k, current.l);
        }
        if (current.i < current.k) {
            swap(current.i, current.k);
            swap(current.j, current.l);
        }
    }
    else{
        if (current.i < current.j) {
            swap(current.i, current.j);
        }
        if (current.k < current.l) {
            swap(current.k, current.l);
        }
        if ((current.i < current.k) || (current.i == current.k && current.j < current.l)) {
            swap(current.i, current.k);
            swap(current.j, current.l);
        }
    }
}


void IntegralsIterator::next(){
    if (usi == usj && usk == usl && usi == usk) {
        ++ll;
        if(ll > llmax){
            ++kk;
            ll = 0;
            if(kk > kkmax){
                kk = 0;
                ++jj;
                if(jj > jjmax){
                    jj = 0;
                    ++ii;
                    if(ii > iimax){
                        done = true;
                    }
                    jjmax = ii;
                }
                kkmax = ii;

            }
            llmax = (kk==ii) ? jj : kk;
        }
        current.i = ii + fii;
        current.j = jj + fij;
        current.k = kk + fik;
        current.l = ll + fil;
        current.index = ll+nl*(kk+nk*(jj+nj*ii));

    }
    else if(usi == usk && usj == usl){ //(ab|ab)
        ++ll;
        if(ll > llmax){
            ++kk;
            ll = 0;
            if(kk > kkmax){
                kk = 0;
                ++jj;
                if(jj > jjmax){
                    jj = 0;
                    ++ii;
                    if(ii > iimax){
                        done = true;
                    }
                }
                kkmax = ii;
            }
            llmax = (kk == ii) ? jj : nl - 1;
        }
        current.i = ii + fii;
        current.j = jj + fij;
        current.k = kk + fik;
        current.l = ll + fil;
        current.index = ll+nl*(kk+nk*(jj+nj*ii));
        if (current.i < current.j) {
            swap(current.i, current.j);
            swap(current.k, current.l);
        }
        if (current.i < current.k) {
            swap(current.i, current.k);
            swap(current.j, current.l);
        }
    }
    else{
        ++ll;
        if(ll > llmax){
            ++kk;
            ll = 0;
            if(kk > kkmax){
                kk = 0;
                ++jj;
                if(jj > jjmax){
                    jj = 0;
                    ++ii;
                    if(ii > iimax){
                        done = true;
                    }
                    jjmax = (usi == usj) ? ii : nj - 1;
                }
            }
            llmax = (usk==usl) ? kk : nl - 1;
        }
        current.i = ii + fii;
        current.j = jj + fij;
        current.k = kk + fik;
        current.l = ll + fil;
        current.index = ll+nl*(kk+nk*(jj+nj*ii));
        if (current.i < current.j) {
            swap(current.i, current.j);
        }
        if (current.k < current.l) {
            swap(current.k, current.l);
        }
        if ((current.i < current.k) || (current.i == current.k && current.j < current.l)) {
            swap(current.i, current.k);
            swap(current.j, current.l);
        }
    }

}

