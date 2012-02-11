#include "vector3.h"
#include "integral.h"
#include "gshell.h"
#include "basisset.h"
#include "sobasis.h"
#include "sointegral.h"
#include <boost/shared_ptr.hpp>
#include <algorithm>

using namespace boost;
using namespace psi;

AOIntegralsIterator::AOIntegralsIterator(const GaussianShell& s1, const GaussianShell& s2,
                                     const GaussianShell& s3, const GaussianShell& s4)
    : usi(s1), usj(s2), usk(s3), usl(s4)
{
    done = false;
    ni = usi.nfunction();
    nj = usj.nfunction();
    nk = usk.nfunction();
    nl = usl.nfunction();

    fii = usi.function_index();
    fij = usj.function_index();
    fik = usk.function_index();
    fil = usl.function_index();

    iimax = ni - 1;
    if (&usi == &usj && &usk == &usl && &usi == &usk) {
        kkmax = 0;
        llmax = 0;
        jjmax = 0;
    }
    else if(&usi == &usk && &usj == &usl){
        kkmax = 0;
        llmax = 0;
        jjmax = nj - 1;
    }
    else{
        kkmax = nk - 1;
        jjmax = (&usi == &usj) ? 0 : nj - 1;
        llmax = (&usk == &usl) ? 0 : nl - 1;
    }

    ii = 0;
    jj = 0;
    kk = 0;
    ll = 0;
}

void AOIntegralsIterator::first()
{
    current.i = 0 + fii;
    current.j = 0 + fij;
    current.k = 0 + fik;
    current.l = 0 + fil;
    current.index = 0;
    if (&usi == &usj && &usk == &usl && &usi == &usk) {     // (aa|aa) case
    }
    else if(&usi== &usk && &usj == &usl){
        if (current.i < current.j) {
            std::swap(current.i, current.j);
            std::swap(current.k, current.l);
        }
        if (current.i < current.k) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
        }
    }
    else{
        if (current.i < current.j) {
            std::swap(current.i, current.j);
        }
        if (current.k < current.l) {
            std::swap(current.k, current.l);
        }
        if ((current.i < current.k) || (current.i == current.k && current.j < current.l)) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
        }
    }
}

void AOIntegralsIterator::next()
{
    if (&usi == &usj && &usk == &usl && &usi == &usk) {
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
    else if(&usi == &usk && &usj == &usl){ //(ab|ab)
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
            std::swap(current.i, current.j);
            std::swap(current.k, current.l);
        }
        if (current.i < current.k) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
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
                    jjmax = (&usi == &usj) ? ii : nj - 1;
                }
            }
            llmax = (&usk == &usl) ? kk : nl - 1;
        }
        current.i = ii + fii;
        current.j = jj + fij;
        current.k = kk + fik;
        current.l = ll + fil;
        current.index = ll+nl*(kk+nk*(jj+nj*ii));
        if (current.i < current.j) {
            std::swap(current.i, current.j);
        }
        if (current.k < current.l) {
            std::swap(current.k, current.l);
        }
        if ((current.i < current.k) || (current.i == current.k && current.j < current.l)) {
            std::swap(current.i, current.k);
            std::swap(current.j, current.l);
        }
    }

}

// ===========================================================================
//  AOShellCombinationsIterator
// ===========================================================================
AOShellCombinationsIterator::AOShellCombinationsIterator(boost::shared_ptr<BasisSet>bs1, boost::shared_ptr<BasisSet>bs2,
                                                         boost::shared_ptr<BasisSet>bs3, boost::shared_ptr<BasisSet>bs4) :
    bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4)
{

}

AOShellCombinationsIterator::AOShellCombinationsIterator()
{

}

AOIntegralsIterator AOShellCombinationsIterator::integrals_iterator()
{
    return AOIntegralsIterator(bs1_->shell(p()), bs2_->shell(q()), bs3_->shell(r()), bs4_->shell(s()));
}

void AOShellCombinationsIterator::init(boost::shared_ptr<BasisSet>bs1, boost::shared_ptr<BasisSet>bs2,
                                     boost::shared_ptr<BasisSet>bs3, boost::shared_ptr<BasisSet>bs4)
{
    bs1_=bs1;
    bs2_=bs2;
    bs3_=bs3;
    bs4_=bs4;
}

void AOShellCombinationsIterator::first()
{
    usii = usjj = uskk = usll = upk = 0;
    done = false;

    num_unique_pk = 1;
    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;

    int usi, usj, usk, usl;
    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];

    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->shell(usi).am() < bs2_->shell(usj).am()) {
        std::swap(usi, usj);
    }
    if (bs3_->shell(usk).am() < bs4_->shell(usl).am()) {
        std::swap(usk, usl);
    }
    if (bs1_->shell(usi).am() + bs2_->shell(usj).am() >
            bs3_->shell(usk).am() + bs4_->shell(usl).am()) {
        std::swap(usi, usk);
        std::swap(usj, usl);
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

void AOShellCombinationsIterator::next()
{
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
        if ((usii == usjj && usii == uskk) || (usjj == uskk && usjj == usll))
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
    if (bs1_->shell(usi).am() < bs2_->shell(usj).am()) {
        std::swap(usi, usj);
    }
    if (bs3_->shell(usk).am() < bs4_->shell(usl).am()) {
        std::swap(usk, usl);
    }
    if (bs1_->shell(usi).am() + bs2_->shell(usj).am() >
            bs3_->shell(usk).am() + bs4_->shell(usl).am()) {
        std::swap(usi, usk);
        std::swap(usj, usl);
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

// ===========================================================================
//  SOShellCombinationsIterator
// ===========================================================================
SOShellCombinationsIterator::SOShellCombinationsIterator(boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
                                                         boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4) :
    bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4)
{

}

SOShellCombinationsIterator::SOShellCombinationsIterator()
{

}

void SOShellCombinationsIterator::init(boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
                                       boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4)
{
    bs1_=bs1;
    bs2_=bs2;
    bs3_=bs3;
    bs4_=bs4;
}

void SOShellCombinationsIterator::first()
{
    usii = usjj = uskk = usll = upk = 0;
    done = false;

    num_unique_pk = 1;
    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;

    int usi, usj, usk, usl;
    usi = usi_arr[upk]; usj = usj_arr[upk]; usk = usk_arr[upk]; usl = usl_arr[upk];

    current.P = usi; current.Q = usj; current.R = usk; current.S = usl; current.end_of_PK = false;

    if (upk == num_unique_pk - 1) {
        // If this is the last unique shell flag it as end of a pk block.
        current.end_of_PK = true;
    }
    else{
        current.end_of_PK = false;
    }
}

void SOShellCombinationsIterator::next()
{
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
//        fprintf(outfile, ">usii %d usjj %d uskk %d usll %d\n", usii, usjj, uskk, usll);

        usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
        if ((usii == usjj && usii == uskk) || (usjj == uskk && usjj == usll))
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

//    fprintf(outfile, ">si %d usj %d usk %d usl %d\n", usi, usj, usk, usl);

    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->am(usi) < bs2_->am(usj)) {
        std::swap(usi, usj);
    }
    if (bs3_->am(usk) < bs4_->am(usl)) {
        std::swap(usk, usl);
    }
    if (bs1_->am(usi) + bs2_->am(usj) >
            bs3_->am(usk) + bs4_->am(usl)) {
        std::swap(usi, usk);
        std::swap(usj, usl);
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

// ===========================================================================
//  SO_PQ_Iterator
// ===========================================================================
SO_PQ_Iterator::SO_PQ_Iterator(boost::shared_ptr<SOBasisSet>bs1) :
    bs1_(bs1)
{

}

SO_PQ_Iterator::SO_PQ_Iterator()
{

}

void SO_PQ_Iterator::first()
{
    ii = jj = bs1_->nshell()-1;
    done = false;

    current.P = ii;
    current.Q = jj;
}

void SO_PQ_Iterator::next()
{
    if (jj > 0) {
        jj--;
    }
    else {
        ii--;
        jj = ii;
        if (ii < 0) {
            done = true;
            return;
        }
    }

    current.P = ii;
    current.Q = jj;
}

// ===========================================================================
//  SO_RS_Iterator
// ===========================================================================
SO_RS_Iterator::SO_RS_Iterator(const int &P, const int &Q,
                               boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
                               boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4) :
    usii(P), usjj(Q), bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4)
{
}

SO_RS_Iterator::SO_RS_Iterator(boost::shared_ptr<SOBasisSet>bs1, boost::shared_ptr<SOBasisSet>bs2,
                               boost::shared_ptr<SOBasisSet>bs3, boost::shared_ptr<SOBasisSet>bs4) :
     bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4)
{
}

SO_RS_Iterator::SO_RS_Iterator() : usii(0), usjj(0)
{

}

void SO_RS_Iterator::first()
{
    uskk = usll = upk = 0;
    done = false;

    int usi, usj, usk, usl;

//    fprintf(outfile, ">usii %d usjj %d uskk %d usll %d\n", usii, usjj, uskk, usll);

    usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
    if ((usii == usjj && usii == uskk) || (usjj == uskk && usjj == usll))
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

    usi = usii; usj = usjj; usk = uskk; usl = usll;

//    fprintf(outfile, ">si %d usj %d usk %d usl %d\n", usi, usj, usk, usl);

    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->am(usi) < bs2_->am(usj)) {
        std::swap(usi, usj);
    }
    if (bs3_->am(usk) < bs4_->am(usl)) {
        std::swap(usk, usl);
    }
    if (bs1_->am(usi) + bs2_->am(usj) >
            bs3_->am(usk) + bs4_->am(usl)) {
        std::swap(usi, usk);
        std::swap(usj, usl);
    }

    current.P = usi; current.Q = usj; current.R = usk; current.S = usl; //current.end_of_PK = false;
 }

void SO_RS_Iterator::next()
{
    ++upk;
    if(upk >= num_unique_pk){
        upk = 0;
//        if (usii == 0 && usjj == 0) {
//            done = true;
//            return;
//        }
        ++usll;
        if (usll > uskk){
            ++uskk;
            if ((usll-1) == usjj && (uskk-1) == usjj) {
                done = true;
                return;
            }
            usll = 0;
        }

//        fprintf(outfile, ">usii %d usjj %d uskk %d usll %d\n", usii, usjj, uskk, usll);

        usi_arr[0] = usii; usj_arr[0] = usjj; usk_arr[0] = uskk; usl_arr[0] = usll;
        if ((usii == usjj && usii == uskk) || (usjj == uskk && usjj == usll))
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

//    fprintf(outfile, ">si %d usj %d usk %d usl %d\n", usi, usj, usk, usl);

    // Sort shells based on AM, save ERI some work doing permutation resorting.
    if (bs1_->am(usi) < bs2_->am(usj)) {
        std::swap(usi, usj);
    }
    if (bs3_->am(usk) < bs4_->am(usl)) {
        std::swap(usk, usl);
    }
    if (bs1_->am(usi) + bs2_->am(usj) >
            bs3_->am(usk) + bs4_->am(usl)) {
        std::swap(usi, usk);
        std::swap(usj, usl);
    }

    current.P = usi; current.Q = usj; current.R = usk; current.S = usl; //current.end_of_PK = false;

//    if (upk == num_unique_pk - 1) {
//        // If this is the last unique shell flag it as end of a pk block.
//        current.end_of_PK = true;
//    }
//    else{
//        current.end_of_PK = false;
//    }
}
