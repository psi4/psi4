#include "vector3.h"
#include "integral.h"
#include "gshell.h"

#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace psi;

IntegralsIterator::IntegralsIterator(shared_ptr<GaussianShell> s1, shared_ptr<GaussianShell> s2,
                      shared_ptr<GaussianShell> s3, shared_ptr<GaussianShell> s4) 
{
    done = false;
    usi = s1;
    usj = s2;
    usk = s3;
    usl = s4;
    ni =usi->nfunction();
    nj =usj->nfunction();
    nk =usk->nfunction();
    nl =usl->nfunction();
    
    fii = usi->function_index();
    fij = usj->function_index();
    fik = usk->function_index();
    fil = usl->function_index();
    
    iimax = ni - 1;
    if (usi == usj && usk == usl && usi == usk) {
        kkmax = 0;
        llmax = 0;
        jjmax = 0;
    }
    else if(usi == usk && usj == usl){
        kkmax = 0;
        llmax = 0;
        jjmax = nj - 1;
    }
    else{
        kkmax = nk - 1;
        jjmax = (usi == usj) ? 0 : nj - 1;
        llmax = (usk == usl) ? 0 : nl - 1;
    }
    
    ii = 0;
    jj = 0;
    kk = 0;
    ll = 0;
}

ShellCombinationsIterator::ShellCombinationsIterator(shared_ptr<BasisSet>bs1, shared_ptr<BasisSet>bs2,
                              shared_ptr<BasisSet>bs3, shared_ptr<BasisSet>bs4) :
                       bs1_(bs1), bs2_(bs2), bs3_(bs3), bs4_(bs4)
{

}

ShellCombinationsIterator::ShellCombinationsIterator()
{

}

void ShellCombinationsIterator::init(shared_ptr<BasisSet>bs1, shared_ptr<BasisSet>bs2,
            shared_ptr<BasisSet>bs3, shared_ptr<BasisSet>bs4)

{
    bs1_=bs1; 
    bs2_=bs2;
    bs3_=bs3; 
    bs4_=bs4;
}

