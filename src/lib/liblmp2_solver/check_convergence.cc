/*! \defgroup LMP2 lmp2: LMP2 Evaluation of Energy */

/*!
 ** \file
 ** \ingroup LMP2
 ** \LMP2 evaluation of energy
 */

#include <liblmp2_solver/lmp2.h>
#include <libqt/qt.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{

class BasisSet;
class Options;
class PSIO;
class Chkpt;

namespace lmp2 {

#ifdef HAVE_MADNESS

madness::Future<double> LMP2::T2_rms(const SharedMatrix T2,
                                     const SharedMatrix T2_old,
                                     const int &ij)
{

    int i = ij_i_map_[ij];
    int j = ij_j_map_[ij];

    double drms = 0.0;

    SharedMatrix temp = SharedMatrix(T2->clone());
    temp->copy(T2);
    temp->subtract(T2_old);
    if (i != j)
        drms += 2 * temp->vector_dot(temp);
    else
        drms += temp->vector_dot(temp);

//    for(int a=0; a < pair_domain_len_[ij]; a++) {
//        for(int b=0; b < pair_domain_len_[ij]; b++) {
//            double Tab = T2->get(0,a,b);
//            double Tab_old = T2_old->get(0,a,b);

//            if(fabs(Tab) > 1e-14 && fabs(Tab_old) > 1e-14) {
//                if (i != j)
//                    drms += 2 * (Tab - Tab_old) *
//                            (Tab - Tab_old);
//                else
//                    drms += (Tab - Tab_old) *
//                            (Tab - Tab_old);
//            }
//        }
//    }

    return madness::Future<double>(drms);
}

int LMP2::store_T2(const SharedMatrix T2, const int &ij) {
    T2_amp_[div_][ij]->copy(T2);
    return 0;
}

#endif // Have_madness
}} // End of psi and lmp2 namespaces
