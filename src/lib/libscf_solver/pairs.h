/*! \file
    \ingroup LMP2
    \brief Enter brief description of file here
*/
#ifndef _psi_src_lib_libscf_solver_pairs_
#define _psi_src_lib_libscf_solver_pairs_

#include <cstdio>
#include <libciomr/libciomr.h>
#include <libparallel/parallel.h>
#include <libqt/qt.h>
#include "rhf.h"
#include <psi4-dec.h>

#ifdef HAVE_MADNESS

using namespace madness;

namespace psi{ namespace scf{

class PARALLEL_G_BUILD_INFO {
    private:
        boost::shared_ptr<madness::Spinlock> eri_lock;
        boost::shared_ptr<BasisSet> basis_info;
        boost::shared_ptr<IntegralFactory> integral;
        boost::shared_ptr<TwoBodyAOInt> eri_info;
        boost::shared_ptr<SimpleMatrix> pG;
        boost::shared_ptr<SimpleMatrix> pD;
        int nso;

    public:

        PARALLEL_G_BUILD_INFO() {};
        ~PARALLEL_G_BUILD_INFO() {};

        void initialize(boost::shared_ptr<BasisSet> _basis, SharedMatrix sh_D, const int & _nso) {
            if(basis_info.get() == NULL)
                basis_info = _basis;

            nso = _nso;

            if(integral.get() == NULL)
                integral = boost::shared_ptr<IntegralFactory> (new IntegralFactory(basis_info, basis_info, basis_info, basis_info));

            if (eri_info.get() == NULL)
                eri_info = boost::shared_ptr<TwoBodyAOInt>(integral->eri());

            if(pD.get() == NULL)
                pD = boost::shared_ptr<SimpleMatrix>(sh_D->to_simple_matrix());
            else {
                pD->zero();
                pD = boost::shared_ptr<SimpleMatrix>(sh_D->to_simple_matrix());
            }

            if(pG.get() == NULL) {
                pG = boost::shared_ptr<SimpleMatrix>(new SimpleMatrix(nso, nso));
                pG->zero();
            }
            else {
                pG = boost::shared_ptr<SimpleMatrix>(new SimpleMatrix(nso, nso));
                pG->zero();
            }

            if(eri_lock.get() == NULL)
                eri_lock = boost::shared_ptr<madness::Spinlock> (Communicator::world->mad_mutex());
        }

        boost::shared_ptr<IntegralsIterator> create_int_iter(const int & P, const int & Q,
          const int & R, const int & S) {

            return boost::shared_ptr<IntegralsIterator>(new IntegralsIterator(basis_info->shell(P), basis_info->shell(Q),
               basis_info->shell(R), basis_info->shell(S)));

        }

        boost::shared_ptr<ShellCombinationsIterator> create_shell_iter() {
            return boost::shared_ptr<ShellCombinationsIterator>(new ShellCombinationsIterator(basis_info, basis_info, basis_info, basis_info));
        }


        void compute_integrals(const int &P, const int &Q,
          const int &R, const int &S) {

            eri_info->compute_shell(P, Q, R, S);
        }

        const double* get_buffer() {
            return eri_info->buffer();
        }

        double get_pD_val(const int & _i, const int & _j) {
            return pD->get(_i, _j);
        }

       void add_pG(const int & _i, const int & _j, const double & _temp) {
            pG->add(_i, _j, _temp);
        }

       void sum_G() {
           Communicator::world->sum(pG->ptr(), nso*nso, pG->ptr(), 0);
       }

       double* get_ptr_pG() {
           return pG->ptr();
       }

       void set_G_(SharedMatrix _G) {
           _G->set(pG);
       }

        void lock_mutex() {
            eri_lock->lock();
        }

        void unlock_mutex() {
            eri_lock->unlock();
        }

};

boost::shared_ptr<PARALLEL_G_BUILD_INFO> g_info(new PARALLEL_G_BUILD_INFO());
//boost::shared_ptr<PARALLEL_G_BUILD_INFO> g_info();

inline int integral_type(int i, int j, int k, int l)
{
    int type;

    if (i == j && i == k && i == l)     // (ij|kl)  (11|11)
        type = 1;
    else if (i == j && k == l && i > k) // (ij|kl)  (22|11)
        type = 2;
    else if (i == j && i == k && i > l) // (ij|kl)  (22|21)
        type = 3;
    else if (j == k && j == l && i > j) // (ij|kl)  (21|11)
        type = 4;
    else if (i == k && j == l)          // (ij|kl)  (21|21)
        type = 5;
    else if (i == j)                    // (ij|kl)  (33|21)
        type = 6;
    else if (j >  k && k == l)          // (ij|kl)  (32|11)
        type = 7;
    else if (k == l)                    // (ij|kl)  (31|22)
        type = 8;
    else if (i == k)                    // (ij|kl)  (32|31)
        type = 9;
    else if (j == k)                    // (ij|kl)  (32|21)
        type = 10;
    else if (j == l)                    // (ij|kl)  (31|21)
        type = 11;
    else if (j >  k)                    // (ij|kl)  (43|21)
        type = 12;
    else if (j >  l)                    // (ij|kl)  (42|31)
        type = 13;
    else                                // (ij|kl)  (41|32)
        type = 14;

    return type;
}

class Pair {
    private:
        int ij;

    public:
        Pair() {};

        Pair(int ij)
        : ij(ij){}

        ~Pair() {}

        bool operator==(const Pair& b) const {
            return ij==b.ij;
        }

        hashT hash() const {
            return ij;
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & ij;
        }
};


class G_MAT {
    private:

    public:
        G_MAT() {}

        ~G_MAT() {}

        Void build_G(const int & _P, const int _Q, const int & _R, const int & _S) {

            g_info->lock_mutex();

            double temp1, temp2, temp3, temp4, temp5, temp6, value;
            int itype;

            // Get the storage buffer from the eri object

            const double *buffer = g_info->get_buffer();

            g_info->compute_integrals(_P, _Q, _R, _S);

            boost::shared_ptr<IntegralsIterator> integral_iter = g_info->create_int_iter(_P, _Q, _R, _S);

            for(integral_iter->first(); !integral_iter->is_done(); integral_iter->next()) {
                int i = integral_iter->i();
                int j = integral_iter->j();
                int k = integral_iter->k();
                int l = integral_iter->l();
                int index = integral_iter->index();

                value = buffer[index];

                if (fabs(value) > 1.0e-14) {

                    itype = integral_type(i, j, k, l);
                    switch(itype) {
                        case 1:
                        temp1 = g_info->get_pD_val(i, i) * value;

                        g_info->add_pG(i, i, temp1);
                        break;

                        case 2:
                        temp1 = g_info->get_pD_val(k, k) * 2.0 * value;
                        temp2 = g_info->get_pD_val(i, k) * value;
                        temp3 = g_info->get_pD_val(i, i) * 2.0 * value;

                        g_info->add_pG(i, i, temp1);
                        g_info->add_pG(k, k, temp3);
                        g_info->add_pG(i, k, -temp2);
                        g_info->add_pG(k, i, -temp2);

                        break;

                        case 3:
                        temp1 = g_info->get_pD_val(i, i) * value;
                        temp2 = g_info->get_pD_val(i, l) * value * 2.0;

                        g_info->add_pG(i, l, temp1);
                        g_info->add_pG(l, i, temp1);
                        g_info->add_pG(i, i, temp2);

                        break;

                        case 4:
                        temp1 = g_info->get_pD_val(j, j) * value;
                        temp2 = g_info->get_pD_val(i, j) * value * 2.0;

                        g_info->add_pG(i, j, temp1);
                        g_info->add_pG(j, i, temp1);
                        g_info->add_pG(j, j, temp2);

                        break;

                        case 5:
                        temp1 = g_info->get_pD_val(i, j) * value * 3.0;
                        temp2 = g_info->get_pD_val(i, i) * value;
                        temp3 = g_info->get_pD_val(j, j) * value;

                        g_info->add_pG(i, j, temp1);
                        g_info->add_pG(j, i, temp1);
                        g_info->add_pG(j, j, -temp2);
                        g_info->add_pG(i, i, -temp3);

                        break;

                        case 6:
                        temp1 = g_info->get_pD_val(k, l) * value * 4.0;
                        temp2 = g_info->get_pD_val(i, l) * value;
                        temp3 = g_info->get_pD_val(i, i) * value * 2.0;
                        temp4 = g_info->get_pD_val(i, k) * value;

                        g_info->add_pG(i, i, temp1);
                        g_info->add_pG(i, k, -temp2);
                        g_info->add_pG(k, i, -temp2);
                        g_info->add_pG(k, l, temp3);
                        g_info->add_pG(l, k, temp3);
                        g_info->add_pG(i, l, -temp4);
                        g_info->add_pG(l, i, -temp4);

                        break;

                        case 7:
                        temp1 = g_info->get_pD_val(i, j) * value * 4.0;
                        temp2 = g_info->get_pD_val(j, k) * value;
                        temp3 = g_info->get_pD_val(i, k) * value;
                        temp4 = g_info->get_pD_val(k, k) * value * 2.0;

                        g_info->add_pG(k, k,  temp1);
                        g_info->add_pG(i, k, -temp2);
                        g_info->add_pG(k, i, -temp2);
                        g_info->add_pG(j, k, -temp3);
                        g_info->add_pG(k, j, -temp3);
                        g_info->add_pG(i, j,  temp4);
                        g_info->add_pG(j, i,  temp4);

                        break;

                        case 8:
                        temp1 = g_info->get_pD_val(k, k) * value * 2.0;
                        temp2 = g_info->get_pD_val(i, j) * value * 4.0;
                        temp3 = g_info->get_pD_val(j, k) * value;
                        temp4 = g_info->get_pD_val(i, k) * value;

                        g_info->add_pG(i, j, temp1);
                        g_info->add_pG(j, i, temp1);
                        g_info->add_pG(k, k, temp2);
                        g_info->add_pG(i, k, -temp3);
                        g_info->add_pG(k, i, -temp3);
                        g_info->add_pG(j, k, -temp4);
                        g_info->add_pG(k, j, -temp4);

                        break;

                        case 9:
                        temp1 = g_info->get_pD_val(i, l) * value * 3.0;
                        temp2 = g_info->get_pD_val(i, j) * value * 3.0;
                        temp3 = g_info->get_pD_val(j, l) * value * 2.0;
                        temp4 = g_info->get_pD_val(i, i) * value;

                        g_info->add_pG(i, j, temp1);
                        g_info->add_pG(j, i, temp1);
                        g_info->add_pG(i, l, temp2);
                        g_info->add_pG(l, i, temp2);
                        g_info->add_pG(i, i, -temp3);
                        g_info->add_pG(j, l, -temp4);
                        g_info->add_pG(l, j, -temp4);

                        break;

                        case 10:
                        temp1 = g_info->get_pD_val(j, l) * value * 3.0;
                        temp2 = g_info->get_pD_val(i, j) * value * 3.0;
                        temp3 = g_info->get_pD_val(j, j) * value;
                        temp4 = g_info->get_pD_val(i, l) * value * 2.0;

                        g_info->add_pG(i, j, temp1);
                        g_info->add_pG(j, i, temp1);
                        g_info->add_pG(j, l, temp2);
                        g_info->add_pG(l, j, temp2);
                        g_info->add_pG(i, l, -temp3);
                        g_info->add_pG(l, i, -temp3);
                        g_info->add_pG(j, j, -temp4);

                        break;

                        case 11:
                        temp1 = g_info->get_pD_val(k, j) * value * 3.0;
                        temp2 = g_info->get_pD_val(i, j) * value * 3.0;
                        temp3 = g_info->get_pD_val(j, j) * value;
                        temp4 = g_info->get_pD_val(i, k) * value * 2.0;

                        g_info->add_pG(i, j, temp1);
                        g_info->add_pG(j, i, temp1);
                        g_info->add_pG(k, j, temp2);
                        g_info->add_pG(j, k, temp2);
                        g_info->add_pG(i, k, -temp3);
                        g_info->add_pG(k, i, -temp3);
                        g_info->add_pG(j, j, -temp4);

                        break;

                        case 12:
                        case 13:
                        case 14:
                        temp1 = g_info->get_pD_val(k, l) * value * 4.0;
                        temp2 = g_info->get_pD_val(i, j) * value * 4.0;
                        temp3 = g_info->get_pD_val(j, l) * value;
                        temp4 = g_info->get_pD_val(i, k) * value;
                        temp5 = g_info->get_pD_val(j, k) * value;
                        temp6 = g_info->get_pD_val(i, l) * value;

                        g_info->add_pG(i, j, temp1);
                        g_info->add_pG(j, i, temp1);
                        g_info->add_pG(k, l, temp2);
                        g_info->add_pG(l, k, temp2);
                        g_info->add_pG(i, k, -temp3);
                        g_info->add_pG(k, i, -temp3);
                        g_info->add_pG(j, l, -temp4);
                        g_info->add_pG(l, j, -temp4);
                        g_info->add_pG(i, l, -temp5);
                        g_info->add_pG(l, i, -temp5);
                        g_info->add_pG(j, k, -temp6);
                        g_info->add_pG(k, j, -temp6);

                        break;
                    };
                }
            }

            g_info->unlock_mutex();

            return None;
        }


        template <typename Archive>
        void serialize(const Archive& ar) {
            ar;
        }
};

typedef madness::WorldContainer<Pair,G_MAT> mad_G;

}}

#endif /* end of have_madness */

#endif /* Header guard */


