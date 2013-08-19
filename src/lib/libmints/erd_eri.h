#ifndef psi4_libmints_erd_eri_h_
#define psi4_libmints_erd_eri_h_

#include "psiconfig.h"
#ifdef HAVE_ERD

#include "libmints/twobody.h"

namespace psi{


class IntegralFactory;
class AOShellCombinationsIterator;

class ERDTwoElectronInt : public TwoBodyAOInt
{

typedef int F_INT;

protected:
    /// The list of contraction coefficients
    double *cc_;
    /// The list of renormalized contraction coefficients for center 1
    double *new_cc_1_;
    /// The list of renormalized contraction coefficients for center 2
    double *new_cc_2_;
    /// The list of renormalized contraction coefficients for center 3
    double *new_cc_3_;
    /// The list of renormalized contraction coefficients for center 4
    double *new_cc_4_;
    /// The list of exponents
    double *alpha_;
    /// The current size of the integral buffer
    size_t d_buffer_size_;
    /// The current size of the integer scratch space
    size_t i_buffer_size_;
    /// The address of the first contraction coefficient for each shell on center 1
    int *cc_shell_offsets_1_;
    /// The address of the first contraction coefficient for each shell on center 2
    int *cc_shell_offsets_2_;
    /// The address of the first contraction coefficient for each shell on center 3
    int *cc_shell_offsets_3_;
    /// The address of the first contraction coefficient for each shell on center 4
    int *cc_shell_offsets_4_;
    /// The integer scratch space
    F_INT *iscratch_;
    /// The double scratch space, which has junk at the start, and integrals at the end
    double *dscratch_;
    /// The start address in the target integral buffer
    F_INT buffer_offset_;

    void normalize_basis();
public:
    ERDTwoElectronInt(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~ERDTwoElectronInt();
    void compute_scratch_size();
    virtual size_t compute_shell(const psi::AOShellCombinationsIterator&);
    virtual size_t compute_shell(int, int, int, int);
    virtual size_t compute_shell_deriv1(int, int, int, int);
    virtual size_t compute_shell_deriv2(int, int, int, int);
};

class ERDERI : public ERDTwoElectronInt
{
public:
    ERDERI(const IntegralFactory* integral, int deriv=0, bool use_shell_pairs=false);
    virtual ~ERDERI();
};

}//Namespace
#endif // HAVE_ERD
#endif // header guard
