#ifndef _psi_src_lib_libmints_eri_h
#define _psi_src_lib_libmints_eri_h

#include <libutil/ref.h>

#include <libmints/basisset.h>
#include <libmints/gshell.h>
#include <libmints/osrecur.h>
#include <libmints/onebody.h>
#include <libmints/twobody.h>
#include <libmints/integral.h>

#include <libint/libint.h>
#include <libderiv/libderiv.h>

namespace psi {

    /**
      * \ingroup MINTS
      * Structure to hold precomputed shell pair information
      */
    typedef struct ShellPair_typ {
        //! Shells for this information.
        int i, j;
        //! Matrix over primitives with x, y, z coordinate of average Gaussian
        double ***P;
        //! Distance between shell i and shell j centers
        double AB[3];
        //! Distance between P and shell i center
        double ***PA;
        //! Distance between P and shell j center
        double ***PB;
        //! Array of alphas for both centers
        double *ai, *aj;
        //! Array of the gammas (ai + aj)
        double **gamma;
        //! Contraction coefficients
        double *ci, *cj;
        //! Overlap between primitives on i and j
        double **overlap;
    } ShellPair;

/*! \ingroup MINTS
 *  \class ERI
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class ERI : public TwoBodyInt
{
    //! Libint object.
    Libint_t libint_;
    //! Libderiv object
    Libderiv_t libderiv_;

    //! Maximum cartesian class size.
    int max_cart_;
    double **d_;
    double *denom_;
    double wval_infinity_;
    int itable_infinity_;

    int screen_; //to screen or not to screen, that is the question
    double schwarz2_; //square of schwarz cutoff value;
    double *schwarz_norm_;

    void init_fjt(int);
    void int_fjt(double *, int, double);

    //! Computes the ERIs between four shells.
    void compute_quartet(int, int, int, int);

    //! Computes the ERI derivatives between four shells.
    void compute_quartet_deriv1(int, int, int, int);

    //! Form shell pair information. Must be smart enough to handle arbitrary basis sets
    void init_shell_pairs12();
    void init_shell_pairs34();

    //! Free shell pair information
    void free_shell_pairs12();
    void free_shell_pairs34();

    //! Stack memory pointer, used in init_shell_pairs, freed in destructor
    double *stack12_, *stack34_;

    //! Shell pair information
    ShellPair **pairs12_, **pairs34_;

    //! Evaluates how much memory (in doubles) is needed to store shell pair data
    size_t memory_to_store_shell_pairs(const shared_ptr<BasisSet>&, const shared_ptr<BasisSet>&);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;

    //! Were the indices permuted?
    bool p13p24_, p12_, p34_;

public:
    //! Constructor. Use an IntegralFactory to create this object.
    ERI(shared_ptr<BasisSet>, shared_ptr<BasisSet>, shared_ptr<BasisSet>, shared_ptr<BasisSet>, int deriv=0, double schwarz = 0.0);

    ~ERI();

    //! Performs integral screening calculations
    void form_sieve();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    void compute_shell(int, int, int, int);

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    void compute_shell_deriv1(int, int, int, int);

    //! Determine if a shell is zero based on schwarz sieve
    //Case No Sieve: false
    //Case Sieve, non-negligible integrals: false
    //Case Sieve, negligible integrals: true
    int shell_is_zero(int,int,int,int);
};

}

#endif
