#ifndef _psi_src_lib_libmints_eri_h
#define _psi_src_lib_libmints_eri_h

#include <libint/libint.h>
#include <libderiv/libderiv.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class BasisSet;
class GaussianShell;
class TwoBodyAOInt;
class IntegralFactory;
class SphericalTransform;
class Fjt;
class AOShellCombinationsIterator;
class CorrelationFactor;

/**
  * \ingroup MINTS
  * Structure to hold precomputed shell pair information
  */
typedef struct ShellPair_typ {
    //! Shells for this information.
    int i, j;
    //! Matrix over primitives with x, y, z coordinate of average Gaussian
    double *** restrict P;
    //! Distance between shell i and shell j centers
    double AB[3];
    //! Distance between P and shell i center
    double *** restrict PA;
    //! Distance between P and shell j center
    double *** restrict PB;
    //! Array of alphas for both centers
    double * restrict ai, * restrict aj;
    //! Array of the gammas (ai + aj)
    double ** restrict gamma;
    //! Contraction coefficients
    double * restrict ci, * restrict cj;
    //! Overlap between primitives on i and j
    double ** restrict overlap;
} ShellPair;

/*! \ingroup MINTS
 *  \class ERI
 *  \brief Capable of computing two-electron repulsion integrals.
 */
class TwoElectronInt : public TwoBodyAOInt
{
protected:
    //! Libint object.
    Libint_t libint_;
    //! Libderiv object
    Libderiv_t libderiv_;

    //! Maximum cartesian class size.
    int max_cart_;

    int screen_; //to screen or not to screen, that is the question
    double schwarz2_; //square of schwarz cutoff value;
    double *schwarz_norm_;

    //! Computes the fundamental
    Fjt *fjt_;

    //! Computes the ERIs between four shells.
    void compute_quartet(int, int, int, int);

    //! Computes the ERI derivatives between four shells.
    // TODO: Add integral derivative functor for erf integrals (and others)
    void compute_quartet_deriv1(int, int, int, int);

    //! Form shell pair information. Must be smart enough to handle arbitrary basis sets
    void init_shell_pairs12();
    void init_shell_pairs34();

    //! Free shell pair information
    void free_shell_pairs12();
    void free_shell_pairs34();

    //! Should we use shell pair information?
    bool use_shell_pairs_;

    //! Stack memory pointer, used in init_shell_pairs, freed in destructor
    double *stack12_, *stack34_;

    //! Shell pair information
    ShellPair **pairs12_, **pairs34_;

    //! Evaluates how much memory (in doubles) is needed to store shell pair data
    size_t memory_to_store_shell_pairs(const boost::shared_ptr<BasisSet>&, const boost::shared_ptr<BasisSet>&);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;

    //! Were the indices permuted?
    bool p13p24_, p12_, p34_;


public:
    //! Constructor. Use an IntegralFactory to create this object.
    TwoElectronInt(const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);

    virtual ~TwoElectronInt();

    //! Performs integral screening calculations
    void form_sieve();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    void compute_shell(const AOShellCombinationsIterator&);

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual void compute_shell(int, int, int, int);

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    virtual void compute_shell_deriv1(int, int, int, int);

    /// Get a derivative integral from the buffer
    virtual double get_derivative_integral(int center, int xyz, size_t index);

    //! Determine if a shell is zero based on schwarz sieve
    //Case No Sieve: false
    //Case Sieve, non-negligible integrals: false
    //Case Sieve, negligible integrals: true
    int shell_is_zero(int,int,int,int);
};

class ERI : public TwoElectronInt
{
public:
    ERI(const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);
    virtual ~ERI();
};

class F12 : public TwoElectronInt
{
public:
    F12(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);
    virtual ~F12();
};

class F12Squared : public TwoElectronInt
{
public:
    F12Squared(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);
    virtual ~F12Squared();
};

class F12G12 : public TwoElectronInt
{
public:
    F12G12(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);
    virtual ~F12G12();
};

class F12DoubleCommutator : public TwoElectronInt
{
public:
    F12DoubleCommutator(boost::shared_ptr<CorrelationFactor> cf, const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);
    virtual ~F12DoubleCommutator();
};

class ErfERI : public TwoElectronInt
{
public:
    ErfERI(double omega, const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);
    virtual ~ErfERI();

    void setOmega(double omega);
};

class ErfComplementERI : public TwoElectronInt
{
public:
    ErfComplementERI(double omega, const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);
    virtual ~ErfComplementERI();

    void setOmega(double omega);
};

}

#endif
