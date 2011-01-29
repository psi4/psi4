#ifndef _psi_src_lib_libmints_eri_h
#define _psi_src_lib_libmints_eri_h

#include <libint/libint.h>
#include <libderiv/libderiv.h>

namespace psi {

    class BasisSet;
    class GaussianShell;
    class TwoBodyAOInt;
    class IntegralFactory;
    class SphericalTransform;
    class SimpleMatrix;
    class Fjt;

/**
 * \ingroup MINTS
 * Functor for fundamental ERIs
 * Provides the effective contraction-aware incomplete gamma function
 * Inlined by the compiler
 */
class ERIFundamentalFunctor
{
public:
    /// The fundamental ERI integral (effective contracted incomplete gamma function)
    virtual void operator()(Libint_t &libint, Fjt* fjt, int nprim, double coef1, int max_am, double PQ2, double rho);
}; 
/**
 * \ingroup MINTS
 * Functor for fundamental erf ERIs
 * Provides the effective contraction-aware incomplete gamma function
 * Inlined by the compiler
 */
class ErfERIFundamentalFunctor : public ERIFundamentalFunctor
{
private:
    /// Weight of the 1/r operator 
    double alpha_;
    /// Weight of the erf/r operator
    double beta_;
    /// The omega parameter
    double omega_;
    /// The square of the omega parameter
    double omega2_;
public:
    /** 
    * Constructor for the ErfcERIFundamentalFunctor 
    * 
    * The coefficients provided determine the type of integral evaluated
    * \param omega The $\omega$ parameter of the error function
    * \param alpha The coefficent of the $1/r$ operator
    * \param beta  The coefficient of the $\erf(\omega r)/r$ operator
    * 
    * Integrals returned are of type:
    *  \alpha * 1/r + \beta * \erf(\omega r) / r
    * Which is equivalent to:
    *  (1 - alpha) * 1/r - beta * \erfc(\omega r) / r
    *
    * Some useful parameter pairs:
    *  - $\alpha = 0.0, \beta = 1.0$:             $\erf(\omega r)/r$  (Long range)
    *  - $\alpha = 1.0, \beta = -1.0$:            $\erfc(\omega r)/r$ (Short range)
    *  - $\alpha = 1.0, \beta = \sqrt{2.0} -1.0$:                     (MHG's MOS-MP2 integrals)
    *  - $\alpha = 1.0, \beta = 0.0$:             $1/r$               (ERI)
    */ 
    ErfERIFundamentalFunctor(double omega, double alpha, double beta) : omega_(omega), omega2_(omega*omega), alpha_(alpha), beta_(beta) {}
    // Set the adjustable omega parameter (see above) 
    void set_omega(double omega) { omega_ = omega; omega2_ = omega*omega;}
    // Set the adjustable alpha parameter (see above) 
    void set_alpha(double alpha) { alpha_ = alpha; }
    // Set the adjustable beta parameter (see above) 
    void set_beta(double beta) { beta_ = beta; }
    /// The fundamental ErfcERI integral (effective contracted incomplete gamma function)
    virtual void operator()(Libint_t &libint, Fjt* fjt, int nprim, double coef1, int max_am, double PQ2, double rho);
}; 

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
class ERI : public TwoBodyAOInt
{
protected:
    //! Libint object.
    Libint_t libint_;
    //! Libderiv object
    Libderiv_t libderiv_;

    //! ERIFundamentalFunctor
    ERIFundamentalFunctor *eri_functor_; 

    //! Maximum cartesian class size.
    int max_cart_;

    int screen_; //to screen or not to screen, that is the question
    double schwarz2_; //square of schwarz cutoff value;
    double *schwarz_norm_;

    //! Fj(T)
    Fjt *fjt_;

    //! Computes the ERIs between four shells.
    template <typename FundamentalFunctor>
    void compute_quartet(int, int, int, int, FundamentalFunctor &functor);

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
    size_t memory_to_store_shell_pairs(const shared_ptr<BasisSet>&, const shared_ptr<BasisSet>&);

    //! Original shell index requested
    int osh1_, osh2_, osh3_, osh4_;

    //! Were the indices permuted?
    bool p13p24_, p12_, p34_;

public:
    //! Constructor. Use an IntegralFactory to create this object.
    ERI(const IntegralFactory* integral, int deriv=0, double schwarz = 0.0);

    virtual ~ERI();

    //! Performs integral screening calculations
    void form_sieve();

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    void compute_shell(const ShellCombinationsIterator&);

    /// Compute ERIs between 4 shells. Result is stored in buffer.
    virtual void compute_shell(int, int, int, int);

    /// Compute ERI derivatives between 4 shells. Result is stored in buffer.
    virtual void compute_shell_deriv1(int, int, int, int);

    //! Determine if a shell is zero based on schwarz sieve
    //Case No Sieve: false
    //Case Sieve, non-negligible integrals: false
    //Case Sieve, negligible integrals: true
    int shell_is_zero(int,int,int,int);
};

/** \ingroup MINTS
* ErfERI - two-electron integrals involving error functions 
*
* Integrals returned are of type:
*  \alpha * 1/r + \beta * \erf(\omega r) / r
* Which is equivalent to:
*  (1 - alpha) * 1/r - beta * \erfc(\omega r) / r
*
* Some useful parameter pairs:
*  - $\alpha = 0.0, \beta = 1.0$:             $\erf(\omega r)/r$  (Long range)
*  - $\alpha = 1.0, \beta = -1.0$:            $\erfc(\omega r)/r$ (Short range)
*  - $\alpha = 1.0, \beta = \sqrt{2.0} -1.0$:                     (MHG's MOS-MP2 integrals)
*  - $\alpha = 1.0, \beta = 0.0$:             $1/r$               (ERI)
*/ 
class ErfERI : public ERI {

public:
    //! Constructor. Use an IntegralFactory to create this object.
    ErfERI(const IntegralFactory* integral, double omega, double alpha, double beta, int deriv=0, double schwarz = 0.0);

    virtual ~ErfERI();

    /// Set omega
    void set_omega(double omega) { (static_cast<ErfERIFundamentalFunctor*> (eri_functor_))->set_omega(omega); }

    /// Set alpha
    void set_alpha(double alpha) { (static_cast<ErfERIFundamentalFunctor*> (eri_functor_))->set_alpha(alpha); }

    /// Set beta
    void set_beta(double beta) { (static_cast<ErfERIFundamentalFunctor*> (eri_functor_))->set_beta(beta); }

};

}

#endif
