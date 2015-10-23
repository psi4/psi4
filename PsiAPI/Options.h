/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */
#ifndef PSIAPI_OPTIONS_H_
#define PSIAPI_OPTIONS_H_
#include <boost/any.hpp>
#include "Utils/Exception.h"
namespace PsiAPI{

/** \brief A class to hold all of the options
 *
 *  Given the myriad of methods it is no surprise that there is almost a
 *  countable infinity of options available in each package.  This class
 *  provides a portable, uniform way of accessing these options.
 *  We define an option as a tag and a value, the latter can be any type.
 *  The tag is meant to be descriptive of what the option contains and
 *  is set by the method when it is initially passed an options object.
 *
 *  Typical usage cases for this class are:
 *  \code
 *  //Get an Options object (there's one in your Reference)
 *  Reference MyRef;
 *  Options& MyOptions=MyRef.Options();
 *
 *  //See if anyone set the integral thresholds
 *  bool HasThresh=MyOptions.HasOption("INTS_TOLERANCE");
 *
 *  //Of course they did, because you set it when the driver called
 *  //your method's InitializeOptions() member function
 *  double Tol=0.0;
 *  MyOptions.GetOption("INTS_TOLERANCE",Tol);
 *
 *  //Set your super-secret option to 4
 *  MyOptions.SetOption("MY_SUPER_SECRET_OPTION",4);
 *  \endcode
 *
 *
 *  If this class is going to work with many programs there has to be
 *  a standardization of the option names.  Currently we adopt those of
 *  Psi4 (copy/pasted here from Psi4Public Wiki for convenience):

    - All option names should be all caps and separated by underscores.
    - In deciding how to arrange words in an option name, place the
    context first (e.g., MP2_AMPS_PRINT, TRIPLES_DIIS). This means PRINT
    will generally be at the end of an option name.
    - You're not welcome to add CHARGE or MULTP options. Plan to get these
    quantities from the molecule object.
    - Convergence of a method should be governed by an E_CONVERGENCE for energy
    and either a D_CONVERGENCE for density or a R_CONVERGENCE for
    residual/amplitudes. All of these should be doubles- let the input
    parser handle the flexible input format.
    - Diis should have a boolean DIIS (not do_diis, not use_diis) to turn
    on/off diis extrapolation, a DIIS_MIN_VECS and DIIS_MAX_VECS for minimum
    and maximum number of diis vectors to use, and a DIIS_START which is
    the iteration at which to start saving vectors for diis.
    - INTS (not integrals), also OEI (not oe_integrals) for one-electron
    integrals and TEI (not te_integrals) for two-electron integrals
    - Use PRINT options to indicate printing to output file. Use WRITE
    options to indicate printing to another file. This probably isn't
    entirely valid now but should be observed in future. The complement
    to WRITE is READ. PRINT, READ, and WRITE will usually be the last
    words in an option name.
    - If you have a quantity you'd like to call a cutoff, a threshold, a
    tolerance, or a convergence, consider the following guidelines in
    naming it:
      - If its value is typically greater than ~0.001, give it a name
        with CUTOFF.
      - If its value is typically less than ~0.001 and quantities being
        tested against the option are more valuable with larger values
        (e.g., integrals, occupations, eigenvectors), give it a name with
        TOLERANCE.
      - If its value is typically less than ~0.001 and quantities being
        tested against the option are more valuable with smaller values
        (e.g., energy changes, residual errors, gradients), give it a
        name with CONVERGENCE.
    - H in an option name is reserved for Hamiltonian (or hydrogen).
    Hessian should be HESS.
    - If you have an option that instructs your module to do something
    not too computationally intensive and then quit, append _EXIT to the
    option name.
    - Scaling terms (like for scs) should follow the pattern MP2_SS_SCALE
    and SAPT_OS_SCALE.
    - For level-shifting, let's try to have it governed by (double)
    LEVEL_SHIFT only and not a boolean/double combo since the procedure
    can be turned on (role of boolean) if the value (role of double) has
    changed.
    - Use AO and MO for atomic and molecular orbitals. When 'O' for
    orbitals is too obscure or would make for too short a keyword, as in
    "bool NO" for "Do use natural orbitals", use ORBS for orbitals.
     So natural orbitals are NAT_ORBS and Brueckner orbitals are
     BRUECKNER_ORBS.
    - WRITE/READ for info transfer across jobs. SAVE/RESTART for same in
    context of restart.
    - Damping should interface through option (double) DAMPING_PERCENTAGE,
    where a value of 0.0 indicates no damping.
    - Try to avoid COMPUTE or CALC in an option name. If it's a boolean
    like "opdm_compute" for "Do compute the one-particle density matrix",
    just use OPDM.
    - Properties should be governed by a PROPERTIES array for the
    root of interest or by a PROPERTIES_ALL array for all roots in
    a multi-root calc. Since no module conforms to this right now,
    use PROPERTY alone and PROP in multi-part option as PROP_ROOT,
    PROP_ALL, PROP_SYM to conform.
    - Use DF (not ri) for density-fitting and resolution-of-the-identity
    option names. Only the basis sets are staying as -RI since that's
    what EMSL uses.
    - AMPS (not amplitude, not amp) for amplitudes
    - NUM_ (not n) for number (e.g., NUM_AMPS_PRINT, MAX_NUM_VECS, NUM_THREADS)
    - Use MAXITER, CACHELEVEL, PUREAM, DERTYPE.
    - PERTURB (not pert) for perturbation
    - Use FOLLOW_ROOT for the state to be followed in geometry optimizations
    - WFN (not wavefunction)
    - Use INTS_TOLERANCE (not schwarz_cutoff)
    - TRIPLES (not trip), TRIPLETS (not trip), SINGLES (not sing), SINGLETS (not sing)
    - CONVERGENCE (not conv, not converge) and TOLERANCE (not tol)
    - For Tikhonow regularization, use TIKONOW_OMEGA, not regularizer.
    - OCC for occupied/occupation (e.g., DOCC, LOCK_OCC, OCC_TOLERANCE).
    - COND for condition and CONDITIONER for conditioner.
    - LOCAL (not localize).
    - FRAG for fragment.
    - AVG for average.
    - LEVEL (not LVL, not LEV).
    - EX for excitation.
    - VAL for valence.
    - GEOM (not geo, not geometry).
    - SYM (not symm, not symmetry).
    - FILE (unless truly multiple FILES).

 */
class Options{
   private:
      ///The actual options
      PsiMap<std::string,boost::any> Options_;
   public:
      /** \brief The call to insert another option
       *
       *  Don't let the template scare you, this call is easy. Just give
       *  the options object a string describing what the option is and
       *  a value for the Option.  The value can be any valid type: int,
       *  double, std::vector<double>, or even a user defined type.
       *
       *
       *  \param[in] Name An ASCII name for your option
       *  \param[in] Value The value of your option
       *  \param[in] T The type of the value (determined by compiler and
       *               can be ignored)
       *
       */
      template<typename T>
      void AddOption(const std::string& Name,const T& Value){
         Options_[Name]=Value;
      }

      /** \brief The call to obtain an option value
       *
       *
       *
       */
      template<typename T>
      void GetOption(const std::string& Name, T& Option){
         if(Options_[Name].type()!=typeid(T))
            PSIERROR("Types are not compatible");
         Option=boost::any_cast<T>(Options_[Name]);
      }
      ///Returns true if this options object contains an
      bool HasOption(const std::string& Name)const{
         return Options_.count(Name)==1;
      }
};



}



#endif /* PSIAPI_OPTIONS_H_ */
