#ifndef _psi_src_lib_libmints_sobasis_h_
#define _psi_src_lib_libmints_sobasis_h_

#include <vector>

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi {

class BasisSet;

/*! \ingroup MINTS */
class SOTransformComponent
{
protected:
    /// Coefficient of the AO
    double coef_;
    /// AO function number
    int aofunc_;
    /// SO function number
    int sofunc_;
    /// SO function's irrep
    int irrep_;
    /// SO function in irrep
    int sofuncirrep_;

public:
    SOTransformComponent();

    /// Returns the coefficient of the AO
    double coef() const { return coef_; }
    /// Returns the AO function number
    int aofunc() const  { return aofunc_; }
    /// Returns the SO function number
    int sofunc() const  { return sofunc_; }
    /// Returns the SO function's irrep
    int irrep() const   { return irrep_; }
    /// Return the SO function number in its irrep
    int sofuncirrep() const { return sofuncirrep_; }

    void init(int aofunc, int sofunc, int sofuncirrep, int irrep, double coef);
};

/*! \ingroup MINTS */
class SOTransformShell
{
protected:
    /// The number of the AO shell from which these functions come.
    int aoshell_;
    /// Array of SOTransformComponent objects describing the transform
    std::vector<SOTransformComponent> funcs_;

public:
    SOTransformShell();

    void add_function(int irrep, double coef, int aofunc, int sofunc, int sofuncirrep);

    void aoshell(int i) { aoshell_ = i;    }
    int aoshell() const { return aoshell_; }
    int nfunc() const   { return funcs_.size();    }
    SOTransformComponent* func(int i) { return &(funcs_[i]); }
};

/*! \ingroup MINTS */
class SOTransformIter
{
private:
    SOTransformShell *trans_;
    int i_;

public:
    SOTransformIter(SOTransformShell* trans) { trans_ = trans; i_ = 0; }

    void first() { i_ = 0; }
    void next()  { i_++;   }
    bool is_done() { return i_ < trans_->nfunc() ? false : true; }

    /// Returns the coefficient of component i
    double coef() const { return trans_->func(i_)->coef(); }
    /// Returns the AO function number of component i
    int aofunc() const  { return trans_->func(i_)->aofunc(); }
    /// Returns the SO function number of component i
    int sofunc() const  { return trans_->func(i_)->sofunc(); }
    /// Returns the SO function's irrep of component i
    int irrep() const   { return trans_->func(i_)->irrep(); }
    /// Returns the SO function number in its irrep for component i
    int sofuncirrep() const { return trans_->func(i_)->sofuncirrep(); }
};

/*! \ingroup MINTS */
class SOTransform
{
protected:
    /// The SOTransformShell object for each AO
    std::vector<SOTransformShell> aoshell_;

public:
    SOTransform() {};

    /// Initialize
    void init(int nshells);

    /// Add another term to the transform.
    void add_transform(int aoshell, int irrep, int sofuncirrep, double coef, int aofunc, int sofunc);

    /// Returns the number of ao shells
    int naoshell() const { return aoshell_.size(); }

    /// Return the i'th ao shell
    SOTransformShell* aoshell(int i) { return &(aoshell_[i]); }
};

/** An SOBasis object describes the transformation from an atomic orbital basis
    to a symmetry orbital basis. */
class SOBasis
{
protected:
    boost::shared_ptr<BasisSet> basis_;

    int nshell_;
    int nirrep_;
    int *ncomp_;
    int **nfunc_;
    int *naofunc_;
    int **funcoff_;

    int *nfunc_in_irrep_;
    int *func_;
    int *irrep_;
    int *func_within_irrep_;

    SOTransform *trans_;

public:
    /// Create an SOBasis object given a BasisSet and Integral objects.
    SOBasis(const boost::shared_ptr<BasisSet>&, const boost::shared_ptr<IntegralFactory>&);
    ~SOBasis();
    /// Return the number of shells.
    int nshell() const { return nshell_; }
    /// Return the number of irreps.
    int nirrep() const { return nirrep_; }
    int ncomponent(int iirrep) const { return ncomp_[iirrep]; }
    /// Return the number of functions in the given irrep.
    int nfunction_in_irrep(int irrep) const { return nfunc_in_irrep_[irrep]; }
    /// Return the offset for the first function of the given irrep.
    int function_offset_for_irrep(int irrep) const;
    /// Return the number of functions in the given shell.
    int nfunction(int ishell) const;
    /** Return the number of functions in the AO shell that make up
        the given SO shell. */
    int naofunction(int ishell) const { return naofunc_[ishell]; }
    /// Returns the number of functions in the shell in a given irrep.
    int nfunction(int ishell, int iirrep) const;
    /** Returns the maximum number of functions in a shell (summed over all
        irreps) */
    int max_nfunction_in_shell() const;
    /** Normally, SO shell numbering starts at zero within each irrep.
        This returns an offset to make SO shell numbers unique within the
        shell. */
    int function_offset_within_shell(int ishell, int iirrep) const;

    /** Convert the SO shell number to the overall number of the first
        function within that shell. */
    int function(int ishell);

    /// Convert SO shell and function number within shell to irrep.
    int irrep(int ishell, int ifunc) const;
    /// Convert SO shell and function number to number within irrep.
    int function_within_irrep(int ishell, int ifunc) const;

    /// Return the SOTransform object for the given shell.
    const SOTransform &trans(int i) const { return trans_[i]; }

    //void print() const;
};

inline int SOBasis::function(int ishell)
{
  return func_[ishell];
}

inline int SOBasis::irrep(int ishell, int ifunc) const
{
  return irrep_[func_[ishell]+ifunc];
}

inline int SOBasis::function_offset_for_irrep(int irrep) const
{
  int r = 0;
  for (int i=0; i<irrep; i++) {
      r += nfunc_in_irrep_[i];
    }
  return r;
}

inline int SOBasis::function_within_irrep(int ishell, int ifunc) const
{
  return func_within_irrep_[func_[ishell]+ifunc];
}

inline int SOBasis::nfunction(int ishell, int iirrep) const
{
  return nfunc_[ishell][iirrep];
}

inline int SOBasis::function_offset_within_shell(int ishell, int iirrep) const
{
  return funcoff_[ishell][iirrep];
}

}

#endif
