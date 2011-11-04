#ifndef _psi_src_lib_libmints_sobasis_h_
#define _psi_src_lib_libmints_sobasis_h_

#include <vector>
#include <cstdio>

#include "dimension.h"

namespace boost {
template<class T>
class shared_ptr;
}

namespace psi {

extern FILE *outfile;

class BasisSet;
class IntegralFactory;
class Dimension;
class PetiteList;

/*! \ingroup MINTS */
/** SOTransformFunction describes how an AO function contributes to an SO
    function in a particular SO shell. */
class SOTransformFunction {
public:
    /// The coefficient of the AO.
    double coef;
    /// The AO function number.
    int aofunc;
    /// The SO function number.
    int sofunc;
    /// The SO function's irrep.
    int irrep;
};

class AOTransformFunction {
public:
    double coef;
    int aofunc;
    int sofunc;
    int irrep;

    AOTransformFunction(double cf, int af, int sf, int irr)
        : coef(cf), aofunc(af), sofunc(sf), irrep(irr) {}
};

/*! \ingroup MINTS */
/** SOTransformShell maintains a list of AO functions contribute to an SO
    function in a particular SO shell.  The information is stored in
    objects of type SOTransformFunction. */
class SOTransformShell {
public:
    /// The number of the AO shell from which these functions come.
    int aoshell;
    /// The number of AO/SO function pairs contributing.
    int nfunc;
    /// The array of SOTransformFunction objects describing the transform.
    SOTransformFunction *func;
    SOTransformShell();
    ~SOTransformShell();
    /// Add another function to the transform.
    void add_func(int irrep, double coef, int aofunc, int sofunc);
};

/*! \ingroup MINTS */
/** SOTransform maintains a list of AO shells that are be used
    to compute the SO.  The information is stored in objects of
    type SOTransformShell. */
class SOTransform {
public:
    int naoshell_allocated;
    /// The number of AO shells that make up this SO shell.
    int naoshell;
    /// The SOTransformShell object for each AO.
    SOTransformShell *aoshell;
    SOTransform();
    ~SOTransform();
    void set_naoshell(int n);
    /// Adds another term to the transform.
    void add_transform(int aoshell, int irrep,
                       double coef, int aofunc, int sofunc);
};

class AOTransform {
public:
    std::vector<AOTransformFunction> soshell;

    AOTransform();
    ~AOTransform();
    void add_transform(int irrep, double coef, int aofunc, int sofunc);
};

/** An SOBasis object describes the transformation from an atomic orbital basis
    to a symmetry orbital basis. */
class SOBasisSet
{
protected:
    boost::shared_ptr<BasisSet> basis_;
    const IntegralFactory* integral_;
    boost::shared_ptr<PetiteList> petite_;

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

    int *ushell_am_;

    SOTransform *sotrans_;
    AOTransform *aotrans_;

    //! vector of so shells numbers sorted in acending AM order.
    std::vector< int > sorted_so_shell_list_;

    /// Handles initializing SOBasis
    void init();

public:
    /// Create an SOBasis object given a BasisSet and Integral objects.
    SOBasisSet(const boost::shared_ptr<BasisSet>&, const boost::shared_ptr<IntegralFactory>&);
    SOBasisSet(const boost::shared_ptr<BasisSet>&, const IntegralFactory*);
    ~SOBasisSet();

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
    /// Return tehe angular momentum of the unique shell
    int am(int ishell) const { return ushell_am_[ishell]; }
    /** Return the number of functions in the AO shell that make up
        the given SO shell. */
    int naofunction(int ishell) const { return naofunc_[ishell]; }
    /// Returns the number of functions in the shell in a given irrep.
    int nfunction(int ishell, int iirrep) const;
    /** Returns the maximum number of functions in a shell (summed over all
        irreps) */
    int max_nfunction_in_shell() const;
    int *function_offset_within_shell(int shell) const { return funcoff_[shell]; }

    /** Normally, SO shell numbering starts at zero within each irrep.
        This returns an offset to make SO shell numbers unique within the
        shell. */
    int function_offset_within_shell(int ishell, int iirrep) const;

    /** Convert the SO shell number to the overall number of the first
        function within that shell. */
    int function(int ishell);

    /// Convert SO shell and function number within shell to irrep.
    int irrep(int ishell, int ifunc) const;
    /// Convert SO shell and function number within shell to irrep.
    int irrep(int ifunc) const;
    /// Convert SO shell and function number to number within irrep.
    int function_within_irrep(int ishell, int ifunc) const;
    /// Convert SO shell and function number to number within irrep.
    int function_within_irrep(int ifunc) const;

    /// Return the SOTransform object for the given shell.
    const SOTransform &sotrans(int i) const { return sotrans_[i]; }
    /// Return the AOTransform object for the given shell.
    const AOTransform &aotrans(int i) const { return aotrans_[i]; }

    /// Return the PetiteList object used in creating this SOBasis.
    const boost::shared_ptr<PetiteList> petitelist() const;

    /// Returns the dimension for each irrep.
    Dimension dimension() const;

    //! Returns the value of the sorted shell list.
    int get_so_sorted_shell(const int &i) { return sorted_so_shell_list_[i]; }
    //! Returns the vector of sorted shell list.
    std::vector<int> get_so_sorted_list() { return sorted_so_shell_list_; }

    void print(FILE *out = outfile) const;
};

inline int SOBasisSet::function(int ishell)
{
    return func_[ishell];
}

inline int SOBasisSet::irrep(int ishell, int ifunc) const
{
    return irrep_[func_[ishell]+ifunc];
}

inline int SOBasisSet::irrep(int ifunc) const
{
    return irrep_[ifunc];
}

inline int SOBasisSet::function_offset_for_irrep(int irrep) const
{
    int r = 0;
    for (int i=0; i<irrep; i++) {
        r += nfunc_in_irrep_[i];
    }
    return r;
}

inline int SOBasisSet::function_within_irrep(int ishell, int ifunc) const
{
    return func_within_irrep_[func_[ishell]+ifunc];
}

inline int SOBasisSet::function_within_irrep(int ifunc) const
{
    return func_within_irrep_[ifunc];
}

inline int SOBasisSet::nfunction(int ishell, int iirrep) const
{
    return nfunc_[ishell][iirrep];
}

inline int SOBasisSet::function_offset_within_shell(int ishell, int iirrep) const
{
    return funcoff_[ishell][iirrep];
}

}

#endif
