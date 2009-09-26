#ifndef _psi_src_lib_libmints_sobasis_h_
#define _psi_src_lib_libmints_sobasis_h_

/*!
    \file libmints/sobasis.h
    \ingroup MINTS
*/

#include <vector>

namespace psi {
    
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

class SOTransformIter
{
private:
    SOTransformShell *trans_;
    int i_;
    
public:
    SOTransformIter(SOTransformShell* trans) { trans_ = trans; i_ = 0; }
    
    void first() { i_ = 0; }
    void next()  { i_++;   }
    bool is_done() { return i_ < trans_->nfunc() ? true : false; }
    
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

}

#endif
