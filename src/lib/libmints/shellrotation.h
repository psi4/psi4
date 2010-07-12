#ifndef _psi_src_lib_libmints_shellrotation_h_
#define _psi_src_lib_libmints_shellrotation_h_

namespace psi {

class SymmetryOperation;
class IntegralFactory;

class ShellRotation
{
    int n_;
    int am_;
    double **r_;
    
    void done();

public:
    /// Initialize this ShellRotation to hold a n by n transformation.
    ShellRotation(int n);
    /// Initialize this from another ShellRotation.
    ShellRotation(const ShellRotation&);
    /// Initialize using init(...) or, if pure is nonzero, init_pure(...).
    ShellRotation(int a, SymmetryOperation&, const shared_ptr<IntegralFactory>&);
    virtual ~ShellRotation();

    /// Assign this to another shell rotation.
    ShellRotation& operator=(const ShellRotation&);
    
    /** Initialize the ShellRotation for Cartesian functions, given the
        angular momentum, a symmetry operation, and an Integral object. */
    void init(int a, SymmetryOperation&, const shared_ptr<IntegralFactory>&);

    /// Return the angular momentum.
    int am() const { return am_; }
    /// Return the number of functions in a shell.
    int dim() const { return n_; }
    
    /// Return an element of the transform matrix.
    double& operator()(int i, int j) { return r_[i][j]; }
    /// Return a row of the transform matrix.
    double* operator[](int i) { return r_[i]; }
    
    /// Returns the result of rot*this.
    ShellRotation operate(const ShellRotation&rot) const;
    /// Returns the result of rot*this*transpose(rot).
    ShellRotation transform(const ShellRotation&rot) const;
    
    /// Return the trace of the transformation.
    double trace() const;
    
    /// Print the object.
    void print() const;
};

}

#endif // _psi_src_lib_libmints_shellrotation_h_
