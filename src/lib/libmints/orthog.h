#ifndef ORTHOG_H
#define ORTHOG_H

#include "typedefs.h"

namespace psi {

class OverlapOrthog
{
public:
    /// An enum for the types of orthogonalization.
    enum OrthogMethod { Symmetric, Canonical, GramSchmidt };

private:
    int debug_;

    Dimension dim_;
    Dimension orthog_dim_;

    /** The tolerance for linearly independent basis functions.
     * The interpretation depends on the orthogonalization
     * method.
     */
    double lindep_tol_;
    /// The number of linearly dependent functions
    int nlindep_;
    /// The orthogonalization method
    OrthogMethod orthog_method_;
    // The orthogonalization matrices
    SharedMatrix orthog_trans_;
    SharedMatrix orthog_trans_inverse_;

    /// @{
    /** The minimum and maximum residual from the
     * orthogonalizationprocedure. The interpretation depends
     * on the method used. For symmetry and canonical, these
     * are the min and max overlap eigenvalues. These are the
     * residuals for the basis functions that actually end
     * up being used.
     */
    double min_orthog_res_;
    double max_orthog_res_;
    /// @}

    void compute_overlap_eig(SharedMatrix overlap_eigvec,
                             SharedVector isqrt_eigval,
                             SharedVector sqrt_eigval);
    void compute_symmetric_orthog();
    void compute_canonical_orthog();
    void compute_gs_orthog();
    void compute_orthog_trans();

    SharedMatrix overlap_;

public:
    OverlapOrthog(OrthogMethod method,
                  SharedMatrix overlap,
                  double lindep_tolerance,
                  int debug = 0);

    double min_orthog_res() const { return min_orthog_res_; }
    double max_orthog_res() const { return max_orthog_res_; }

    OrthogMethod orthog_method() const { return orthog_method_; }

    double lindep_tol() const { return lindep_tol_; }

    /** Returns a matrix which does the requested transform
        from a basis to an orthogonal basis.  This could be
        either the symmetric or canonical orthogonalization
        matrix.  The row dimension is the basisdimension and
        the column dimension is orthogonal basis dimension.An
        operator \f$O\f$ in the orthogonal basis is given by
        \f$X OX^T\f$ where \f$X\f$ is the return value of this
        function.
    */
    SharedMatrix basis_to_orthog_basis();

    /** Returns the inverse of the transformation returned by
     * basis_to_orthog_basis().
     */
    SharedMatrix basis_to_orthog_basis_inverse();

    /// Return an $S^{-1}$.
    SharedMatrix overlap_inverse();

    Dimension dim();
    Dimension orthog_dim();

    int nlindep();
};

} // namespace psi

#endif // ORTHOG_H
