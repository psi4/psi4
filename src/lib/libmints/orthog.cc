#include <cstdlib>
#include <cmath>
#include <iterator>

#include "matrix.h"
#include "vector.h"
#include "orthog.h"

namespace psi {

namespace {

struct max_abs {
  bool operator() (const double& i, const double& j)
    { return std::abs(i)<std::abs(j); }
};

}

OverlapOrthog::OverlapOrthog(OrthogMethod method,
                             SharedMatrix overlap,
                             double lindep_tolerance,
                             int debug)
    : nlindep_(0), min_orthog_res_(0.0), max_orthog_res_(0.0)
{
    orthog_method_ = method;
    overlap_ = overlap;
    lindep_tol_ = lindep_tolerance;
    debug_ = debug;
    dim_ = overlap_->rowspi();
    orthog_dim_.init("Orthogonal Dimension", dim_.n());
}

void OverlapOrthog::compute_overlap_eig(Matrix& overlap_eigvec,
                                        Vector& isqrt_eigval,
                                        Vector& sqrt_eigval)
{
    SharedMatrix U(new Matrix("U", overlap_->rowspi(), overlap_->colspi()));
    SharedVector m(new Vector(overlap_->colspi()));

    overlap_->diagonalize(U, m);

    double maxabs = *std::max_element(m->begin(), m->end(), max_abs());
    double s_tol = lindep_tol_ * maxabs;
    double minabs = maxabs;

    std::vector<double> m_sqrt(dim_.sum());
    std::vector<double> m_isqrt(dim_.sum());
    std::vector<int> m_index(dim_.sum());
    std::vector<int> nfunc(dim_.n());
    int nfunctotal = 0;

    nlindep_ = 0;
    for (int h=0; h<m->nirrep(); ++h) {
        for (Vector::iterator iter=m->begin_irrep(h);
             iter != m->end_irrep(h); ++iter) {
            if (*iter > s_tol) {
                if (*iter < minabs)
                    minabs = *iter;
                m_sqrt[nfunctotal] = sqrt(*iter);
                m_isqrt[nfunctotal] = 1.0/m_sqrt[nfunctotal];
                m_index[nfunctotal] = std::distance(m->begin(), iter);
                nfunc[h]++;
                nfunctotal++;
            }
            else if (orthog_method_ == Symmetric) {
                m_sqrt[nfunctotal] = 0.0;
                m_isqrt[nfunctotal] = 0.0;
                m_index[nfunctotal] = std::distance(m->begin(), iter);;
                nfunc[h]++;
                nfunctotal++;
                nlindep_++;
            }
            else
                nlindep_++;
        }
    }

    if (nlindep_ > 0 && orthog_method_ == Symmetric) {
        fprintf(outfile, "    WARNING: %d basis function%s ignored in symmetric orthogonalization.\n", nlindep_, (dim_.sum()-orthog_dim_.sum()>1)?"s":"");
    }

    if (orthog_method_ == Symmetric) {
        orthog_dim_.init("ortho basis (symmetric)", m->nirrep());
        orthog_dim_ = m->dimpi();
    }
    else {
        orthog_dim_.init("ortho basis (canonical)", m->nirrep());
        orthog_dim_ = nfunc;
    }

    overlap_eigvec.init(dim_, orthog_dim_, "Overlap Eigenvectors");
    if (orthog_method_ == Symmetric)
        overlap_eigvec.copy(U);
    else {
        int jfunc=0;
        // Copy over the vectors we need
        for (int h=0; h<m->nirrep(); ++h) {
            for (int j=0; j<nfunc[h]; ++j) {
                for (int i=0; i<dim_[h]; ++i) {
                    overlap_eigvec.set(h, i, j, U->get(h, i, m_index[jfunc]));
                }
                jfunc++;
            }
        }
    }

    sqrt_eigval.init(orthog_dim_);
    isqrt_eigval.init(orthog_dim_);

    std::copy(m_sqrt.begin(), m_sqrt.end(), sqrt_eigval.begin());
    std::copy(m_isqrt.begin(), m_isqrt.end(), isqrt_eigval.begin());

    max_orthog_res_ = maxabs;
    min_orthog_res_ = minabs;

    if (debug_ > 1) {
        overlap_->print();
        overlap_eigvec.print();
        isqrt_eigval.print();
        sqrt_eigval.print();
        fflush(outfile);
    }
}

void OverlapOrthog::compute_symmetric_orthog()
{
    Matrix overlap_eigvec;
    Vector overlap_isqrt_eigval;
    Vector overlap_sqrt_eigval;

    compute_overlap_eig(overlap_eigvec,
                        overlap_isqrt_eigval,
                        overlap_sqrt_eigval);

    SharedMatrix overlap_isqrt_eigval_mat(new Matrix(orthog_dim_, orthog_dim_));
    overlap_isqrt_eigval_mat->set_diagonal(overlap_isqrt_eigval);
    SharedMatrix overlap_sqrt_eigval_mat(new Matrix(orthog_dim_, orthog_dim_));
    overlap_sqrt_eigval_mat->set_diagonal(overlap_sqrt_eigval);

    orthog_trans_ = SharedMatrix(new Matrix("Orthogonal Transformation", dim_, dim_));
    orthog_trans_->transform(*overlap_isqrt_eigval_mat.get(), overlap_eigvec);
    orthog_trans_inverse_ = SharedMatrix(new Matrix("Orthogonal Inverse Transformation", dim_, dim_));
    orthog_trans_inverse_->transform(*overlap_sqrt_eigval_mat.get(), overlap_eigvec);
}

void OverlapOrthog::compute_canonical_orthog()
{
    Matrix overlap_eigvec;
    Vector overlap_isqrt_eigval;
    Vector overlap_sqrt_eigval;

    compute_overlap_eig(overlap_eigvec,
                        overlap_isqrt_eigval,
                        overlap_sqrt_eigval);

    SharedMatrix overlap_isqrt_eigval_mat(new Matrix(orthog_dim_, orthog_dim_));
    overlap_isqrt_eigval_mat->set_diagonal(overlap_isqrt_eigval);
    SharedMatrix overlap_sqrt_eigval_mat(new Matrix(orthog_dim_, orthog_dim_));
    overlap_sqrt_eigval_mat->set_diagonal(overlap_sqrt_eigval);

    orthog_trans_ = SharedMatrix(new Matrix("Orthogonal Transformation", orthog_dim_, dim_));
    orthog_trans_->gemm(false, true, 1.0, overlap_isqrt_eigval_mat, overlap_eigvec, 0.0);
    orthog_trans_inverse_ = SharedMatrix(new Matrix("Orthogonal Inverse Transformation", dim_, orthog_dim_));
    orthog_trans_inverse_->gemm(false, false, 1.0, overlap_eigvec, overlap_sqrt_eigval_mat, 0.0);
}

void OverlapOrthog::compute_gs_orthog()
{
    throw NotImplementedException("OverlapOrthog::compute_gs_orthog");
}

void OverlapOrthog::compute_orthog_trans()
{
    switch(orthog_method_) {
    case GramSchmidt:
        fprintf(outfile, "    Using Gram-Schmidt orthogonalization.\n");
        compute_gs_orthog();
        break;
    case Symmetric:
        fprintf(outfile, "    Using symmetric orthogonalization.\n");
        compute_symmetric_orthog();
        break;
    case Canonical:
        fprintf(outfile, "    Using canonical orthogonalization.\n");
        compute_canonical_orthog();
        break;
    default:
        throw PSIEXCEPTION("OverlapOrthog::compute_orthog_tarns: bad value.");
    }
}

SharedMatrix OverlapOrthog::basis_to_orthog_basis()
{
    if (!orthog_trans_)
        compute_orthog_trans();

    return orthog_trans_;
}

SharedMatrix OverlapOrthog::basis_to_orthog_basis_inverse()
{
    if (!orthog_trans_inverse_)
        compute_orthog_trans();

    return orthog_trans_inverse_;
}

Dimension OverlapOrthog::dim()
{
    return dim_;
}

Dimension OverlapOrthog::orthog_dim()
{
    if (!orthog_trans_)
        compute_orthog_trans();

    return orthog_dim_;
}

int OverlapOrthog::nlindep()
{
    if (!orthog_trans_)
        compute_orthog_trans();

    return nlindep_;
}

}
