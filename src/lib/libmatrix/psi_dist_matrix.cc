#include <cstring>
#include <libmatrix/matrixwrapper.h>

using namespace psi;
using namespace boost;

namespace psi {

#if parallel

PsiDistMatrix::PsiDistMatrix(int rows, int cols, int tile_sz)
    : MatrixWrapper()
{
    std::vector<int> dims(2);
    dims[0] = 2; dims[1] = 2;
    process_grid<2> pgrid(dims);

    psidistmat_ = boost::shared_ptr<Distributed_Matrix> (new Distributed_Matrix(pgrid, rows, cols, tile_sz, "dist_matrix"));
}

PsiDistMatrix::~PsiDistMatrix()
{ }

void PsiDistMatrix::print()
{
    psidistmat_->print();
}

void PsiDistMatrix::fill(double val)
{
    psidistmat_->fill(val);
}

void PsiDistMatrix::add(boost::shared_ptr<MatrixWrapper> rhs)
{
    psidistmat_->add(rhs->psidistmat_);
}

#endif // End of HAVE_MADNESS

} // End of namespace PSI

