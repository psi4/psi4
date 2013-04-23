/*
 *  dist_mat_set.cc
 *  part of distributed matrix
 *
 *  Created by Ben Mintz on 12/14/11.
 *
 */

#include "dist_mat.h"

#ifdef HAVE_MADNESS

namespace psi {

bool Distributed_Matrix::operator ==(const Distributed_Matrix &comp) const
{
    if (this->nelements_ != comp.nelements_) return false;
    else if (this->nrows_ != comp.nrows_) return false;
    else if (this->ncols_ != comp.ncols_) return false;
    else if (this->tile_sz_ != comp.tile_sz_) return false;
    else if (this->pgrid_ != comp.pgrid_) return false;
    else if (this->ntiles_ != comp.ntiles_) return false;
    else if (this->tile_nrows_ != comp.tile_nrows_) return false;
    else if (this->tile_ncols_ != comp.tile_ncols_) return false;
    else return true;
}

bool Distributed_Matrix::operator ==(const Distributed_Matrix *comp) const
{
    if (this->nelements_ != comp->nelements_) return false;
    else if (this->nrows_ != comp->nrows_) return false;
    else if (this->ncols_ != comp->ncols_) return false;
    else if (this->tile_sz_ != comp->tile_sz_) return false;
    else if (this->pgrid_ != comp->pgrid_) return false;
    else if (this->ntiles_ != comp->ntiles_) return false;
    else if (this->tile_nrows_ != comp->tile_nrows_) return false;
    else if (this->tile_ncols_ != comp->tile_ncols_) return false;
    else return true;
}

bool Distributed_Matrix::operator ==(boost::shared_ptr<Distributed_Matrix> comp) const
{
    if (this->nelements_ != comp->nelements_) return false;
    else if (this->nrows_ != comp->nrows_) return false;
    else if (this->ncols_ != comp->ncols_) return false;
    else if (this->tile_sz_ != comp->tile_sz_) return false;
    else if (this->pgrid_ != comp->pgrid_) return false;
    else if (this->ntiles_ != comp->ntiles_) return false;
    else if (this->tile_nrows_ != comp->tile_nrows_) return false;
    else if (this->tile_ncols_ != comp->tile_ncols_) return false;
    else return true;
}


bool Distributed_Matrix::operator !=(const Distributed_Matrix &rhs) const
{
    if (*this == rhs) return false;
    else return true;
}

bool Distributed_Matrix::operator !=(const Distributed_Matrix *rhs) const
{
    if (*this == rhs) return false;
    else return true;
}

} // End of namespace psi


#endif
