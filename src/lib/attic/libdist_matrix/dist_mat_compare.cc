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
