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

using namespace psi;
using namespace std;

namespace psi {

// anonymous namespace, only visible in this file.
//namespace {
//string to_string(const int val)
//{
//    stringstream strm;
//    strm <<  val;
//    return strm.str();
//}
//}

madness::Void Distributed_Matrix::print_tile(const int &ti, const int &tj,
                                        const madness::Tensor<double> &tile) const
{

    int tij = ti * tile_ncols_ + tj;
    print_tile_tij(tij, tile);

}

madness::Void Distributed_Matrix::print_tile_tij(const int &tij, const madness::Tensor<double> &tile) const
{
    print_mutex_->lock();

    int nrows = tile.dim(0);
    int ncols = tile.dim(1);

    std::string fname;
    if (name_.size() != 0) fname = name_ + ": Owner " + to_string(owner(tij)) +
            ": Tile " + to_string(tij);
    else fname = ": Owner " + to_string(owner(tij)) + ": Tile " + to_string(tij);

    psi::fprintf(outfile, "\n  ## %s ##\n", fname.c_str());

    if (tile.size() == 0) {
        psi::fprintf(outfile, "\n\t## %s ## (empty)\n", fname.c_str());
    }
    else {
        const int print_ncol = 5;
        int num_frames = int(ncols/print_ncol);
        int num_frames_rem = ncols%print_ncol; //adding one for changing 0->1 start
        int num_frame_counter = 0;
        //for each frame
        for(num_frame_counter=0;num_frame_counter<num_frames;num_frame_counter++){
            psi::fprintf(outfile,"\n");
            for(int j=print_ncol*num_frame_counter+1;j<print_ncol*num_frame_counter+print_ncol+1;j++){
                if(j==print_ncol*num_frame_counter+1){ psi::fprintf(outfile,"%18d",j); }
                else{ psi::fprintf(outfile,"        %5d",j); }
            }
            psi::fprintf(outfile,"\n\n");

            for(int k=1; k<=nrows; ++k){
                for(int j=print_ncol*num_frame_counter+1;j<print_ncol*num_frame_counter+print_ncol+2;j++){
                    if(j==print_ncol*num_frame_counter+1){ psi::fprintf(outfile,"%5d",k);}
                    else{ psi::fprintf(outfile," %12.7f", tile(k-1, j-2)); } //[(k-1)*ncols + (j-2)]); }
                }
                psi::fprintf(outfile,"\n");
            }
        }

        // ALREADY DID THE FULL FRAMES BY THIS POINT
        // NEED TO TAKE CARE OF THE REMAINDER
        if(num_frames_rem != 0){
            psi::fprintf(outfile,"\n");
            for(int j=print_ncol*num_frame_counter+1;j<=ncols;j++){
                if(j==print_ncol*num_frame_counter+1){ psi::fprintf(outfile,"%18d",j); }
                else{ psi::fprintf(outfile,"        %5d",j); }
            }
            psi::fprintf(outfile,"\n\n");

            for(int k=1; k<=nrows; ++k){
                for(int j=print_ncol*num_frame_counter+1;j<ncols+2;j++){
                    if(j==print_ncol*num_frame_counter+1){ psi::fprintf(outfile,"%5d",k); }
                    else{ psi::fprintf(outfile," %12.7f", tile(k-1, j-2)); } //tile_.get()[(k-1)*ncols + (j-2)]); }
                }
                psi::fprintf(outfile,"\n");
            }
        }
        psi::fprintf(outfile,"\n\n");        }

    fflush(outfile);
    print_mutex_->unlock();

}


madness::Void Distributed_Matrix::print(const std::string str) const
{
    WorldComm->sync();
    if (me_ == 0) {
        if (nelements_) {
            psi::fprintf(outfile, "\n\n%s\n", str.c_str());
            for (int tij=0; tij < ntiles_; tij++) {
                madness::Future<madness::Tensor<double> > tile = task(owner(tij), &Distributed_Matrix::get_remote_tile_tij, tij, 0, 0);
                task(me_, &Distributed_Matrix::print_tile_tij, tij, tile);
            }
        }
    }
    WorldComm->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::print(const int &ti, const int &tj) const
{
    WorldComm->sync();
    if (me_ == 0) {
        madness::Future<madness::Tensor<double> > tile = task(owner(ti,tj), &Distributed_Matrix::get_tile, ti, tj);
        task(me_, &Distributed_Matrix::print_tile, ti, tj, tile);
    }
    WorldComm->sync();
    return madness::None;
}

madness::Void Distributed_Matrix::print(const int &tij) const
{
    WorldComm->sync();
    if (me_ == 1) {
        madness::Future<madness::Tensor<double> > tile = task(owner(tij), &Distributed_Matrix::get_tile_tij, tij);
        task(me_, &Distributed_Matrix::print_tile_tij, tij, tile);
    }
    WorldComm->sync();
    return madness::None;

}

void Distributed_Matrix::print_matrix_info()
{
    if (me_ == 0) {
        print_mutex_->lock();
        std::cout << "global nrows = " << nrows_ << std::endl;
        std::cout << "global ncols = " << ncols_ << std::endl;
        std::cout << "global_nelements = " << nelements_ << std::endl;
        std::cout << "global_tile_sz = " << tile_sz_ << std::endl;
        std::cout << "global tile nrows = " << tile_nrows_ << std::endl;
        std::cout << "global tile ncols = " << tile_ncols_ << std::endl;
        std::cout << "global ntiles = " << ntiles_ << std::endl;

        for (int ti=0; ti < tile_nrows_; ti++) {
            for (int tj=0; tj < tile_ncols_; tj++) {
                std::cout << "owner(" << ti << ", " << tj << ") = " << owner(ti,tj) << std::endl;
            }
        }
        for (int tij=0; tij < ntiles_; tij++) {
            std::cout << "owner(" << tij << ") = " << owner(tij) << std::endl;
        }

        print_mutex_->unlock();
    }
    print_mutex_->lock();

    for (int ti=0, tij=0; ti < tile_nrows_; ti++) {
        for (int tj=0; tj < tile_ncols_; tj++, tij++) {
            if (me_ == owner(tij)) {
                std::cout << "proc " << me_ << ": global_local_map[" << tij << "] = " <<
                             local(tij) << std::endl;
                //                    std::cout << "proc " << me_ << ": Tile[" << ti << "][" << tj << "] = " << data_[local(tij)].dim(0) << " x " << data_[local(tij)].dim(1) << std::endl;
            }
            WorldComm->sync();
        }
    }
    print_mutex_->unlock();

}


} // End of namespace psi

#endif
