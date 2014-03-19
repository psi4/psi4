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

#include "gridprop.h"
#include <libmints/mints.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <libqt/qt.h>
#include <psi4-dec.h>
#include <physconst.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

CubicScalarGrid::CubicScalarGrid(
        boost::shared_ptr<BasisSet> primary) :
        primary_(primary),
        mol_(primary->molecule()),
        options_(Process::environment.options)
{
    filepath_ = "";
    npoints_ = 0L;
    x_ = NULL;
    y_ = NULL;
    z_ = NULL;
    v_ = NULL;
    N_ = new int[3];
    D_ = new double[3];
    O_ = new double[3];

    build_grid(); // Defaults from Options
}
CubicScalarGrid::~CubicScalarGrid()
{
    if (x_) delete x_;
    if (y_) delete y_;
    if (z_) delete z_;
    if (v_) delete v_;
    delete[] N_;
    delete[] D_;
    delete[] O_;
}
void CubicScalarGrid::build_grid(const std::string filepath, int* N, double* D, double* O)
{
    filepath_ = filepath;    

    for (int k = 0; k < 3; k++) {
        N_[k] = N[k];
        O_[k] = O[k];
        D_[k] = D[k];
    }
    
    populate_grid();
}
void CubicScalarGrid::build_grid()
{
    filepath_ = options_.get_str("CUBIC_GRID_FILEPATH");

    double L[3];
    if (options_["CUBIC_GRID_OVERAGE"].size() != 3) { 
        L[0] = 2.0;
        L[1] = 2.0;
        L[2] = 2.0;
    } else {
        L[0] = options_["CUBIC_GRID_OVERAGE"][0].to_double();
        L[1] = options_["CUBIC_GRID_OVERAGE"][1].to_double();
        L[2] = options_["CUBIC_GRID_OVERAGE"][2].to_double();
    }

    double D[3];
    if (options_["CUBIC_GRID_SPACING"].size() != 3) {
        D[0] = 0.2;
        D[1] = 0.2;
        D[2] = 0.2;
    } else {
        D[0] = options_["CUBIC_GRID_SPACING"][0].to_double();
        D[1] = options_["CUBIC_GRID_SPACING"][1].to_double();
        D[2] = options_["CUBIC_GRID_SPACING"][2].to_double();
    }

    double Xmin[3];
    double Xmax[3];
    double Xdel[3];
    int N[3];
    double O[3];

    for (int k = 0; k < 3; k++) {
        Xmin[k] = Xmax[k] = mol_->xyz(0,k);
        for (int A = 0; A < mol_->natom(); A++) {
            Xmin[k] = (Xmin[k] > mol_->xyz(A,k) ? mol_->xyz(A,k) : Xmin[k]);
            Xmax[k] = (Xmax[k] < mol_->xyz(A,k) ? mol_->xyz(A,k) : Xmax[k]);
        }
        Xdel[k] = Xmax[k] - Xmin[k];
        N[k] = (int) ((Xmax[k] - Xmin[k] + 2.0 * L[k]) / (D[k]));
        if (D[k] * (double) N[k] < (Xmax[k] - Xmin[k] + 2.0 * L[k])) N[k]++;
        O[k] = Xmin[k] - (D[k] * (double) N[k] - (Xmax[k] - Xmin[k])) / 2.0;
        N_[k] = N[k];
        O_[k] = O[k];
        D_[k] = D[k];
    }

    populate_grid();
}
void CubicScalarGrid::populate_grid()
{
    if (x_) delete x_;
    if (y_) delete y_;
    if (z_) delete z_;
    if (v_) delete v_;

    npoints_ = (N_[0] + 1L) * (N_[1] + 1L) * (N_[2] + 1L);
    x_ = new double[npoints_];
    y_ = new double[npoints_];
    z_ = new double[npoints_];
    v_ = new double[npoints_];

    double epsilon = options_.get_double("CUBIC_BASIS_TOLERANCE");
    extents_ = boost::shared_ptr<BasisExtents> (new BasisExtents(primary_, epsilon));

    int max_points = options_.get_int("CUBIC_BLOCK_MAX_POINTS");
    double xyz = pow((double) max_points, 1.0/3.0);
    nxyz_ = (size_t) pow((double) max_points, 1.0/3.0);

    size_t offset = 0L;
    for (int istart = 0L; istart <= N_[0]; istart+=nxyz_) {
        int ni = (istart + nxyz_ > N_[0] ? (N_[0] + 1) - istart : nxyz_);
        for (int jstart = 0L; jstart <= N_[1]; jstart+=nxyz_) {
            int nj = (jstart + nxyz_ > N_[1] ? (N_[1] + 1) - jstart : nxyz_);
            for (int kstart = 0L; kstart <= N_[2]; kstart+=nxyz_) {
                int nk = (kstart + nxyz_ > N_[2] ? (N_[2] + 1) - kstart : nxyz_);

                double* xp = &x_[offset];
                double* yp = &y_[offset];
                double* zp = &z_[offset];
                double* vp = &v_[offset];
                
                size_t block_size = 0L;
                for (int i = istart; i < istart + ni; i++) {
                    for (int j = jstart; j < jstart + nj; j++) {
                        for (int k = kstart; k < kstart + nk; k++) {
                            x_[offset] = O_[0] + i * D_[0];
                            y_[offset] = O_[1] + j * D_[1];
                            z_[offset] = O_[2] + k * D_[2];
                            offset++;
                            block_size++;
                        }
                    }
                }
                blocks_.push_back(boost::shared_ptr<BlockOPoints>(new BlockOPoints(block_size,xp,yp,zp,vp,extents_)));
            }
        }
    }
            
    int max_functions = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        max_functions = (max_functions >= blocks_[ind]->functions_local_to_global().size() ? 
            max_functions : blocks_[ind]->functions_local_to_global().size());
    }
 
    points_ = boost::shared_ptr<RKSFunctions>(new RKSFunctions(primary_,max_points,max_functions));
    points_->set_ansatz(0);
}
void CubicScalarGrid::print_header()
{
    fprintf(outfile,"  ==> CubicScalarGrid <==\n\n");

    fprintf(outfile,"    Filepath     = %s\n", filepath_.c_str());
    fprintf(outfile,"    Total Points = %16zu\n", npoints_);
    fprintf(outfile,"    XYZ Blocking = %16zu\n", nxyz_);
    fprintf(outfile,"    X Points     = %16zu\n", N_[0] + 1L);
    fprintf(outfile,"    Y Points     = %16zu\n", N_[1] + 1L);
    fprintf(outfile,"    Z Points     = %16zu\n", N_[2] + 1L);
    fprintf(outfile,"    X Spacing    = %16.3E\n", D_[0]);
    fprintf(outfile,"    Y Spacing    = %16.3E\n", D_[1]);
    fprintf(outfile,"    Z Spacing    = %16.3E\n", D_[2]);
    fprintf(outfile,"    X Minimum    = %16.3E\n", O_[0]);
    fprintf(outfile,"    Y Minimum    = %16.3E\n", O_[1]);
    fprintf(outfile,"    Z Minimum    = %16.3E\n", O_[2]);
    fprintf(outfile,"    X Maximum    = %16.3E\n", O_[0] + D_[0] * N_[0]);
    fprintf(outfile,"    Y Maximum    = %16.3E\n", O_[1] + D_[1] * N_[1]);
    fprintf(outfile,"    Z Maximum    = %16.3E\n", O_[2] + D_[2] * N_[2]);
    fprintf(outfile,"\n");

    fflush(outfile);
}
void CubicScalarGrid::zero()
{
    ::memset(v_, '\0', sizeof(double) * npoints_);
}
void CubicScalarGrid::compute_density(boost::shared_ptr<Matrix> D, double scale)
{
    points_->set_pointers(D);
    boost::shared_ptr<Vector> rho = points_->point_value("RHO_A");
    double* rhop = rho->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_points(blocks_[ind]);
        size_t npoints = blocks_[ind]->npoints();
        C_DAXPY(npoints,scale,rhop,1,&v_[offset],1);        
        offset += npoints;
    }
}
void CubicScalarGrid::compute_orbital(boost::shared_ptr<Vector> C, double scale)
{
    throw PSIEXCEPTION("Not yet implemented"); 
}
void CubicScalarGrid::compute_basis(int function, double scale)
{
    throw PSIEXCEPTION("Not yet implemented"); 
}
void CubicScalarGrid::write_cube_file(const std::string& name, double* v)
{
    if (v == NULL) v = v_;

    // => Reorder the grid <= //

    double* v2 = new double[npoints_];
    size_t offset = 0L;
    for (int istart = 0L; istart <= N_[0]; istart+=nxyz_) {
        int ni = (istart + nxyz_ > N_[0] ? (N_[0] + 1) - istart : nxyz_);
        for (int jstart = 0L; jstart <= N_[1]; jstart+=nxyz_) {
            int nj = (jstart + nxyz_ > N_[1] ? (N_[1] + 1) - jstart : nxyz_);
            for (int kstart = 0L; kstart <= N_[2]; kstart+=nxyz_) {
                int nk = (kstart + nxyz_ > N_[2] ? (N_[2] + 1) - kstart : nxyz_);
                for (int i = istart; i < istart + ni; i++) {
                    for (int j = jstart; j < jstart + nj; j++) {
                        for (int k = kstart; k < kstart + nk; k++) {
                            size_t index = i * (N_[1] + 1L) * (N_[2] + 1L) + j * (N_[2] + 1L) + k;
                            v2[index] = v[offset];
                            offset++;
                        }
                    }
                }
            }
        }
    }

    // => Drop the grid out <= //

    std::stringstream ss;
    ss << filepath_ << "/" << name << ".cube";
    FILE* fh = fopen(ss.str().c_str(), "w");

    // Two comment lines
    fprintf(fh, "PSI4 Gaussian Cube File.\n");
    fprintf(fh, "Property: %s\n", name.c_str());
   
    // Number of atoms plus origin of data 
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", mol_->natom(), O_[0], O_[1], O_[2]);
    
    // Number of points along axis, displacement along x,y,z
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[0] + 1, D_[0], 0.0, 0.0);
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[1] + 1, 0.0, D_[1], 0.0);
    fprintf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[2] + 1, 0.0, 0.0, D_[2]);

    // Atoms of molecule (Z, Q?, x, y, z)
    for (int A = 0; A < mol_->natom(); A++) {
        fprintf(fh, "%3d %10.6f %10.6f %10.6f %10.6f\n", (int) mol_->Z(A), 0.0, mol_->x(A), mol_->y(A), mol_->z(A));
    }

    // Data, striped (x, y, z)
    for (size_t ind = 0; ind < npoints_; ind++) {
        fprintf(fh, "%12.5E ", v2[ind]);
        if (ind % 6 == 5) fprintf(fh,"\n");
    }

    fclose(fh);   
}
void CubicScalarGrid::compute_density_cube(boost::shared_ptr<Matrix> D, const std::string name)
{
    zero();
    compute_density(D);
    write_cube_file(name);
}
void CubicScalarGrid::compute_orbital_cubes(boost::shared_ptr<Matrix> C, const std::string name, const std::vector<int>& indices)
{
    throw PSIEXCEPTION("Not yet implemented"); 
}
void CubicScalarGrid::compute_basis_cubes(const std::string name, const std::vector<int>& indices)
{
    
    throw PSIEXCEPTION("Not yet implemented"); 
}

//void CubicScalarGrid::compute_orbitals(boost::shared_ptr<Matrix> C, double clamp, const::std::string& label)
//{
//    if (!npoints_) throw PSIEXCEPTION("CubicScalarGrid::compute: call build_grid first");
//   
//    int na = C->colspi()[0];
//
//    points_->set_Cs(C);
//    boost::shared_ptr<Matrix> psi = points_->orbital_value("PSI_A");
//    double** psip = psi->pointer();
//
//    int max_points = points_->max_points();
//
//    boost::shared_ptr<Matrix> W(new Matrix("W", na, npoints_));
//    double** Wp = W->pointer();
//
//    size_t offset = 0L;
//    for (int ind = 0; ind < blocks_.size(); ind++) {
//        size_t npoints = blocks_[ind]->npoints();
//        points_->compute_orbitals(blocks_[ind]);
//        for (int A = 0; A < na; A++) {
//            ::memcpy(&Wp[A][offset],psip[A],sizeof(double)*npoints);
//        }
//        offset += npoints;
//    }
//
//    for (int A = 0; A < na; A++) {
//        fprintf(outfile,"    Saving %4d Orbital Voxel Partition.\n", A+1);
//        std::stringstream ss;
//        ss << label << "f" << A+1 << ".raw";
//        drop_raw(ss.str(), clamp, Wp[A]);
//    }
//}
//void CubicScalarGrid::drop_raw(const std::string& file, double clamp, double* v)
//{
//    if (!npoints_) throw PSIEXCEPTION("CubicScalarGrid::drop_raw: call build_grid first");
//
//    //ASCII Variant
//    //std::stringstream ss2;
//    //ss2 << V_->name() << ".debug";
//    //FILE* fh2 = fopen(ss2.str().c_str(), "w");
//    //for (size_t ind = 0; ind < npoints_; ind++) {
//    //    fprintf(fh2,"  %16zu %24.16E %24.16E %24.16E %24.16E\n", ind, x_[ind], y_[ind], z_[ind], v_[ind]);
//    //}
//    //fclose(fh2);
//
//    if (v == NULL) { 
//        v = v_;
//    }
//
//    double s = 1.0 / clamp;
//    double maxval = 0.0;
//
//    //double total = 0.0;
//
//    float* v2 = new float[npoints_];
//    size_t offset = 0L;
//    for (int istart = 0L; istart <= N_[0]; istart+=nxyz_) {
//        int ni = (istart + nxyz_ > N_[0] ? (N_[0] + 1) - istart : nxyz_);
//        for (int jstart = 0L; jstart <= N_[1]; jstart+=nxyz_) {
//            int nj = (jstart + nxyz_ > N_[1] ? (N_[1] + 1) - jstart : nxyz_);
//            for (int kstart = 0L; kstart <= N_[2]; kstart+=nxyz_) {
//                int nk = (kstart + nxyz_ > N_[2] ? (N_[2] + 1) - kstart : nxyz_);
//                for (int i = istart; i < istart + ni; i++) {
//                    for (int j = jstart; j < jstart + nj; j++) {
//                        for (int k = kstart; k < kstart + nk; k++) {
//                            size_t index = i * (N_[1] + 1L) * (N_[2] + 1L) + j * (N_[2] + 1L) + k;
//                            double val = v[offset];
//                            maxval = (maxval >= fabs(val) ? maxval : fabs(val));
//                            //total += val;
//                            val = (val <= clamp ? val : clamp);
//                            val = (val >= -clamp ? val : -clamp);
//                            val *= s;
//                            v2[index] = val;
//                            offset++;
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    fprintf(outfile,"    Max val = %11.3E out of %11.3E: %s.\n", maxval, clamp, (clamp >= maxval ? "No clamping" : "Clamped")); 
//
//    //total *= D_[0] * D_[1] * D_[2];
//    //fprintf(outfile,"    Integral value is %24.16E\n", total);
//
//    //Dirty Hack: I love it!
//    v2[npoints_-1L] =  0.0;
//    v2[npoints_-2L] =  0.0;
//    v2[npoints_-3L] =  1.0;
//    v2[npoints_-4L] = -1.0;
//    
//    FILE* fh = fopen(file.c_str(), "wb");
//    fwrite(v2,sizeof(float),npoints_,fh);
//    fclose(fh);
//}
//void CubicScalarGrid::drop_uvf(const std::string& file, double clamp, double* v)
//{
//    if (!npoints_) throw PSIEXCEPTION("CubicScalarGrid::drop_raw: call build_grid first");
//
//    if (v == NULL) { 
//        v = v_;
//    }
//
//    double s = 1.0 / clamp;
//    double maxval = 0.0;
//
//    float* v2 = new float[npoints_];
//    size_t offset = 0L;
//    for (int istart = 0L; istart <= N_[0]; istart+=nxyz_) {
//        int ni = (istart + nxyz_ > N_[0] ? (N_[0] + 1) - istart : nxyz_);
//        for (int jstart = 0L; jstart <= N_[1]; jstart+=nxyz_) {
//            int nj = (jstart + nxyz_ > N_[1] ? (N_[1] + 1) - jstart : nxyz_);
//            for (int kstart = 0L; kstart <= N_[2]; kstart+=nxyz_) {
//                int nk = (kstart + nxyz_ > N_[2] ? (N_[2] + 1) - kstart : nxyz_);
//                for (int i = istart; i < istart + ni; i++) {
//                    for (int j = jstart; j < jstart + nj; j++) {
//                        for (int k = kstart; k < kstart + nk; k++) {
//                            size_t index = i * (N_[1] + 1L) * (N_[2] + 1L) + j * (N_[2] + 1L) + k;
//                            double val = v[offset];
//                            maxval = (maxval >= fabs(val) ? maxval : fabs(val));
//                            val = (val <= clamp ? val : clamp);
//                            val = (val >= -clamp ? val : -clamp);
//                            val *= s;
//                            v2[index] = val;
//                            offset++;
//                        }
//                    }
//                }
//            }
//        }
//    }
//
//    fprintf(outfile,"    Max val = %11.3E out of %11.3E: %s.\n", maxval, clamp, (clamp >= maxval ? "No clamping" : "Clamped")); 
//
//    //Dirty Hack: I love it!
//    v2[npoints_-1L] =  0.0;
//    v2[npoints_-2L] =  0.0;
//    v2[npoints_-3L] =  1.0;
//    v2[npoints_-4L] = -1.0;
//    
//    FILE* fh = fopen(file.c_str(), "wb");
//    const char* header = "UVF-DATA";
//    fwrite(header,sizeof(char),8,fh);
//    const char* endian = "\0";
//    fwrite(endian,sizeof(char),1,fh);
//    unsigned long int version = 3L;
//    fwrite(&version,sizeof(unsigned long int),1,fh);
//    unsigned long int nchecksum = 0L;   
//    fwrite(&nchecksum,sizeof(unsigned long int),1,fh);
//    unsigned long int offset2 = 0L;
//    fwrite(&offset2,sizeof(unsigned long int),1,fh);
//    
//    throw PSIEXCEPTION("Not implemented");
// 
//    fclose(fh);
//}

}
