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

#include "vis.h"
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

double GetBSRadius(unsigned Z)
{
    // Not sure where these numbers come from...
    static const double BSRadii[55] = {
        1.000,
        1.001,                                                                                                                 1.012,
        0.825, 1.408,                                                                       1.485, 1.452, 1.397, 1.342, 1.287, 1.243,
        1.144, 1.364,                                                                       1.639, 1.716, 1.705, 1.683, 1.639, 1.595,
        1.485, 1.474, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.562, 1.650, 1.727, 1.760, 1.771, 1.749, 1.727,
        1.628, 1.606, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.639, 1.672, 1.804, 1.881, 1.892, 1.892, 1.881,
    };
    if (Z < sizeof BSRadii/sizeof BSRadii[0])
        return BSRadii[Z];
    else
        return 1.881;
};

namespace psi {
namespace dftsapt {

ASAPTVis::ASAPTVis(
        boost::shared_ptr<BasisSet> primary, 
        boost::shared_ptr<Molecule> monomer_A,
        boost::shared_ptr<Molecule> monomer_B,
        boost::shared_ptr<Matrix> Locc_A, 
        boost::shared_ptr<Matrix> Locc_B, 
        boost::shared_ptr<Matrix> Q_A, 
        boost::shared_ptr<Matrix> Q_B
        ) : 
        primary_(primary),
        monomer_A_(monomer_A),
        monomer_B_(monomer_B),
        Locc_A_(Locc_A),
        Locc_B_(Locc_B),
        Q_A_(Q_A),
        Q_B_(Q_B),
        options_(Process::environment.options)
{
    for (int i = 0; i < options_["ASAPT_TASKS"].size(); i++) {
        tasks_.insert(options_["ASAPT_TASKS"][i].to_string());
    } 
}
ASAPTVis::~ASAPTVis()
{
}
void ASAPTVis::analyze()
{
    fprintf(outfile, "  ANALYSIS:\n\n");

    monomer_A_->save_xyz("mA.xyz",false);
    monomer_B_->save_xyz("mB.xyz",false);

    summations();
    if (tasks_.count("ATOMIC1")) {
        drop_atomic_1();
    }
    if (tasks_.count("ATOMIC2")) {
        drop_atomic_2();
    }
    if (tasks_.count("ORBITAL1")) {
        drop_orbital_1();
    }
    if (tasks_.count("ORBITAL2")) {
        drop_orbital_2();
    }
    if (tasks_.count("VOXEL")) {
        drop_voxel();
    }
}
void ASAPTVis::summations()
{
    if (!vars_.count("Elst_AB"))   throw PSIEXCEPTION("Missing Elst_AB");
    if (!vars_.count("Exch_ab"))   throw PSIEXCEPTION("Missing Exch_ab");
    if (!vars_.count("IndAB_aB"))  throw PSIEXCEPTION("Missing Ind_aB");
    if (!vars_.count("IndBA_Ab"))  throw PSIEXCEPTION("Missing Ind_Ab");
    if (!vars_.count("Disp_ab"))   throw PSIEXCEPTION("Missing Disp_ab");

    vars_["Exch_AB"]  = Matrix::triplet(Q_A_,vars_["Exch_ab"],Q_B_,false,false,true);
    vars_["IndAB_AB"] = Matrix::doublet(Q_A_,vars_["IndAB_aB"],false,false);
    vars_["IndBA_AB"] = Matrix::doublet(vars_["IndBA_Ab"],Q_B_,false,true);
    vars_["Disp_AB"]  = Matrix::triplet(Q_A_,vars_["Disp_ab"],Q_B_,false,false,true);

    vars_["Elst_A"]  = vars_["Elst_AB"]->collapse(1);
    vars_["Elst_B"]  = vars_["Elst_AB"]->collapse(0);
    vars_["Exch_A"]  = vars_["Exch_AB"]->collapse(1);
    vars_["Exch_B"]  = vars_["Exch_AB"]->collapse(0);
    vars_["Exch_a"]  = vars_["Exch_ab"]->collapse(1);
    vars_["Exch_b"]  = vars_["Exch_ab"]->collapse(0);
    vars_["IndAB_A"] = vars_["IndAB_AB"]->collapse(1);
    vars_["IndAB_B"] = vars_["IndAB_AB"]->collapse(0);
    vars_["IndAB_a"] = vars_["IndAB_aB"]->collapse(1);
    vars_["IndBA_A"] = vars_["IndBA_AB"]->collapse(1);
    vars_["IndBA_B"] = vars_["IndBA_AB"]->collapse(0);
    vars_["IndBA_b"] = vars_["IndBA_Ab"]->collapse(0);
    vars_["Disp_A"]  = vars_["Disp_AB"]->collapse(1);
    vars_["Disp_B"]  = vars_["Disp_AB"]->collapse(0);
    vars_["Disp_a"]  = vars_["Disp_ab"]->collapse(1);
    vars_["Disp_b"]  = vars_["Disp_ab"]->collapse(0);

    // Make sure the names line up
    for (std::map<std::string, boost::shared_ptr<Matrix> >::iterator it = vars_.begin();
        it != vars_.end(); ++it) {
        (*it).second->set_name((*it).first);
    }
}
void ASAPTVis::drop_atomic_1()
{
    fprintf(outfile,"    Saving Order-1 Atomic Partition.\n");

    drop("Elst_A");
    drop("Elst_B");
    drop("Exch_A");
    drop("Exch_B");
    drop("IndAB_A");
    drop("IndAB_B");
    drop("IndBA_A");
    drop("IndBA_B");
    drop("Disp_A");
    drop("Disp_B");
}
void ASAPTVis::drop_atomic_2()
{
    fprintf(outfile,"    Saving Order-2 Atomic Partition.\n");

    drop("Elst_AB");
    drop("Exch_AB");
    drop("IndAB_AB");
    drop("IndBA_AB");
    drop("Disp_AB");
}
void ASAPTVis::drop_orbital_1()
{
    fprintf(outfile,"    Saving Order-1 Orbital Partition.\n");

    drop("Exch_a");
    drop("Exch_b");
    drop("IndAB_a");
    drop("IndBA_b");
    drop("Disp_a");
    drop("Disp_b");
}
void ASAPTVis::drop_orbital_2()
{
    fprintf(outfile,"    Saving Order-2 Orbital Partition.\n");

    drop("Exch_ab");
    drop("IndAB_aB");
    drop("IndBA_Ab");
    drop("Disp_ab");
}
void ASAPTVis::drop(const std::string& var) 
{
    boost::shared_ptr<Matrix> T = vars_[var];
    double** Tp = T->pointer();
    int nrow = T->rowspi()[0];
    int ncol = T->colspi()[0];
    
    std::stringstream ss;
    ss << var << ".dat";
    FILE* fh = fopen(ss.str().c_str(),"w");
    for (int i = 0; i < nrow; i++) {
        for (int j = 0; j < ncol; j++) {
            fprintf(fh,"%24.16E%s", Tp[i][j], (j + 1 == ncol ? "" : " "));
        }
        fprintf(fh,"\n");
    }
    fclose(fh);
}
void ASAPTVis::drop_voxel()
{
    fprintf(outfile,"    Saving Order-1 Voxel Partition:\n\n");

    // ==> Sizing <== // 

    int nn = Locc_A_->rowspi()[0];
    int na = Locc_A_->colspi()[0];
    int nb = Locc_B_->colspi()[0];

    // ==> Saturation Scales <== //

    double dens_scale = options_.get_double("ASAPT_DENSITY_CLAMP");
    double int_scale  = options_.get_double("ASAPT_ENERGY_CLAMP");
    double gauss_scale = options_.get_double("ASAPT_GAUSSIAN_SCALE");
    
    // ==> Temporaries <== //

    boost::shared_ptr<Matrix> TA(Locc_A_->clone());
    boost::shared_ptr<Matrix> TB(Locc_B_->clone());

    double** LAp = Locc_A_->pointer(); 
    double** LBp = Locc_B_->pointer(); 
    double** TAp = TA->pointer(); 
    double** TBp = TB->pointer(); 
    
    // ==> Floating Temporaries <== //

    double** EAp;
    double** EBp;

    boost::shared_ptr<Matrix> DA;
    boost::shared_ptr<Matrix> DB;

    boost::shared_ptr<Matrix> QA(new Matrix("QA",vars_["Elst_A"]->nrow(),5));
    boost::shared_ptr<Matrix> QB(new Matrix("QA",vars_["Elst_B"]->nrow(),5));
    double** QAp = QA->pointer();
    double** QBp = QB->pointer();
    int offsetA = 0;
    for (int A = 0; A < monomer_A_->natom(); A++) {
        if (monomer_A_->Z(A) == 0.0) continue;
        QAp[offsetA][0] = monomer_A_->x(A);
        QAp[offsetA][1] = monomer_A_->y(A);
        QAp[offsetA][2] = monomer_A_->z(A);
        double RBS = GetBSRadius(monomer_A_->true_atomic_number(A));
        QAp[offsetA][4] = gauss_scale / (RBS * RBS);
        offsetA++;
    }
    int offsetB = 0;
    for (int B = 0; B < monomer_B_->natom(); B++) {
        if (monomer_B_->Z(B) == 0.0) continue;
        QBp[offsetB][0] = monomer_B_->x(B);
        QBp[offsetB][1] = monomer_B_->y(B);
        QBp[offsetB][2] = monomer_B_->z(B);
        double RBS = GetBSRadius(monomer_B_->true_atomic_number(B));
        QBp[offsetB][4] = gauss_scale / (RBS * RBS);
        offsetB++;
    }

    // ==> The Grid: A Digital Frontier to Reshape the Chemist's Condition <== // 

    boost::shared_ptr<CubicDensityGrid> grid = boost::shared_ptr<CubicDensityGrid>(new CubicDensityGrid(primary_));
    grid->build_grid();
    grid->print_header();
    
    // ==> Elst <== //

    fprintf(outfile, "    Saving Elst Voxel Partition.\n");

    EAp = vars_["Elst_A"]->pointer();
    EBp = vars_["Elst_B"]->pointer();
    
    C_DCOPY(QA->nrow(),EAp[0],1,&QAp[0][3],5);
    C_DCOPY(QB->nrow(),EBp[0],1,&QBp[0][3],5);

    grid->zero();
    grid->compute_atomic(QA);
    grid->compute_atomic(QB);
    grid->drop_raw("Elst.raw",int_scale);
    
    // ==> Exch <== //

    fprintf(outfile, "    Saving Exch Voxel Partition.\n");

    EAp = vars_["Exch_a"]->pointer();
    EBp = vars_["Exch_b"]->pointer();

    TA->copy(Locc_A_);
    for (int a = 0; a < na; a++) {
        C_DSCAL(nn,EAp[a][0],&TAp[0][a],na);
    }
    DA = Matrix::doublet(TA,Locc_A_,false,true);

    TB->copy(Locc_B_);
    for (int b = 0; b < nb; b++) {
        C_DSCAL(nn,EBp[b][0],&TBp[0][b],nb);
    }
    DB = Matrix::doublet(TB,Locc_B_,false,true);

    DA->add(DB); 
    DA->set_name("Exch");

    grid->zero();
    grid->compute_electronic(DA);
    grid->drop_raw("Exch.raw",int_scale);
    
    // ==> IndAB <== //

    fprintf(outfile, "    Saving IndAB Voxel Partition.\n");

    EAp = vars_["IndAB_a"]->pointer();
    EBp = vars_["IndAB_B"]->pointer();
    
    C_DCOPY(QB->nrow(),EBp[0],1,&QBp[0][3],5);

    TA->copy(Locc_A_);
    for (int a = 0; a < na; a++) {
        C_DSCAL(nn,EAp[a][0],&TAp[0][a],na);
    }
    DA = Matrix::doublet(TA,Locc_A_,false,true);

    DA->set_name("IndAB");

    grid->zero();
    grid->compute_electronic(DA);
    grid->compute_atomic(QB);
    grid->drop_raw("IndAB.raw",int_scale);

    // ==> IndBA <== //

    fprintf(outfile, "    Saving IndBA Voxel Partition.\n");

    EBp = vars_["IndBA_b"]->pointer();
    EAp = vars_["IndBA_A"]->pointer();
    
    C_DCOPY(QA->nrow(),EAp[0],1,&QAp[0][3],5);

    TB->copy(Locc_B_);
    for (int b = 0; b < nb; b++) {
        C_DSCAL(nn,EBp[b][0],&TBp[0][b],nb);
    }
    DB = Matrix::doublet(TB,Locc_B_,false,true);

    DB->set_name("IndBA");

    grid->zero();
    grid->compute_electronic(DB);
    grid->compute_atomic(QA);
    grid->drop_raw("IndBA.raw",int_scale);

    // ==> Disp <== //

    fprintf(outfile, "    Saving Disp Voxel Partition.\n");

    EAp = vars_["Disp_a"]->pointer();
    EBp = vars_["Disp_b"]->pointer();

    TA->copy(Locc_A_);
    for (int a = 0; a < na; a++) {
        C_DSCAL(nn,EAp[a][0],&TAp[0][a],na);
    }
    DA = Matrix::doublet(TA,Locc_A_,false,true);

    TB->copy(Locc_B_);
    for (int b = 0; b < nb; b++) {
        C_DSCAL(nn,EBp[b][0],&TBp[0][b],nb);
    }
    DB = Matrix::doublet(TB,Locc_B_,false,true);

    DA->add(DB); 
    DA->set_name("Disp");

    grid->zero();
    grid->compute_electronic(DA);
    grid->drop_raw("Disp.raw",int_scale);

    // ==> Dens <== //

    fprintf(outfile, "    Saving Dens Voxel Partition.\n");

    DA = Matrix::doublet(Locc_A_,Locc_A_,false,true);
    DB = Matrix::doublet(Locc_B_,Locc_B_,false,true);

    DA->add(DB); 
    DA->set_name("Dens");

    grid->zero();
    grid->compute_electronic(DA);
    grid->drop_raw("Dens.raw",dens_scale);
}

CubicDensityGrid::CubicDensityGrid(
        boost::shared_ptr<BasisSet> primary) :
        primary_(primary),
        mol_(primary->molecule()),
        options_(Process::environment.options)
{
    npoints_ = 0L;
    x_ = NULL;
    y_ = NULL;
    z_ = NULL;
    v_ = NULL;
}
CubicDensityGrid::~CubicDensityGrid()
{
    if (x_) delete x_;
    if (y_) delete y_;
    if (z_) delete z_;
    if (v_) delete v_;
}
void CubicDensityGrid::build_grid()
{
    if (options_["CUBIC_GRID_OVERAGE"].size() != 3) throw PSIEXCEPTION("CubicDensityGrid: CUBIC_GRID_OVERAGE not set.");
    double L[3];
    L[0] = options_["CUBIC_GRID_OVERAGE"][0].to_double();
    L[1] = options_["CUBIC_GRID_OVERAGE"][1].to_double();
    L[2] = options_["CUBIC_GRID_OVERAGE"][2].to_double();

    if (options_["CUBIC_GRID_SPACING"].size() != 3) throw PSIEXCEPTION("CubicDensityGrid: CUBIC_GRID_SPACING not set.");
    double D[3];
    D[0] = options_["CUBIC_GRID_SPACING"][0].to_double();
    D[1] = options_["CUBIC_GRID_SPACING"][1].to_double();
    D[2] = options_["CUBIC_GRID_SPACING"][2].to_double();

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

    npoints_ = (N[0] + 1L) * (N[1] + 1L) * (N[2] + 1L);
    x_ = new double[npoints_];
    y_ = new double[npoints_];
    z_ = new double[npoints_];
    v_ = new double[npoints_];

    double epsilon = options_.get_double("DFT_BASIS_TOLERANCE");
    extents_ = boost::shared_ptr<BasisExtents> (new BasisExtents(primary_, epsilon));

    int max_points = options_.get_int("DFT_BLOCK_MAX_POINTS");
    double xyz = pow((double) max_points, 1.0/3.0);
    nxyz_ = (size_t) pow((double) max_points, 1.0/3.0);

    size_t offset = 0L;
    for (int istart = 0L; istart <= N[0]; istart+=nxyz_) {
        int ni = (istart + nxyz_ > N[0] ? (N[0] + 1) - istart : nxyz_);
        for (int jstart = 0L; jstart <= N[1]; jstart+=nxyz_) {
            int nj = (jstart + nxyz_ > N[1] ? (N[1] + 1) - jstart : nxyz_);
            for (int kstart = 0L; kstart <= N[2]; kstart+=nxyz_) {
                int nk = (kstart + nxyz_ > N[2] ? (N[2] + 1) - kstart : nxyz_);

                double* xp = &x_[offset];
                double* yp = &y_[offset];
                double* zp = &z_[offset];
                double* vp = &v_[offset];
                
                size_t block_size = 0L;
                for (int i = istart; i < istart + ni; i++) {
                    for (int j = jstart; j < jstart + nj; j++) {
                        for (int k = kstart; k < kstart + nk; k++) {
                            x_[offset] = O[0] + i * D[0];
                            y_[offset] = O[1] + j * D[1];
                            z_[offset] = O[2] + k * D[2];
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
void CubicDensityGrid::print_header()
{
    fprintf(outfile,"  ==> CubicDensityGrid <==\n\n");

    fprintf(outfile,"    Total Points = %16zu\n", npoints_);
    fprintf(outfile,"    XYZ Blocking = %16zu\n", nxyz_);
    fprintf(outfile,"    X Points     = %16zu\n", N_[0] + 1L);
    fprintf(outfile,"    Y Points     = %16zu\n", N_[1] + 1L);
    fprintf(outfile,"    Z Points     = %16zu\n", N_[2] + 1L);
    fprintf(outfile,"    X Spacing    = %16.3E\n", D_[0]);
    fprintf(outfile,"    Y Spacing    = %16.3E\n", D_[1]);
    fprintf(outfile,"    Z Spacing    = %16.3E\n", D_[2]);
    fprintf(outfile,"    X Maximum    = %16.3E\n", O_[0]);
    fprintf(outfile,"    Y Maximum    = %16.3E\n", O_[1]);
    fprintf(outfile,"    Z Maximum    = %16.3E\n", O_[2]);
    fprintf(outfile,"    X Minimum    = %16.3E\n", O_[0] + D_[0] * N_[0]);
    fprintf(outfile,"    Y Minimum    = %16.3E\n", O_[1] + D_[1] * N_[1]);
    fprintf(outfile,"    Z Minimum    = %16.3E\n", O_[2] + D_[2] * N_[2]);
    fprintf(outfile,"\n");

    fflush(outfile);
}
void CubicDensityGrid::zero()
{
    ::memset(v_, '\0', sizeof(double) * npoints_);
}
void CubicDensityGrid::compute_electronic(boost::shared_ptr<Matrix> D)
{
    if (!npoints_) throw PSIEXCEPTION("CubicDensityGrid::compute: call build_grid first");
   
    points_->set_pointers(D);
    boost::shared_ptr<Vector> rho = points_->point_value("RHO_A");
    double* rhop = rho->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_points(blocks_[ind]);
        size_t npoints = blocks_[ind]->npoints();
        C_DAXPY(npoints,1.0,rhop,1,&v_[offset],1);        
        offset += npoints;
    }
}
void CubicDensityGrid::compute_atomic(boost::shared_ptr<Matrix> Q)
{
    if (Q->ncol() != 5) throw PSIEXCEPTION("CubicDensityGrid::compute_atomic: Q should be (x,y,z,V,\\alpha) in a.u.");
    int natom = Q->nrow();
    double** Qp = Q->pointer();

    double* a = new double[natom];
    double* N = new double[natom];
    double* x = new double[natom];
    double* y = new double[natom];
    double* z = new double[natom];

    for (int A = 0; A < natom; A++) {
        x[A] = Qp[A][0];
        y[A] = Qp[A][1];
        z[A] = Qp[A][2];
        a[A] = Qp[A][4];
        N[A] = 0.5 * pow(M_PI / a[A], 3.0/2.0) * Qp[A][3]; 
    }

    #pragma omp parallel for
    for (size_t ind = 0; ind < npoints_; ind++) {
        double xp = x_[ind];
        double yp = y_[ind];
        double zp = z_[ind];
        double val = 0.0; 
        for (int A = 0; A < natom; A++) {
            double R2 = (xp - x[A]) * (xp - x[A]) + 
                        (yp - y[A]) * (yp - y[A]) + 
                        (zp - z[A]) * (zp - z[A]);
            val += N[A] * exp(-a[A] * R2);
        }
        v_[ind] += val;
    }

    delete[] a;
    delete[] N;
    delete[] x;
    delete[] y;
    delete[] z;
}
void CubicDensityGrid::drop_raw(const std::string& file, double clamp)
{
    if (!npoints_) throw PSIEXCEPTION("CubicDensityGrid::drop_raw: call build_grid first");

    //ASCII Variant
    //std::stringstream ss2;
    //ss2 << V_->name() << ".debug";
    //FILE* fh2 = fopen(ss2.str().c_str(), "w");
    //for (size_t ind = 0; ind < npoints_; ind++) {
    //    fprintf(fh2,"  %16zu %24.16E %24.16E %24.16E %24.16E\n", ind, x_[ind], y_[ind], z_[ind], v_[ind]);
    //}
    //fclose(fh2);

    double s = 1.0 / clamp;
    double maxval = 0.0;

    float* v2 = new float[npoints_];
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
                            double val = v_[offset];
                            maxval = (maxval >= fabs(val) ? maxval : fabs(val));
                            val = (val <= clamp ? val : clamp);
                            val = (val >= -clamp ? val : -clamp);
                            val *= s;
                            v2[index] = val;
                            offset++;
                        }
                    }
                }
            }
        }
    }

    fprintf(outfile,"    Max val = %11.3E out of %11.3E: %s.\n", maxval, clamp, (clamp >= maxval ? "No clamping" : "Clamped")); 

    //Dirty Hack: I love it!
    v2[npoints_-1L] =  0.0;
    v2[npoints_-2L] =  0.0;
    v2[npoints_-3L] =  1.0;
    v2[npoints_-4L] = -1.0;
    
    FILE* fh = fopen(file.c_str(), "wb");
    fwrite(v2,sizeof(float),npoints_,fh);
    fclose(fh);
}
void CubicDensityGrid::drop_uvf(const std::string& file, double clamp)
{
    if (!npoints_) throw PSIEXCEPTION("CubicDensityGrid::drop_raw: call build_grid first");

    double s = 1.0 / clamp;
    double maxval = 0.0;

    float* v2 = new float[npoints_];
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
                            double val = v_[offset];
                            maxval = (maxval >= fabs(val) ? maxval : fabs(val));
                            val = (val <= clamp ? val : clamp);
                            val = (val >= -clamp ? val : -clamp);
                            val *= s;
                            v2[index] = val;
                            offset++;
                        }
                    }
                }
            }
        }
    }

    fprintf(outfile,"    Max val = %11.3E out of %11.3E: %s.\n", maxval, clamp, (clamp >= maxval ? "No clamping" : "Clamped")); 

    //Dirty Hack: I love it!
    v2[npoints_-1L] =  0.0;
    v2[npoints_-2L] =  0.0;
    v2[npoints_-3L] =  1.0;
    v2[npoints_-4L] = -1.0;
    
    FILE* fh = fopen(file.c_str(), "wb");
    const char* header = "UVF-DATA";
    fwrite(header,sizeof(char),8,fh);
    const char* endian = "\0";
    fwrite(endian,sizeof(char),1,fh);
    unsigned long int version = 3L;
    fwrite(&version,sizeof(unsigned long int),1,fh);
    unsigned long int nchecksum = 0L;   
    fwrite(&nchecksum,sizeof(unsigned long int),1,fh);
    unsigned long int offset = 0L;
    fwrite(&offset,sizeof(unsigned long int),1,fh);
    
    throw PSIEXCEPTION("Not implemented");
 
    fclose(fh);
}

}}
