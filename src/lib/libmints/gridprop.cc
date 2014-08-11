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
#include <libmints/sieve.h>
#include <libmints/mints.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
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
    w_ = NULL;
    N_ = new int[3];
    D_ = new double[3];
    O_ = new double[3];

    build_grid(); // Defaults from Options
}
CubicScalarGrid::~CubicScalarGrid()
{
    if (x_) delete[] x_;
    if (y_) delete[] y_;
    if (z_) delete[] z_;
    if (w_) delete[] w_;
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
void CubicScalarGrid::build_grid(boost::shared_ptr<CubicScalarGrid> other)
{
    filepath_ = other->filepath();

    for (int k = 0; k < 3; k++) {
        N_[k] = other->N()[k];
        O_[k] = other->O()[k];
        D_[k] = other->D()[k];
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
    if (x_) delete[] x_;
    if (y_) delete[] y_;
    if (z_) delete[] z_;
    if (w_) delete[] w_;

    npoints_ = (N_[0] + 1L) * (N_[1] + 1L) * (N_[2] + 1L);
    x_ = new double[npoints_];
    y_ = new double[npoints_];
    z_ = new double[npoints_];
    w_ = new double[npoints_];

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
                double* wp = &w_[offset];
                
                size_t block_size = 0L;
                for (int i = istart; i < istart + ni; i++) {
                    for (int j = jstart; j < jstart + nj; j++) {
                        for (int k = kstart; k < kstart + nk; k++) {
                            x_[offset] = O_[0] + i * D_[0];
                            y_[offset] = O_[1] + j * D_[1];
                            z_[offset] = O_[2] + k * D_[2];
                            w_[offset] = D_[0] * D_[1] * D_[2];
                            offset++;
                            block_size++;
                        }
                    }
                }
                blocks_.push_back(boost::shared_ptr<BlockOPoints>(new BlockOPoints(block_size,xp,yp,zp,wp,extents_)));
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
    outfile->Printf("  ==> CubicScalarGrid <==\n\n");

    outfile->Printf("    Filepath     = %s\n", filepath_.c_str());
    outfile->Printf("    Total Points = %16zu\n", npoints_);
    outfile->Printf("    XYZ Blocking = %16zu\n", nxyz_);
    outfile->Printf("    X Points     = %16zu\n", N_[0] + 1L);
    outfile->Printf("    Y Points     = %16zu\n", N_[1] + 1L);
    outfile->Printf("    Z Points     = %16zu\n", N_[2] + 1L);
    outfile->Printf("    X Spacing    = %16.3E\n", D_[0]);
    outfile->Printf("    Y Spacing    = %16.3E\n", D_[1]);
    outfile->Printf("    Z Spacing    = %16.3E\n", D_[2]);
    outfile->Printf("    X Minimum    = %16.3E\n", O_[0]);
    outfile->Printf("    Y Minimum    = %16.3E\n", O_[1]);
    outfile->Printf("    Z Minimum    = %16.3E\n", O_[2]);
    outfile->Printf("    X Maximum    = %16.3E\n", O_[0] + D_[0] * N_[0]);
    outfile->Printf("    Y Maximum    = %16.3E\n", O_[1] + D_[1] * N_[1]);
    outfile->Printf("    Z Maximum    = %16.3E\n", O_[2] + D_[2] * N_[2]);
    outfile->Printf("\n");

    
}
void CubicScalarGrid::write_gen_file(double* v, const std::string& name, const std::string& type)
{
    if (type == "CUBE") {
        write_cube_file(v, name);
    } else {
        throw PSIEXCEPTION("CubicScalarGrid: Unrecognized output file type");
    }
}
void CubicScalarGrid::write_cube_file(double* v, const std::string& name)
{
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
    outfile->Printf(fh, "PSI4 Gaussian Cube File.\n");
    outfile->Printf(fh, "Property: %s\n", name.c_str());
   
    // Number of atoms plus origin of data 
    outfile->Printf(fh, "%6d %10.6f %10.6f %10.6f\n", mol_->natom(), O_[0], O_[1], O_[2]);
    
    // Number of points along axis, displacement along x,y,z
    outfile->Printf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[0] + 1, D_[0], 0.0, 0.0);
    outfile->Printf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[1] + 1, 0.0, D_[1], 0.0);
    outfile->Printf(fh, "%6d %10.6f %10.6f %10.6f\n", N_[2] + 1, 0.0, 0.0, D_[2]);

    // Atoms of molecule (Z, Q?, x, y, z)
    for (int A = 0; A < mol_->natom(); A++) {
        outfile->Printf(fh, "%3d %10.6f %10.6f %10.6f %10.6f\n", (int) mol_->Z(A), 0.0, mol_->x(A), mol_->y(A), mol_->z(A));
    }

    // Data, striped (x, y, z)
    for (size_t ind = 0; ind < npoints_; ind++) {
        outfile->Printf(fh, "%12.5E ", v2[ind]);
        if (ind % 6 == 5) outfile->Printf(fh,"\n");
    }

    fclose(fh);   
}
void CubicScalarGrid::add_density(double* v, boost::shared_ptr<Matrix> D)
{
    points_->set_pointers(D);
    boost::shared_ptr<Vector> rho = points_->point_value("RHO_A");
    double* rhop = rho->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_points(blocks_[ind]);
        size_t npoints = blocks_[ind]->npoints();
        C_DAXPY(npoints,1.0,rhop,1,&v[offset],1);        
        offset += npoints;
    }
}
void CubicScalarGrid::add_esp(double* v, boost::shared_ptr<Matrix> D)
{
    // => Auxiliary Basis Set (TODO: Get appropriate default) <= //

    boost::shared_ptr<BasisSetParser> parser(new Gaussian94BasisSetParser());
    boost::shared_ptr<BasisSet> auxiliary = BasisSet::construct(parser, primary_->molecule(), "DF_BASIS_SCF");
    
    // => DF Options (TODO: Should these be in here?) <= //

    double cutoff    = options_.get_double("INTS_TOLERANCE");
    double condition = options_.get_double("DF_FITTING_CONDITION");

    // => Sizing <= //

    int nbf  = primary_->nbf();
    int naux = auxiliary->nbf();
    int maxP = auxiliary->max_function_per_shell();

    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    // => Density Fitting (TODO: Could be sped up) <= //

    boost::shared_ptr<IntegralFactory> Ifact(new IntegralFactory(auxiliary,BasisSet::zero_ao_basis_set(), primary_, primary_));
    std::vector<boost::shared_ptr<TwoBodyAOInt> > ints;
    for (int thread = 0; thread < nthreads; thread++) {
        ints.push_back(boost::shared_ptr<TwoBodyAOInt>(Ifact->eri()));
    }

    boost::shared_ptr<ERISieve> sieve(new ERISieve(primary_, cutoff));
    const std::vector<std::pair<int,int> >& pairs = sieve->shell_pairs();  

    boost::shared_ptr<Vector> c(new Vector("c", naux));
    double* cp = c->pointer();

    boost::shared_ptr<Matrix> Amn(new Matrix("Amn", maxP, nbf*nbf));
    double** Amnp = Amn->pointer();

    double** Dp = D->pointer();

    for (int P = 0; P < auxiliary->nshell(); P++) {

        int nP = auxiliary->shell(P).nfunction();
        int oP = auxiliary->shell(P).function_index();

        Amn->zero();        

        // Integrals
        #pragma omp parallel for schedule(dynamic)
        for (int task = 0; task < pairs.size(); task++) {
            
            int thread = 0;
            #ifdef _OPENMP
                thread = omp_get_thread_num();
            #endif

            int M = pairs[task].first;
            int N = pairs[task].second;

            ints[thread]->compute_shell(P,0,M,N);
            const double* buffer = ints[thread]->buffer();

            int nM = primary_->shell(M).nfunction();
            int oM = primary_->shell(M).function_index();
            int nN = primary_->shell(N).nfunction();
            int oN = primary_->shell(N).function_index();

            int index = 0;
            for (int p = 0; p < nP; p++) {
                for (int m = 0; m < nM; m++) {
                    for (int n = 0; n < nN; n++) {
                        Amnp[p][(m + oM) * nbf + (n + oN)] =
                        Amnp[p][(n + oN) * nbf + (m + oM)] =
                        buffer[index++];
                    }
                }
            }
            
        }
        
        // Contraction
        C_DGEMV('N', nP, nbf*nbf, 1.0, Amnp[0], nbf*nbf, Dp[0], 1, 0.0, cp + oP, 1);

    }

    Amn.reset();
    Ifact.reset();
    ints.clear();

    boost::shared_ptr<Matrix> J(new Matrix("J", naux, naux));
    double** Jp = J->pointer();

    boost::shared_ptr<IntegralFactory> Jfact(new IntegralFactory(
        auxiliary,BasisSet::zero_ao_basis_set(), 
        auxiliary,BasisSet::zero_ao_basis_set()));

    boost::shared_ptr<TwoBodyAOInt> Jints(Jfact->eri());
    const double* Jbuffer = Jints->buffer();

    for (int P = 0; P < auxiliary->nshell(); P++) {
        int nP = auxiliary->shell(P).nfunction();
        int oP = auxiliary->shell(P).function_index();

        for (int Q = 0; Q <= P; Q++) {
            int nQ = auxiliary->shell(Q).nfunction();
            int oQ = auxiliary->shell(Q).function_index();

            Jints->compute_shell(P,0,Q,0);

            int index = 0;
            for (int p = 0; p < nP; p++) {
                for (int q = 0; q < nQ; q++) {
                    Jp[p + oP][q + oQ] = 
                    Jp[q + oQ][p + oP] = 
                    Jbuffer[index++];
                }
            }
        }
    }

    Jfact.reset();
    Jints.reset(); 
   
    J->power(-1.0, condition); 

    boost::shared_ptr<Vector> d(new Vector("d", naux));
    double* dp = d->pointer();

    C_DGEMV('N',naux,naux,1.0,Jp[0],naux,cp,1,0.0,dp,1);

    //c->print();
    //d->print();

    J.reset();

    // => Electronic Part <= //

    boost::shared_ptr<IntegralFactory> Vfact(new IntegralFactory(auxiliary,BasisSet::zero_ao_basis_set()));
    std::vector<boost::shared_ptr<Matrix> > ZxyzT;
    std::vector<boost::shared_ptr<Matrix> > VtempT;
    std::vector<boost::shared_ptr<PotentialInt> > VintT;
    for (int thread = 0; thread < nthreads; thread++) {
        ZxyzT.push_back(boost::shared_ptr<Matrix>(new Matrix("Zxyz",1,4)));
        VtempT.push_back(boost::shared_ptr<Matrix>(new Matrix("Vtemp",naux,1)));
        VintT.push_back(boost::shared_ptr<PotentialInt>(static_cast<PotentialInt*>(Vfact->ao_potential())));
        VintT[thread]->set_charge_field(ZxyzT[thread]);
    }

    #pragma omp parallel for schedule(dynamic)
    for (int P = 0; P < npoints_; P++) {

        // Thread info
        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif

        // Pointers
        double** ZxyzTp = ZxyzT[thread]->pointer();
        double** VtempTp = VtempT[thread]->pointer();

        // Integrals
        VtempT[thread]->zero();
        ZxyzTp[0][0] = 1.0;
        ZxyzTp[0][1] = x_[P];
        ZxyzTp[0][2] = y_[P];
        ZxyzTp[0][3] = z_[P]; 
        VintT[thread]->compute(VtempT[thread]);

        // Contraction 
        v[P] += C_DDOT(naux,dp,1,VtempTp[0],1);
    }
    
    // => Nuclear Part <= //

    for (int A = 0; A < mol_->natom(); A++) {
        double Z = mol_->Z(A);
        double x = mol_->x(A);
        double y = mol_->y(A);
        double z = mol_->z(A);
        for (int P = 0; P < npoints_; P++) {
            double R = sqrt(
                (x - x_[P]) * (x - x_[P]) +
                (y - y_[P]) * (y - y_[P]) +
                (z - z_[P]) * (z - z_[P]));
            v[P] += (R >= 1.0E-15 ? Z / R : 0.0);
        }
    }
}
void CubicScalarGrid::add_basis_functions(double** v, const std::vector<int>& indices)
{
    boost::shared_ptr<Matrix> phi = points_->basis_value("PHI");
    double** phip = phi->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_functions(blocks_[ind]);

        size_t npoints = blocks_[ind]->npoints();
        const std::vector<int>& function_map = blocks_[ind]->functions_local_to_global();
        int nlocal  = function_map.size();
        int nglobal = points_->max_functions();

        for (int ind1 = 0; ind1 < indices.size(); ind1++) {
            for (int ind2 = 0; ind2 < function_map.size(); ind2++) {
                if (indices[ind1] == function_map[ind2]) {
                    C_DAXPY(npoints,1.0,&phip[0][ind2],nglobal,&v[ind1][offset],1);
                }                
            }
        }

        offset += npoints;
    }
}
void CubicScalarGrid::add_orbitals(double** v, boost::shared_ptr<Matrix> C)
{
    int na = C->colspi()[0];    

    points_->set_Cs(C);
    boost::shared_ptr<Matrix> psi = points_->orbital_value("PSI_A");
    double** psip = psi->pointer();

    size_t offset = 0L;
    for (int ind = 0; ind < blocks_.size(); ind++) {
        points_->compute_orbitals(blocks_[ind]);

        size_t npoints = blocks_[ind]->npoints();
        for (int a = 0; a < na; a++) {
            C_DAXPY(npoints,1.0,psip[a],1,&v[a][offset],1);
        }    

        offset += npoints;
    }
}
void CubicScalarGrid::compute_density(boost::shared_ptr<Matrix> D, const std::string& name, const std::string& type)
{
    double* v = new double[npoints_];
    memset(v,'\0',npoints_*sizeof(double));
    add_density(v,D);
    write_gen_file(v,name,type);
    delete[] v;
}
void CubicScalarGrid::compute_esp(boost::shared_ptr<Matrix> D, const std::string& name, const std::string& type)
{
    double* v = new double[npoints_];
    memset(v,'\0',npoints_*sizeof(double));
    add_esp(v,D);
    write_gen_file(v,name,type);
    delete[] v;
}
void CubicScalarGrid::compute_basis_functions(const std::vector<int>& indices, const std::string& name, const std::string& type)
{
    double** v = block_matrix(indices.size(), npoints_);
    memset(v[0],'\0',indices.size()*npoints_*sizeof(double));
    add_basis_functions(v, indices);
    for (int k = 0; k < indices.size(); k++) {
        std::stringstream ss; 
        ss << name << "_" << (indices[k] + 1);
        write_gen_file(v[k],ss.str(),type);
    }   
    free_block(v);
}
void CubicScalarGrid::compute_orbitals(boost::shared_ptr<Matrix> C, const std::vector<int>& indices, const std::string& name, const std::string& type)
{
    boost::shared_ptr<Matrix> C2(new Matrix(primary_->nbf(), indices.size()));  
    double** Cp  = C->pointer();
    double** C2p = C2->pointer();
    for (int k = 0; k < indices.size(); k++) {
        C_DCOPY(primary_->nbf(), &Cp[0][indices[k]], C->colspi()[0], &C2p[0][k], C2->colspi()[0]); 
    }
    double** v = block_matrix(indices.size(), npoints_);
    memset(v[0],'\0',indices.size()*npoints_*sizeof(double));
    add_orbitals(v, C2);
    for (int k = 0; k < indices.size(); k++) {
        std::stringstream ss; 
        ss << name << "_" << (indices[k] + 1);
        write_gen_file(v[k],ss.str(),type);
    }   
    free_block(v);
}

}
