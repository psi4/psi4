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

#include "atomic.h"
#include <libmints/mints.h>
#include <libfock/cubature.h>
#include <libfock/points.h>
#include <libdiis/diismanager.h>
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

AtomicDensity::AtomicDensity() : 
    print_(1),
    debug_(0),
    bench_(0)
{
}
AtomicDensity::~AtomicDensity() 
{
}
boost::shared_ptr<AtomicDensity> AtomicDensity::build(const std::string& type, boost::shared_ptr<BasisSet> basis, Options& options)
{
    AtomicDensity* target;
    if (type == "STOCKHOLDER") {
        StockholderDensity* isa = new StockholderDensity();
        isa->convergence_ = options.get_double("ISA_CONVERGENCE");
        isa->maxiter_ = options.get_int("ISA_MAXITER");
        isa->diis_ = options.get_bool("ISA_DIIS");
        isa->diis_min_vecs_ = options.get_int("ISA_DIIS_MIN_VECS");
        isa->diis_max_vecs_ = options.get_int("ISA_DIIS_MAX_VECS");
        isa->diis_flush_ = options.get_int("ISA_DIIS_FLUSH");
        target = static_cast<AtomicDensity*>(isa); 
    } else {
        throw PSIEXCEPTION("AtomicDensity::build: Unrecognized Atomic Density Type");
    }

    target->primary_ = basis;
    target->molecule_ = basis->molecule();
    target->grid_ = boost::shared_ptr<MolecularGrid>(new DFTGrid(basis->molecule(),basis,options));
    return boost::shared_ptr<AtomicDensity>(target);
}
void AtomicDensity::compute_total_density()
{
    // Fast grid sizing
    int npoints2 = grid_->npoints();
    double* x2p = grid_->x();
    double* y2p = grid_->y();
    double* z2p = grid_->z();
    double* w2p = grid_->w();
    int* index = grid_->index(); 

    // Density computer (in blocks)
    boost::shared_ptr<RKSFunctions> points = boost::shared_ptr<RKSFunctions>(new RKSFunctions(primary_,grid_->max_points(),grid_->max_functions()));
    points->set_ansatz(0);
    points->set_pointers(D_);
    boost::shared_ptr<Vector> rho3 = points->point_value("RHO_A");
    double* rho3p = rho3->pointer();

    // Build density (in fast ordering)
    boost::shared_ptr<Vector> rho2(new Vector("rho2", npoints2));
    double* rho2p = rho2->pointer();
    const std::vector<boost::shared_ptr<BlockOPoints> >& blocks = grid_->blocks();    
    size_t offset = 0L;
    for (int ind = 0; ind < blocks.size(); ind++) {
        int npoints = blocks[ind]->npoints();
        points->compute_points(blocks[ind]);
        C_DAXPY(npoints,1.0,rho3p,1,&rho2p[offset],1);        
        offset += npoints;
    }

    // How many clean points?
    const std::vector<std::vector<boost::shared_ptr<SphericalGrid> > >& spheres = grid_->spherical_grids();
    int npoints = 0;
    for (int A = 0; A < spheres.size(); A++) {
        for (int R = 0; R < spheres[A].size(); R++) {
            npoints += spheres[A][R]->npoints();
        }
    }

    // Targets
    x_ = boost::shared_ptr<Vector>(new Vector("x",npoints)); 
    y_ = boost::shared_ptr<Vector>(new Vector("y",npoints)); 
    z_ = boost::shared_ptr<Vector>(new Vector("z",npoints)); 
    w_ = boost::shared_ptr<Vector>(new Vector("w",npoints)); 
    rho_ = boost::shared_ptr<Vector>(new Vector("rho",npoints)); 
    double* xp = x_->pointer();
    double* yp = y_->pointer();
    double* zp = z_->pointer();
    double* wp = w_->pointer();
    double* rhop = rho_->pointer();

    for (int i = 0; i < npoints2; i++) {
        xp[index[i]] = x2p[i];
        yp[index[i]] = y2p[i];
        zp[index[i]] = z2p[i];
        wp[index[i]] = w2p[i];
        rhop[index[i]] = rho2p[i];
    }

    // Debug stuff 
    //x_->print();
    //y_->print();
    //z_->print();
    //w_->print();
    //rho_->print();
}

StockholderDensity::StockholderDensity() : AtomicDensity() 
{
}
StockholderDensity::~StockholderDensity()
{
}
void StockholderDensity::print_header() const 
{
    fprintf(outfile,"  ==> Stockholder Atomic Densities <==\n\n");
    molecule_->print();
    primary_->print();
    grid_->print();

    fflush(outfile);
}
void StockholderDensity::compute(boost::shared_ptr<Matrix> D)
{
    D_ = D;

    print_header();

    compute_total_density();

    // ==> Main ISA Algorithm <== //

    // Where are the true atoms?
    std::vector<int> Aind;
    std::vector<int> Aind2;
    for (int A = 0; A < molecule_->natom(); A++) {
        if (molecule_->Z(A) != 0) {
            Aind2.push_back(Aind.size());
            Aind.push_back(A);
        } else {
            Aind2.push_back(-1);
        }
    }
    int nA = Aind.size();
    int nA2 = molecule_->natom();

    // Grid indexing
    int nP = x_->dimpi()[0];
    double* xp = x_->pointer();
    double* yp = y_->pointer();
    double* zp = z_->pointer();
    double* wp = w_->pointer();
    double* rhop = rho_->pointer();

    //fprintf(outfile,"  Electron count is %24.16E\n", C_DDOT(nP,wp,1,rhop,1));

    const std::vector<boost::shared_ptr<RadialGrid> >& rads = grid_->radial_grids();
    const std::vector<std::vector<boost::shared_ptr<SphericalGrid> > >& spheres = grid_->spherical_grids();

    // True atom spherical quadratures 
    std::vector<std::vector<int> > orders;
    rs_.resize(nA);
    ws_.resize(nA);
    orders.resize(nA);
    int nstate = 0;
    for (int A = 0; A < nA; A++) {
        int Aabs = Aind[A];
        std::vector<double> rs2;
        std::vector<double> ws2;
        std::vector<std::pair<double, int> > index;
        for (int R = 0; R < rads[Aabs]->npoints(); R++) {
            rs2.push_back(rads[Aabs]->r()[R]);
            ws2.push_back(1.0); // Initial guess
            index.push_back(std::pair<double, int>(rs2[R],R));
            nstate++;
        } 
        std::sort(index.begin(), index.end());
        for (int R = 0; R < rads[Aabs]->npoints(); R++) {
            int Rold = index[R].second;
            rs_[A].push_back(rs2[Rold]);
            ws_[A].push_back(ws2[Rold]);
            orders[A].push_back(Rold);
        } 
    }

    std::vector<std::vector<int> > orders2;
    orders2.resize(nA);
    for (int A = 0; A < nA; A++) {
        orders2[A].resize(orders[A].size());
        for (int R = 0; R < orders[A].size(); R++) {
            orders2[A][orders[A][R]] = R;
        }
    }    

    // Low-memory atomic charges target
    int max_points = 0;
    std::vector<std::vector<std::vector<double> > > Q;
    std::vector<int> atomic_points;
    Q.resize(nA);
    atomic_points.resize(nA);
    for (int A = 0; A < nA; A++) {
        int Aabs = Aind[A];
        Q[A].resize(orders[A].size());
        int atom_points = 0;
        for (int R = 0; R < orders[A].size(); R++) {
            int Rabs = orders[A][R];
            Q[A][R].resize(spheres[Aabs][Rabs]->npoints());
            atom_points += spheres[Aabs][Rabs]->npoints();
        }
        atomic_points[A] = atom_points;
        max_points = (max_points >= atom_points ? max_points : atom_points);
    }
    
    // Temps 
    boost::shared_ptr<Matrix> Q2(new Matrix("Q2", 1, max_points));
    double** Q2p = Q2->pointer();

    // I like to work in log space for interpolation window root finding
    std::vector<std::vector<double> > ls;
    ls.resize(nA);
    for (int A = 0; A < rs_.size(); A++) {
        for (int R = 0; R < rs_[A].size(); R++) {
            ls[A].push_back(log(rs_[A][R]));
        }
    }

    // DIIS Setup
    boost::shared_ptr<Matrix> state(new Matrix("State", nstate, 1));
    boost::shared_ptr<Matrix> error(new Matrix("Error", nstate, 1));
    double** statep = state->pointer();
    double** errorp = error->pointer();

    boost::shared_ptr<DIISManager> diis_manager(new DIISManager(diis_max_vecs_, "ISA DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
    diis_manager->set_error_vector_size(1, DIISEntry::Matrix, error.get());
    diis_manager->set_vector_size(1, DIISEntry::Matrix, state.get());

    // Store last iteration
    std::vector<std::vector<double> > ws_old = ws_;

    // => Master Loop <= //

    fprintf(outfile, "   > ISA Iterations <\n\n");
    fprintf(outfile, "    Convergence:        %11.3E\n", convergence_);
    fprintf(outfile, "    Maximum iterations: %11d\n", maxiter_);   
    fprintf(outfile, "    DIIS:               %11s\n", (diis_ ? "Yes" : "No"));
    fprintf(outfile, "    DIIS Min Vecs:      %11d\n", diis_min_vecs_);
    fprintf(outfile, "    DIIS Max Vecs:      %11d\n", diis_max_vecs_);
    fprintf(outfile, "    DIIS Flush Vecs:    %11d\n", diis_flush_);
    fprintf(outfile, "\n");
    fflush(outfile);

    bool converged = false;
    for (int iter = 0, diis_iter = 0; iter <= maxiter_; iter++) {

        // Checkpoint last iteration
        ws_old = ws_;
    
        // Needed atomic weights
        int offset = 0;
        for (int A = 0; A < nA2; A++) {
            if (molecule_->Z(A) != 0) {
                int Arel = Aind2[A]; 
                compute_weights(atomic_points[Arel],&xp[offset],&yp[offset],&zp[offset],Q2p,&rhop[offset],Arel);
                int offset2 = 0;
                for (int R = 0; R < orders[Arel].size(); R++) {
                    for (int k = 0; k < spheres[A][orders[Arel][R]]->npoints(); k++) {
                        Q[Arel][orders[Arel][R]][k] = Q2p[0][offset2];
                        offset2++; 
                    }
                }
            }
            for (int R = 0; R < spheres[A].size(); R++) {
                offset += spheres[A][R]->npoints();
            }
        }

        // Spherical averaging
        for (int A = 0; A < Q.size(); A++) {
            for (int R = 0; R < Q[A].size(); R++) { 
                double val = 0.0;
                for (int k = 0; k < Q[A][R].size(); k++) {
                    val += Q[A][R][k];
                }
                val /= Q[A][R].size();
                ws_[A][R] = val;
            }
        }

        // Residual and error norm
        double norm = 0.0;
        std::vector<std::vector<double> > dws = ws_;
        for (int A = 0; A < ws_.size(); A++) {
            for (int R = 0; R < ws_[A].size(); R++) {
                dws[A][R] = (ws_[A][R] - ws_old[A][R]) * rs_[A][R] * rs_[A][R];
                norm += abs(dws[A][R]);
            }
        }

        // Print iterative trace
        fprintf(outfile,"    @ISA Iter %4d %24.16E ", iter, norm);
        fflush(outfile);

        // Convergence check
        if (norm < convergence_) { 
            fprintf(outfile,"\n");
            converged = true;
            break; 
        }

        // Periodically flush the DIIS subspace
        if (diis_ && iter > 0 && iter % diis_flush_ == 0) {
            diis_manager->reset_subspace();
            diis_iter = 0;
        }

        // DIIS Add
        if (diis_ && iter > 0) {
            int offset2 = 0;
            for (int A = 0; A < ws_.size(); A++) {
                for (int R = 0; R < ws_[A].size(); R++) {
                    statep[0][offset2] = (ws_[A][R] == 0.0 ? -std::numeric_limits<double>::infinity() : log(ws_[A][R]));
                    //statep[0][offset2] = ws_[A][R];
                    errorp[0][offset2] = dws[A][R];
                    offset2++;
                }
            }
            diis_manager->add_entry(2,error.get(),state.get());
            diis_iter++;
        }

        // DIIS Extrapolate
        if (diis_ && diis_iter >= diis_min_vecs_) {
            diis_manager->extrapolate(1,state.get());
            int offset2 = 0;
            for (int A = 0; A < ws_.size(); A++) {
                for (int R = 0; R < ws_[A].size(); R++) {
                    ws_[A][R] = exp(statep[0][offset2]);
                    //ws_[A][R] = statep[0][offset2];
                    offset2++;
                }
            }
            fprintf(outfile,"DIIS");
        }

        fprintf(outfile,"\n");
        fflush(outfile);
 
    }

    diis_manager->delete_diis_file();

    fprintf(outfile,"\n");
    if (converged) { 
        fprintf(outfile,"    ISA Converged.\n\n"); 
    } else {
        fprintf(outfile,"    ISA Failed.\n\n"); 
    }
    fflush(outfile);

    // => Compute normalizations for later <= //
    
    // Target
    N_ = boost::shared_ptr<Vector>(new Vector("N", nA));
    double* Np = N_->pointer();
    
    boost::shared_ptr<Matrix> Q3(new Matrix("Q3", nA, max_points));
    double** Q3p = Q3->pointer();

    for (int index = 0; index < nP; index+=max_points) {
        int points = (index + max_points >= nP ? nP - index : max_points);
        compute_weights(points,&xp[index],&yp[index],&zp[index],Q3p,&rhop[index]);
        for (int A = 0; A < nA; A++) {
            Np[A] += C_DDOT(points,&wp[index],1,Q3p[A],1);
        }
    }
}
void StockholderDensity::compute_weights(int nP, double* xp, double* yp, double* zp, double** wp, double* rhop, int atom)
{
    // Where are the true atoms?
    std::vector<int> Aind;
    std::vector<int> Aind2;
    for (int A = 0; A < molecule_->natom(); A++) {
        if (molecule_->Z(A) != 0) {
            Aind2.push_back(Aind.size());
            Aind.push_back(A);
        } else {
            Aind2.push_back(-1);
        }
    }
    int nA = Aind.size();
    int nA2 = molecule_->natom();
    
    // I like to work in log space for interpolation window root finding
    std::vector<std::vector<double> > ls;
    ls.resize(nA);
    for (int A = 0; A < rs_.size(); A++) {
        for (int R = 0; R < rs_[A].size(); R++) {
            ls[A].push_back(log(rs_[A][R]));
        }
    }

    // Doesn't like being on the stack?
    int nthreads = 1;
    #ifdef _OPENMP
        nthreads = omp_get_max_threads();
    #endif

    std::vector<double* > wA2;
    for (int thread = 0; thread < nthreads; thread++) {
        wA2.push_back(new double[nA]);
    }
    
    // Compute Q_A^P via W_A^P
    #pragma omp parallel for schedule(dynamic) 
    for (int P = 0; P < nP; P++) {
        double xc = xp[P];
        double yc = yp[P];
        double zc = zp[P];

        int thread = 0;
        #ifdef _OPENMP
            thread = omp_get_thread_num();
        #endif
        double* wA = wA2[thread];
        double wT = 0.0;
    
        for (int A = 0; A < nA; A++) {
            int Aabs = Aind[A];

            // Default value
            wA[A] = 0.0;

            // What is R_A^P?
            double xA = molecule_->x(Aabs);
            double yA = molecule_->y(Aabs);
            double zA = molecule_->z(Aabs);
            double R = sqrt((xA - xc) * (xA - xc) +
                            (yA - yc) * (yA - yc) +
                            (zA - zc) * (zA - zc));

            // Interpolation problem
            const std::vector<double>& rv = rs_[A];
            const std::vector<double>& lv = ls[A];
            const std::vector<double>& wv = ws_[A];
            int nR = rv.size();

            // Root solving to determine window (log-transformed Regula Falsi/Illinois, like a ninja)
            int ind = 0;
            double lR = log(R);
            if (R <= rv[0]) {
                if (fabs(R - rv[0]) < 1.0E-12 * rv[0]) {
                    ind = 0; // Within epsilon of the inside shell
                } else {
                    continue; // Too close inside
                }
            } else if (R >= rv[nR - 1]) {
                if (fabs(R - rv[nR - 1]) < 1.0E-12 * rv[nR - 1]) {
                    ind = nR - 2; // Within epsilon of the outside shell
                } else {
                    continue; // Too far outside
                }
            } else {
                double dl = lv[0] - lR; // Negative
                double dr = lv[nR - 1] - lR; // Positive
                double xl = 0;
                double xr = nR-1;
                int retl = 0;
                int retr = 0;
                do {
                    double xc = xr - dr * (xr - xl) / (dr - dl);
                    int xind = (int) xc;
                    if (xind == nR - 1 && lv[xind - 1] - lR <= 0.0) {
                        ind = nR - 2;
                        break;
                    } else if (lv[xind] - lR <= 0.0 && lv[xind + 1] - lR > 0.0) {
                        ind = xind;
                        break;
                    }
                    double dc = lv[xind] - lR;
                    if (dc <= 0.0) {
                        xl = (double) xind;
                        dl = dc;
                        retl = 0;
                        retr++;
                    } else {
                        xr = (double) xind;
                        dr = dc;
                        retr = 0;
                        retl++;
                    }
                    if (retr > 1) {
                        dr *= 0.5;
                    } 
                    if (retl > 1) {
                        dl *= 0.5;
                    } 
                } while (true);
            }

            // Exponential interpolation pair
            double Rl = rs_[A][ind];
            double Rr = rs_[A][ind+1];
            double wl = ws_[A][ind];
            double wr = ws_[A][ind+1];

            if (wl == 0.0 || wr == 0.0) continue;

            double lwl = log(wl);
            double lwr = log(wr);

            wA[A] = exp(lwl + (lwr - lwl)*(R - Rl)/(Rr - Rl));

            wT += wA[A];
        }
        if (wT == 0.0) wT = 1.0; // Don't NaN due to no density, yo?
        
        // => Output <= //

        double scale = (rhop == NULL ? 1.0 / wT : rhop[P] / wT);
        if (atom == -1) {
            for (int A = 0; A < nA; A++) {
                wp[A][P] = wA[A] * scale;
            }
        } else {
            wp[0][P] = wA[atom] * scale;
        }
    }

    for (int thread = 0; thread < nthreads; thread++) {
        delete[] wA2[thread];
    }
}
void StockholderDensity::compute_charges(double scale) 
{
    // Where are the true atoms?
    std::vector<int> Aind;
    std::vector<int> Aind2;
    for (int A = 0; A < molecule_->natom(); A++) {
        if (molecule_->Z(A) != 0) {
            Aind2.push_back(Aind.size());
            Aind.push_back(A);
        } else {
            Aind2.push_back(-1);
        }
    }
    int nA = Aind.size();
    int nA2 = molecule_->natom();

    // Where are the atomic charges?
    double* Np = N_->pointer();

    // Print    
    fprintf(outfile,"   > Atomic Charges <\n\n");
    fprintf(outfile,"    %4s %3s %11s %11s %11s\n", 
        "N", "Z", "Nuclear", "Electronic", "Atomic");
    double Ztot;
    double Qtot;
    for (int A = 0; A < nA; A++) {
        int Aabs = Aind[A];
        double Z = molecule_->Z(Aabs);
        double Q = -scale * Np[A];
        fprintf(outfile,"    %4d %3s %11.3E %11.3E %11.3E\n", 
            Aabs+1, molecule_->symbol(Aabs).c_str(), Z, Q, Z + Q);
        Ztot += Z;
        Qtot += Q;
    }
    fprintf(outfile,"    %8s %11.3E %11.3E %11.3E\n", 
            "Total", Ztot, Qtot, Ztot + Qtot);
    fprintf(outfile,"\n");

    fprintf(outfile,"    True Molecular Charge: %11.3E\n", (double) molecule_->molecular_charge());
    fprintf(outfile,"    Grid Molecular Charge: %11.3E\n", Ztot + Qtot);
    fprintf(outfile,"    Grid Error:            %11.3E\n", Ztot + Qtot - (double) molecule_->molecular_charge());
    fprintf(outfile,"\n");

    fflush(outfile);
}
boost::shared_ptr<Matrix> StockholderDensity::charges(double scale) 
{
    // Where are the true atoms?
    std::vector<int> Aind;
    std::vector<int> Aind2;
    for (int A = 0; A < molecule_->natom(); A++) {
        if (molecule_->Z(A) != 0) {
            Aind2.push_back(Aind.size());
            Aind.push_back(A);
        } else {
            Aind2.push_back(-1);
        }
    }
    int nA = Aind.size();
    int nA2 = molecule_->natom();

    // Where are the atomic charges?
    double* Np = N_->pointer();

    boost::shared_ptr<Matrix> T(new Matrix("Q", nA, 1));
    double* Tp = T->pointer()[0];

    for (int A = 0; A < nA; A++) {
        int Aabs = Aind[A];
        double Z = molecule_->Z(Aabs);
        double Q = -scale * Np[A];
        Tp[A] = Z + Q;
    }
    
    return T;
}

}
