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
        points->compute_points(blocks[ind]);
        int npoints = blocks[ind]->npoints();
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
    //grid_->print_details();

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
    std::vector<std::vector<double> > rs;
    std::vector<std::vector<double> > ws;
    rs.resize(nA);
    ws.resize(nA);
    int nstate = 0;
    for (int A = 0; A < nA; A++) {
        int Aabs = Aind[A];
        for (int R = 0; R < rads[Aabs]->npoints(); R++) {
            rs[A].push_back(rads[Aabs]->r()[R]);
            ws[A].push_back(1.0); // Initial guess
            nstate++;
        } 
    }

    // Target
    Q_ = boost::shared_ptr<Matrix>(new Matrix("Q", nA, nP));
    double** Qp = Q_->pointer();

    // Temps 
    double wT;
    double* wA = new double[nA];

    // DIIS Setup
    boost::shared_ptr<Matrix> state(new Matrix("State", nstate, 1));
    boost::shared_ptr<Matrix> error(new Matrix("Error", nstate, 1));
    double** statep = state->pointer();
    double** errorp = error->pointer();

    boost::shared_ptr<DIISManager> diis_manager(new DIISManager(diis_max_vecs_, "ISA DIIS vector", DIISManager::LargestError, DIISManager::OnDisk));
    diis_manager->set_error_vector_size(1, DIISEntry::Matrix, error.get());
    diis_manager->set_vector_size(1, DIISEntry::Matrix, state.get());

    // Store last iteration
    std::vector<std::vector<double> > ws_old = ws;

    // => Master Loop <= //

    fprintf(outfile, "   > ISA Iterations <\n\n");
    fprintf(outfile, "    Convergence:        %11.3E\n", convergence_);
    fprintf(outfile, "    Maximum iterations: %11d\n", maxiter_);   
    fprintf(outfile, "    DIIS:               %11s\n", (diis_ ? "Yes" : "No"));
    fprintf(outfile, "    DIIS Min Vecs:      %11d\n", diis_min_vecs_);
    fprintf(outfile, "    DIIS Max Vecs:      %11d\n", diis_max_vecs_);
    fprintf(outfile, "\n");
    fflush(outfile);

    bool converged = false;
    for (int iter = 0; iter < maxiter_; iter++) {

        // Compute Q_A^P via W_A^P
        for (int P = 0; P < nP; P++) {
            double xc = xp[P];
            double yc = yp[P];
            double zc = zp[P];
            double rhoc = rhop[P];
            wT = 0.0;
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

                // What is the closest radial node? (could be faster)
                double dR = std::numeric_limits<double>::max();
                int ind = 0;
                for (int k = 0; k < rs[A].size(); k++) {
                    double d2R = fabs(rs[A][k] - R);
                    if (d2R < dR) {
                        dR = d2R;
                        ind = k;
                    }
                }

                // Edge cases (and unique ind, ind+1 terp pair)
                if (ind == 0 && (R - rs[A][0]) * (R - rs[A][1]) > 0.0) continue; // Too close inside
                if (ind == rs[A].size() - 1) {
                    if ((R - rs[A][rs[A].size()-1]) * (R - rs[A][rs[A].size()-2]) > 0.0) { 
                        continue; // Too far outside
                    } else {
                        ind--; // Decrement to use last two indices
                    }
                }

                // Exponential interpolation pair
                double Rl = rs[A][ind];
                double Rr = rs[A][ind+1];
                double wl = ws[A][ind];
                double wr = ws[A][ind+1];

                if (wl == 0.0 || wr == 0.0) continue;

                double lwl = log(wl);
                double lwr = log(wr);

                wA[A] = exp(lwl + (lwr - lwl)*(R - Rl)/(Rr - Rl));

                wT += wA[A];
            }
            if (wT == 0.0) wT = 1.0; // Don't NaN due to no density, yo?
            for (int A = 0; A < nA; A++) {
                Qp[A][P] = wA[A] / wT * rhoc;
            }
        }

        // Debug printing (only for isotropic radial counts)
        //boost::shared_ptr<Matrix> wtemp(new Matrix("W", ws.size(), ws[0].size()));
        //double** wtempp = wtemp->pointer();
        //for (int A = 0; A < ws.size(); A++) {
        //    for (int R = 0; R < ws[0].size(); R++) {    
        //        wtempp[A][R] = ws[A][R];
        //    }
        //}
        //wtemp->print();
        
        // Debug printing
        //Q_->print();

        // Perform spherical averaging
        int offset = 0;
        for (int A = 0; A < nA2; A++) {
            for (int R = 0; R < spheres[A].size(); R++) {
                if (molecule_->Z(A) != 0) {
                    int Arel = Aind2[A]; 
                    double val = 0.0;
                    for (int i = 0; i < spheres[A][R]->npoints(); i++) {
                        val += Qp[Arel][i + offset];
                    }
                    val /= spheres[A][R]->npoints();
                    ws[Arel][R] = val; 
                }
                offset += spheres[A][R]->npoints();
            }
        }

        // Residual and error norm
        double norm = 0.0;
        std::vector<std::vector<double> > dws = ws;
        for (int A = 0; A < ws.size(); A++) {
            for (int R = 0; R < ws[A].size(); R++) {
                dws[A][R] = ws[A][R] - ws_old[A][R];
                norm += abs(dws[A][R]) * rs[A][R] * rs[A][R];
            }
        }

        // Checkpoint last iteration
        ws_old = ws;
    
        // Print iterative trace
        fprintf(outfile,"    @ISA Iter %4d %24.16E ", iter, norm);
        fflush(outfile);

        // Convergence check
        if (norm < convergence_) { 
            fprintf(outfile,"\n");
            converged = true;
            break; 
        }

        // DIIS Add
        if (diis_ && iter > 0) {
            int offset2 = 0;
            for (int A = 0; A < ws.size(); A++) {
                for (int R = 0; R < ws[A].size(); R++) {
                    statep[0][offset2] = ws[A][R];
                    errorp[0][offset2] = dws[A][R];
                    offset2++;
                }
            }
            diis_manager->add_entry(2,error.get(),state.get());
        }

        // DIIS Extrapolate
        if (diis_ && iter >= diis_min_vecs_) {
            diis_manager->extrapolate(1,state.get());
            int offset2 = 0;
            for (int A = 0; A < ws.size(); A++) {
                for (int R = 0; R < ws[A].size(); R++) {
                    ws[A][R] = statep[0][offset2];
                    offset2++;
                }
            }
            fprintf(outfile,"DIIS");
        }

        fprintf(outfile,"\n");
        fflush(outfile);
 
    }

    delete[] wA;

    // Store targets
    rs_ = rs;
    ws_ = ws;

    fprintf(outfile,"\n");
    if (converged) { 
        fprintf(outfile,"    ISA Converged.\n\n"); 
    } else {
        fprintf(outfile,"    ISA Failed.\n\n"); 
    }
    fflush(outfile);
}
void StockholderDensity::compute_weights(int nP, double* xp, double* yp, double* zp, double** wp, double* rhop)
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
    
    // Temps 
    double wT;
    double* wA = new double[nA];

    // Compute Q_A^P via W_A^P
    for (int P = 0; P < nP; P++) {
        double xc = xp[P];
        double yc = yp[P];
        double zc = zp[P];
        wT = 0.0;
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

            // What is the closest radial node? (could be faster)
            double dR = std::numeric_limits<double>::max();
            int ind = 0;
            for (int k = 0; k < rs_[A].size(); k++) {
                double d2R = fabs(rs_[A][k] - R);
                if (d2R < dR) {
                    dR = d2R;
                    ind = k;
                }
            }

            // Edge cases (and unique ind, ind+1 terp pair)
            if (ind == 0 && (R - rs_[A][0]) * (R - rs_[A][1]) > 0.0) continue; // Too close inside
            if (ind == rs_[A].size() - 1) {
                if ((R - rs_[A][rs_[A].size()-1]) * (R - rs_[A][rs_[A].size()-2]) > 0.0) { 
                    continue; // Too far outside
                } else {
                    ind--; // Decrement to use last two indices
                }
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
        if (rhop != NULL) {
            for (int A = 0; A < nA; A++) {
                wp[A][P] = wA[A] / wT * rhop[P];
            }
        } else {
            for (int A = 0; A < nA; A++) {
                wp[A][P] = wA[A] / wT;
            }
        }
    }
    delete[] wA;
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

    // Target
    boost::shared_ptr<Vector> Q2(new Vector("Q2", nA));
    double* Q2p = Q2->pointer();

    int max_points = grid_->max_points();
    boost::shared_ptr<Matrix> v(new Matrix("v", nA, max_points));    
    double** vp = v->pointer();
   
    int npoints = rho_->dimpi()[0];
    double* rhop = rho_->pointer();
    double* xp = x_->pointer(); 
    double* yp = y_->pointer(); 
    double* zp = z_->pointer(); 
    double* wp = w_->pointer(); 

    for (int index = 0; index < npoints; index+=max_points) {
        int points = (index + max_points >= npoints ? npoints - index : max_points);
        compute_weights(points,&xp[index],&yp[index],&zp[index],vp,&rhop[index]);
        for (int A = 0; A < nA; A++) {
            Q2p[A] += C_DDOT(points,&wp[index],1,vp[A],1);
        }
    }
    Q2->scale(-scale);

    // Print    
    fprintf(outfile,"   > Atomic Charges <\n\n");
    fprintf(outfile,"    %4s %3s %11s %11s %11s\n", 
        "N", "Z", "Nuclear", "Electronic", "Atomic");
    double Ztot;
    double Qtot;
    for (int A = 0; A < nA; A++) {
        int Aabs = Aind[A];
        double Z = molecule_->Z(Aabs);
        double Q = Q2p[A];
        fprintf(outfile,"    %4d %3s %11.3E %11.3E %11.3E\n", 
            Aabs, molecule_->symbol(Aabs).c_str(), Z, Q, Z + Q);
        Ztot += Z;
        Qtot += Q;
    }
    fprintf(outfile,"    %8s %11.3E %11.3E %11.3E\n", 
            "Total", Ztot, Qtot, Ztot + Qtot);
    fprintf(outfile,"\n");
    fflush(outfile);
}

}
