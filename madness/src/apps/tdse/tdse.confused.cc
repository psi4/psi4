/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id$
*/
/// \file tdse/tdse.confused.cc
/// \brief Evolves the hydrogen atom in imaginary and also real time


#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/qmprop.h>
#include <mra/operator.h>
#include <constants.h>
#include <tensor/vmath.h>

#include <mra/lbdeux.h>

using namespace madness;

struct InputParameters {
  static const int MAXNATOM=99;

    // IF YOU ADD A NEW PARAMETER DON'T FORGET TO INCLUDE IT IN
    // a) read()
    // b) serialize()
    // c) operator<<()

  double L;           // Box size for the simulation
  double Lsmall;      // Box size for small (near nucleus) plots
  double Llarge;      // Box size for large (far from nucleus) plots
  double F;           // Laser field strength
  double omega;       // Laser frequency
  double ncycle;      // Number of laser cycles in envelope
  int natom;          // Number of atoms
  double Z[MAXNATOM]; // Nuclear charge of atoms
  double R[MAXNATOM][3]; // Coordinates of atoms
  int k;              // wavelet order
  double thresh;      // precision for truncating wave function
  double safety;      // additional precision (thresh*safety) for operators and potential
  double cut;         // smoothing parameter for 1/r (same for all atoms for now)
  string prefix;      // Prefix for filenames
  int ndump;          // dump wave function to disk every ndump steps
  int nplot;          // dump opendx plot to disk every nplot steps
  int nprint;         // print stats every nprint steps
  int nloadbal;       // load balance every nloadbal steps
  int nio;            // Number of IO nodes

  double tScale;      // Scaling parameter for optimization

  double target_time;// Target end-time for the simulation

  void read(const char* filename) {
    ifstream f(filename);
    string tag;
    printf("\n");
    printf("       Simulation parameters\n");
    printf("       ---------------------\n");
    while(f >> tag) {
        if (tag[0] == '#') {
            char ch;
            printf("    comment  %s ",tag.c_str());
            while (f.get(ch)) {
                printf("%c",ch);
                if (ch == '\n') break;
            }
        }
        else if (tag == "L") {
            f >> L;
            printf("             L = %.1f\n", L);
        }
        else if (tag == "Lsmall") {
            f >> Lsmall;
            printf("        Lsmall = %.1f\n", Lsmall);
        }
        else if (tag == "Llarge") {
            f >> Llarge;
            printf("        Llarge = %.1f\n", Llarge);
        }
        else if (tag == "F") {
            f >> F;
            printf("             F = %.6f\n", F);
        }
        else if (tag == "omega") {
            f >> omega;
            printf("         omega = %.6f\n", omega);
        }
        else if (tag == "ncycle") {
            f >> ncycle;
            printf("         ncycle = %.6f\n", ncycle);
        }
        else if (tag == "natom") {
            f >> natom;
            printf("         natom = %d\n", natom);
            for (int i=0; i<natom; i++) {
                f >> Z[i] >> R[i][0] >> R[i][1] >> R[i][2];
                printf("           atom %2d   %.1f  %10.6f  %10.6f  %10.6f\n", i, Z[i], R[i][0], R[i][1], R[i][2]);
            }
        }
        else if (tag == "k") {
            f >> k;
            printf("             k = %d\n", k);
        }
        else if (tag == "thresh") {
            f >> thresh;
            printf("        thresh = %.1e\n", thresh);
        }
        else if (tag == "safety") {
            f >> safety;
            printf("        safety = %.1e\n", safety);
        }
        else if (tag == "cut") {
            f >> cut;
            printf("           cut = %.2f\n", cut);
        }
        else if (tag == "prefix") {
            f >> prefix;
            printf("        prefix = %s\n", prefix.c_str());
        }
        else if (tag == "ndump") {
            f >> ndump;
            printf("         ndump = %d\n", ndump);
        }
        else if (tag == "nplot") {
            f >> nplot;
            printf("         nplot = %d\n", nplot);
        }
        else if (tag == "nprint") {
            f >> nprint;
            printf("         nprint = %d\n", nprint);
        }
        else if (tag == "nloadbal") {
            f >> nloadbal;
            printf("         nloadbal = %d\n", nloadbal);
        }
        else if (tag == "nio") {
            f >> nio;
            printf("           nio = %d\n", nio);
        }
        else if (tag == "target_time") {
            f >> target_time;
            printf("   target_time = %.3f\n", target_time);
        }
        else if (tag == "tScale") {
            f >> tScale;
            printf("           tScale = %.5f\n", tScale);
        }
        else {
            MADNESS_EXCEPTION("unknown input option", 0);
        }
    }
  }

  template <typename Archive>
  void serialize(Archive & ar) {
    ar & L & Lsmall & Llarge & F & omega & ncycle & natom & Z;
    ar & archive::wrap(&(R[0][0]), 3*MAXNATOM);
    ar & k & thresh & safety & cut & prefix & ndump & nplot & nprint & nloadbal & nio & target_time &
		tScale;
  }
};

ostream& operator<<(ostream& s, const InputParameters& p) {
    s << p.L<< " " << p.Lsmall<< " " << p.Llarge<< " " << p.F<< " " << p.omega<<
        " " << p.ncycle << " " << p.Z<< " " << p.R[0]<< " " << p.k<< " " <<
        p.thresh<< " " << p.cut<< " " << p.prefix<< " " << p.ndump<< " " <<
        p.nplot << " " << p.nprint << " "  << p.nloadbal << " " << p.nio << p.tScale << std::endl;
return s;
}

InputParameters param;

static double zero_field_time;      // Laser actually switches on after this time (set by propagate)
                                    // Delay provides for several steps with no field before start

// typedefs to make life less verbose
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr< FunctionFunctorInterface<double_complex,3> > complex_functorT;
typedef Function<double_complex,3> complex_functionT;
typedef FunctionFactory<double_complex,3> complex_factoryT;
typedef SeparatedConvolution<double_complex,3> complex_operatorT;
typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;

// This controls the distribution of data across the machine
class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};

    LevelPmap(World& world) : nproc(world.nproc()) {}

    // Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;

        // This randomly hashes levels 0-2 and then
        // hashes nodes by their grand-parent key so as
        // to increase locality separately on each level.
        //if (n <= 2) hash = key.hash();
        //else hash = key.parent(2).hash();

        // This randomly hashes levels 0-3 and then
        // maps nodes on even levels to the same
        // random node as their parent.
        // if (n <= 3 || (n&0x1)) hash = key.hash();
        // else hash = key.parent().hash();

        // This randomly hashes each key
        hash = key.hash();

        return hash%nproc;
    }
};

// Derivative of the smoothed 1/r approximation

// Invoke as \c du(r/c)/(c*c) where \c c is the radius of the smoothed volume.
static double d_smoothed_potential(double r) {
    double r2 = r*r;

    if (r > 6.5) {
        return -1.0/r2;
    }
    else if (r > 1e-2) {
        return -(1.1283791670955126*(0.88622692545275800*erf(r)-exp(-r2)*r*(1.0-r2)))/r2;
    }
    else {
        return (-1.880631945159187623160265+(1.579730833933717603454623-0.7253866074185437975046736*r2)*r2)*r;
    }
}

// Smoothed 1/r potential.

// Invoke as \c u(r/c)/c where \c c is the radius of the smoothed volume.
static double smoothed_potential(double r) {
    double r2 = r*r, pot;
    if (r > 6.5){
        pot = 1.0/r;
    } else if (r > 1e-2) {
        pot = erf(r)/r + exp(-r2)*0.56418958354775630;
    } else{
        pot = 1.6925687506432689-r2*(0.94031597257959381-r2*(0.39493270848342941-0.12089776790309064*r2));
    }

    return pot;
}

// Nuclear attraction potential
static double V(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    double sum = 0.0;
    for (int i=0; i<param.natom; i++) {
      double xx = x-param.R[i][0];
      double yy = y-param.R[i][1];
      double zz = z-param.R[i][2];
      double rr = sqrt(xx*xx+yy*yy+zz*zz);

      sum +=  -param.Z[i]*smoothed_potential(rr/param.cut)/param.cut;
    }
    return sum;
}

// dV/dx
static double dVdx(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    double sum = 0.0;
    for (int i=0; i<param.natom; i++) {
      double xx = x-param.R[i][0];
      double yy = y-param.R[i][1];
      double zz = z-param.R[i][2];
      if (xx) {
	double rr = sqrt(xx*xx+yy*yy+zz*zz);
	sum += -param.Z[i]*xx*d_smoothed_potential(rr/param.cut)/(rr*param.cut*param.cut);
      }
    }
    return sum;
}

// dV/dy
static double dVdy(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    double sum = 0.0;
    for (int i=0; i<param.natom; i++) {
      double xx = x-param.R[i][0];
      double yy = y-param.R[i][1];
      double zz = z-param.R[i][2];
      if (yy) {
	double rr = sqrt(xx*xx+yy*yy+zz*zz);
	sum += -param.Z[i]*yy*d_smoothed_potential(rr/param.cut)/(rr*param.cut*param.cut);
      }
    }
    return sum;
}

// dV/dz
static double dVdz(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    double sum = 0.0;
    for (int i=0; i<param.natom; i++) {
      double xx = x-param.R[i][0];
      double yy = y-param.R[i][1];
      double zz = z-param.R[i][2];
      if (zz) {
	double rr = sqrt(xx*xx+yy*yy+zz*zz);
	sum += -param.Z[i]*zz*d_smoothed_potential(rr/param.cut)/(rr*param.cut*param.cut);
      }
    }
    return sum;
}

// Initial guess wave function using symmetric superposition of 1s orbital on atoms
static double guess(const coordT& r) {
    const double x=r[0], y=r[1], z=r[2];
    double sum = 0.0;
    for (int i=0; i<param.natom; i++) {
      double xx = x-param.R[i][0];
      double yy = y-param.R[i][1];
      double zz = z-param.R[i][2];
      double rr = sqrt(xx*xx+yy*yy+zz*zz+param.cut*param.cut);
      sum += sqrt(3.14*param.Z[i]*param.Z[i]*param.Z[i])*exp(-param.Z[i]*rr);
    }
    return sum;
}

// x-dipole
double xdipole(const coordT& r) {
    return r[0];
}

// y-dipole
double ydipole(const coordT& r) {
    return r[1];
}

// z-dipole
double zdipole(const coordT& r) {
    return r[2];
}

// Strength of the laser field at time t
double laser(double t) {
    double omegat = param.omega*t;

    if (omegat < 0.0 || omegat/(2*param.ncycle) > constants::pi) return 0.0;

    double envelope = sin(omegat/(2*param.ncycle));
    envelope *= envelope;
    return param.F*envelope*sin(omegat);
}

double myreal(double t) {return t;}

double myreal(const double_complex& t) {return real(t);}

// Given psi and V evaluate the energy ... leaves psi compressed, potn reconstructed
template <typename T>
double energy(World& world, const Function<T,3>& psi, const functionT& potn) {
    // First do all work in the scaling function basis
    bool DOFENCE = false;
    psi.reconstruct();
    Function<T,3> dx = diff(psi,0,DOFENCE);
    Function<T,3> dy = diff(psi,1,DOFENCE);
    Function<T,3> dz = diff(psi,2,DOFENCE);
    Function<T,3> Vpsi = psi*potn;

    // Now do all work in the wavelet basis
    psi.compress(DOFENCE); Vpsi.compress(DOFENCE); dx.compress(DOFENCE); dy.compress(DOFENCE); dz.compress(true);
    T S = psi.inner(psi);
    T PE = psi.inner(Vpsi);
    T KE = 0.5*(inner(dx,dx) + inner(dy,dy) + inner(dz,dz));
    T E = (KE+PE)/S;

    dx.clear(); dy.clear(); dz.clear(); Vpsi.clear(); // To free memory on return
    world.gop.fence();
//     if (world.rank() == 0) {
//         print("the overlap integral is",S);
//         print("the kinetic energy integral",KE);
//         print("the potential energy integral",PE);
//         print("the total energy",E);
//     }

    return myreal(E);
}

void converge(World& world, functionT& potn, functionT& psi, double& eps) {

    for (int iter=0; iter<30; iter++) {
        if (world.rank() == 0) print("beginning iter", iter, wall_time());
        operatorT op = BSHOperator3D<double>(world, sqrt(-2*eps), param.k, param.cut, param.thresh);
        functionT Vpsi = (potn*psi);
        if (world.rank() == 0) print("made V*psi", wall_time());
        Vpsi.scale(-2.0).truncate();
        if (world.rank() == 0) print("tryuncated V*psi", wall_time());
        functionT tmp = apply(op,Vpsi).truncate(param.thresh);
        if (world.rank() == 0) print("applied operator", wall_time());
        double norm = tmp.norm2();
        functionT r = tmp-psi;
        double rnorm = r.norm2();
        double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
        if (world.rank() == 0) {
            print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        }
        psi = tmp.scale(1.0/norm);
        eps = eps_new;
        if (rnorm < std::max(1e-5,param.thresh)) break;
    }
}

complex_functionT chin_chen(const complex_functionT& expV_0,
                            const complex_functionT& expV_tilde,
                            const complex_functionT& expV_1,
                            const complex_operatorT& G,
                            const complex_functionT& psi0) {
    // psi(t) = exp(-i*V(t)*t/6) exp(-i*T*t/2) exp(-i*2*Vtilde(t/2)*t/3) exp(-i*T*t/2) exp(-i*V(0)*t/6)
    // .             expV_1            G               expV_tilde             G             expV_0

    complex_functionT psi1;

    double t0 = wall_time();
    psi1 = expV_0*psi0;     psi1.truncate();
    double t1 = wall_time();
    psi1 = apply(G,psi1);   psi1.truncate();
    double t2 = wall_time();
    psi1 = expV_tilde*psi1; psi1.truncate();
    double t3 = wall_time();
    psi1 = apply(G,psi1);   psi1.truncate();
    double t4 = wall_time();
    psi1 = expV_1*psi1;     psi1.truncate(param.thresh);
    double t5 = wall_time();

    if (psi1.world().rank() == 0) {
        printf("chin-chen: %.2f %.2f %.2f %.2f %.2f\n", t1-t0, t2-t1, t3-t2, t4-t3, t5-t4);
    }

    return psi1;
}

complex_functionT trotter(World& world,
                          const complex_functionT& expV,
                          const complex_operatorT& G,
                          const complex_functionT& psi0) {
    //    psi(t) = exp(-i*T*t/2) exp(-i*V(t/2)*t) exp(-i*T*t/2) psi(0)

    complex_functionT psi1;

    unsigned long size = psi0.size();
    if (world.rank() == 0) print("APPLYING G", size);
    psi1 = apply(G,psi0);  psi1.truncate();  size = psi1.size();
    if (world.rank() == 0) print("APPLYING expV", size);
    psi1 = expV*psi1;      psi1.truncate();  size = psi1.size();
    if (world.rank() == 0) print("APPLYING G again", size);
    psi1 = apply(G,psi1);  psi1.truncate(param.thresh);  size = psi1.size();
    if (world.rank() == 0) print("DONE", size);

    return psi1;
}

template<typename T, int NDIM>
struct unaryexp {
//     void operator()(const Key<NDIM>& key, Tensor<T>& t) const {
//         UNARY_OPTIMIZED_ITERATOR(T, t, *_p0 = exp(*_p0););
//     }
//     template <typename Archive>
//     void serialize(Archive& ar) {}
};


template<int NDIM>
struct unaryexp<double_complex,NDIM> {
    void operator()(const Key<NDIM>& key, Tensor<double_complex>& t) const {
        //vzExp(t.size, t.ptr(), t.ptr());
        UNARY_OPTIMIZED_ITERATOR(double_complex, t, *_p0 = exp(*_p0););
    }
    template <typename Archive>
    void serialize(Archive& ar) {}
};


// Returns exp(-I*t*V)
complex_functionT make_exp(double t, const functionT& v) {
    v.reconstruct();
    complex_functionT expV = double_complex(0.0,-t)*v;
    expV.unaryop(unaryexp<double_complex,3>());
    return expV;
}

void print_stats_header(World& world) {
    if (world.rank() == 0) {
        printf("  step       time            field           energy            norm           overlap0         x-dipole         y-dipole         z-dipole           accel      wall-time(s)\n");
        printf("------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ---------------- ------------\n");
    }
}

void print_stats(World& world, int step, double t, const functionT& v,
                 const functionT& x, const functionT& y, const functionT& z,
                 const functionT& dV_dz,
                 const complex_functionT& psi0, const complex_functionT& psi) {
    double norm = psi.norm2();
    double current_energy = energy(world, psi, v);
    double xdip = real(inner(psi, x*psi))/(norm*norm);
    double ydip = real(inner(psi, y*psi))/(norm*norm);
    double zdip = real(inner(psi, z*psi))/(norm*norm);
    double overlap0 = std::abs(psi.inner(psi0))/norm;
    double accel = std::abs(psi.inner(psi*dV_dz))/(norm*norm);
    if (world.rank() == 0) {
        printf("%7d %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %16.8e %9.1f\n", step, t, laser(t), current_energy, norm, overlap0, xdip, ydip, zdip, accel, wall_time());
    }
}

const char* wave_function_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5d", param.prefix.c_str(), step);
    return fname;
}

const char* wave_function_small_plot_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5dS.dx", param.prefix.c_str(), step);
    return fname;
}

const char* wave_function_large_plot_filename(int step) {
    static char fname[1024];
    sprintf(fname, "%s-%5.5dL.dx", param.prefix.c_str(), step);
    return fname;
}

complex_functionT wave_function_load(World& world, int step) {
    complex_functionT psi;
    ParallelInputArchive ar(world, wave_function_filename(step));
    ar & psi;
    return psi;
}

void wave_function_store(World& world, int step, const complex_functionT& psi) {
    ParallelOutputArchive ar(world, wave_function_filename(step), param.nio);
    ar & psi;
}

bool wave_function_exists(World& world, int step) {
    return ParallelInputArchive::exists(world, wave_function_filename(step));
}

void doplot(World& world, int step, const complex_functionT& psi, double Lplot, long numpt, const char* fname) {
    double start = wall_time();
    Tensor<double> cell(3,2);
    std::vector<long> npt(3, numpt);
    cell(_,0) = -Lplot;
    cell(_,1) =  Lplot;
    //plotdx(psi, fname, cell, npt);
    if (world.rank() == 0) print("plotting used", wall_time()-start);
}


void line_plot(World& world, int step, complex_functionT& psi) {
    static const int npt = 10001;
    double_complex v[10001];
    psi.reconstruct();
    for (int i=0; i<npt; i++)
        v[i] = 0.0;
    for (int i=world.rank(); i<npt; i+=world.size()) {
        double z = -param.Llarge + 2.0*i*param.Llarge/(npt-1);
        coordT r(0.0);
        r[2] = z;
        v[i] = psi.eval(r);
    }
    world.gop.fence();
    world.gop.sum(v, npt);
    if (world.rank() == 0) {
        char buf[256];
        sprintf(buf, "%s.lineplot", wave_function_filename(step));
        ofstream f(buf);
        f.precision(10);
        for (int i=0; i<npt; i++) {
            double z = -param.Llarge + 2.0*i*param.Llarge/(npt-1);
            f << z << " " << v[i] << "\n";
        }
    }
    world.gop.fence();
}

template <typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=1.0) : leaf_value(leaf_value), parent_value(parent_value) {}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (key.level() <= 1) return 1000.0;
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
        //return key.level()+1.0;
    }
};

void loadbal(World& world,
             functionT& potn, functionT& vt,
             functionT& dpotn_dx, functionT& dpotn_dy, functionT& dpotn_dz,
             functionT& dpotn_dx_sq, functionT& dpotn_dy_sq,
             complex_functionT& psi, complex_functionT& psi0,
             functionT& x, functionT& y, functionT& z) {
    if (world.size() < 2) return;
    if (world.rank() == 0) print("starting LB");
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(potn, lbcost<double,3>(1.0,1.0));
    lb.add_tree(psi, lbcost<double_complex,3>(10.0,10.0));
    FunctionDefaults<3>::set_pmap(lb.load_balance());
    world.gop.fence();
    if (world.rank() == 0) print("starting LB copies");
    world.gop.fence();
    potn = copy(potn, FunctionDefaults<3>::get_pmap(), true);
    vt = copy(vt, FunctionDefaults<3>::get_pmap(), true);
    dpotn_dx = copy(dpotn_dx, FunctionDefaults<3>::get_pmap(), true);
    dpotn_dy = copy(dpotn_dy, FunctionDefaults<3>::get_pmap(), true);
    dpotn_dz = copy(dpotn_dz, FunctionDefaults<3>::get_pmap(), true);
    dpotn_dx_sq = copy(dpotn_dx_sq, FunctionDefaults<3>::get_pmap(), true);
    dpotn_dy_sq = copy(dpotn_dy_sq, FunctionDefaults<3>::get_pmap(), true);
    psi = copy(psi, FunctionDefaults<3>::get_pmap(), true);
    psi0 = copy(psi0, FunctionDefaults<3>::get_pmap(), true);
    x = copy(x, FunctionDefaults<3>::get_pmap(), true);
    y = copy(y, FunctionDefaults<3>::get_pmap(), true);
    z = copy(z, FunctionDefaults<3>::get_pmap(), true);
    world.gop.fence();
    if (world.rank() == 0) print("done with LB");
    world.gop.fence();
}

// Evolve the wave function in real time starting from given time step on disk
void propagate(World& world, int step0) {
    //double ctarget = 10.0/param.cut;                // From Fourier analysis of the potential
    double ctarget = 5.0/param.cut;
    //double c = 1.86*ctarget; // This for 10^5 steps
    double c = 1.72*ctarget;   // This for 10^4 steps
    double tcrit = 2*constants::pi/(c*c);

    double time_step = tcrit * param.tScale;

    zero_field_time = 10.0*time_step;

    int nstep = (param.target_time + zero_field_time)/time_step + 1;

    // Ensure everyone has the same data
    world.gop.broadcast(c);
    world.gop.broadcast(time_step);
    world.gop.broadcast(nstep);

    // Free particle propagator for both Trotter and Chin-Chen --- exp(-I*T*time_step/2)
    SeparatedConvolution<double_complex,3> G = qm_free_particle_propagator<3>(world, param.k, c, 0.5*time_step, 2*param.L);
    //G.doleaves = true;

    // The time-independent part of the potential plus derivatives for
    // Chin-Chen and also for computing the power spectrum ... compute
    // derivatives analytically to reduce numerical noise
    functionT potn = factoryT(world).f(V);         potn.truncate();
    functionT dpotn_dx = factoryT(world).f(dVdx);
    functionT dpotn_dy = factoryT(world).f(dVdy);
    functionT dpotn_dz = factoryT(world).f(dVdz);

    functionT dpotn_dx_sq = dpotn_dx*dpotn_dx;
    functionT dpotn_dy_sq = dpotn_dy*dpotn_dy;

    // Dipole moment functions for laser field and for printing statistics
    functionT x = factoryT(world).f(xdipole);
    functionT y = factoryT(world).f(ydipole);
    functionT z = factoryT(world).f(zdipole);

    // Wave function at time t=0 for printing statistics
    complex_functionT psi0 = wave_function_load(world, 0);

    int step = step0;  // The current step
    double t = step0 * time_step - zero_field_time;        // The current time
    complex_functionT psi = wave_function_load(world, step); // The wave function at time t
    functionT vt = potn+laser(t)*z; // The total potential at time t

    if (world.rank() == 0) {
        printf("\n");
        printf("        Evolution parameters\n");
        printf("       --------------------\n");
        printf("     bandlimit = %.2f\n", ctarget);
        printf(" eff-bandlimit = %.2f\n", c);
        printf("         tcrit = %.6f\n", tcrit);
        printf("     time step = %.6f\n", time_step);
        printf(" no field time = %.6f\n", zero_field_time);
        printf("   target time = %.2f\n", param.target_time);
        printf("         nstep = %d\n", nstep);
        printf("\n");
        printf("  restart step = %d\n", step0);
        printf("  restart time = %.6f\n", t);
        printf("\n");
    }

    print_stats_header(world);
    {
        // Make gradient of potential at time t in z direction to compute HHG
        functionT dV_dz = copy(dpotn_dz);
        dV_dz.add_scalar(laser(t));

        print_stats(world, step, t, vt, x, y, z, dV_dz, psi0, psi);
        line_plot(world, step, psi);
    }
    world.gop.fence();

    psi.truncate();

    bool use_trotter = false;
    while (step < nstep) {
        double t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;
        t0 = wall_time();
        if (step < 2 || (step%param.nloadbal) == 0)
            loadbal(world, potn, vt, dpotn_dx, dpotn_dy, dpotn_dz, dpotn_dx_sq, dpotn_dy_sq, psi, psi0, x, y, z);
        t1 = wall_time();

        long depth = psi.max_depth(); long size=psi.size();
        if (use_trotter) {
            // Make the potential at time t + step/2
            functionT vhalf = potn + laser(t+0.5*time_step)*z;

            // Apply Trotter to advance from time t to time t+step
            complex_functionT expV = make_exp(time_step, vhalf);
            psi = trotter(world, expV, G, psi);
        }
        else { // Chin-Chen
            // Make z-component of del V at time tstep/2

            functionT dV_dz = copy(dpotn_dz);
            t2 = wall_time();
            dV_dz.add_scalar(laser(t+0.5*time_step));

            // Make Vtilde at time tstep/2
            functionT Vtilde = potn + laser(t+0.5*time_step)*z;
            t3 = wall_time();

            //functionT dvsq = dpotn_dx_sq + dpotn_dy_sq + dV_dz*dV_dz;
            world.gop.fence();
            functionT dV_dz_sq = dV_dz*dV_dz;
            world.gop.fence();
            dV_dz_sq.compress();
            world.gop.fence();
            functionT dvsq = dpotn_dx_sq + dpotn_dy_sq + dV_dz_sq;
            world.gop.fence();
            dV_dz_sq.clear();

            t4 = wall_time();
            Vtilde.gaxpy(1.0, dvsq, -time_step*time_step/48.0);
            t5 = wall_time();

            // Exponentiate potentials
            complex_functionT expv_0     = make_exp(time_step/6.0, vt);
            t6 = wall_time();
            complex_functionT expv_tilde = make_exp(2.0*time_step/3.0, Vtilde);
            t7 = wall_time();
            complex_functionT expv_1     = make_exp(time_step/6.0, potn + laser(t+time_step)*z);
            t8 = wall_time();

            // Free up some memory
            dV_dz.clear();
            Vtilde.clear();
            dvsq.clear();
            world.gop.fence();
            t9 = wall_time();

            // Apply Chin-Chen
            psi = chin_chen(expv_0, expv_tilde, expv_1, G, psi);
            t10 = wall_time();
        }

        // Update counters, print info, dump/plot as necessary
        step++;
        t += time_step;
        vt = potn+laser(t)*z;


        {
            // Make gradient of potential at time t in z direction to compute HHG
            functionT dV_dz = copy(dpotn_dz);
            t11 = wall_time();
            dV_dz.add_scalar(laser(t));
            t12 = wall_time();

            if ((step%param.nprint) == 0 || step==nstep) {
                print_stats(world, step, t, vt, x, y, z, dV_dz, psi0, psi);
                line_plot(world, step, psi);
                if (world.rank() == 0) print(step, "depth", depth, "size", size);
            }

            t13 = wall_time();
            if (world.rank() == 0)
                printf("loadbal=%.2f copy=%.2f Vtil=%.2f dvsq =%.2f gaxpy=%.2f exp1=%.2f exp2=%.2f exp3=%.2f clear=%.2f CC=%.2f copy=%.2f addscl=%.2f prnt=%.2f\n",
                       t1-t0, t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6, t8-t7, t9-t8, t10-t9, t11-t10, t12-t11, t13-t12);
        }

        if ((step%param.ndump) == 0 || step==nstep) {
            double start = wall_time();
            wave_function_store(world, step, psi);
            // Update the restart file for automatic restarting
            if (world.rank() == 0) {
                std::ofstream("restart") << step << std::endl;
                print("dumping took", wall_time()-start);
            }
            world.gop.fence();
        }

        if ((step%param.nplot) == 0 || step==nstep) {
            doplot(world, step, psi, param.Lsmall, 201, wave_function_small_plot_filename(step));
            doplot(world, step, psi, param.Llarge, 201, wave_function_large_plot_filename(step));
        }
    }
}

void doit(World& world) {
    cout.precision(8);

    if (world.rank() == 0) param.read("input");
    world.gop.broadcast_serializable(param, 0);

    FunctionDefaults<3>::set_k(param.k);                        // Wavelet order
    FunctionDefaults<3>::set_thresh(param.thresh*param.safety);       // Accuracy
    FunctionDefaults<3>::set_initial_level(4);
    FunctionDefaults<3>::set_cubic_cell(-param.L,param.L);
    FunctionDefaults<3>::set_apply_randomize(true);
    FunctionDefaults<3>::set_autorefine(false);
    FunctionDefaults<3>::set_truncate_mode(2);
    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));

    // Make the potential
    functionT potn = factoryT(world).f(V);  potn.truncate();

    // Before doing anything else use the potential to guide future
    // data distribution
    LoadBalanceDeux<3> lb(world);
    lb.add_tree(potn, lbcost<double,3>());
    FunctionDefaults<3>::set_pmap(lb.load_balance());
    world.gop.fence();
    potn = copy(potn, FunctionDefaults<3>::get_pmap());

    // Read restart information
    int step0;               // Initial time step ... filenames are <prefix>-<step0>
    if (world.rank() == 0) std::ifstream("restart") >> step0;
    world.gop.broadcast(step0);

    bool exists = wave_function_exists(world, step0);

    if (!exists) {
        if (step0 == 0) {
            if (world.rank() == 0) print("Computing initial ground state wavefunction", wall_time());
            functionT psi = factoryT(world).f(guess);
            psi.scale(1.0/psi.norm2());
            psi.truncate();
            psi.scale(1.0/psi.norm2());
            if (world.rank() == 0) print("got psi", wall_time());

            double eps = energy(world, psi, potn);
            if (world.rank() == 0) print("guess energy", eps, wall_time());
            converge(world, potn, psi, eps);

            psi.truncate(param.thresh);

            complex_functionT psic = double_complex(1.0,0.0)*psi;
            wave_function_store(world, 0, psic);
        }
        else {
            if (world.rank() == 0) {
                print("The requested restart was not found ---", step0);
                error("restart failed", 0);
            }
            world.gop.fence();
        }
    }

    propagate(world, step0);
}

int main(int argc, char** argv) {
    initialize(argc,argv);
    World world(MPI::COMM_WORLD);

    startup(world,argc,argv);

    try {
        doit(world);
    } catch (const MPI::Exception& e) {
        //print(e); std::cout.flush();
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e); std::cout.flush();
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e); std::cout.flush();
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (char* s) {
        print(s); std::cout.flush();
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s); std::cout.flush();
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what()); std::cout.flush();
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }


    world.gop.fence();

    ThreadPool::end();
    print_stats(world);
    finalize();
    return 0;
}


