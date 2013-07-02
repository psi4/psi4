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

/// \file moldft.cc
/// \brief Molecular HF and DFT code

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <mra/mra.h>
#include <mra/lbdeux.h>
#include <misc/ran.h>
#include <linalg/solvers.h>
#include <ctime>
#include <list>
using namespace madness;

extern int x_rks_s__(const double *rho, double *f, double *dfdra);
extern int c_rks_vwn5__(const double *rho, double *f, double *dfdra);
// extern int x_uks_s__(const double *rho, double *f, double *dfdra);
// extern int c_uks_vwn5__(const double *rho, double *f, double *dfdra);

#include <moldft/molecule.h>
#include <moldft/molecularbasis.h>
#include <moldft/corepotential.h>

/// Simple (?) version of BLAS-1 DROT(N, DX, INCX, DY, INCY, DC, DS)
void drot(long n, double* restrict a, double* restrict b, double s, double c, long inc) {
    if (inc == 1) {
        for (long i=0; i<n; ++i) {
            double aa = a[i]*c - b[i]*s;
            double bb = b[i]*c + a[i]*s;
            a[i] = aa;
            b[i] = bb;
        }
    }
    else {
        for (long i=0; i<(n*inc); i+=inc) {
            double aa = a[i]*c - b[i]*s;
            double bb = b[i]*c + a[i]*s;
            a[i] = aa;
            b[i] = bb;
        }
    }
}

void drot3(long n, double* restrict a, double* restrict b, double s, double c, long inc) {
    if (inc == 1) {
        n*=3;
        for (long i=0; i<n; i+=3) {
            double aa0 = a[i  ]*c - b[i  ]*s;
            double bb0 = b[i  ]*c + a[i  ]*s;
            double aa1 = a[i+1]*c - b[i+1]*s;
            double bb1 = b[i+1]*c + a[i+1]*s;
            double aa2 = a[i+2]*c - b[i+2]*s;
            double bb2 = b[i+2]*c + a[i+2]*s;
            a[i  ] = aa0;
            b[i  ] = bb0;
            a[i+1] = aa1;
            b[i+1] = bb1;
            a[i+2] = aa2;
            b[i+2] = bb2;
        }
    }
    else {
        inc*=3;
        n*=inc;
        for (long i=0; i<n; i+=inc) {
            double aa0 = a[i  ]*c - b[i  ]*s;
            double bb0 = b[i  ]*c + a[i  ]*s;
            double aa1 = a[i+1]*c - b[i+1]*s;
            double bb1 = b[i+1]*c + a[i+1]*s;
            double aa2 = a[i+2]*c - b[i+2]*s;
            double bb2 = b[i+2]*c + a[i+2]*s;
            a[i  ] = aa0;
            b[i  ] = bb0;
            a[i+1] = aa1;
            b[i+1] = bb1;
            a[i+2] = aa2;
            b[i+2] = bb2;
        }
    }
}

class LevelPmap : public WorldDCPmapInterface< Key<3> > {
private:
    const int nproc;
public:
    LevelPmap() : nproc(0) {};

    LevelPmap(World& world) : nproc(world.nproc()) {}

    /// Find the owner of a given key
    ProcessID owner(const Key<3>& key) const {
        Level n = key.level();
        if (n == 0) return 0;
        hashT hash;
        if (n <= 3 || (n&0x1)) hash = key.hash();
        else hash = key.parent().hash();
        //hashT hash = key.hash();
        return hash%nproc;
    }
};

typedef std::shared_ptr< WorldDCPmapInterface< Key<3> > > pmapT;
typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef std::vector<functionT> vecfuncT;
typedef std::pair<vecfuncT,vecfuncT> pairvecfuncT;
typedef std::vector<pairvecfuncT> subspaceT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;

double ttt, sss;
void START_TIMER(World& world) {
    world.gop.fence(); ttt=wall_time(); sss=cpu_time();
}

void END_TIMER(World& world, const char* msg) {
    ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %20.20s %8.2fs %8.2fs\n", msg, sss, ttt);
}


inline double mask1(double x) {
    /* Iterated first beta function to switch smoothly
       from 0->1 in [0,1].  n iterations produce 2*n-1
       zero derivatives at the end points. Order of polyn
       is 3^n.

       Currently use one iteration so that first deriv.
       is zero at interior boundary and is exactly representable
       by low order multiwavelet without refinement */

    x = (x*x*(3.-2.*x));
    return x;
}

double mask3(const coordT& ruser) {
    coordT rsim;
    user_to_sim(ruser, rsim);
    double x= rsim[0], y=rsim[1], z=rsim[2];
    double lo = 0.0625, hi = 1.0-lo, result = 1.0;
    double rlo = 1.0/lo;

    if (x<lo)
        result *= mask1(x*rlo);
    else if (x>hi)
        result *= mask1((1.0-x)*rlo);
    if (y<lo)
        result *= mask1(y*rlo);
    else if (y>hi)
        result *= mask1((1.0-y)*rlo);
    if (z<lo)
        result *= mask1(z*rlo);
    else if (z>hi)
        result *= mask1((1.0-z)*rlo);

    return result;
}

class MolecularPotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
public:
    MolecularPotentialFunctor(const Molecule& molecule)
        : molecule(molecule) {}

    double operator()(const coordT& x) const {
        return molecule.nuclear_attraction_potential(x[0], x[1], x[2]);
    }
};

class MolecularCorePotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
public:
    MolecularCorePotentialFunctor(const Molecule& molecule)
        : molecule(molecule) {}

    double operator()(const coordT& x) const {
        return molecule.molecular_core_potential(x[0], x[1], x[2]);
    }
};

class MolecularGuessDensityFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const AtomicBasisSet& aobasis;
public:
    MolecularGuessDensityFunctor(const Molecule& molecule, const AtomicBasisSet& aobasis)
        : molecule(molecule), aobasis(aobasis) {}

    double operator()(const coordT& x) const {
        return aobasis.eval_guess_density(molecule, x[0], x[1], x[2]);
    }
};


class AtomicBasisFunctor : public FunctionFunctorInterface<double,3> {
private:
    const AtomicBasisFunction aofunc;
    std::vector<coordT> specialpt;
public:
    AtomicBasisFunctor(const AtomicBasisFunction& aofunc)
        : aofunc(aofunc)
    {
    	double x, y, z;
        aofunc.get_coords(x,y,z);
        coordT r;
        r[0]=x; r[1]=y; r[2]=z;
        specialpt=std::vector<coordT>(1,r);
    }

    double operator()(const coordT& x) const {
        return aofunc(x[0], x[1], x[2]);
    }

    std::vector<coordT> special_points() const {return specialpt;}
};

class MolecularDerivativeFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const int atom;
    const int axis;
public:
    MolecularDerivativeFunctor(const Molecule& molecule, int atom, int axis)
        : molecule(molecule), atom(atom), axis(axis) {}

    double operator()(const coordT& x) const {
        return molecule.nuclear_attraction_potential_derivative(atom, axis, x[0], x[1], x[2]);
    }
};

class CorePotentialDerivativeFunctor : public FunctionFunctorInterface<double,3> {
private:
    const Molecule& molecule;
    const int atom;
    const int axis;
public:
    CorePotentialDerivativeFunctor(const Molecule& molecule, int atom, int axis)
        : molecule(molecule), atom(atom), axis(axis) {}

    double operator()(const coordT& r) const {
        return molecule.core_potential_derivative(atom, axis, r[0], r[1], r[2]);
    }
};

class CoreOrbitalFunctor : public FunctionFunctorInterface<double,3> {
    const Molecule molecule;
    const int atom;
    const unsigned int core;
    const int m;
public:
    CoreOrbitalFunctor(Molecule& molecule, int atom, unsigned int core, int m)
        : molecule(molecule), atom(atom), core(core), m(m) {};
    double operator()(const coordT& r) const {
        return molecule.core_eval(atom, core, m, r[0], r[1], r[2]);
    };
};

class CoreOrbitalDerivativeFunctor : public FunctionFunctorInterface<double,3> {
    const Molecule molecule;
    const int atom, axis;
    const unsigned int core;
    const int m;
public:
    CoreOrbitalDerivativeFunctor(Molecule& molecule, int atom, int axis, unsigned int core, int m)
        : molecule(molecule), atom(atom), axis(axis), core(core), m(m) {};
    double operator()(const coordT& r) const {
        return molecule.core_derivative(atom, axis, core, m, r[0], r[1], r[2]);
    };
};


/// A MADNESS functor to compute either x, y, or z
class DipoleFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int axis;
public:
    DipoleFunctor(int axis) : axis(axis) {}
    double operator()(const coordT& x) const {
        return x[axis];
    }
};

double rsquared(const coordT& r) {
    return r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
}

/// A MADNESS functor to compute the cartesian moment x^i * y^j * z^k (i, j, k integer and >= 0)
class MomentFunctor : public FunctionFunctorInterface<double,3> {
private:
    const int i, j, k;
public:
    MomentFunctor(int i, int j, int k) : i(i), j(j), k(k) {}
    MomentFunctor(const std::vector<int>& x) : i(x[0]), j(x[1]), k(x[2]) {}
    double operator()(const coordT& r) const {
        double xi=1.0, yj=1.0, zk=1.0;
        for (int p=0; p<i; ++p) xi *= r[0];
        for (int p=0; p<j; ++p) yj *= r[1];
        for (int p=0; p<k; ++p) zk *= r[2];
        return xi*yj*zk;
    }
};


/// Given overlap matrix, return rotation with 3rd order error to orthonormalize the vectors
tensorT Q3(const tensorT& s) {
    tensorT Q = inner(s,s);
    Q.gaxpy(0.2,s,-2.0/3.0);
    for (int i=0; i<s.dim(0); ++i) Q(i,i) += 1.0;
    return Q.scale(15.0/8.0);
}

/// Computes matrix square root (not used any more?)
tensorT sqrt(const tensorT& s, double tol=1e-8) {
    int n=s.dim(0), m=s.dim(1);
    MADNESS_ASSERT(n==m);
    tensorT c, e;
    //s.gaxpy(0.5,transpose(s),0.5); // Ensure exact symmetry
    syev(s, c, e);
    for (int i=0; i<n; ++i) {
        if (e(i) < -tol) {
            MADNESS_EXCEPTION("Matrix square root: negative eigenvalue",i);
        }
        else if (e(i) < tol) { // Ugh ..
            print("Matrix square root: Warning: small eigenvalue ", i, e(i));
            e(i) = tol;
        }
        e(i) = 1.0/sqrt(e(i));
    }
    for (int j=0; j<n; ++j) {
        for (int i=0; i<n; ++i) {
            c(j,i) *= e(i);
        }
    }
    return c;
}

template <typename T, int NDIM>
struct lbcost {
    double leaf_value;
    double parent_value;
    lbcost(double leaf_value=1.0, double parent_value=0.0) : leaf_value(leaf_value), parent_value(parent_value) {}
    double operator()(const Key<NDIM>& key, const FunctionNode<T,NDIM>& node) const {
        if (key.level() <= 1) {
            return 100.0*(leaf_value+parent_value);
        }
        else if (node.is_leaf()) {
            return leaf_value;
        }
        else {
            return parent_value;
        }
    }
};


struct CalculationParameters {
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // !!!                                                                   !!!
    // !!! If you add more data don't forget to add them to serialize method !!!
    // !!!                                                                   !!!
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    // First list input parameters
    double charge;              ///< Total molecular charge
    double smear;               ///< Smearing parameter
    double econv;               ///< Energy convergence
    double dconv;               ///< Density convergence
    double L;                   ///< User coordinates box size
    double maxrotn;             ///< Step restriction used in autoshift algorithm
    int nvalpha;                ///< Number of alpha virtuals to compute
    int nvbeta;                 ///< Number of beta virtuals to compute
    int nopen;                  ///< Number of unpaired electrons = napha-nbeta
    int maxiter;                ///< Maximum number of iterations
    int nio;                    ///< No. of io servers to use
    bool spin_restricted;       ///< True if spin restricted
    bool lda;                   ///< True if LDA (HF if false)
    int plotlo,plothi;          ///< Range of MOs to print (for both spins if polarized)
    bool plotdens;              ///< If true print the density at convergence
    bool plotcoul;              ///< If true plot the total coulomb potential at convergence
    bool localize;              ///< If true solve for localized orbitals
    bool localize_pm;           ///< If true use PM for localization
    bool restart;               ///< If true restart from orbitals on disk
    unsigned int maxsub;        ///< Size of iterative subspace ... set to 0 or 1 to disable
    int npt_plot;               ///< No. of points to use in each dim for plots
    tensorT plot_cell;          ///< lo hi in each dimension for plotting (default is all space)
    std::string aobasis;        ///< AO basis used for initial guess (6-31g or sto-3g)
    std::string core_type;      ///< core potential type ("" or "mcp")
    bool derivatives;           ///< If true calculate derivatives
    bool dipole;                ///< If true calculatio dipole moment
    bool conv_only_dens;        ///< If true remove bsh_residual from convergence criteria   how ugly name is...
    // Next list inferred parameters
    int nalpha;                 ///< Number of alpha spin electrons
    int nbeta;                  ///< Number of beta  spin electrons
    int nmo_alpha;              ///< Number of alpha spin molecular orbitals
    int nmo_beta;               ///< Number of beta  spin molecular orbitals
    double lo;                  ///< Smallest length scale we need to resolve

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & charge & smear & econv & dconv & L & maxrotn & nvalpha & nvbeta & nopen & maxiter & nio & spin_restricted & lda;
        ar & plotlo & plothi & plotdens & plotcoul & localize & localize_pm & restart & maxsub & npt_plot & plot_cell & aobasis;
        ar & nalpha & nbeta & nmo_alpha & nmo_beta & lo;
        ar & core_type & derivatives & conv_only_dens & dipole;
    }

    CalculationParameters()
        : charge(0.0)
        , smear(0.0)
        , econv(1e-5)
        , dconv(1e-4)
        , L(0.0)
        , maxrotn(0.25)
        , nvalpha(0)
        , nvbeta(0)
        , nopen(0)
        , maxiter(20)
        , nio(1)
        , spin_restricted(true)
        , lda(true)
        , plotlo(0)
        , plothi(-1)
        , plotdens(false)
        , plotcoul(false)
        , localize(true)
        , localize_pm(false)
        , restart(false)
        , maxsub(8)
        , npt_plot(101)
        , aobasis("6-31g")
        , core_type("")
        , derivatives(false)
        , dipole(false)
        , conv_only_dens(false)
        , nalpha(0)
        , nbeta(0)
        , nmo_alpha(0)
        , nmo_beta(0)
        , lo(1e-10) {}

    void read_file(const std::string& filename) {
        std::ifstream f(filename.c_str());
        position_stream(f, "dft");
        std::string s;
        while (f >> s) {
            if (s == "end") {
                break;
            }
            else if (s == "charge") {
                f >> charge;
            }
            else if (s == "smear") {
                f >> smear;
            }
            else if (s == "econv") {
                f >> econv;
            }
            else if (s == "dconv") {
                f >> dconv;
            }
            else if (s == "L") {
                f >> L;
            }
            else if (s == "maxrotn") {
                f >> maxrotn;
            }
            else if (s == "nvalpha") {
                f >> nvalpha;
            }
            else if (s == "nvbeta") {
                f >> nvbeta;
            }
            else if (s == "nopen") {
                f >> nopen;
            }
            else if (s == "unrestricted") {
                spin_restricted = false;
            }
            else if (s == "restricted") {
                spin_restricted = true;
            }
            else if (s == "maxiter") {
                f >> maxiter;
            }
            else if (s == "nio") {
                f >> nio;
            }
            else if (s == "lda") {
                lda = true;
            }
            else if (s == "hf") {
                lda = false;
            }
            else if (s == "plotmos") {
                f >> plotlo >> plothi;
            }
            else if (s == "plotdens") {
                plotdens = true;
            }
            else if (s == "plotcoul") {
                plotcoul = true;
            }
            else if (s == "plotnpt") {
                f >> npt_plot;
            }
            else if (s == "plotcell") {
                plot_cell = tensorT(3L,2L);
                f >> plot_cell(0,0) >> plot_cell(0,1) >> plot_cell(1,0) >> plot_cell(1,1) >> plot_cell(2,0) >> plot_cell(2,1);
            }
            else if (s == "aobasis") {
                f >> aobasis;
                if (aobasis!="sto-3g" && aobasis!="6-31g") {
                    std::cout << "moldft: unrecognized aobasis (sto-3g or 6-31g only): " << aobasis << std::endl;
                    MADNESS_EXCEPTION("input_error", 0);
                }
            }
            else if (s == "canon") {
                localize = false;
            }
            else if (s == "local") {
                localize = true;
            }
            else if (s == "pm") {
                localize_pm = true;
            }
            else if (s == "boys") {
                localize_pm = false;
            }
            else if (s == "restart") {
                restart = true;
            }
            else if (s == "maxsub") {
                f >> maxsub;
                if (maxsub <= 0) maxsub = 1;
                if (maxsub > 20) maxsub = 20;
            }
            else if (s == "core_type") {
                f >> core_type;
            }
            else if (s == "derivatives") {
                derivatives = true;
            }
            else if (s == "dipole") {
                dipole = true;
            }
            else if (s == "convonlydens") {
                conv_only_dens = true;
            }
            else {
                std::cout << "moldft: unrecognized input keyword " << s << std::endl;
                MADNESS_EXCEPTION("input error",0);
            }
            if (nopen != 0) spin_restricted = false;
        }
    }

    void set_molecular_info(const Molecule& molecule, const AtomicBasisSet& aobasis, unsigned int n_core) {
        double z = molecule.total_nuclear_charge();
        int nelec = int(z - charge - n_core*2);
        if (fabs(nelec+charge+n_core*2-z) > 1e-6) {
            error("non-integer number of electrons?", nelec+charge+n_core*2-z);
        }
        nalpha = (nelec + nopen)/2;
        nbeta  = (nelec - nopen)/2;
        if (nalpha < 0) error("negative number of alpha electrons?", nalpha);
        if (nbeta < 0) error("negative number of beta electrons?", nbeta);
        if ((nalpha+nbeta) != nelec) error("nalpha+nbeta != nelec", nalpha+nbeta);
        nmo_alpha = nalpha + nvalpha;
        nmo_beta = nbeta + nvbeta;
        if (nalpha != nbeta) spin_restricted = false;

        // Ensure we have enough basis functions to guess the requested
        // number of states ... a minimal basis for a closed-shell atom
        // might not have any functions for virtuals.
        int nbf = aobasis.nbf(molecule);
        nmo_alpha = std::min(nbf,nmo_alpha);
        nmo_beta = std::min(nbf,nmo_beta);
        if (nalpha>nbf || nbeta>nbf) error("too few basis functions?", nbf);
        nvalpha = nmo_alpha - nalpha;
        nvbeta = nmo_beta - nbeta;

        // Unless overridden by the user use a cell big enough to
        // have exp(-sqrt(2*I)*r) decay to 1e-6 with I=1ev=0.037Eh
        // --> need 50 a.u. either side of the molecule

        if (L == 0.0) {
            L = molecule.bounding_cube() + 50.0;
        }

        lo = molecule.smallest_length_scale();
    }

    void print(World& world) const {
        //time_t t = time((time_t *) 0);
        //char *tmp = ctime(&t);
        //tmp[strlen(tmp)-1] = 0; // lose the trailing newline
        const char* calctype[2] = {"Hartree-Fock","LDA"};
        //madness::print(" date of calculation ", tmp);
        madness::print("             restart ", restart);
        madness::print(" number of processes ", world.size());
        madness::print("   no. of io servers ", nio);
        madness::print("        total charge ", charge);
        madness::print("            smearing ", smear);
        madness::print(" number of electrons ", nalpha, nbeta);
        madness::print("  number of orbitals ", nmo_alpha, nmo_beta);
        madness::print("     spin restricted ", spin_restricted);
        madness::print("  energy convergence ", econv);
        madness::print(" density convergence ", dconv);
        if (conv_only_dens)
            madness::print(" Convergence criterion is only density delta.");
        else
            madness::print(" Convergence criteria are density delta & BSH residual.");
        madness::print("    maximum rotation ", maxrotn);
        if (core_type != "")
            madness::print("           core type ", core_type);
        madness::print(" initial guess basis ", aobasis);
        madness::print(" max krylov subspace ", maxsub);
        madness::print("    calculation type ", calctype[int(lda)]);
        madness::print("     simulation cube ", -L, L);
        madness::print("        plot density ", plotdens);
        madness::print("        plot coulomb ", plotcoul);
        madness::print("        plot orbital ", plotlo, plothi);
        madness::print("        plot npoints ", npt_plot);
        if (plot_cell.size() > 0)
            madness::print("        plot  volume ", plot_cell(0,0), plot_cell(0,1),
                           plot_cell(1,0), plot_cell(1,1), plot_cell(2,0), plot_cell(2,1));
        else
            madness::print("        plot  volume ", "default");

        std::string loctype = "boys";
        if (localize_pm) loctype = "pm";
        if (localize)
            madness::print("  localized orbitals ", loctype);
        else
            madness::print("  canonical orbitals ");
        if (derivatives)
            madness::print("    calc derivatives ");
        if (dipole)
            madness::print("         calc dipole ");
    }
};

struct Calculation {
    Molecule molecule;
    CalculationParameters param;
    AtomicBasisSet aobasis;
    functionT vnuc;
    functionT mask;
    vecfuncT amo, bmo;
    std::vector<int> aset, bset;
    vecfuncT ao;
    std::vector<int> at_to_bf, at_nbf;
    tensorT aocc, bocc;
    tensorT aeps, beps;
    poperatorT coulop;
    std::vector< std::shared_ptr<real_derivative_3d> > gradop;
    double vtol;
    double current_energy;
    Calculation(World & world, const char *filename)
    {
        if(world.rank() == 0) {
            molecule.read_file(filename);
            param.read_file(filename);
            unsigned int n_core = 0;
            if (param.core_type != "") {
                molecule.read_core_file(param.core_type);
                param.aobasis = molecule.guess_file();
                n_core = molecule.n_core_orb_all();
            }
            molecule.orient();
            aobasis.read_file(param.aobasis);
            param.set_molecular_info(molecule, aobasis, n_core);
        }
        world.gop.broadcast_serializable(molecule, 0);
        world.gop.broadcast_serializable(param, 0);
        world.gop.broadcast_serializable(aobasis, 0);
        FunctionDefaults<3>::set_cubic_cell(-param.L, param.L);
        set_protocol(world, 1e-4);
    }

    void set_protocol(World & world, double thresh)
    {
        int k;
        if(thresh >= 1e-2)
            k = 4;
        else if(thresh >= 1e-4)
            k = 6;
        else if(thresh >= 1e-6)
            k = 8;
        else if(thresh >= 1e-8)
            k = 10;
        else
            k = 12;

        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(thresh);
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_initial_level(2);
        FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_apply_randomize(false);
        FunctionDefaults<3>::set_project_randomize(false);
        GaussianConvolution1DCache<double>::map.clear();
        double safety = 0.1;
        vtol = FunctionDefaults<3>::get_thresh() * safety;
        coulop = poperatorT(CoulombOperatorPtr(world, param.lo, thresh));
        gradop = gradient_operator<double,3>(world);
        mask = functionT(factoryT(world).f(mask3).initial_level(4).norefine());
        if(world.rank() == 0){
            print("\nSolving with thresh", thresh, "    k", FunctionDefaults<3>::get_k(), "   conv", std::max(thresh, param.dconv), "\n");
        }
    }

    void save_mos(World& world) {
        archive::ParallelOutputArchive ar(world, "restartdata", param.nio);
        ar & param.spin_restricted;
        ar & (unsigned int)(amo.size());
        ar & aeps & aocc & aset;
        for (unsigned int i=0; i<amo.size(); ++i) ar & amo[i];
        if (!param.spin_restricted) {
            ar & (unsigned int)(bmo.size());
            ar & beps & bocc & bset;
            for (unsigned int i=0; i<bmo.size(); ++i) ar & bmo[i];
        }
    }

    void load_mos(World& world) {
        const double trantol = vtol / std::min(30.0, double(param.nalpha));
        const double thresh = FunctionDefaults<3>::get_thresh();
        const int k = FunctionDefaults<3>::get_k();
        unsigned int nmo;
        bool spinrest;
        amo.clear(); bmo.clear();

        archive::ParallelInputArchive ar(world, "restartdata");

        /*
          File format:

          bool spinrestricted --> if true only alpha orbitals are present

          unsigned int nmo_alpha;
          Tensor<double> aeps;
          Tensor<double> aocc;
          vector<int> aset;
          for i from 0 to nalpha-1:
          .   Function<double,3> amo[i]

          repeat for beta if !spinrestricted

        */

        // LOTS OF LOGIC MISSING HERE TO CHANGE OCCUPATION NO., SET,
        // EPS, SWAP, ... sigh

        ar & spinrest;

        ar & nmo;
        MADNESS_ASSERT(nmo >= unsigned(param.nmo_alpha));
        ar & aeps & aocc & aset;
        amo.resize(nmo);
        for (unsigned int i=0; i<amo.size(); ++i) ar & amo[i];
        unsigned int n_core = molecule.n_core_orb_all();
        if (nmo > unsigned(param.nmo_alpha)) {
            aset = vector<int>(aset.begin()+n_core, aset.begin()+n_core+param.nmo_alpha);
            amo = vecfuncT(amo.begin()+n_core, amo.begin()+n_core+param.nmo_alpha);
            aeps = copy(aeps(Slice(n_core, n_core+param.nmo_alpha-1)));
            aocc = copy(aocc(Slice(n_core, n_core+param.nmo_alpha-1)));
        }

        if (amo[0].k() != k) {
            reconstruct(world,amo);
            for(unsigned int i = 0;i < amo.size();++i) amo[i] = madness::project(amo[i], k, thresh, false);
            world.gop.fence();
        }
        normalize(world, amo);
        amo = transform(world, amo, Q3(matrix_inner(world, amo, amo)), trantol, true);
        truncate(world, amo);
        normalize(world, amo);

        if (!param.spin_restricted) {

            if (spinrest) { // Only alpha spin orbitals were on disk
                MADNESS_ASSERT(param.nmo_alpha >= param.nmo_beta);
                bmo.resize(param.nmo_beta);
                bset.resize(param.nmo_beta);
                beps = copy(aeps(Slice(0,param.nmo_beta-1)));
                bocc = copy(aocc(Slice(0,param.nmo_beta-1)));
                for (int i=0; i<param.nmo_beta; ++i) bmo[i] = copy(amo[i]);
            }
            else {
                ar & nmo;
                ar & beps & bocc & bset;

                bmo.resize(nmo);
                for (unsigned int i=0; i<bmo.size(); ++i) ar & bmo[i];

                if (nmo > unsigned(param.nmo_beta)) {
                    bset = vector<int>(bset.begin()+n_core, bset.begin()+n_core+param.nmo_beta);
                    bmo = vecfuncT(bmo.begin()+n_core, bmo.begin()+n_core+param.nmo_beta);
                    beps = copy(beps(Slice(n_core, n_core+param.nmo_beta-1)));
                    bocc = copy(bocc(Slice(n_core, n_core+param.nmo_beta-1)));
                }

                if (bmo[0].k() != k) {
                    reconstruct(world,bmo);
                    for(unsigned int i = 0;i < bmo.size();++i) bmo[i] = madness::project(bmo[i], k, thresh, false);
                    world.gop.fence();
                }

                normalize(world, bmo);
                bmo = transform(world, bmo, Q3(matrix_inner(world, bmo, bmo)), trantol, true);
                truncate(world, bmo);
                normalize(world, bmo);

            }
        }
    }

    void do_plots(World& world) {
        START_TIMER(world);

        std::vector<long> npt(3,param.npt_plot);

        if (param.plot_cell.size() == 0)
            param.plot_cell = copy(FunctionDefaults<3>::get_cell());

        if (param.plotdens || param.plotcoul) {
            functionT rho;
            rho = make_density(world, aocc, amo);

            if (param.spin_restricted) {
                rho.scale(2.0);
            }
            else {
                functionT rhob = make_density(world, bocc, bmo);
                functionT rho_spin = rho - rhob;
                rho += rhob;
                plotdx(rho_spin, "spin_density.dx", param.plot_cell, npt, true);

            }
            plotdx(rho, "total_density.dx", param.plot_cell, npt, true);
            if (param.plotcoul) {
                functionT vlocl = vnuc + apply(*coulop, rho);
                vlocl.truncate();
                vlocl.reconstruct();
                plotdx(vlocl, "coulomb.dx", param.plot_cell, npt, true);
            }
        }

        for (int i=param.plotlo; i<=param.plothi; ++i) {
            char fname[256];
            if (i < param.nalpha) {
                sprintf(fname, "amo-%5.5d.dx", i);
                plotdx(amo[i], fname, param.plot_cell, npt, true);
            }
            if (!param.spin_restricted && i < param.nbeta) {
                sprintf(fname, "bmo-%5.5d.dx", i);
                plotdx(bmo[i], fname, param.plot_cell, npt, true);
            }
        }
        END_TIMER(world, "plotting");
    }

    void project(World & world)
    {
        reconstruct(world, amo);
        for(unsigned int i = 0;i < amo.size();++i){
            amo[i] = madness::project(amo[i], FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(), false);
        }
        world.gop.fence();
        truncate(world, amo);
        normalize(world, amo);
        if(param.nbeta && !param.spin_restricted){
            reconstruct(world, bmo);
            for(unsigned int i = 0;i < bmo.size();++i){
                bmo[i] = madness::project(bmo[i], FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(), false);
            }
            world.gop.fence();
            truncate(world, bmo);
            normalize(world, bmo);
        }
    }

    void make_nuclear_potential(World & world)
    {
        START_TIMER(world);
        vnuc = factoryT(world).functor(functorT(new MolecularPotentialFunctor(molecule))).thresh(vtol).truncate_on_project();
        vnuc.set_thresh(FunctionDefaults<3>::get_thresh());
        vnuc.reconstruct();
        END_TIMER(world, "Project vnuclear");
        if (param.core_type != "") {
            START_TIMER(world);
            functionT c_pot = factoryT(world).functor(functorT(new MolecularCorePotentialFunctor(molecule))).thresh(vtol).initial_level(4);
            c_pot.set_thresh(FunctionDefaults<3>::get_thresh());
            c_pot.reconstruct();
            END_TIMER(world, "Project Core Pot.");
            vnuc += c_pot;
            vnuc.truncate();
        }
    }

    void project_ao_basis(World & world)
    {
        // Make at_to_bf, at_nbf ... map from atom to first bf on atom, and nbf/atom
        aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

        START_TIMER(world);
        ao = vecfuncT(aobasis.nbf(molecule));
        Level initial_level = 2;
        for(int i = 0;i < aobasis.nbf(molecule);++i){
            functorT aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule, i)));
            ao[i] = factoryT(world).functor(aofunc).initial_level(initial_level).truncate_on_project().nofence();
        }
        world.gop.fence();
        std::vector<double> norms;
        while(1){
            norms = norm2s(world, ao);
            initial_level += 2;
            if(initial_level >= 11)
                throw "project_ao_basis: projection failed?";

            int nredone = 0;
            for(int i = 0;i < aobasis.nbf(molecule);++i){
                if(norms[i] < 0.25){
                    ++nredone;
                    functorT aofunc(new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule, i)));
                    ao[i] = factoryT(world).functor(aofunc).initial_level(6).truncate_on_project().nofence();
                }
            }

            world.gop.fence();
            if(nredone == 0)
                break;

        }
        truncate(world, ao);
        // Don't want to renorm since some d and f functions won't be normalized to unity
        // for(int i = 0;i < aobasis.nbf(molecule);++i){
        //     if(world.rank() == 0 && fabs(norms[i] - 1.0) > 1e-2)
        //         print(i, " bad ao norm?", norms[i]);

        //     norms[i] = 1.0 / norms[i];
        // }
        // scale(world, ao, norms);
        END_TIMER(world, "project ao basis");
    }

    double PM_q(const tensorT & S, const tensorT & C, int i, int j, int lo, int nbf)
    {
        double qij = 0.0;
        if (nbf == 1) { // H atom in STO-3G ... often lots of these!
            qij = C(i,lo)*S(0,0)*C(j,lo);
        }
        else {
            for(int mu = 0;mu < nbf;++mu){
                double Smuj = 0.0;
                for(int nu = 0;nu < nbf;++nu){
                    Smuj += S(mu, nu) * C(j, nu + lo);
                }
                qij += C(i, mu + lo) * Smuj;
            }
        }

        return qij;
    }

    void localize_PM_task_kernel(tensorT & Q, std::vector<tensorT> & Svec, tensorT & C,
                                 const bool & doprint, const std::vector<int> & set,
                                 const double thetamax, tensorT & U, const double thresh)
    {
        long nmo = C.dim(0);
        long nao = C.dim(1);
        long natom = molecule.natom();

        for(long i = 0;i < nmo;++i){
            for(long a = 0;a < natom;++a){
                Q(i, a) = PM_q(Svec[a], C, i, i, at_to_bf[a], at_nbf[a]);
            }
        }

        double tol = 0.1;
        long ndone = 0;
        bool converged = false;
        for(long iter = 0;iter < 100;++iter){
            double sum = 0.0;
            for(long i = 0;i < nmo;++i){
                for(long a = 0;a < natom;++a){
                    double qiia = Q(i, a);
                    sum += qiia * qiia;
                }
            }

            long ndone_iter = 0;
            double maxtheta = 0.0;
            // if(doprint)
            //     printf("iteration %ld sum=%.4f ndone=%ld tol=%.2e\n", iter, sum, ndone, tol);

            for(long i = 0;i < nmo;++i){
                for(long j = 0;j < i;++j){
                    if(set[i] == set[j]){
                        double ovij = 0.0;
                        for(long a = 0;a < natom;++a)
                            ovij += Q(i, a) * Q(j, a);

                        if(fabs(ovij) > tol * tol){
                            double aij = 0.0;
                            double bij = 0.0;
                            for(long a = 0;a < natom;++a){
                                double qiia = Q(i, a);
                                double qija = PM_q(Svec[a], C, i, j, at_to_bf[a], at_nbf[a]);
                                double qjja = Q(j, a);
                                double d = qiia - qjja;
                                aij += qija * qija - 0.25 * d * d;
                                bij += qija * d;
                            }
                            double theta = 0.25 * acos(-aij / sqrt(aij * aij + bij * bij));
                            if(bij > 0.0)
                                theta = -theta;

                            if(theta > thetamax)
                                theta = thetamax;

                            else
                                if(theta < -thetamax)
                                    theta = -thetamax;


                            maxtheta = std::max(fabs(theta), maxtheta);
                            if(fabs(theta) >= tol){
                                ++ndone_iter;
                                double c = cos(theta);
                                double s = sin(theta);
                                drot(nao, &C(i, 0), &C(j, 0), s, c, 1);
                                drot(nmo, &U(i, 0), &U(j, 0), s, c, 1);
                                for(long a = 0;a < natom;++a){
                                    Q(i, a) = PM_q(Svec[a], C, i, i, at_to_bf[a], at_nbf[a]);
                                    Q(j, a) = PM_q(Svec[a], C, j, j, at_to_bf[a], at_nbf[a]);
                                }
                            }

                        }

                    }

                }

            }

            ndone += ndone_iter;
            if(ndone_iter == 0 && tol == thresh){
                if(doprint)
                    print("PM localization converged in", ndone,"steps");

                converged = true;
                break;
            }
            tol = std::max(0.1 * std::min(maxtheta, tol), thresh);
        }

    }

    tensorT localize_PM(World & world, const vecfuncT & mo, const std::vector<int> & set, const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true, const bool doprint = true)
    {
        START_TIMER(world);
        long nmo = mo.size();
        long natom = molecule.natom();

        tensorT S = matrix_inner(world, ao, ao, true);
        std::vector<tensorT> Svec(natom);
        for(long a = 0;a < natom;++a){
            Slice as(at_to_bf[a], at_to_bf[a] + at_nbf[a] - 1);
            Svec[a] = copy(S(as, as));
        }
        S = tensorT();
        tensorT C = matrix_inner(world, mo, ao);
        tensorT U(nmo, nmo);
        tensorT Q(nmo, natom);
        if(world.rank() == 0){
            for(long i = 0;i < nmo;++i)
                U(i, i) = 1.0;

            localize_PM_task_kernel(Q, Svec, C, doprint, set, thetamax, U, thresh);
            U = transpose(U);
        }
        world.gop.broadcast(U.ptr(), U.size(), 0);
        END_TIMER(world, "Pipek-Mezy localize");
        return U;
    }

    void analyze_vectors(World & world, const vecfuncT & mo, const tensorT & occ = tensorT(), const tensorT & energy = tensorT(), const std::vector<int> & set = std::vector<int>())
    {
        tensorT Saomo = matrix_inner(world, ao, mo);
        tensorT Saoao = matrix_inner(world, ao, ao, true);
        int nmo = mo.size();
        tensorT rsq, dip(3, nmo);
        {
            functionT frsq = factoryT(world).f(rsquared).initial_level(4);
            rsq = inner(world, mo, mul_sparse(world, frsq, mo, vtol));
            for(int axis = 0;axis < 3;++axis){
                functionT fdip = factoryT(world).functor(functorT(new DipoleFunctor(axis))).initial_level(4);
                dip(axis, _) = inner(world, mo, mul_sparse(world, fdip, mo, vtol));
                for(int i = 0;i < nmo;++i)
                    rsq(i) -= dip(axis, i) * dip(axis, i);

            }
        }
        if(world.rank() == 0){
            tensorT C;
            gesv(Saoao, Saomo, C);
            C = transpose(C);
            long nmo = mo.size();
            for(long i = 0;i < nmo;++i){
                printf("  MO%4ld : ", i);
                if(set.size())
                    printf("set=%d : ", set[i]);

                if(occ.size())
                    printf("occ=%.2f : ", occ(i));

                if(energy.size())
                    printf("energy=%11.6f : ", energy(i));

                printf("center=(%.2f,%.2f,%.2f) : radius=%.2f\n", dip(0, i), dip(1, i), dip(2, i), sqrt(rsq(i)));
                aobasis.print_anal(molecule, C(i, _));
            }
        }

    }

    inline double DIP(const tensorT & dip, int i, int j, int k, int l)
    {
        return dip(i, j, 0) * dip(k, l, 0) + dip(i, j, 1) * dip(k, l, 1) + dip(i, j, 2) * dip(k, l, 2);
    }

    tensorT localize_boys(World & world, const vecfuncT & mo, const std::vector<int> & set, const double thresh = 1e-9, const double thetamax = 0.5, const bool randomize = true)
    {
        START_TIMER(world);
        const bool doprint = false;
        long nmo = mo.size();
        tensorT dip(nmo, nmo, 3);
        for(int axis = 0;axis < 3;++axis){
            functionT fdip = factoryT(world).functor(functorT(new DipoleFunctor(axis))).initial_level(4);
            dip(_, _, axis) = matrix_inner(world, mo, mul_sparse(world, fdip, mo, vtol), true);
        }
        tensorT U(nmo, nmo);
        if(world.rank() == 0){
            for(long i = 0;i < nmo;++i)
                U(i, i) = 1.0;

            double tol = thetamax;
            long ndone = 0;
            bool converged = false;
            for(long iter = 0;iter < 300;++iter){
                double sum = 0.0;
                for(long i = 0;i < nmo;++i){
                    sum += DIP(dip, i, i, i, i);
                }
                long ndone_iter = 0;
                double maxtheta = 0.0;
                if(doprint)
                    printf("iteration %ld sum=%.4f ndone=%ld tol=%.2e\n", iter, sum, ndone, tol);

                for(long i = 0;i < nmo;++i){
                    for(long j = 0;j < i;++j){
                        if (set[i] == set[j]) {
                            double g = DIP(dip, i, j, j, j) - DIP(dip, i, j, i, i);
                            double h = 4.0 * DIP(dip, i, j, i, j) + 2.0 * DIP(dip, i, i, j, j) - DIP(dip, i, i, i, i) - DIP(dip, j, j, j, j);
                            double sij = DIP(dip, i, j, i, j);
                            bool doit = false;
                            if(h >= 0.0){
                                doit = true;
                                if(doprint)
                                    print("             forcing negative h", i, j, h);

                                h = -1.0;
                            }
                            double theta = -g / h;
                            maxtheta = std::max<double>(std::abs(theta), maxtheta);
                            if(fabs(theta) > thetamax){
                                doit = true;
                                if(doprint)
                                    print("             restricting", i, j);

                                if(g < 0)
                                    theta = -thetamax;

                                else
                                    theta = thetamax * 0.8;

                            }
                            bool randomized = false;
                            if(randomize && iter == 0 && sij > 0.01 && fabs(theta) < 0.01){
                                randomized = true;
                                if(doprint)
                                    print("             randomizing", i, j);

                                theta += (RandomValue<double>() - 0.5);
                            }
                            if(fabs(theta) >= tol || randomized || doit){
                                ++ndone_iter;
                                if(doprint)
                                    print("     rotating", i, j, theta);

                                double c = cos(theta);
                                double s = sin(theta);
                                drot3(nmo, &dip(i, 0, 0), &dip(j, 0, 0), s, c, 1);
                                drot3(nmo, &dip(0, i, 0), &dip(0, j, 0), s, c, nmo);
                                drot(nmo, &U(i, 0), &U(j, 0), s, c, 1);
                            }
                        }
                    }
                }

                ndone += ndone_iter;
                if(ndone_iter == 0 && tol == thresh){
                    if(doprint)
                        print("Boys localization converged in", ndone,"steps");

                    converged = true;
                    break;
                }
                tol = std::max(0.1 * maxtheta, thresh);
            }

            if(!converged){
                print("warning: boys localization did not fully converge: ", ndone);
            }
            U = transpose(U);
        }

        world.gop.broadcast(U.ptr(), U.size(), 0);
        END_TIMER(world, "Boys localize");
        return U;
    }

    tensorT kinetic_energy_matrix(World & world, const vecfuncT & v)
    {
        reconstruct(world, v);
        int n = v.size();
        tensorT r(n, n);
        for(int axis = 0;axis < 3;++axis){
            vecfuncT dv = apply(world, *(gradop[axis]), v);
            r += matrix_inner(world, dv, dv, true);
            dv.clear();
        }
        return r.scale(0.5);
    }

    struct GuessDensity : public FunctionFunctorInterface<double,3>
    {
        const Molecule & molecule;
        const AtomicBasisSet & aobasis;
        double operator ()(const coordT & x) const
        {
            return aobasis.eval_guess_density(molecule, x[0], x[1], x[2]);
        }

        GuessDensity(const Molecule & molecule, const AtomicBasisSet & aobasis)
            :molecule(molecule), aobasis(aobasis)
        {
        }

    };

    vecfuncT core_projection(World & world, const vecfuncT & psi, const bool include_Bc = true)
    {
        int npsi = psi.size();
        if (npsi == 0) return psi;
        int natom = molecule.natom();
        vecfuncT proj = zero_functions<double,3>(world, npsi);
        tensorT overlap_sum(static_cast<long>(npsi));

        for (int i=0; i<natom; ++i) {
            Atom at = molecule.get_atom(i);
            unsigned int atn = at.atomic_number;
            unsigned int nshell = molecule.n_core_orb(atn);
            if (nshell == 0) continue;
            for (unsigned int c=0; c<nshell; ++c) {
                unsigned int l = molecule.get_core_l(atn, c);
                int max_m = (l+1)*(l+2)/2;
                nshell -= max_m - 1;
                for (int m=0; m<max_m; ++m) {
                    functionT core = factoryT(world).functor(functorT(new CoreOrbitalFunctor(molecule, i, c, m)));
                    tensorT overlap = inner(world, core, psi);
                    overlap_sum += overlap;
                    for (int j=0; j<npsi; ++j) {
                        if (include_Bc) overlap[j] *= molecule.get_core_bc(atn, c);
                        proj[j] += core.scale(overlap[j]);
                    }
                }
            }
            world.gop.fence();
        }
        if (world.rank() == 0) print("sum_k <core_k|psi_i>:", overlap_sum);
        return proj;
    }

    double core_projector_derivative(World & world, const vecfuncT & mo, const tensorT & occ, int atom, int axis)
    {
        vecfuncT cores, dcores;
        std::vector<double> bc;
        unsigned int atn = molecule.get_atom(atom).atomic_number;
        unsigned int ncore = molecule.n_core_orb(atn);

        // projecting core & d/dx core
        for (unsigned int c=0; c<ncore; ++c) {
            unsigned int l = molecule.get_core_l(atn, c);
            int max_m = (l+1)*(l+2)/2;
            for (int m=0; m<max_m; ++m) {
                functorT func = functorT(new CoreOrbitalFunctor(molecule, atom, c, m));
                cores.push_back(functionT(factoryT(world).functor(func).truncate_on_project()));
                func = functorT(new CoreOrbitalDerivativeFunctor(molecule, atom, axis, c, m));
                dcores.push_back(functionT(factoryT(world).functor(func).truncate_on_project()));
                bc.push_back(molecule.get_core_bc(atn, c));
            }
        }

        // calc \sum_i occ_i <psi_i|(\sum_c Bc d/dx |core><core|)|psi_i>
        double r = 0.0;
        for (unsigned int c=0; c<cores.size(); ++c) {
            double rcore= 0.0;
            tensorT rcores = inner(world, cores[c], mo);
            tensorT rdcores = inner(world, dcores[c], mo);
            for (unsigned int i=0; i<mo.size(); ++i) {
                rcore += rdcores[i] * rcores[i] * occ[i];
            }
            r += 2.0 * bc[c] * rcore;
        }

        return r;
    }

    void initial_guess(World & world)
    {
        START_TIMER(world);
        if (param.restart) {
            load_mos(world);
        }
        else {
            // Use the initial density and potential to generate a better process map
            functionT rho = factoryT(world).functor(functorT(new GuessDensity(molecule, aobasis))).truncate_on_project();
            END_TIMER(world, "guess density");
            double nel = rho.trace();
            if(world.rank() == 0)
                print("guess dens trace", nel);

            if(world.size() > 1) {
                START_TIMER(world);
                LoadBalanceDeux<3> lb(world);
                lb.add_tree(vnuc, lbcost<double,3>(1.0, 0.0), false);
                lb.add_tree(rho, lbcost<double,3>(1.0, 1.0), true);

                FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));

//                 FunctionDefaults<3>::set_pmap(lb.load_balance(6.0));
//                 rho = copy(rho, FunctionDefaults<3>::get_pmap(), false);
//                 vnuc = copy(vnuc, FunctionDefaults<3>::get_pmap(), false);
//                 world.gop.fence();
                END_TIMER(world, "guess loadbal");
            }

            // Diag approximate fock matrix to get initial mos
            functionT vlocal;
            if(param.nalpha + param.nbeta > 1){
                START_TIMER(world);
                vlocal = vnuc + apply(*coulop, rho);
                END_TIMER(world, "guess Coulomb potn");
                bool save = param.spin_restricted;
                param.spin_restricted = true;
                rho.scale(0.5);
                vlocal = vlocal + make_lda_potential(world, rho, rho, functionT(), functionT());
                vlocal.truncate();
                param.spin_restricted = save;
            } else {
                vlocal = vnuc;
            }
            rho.clear();
            vlocal.reconstruct();
            if(world.size() > 1){
                LoadBalanceDeux<3> lb(world);
                lb.add_tree(vnuc, lbcost<double,3>(1.0, 1.0), false);
                for(unsigned int i = 0;i < ao.size();++i){
                    lb.add_tree(ao[i], lbcost<double,3>(1.0, 1.0), false);
                }

                FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));

//                 FunctionDefaults<3>::set_pmap(lb.load_balance(6.0));
//                 vnuc = copy(vnuc, FunctionDefaults<3>::get_pmap(), false);
//                 vlocal = copy(vlocal, FunctionDefaults<3>::get_pmap(), false);
//                 for(unsigned int i = 0;i < ao.size();++i){
//                     ao[i] = copy(ao[i], FunctionDefaults<3>::get_pmap(), false);
//                 }
//                 world.gop.fence();
            }

            tensorT overlap = matrix_inner(world, ao, ao, true);
            tensorT kinetic = kinetic_energy_matrix(world, ao);
            reconstruct(world, ao);
            vlocal.reconstruct();
            vecfuncT vpsi = mul_sparse(world, vlocal, ao, vtol);
            compress(world, vpsi);
            truncate(world, vpsi);
            compress(world, ao);
            tensorT potential = matrix_inner(world, vpsi, ao, true);
            vpsi.clear();
            tensorT fock = kinetic + potential;
            fock = 0.5 * (fock + transpose(fock));
            tensorT c, e;
            sygv(fock, overlap, 1, c, e);
            world.gop.broadcast(c.ptr(), c.size(), 0);
            world.gop.broadcast(e.ptr(), e.size(), 0);
            if(world.rank() == 0){
                print("initial eigenvalues");
                print(e);
            }
            compress(world, ao);

            unsigned int ncore = 0;
            if (param.core_type != "") {
                ncore = molecule.n_core_orb_all();
            }
            amo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_alpha - 1)), 0.0, true);
            truncate(world, amo);
            normalize(world, amo);
            aeps = e(Slice(ncore, ncore + param.nmo_alpha - 1));

            aocc = tensorT(param.nmo_alpha);
            for(int i = 0;i < param.nalpha;++i)
                aocc[i] = 1.0;

            aset = std::vector<int>(param.nmo_alpha,0);
            //if (param.localize_pm) {
                aset[0] = 0;
                if(world.rank() == 0)
                    std::cout << "alpha set " << 0 << " " << 0 << "-";

                for(int i = 1;i < param.nmo_alpha;++i) {
                    aset[i] = aset[i - 1];
                    if(aeps[i] - aeps[i - 1] > 1.5 || aocc[i] != 1.0){
                        ++(aset[i]);
                        if(world.rank() == 0){
                            std::cout << i - 1 << std::endl;
                            std::cout << "alpha set " << aset[i] << " " << i << "-";
                        }
                    }
                }
                if(world.rank() == 0)
                    std::cout << param.nmo_alpha - 1 << std::endl;
            //}

            if(param.nbeta && !param.spin_restricted){
                bmo = transform(world, ao, c(_, Slice(ncore, ncore + param.nmo_beta - 1)), 0.0, true);
                truncate(world, bmo);
                normalize(world, bmo);
                beps = e(Slice(ncore, ncore + param.nmo_beta - 1));
                bocc = tensorT(param.nmo_beta);
                for(int i = 0;i < param.nbeta;++i)
                    bocc[i] = 1.0;

                bset = std::vector<int>(param.nmo_beta,0);
                //if (param.localize_pm) {
                    bset[0] = 0;
                    if(world.rank() == 0)
                        std::cout << " beta set " << 0 << " " << 0 << "-";

                    for(int i = 1;i < param.nmo_beta;++i) {
                        bset[i] = bset[i - 1];
                        if(beps[i] - beps[i - 1] > 1.5 || bocc[i] != 1.0){
                            ++(bset[i]);
                            if(world.rank() == 0){
                                std::cout << i - 1 << std::endl;
                                std::cout << " beta set " << bset[i] << " " << i << "-";
                            }
                        }
                    }
                    if(world.rank() == 0)
                        std::cout << param.nmo_beta - 1 << std::endl;
                //}

            }
        }
    }

    void initial_load_bal(World & world)
    {
        LoadBalanceDeux<3> lb(world);
        lb.add_tree(vnuc, lbcost<double,3>(1.0, 0.0));

        FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));

//         FunctionDefaults<3>::set_pmap(lb.load_balance(6.0));
//         world.gop.fence();
    }

    functionT make_density(World & world, const tensorT & occ, const vecfuncT & v)
    {
        vecfuncT vsq = square(world, v);
        compress(world, vsq);
        functionT rho = factoryT(world);
        rho.compress();
        for(unsigned int i = 0;i < vsq.size();++i){
            if(occ[i])
                rho.gaxpy(1.0, vsq[i], occ[i], false);

        }
        world.gop.fence();
        vsq.clear();
        return rho;
    }

    std::vector<poperatorT> make_bsh_operators(World & world, const tensorT & evals)
    {
        int nmo = evals.dim(0);
        std::vector<poperatorT> ops(nmo);
        double tol = FunctionDefaults<3>::get_thresh();
        for(int i = 0;i < nmo;++i){
            double eps = evals(i);
            if(eps > 0){
                if(world.rank() == 0){
                    print("bsh: warning: positive eigenvalue", i, eps);
                }
                eps = -0.1;
            }

            ops[i] = poperatorT(BSHOperatorPtr3D(world, sqrt(-2.0 * eps),  param.lo, tol));
        }

        return ops;
    }

    vecfuncT apply_hf_exchange(World & world, const tensorT & occ, const vecfuncT & psi, const vecfuncT & f)
    {
        const bool same = (&psi == &f);
        int nocc = psi.size();
        int nf = f.size();
        double tol = FunctionDefaults<3>::get_thresh(); /// Important this is consistent with Coulomb
        vecfuncT Kf = zero_functions<double,3>(world, nf);
        compress(world, Kf);
        reconstruct(world, psi);
        norm_tree(world, psi);
        if (!same) {
            reconstruct(world, f);
            norm_tree(world, f);
        }

//         // Smaller memory algorithm ... possible 2x saving using i-j sym
//         for(int i=0; i<nocc; ++i){
//             if(occ[i] > 0.0){
//                 vecfuncT psif = mul_sparse(world, psi[i], f, tol); /// was vtol
//                 truncate(world, psif);
//                 psif = apply(world, *coulop, psif);
//                 truncate(world, psif);
//                 psif = mul_sparse(world, psi[i], psif, tol); /// was vtol
//                 gaxpy(world, 1.0, Kf, occ[i], psif);
//             }
//         }

        // Larger memory algorithm ... use i-j sym if psi==f
        vecfuncT psif;
        for (int i=0; i<nocc; ++i) {
            int jtop = nf;
            if (same) jtop = i+1;
            for (int j=0; j<jtop; ++j) {
                psif.push_back(mul_sparse(psi[i], f[j], tol, false));
            }
        }

        world.gop.fence();
        truncate(world, psif);
        psif = apply(world, *coulop, psif);
        truncate(world, psif, tol);
        reconstruct(world, psif);
        norm_tree(world, psif);
        vecfuncT psipsif = zero_functions<double,3>(world, nf*nocc);
        int ij = 0;
        for (int i=0; i<nocc; ++i) {
            int jtop = nf;
            if (same) jtop = i+1;
            for (int j=0; j<jtop; ++j,++ij) {
                psipsif[i*nf+j] = mul_sparse(psif[ij],psi[i],false);
                if (same && i!=j) {
                    psipsif[j*nf+i] = mul_sparse(psif[ij],psi[j],false);
                }
            }
        }
        world.gop.fence();
        psif.clear();
        world.gop.fence();
        compress(world, psipsif);
        for (int i=0; i<nocc; ++i) {
            for (int j=0; j<nf; ++j) {
                Kf[j].gaxpy(1.0,psipsif[i*nf+j],occ[i],false);
            }
        }
        world.gop.fence();
        psipsif.clear();
        world.gop.fence();

        truncate(world, Kf, tol);
        return Kf;
    }

    static double munge(double r)
    {
        if(r < 1e-12)
            r = 1e-12;

        return r;
    }

    template <typename T>
    static void ldaop(const Key<3> & key, Tensor<T>& t)
    {
        UNARY_OPTIMIZED_ITERATOR(T, t, double r=munge(2.0* *_p0); double q; double dq1; double dq2;x_rks_s__(&r, &q, &dq1);c_rks_vwn5__(&r, &q, &dq2); *_p0 = dq1+dq2);
    }

    template <typename T>
    static void ldaeop(const Key<3> & key, Tensor<T>& t)
    {
        UNARY_OPTIMIZED_ITERATOR(T, t, double r=munge(2.0* *_p0); double q1; double q2; double dq;x_rks_s__(&r, &q1, &dq);c_rks_vwn5__(&r, &q2, &dq); *_p0 = q1+q2);
    }

    functionT make_lda_potential(World & world, const functionT & arho, const functionT & brho, const functionT & adelrhosq, const functionT & bdelrhosq)
    {
        MADNESS_ASSERT(param.spin_restricted);
        functionT vlda = copy(arho);
        vlda.reconstruct();
        vlda.unaryop(&ldaop<double>);
        return vlda;
    }

    double make_lda_energy(World & world, const functionT & arho, const functionT & brho, const functionT & adelrhosq, const functionT & bdelrhosq)
    {
        MADNESS_ASSERT(param.spin_restricted);
        functionT vlda = copy(arho);
        vlda.reconstruct();
        vlda.unaryop(&ldaeop<double>);
        return vlda.trace();
    }

    vecfuncT apply_potential(World & world, const tensorT & occ, const vecfuncT & amo, const functionT & arho, const functionT & brho, const functionT & adelrhosq, const functionT & bdelrhosq, const functionT & vlocal, double & exc)
    {
        functionT vloc = vlocal;
        if(param.lda){
            START_TIMER(world);
            exc = make_lda_energy(world, arho, brho, adelrhosq, bdelrhosq);
            vloc = vloc + make_lda_potential(world, arho, brho, adelrhosq, bdelrhosq);
            END_TIMER(world, "LDA potential");
        }
        START_TIMER(world);
        vecfuncT Vpsi = mul_sparse(world, vloc, amo, vtol);
        END_TIMER(world, "V*psi");
        if(!param.lda){
            START_TIMER(world);
            vecfuncT Kamo = apply_hf_exchange(world, occ, amo, amo);
            tensorT excv = inner(world, Kamo, amo);
            exc = 0.0;
            for(unsigned long i = 0;i < amo.size();++i){
                exc -= 0.5 * excv[i] * occ[i];
            }
            gaxpy(world, 1.0, Vpsi, -1.0, Kamo);
            Kamo.clear();
            END_TIMER(world, "HF exchange");
        }

        if (param.core_type.substr(0,3) == "mcp") {
            START_TIMER(world);
            gaxpy(world, 1.0, Vpsi, 1.0, core_projection(world, amo));
            END_TIMER(world, "MCP Core Projector");
        }

        START_TIMER(world);
        truncate(world, Vpsi);
        END_TIMER(world, "Truncate Vpsi");
        world.gop.fence();
        return Vpsi;
    }

    tensorT derivatives(World & world)
    {
        START_TIMER(world);

        functionT rho = make_density(world, aocc, amo);
        functionT brho = rho;
        if (!param.spin_restricted) brho = make_density(world, bocc, bmo);
        rho.gaxpy(1.0, brho, 1.0);

        vecfuncT dv(molecule.natom() * 3);
        vecfuncT du = zero_functions<double,3>(world, molecule.natom() * 3);
        tensorT rc(molecule.natom() * 3);
        for(int atom = 0;atom < molecule.natom();++atom){
            for(int axis = 0;axis < 3;++axis){
                functorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
                dv[atom * 3 + axis] = functionT(factoryT(world).functor(func).nofence().truncate_on_project());
                if (param.core_type != "" && molecule.is_potential_defined_atom(atom)) {
                    // core potential contribution
                    func = functorT(new CorePotentialDerivativeFunctor(molecule, atom, axis));
                    du[atom * 3 + axis] = functionT(factoryT(world).functor(func).truncate_on_project());

                    // core projector contribution
                    rc[atom * 3 + axis] = core_projector_derivative(world, amo, aocc, atom, axis);
                    if (!param.spin_restricted) {
                        if (param.nbeta) rc[atom * 3 + axis] += core_projector_derivative(world, bmo, bocc, atom, axis);
                    }
                    else {
                        rc[atom * 3 + axis] *= 2 * 2;
                            // because of 2 electrons in each valence orbital bra+ket
                    }
                }
            }
        }

        world.gop.fence();
        tensorT r = inner(world, rho, dv);
        world.gop.fence();
        tensorT ru = inner(world, rho, du);
        dv.clear();
        du.clear();
        world.gop.fence();
        tensorT ra(r.size());
        for(int atom = 0;atom < molecule.natom();++atom){
            for(int axis = 0;axis < 3;++axis){
                ra[atom * 3 + axis] = molecule.nuclear_repulsion_derivative(atom, axis);
            }
        }
        //if (world.rank() == 0) print("derivatives:\n", r, ru, rc, ra);
        r +=  ra + ru + rc;
        END_TIMER(world,"derivatives");

        if (world.rank() == 0) {
            print("\n Derivatives (a.u.)\n -----------\n");
            print("  atom        x            y            z          dE/dx        dE/dy        dE/dz");
            print(" ------ ------------ ------------ ------------ ------------ ------------ ------------");
            for (int i=0; i<molecule.natom(); ++i) {
                const Atom& atom = molecule.get_atom(i);
                printf(" %5d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f\n",
                       i, atom.x, atom.y, atom.z,
                       r[i*3+0], r[i*3+1], r[i*3+2]);
            }
        }
        return r;
    }

    tensorT dipole(World & world)
    {
        START_TIMER(world);
        tensorT mu(3);
        for (unsigned int axis=0; axis<3; ++axis) {
            std::vector<int> x(3, 0);
            x[axis] = true;
            functionT dipolefunc = factoryT(world).functor(functorT(new MomentFunctor(x)));
            functionT rho = make_density(world, aocc, amo);
            if (!param.spin_restricted) {
                if (param.nbeta) rho += make_density(world, bocc, bmo);
            }
            else {
                rho.scale(2.0);
            }
            mu[axis] = -dipolefunc.inner(rho);
            mu[axis] += molecule.nuclear_dipole(axis);
        }

        if (world.rank() == 0) {
            print("\n Dipole Moment (a.u.)\n -----------\n");
            print("     x: ", mu[0]);
            print("     y: ", mu[1]);
            print("     z: ", mu[2]);
            print(" Total Dipole Moment: ", mu.normf());
        }
        END_TIMER(world, "dipole");

        return mu;
    }

    void vector_stats(const std::vector<double> & v, double & rms, double & maxabsval)
    {
        rms = 0.0;
        maxabsval = v[0];
        for(unsigned int i = 0;i < v.size();++i){
            rms += v[i] * v[i];
            maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
        }
        rms = sqrt(rms / v.size());
    }

    vecfuncT compute_residual(World & world, tensorT & occ, tensorT & fock, const vecfuncT & psi, vecfuncT & Vpsi, double & err)
    {
        double trantol = vtol / std::min(30.0, double(psi.size()));
        int nmo = psi.size();

        tensorT eps(nmo);
        for(int i = 0;i < nmo;++i){
            eps(i) = std::min(-0.05, fock(i, i));
            fock(i, i) -= eps(i);
        }
        vecfuncT fpsi = transform(world, psi, fock, trantol, true);

        for(int i = 0;i < nmo;++i){ // Undo the damage
            fock(i, i) += eps(i);
        }

        gaxpy(world, 1.0, Vpsi, -1.0, fpsi);
        fpsi.clear();
        std::vector<double> fac(nmo, -2.0);
        scale(world, Vpsi, fac);
        std::vector<poperatorT> ops = make_bsh_operators(world, eps);
        set_thresh(world, Vpsi, FunctionDefaults<3>::get_thresh());
        if(world.rank() == 0)
            std::cout << "entering apply\n";

        START_TIMER(world);
        vecfuncT new_psi = apply(world, ops, Vpsi);
        END_TIMER(world, "Apply BSH");
        ops.clear();
        Vpsi.clear();
        world.gop.fence();
        START_TIMER(world);
        truncate(world, new_psi);
        END_TIMER(world, "Truncate new psi");
        vecfuncT r = sub(world, psi, new_psi);
        std::vector<double> rnorm = norm2s(world, r);
        if (world.rank() == 0) print("residuals", rnorm);
        double rms, maxval;
        vector_stats(rnorm, rms, maxval);
        err = maxval;
        if(world.rank() == 0)
            print("BSH residual: rms", rms, "   max", maxval);

        return r;
    }

    tensorT make_fock_matrix(World & world, const vecfuncT & psi, const vecfuncT & Vpsi, const tensorT & occ, double & ekinetic)
    {
        START_TIMER(world);
        tensorT pe = matrix_inner(world, Vpsi, psi, true);
        END_TIMER(world, "PE matrix");
        START_TIMER(world);
        tensorT ke = kinetic_energy_matrix(world, psi);
        END_TIMER(world, "KE matrix");
        int nocc = occ.size();
        ekinetic = 0.0;
        for(int i = 0;i < nocc;++i){
            ekinetic += occ[i] * ke(i, i);
        }
        ke += pe;
        pe = tensorT();
        ke.gaxpy(0.5, transpose(ke), 0.5);
        return ke;
    }

    tensorT matrix_exponential(const tensorT& A) {
        const double tol = 1e-13;
        MADNESS_ASSERT(A.dim((0) == A.dim(1)));

        // Scale A by a power of 2 until it is "small"
        double anorm = A.normf();
        int n = 0;
        double scale = 1.0;
        while (anorm*scale > 0.1) {
            ++n;
            scale *= 0.5;
        }
        tensorT B = scale*A;    // B = A*2^-n

        // Compute exp(B) using Taylor series
        tensorT expB = tensorT(2, B.dims());
        for (int i=0; i<expB.dim(0); ++i) expB(i,i) = 1.0;

        int k = 1;
        tensorT term = B;
        while (term.normf() > tol) {
            expB += term;
            term = inner(term,B);
            ++k;
            term.scale(1.0/k);
        }

        // Repeatedly square to recover exp(A)
        while (n--) {
            expB = inner(expB,expB);
        }

        return expB;
    }

    tensorT diag_fock_matrix(World & world, tensorT& fock, vecfuncT & psi, vecfuncT & Vpsi, tensorT & evals, double thresh)
    {
        long nmo = psi.size();
        tensorT overlap = matrix_inner(world, psi, psi, true);
        START_TIMER(world);
        tensorT U;

        sygv(fock, overlap, 1, U, evals);
        END_TIMER(world, "Diagonalization");

        // Fix phases.
        long j;
        for (long i=0; i<nmo; ++i) {
            U(_,i).absmax(&j);
            if (U(j,i) < 0) U(_,i).scale(-1.0);
        }

        // Within blocks with the same occupation number attempt to
        // keep orbitals in the same order (to avoid confusing the
        // non-linear solver).  Have to run the reordering multiple
        // times to handle multiple degeneracies.
        for (int pass=0; pass<5; ++pass) {
            for (long i=0; i<nmo; ++i) {
                U(_,i).absmax(&j);
                if (i != j) {
                    tensorT tmp = copy(U(_,i));
                    U(_,i) = U(_,j);
                    U(_,j) = tmp;
                    std::swap(evals[i],evals[j]);
                }
            }
        }

        // Rotations between effectively degenerate states confound
        // the non-linear equation solver ... undo these rotations
        long ilo = 0; // first element of cluster
        while (ilo < nmo-1) {
            long ihi = ilo;
            while (fabs(evals[ilo]-evals[ihi+1]) < thresh*10.0*std::max(fabs(evals[ilo]),1.0)) {
                ++ihi;
                if (ihi == nmo-1) break;
            }
            long nclus = ihi - ilo + 1;
            if (nclus > 1) {
                //print("   found cluster", ilo, ihi);
                tensorT q = copy(U(Slice(ilo,ihi),Slice(ilo,ihi)));
                //print(q);
                // Special code just for nclus=2
                // double c = 0.5*(q(0,0) + q(1,1));
                // double s = 0.5*(q(0,1) - q(1,0));
                // double r = sqrt(c*c + s*s);
                // c /= r;
                // s /= r;
                // q(0,0) = q(1,1) = c;
                // q(0,1) = -s;
                // q(1,0) = s;

                // Iteratively construct unitary rotation by
                // exponentiating the antisymmetric part of the matrix
                // ... is quadratically convergent so just do 3
                // iterations
                tensorT rot = matrix_exponential(-0.5*(q - transpose(q)));
                q = inner(q,rot);
                tensorT rot2 = matrix_exponential(-0.5*(q - transpose(q)));
                q = inner(q,rot2);
                tensorT rot3 = matrix_exponential(-0.5*(q - transpose(q)));
                q = inner(rot,inner(rot2,rot3));
                U(_,Slice(ilo,ihi)) = inner(U(_,Slice(ilo,ihi)),q);
            }
            ilo = ihi+1;
        }

        // print("Fock");
        // print(fock);
        // print("Evec");
        // print(U);;
        // print("Eval");
        // print(evals);

        world.gop.broadcast(U.ptr(), U.size(), 0);
        world.gop.broadcast(evals.ptr(), evals.size(), 0);

        fock = 0;
        for (unsigned int i=0; i<psi.size(); ++i) fock(i,i) = evals(i);

        Vpsi = transform(world, Vpsi, U, vtol / std::min(30.0, double(psi.size())), false);
        psi = transform(world, psi, U, FunctionDefaults<3>::get_thresh() / std::min(30.0, double(psi.size())), true);
        truncate(world, Vpsi, vtol, false);
        truncate(world, psi);
        normalize(world, psi);

        return U;
    }

    void loadbal(World & world, functionT & arho, functionT & brho, functionT & arho_old, functionT & brho_old, subspaceT & subspace)
    {
        if(world.size() == 1)
            return;

        LoadBalanceDeux<3> lb(world);
        lb.add_tree(vnuc, lbcost<double,3>(1.0, 0.0), false);
        lb.add_tree(arho, lbcost<double,3>(1.0, 1.0), false);
        for(unsigned int i = 0;i < amo.size();++i){
            lb.add_tree(amo[i], lbcost<double,3>(1.0, 1.0), false);
        }
        if(param.nbeta && !param.spin_restricted){
            lb.add_tree(brho, lbcost<double,3>(1.0, 1.0), false);
            for(unsigned int i = 0;i < bmo.size();++i){
                lb.add_tree(bmo[i], lbcost<double,3>(1.0, 1.0), false);
            }
        }

        FunctionDefaults<3>::redistribute(world, lb.load_balance(6.0));

//         FunctionDefaults<3>::set_pmap(lb.load_balance(6.0));
//         vnuc = copy(vnuc, FunctionDefaults<3>::get_pmap(), false);
//         arho = copy(arho, FunctionDefaults<3>::get_pmap(), false);
//         if(arho_old.is_initialized())
//             arho_old = copy(arho_old, FunctionDefaults<3>::get_pmap(), false);

//         for(unsigned int i = 0;i < ao.size();++i){
//             ao[i] = copy(ao[i], FunctionDefaults<3>::get_pmap(), false);
//         }
//         for(unsigned int i = 0;i < amo.size();++i){
//             amo[i] = copy(amo[i], FunctionDefaults<3>::get_pmap(), false);
//         }
//         if(param.nbeta && !param.spin_restricted){
//             brho = copy(brho, FunctionDefaults<3>::get_pmap(), false);
//             if(brho_old.is_initialized())
//                 brho_old = copy(brho_old, FunctionDefaults<3>::get_pmap(), false);

//             for(unsigned int i = 0;i < bmo.size();++i){
//                 bmo[i] = copy(bmo[i], FunctionDefaults<3>::get_pmap(), false);
//             }
//         }

//         for(unsigned int i = 0;i < subspace.size();++i){
//             vecfuncT & v = subspace[i].first;
//             vecfuncT & r = subspace[i].second;
//             for(unsigned int j = 0;j < v.size();++j){
//                 v[j] = copy(v[j], FunctionDefaults<3>::get_pmap(), false);
//                 r[j] = copy(r[j], FunctionDefaults<3>::get_pmap(), false);
//             }
//         }

//         world.gop.fence();
    }

    void rotate_subspace(World& world, const tensorT& U, subspaceT& subspace, int lo, int nfunc, double trantol) {
        for (unsigned int iter=0; iter<subspace.size(); ++iter) {
            vecfuncT& v = subspace[iter].first;
            vecfuncT& r = subspace[iter].second;
            transform(world, vecfuncT(&v[lo],&v[lo+nfunc]), U, trantol, false);
            transform(world, vecfuncT(&r[lo],&r[lo+nfunc]), U, trantol, true);
        }
    }

    void update_subspace(World & world,
                         vecfuncT & Vpsia, vecfuncT & Vpsib,
                         tensorT & focka, tensorT & fockb,
                         subspaceT & subspace, tensorT & Q,
                         double & bsh_residual, double & update_residual)
    {
        double aerr = 0.0, berr = 0.0;
        vecfuncT vm = amo;
        vecfuncT rm = compute_residual(world, aocc, focka, amo, Vpsia, aerr);
        if(param.nbeta && !param.spin_restricted){
            vecfuncT br = compute_residual(world, bocc, fockb, bmo, Vpsib, berr);
            vm.insert(vm.end(), bmo.begin(), bmo.end());
            rm.insert(rm.end(), br.begin(), br.end());
        }
        bsh_residual = std::max(aerr, berr);
        world.gop.broadcast(bsh_residual, 0);
        compress(world, vm, false);
        compress(world, rm, false);
        world.gop.fence();
        subspace.push_back(pairvecfuncT(vm, rm));
        int m = subspace.size();
        tensorT ms(m);
        tensorT sm(m);
        for(int s = 0;s < m;++s){
            const vecfuncT & vs = subspace[s].first;
            const vecfuncT & rs = subspace[s].second;
            for(unsigned int i = 0;i < vm.size();++i){
                ms[s] += vm[i].inner_local(rs[i]);
                sm[s] += vs[i].inner_local(rm[i]);
            }
        }

        world.gop.sum(ms.ptr(), m);
        world.gop.sum(sm.ptr(), m);
        tensorT newQ(m, m);
        if(m > 1)
            newQ(Slice(0, -2), Slice(0, -2)) = Q;

        newQ(m - 1, _) = ms;
        newQ(_, m - 1) = sm;
        Q = newQ;
        //if (world.rank() == 0) { print("kain Q"); print(Q); }
        tensorT c;
        if(world.rank() == 0){
            double rcond = 1e-12;
            while(1){
                c = KAIN(Q, rcond);
                //if (world.rank() == 0) print("kain c:", c);
                if(std::abs(c[m - 1]) < 3.0){
                    break;
                } else  if(rcond < 0.01){
                    print("Increasing subspace singular value threshold ", c[m - 1], rcond);
                    rcond *= 100;
                } else {
                    print("Forcing full step due to subspace malfunction");
                    c = 0.0;
                    c[m - 1] = 1.0;
                    break;
                }
            }
        }

        world.gop.broadcast_serializable(c, 0);
        if(world.rank() == 0){
            print("Subspace solution", c);
        }
        START_TIMER(world);
        vecfuncT amo_new = zero_functions<double,3>(world, amo.size());
        vecfuncT bmo_new = zero_functions<double,3>(world, bmo.size());
        compress(world, amo_new, false);
        compress(world, bmo_new, false);
        world.gop.fence();
        for(unsigned int m = 0;m < subspace.size();++m){
            const vecfuncT & vm = subspace[m].first;
            const vecfuncT & rm = subspace[m].second;
            const vecfuncT vma(vm.begin(), vm.begin() + amo.size());
            const vecfuncT rma(rm.begin(), rm.begin() + amo.size());
            const vecfuncT vmb(vm.end() - bmo.size(), vm.end());
            const vecfuncT rmb(rm.end() - bmo.size(), rm.end());
            gaxpy(world, 1.0, amo_new, c(m), vma, false);
            gaxpy(world, 1.0, amo_new, -c(m), rma, false);
            gaxpy(world, 1.0, bmo_new, c(m), vmb, false);
            gaxpy(world, 1.0, bmo_new, -c(m), rmb, false);
        }
        world.gop.fence();
        END_TIMER(world, "Subspace transform");
        if(param.maxsub <= 1){
            subspace.clear();
        } else  if(subspace.size() == param.maxsub){
            subspace.erase(subspace.begin());
            Q = Q(Slice(1, -1), Slice(1, -1));
        }

        std::vector<double> anorm = norm2s(world, sub(world, amo, amo_new));
        std::vector<double> bnorm = norm2s(world, sub(world, bmo, bmo_new));
        int nres = 0;
        for(unsigned int i = 0;i < amo.size();++i){
            if(anorm[i] > param.maxrotn){
                double s = param.maxrotn / anorm[i];
                ++nres;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for alpha orbitals:");

                    printf(" %d", i);
                }
                amo_new[i].gaxpy(s, amo[i], 1.0 - s, false);
            }

        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");

        nres = 0;
        for(unsigned int i = 0;i < bmo.size();++i){
            if(bnorm[i] > param.maxrotn){
                double s = param.maxrotn / bnorm[i];
                ++nres;
                if(world.rank() == 0){
                    if(nres == 1)
                        printf("  restricting step for  beta orbitals:");

                    printf(" %d", i);
                }
                bmo_new[i].gaxpy(s, bmo[i], 1.0 - s, false);
            }

        }
        if(nres > 0 && world.rank() == 0)
            printf("\n");

        world.gop.fence();
        double rms, maxval;
        vector_stats(anorm, rms, maxval);
        if(world.rank() == 0)
            print("Norm of vector changes alpha: rms", rms, "   max", maxval);

        update_residual = maxval;
        if(bnorm.size()){
            vector_stats(bnorm, rms, maxval);
            if(world.rank() == 0)
                print("Norm of vector changes  beta: rms", rms, "   max", maxval);

            update_residual = std::max(update_residual, maxval);
        }
        START_TIMER(world);
        double trantol = vtol / std::min(30.0, double(amo.size()));
        normalize(world, amo_new);
        amo_new = transform(world, amo_new, Q3(matrix_inner(world, amo_new, amo_new)), trantol, true);
        truncate(world, amo_new);
        normalize(world, amo_new);
        if(param.nbeta && !param.spin_restricted){
            normalize(world, bmo_new);
            bmo_new = transform(world, bmo_new, Q3(matrix_inner(world, bmo_new, bmo_new)), trantol, true);
            truncate(world, bmo_new);
            normalize(world, bmo_new);
        }
        END_TIMER(world, "Orthonormalize");
        amo = amo_new;
        bmo = bmo_new;
    }

    // For given protocol, solve the DFT/HF/response equations
    void solve(World & world)
    {
        functionT arho_old, brho_old;
        functionT adelrhosq, bdelrhosq;
        const double dconv = std::max(FunctionDefaults<3>::get_thresh(), param.dconv);
        const double trantol = vtol / std::min(30.0, double(amo.size()));
        const double tolloc = 1e-3;
        double update_residual = 0.0, bsh_residual = 0.0;
        subspaceT subspace;
        tensorT Q;
        bool do_this_iter = true;
        // Shrink subspace until stop localizing/canonicalizing
        int maxsub_save = param.maxsub;
        param.maxsub = 2;

        for(int iter = 0;iter < param.maxiter;++iter){
            if(world.rank() == 0)
                printf("\nIteration %d at time %.1fs\n\n", iter, wall_time());

            if (iter > 0 && update_residual < 0.1) {
                //do_this_iter = false;
                param.maxsub = maxsub_save;
            }

            if(param.localize && do_this_iter) {
                tensorT U;
                if (param.localize_pm) {
                    U = localize_PM(world, amo, aset, tolloc, 0.25, iter == 0);
                }
                else {
                    U = localize_boys(world, amo, aset, tolloc, 0.25, iter==0);
                }
                amo = transform(world, amo, U, trantol, true);
                truncate(world, amo);
                normalize(world, amo);
                rotate_subspace(world, U, subspace, 0, amo.size(), trantol);
                if(!param.spin_restricted && param.nbeta){
                    if (param.localize_pm) {
                        U = localize_PM(world, bmo, bset, tolloc, 0.25, iter == 0);
                    }
                    else {
                        U = localize_boys(world, bmo, bset, tolloc, 0.25, iter==0);
                    }
                    bmo = transform(world, bmo, U, trantol, true);
                    truncate(world, bmo);
                    normalize(world, bmo);
                    rotate_subspace(world, U, subspace, amo.size(), bmo.size(), trantol);
                }
            }

            START_TIMER(world);
            functionT arho = make_density(world, aocc, amo);
            functionT brho;
            //if(!param.spin_restricted && param.nbeta)
            //    brho = make_density(world, bocc, bmo);
            //else
            //    brho = arho; // wrong with 1-electron system
            if (param.nbeta) {
                if (param.spin_restricted) brho = arho;
                else brho = make_density(world, bocc, bmo);
            }
            else {
                brho = functionT(world); // zero
            }


            END_TIMER(world, "Make densities");
            if(iter < 2 || (iter % 10) == 0){
                START_TIMER(world);
                loadbal(world, arho, brho, arho_old, brho_old, subspace);
                END_TIMER(world, "Load balancing");
            }
            double da = 0.0, db = 0.0;
            if(iter > 0){
                da = (arho - arho_old).norm2();
                db = (brho - brho_old).norm2();
                if(world.rank() == 0)
                    print("delta rho", da, db, "residuals", bsh_residual, update_residual);

            }
            arho_old = arho;
            brho_old = brho;
            functionT rho = arho + brho;
            rho.truncate();
            double enuclear = inner(rho, vnuc);
            START_TIMER(world);
            functionT vcoul = apply(*coulop, rho);
            END_TIMER(world, "Coulomb");
            double ecoulomb = 0.5 * inner(rho, vcoul);
            rho.clear(false);
            functionT vlocal = vcoul + vnuc;
            vcoul.clear(false);
            vlocal.truncate();
            double exca = 0.0, excb = 0.0;
            vecfuncT Vpsia = apply_potential(world, aocc, amo, arho, brho, adelrhosq, bdelrhosq, vlocal, exca);
            vecfuncT Vpsib;
            if(param.spin_restricted){
                if(!param.lda)  excb = exca;
            }
            else if(param.nbeta) {
                Vpsib = apply_potential(world, bocc, bmo, brho, arho, bdelrhosq, adelrhosq, vlocal, excb);
            }

            double ekina = 0.0, ekinb = 0.0;
            tensorT focka = make_fock_matrix(world, amo, Vpsia, aocc, ekina);
            tensorT fockb = focka;

            if (!param.spin_restricted && param.nbeta)
                fockb = make_fock_matrix(world, bmo, Vpsib, bocc, ekinb);
            else
                ekinb = ekina;

            if (!param.localize && do_this_iter) {
                tensorT U = diag_fock_matrix(world, focka, amo, Vpsia, aeps, dconv);
                rotate_subspace(world, U, subspace, 0, amo.size(), trantol);
                if (!param.spin_restricted && param.nbeta) {
                    U = diag_fock_matrix(world, fockb, bmo, Vpsib, beps, dconv);
                    rotate_subspace(world, U, subspace, amo.size(), bmo.size(), trantol);
                }
            }

            double enrep = molecule.nuclear_repulsion_energy();
            double ekinetic = ekina + ekinb;
            double exc = exca + excb;
            double etot = ekinetic + enuclear + ecoulomb + exc + enrep;
            current_energy = etot;

            if(world.rank() == 0){
                printf("\n              kinetic %16.8f\n", ekinetic);
                printf("   nuclear attraction %16.8f\n", enuclear);
                printf("              coulomb %16.8f\n", ecoulomb);
                printf(" exchange-correlation %16.8f\n", exc);
                printf("    nuclear-repulsion %16.8f\n", enrep);
                printf("                total %16.8f\n\n", etot);
            }

            if(iter > 0){
                //print("##convergence criteria: density delta=", da < dconv * molecule.natom() && db < dconv * molecule.natom(), ", bsh_residual=", (param.conv_only_dens || bsh_residual < 5.0*dconv));
                if(da < dconv * molecule.natom() && db < dconv * molecule.natom() && (param.conv_only_dens || bsh_residual < 5.0*dconv)){
                    if(world.rank() == 0) {
                        print("\nConverged!\n");
                    }

                    // Diagonalize to get the eigenvalues and if desired the final eigenvectors
                    tensorT U;
                    tensorT overlap = matrix_inner(world, amo, amo, true);
                    sygv(focka, overlap, 1, U, aeps);
                    if (!param.localize) {
                        amo = transform(world, amo, U, trantol, true);
                        truncate(world, amo);
                        normalize(world, amo);
                    }
                    if(param.nbeta && !param.spin_restricted){
                        overlap = matrix_inner(world, bmo, bmo, true);
                        sygv(fockb, overlap, 1, U, beps);
                        if (!param.localize) {
                            bmo = transform(world, bmo, U, trantol, true);
                            truncate(world, bmo);
                            normalize(world, bmo);
                        }
                    }

                    if(world.rank() == 0) {
                        print(" ");
                        print("alpha eigenvalues");
                        print(aeps);
                        if(param.nbeta && !param.spin_restricted){
                            print("beta eigenvalues");
                            print(beps);
                        }
                    }

                    if (param.localize) {
                        // Restore the diagonal elements for the analysis
                        for (unsigned int i=0; i<amo.size(); ++i) aeps[i] = focka(i,i);
                        for (unsigned int i=0; i<bmo.size(); ++i) beps[i] = fockb(i,i);
                    }

                    break;
                }

            }

            update_subspace(world, Vpsia, Vpsib, focka, fockb, subspace, Q, bsh_residual, update_residual);
        }

        if(world.rank() == 0) {
            if (param.localize) print("Orbitals are localized - energies are diagonal Fock matrix elements\n");
            else print("Orbitals are eigenvectors - energies are eigenvalues\n");
            print("Analysis of alpha MO vectors");
        }

        analyze_vectors(world, amo, aocc, aeps);
        if(param.nbeta && !param.spin_restricted){
            if(world.rank() == 0)
                print("Analysis of beta MO vectors");

            analyze_vectors(world, bmo, bocc, beps);
        }
    }
};


// Computes molecular energy as a function of the geometry
// This is cludgy ... need better factorization of functionality
// between calculation, main program and this ... or just merge it all.
class MolecularEnergy : public OptimizationTargetInterface {
    World& world;
    Calculation& calc;
    mutable double coords_sum;     // sum of square of coords at last solved geometry
    mutable double E; //< Current energy

public:
    MolecularEnergy(World& world, Calculation& calc)
        : world(world)
        , calc(calc)
        , coords_sum(-1.0)
    {}

    bool provides_gradient() const {return true;}

    double value(const Tensor<double>& x) {
        double xsq = x.sumsq();
        if (xsq == coords_sum) {
            return calc.current_energy;
        }
        calc.molecule.set_all_coords(x.reshape(calc.molecule.natom(),3));
        coords_sum = xsq;

        // The below is missing convergence test logic, etc.

        // Make the nuclear potential, initial orbitals, etc.
        calc.set_protocol(world,1e-4);
        calc.make_nuclear_potential(world);
        calc.project_ao_basis(world);

        //calc.project(world);
        if (calc.param.restart) {
            calc.load_mos(world);
        }
        else {
            calc.initial_guess(world);
            calc.param.restart = true;
        }

        // If the basis for the inital guess was not sto-3g
        // switch to sto-3g since this is needed for analysis
        // of the MOs and orbital localization
        if (calc.param.aobasis != "sto-3g") {
            calc.param.aobasis = "sto-3g";
            calc.project_ao_basis(world);
        }

        calc.solve(world);
        calc.save_mos(world);

        calc.set_protocol(world,1e-6);
        calc.make_nuclear_potential(world);
        calc.project_ao_basis(world);
        calc.project(world);
        calc.solve(world);
        calc.save_mos(world);

        //         calc.set_protocol(world,1e-8);
        //         calc.make_nuclear_potential(world);
        //         calc.project(world);
        //         calc.solve(world);
        //         calc.save_mos(world);

        return calc.current_energy;
    }

    madness::Tensor<double> gradient(const Tensor<double>& x) {
        value(x); // Ensures DFT equations are solved at this geometry

        return calc.derivatives(world);
    }
};



int main(int argc, char** argv) {
    initialize(argc, argv);

    { // limit lifetime of world so that finalize() can execute cleanly
      World world(MPI::COMM_WORLD);

      try {
        // Load info for MADNESS numerical routines
        startup(world,argc,argv);
        FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));

        std::cout.precision(6);

        // Process 0 reads input information and broadcasts
        Calculation calc(world, "input");

        // Warm and fuzzy for the user
        if (world.rank() == 0) {
          print("\n\n");
          print(" MADNESS Hartree-Fock and Density Functional Theory Program");
          print(" ----------------------------------------------------------\n");
          print("\n");
          calc.molecule.print();
          print("\n");
          calc.param.print(world);
        }

        // Come up with an initial OK data map
        if (world.size() > 1) {
          calc.set_protocol(world,1e-4);
          calc.make_nuclear_potential(world);
          calc.initial_load_bal(world);
        }

        MolecularEnergy E(world, calc);
        E.value(calc.molecule.get_all_coords().flat()); // ugh!
        if (calc.param.derivatives) calc.derivatives(world);
        if (calc.param.dipole) calc.dipole(world);
        calc.do_plots(world);

      }
      catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
      }
      catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
      }
      catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
      }
      catch (char* s) {
        print(s);
        error("caught a string exception");
      }
      catch (const char* s) {
        print(s);
        error("caught a string exception");
      }
      catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
      }
      catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
      }
      catch (...) {
        error("caught unhandled exception");
      }

      // Nearly all memory will be freed at this point
      world.gop.fence();
      world.gop.fence();
      ThreadPool::end();
      print_stats(world);
    } // world is dead -- ready to finalize
    finalize();

    return 0;
}
