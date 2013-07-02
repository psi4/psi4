/*
  This file is part of MADNESS.

  Copyright (C) 2007,2011 Oak Ridge National Laboratory

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

/// \file mcpfit.cc
/// \brief fitting parameters of Model Core Potential

#include <mra/mra.h>
#include <linalg/solvers.h>
#include <moldft/corepotential.h>
#include <iostream>
#include <iomanip>
#include <set>
#include <cstdio>
using namespace madness;

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
typedef vector<functionT> vecfuncT;
typedef Tensor<double> tensorT;
typedef FunctionFactory<double,3> factoryT;
typedef SeparatedConvolution<double,3> operatorT;
typedef std::shared_ptr<operatorT> poperatorT;

static double ttt, sss;
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

/// Given overlap matrix, return rotation with 3rd order error to orthonormalize the vectors
tensorT Q3(const tensorT& s) {
    tensorT Q = inner(s,s);
    Q.gaxpy(0.2,s,-2.0/3.0);
    for (int i=0; i<s.dim(0); i++) Q(i,i) += 1.0;
    return Q.scale(15.0/8.0);
}

/*
static double dot_product (const tensorT & t1, const tensorT & t2) {
    double s = 0.0;
    for (int i=0; i<t1.dim(0); i++) {
        s += t1[i] * t2[i];
    }
    return s;
}
 */

static tensorT vec2tensor (const vector<double> & v) {
    size_t s = v.size();
    tensorT t(static_cast<long>(s));
    for (unsigned int i=0; i<s; i++) {
        t[i] = v[i];
    }
    return t;
}

static vector<double> tensor2vec (const tensorT & t) {
    size_t s = t.size();
    vector<double> v(s);
    for (unsigned int i=0; i<s; i++) {
        v[i] = t[i];
    }
    return v;
}

void print (const tensorT & t) {
    std::cout << std::scientific << std::setprecision(12);
    std::cout << "[*]" << std::endl;
    for (unsigned int i=0; i<t.dim(0); ++i) {
        printf("%02d: %.12e\n", i, t[i]);
    }
}

class CorePotentialFunctor : public FunctionFunctorInterface<double,3> {
private:
    const CorePotential& pot;
public:
    CorePotentialFunctor(const CorePotential& pot) : pot(pot) {}
    double operator()(const coordT& x) const {
        double r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
        return pot.eval(r);
    }
};

inline static double radius_function (double rsq) {
    //return sqrt(rsq); // r
    return rsq; // r^2
}

class RadiusFunctor : public FunctionFunctorInterface<double,3> {
public:
    RadiusFunctor () {}
    double operator() (const coordT& x) const {
        return radius_function(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    }
};

class RadiusSquareFunctor : public FunctionFunctorInterface<double,3> {
public:
    RadiusSquareFunctor () {}
    double operator() (const coordT& x) const {
        return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
    }
};

class GaussianFunctor : public FunctionFunctorInterface<double,3> {
    double alpha;
public:
    GaussianFunctor (double alpha) : alpha(alpha) {}
    double operator() (const coordT& x) const {
        double rsq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        return exp(-alpha*rsq);
    }
};

class PotentialBasisFunctor : public FunctionFunctorInterface<double,3> {
    int n;
    double alpha, rcut;
    double rn;
public:
    PotentialBasisFunctor (int n, double alpha, double rcut) : n(n), alpha(alpha), rcut(rcut) {}
    double operator() (const coordT& x) const {
        double rsq = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
        double r = sqrt(rsq);
        double sp = smoothed_potential(r*rcut)*rcut;
        double rn = 1.0;
        switch (n) {
            case 0: rn = sp*sp; break;
            case 1: rn = sp; break;
            case 2: rn = 1.0; break;
            case 3: rn = r; break;
            case 4: rn = r*r; break;
            default: rn = pow(r, n - 2);
        }
        return rn * exp(-alpha*rsq);
    }
};

class NcOverR : public FunctionFunctorInterface<double,3> {
    int Nc;
    double c;
public:
    NcOverR (int Nc, double c) : Nc(Nc), c(c) {}
    double operator() (const coordT& x) const {
        return Nc * smoothed_potential(sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]) / c) / c;
    }
};

class CoreOrbitalFunctor : public FunctionFunctorInterface<double,3> {
    CorePotentialManager & cpm;
    unsigned int atn;
    unsigned int core;
    int m;
public:
    CoreOrbitalFunctor (CorePotentialManager & cpm, unsigned int atn, unsigned int core, int m)
        : cpm(cpm), atn(atn), core(core), m(m) {}
    double operator() (const coordT& x) const {
        double xx = x[0];
        double yy = x[1];
        double zz = x[2];
        double rsq = xx*xx + yy*yy + zz*zz;
        return cpm.core_eval(atn, core, m, rsq, xx, yy, zz);
    }
};

struct CalculationParameters {
    double dconv;   ///< conversion criteria
    double L;       ///< box size
    int maxiter;    ///< max number of iteration
    std::string symbol;  ///< symbol of target atom
    double lo;      ///< smallest length scale we need to resolve
    double thresh;  ///< truncation threshold
    double eprec;   ///< precision of smoothing parameter
    double delta;   ///< step size to calc numerical derivative of residual
    int nio;        ///< number of I/O node
    bool nonlinear; ///< if true do non linear optimization
    bool plot;      ///< if true plot vmo or umo

    CalculationParameters()
        : dconv(1e-6)
        , L(50.0)
        , maxiter(50)
        , symbol("Li")
        , lo(1e-10)
        , thresh(1e-6)
        , eprec(1e-4)
        , delta(0.1)
        , nio(1)
        , nonlinear(false)
        , plot(false)
        {}

    void print_info () {
        print(" ** Parameter settings ** ");
        print("            target atom", symbol);
        print("               box size", L);
        print("              lo of mra", lo);
        print("          thresh of mra", thresh);
        print("                  dconv", dconv);
        print("               max iter", maxiter);
        print("eprec for smoothed pot.", eprec);
        print("   delta for derivative", delta);
        print("      number of io node", nio);
        print("non linear optimization", nonlinear);
        print("                   plot", plot);
    }

    void read_file (const std::string filename) {
        std::ifstream f(filename.c_str());
        std::string s;
        while (f >> s) {
            if (s == "atom") {
                f >> symbol;
            }
            else if (s == "conv") {
                f >> dconv;
            }
            else if (s == "maxiter") {
                f >> maxiter;
            }
            else if (s == "L") {
                f >> L;
            }
            else if (s == "thresh") {
                f >> thresh;
            }
            else if (s == "delta") {
                f >> delta;
            }
            else if (s == "nio") {
                f >> nio;
            }
            else if (s == "nonlinear") {
                nonlinear = true;
            }
            else if (s == "plot") {
                plot = true;
            }
            else {
                throw "Unknown input parameter:" + s;
            }
        }
        print_info();
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & dconv & L & maxiter & symbol & lo & thresh & eprec & delta & nio;
    }
};

struct Calculation {
    CalculationParameters param;
    CorePotentialManager corepot;
    bool spin_restricted;
    tensorT aeps, beps;
    tensorT aocc, bocc;
    vector<int> aset, bset;
    vecfuncT amo, bmo;
    vecfuncT vamo, vbmo;
    poperatorT coulop;
    unsigned int atn;               ///< atomic number of target atom
    unsigned int ncore;             ///< number of core orbitals
    unsigned int nalpha, nbeta;     ///< number of valence orbitals alpha/beta
    double vtol;                    ///< multiplication tolerance
    functionT mask;

    Calculation (World & world)
    {
        if (world.rank() == 0) {
            param.read_file("fitinput");
            std::set<unsigned int> atomset;
            atn = symbol_to_atomic_number(param.symbol);
            atomset.insert(atn);
            corepot.read_file("mcp", atomset, param.eprec);
            MADNESS_ASSERT(corepot.is_defined(atn));
        }
        world.gop.broadcast_serializable(param, 0);
        world.gop.broadcast_serializable(corepot, 0);
        FunctionDefaults<3>::set_cubic_cell(-param.L, param.L);
        atn = symbol_to_atomic_number(param.symbol);
        ncore = corepot.n_core_orb(atn);
        set_protocol(world, param.thresh);
        load_mos(world, param.symbol);
        print_info();
        vamo = make_reference(world, amo, aocc);
        if (!spin_restricted && nbeta)
            vbmo = make_reference(world, bmo, bocc);

        // erase core vector
        vecfuncT val_amo = vecfuncT(amo.begin()+ncore, amo.end());
        amo = val_amo;
        val_amo.clear();
        if (!spin_restricted && nbeta) {
            vecfuncT val_bmo = vecfuncT(bmo.begin()+ncore, bmo.end());
            bmo = val_bmo;
            val_bmo.clear();
        }
        print("norm2 of references");
        print(norm2(world, vamo));
        if (!spin_restricted && nbeta)
            print(norm2(world, vbmo));
    }

    void print_info () {
        print("   atomic number :", atn);
        print("  number of core :", ncore);
        print("   nalpha, nbeta :", nalpha, nbeta);
        print(" spin_restricted :", spin_restricted);
        print("            aeps :", aeps);
        print("            aocc :", aocc);
        if (!spin_restricted && nbeta) {
            print("            beps :", beps);
            print("            bocc :", bocc);
        }
    }

    void set_protocol (World & world, double thresh)
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
        mask = functionT(factoryT(world).f(mask3).initial_level(4).norefine());
        if(world.rank() == 0){
            print("\nSolving with thresh", thresh, "    k", FunctionDefaults<3>::get_k(), "   conv", std::max(thresh, param.dconv), "\n");
        }
    }

    void load_mos (World& world, std::string atom)
    {
        const double thresh = FunctionDefaults<3>::get_thresh();
        const int k = FunctionDefaults<3>::get_k();
        unsigned int nmo;
        amo.clear(); bmo.clear();

        std::string refdir = "references/";
        std::string filename = refdir + atom;
        try {
            archive::ParallelInputArchive ar(world, filename.c_str());
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

            ar & spin_restricted;

            ar & nmo;
            nalpha = nmo - ncore;
            nbeta = 0;
            const double trantol = vtol / std::min(30.0, double(nalpha));
            ar & aeps & aocc & aset;
            amo.resize(nmo);
            for (unsigned int i=0; i<amo.size(); i++) ar & amo[i];

            if (amo[0].k() != k) {
                reconstruct(world,amo);
                for(unsigned int i = 0;i < amo.size();i++) amo[i] = madness::project(amo[i], k, thresh, false);
                world.gop.fence();
            }
            normalize(world, amo);
            amo = transform(world, amo, Q3(matrix_inner(world, amo, amo)), trantol, true);
            truncate(world, amo);
            normalize(world, amo);

            if (!spin_restricted) {
                ar & nmo;
                nbeta = nmo - ncore;
                ar & beps & bocc & bset;

                bmo.resize(nmo);
                for (unsigned int i=0; i<bmo.size(); i++) ar & bmo[i];

                if (bmo[0].k() != k) {
                    reconstruct(world,bmo);
                    for(unsigned int i = 0;i < bmo.size();i++) bmo[i] = madness::project(bmo[i], k, thresh, false);
                    world.gop.fence();
                }

                normalize(world, bmo);
                bmo = transform(world, bmo, Q3(matrix_inner(world, bmo, bmo)), trantol, true);
                truncate(world, bmo);
                normalize(world, bmo);
            }
        }
        catch (const madness::MadnessException & e) {
            print(e);
            std::string err_msg("The `mcpfit' requires ");
            err_msg += filename;
            err_msg += ".00000 (and additional numbers of nio) to make reference potential.\n";
            err_msg += "It can be produced by atom SCF using moldft.";
            error(err_msg.c_str());
            throw e;
        }

    }

    vecfuncT make_reference (World & world, vecfuncT & mo, tensorT & occ)
    {
        int nv = mo.size() - ncore;
        vecfuncT vmo = zero_functions<double,3>(world, nv);
        compress(world, vmo);
        reconstruct(world, amo);
        norm_tree(world, amo);
        if (!spin_restricted && bmo.size()) {
            reconstruct(world, bmo);
            norm_tree(world, bmo);
        }

        // coulomb potential
        START_TIMER(world);
        vecfuncT psi_c(amo.begin(), amo.begin()+ncore);
        vecfuncT vsq = square(world, psi_c);
        compress(world, vsq);
        functionT rho_c = factoryT(world);
        rho_c.compress();
        for (unsigned int i=0; i<vsq.size(); i++) {
            if (occ[i])
                rho_c.gaxpy(1.0, vsq[i], aocc[i], false);
        }
        world.gop.fence();
        vsq.clear();
        if (!spin_restricted && bmo.size()) {
            psi_c = vecfuncT(bmo.begin(), bmo.begin()+ncore);
            vsq = square(world, psi_c);
            compress(world, vsq);
            for (unsigned int i=0; i<vsq.size(); i++) {
                if (occ[i])
                    rho_c.gaxpy(1.0, vsq[i], bocc[i], false);
            }
            world.gop.fence();
            vsq.clear();
        }
        else {
            rho_c *= 2.0;
        }
        rho_c.truncate();
        functionT coul = apply(*coulop, rho_c);

        // coulomb - Nc/r (long range term)
        functionT nc_over_r = factoryT(world).functor(functorT(new NcOverR(2*ncore, smoothing_parameter(2*ncore, param.eprec)))).thresh(vtol).initial_level(4);
        nc_over_r.reconstruct();
        coul -= nc_over_r;
        nc_over_r.clear();

        vecfuncT valencemo(mo.begin()+ncore, mo.end());
        gaxpy(world, 1.0, vmo, 1.0, mul_sparse(world, coul, valencemo, vtol));
        truncate(world, vmo);
        END_TIMER(world, "coulomb pot.");

        // exchange potential
        START_TIMER(world);
        psi_c = vecfuncT(mo.begin(), mo.begin()+ncore);
        vecfuncT psif;
        for (unsigned int c=0; c<ncore; c++) {
            for (unsigned int i=ncore; i<mo.size(); i++) {
                psif.push_back(mul_sparse(mo[i], mo[c], vtol, false));
            }
        }
        world.gop.fence();
        truncate(world, psif, vtol);
        psif = apply(world, *coulop, psif);
        truncate(world, psif, vtol);
        reconstruct(world, psif);
        norm_tree(world, psif);
        for (unsigned int c=0; c<ncore; c++) {
            for (unsigned int i=ncore; i<mo.size(); i++) {
                functionT psipsif = mul_sparse(psif[c*nv+i-ncore], mo[c], vtol, false);
                world.gop.fence();
                psipsif.compress();
                vmo[i-ncore].gaxpy(1.0, psipsif, -1.0, false);
            }
        }
        truncate(world, vmo);
        END_TIMER(world, "exchange pot.");

        // shift operator
        START_TIMER(world);
        unsigned int nshell = ncore;
        for (unsigned int c = 0; c < nshell; ++c) {
            unsigned int l = corepot.get_core_l(atn, c);
            unsigned int max_m = (l+1)*(l+2)/2;
            nshell -= max_m - 1;
            double bc = corepot.get_core_bc(atn, c);
            for (unsigned int m = 0; m < max_m; ++m) {
                functionT core = factoryT(world).functor(functorT(new CoreOrbitalFunctor(corepot, atn, c, m)));
                tensorT overlap = inner(world, core, mo);
                print("for m=", m, " core orbital of ", c, "'th shell");
                print("  <core|mo> = ", overlap);
                overlap *= bc;
                for (unsigned int i = ncore; i < mo.size(); ++i) {
                    vmo[i-ncore] -= overlap[i] * core;
                }
            }
        }
        truncate(world, vmo);
        END_TIMER(world, "shift operator");

        return vmo;
    }

    functionT project_potential_basis (World & world, CorePotential & cp, int i)
    {
        int level = 4;
        if (cp.alpha[i] > 100) level += 2;
        //functionT u1 = factoryT(world).functor(functorT(new PotentialBasisFunctor(cp.n[i], cp.alpha[i], cp.rcut))).initial_level(level).truncate_on_project().autorefine();
        functionT u1;
        double norm = 0.0;
        //print("  norm2 of potential basis", i, ":", u1.norm2());
        while (norm < 1e-3) {
            //u1.clear();
            //world.gop.fence();
            u1 = functionT(factoryT(world).functor(functorT(new PotentialBasisFunctor(cp.n[i], cp.alpha[i], cp.rcut))).initial_level(++level).truncate_on_project().noautorefine());
            print("  norm2 of potential basis", i, ":", u1.norm2());
            norm = u1.norm2();
        }
        u1.truncate();

        return u1;
    }

    vecfuncT make_Upsi (World & world, CorePotential & cp, vecfuncT & mo)
    {
        START_TIMER(world);
        functionT umcp = factoryT(world).functor(functorT(new CorePotentialFunctor(cp))).thresh(FunctionDefaults<3>::get_thresh()).initial_level(5);
        umcp.truncate();
        umcp.reconstruct();
        vecfuncT umo = mul_sparse(world, umcp, mo, vtol);
        umcp.clear();
        truncate(world, umo);
        world.gop.fence();
        END_TIMER(world, "make Upsi");
        return umo;
    }

    void plot_z (const char* filename, const functionT & f)
    {
        const int npt = 2049;
        Vector<double,3> lo(vector<double>(3,0.0)), hi(vector<double>(3,0.0));
        hi[2] = 5.0;
        f.reconstruct();
        plot_line(filename, npt, lo, hi, f);
    }

    void plot_p (const char* filename, const vecfuncT& f)
    {
        const int npt = 2049;
        Vector<double,3> lo(vector<double>(3,0.0)), hi(vector<double>(3,0.0));
        hi[2] = 5.0;
        Vector<double,3> h = (hi - lo)*(1.0/(npt-1));

        double sum = 0.0;
        for (int i=0; i<3; i++) sum += h[i]*h[i];
        sum = sqrt(sum);

        World& world = f[0].world();
        reconstruct(world, f);
        if (world.rank() == 0) {
            FILE* file = fopen(filename,"w");
            for (int i=0; i<npt; i++) {
                coordT r = lo + h*double(i);
                fprintf(file, "%.14e ", i*sum);
                plot_line_print_value(file, f[0].eval(r));
                double f0 = f[0].eval(r);
                double f1 = f[1].eval(r);
                double f2 = f[2].eval(r);
                double v = f0*f0+f1*f1+f2*f2;
                fprintf(file, " %.14e", v);
                fprintf(file, "\n");
            }
            fclose(file);
        }
        world.gop.fence();
    }

    CorePotential calc_optimal_coeffs (World & world, CorePotential & cp)
    {
        long nparam = cp.alpha.size() - 1;
        vecfuncT pbs(nparam);
        functionT rf = factoryT(world).functor(functorT(new RadiusFunctor()));
        rf.reconstruct();

        print("alpha: ", cp.alpha);
        print("n: ", cp.n);
        for (unsigned int i=1; i<cp.alpha.size(); i++) {
            printf("making pb for alpha[%d]=%e\n", i, cp.alpha[i]);
            pbs[i-1] = project_potential_basis(world, cp, i);
            pbs[i-1].reconstruct();
            pbs[i-1] = rf * pbs[i-1];
            pbs[i-1].truncate();
        }

        tensorT S(nparam, nparam);
        tensorT b(nparam);
        print("making S, b");
        print("nparam=", nparam, "nalpha,nbeta=", nalpha, nbeta);
        for (unsigned int j=0; j<nalpha; j++) {
            print("for alpha elec", j);
            vecfuncT pb_amo = mul_sparse(world, amo[j], pbs, vtol);
            truncate(world, pb_amo);
            S += matrix_inner(world, pb_amo, pb_amo);
            b += inner(world, rf*vamo[j], pb_amo);
            pb_amo.clear();
            world.gop.fence();
        }
        if (!spin_restricted && nbeta) {
            for (unsigned int j=0; j<nbeta; j++) {
                print("for beta elec", j);
                vecfuncT pb_bmo = mul_sparse(world, bmo[j], pbs, vtol);
                truncate(world, pb_bmo);
                S += matrix_inner(world, pb_bmo, pb_bmo);
                b += inner(world, rf*vbmo[j], pb_bmo);
                pb_bmo.clear();
                world.gop.fence();
            }
        }
        print("S,b=\n", S, b);

        tensorT c, s, sumsq;
        long rank;
        if (world.rank() == 0) gelss(S, b, -1, c, s, rank, sumsq); // solve linear equation
        world.gop.broadcast_serializable(c, 0);
        world.gop.broadcast_serializable(s, 0);
        world.gop.broadcast_serializable(sumsq, 0);
        world.gop.broadcast_serializable(rank, 0);
        print("s,sumsq=", s, sumsq);
        print("rank=", rank);
        print("result potential:");
        printf("        %d    %d    %.6e    %.12e\n", cp.l[0], cp.n[0], cp.alpha[0], cp.A[0]);
        for (unsigned int i=0; i<cp.alpha.size()-1; ++i) {
            printf("        %d    %d    %.6e    %.12e\n", cp.l[i+1], cp.n[i+1], cp.alpha[i+1], c[i]);
        }

        CorePotential result = cp;
        vector<double> A = tensor2vec(c);
        A.insert(A.begin(), cp.A[0]);
        result.A = A;

        return result;
    }

    double compute_residuals (World & world, CorePotential & cp, vecfuncT & mo, vecfuncT & vmo)
    {
        //CorePotential cp_opt = calc_optimal_coeffs(world, cp);
        //cp = cp_opt;
        vecfuncT umo = make_Upsi(world, cp, mo);
        for (unsigned int i=0; i<vmo.size(); i++) {
            print("norm(vmo)", vmo[i].norm2());
            print("norm(umo)", umo[i].norm2());
            //double normratio = vmo[i].norm2() / umo[i].norm2();
            //umo[i] *= normratio;
        }
        vecfuncT rv = sub(world, vmo, umo);
        truncate(world, rv);
        functionT rf = factoryT(world).functor(functorT(new RadiusFunctor()));
        vecfuncT err = mul_sparse(world, rf, rv, vtol);
        if (param.plot) {
            static unsigned int num=0;
            for (unsigned int i=0; i<rv.size(); i++) {
                char fn[256];
                sprintf(fn, "umo%02d.txt", num);
                plot_z(fn, umo[i]);
                sprintf(fn, "vmo%02d.txt", num);
                plot_z(fn, vmo[i]);
                sprintf(fn, "vmo-umo%02d.txt", num);
                plot_z(fn, rv[i]);
                sprintf(fn, "rvmo-rumo%02d.txt", num);
                plot_z(fn, err[i]);
                num++;
            }
        }
        umo.clear();
        rv.clear();
        world.gop.fence();
        //tensorT rnorm = vec2tensor(norm2s(world, err)); // sqrt of norm2's
        tensorT rnorm(static_cast<long>(err.size()));
        for (unsigned int i=0; i<err.size(); i++) rnorm[i] = err[i].norm2();
        print("residual norm:", rnorm);

        double r = 0.0;
        for (int i=0; i<rnorm.dim(0); i++) {
            r += rnorm[i];
        }
        return r / rnorm.dim(0);
    }

    tensorT calc_deriv (World & world, CorePotential & cp, double r)
    {
        size_t nparam = cp.alpha.size();
        tensorT deps(static_cast<long>(nparam));
        for (unsigned int i=1; i<nparam; i++) {
            double a = cp.alpha[i];
            double da = a * param.delta;
            CorePotential cp_mod = cp;
            cp_mod.alpha[i] = a+da;
            cp_mod = calc_optimal_coeffs(world, cp_mod);
            double rmod = compute_residuals(world, cp_mod, amo, vamo);
            if (!spin_restricted && nbeta) rmod += compute_residuals(world, cp_mod, bmo, vbmo);
            deps(i) = (rmod - r) / da;
        }

        return deps;
    }

    CorePotential update (World & world, CorePotential & cp, tensorT & deriv, double & r)
    {
        const double tau = 0.5;

        CorePotential cp_mod = cp;
        double normd = deriv.normf();
        double step = tau / normd;
        if (step > 1) step = 1;
        double rmod;
        while (1) {
            //double rule_value = r - step * normd; // armijo
            double rule_value = r; // simple descent rule
            if (false && rule_value < 0) {
                step *= tau;
                continue;
            }
            if (step < param.dconv) throw "searching vector failed";
            printf("step = %15e, Rule value =%15e\n", step, rule_value);
            vector<double> amod = tensor2vec(vec2tensor(cp.alpha) - step * deriv);
            cp_mod.alpha = amod;
            rmod = compute_residuals(world, cp_mod, amo, vamo);
            if (!spin_restricted && nbeta) rmod += compute_residuals(world, cp_mod, bmo, vbmo);
            printf("r(amod) = %15e\n\n", rmod);

            // wolfe curvature constraint
            //tensorT d2 = calc_deriv(world, cp_mod, rmod);
            //double lhs = dot_product(d2, -deriv);
            //double rhs = -normd;
            //print("deriv");
            //print(d2);
            //printf("lhs:%10e, rhs:%10e\n\n", lhs, rhs);
            printf("rmod:%10e, rule:%10e\n", rmod, rule_value);

            if (rmod < rule_value) break; // found
            //if (rmod < rule_value && lhs > rhs) break; // found

            step *= tau;
        }

        r = rmod;

        return cp_mod;
    }
};

struct CoreFittingTarget : public OptimizationTargetInterface {
    World & world;
    Calculation & calc;
    CorePotential & cp;
    mutable double xsq;
    mutable double r;
    mutable tensorT lastx;

    CoreFittingTarget (World & world, Calculation & calc, CorePotential & cp)
        : world(world), calc(calc), cp(cp), xsq(0.0), r(100.0)
    {}

    bool provides_gradient() const {return true;}

    double value (const tensorT & x) {
        double xsq = x.sumsq();
        if (fabs(xsq - this->xsq) < 1e-10) {
            print("target: using current value (delta xsq =", xsq - this->xsq, ")");
            return this->r;
        }
        this->xsq = xsq;
        lastx = copy(x);
        print("target: calc for alpha = ", x);
        CorePotential cp_mod = cp;
        cp_mod.alpha = tensor2vec(x);
        cp_mod = calc.calc_optimal_coeffs(world, cp_mod);
        double r = calc.compute_residuals(world, cp_mod, calc.amo, calc.vamo);
        if (!calc.spin_restricted && calc.nbeta) r += calc.compute_residuals(world, cp_mod, calc.bmo, calc.vbmo);
        world.gop.fence();
        print("target: value for ", x, " = ", r);
        this->r = r;
        return r;
    }

    tensorT gradient(const tensorT& x) {
        CorePotential cp_mod = cp;
        cp_mod.alpha = tensor2vec(x);
        double r = value(x);
        tensorT g = calc.calc_deriv(world, cp_mod, r);
        world.gop.fence();
        print("target: gradient for ", x, " = ", g);
        return g;
    }
};

class MySteepestDescent : public OptimizerInterface {
    std::shared_ptr<OptimizationTargetInterface> target;
    const double tol;
    const double value_precision;  // Numerical precision of value
    const double gradient_precision; // Numerical precision of each element of residual
    double f;
    double gnorm;

    public:
    MySteepestDescent(const std::shared_ptr<OptimizationTargetInterface>& target,
            double tol = 1e-6,
            double value_precision = 1e-12,
            double gradient_precision = 1e-12)
        : target(target)
          , tol(tol)
          , value_precision(value_precision)
          , gradient_precision(gradient_precision)
          , gnorm(tol*1e16)
    {
        if (!target->provides_gradient()) throw "Steepest descent requires the gradient";
    }

    bool check_positive (tensorT & x)
    {
        for (TensorIterator<double> e = x.unary_iterator(); e != x.end(); ++e) {
            if (*e < 0.0) return false;
        }

        return true;
    }

    bool optimize(Tensor<double>& x)
    {
        double step = 10.0;
        double fnew;
        Tensor<double> g;
        target->value_and_gradient(x,f,g);
        gnorm = g.normf();
        for (int i=0; i<100; i++) {
            while (1) {
                Tensor<double> gnew;
                x.gaxpy(1.0, g, -step);
                if (check_positive(x)) {
                    target->value_and_gradient(x,fnew,gnew);
                    if (fnew < f) {
                        f = fnew;
                        g = gnew;
                        break;
                    }
                }
                x.gaxpy(1.0, g, step);
                step *= 0.5;
                print("reducing step size", f, fnew, "(", f - fnew, ")", step);
            }
            Tensor<double> g = target->gradient(x);
            gnorm = g.normf();
            print("iteration",i,"value",f,"gradient",gnorm);
            if (converged()) break;
        }
        return converged();
    }

    bool converged() const
    {
        return gnorm < tol;
    }

    double gradient_norm() const
    {
        return gnorm;
    }

    double value() const
    {
        return f;
    }

    virtual ~MySteepestDescent(){}
};

int main (int argc, char **argv) {
    initialize(argc, argv);

    {
        World world(MPI::COMM_WORLD);

        try {
            startup(world, argc, argv);
            FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap(world)));
            Calculation calc(world);
            CorePotential cp = calc.corepot.get_potential(calc.atn);

            //double r1 = calc.compute_residuals(world, cp, calc.amo, calc.vamo);
            //if (!calc.spin_restricted) r1 += calc.compute_residuals(world, cp, calc.bmo, calc.vbmo);
            //print("r1=", r1);

            if (!calc.param.nonlinear) {
                print("=== start linear optimization");
                cp = calc.calc_optimal_coeffs(world, cp);
                print("=== linear optimal coeffs prepared");
                double r2 = calc.compute_residuals(world, cp, calc.amo, calc.vamo);
                if (!calc.spin_restricted && calc.nbeta) r2 += calc.compute_residuals(world, cp, calc.bmo, calc.vbmo);
                //print("r1=", r1);
                print("r2=", r2);
            }
            else {
                print("=== start non linear optimization");
                //double tol = 1e-4; // tolerance
                //double prec = calc.param.eprec; // precision
                //double gprec = calc.param.eprec; // gradient precision
                std::shared_ptr<CoreFittingTarget> target_p(new CoreFittingTarget(world, calc, cp));
                //QuasiNewton optimizer(target_p, tol, prec, gprec);
                MySteepestDescent optimizer(target_p);
                tensorT x = vec2tensor(cp.alpha);
                optimizer.optimize(x);
                print("last x = ", target_p->lastx);
                cp.alpha = tensor2vec(target_p->lastx);
                double r = calc.compute_residuals(world, cp, calc.amo, calc.vamo);
                if (!calc.spin_restricted && calc.nbeta) r += calc.compute_residuals(world, cp, calc.bmo, calc.vbmo);
                print("last r=", r);
            }
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
        catch (std::bad_alloc & b) {
            error("operator new failed");
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
    }

    finalize();

    return 0;
}

