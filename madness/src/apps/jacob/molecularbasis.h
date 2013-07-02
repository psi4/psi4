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
#ifndef MOLECULAR_BASIS_H
#define MOLECULAR_BASIS_H

#include <madness_config.h>
#include <constants.h>
#include <moldft/molecule.h>
#include <tinyxml/tinyxml.h>
#include <tensor/tensor.h>
using namespace madness;

#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>

/// Represents a single shell of contracted, Cartesian, Gaussian primitives
class ContractedGaussianShell {
    int type;  ///< Angular momentum = 0, 1, 2, ...
    std::vector<double> coeff;
    std::vector<double> expnt;
    double rsqmax;
    int numbf;  ///< Number of basis functions in shell (type+1)*(type+2)/2

    void normalize() {
        // nwcchem cartesian normalization conventions
        // translation of nmcoeff.F into python and thence to c++
        int np = coeff.size();
        if (np == 1) coeff[0] = 1.0e0;

        double pi32=pow(madness::constants::pi,1.5);
        int l_lim = 2*type - 1;
        double f = 1.0e00;
        for (int n=l_lim; n>1; n-=2) f *= n;

        for (int n=0; n<np; ++n)
            coeff[n] *= pow(2.e0*expnt[n]/madness::constants::pi,0.75e0)*pow(4.e0*expnt[n],0.5E0*type)/sqrt(f);

        double sum = 0.0;
        for (int n1=0; n1<np; ++n1) {
            for (int n2=0; n2<np; ++n2) {
                double S =pi32/pow(expnt[n1]+expnt[n2],1.5e0+type)/pow(2e0,type);
                sum = sum + coeff[n1]*coeff[n2]*S;
            }
        }
        sum *= f;

        f = 1e0/sqrt(sum);
        for (int n=0; n<np; ++n) coeff[n] *= f;
    }

public:
    ContractedGaussianShell()
            : type(-1), coeff(), expnt(), rsqmax(0.0), numbf(0) {};

    ContractedGaussianShell(int type,
                            const std::vector<double>& coeff,
                            const std::vector<double>& expnt,
                            bool donorm=true)
            : type(type), coeff(coeff), expnt(expnt), numbf((type+1)*(type+2)/2) {
        if (donorm) normalize();
        double minexpnt = expnt[0];
        for (unsigned int i=1; i<expnt.size(); ++i)
            minexpnt = std::min(minexpnt,expnt[i]);
        rsqmax = 18.4/minexpnt;  // 18.4 = 8*ln(10)
    }


    /// Returns square of the distance beyond which function is less than 1e-8.
    double rangesq() const {
        return rsqmax;
    }


    /// Evaluates the radial part of the contracted function
    double eval_radial(double rsq) const {
        if (rsq > rsqmax) return 0.0;
        double sum = 0.0;
        for (unsigned int i=0; i<coeff.size(); ++i) {
            double ersq = expnt[i]*rsq;
            if (ersq < 18.4) sum += coeff[i]*exp(-ersq);
        }
        return sum;
    }


    /// Evaluates the entire shell returning the incremented result pointer
    double* eval(double rsq, double x, double y, double z, double* bf) const {
        double R = eval_radial(rsq);
        if (fabs(R) < 1e-8) {
            for (int i=0; i<numbf; ++i) bf[i] = 0.0;

        }
        else {
            switch (type) {
            case 0:
                bf[0] =  R;
                break;
            case 1:
                bf[0] =  R*x;
                bf[1] =  R*y;
                bf[2] =  R*z;
                break;
            case 2:
                static const double fac = 1.0; //sqrt(3.0);
                bf[0] = R*x*x;
                bf[1] = R*x*y*fac;
                bf[2] = R*x*z*fac;
                bf[3] = R*y*y;
                bf[4] = R*y*z*fac;
                bf[5] = R*z*z;
                break;
            case 3:
                bf[0] = R*x*x*x;
                bf[1] = R*x*x*y;
                bf[2] = R*x*x*z;
                bf[3] = R*x*y*y;
                bf[4] = R*x*y*z;
                bf[5] = R*x*z*z;
                bf[6] = R*y*y*y;
                bf[7] = R*y*y*z;
                bf[8] = R*y*z*z;
                bf[9] = R*z*z*z;
                break;

            default:
                throw "UNKNOWN ANGULAR MOMENTUM";
            }
        }
        return bf+numbf;
    }


    /// Returns the shell angular momentum
    int angular_momentum() const {
        return type;
    }

    /// Returns the number of basis functions in the shell
    int nbf() const {
        return numbf;
    }

    /// Returns the number of primitives in the contraction
    int nprim() const {
        return coeff.size();
    }

    /// Returns a const reference to the coefficients
    const std::vector<double>& get_coeff() const {
        return coeff;
    }

    /// Returns a const reference to the exponents
    const std::vector<double>& get_expnt() const {
        return expnt;
    }

    /// Returns a string description of the basis function type
    const char* get_desc(int ibf) const {
        static const char* tags[4][10] = {
            {"s"   ,""    ,""    ,""    ,""    ,""    ,""    ,""    ,""    ,""    } ,
            {"px"  ,"py"  ,"pz"  ,""    ,""    ,""    ,""    ,""    ,""    ,""    } ,
            {"dxx" ,"dxy" ,"dxz" ,"dyy" ,"dyz" ,"dzz" ,""    ,""    ,""    ,""    } ,
            {"fxxx","fxxy","fxxz","fxyy","fxyz","fxzz","fxzz","fyyy","fyzz","fzzz"}
        };
        MADNESS_ASSERT(ibf<numbf && ibf >= 0);
        return tags[type][ibf];
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & type & coeff & expnt & rsqmax & numbf;
    }
};

/// Represents multiple shells of contracted gaussians on a single center
class AtomicBasis {
    std::vector<ContractedGaussianShell> g;
    double rmaxsq;
    int numbf;
    Tensor<double> dmat, avec, bvec;

public:
    AtomicBasis() : g(), rmaxsq(0.0), numbf(0) {};

    AtomicBasis(const std::vector<ContractedGaussianShell>& g)
            : g(g) {
        rmaxsq = 0.0;
        numbf = 0;
        for (unsigned int i=0; i<g.size(); ++i) {
            rmaxsq = std::max(rmaxsq, g[i].rangesq());
            numbf += g[i].nbf();
        }
    }

    void set_guess_info(const Tensor<double>& dmat,
                        const Tensor<double>& avec, const Tensor<double>& bvec) {
        this->dmat = copy(dmat);
        this->avec = copy(avec);
        this->bvec = copy(bvec);
    }

    /// Returns the number of basis functions on the center
    int nbf() const {
        return numbf;
    }

    /// Returns the number of shells on the center
    int nshell() const {
        return g.size();
    }

    /// Returns a const reference to the shells
    const std::vector<ContractedGaussianShell>& get_shells() const {
        return g;
    };

    /// Evaluates the basis functions at point x, y, z relative to atomic center

    /// The array bf[] must be large enough to hold nbf() values.
    ///
    /// Returned is the incremented pointer.
    double* eval(double x, double y, double z, double* bf) const {
        double rsq = x*x + y*y + z*z;
        if (rsq > rmaxsq) {
            for (int i=0; i<numbf; ++i) bf[i] = 0.0;
            return bf+numbf;
        }

        double* bfstart = bf;
        for (unsigned int i=0; i<g.size(); ++i) {
            bf = g[i].eval(rsq, x, y, z, bf);
        }
        // paranoia is good
        MADNESS_ASSERT(bf-bfstart == numbf);
        return bf;
    }

    /// Evaluates the guess atomic density at point x, y, z relative to atomic center
    double eval_guess_density(double x, double y, double z) const {
        MADNESS_ASSERT(has_guess_info());
        double rsq = x*x + y*y + z*z;
        if (rsq > rmaxsq) return 0.0;

        double bf[numbf];
        eval(x, y, z, bf);
        const double* p = dmat.ptr();
        double sum = 0.0;
        for (int i=0; i<numbf; ++i, p+=numbf) {
            double sumj = 0.0;
            for (int j=0; j<numbf; ++j)
                sumj += p[j]*bf[j];
            sum += bf[i]*sumj;
        }
        return sum;
    }

    /// Return shell that contains basis function ibf and also return index of function in the shell
    const ContractedGaussianShell& get_shell_from_basis_function(int ibf, int& ibf_in_shell) const {
        int n=0;
        for (unsigned int i=0; i<g.size(); ++i) {
            int nbf_in_shell = g[i].nbf();
            if (ibf>=n && ibf<(n+nbf_in_shell)) {
                ibf_in_shell = ibf-n;
                return g[i];
            }
            else {
                n += g[i].nbf();
            }
        }
        MADNESS_EXCEPTION("AtomicBasis: get_shell_from_basis_function", ibf*100000 + nbf());
    }

    bool has_guess_info() const {
        return dmat.size()>0;
    }

    const Tensor<double>& get_dmat() const {
        return dmat;
    };

    const Tensor<double>& get_avec() const {
        return avec;
    };

    const Tensor<double>& get_bvec() const {
        return bvec;
    };

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & g & rmaxsq & numbf & dmat & avec & bvec;
    }

};

/// Used to represent one basis function from a shell on a specific center
class AtomicBasisFunction {
private:
    const double xx, yy, zz; // Coordinates of the center
    const ContractedGaussianShell& shell; // Reference to the underlying atomic shell
    const int ibf; // Index of basis function in the shell (0, 1, ...)
    const int nbf; // Number of functions in the shell

public:
    AtomicBasisFunction(double x, double y, double z,
                        const ContractedGaussianShell& shell, int ibf)
            : xx(x), yy(y), zz(z), shell(shell), ibf(ibf), nbf(shell.nbf()) {}


    AtomicBasisFunction(const AtomicBasisFunction& aofunc)
            : xx(aofunc.xx)
            , yy(aofunc.yy)
            , zz(aofunc.zz)
            , shell(aofunc.shell)
            , ibf(aofunc.ibf)
            , nbf(aofunc.nbf) {}

    double operator()(double x, double y, double z) const {
        double bf[nbf];
        x-=xx;
        y-=yy;
        z-=zz;
        double rsq = x*x + y*y + z*z;
        shell.eval(rsq, x, y, z, bf);
        return bf[ibf];
    }

    void print_me(std::ostream& s) const;

    const ContractedGaussianShell& get_shell() const {
        return shell;
    }

    int get_index() const {
        return ibf;
    }

    const char* get_desc() const {
        return shell.get_desc(ibf);
    }

    void get_coords(double& x, double& y, double& z) const {
    	x=xx; y=yy; z=zz;
		return;
    }
};

/// Contracted Gaussian basis
class AtomicBasisSet {
    std::string name;
    std::vector<AtomicBasis> ag;  //< Basis associated by atomic number = 1, 2, ...; 0=Bq.

    template <typename T>
    std::vector<T> load_tixml_vector(TiXmlElement* node, int n, const char* name) {
        TiXmlElement* child = node->FirstChildElement(name);
        MADNESS_ASSERT(child);
        std::istringstream s(child->GetText());
        std::vector<T> r(n);
        for (int i=0; i<n; ++i) {
            MADNESS_ASSERT(s >> r[i]);
        }
        return r;
    }

    template <typename T>
    Tensor<T> load_tixml_matrix(TiXmlElement* node, int n, int m, const char* name) {
        TiXmlElement* child = node->FirstChildElement(name);
        MADNESS_ASSERT(child);
        std::istringstream s(child->GetText());
        Tensor<T> r(n,m);
        for (int i=0; i<n; ++i) {
            for (int j=0; j<m; ++j) {
                MADNESS_ASSERT(s >> r(i,j));
            }
        }
        return r;
    }

public:
    AtomicBasisSet() : name("unknown"), ag(110) {}


    AtomicBasisSet(std::string filename) : name(""), ag(110) {
        read_file(filename);
    }

    void read_file(std::string filename) {
        static const bool debug = false;
        TiXmlDocument doc(filename);
        if (!doc.LoadFile()) {
            std::cout << "AtomicBasisSet: Failed loading from file " << filename
                      << " : ErrorDesc  " << doc.ErrorDesc()
                      << " : Row " << doc.ErrorRow()
                      << " : Col " << doc.ErrorCol() << std::endl;
            MADNESS_EXCEPTION("AtomicBasisSet: Failed loading basis set",0);
        }
        for (TiXmlElement* node=doc.FirstChildElement(); node; node=node->NextSiblingElement()) {
            if (strcmp(node->Value(),"name") == 0) {
                name = node->GetText();
                if (debug) std::cout << "Loading basis set " << name << std::endl;
            }
            else if (strcmp(node->Value(), "basis") == 0) {
                const char* symbol = node->Attribute("symbol");
                if (debug) std::cout << "  found basis set for " << symbol << std::endl;
                int atn = symbol_to_atomic_number(symbol);
                std::vector<ContractedGaussianShell> g;
                for (TiXmlElement* shell=node->FirstChildElement(); shell; shell=shell->NextSiblingElement()) {
                    const char* type = shell->Attribute("type");
                    int nprim=-1;
                    shell->Attribute("nprim",&nprim);
                    if (debug) std::cout << "      found shell " << type << " " << nprim << std::endl;
                    std::vector<double> expnt = load_tixml_vector<double>(shell, nprim, "exponents");
                    if (strcmp(type,"L") == 0) {
                        std::vector<double> scoeff = load_tixml_vector<double>(shell, nprim, "scoefficients");
                        std::vector<double> pcoeff = load_tixml_vector<double>(shell, nprim, "pcoefficients");
                        g.push_back(ContractedGaussianShell(0,scoeff,expnt));
                        g.push_back(ContractedGaussianShell(1,pcoeff,expnt));
                    }
                    else {
                        static const char* tag[] = {"S","P","D","F","G"};
                        int i;
                        for (i=0; i<5; ++i) {
                            if (strcmp(type,tag[i]) == 0) goto foundit;
                        }
                        MADNESS_EXCEPTION("Loading atomic basis set: bad shell type?",0);
foundit:
                        std::vector<double> coeff = load_tixml_vector<double>(shell, nprim, "coefficients");
                        g.push_back(ContractedGaussianShell(i, coeff, expnt));
                    }
                }
                ag[atn] = AtomicBasis(g);
            }
            else if (strcmp(node->Value(), "atomicguess") == 0) {
                const char* symbol = node->Attribute("symbol");
                if (debug) std::cout << "  atomic guess info for " << symbol << std::endl;
                int atn = symbol_to_atomic_number(symbol);
                MADNESS_ASSERT(is_supported(atn));
                int nbf = ag[atn].nbf();
                Tensor<double> dmat = load_tixml_matrix<double>(node, nbf, nbf, "guessdensitymatrix");
                Tensor<double> avec = load_tixml_matrix<double>(node, nbf, nbf, "alphavectors");
                Tensor<double> bvec = load_tixml_matrix<double>(node, nbf, nbf, "betavectors");
                ag[atn].set_guess_info(dmat, avec, bvec);
            }
            else {
                MADNESS_EXCEPTION("Loading atomic basis set: unexpected XML element", 0);
            }
        }

    }


    /// Makes map from atoms to first basis function on atom and number of basis functions on atom
    void atoms_to_bfn(const Molecule& molecule, std::vector<int>& at_to_bf, std::vector<int>& at_nbf) {
        at_to_bf = std::vector<int>(molecule.natom());
        at_nbf   = std::vector<int>(molecule.natom());

        int n = 0;
        for (int i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            at_to_bf[i] = n;
            at_nbf[i] = ag[atn].nbf();
            n += at_nbf[i];
        }
    }


    /// Returns the number of the atom the ibf'th basis function is on
    int basisfn_to_atom(const Molecule& molecule, int ibf) const {
        MADNESS_ASSERT(ibf >= 0);
        int n = 0;
        for (int i=0; i<molecule.natom(); ++i) {
            // Is the desired function on this atom?
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            const int nbf_on_atom = ag[atn].nbf();
            if (ibf >= n  && (n+nbf_on_atom) > ibf) {
                return i;
            }
            else {
                n += nbf_on_atom;
            }
        }
        MADNESS_EXCEPTION("AtomicBasisSet: get_atomic_basis_function: confused?", ibf);
    }

    /// Returns the ibf'th atomic basis function
    AtomicBasisFunction get_atomic_basis_function(const Molecule& molecule, int ibf) const {
        MADNESS_ASSERT(ibf >= 0);
        int n = 0;
        for (int i=0; i<molecule.natom(); ++i) {
            // Is the desired function on this atom?
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            const int nbf_on_atom = ag[atn].nbf();
            if (ibf >= n  && (n+nbf_on_atom) > ibf) {
                int index;
                const ContractedGaussianShell& shell =
                    ag[atn].get_shell_from_basis_function(ibf-n, index);
                return AtomicBasisFunction(atom.x, atom.y, atom.z, shell, index);
            }
            else {
                n += nbf_on_atom;
            }
        }
        MADNESS_EXCEPTION("AtomicBasisSet: get_atomic_basis_function: confused?", ibf);
    }


    /// Given a molecule count the number of basis functions
    int nbf(const Molecule& molecule) const {
        int n = 0;
        for (int i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            MADNESS_ASSERT(is_supported(atn));
            n += ag[atn].nbf();
        }
        return n;
    }

    /// Evaluates the basis functions
    void eval(const Molecule& molecule, double x, double y, double z, double *bf) const {
        for (int i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            bf = ag[atn].eval(x-atom.x, y-atom.y, z-atom.z, bf);
        }
    }


    /// Evaluates the guess density
    double eval_guess_density(const Molecule& molecule, double x, double y, double z) const {
        double sum = 0.0;
        for (int i=0; i<molecule.natom(); ++i) {
            const Atom& atom = molecule.get_atom(i);
            const int atn = atom.atomic_number;
            sum += ag[atn].eval_guess_density(x-atom.x, y-atom.y, z-atom.z);
        }
        return sum;
    }

    bool is_supported(int atomic_number) const {
        return ag[atomic_number].nbf() > 0;
    }

    /// Print basis info for atoms in the molecule (once for each unique atom type)
    void print(const Molecule& molecule) const;

    template <typename T>
    class AnalysisSorter {
        const Tensor<T> v;
    public:
        AnalysisSorter(const Tensor<T>& v) : v(v) {}
        bool operator()(long i, long j) const {
            return std::abs(v[i]) > std::abs(v[j]);
        }
    };

    /// Given a vector of AO coefficients prints an analysis

    /// For each significant coeff it prints
    /// - atomic symbol
    /// - atom number
    /// - basis function type (e.g., dxy)
    /// - basis function number
    /// - MO coeff
    template <typename T>
    void print_anal(const Molecule& molecule, const Tensor<T>& v) {
        const double thresh = 0.2*v.normf();
        if (thresh == 0.0) {
            printf("    zero vector\n");
            return;
        }
        long nbf = int(v.dim(0));
        long list[nbf];
        long ngot=0;
        for (long i=0; i<nbf; ++i) {
            if (std::abs(v(i)) > thresh) {
                list[ngot++] = i;
            }
        }
        std::sort(list,list+ngot,AnalysisSorter<T>(v));

        const char* format;
        if (molecule.natom() < 10) {
            format = "  %2s(%1d)%4s(%2ld)%6.3f  ";
        }
        else if (molecule.natom() < 100) {
            format = "  %2s(%2d)%4s(%3ld)%6.3f  ";
        }
        else if (molecule.natom() < 1000) {
            format = "  %2s(%3d)%4s(%4ld)%6.3f  ";
        }
        else {
            format = "  %2s(%4d)%4s(%5ld)%6.3f  ";
        }
        printf("         ");
        for (long ii=0; ii<ngot; ++ii) {
            long ibf = list[ii];

            const int iat = basisfn_to_atom(molecule, ibf);
            const Atom& atom = molecule.get_atom(iat);
            const AtomicBasisFunction ao = get_atomic_basis_function(molecule, ibf);
            const char* desc = ao.get_desc();
            const char* element = get_atomic_data(atom.atomic_number).symbol;

            // This will need wrapping in a template for a complex MO vector
            printf(format, element, iat, desc, ibf, v[ibf]);
        }
        printf("\n");
    }

    /// Print basis info for all supported atoms
    void print_all() const;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & name & ag;
    }
};



#endif
