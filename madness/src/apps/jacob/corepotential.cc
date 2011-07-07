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

/// \file corepotential.cc
/// \brief Simple management of core potential and orbital information

#include <madness_config.h>
#include <constants.h>
#include <mra/mra.h>
#include <tinyxml/tinyxml.h>
#include <moldft/corepotential.h>
#include <cstdio>
#include <cmath>
#include <string>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <set>
using std::string;
using std::vector;
using namespace madness;

typedef Vector<double,3> coordT;
typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;
typedef Function<double,3> functionT;
typedef FunctionFactory<double,3> factoryT;
typedef vector<functionT> vecfuncT;

const double CorePotentialManager::fc = 2.0;

double CorePotential::eval(double r) const {
    double u = 0.0;
    double rr = r*r;
    double sp_n = smoothed_potential(r*rcut)*rcut;
    for (unsigned int i=0; i<A.size(); ++i) {
        double rn = 1.0;
        double sp = sp_n;
        if (i==0) {
            sp = smoothed_potential(r*rcut0)*rcut0;
        }
        switch (n[i]) {
            //case 0: rn = sp*sp; break;
            case 1: rn = sp; break;
            case 2: rn = 1.0; break;
            //case 3: rn = r; break;
            //case 4: rn = rr; break;
            default: rn = pow(r, n[i] - 2);
        }
        u += A[i] * rn * exp(-alpha[i] * rr);
    }
    return u;
}

double CorePotential::eval_derivative(double xi, double r) const {
    double u = 0.0;
    double rr = r*r;
    double sp_n = smoothed_potential(r*rcut)*rcut;
    double dsp_n = -dsmoothed_potential(r*rcut)*rcut*rcut;
    for (unsigned int i=0; i<A.size(); ++i) {
        double rn2 = 1.0;
        double rn4 = 1.0;
        double sp = sp_n;
        double dsp = dsp_n;
        if (i==0) {
            sp = smoothed_potential(r*rcut0)*rcut0;
            dsp = -dsmoothed_potential(r*rcut0)*rcut0*rcut0;
        }
        switch (n[i]) {
            //case 0: rn2 = sp*sp; rn4 = rn2*rn2; break;
            case 1: rn2 = sp; rn4 = dsp*sp; break;
            case 2: rn2 = 1.0; rn4 = dsp; break;
            //case 3: rn2 = r; rn4 = sp; break;
            //case 4: rn2 = rr; rn4 = 1.0; break;
            default: rn2 = pow(r, n[i] - 2); rn4 = pow(r, n[i] - 4);
        }
        u += A[i] * xi * exp(-alpha[i] * rr) * ((n[i]-2) * rn4 - 2.0 * alpha[i] * rn2);
    }
    return u;
}

string CorePotential::to_string () const {
    std::ostringstream oss;
    for (unsigned int i=0; i<A.size(); ++i) {
        oss.precision(8);
        oss << std::scientific;
        std::string sep = "    ";
        oss << l[i] << sep << n[i] << sep << alpha[i] << sep << A[i] << endl;
    }
    return oss.str();
}

double CoreOrbital::eval_radial(double rsq) const {
    double s=0.0;
    for (unsigned int k=0; k<expnt.size(); ++k) {
        s += coeff[k] * pow((2 * expnt[k] / madness::constants::pi), 0.75) * exp(-1.0 * expnt[k] * rsq);
    }
    return s;
}

double CoreOrbital::eval_radial_derivative(double rsq, double xi) const {
    double s=0.0;
    for (unsigned int k=0; k<expnt.size(); ++k) {
        s += coeff[k] * pow((2 * expnt[k] / madness::constants::pi), 0.75) * exp(-1.0 * expnt[k] * rsq) * (-2.0 * expnt[k] * xi);
    }
    return s;
}

double CoreOrbital::eval_spherical_harmonics(int m, double x, double y, double z, double& dp, int axis=0) const {
    double p = 1.0;
    dp = 0.0;
    switch (type) {
        case 0:
            break;
        case 1:
            switch (m) {
                case 0: p *= x; if (axis == 0) dp = 1.0; break;
                case 1: p *= y; if (axis == 1) dp = 1.0; break;
                case 2: p *= z; if (axis == 2) dp = 1.0; break;
                default: throw "INVALID MAGNETIC QUANTUM NUMBER";
            }
            break;
        case 2:
            static const double fac = 1.0; //sqrt(3.0);
            switch (m) {
                case 0: p *= x*x; if (axis == 0) dp = 2*x; break;
                case 1: p *= x*y*fac;
                        if (axis == 0) dp = y*fac;
                        else if (axis == 1) dp = x*fac;
                        break;
                case 2: p *= x*z*fac;
                        if (axis == 0) dp = z*fac;
                        else if (axis == 2) dp = x*fac;
                        break;
                case 3: p *= y*y; if (axis == 1) dp = 2*y; break;
                case 4: p *= y*z*fac;
                        if (axis == 1) dp = z*fac;
                        else if (axis == 2) dp = y*fac;
                        break;
                case 5: p *= z*z; if (axis == 2) dp = 2*z; break;
                default: throw "INVALID MAGNETIC QUANTUM NUMBER";
            }
            break;
        case 3:
            switch (m) {
                case 0: p *= x*x*x;
                        if (axis == 0) dp = 3*x*x;
                        break;
                case 1: p *= x*x*y;
                        if (axis == 0) dp = 2*x*y;
                        else if (axis == 1) dp = x*x;
                        break;
                case 2: p *= x*x*z;
                        if (axis == 0) dp = 2*x*z;
                        else if (axis == 2) dp = x*x;
                        break;
                case 3: p *= x*y*y;
                        if (axis == 0) dp = y*y;
                        else if (axis == 1) dp = 2*x*y;
                        break;
                case 4: p *= x*y*z;
                        if (axis == 0) dp = y*z;
                        else if (axis == 1) dp = x*z;
                        else dp = x*y;
                        break;
                case 5: p *= x*z*z;
                        if (axis == 0) dp = z*z;
                        else if (axis == 2) dp = 2*x*z;
                        break;
                case 6: p *= y*y*y;
                        if (axis == 1) dp = 3*y*y;
                        break;
                case 7: p *= y*y*z;
                        if (axis == 1) dp = 2*y*z;
                        else if (axis == 2) dp = y*y;
                        break;
                case 8: p *= y*z*z;
                        if (axis == 1) dp = z*z;
                        else if (axis == 2) dp = 2*y*z;
                        break;
                case 9: p *= z*z*z;
                        if (axis == 2) dp = 3*z*z;
                        break;
                default: throw "INVALID MAGNETIC QUANTUM NUMBER";
            }
            break;

        default:
            throw "UNKNOWN ANGULAR MOMENTUM";
    }
    return p;
}

double CoreOrbital::eval(int m, double rsq, double x, double y, double z) const {
    if (m < 0 || m >= (type+1)*(type+2)/2) throw "INVALID MAGNETIC QUANTUM NUMBER";
    double R = eval_radial(rsq);
    if (fabs(R) < 1e-8) {
        return 0.0;
    }
    double dummy;
    double p = eval_spherical_harmonics(m, x, y, z, dummy);
    return R*p;
}

double CoreOrbital::eval_derivative(int m, int axis, double xi, double rsq, double x, double y, double z) const {
    if (m < 0 || m >= (type+1)*(type+2)/2) throw "INVALID MAGNETIC QUANTUM NUMBER";
    double R = eval_radial(rsq);
    double dR = eval_radial_derivative(rsq, xi);
    if (fabs(R) < 1e-8) {
        return 0.0;
    }
    double dp;
    double p = eval_spherical_harmonics(m, x, y, z, dp, axis);
    return dR*p + R*dp;
}

static const string dir = "coredata/";

static bool read_potential(TiXmlElement* elem, AtomCore& ac, double eprec) {
    TiXmlElement* p = elem->FirstChildElement("potential");
    if (!p) return false;

    std::istringstream iss(p->GetText());
    int l, n;
    double e, c;
    vector<int> vl, vn;
    vector<double> ve, vc;
    while (iss >> l) {
        iss >> n >> e >> c;
        if (l<0) continue;
        vl.push_back(l);
        vn.push_back(n);
        ve.push_back(e);
        vc.push_back(c);
    }
    ac.potential.l = vl;
    ac.potential.n = vn;
    ac.potential.A = vc;
    ac.potential.alpha = ve;
    ac.potential.eprec = eprec;
    int atn = ac.atomic_number;
    int ncore = ac.ncore;
    //ac.potential.rcut0 = 1.0/smoothing_parameter(ncore*2, eprec);
    ac.potential.rcut0 = 1.0/smoothing_parameter(ncore*2, 1.0);
    //ac.potential.rcut = 1.0/smoothing_parameter(atn-ncore*2, eprec);
    ac.potential.rcut = 1.0/smoothing_parameter(atn-ncore*2, 1.0);

    return true;
}

static bool read_orbital(TiXmlElement* e, AtomCore& ac) {
    TiXmlElement* p = e->FirstChildElement("core");
    if (!p) return false;

    std::istringstream i_num(p->Attribute("num"));
    i_num >> ac.ncore;

    vector<CoreOrbital> vc;

    for (TiXmlElement* node = p->FirstChildElement("orbital"); node; node = node->NextSiblingElement("orbital")) {
        int type;
        vector<double> coeff, expnt;
        double c, e;
        double bc;
        std::istringstream i_bc(node->Attribute("bc"));
        i_bc >> bc;
        std::istringstream i_type(node->Attribute("type"));
        i_type >> type;
        std::istringstream iss(node->GetText());
        while (iss >> e) {
            iss >> c;
            coeff.push_back(c);
            expnt.push_back(e);
        }
        CoreOrbital co(type, coeff, expnt, bc);
        vc.push_back(co);
    }
    ac.orbital = vc;

    return true;
}

static AtomCore read_atom(TiXmlElement* e, unsigned int atn, double eprec) {
    AtomCore ac;
    ac.atomic_number = atn;

    if (!read_orbital(e, ac)) {
        MADNESS_EXCEPTION("CORE_INFO: read_orbital failed.", -1);
    }
    if (!read_potential(e, ac, eprec)) {
        MADNESS_EXCEPTION("CORE_INFO: read_potential failed.", -1);
    }

    return ac;
}

void CorePotentialManager::read_file(string filename, std::set<unsigned int> atomset, double eprec) {
    core_type = filename;
    guess_filename = dir + filename + "_guess";

    TiXmlDocument doc(dir + core_type);
    if (!doc.LoadFile()) {
        MADNESS_EXCEPTION(("CORE_INFO: Failed to load core_info data file: " + dir + core_type).c_str(), -1);
        return;
    }
    TiXmlElement* core_info = doc.FirstChildElement();
    if (!core_info) {
        MADNESS_EXCEPTION("CORE_INFO: core_info data file is not valid.", -1);
        return;
    }

    for (TiXmlElement* node = core_info->FirstChildElement("atom"); node; node = node->NextSiblingElement("atom")) {
        unsigned int atn = symbol_to_atomic_number(node->Attribute("symbol"));
        if (atomset.find(atn) != atomset.end()) {
            AtomCore ac = read_atom(node, atn, eprec);
            if (ac.n_orbital() == 0) {
                MADNESS_EXCEPTION("CORE_INFO: read_atom Failed.", -1);
                return;
            }
            atom_core.insert(std::pair<unsigned int,AtomCore>(atn, ac));
        }
    }

    vector<unsigned int> atns;
    for (std::map<unsigned int, AtomCore>::iterator it = atom_core.begin(); it != atom_core.end(); ++it) {
        atns.push_back(it->first);
    }
    madness::print("MCP parameters loaded for atomic numbers:", atns);
}

void CorePotentialManager::set_eprec(double value) {
    for (std::map<unsigned int,AtomCore>::iterator it=atom_core.begin(); it != atom_core.end(); ++it) {
        it->second.potential.eprec = value;
        double q0 = it->second.ncore * 2;
        double q = it->first - it->second.ncore * 2;
        it->second.potential.rcut0 = 1.0 / smoothing_parameter(q0, value);
        it->second.potential.rcut = 1.0 / smoothing_parameter(q, value);
    }
}

void CorePotentialManager::set_rcut(double value) {
    for (std::map<unsigned int,AtomCore>::iterator it=atom_core.begin(); it != atom_core.end(); ++it) {
        it->second.potential.rcut0 = (value<=0.0) ? 1.0 : value;
        it->second.potential.rcut = (value<=0.0) ? 1.0 : value;
    }
}
