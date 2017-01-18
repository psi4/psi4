/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

#include <cmath>
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/molecule.h"
#include <sstream>
#include <iomanip>
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/coordentry.h"
#include "psi4/psi4-dec.h"
using namespace psi;

namespace {
  double dzero = 0.0;
}

/**
 * Takes a CoordValue object, and returns a string for printing.
 */
std::string variable_to_string(std::shared_ptr<CoordValue>& val, int precision)
{
    std::string valstr;
    if(val->type() == CoordValue::VariableType){
        // A variable, print its name and the value will be printed later
        VariableValue* pval = dynamic_cast<VariableValue*>(val.get());
        if(pval->negated()) valstr = "-";
        valstr += pval->name();
    }else if(val->type() == CoordValue::NumberType){
        // A number, print its value to the requested precision
        std::stringstream ss;
        ss << std::setw(precision+5) << std::setprecision(precision) << std::fixed << val->compute();
        valstr = ss.str();
    }else{
        throw PSIEXCEPTION("Unknown CoordValue type in variable_to_string()");
    }
    return valstr;
}

double VariableValue::compute()
{
    if(geometryVariables_.count(name_) == 0)
        throw PSIEXCEPTION("Variable " + name_ + " used in geometry specification has not been defined");
    return negate_ ? -geometryVariables_[name_] : geometryVariables_[name_];
}

CoordEntry::CoordEntry(int entry_number, double Z, double charge, double mass, const std::string& symbol, const std::string& label)
    : entry_number_(entry_number), computed_(false), Z_(Z),
      charge_(charge), mass_(mass), symbol_(symbol), label_(label), ghosted_(false)
{
}

CoordEntry::CoordEntry(int entry_number, double Z, double charge, double mass, const std::string& symbol,
           const std::string& label, const std::map<std::string, std::string>& basis,
           const std::map<std::string, std::string>& shells)
    : entry_number_(entry_number), computed_(false), Z_(Z),
      charge_(charge), mass_(mass), symbol_(symbol), label_(label), ghosted_(false),
      basissets_(basis), shells_(shells)
{
}

CoordEntry::~CoordEntry()
{
}

const double& CoordEntry::Z() const
{
    if (ghosted_)
        return dzero;
    else
        return Z_;
}


/**
 * Compares the charge, mass, ghost status, and basis set(s) of the current atom with those of another atom
 *
 * @return Whether the two atoms are the same.
 */
bool CoordEntry::is_equivalent_to(const std::shared_ptr<CoordEntry> &other) const
{
    if(other->Z_ != Z_) return false;
    if(other->mass_ != mass_) return false;
    if(other->ghosted_ != ghosted_) return false;
    std::map<std::string, std::string>::const_iterator iter = shells_.begin();
    std::map<std::string, std::string>::const_iterator stop = shells_.end();
    // TODO do we instead care only about orbital basis?
    for(; iter != stop; ++iter){  // compare basis set contents in shells_ not names in basissets_
        std::map<std::string, std::string>::const_iterator other_it = other->shells_.find(iter->first);
        if(other_it == other->shells_.end()) return false; // This basis was never defined for the other atom
        if(iter->second != other_it->second) return false; // The basis sets are different
    }
    return true;
}

void CoordEntry::set_basisset(const std::string& name, const std::string& type)
{
    basissets_[type] = name;
}

const std::string& CoordEntry::basisset(const std::string& type) const
{
    std::map<std::string, std::string>::const_iterator iter = basissets_.find(type);

    if (iter == basissets_.end())
        throw PSIEXCEPTION("CoordEntry::basisset: Basisset not set for "+label_+" and type of " + type);

    return (*iter).second;
}


void CoordEntry::set_shell(const std::string& name, const std::string& type)
{
    shells_[type] = name;
}

const std::string& CoordEntry::shell(const std::string& type) const
{
    std::map<std::string, std::string>::const_iterator iter = shells_.find(type);

    if (iter == shells_.end())
        throw PSIEXCEPTION("CoordEntry::shell: Basisset not set for "+label_+" and type of " + type);

    return (*iter).second;
}


CartesianEntry::CartesianEntry(int entry_number, double Z, double charge, double mass, const std::string& symbol, const std::string& label,
               std::shared_ptr<CoordValue> x, std::shared_ptr<CoordValue> y, std::shared_ptr<CoordValue> z)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label), x_(x), y_(y), z_(z)
{
}

CartesianEntry::CartesianEntry(int entry_number, double Z, double charge, double mass, const std::string& symbol, const std::string& label,
               std::shared_ptr<CoordValue> x, std::shared_ptr<CoordValue> y, std::shared_ptr<CoordValue> z,
               const std::map<std::string, std::string>& basis,
               const std::map<std::string, std::string>& shells)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label, basis, shells), x_(x), y_(y), z_(z)
{
}

void CartesianEntry::print_in_input_format()
{
    std::string xstr(variable_to_string(x_, 12));
    std::string ystr(variable_to_string(y_, 12));
    std::string zstr(variable_to_string(z_, 12));
    outfile->Printf( "  %17s  %17s  %17s\n", xstr.c_str(), ystr.c_str(), zstr.c_str());
}

std::string CartesianEntry::string_in_input_format()
{
    std::string xstr(variable_to_string(x_, 15));
    std::string ystr(variable_to_string(y_, 15));
    std::string zstr(variable_to_string(z_, 15));

    char buffer[120];
    std::stringstream ss;
    sprintf(buffer, "  %17s  %17s  %17s\n", xstr.c_str(), ystr.c_str(), zstr.c_str());
    ss << buffer;
    return ss.str();
}

const Vector3& CartesianEntry::compute()
{
    if(computed_)
        return coordinates_;

    coordinates_[0] = x_->compute();
    coordinates_[1] = y_->compute();
    coordinates_[2] = z_->compute();
    computed_ = true;
    return coordinates_;
}

void CartesianEntry::set_coordinates(double x, double y, double z)
{
    coordinates_[0] = x;
    coordinates_[1] = y;
    coordinates_[2] = z;

    x_->set(x);
    y_->set(y);
    z_->set(z);

    computed_ = true;
}

ZMatrixEntry::ZMatrixEntry(int entry_number, double Z, double charge, double mass, const std::string& symbol, const std::string& label,
                           std::shared_ptr<CoordEntry> rto, std::shared_ptr<CoordValue> rval,
                           std::shared_ptr<CoordEntry> ato, std::shared_ptr<CoordValue> aval,
                           std::shared_ptr<CoordEntry> dto, std::shared_ptr<CoordValue> dval)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label),
      rto_(rto), rval_(rval),
      ato_(ato), aval_(aval),
      dto_(dto), dval_(dval)
{
}

ZMatrixEntry::ZMatrixEntry(int entry_number, double Z, double charge, double mass, const std::string& symbol, const std::string& label,
             const std::map<std::string, std::string>& basis,
             const std::map<std::string, std::string>& shells,
             std::shared_ptr<CoordEntry> rto,
             std::shared_ptr<CoordValue> rval,
             std::shared_ptr<CoordEntry> ato,
             std::shared_ptr<CoordValue> aval,
             std::shared_ptr<CoordEntry> dto,
             std::shared_ptr<CoordValue> dval)
    : CoordEntry(entry_number, Z, charge, mass, symbol, label, basis, shells),
      rto_(rto), rval_(rval),
      ato_(ato), aval_(aval),
      dto_(dto), dval_(dval)
{
}

ZMatrixEntry::~ZMatrixEntry()
{
}


void
ZMatrixEntry::set_coordinates(double x, double y, double z)
{
    coordinates_[0] = fabs(x) < CLEANUP_THRESH ? 0.0 : x;
    coordinates_[1] = fabs(y) < CLEANUP_THRESH ? 0.0 : y;
    coordinates_[2] = fabs(z) < CLEANUP_THRESH ? 0.0 : z;

    if(rto_ != 0){
        if(!rto_->is_computed())
             throw PSIEXCEPTION("Coordinates have been set in the wrong order");
        rval_->set(r(coordinates_, rto_->compute()));
    }

    if(ato_ != 0){
        if(!ato_->is_computed())
             throw PSIEXCEPTION("Coordinates have been set in the wrong order");
        double aval = a(coordinates_, rto_->compute(), ato_->compute());
        // Noise creeps in for linear molecules. Force linearity, if it is close enough.
        double val = 180.0*aval/M_PI;
        aval_->set(val);
    }

    if(dto_ != 0){
        if(!dto_->is_computed())
             throw PSIEXCEPTION("Coordinates have been set in the wrong order");
        double val = d(coordinates_, rto_->compute(), ato_->compute(), dto_->compute());
        // Check for NaN, and don't update if we find one
        if(val == val)
            dval_->set(180.0*val/M_PI);
    }

    computed_ = true;
}

void ZMatrixEntry::print_in_input_format()
{
    if(rto_ == 0 && ato_ == 0 && dto_ == 0){
        /*
         * The first atom
         */
        outfile->Printf( "\n");
    }else if(ato_ == 0 && dto_ == 0){
        /*
         * The second atom
         */
        int rto = rto_->entry_number() + 1;
        std::string rval = variable_to_string(rval_, 6);
        outfile->Printf( "  %5d %11s\n", rto, rval.c_str());
    }else if(dto_ == 0){
        /*
         * The third atom
         */
        int rto = rto_->entry_number() + 1;
        std::string rval = variable_to_string(rval_, 6);
        int ato = ato_->entry_number() + 1;
        std::string aval = variable_to_string(aval_, 6);
        outfile->Printf( "  %5d %11s  %5d %11s\n", rto, rval.c_str(), ato, aval.c_str());
    }else{
        /*
         * Remaining atoms
         */
        int rto = rto_->entry_number() + 1;
        std::string rval = variable_to_string(rval_, 6);
        int ato = ato_->entry_number() + 1;
        std::string aval = variable_to_string(aval_, 6);
        int dto = dto_->entry_number() + 1;
        std::string dval = variable_to_string(dval_, 6);
        outfile->Printf( "  %5d %11s  %5d %11s  %5d %11s\n",
            rto, rval.c_str(), ato, aval.c_str(), dto, dval.c_str());
    }
}


std::string ZMatrixEntry::string_in_input_format()
{
    char buffer[120];
    std::stringstream ss;

    if(rto_ == 0 && ato_ == 0 && dto_ == 0){
        /*
         * The first atom
         */
        sprintf(buffer, "\n");
        ss << buffer;
    }else if(ato_ == 0 && dto_ == 0){
        /*
         * The second atom
         */
        int rto = rto_->entry_number() + 1;
        std::string rval = variable_to_string(rval_, 10);
        sprintf(buffer, "  %5d %11s\n", rto, rval.c_str());
        ss << buffer;
    }else if(dto_ == 0){
        /*
         * The third atom
         */
        int rto = rto_->entry_number() + 1;
        std::string rval = variable_to_string(rval_, 10);
        int ato = ato_->entry_number() + 1;
        std::string aval = variable_to_string(aval_, 10);
        sprintf(buffer, "  %5d %11s  %5d %11s\n", rto, rval.c_str(), ato, aval.c_str());
        ss << buffer;
    }else{
        /*
         * Remaining atoms
         */
        int rto = rto_->entry_number() + 1;
        std::string rval = variable_to_string(rval_, 10);
        int ato = ato_->entry_number() + 1;
        std::string aval = variable_to_string(aval_, 10);
        int dto = dto_->entry_number() + 1;
        std::string dval = variable_to_string(dval_, 10);
        sprintf(buffer, "  %5d %11s  %5d %11s  %5d %11s\n",
            rto, rval.c_str(), ato, aval.c_str(), dto, dval.c_str());
        ss << buffer;
    }
    return ss.str();
}

/**
 * Computes the coordinates of the current atom's entry
 * @return The Cartesian Coordinates, in Bohr
 */
const Vector3& ZMatrixEntry::compute()
{
    if (computed_){
        return coordinates_;
    }

    if(rto_ == 0 && ato_ == 0 && dto_ == 0){
        /*
         * The first atom
         *
         * Place at the origin
         */
        coordinates_[0] = 0.0;
        coordinates_[1] = 0.0;
        coordinates_[2] = 0.0;
    }else if(ato_ == 0 && dto_ == 0){
        /*
         * The second atom
         *
         * Place directly above the first atom
         */
        coordinates_[0] = 0.0;
        coordinates_[1] = 0.0;
        coordinates_[2] = rval_->compute();
    }else if(dto_ == 0){
        /*
         * The third atom
         *
         * The atom specification is
         *      this       rTo   rVal  aTo  aVal
         *        A         B           C
         * We arbitrarily choose to put A pointing upwards.
         */
        double r = rval_->compute();
        double a = aval_->compute() * M_PI/180.0;
        double cosABC = cos(a);
        double sinABC = sin(a);
        const Vector3& B = rto_->compute();
        const Vector3& C = ato_->compute();

        Vector3 eCB = B - C;
        eCB.normalize();
        Vector3 eX, eY;
        if(fabs(1 - fabs(eCB[0])) < 1.0E-5){
            // CB is collinear with X, start by finding Y
            eY[1] = 1.0;
            eX = eY.perp_unit(eCB);
            eY = eX.perp_unit(eCB);
        }else{
            // CB is not collinear with X, we can safely find X first
            eX[0] = 1.0;
            eY = eX.perp_unit(eCB);
            eX = eY.perp_unit(eCB);
        }
        for(int xyz = 0; xyz < 3; ++xyz) {
           coordinates_[xyz] = B[xyz] + r * (eY[xyz] * sinABC - eCB[xyz] * cosABC );
        }
    }else{
        /*
         * The fourth, or subsequent, atom
         *
         * The atom specification is
         *      this       rTo   rVal  aTo  aVal   dTo   dVal
         *        A         B           C           D
         * which allows us to define the vector from C->B (eCB) as the +z axis, and eDC
         * lies in the xz plane.  Then eX, eY and eZ (=eBC) are the x, y, and z axes, respecively.
         */
        double r = rval_->compute();
        double a = aval_->compute() * M_PI/180.0;
        double d = dval_->compute() * M_PI/180.0;
        const Vector3& B = rto_->compute();
        const Vector3& C = ato_->compute();
        const Vector3& D = dto_->compute();

        Vector3 eDC = C - D;
        Vector3 eCB = B - C;
        eDC.normalize();
        eCB.normalize();
        double cosABC = cos(a);
        double sinABC = sin(a);
        double cosABCD = cos(d);
        double sinABCD = sin(d);
        Vector3 eY = eDC.perp_unit(eCB);
        Vector3 eX = eY.perp_unit(eCB);

        for(int xyz = 0; xyz < 3; ++xyz)
           coordinates_[xyz] = B[xyz] + r * (eX[xyz] * sinABC * cosABCD + eY[xyz] * sinABC * sinABCD - eCB[xyz] * cosABC );

//        outfile->Printf( "%5s r %20.14lf, a %20.14lf, d %20.14lf\n", label_.c_str(),
//                r, a, d);
//        outfile->Printf( "      B %20.14lf    %20.14lf    %20.14lf\n", B[0], B[1], B[2]);
//        outfile->Printf( "      C %20.14lf    %20.14lf    %20.14lf\n", C[0], C[1], C[2]);
//        outfile->Printf( "      D %20.14lf    %20.14lf    %20.14lf\n", D[0], D[1], D[2]);
    }

    computed_ = true;
    return coordinates_;
}
