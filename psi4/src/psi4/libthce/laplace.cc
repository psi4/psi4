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

#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility>
#include <iomanip>
#include <regex>

#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.hpp"
#include "thce.h"
#include "laplace.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

namespace psi {

// the third parameter of from_string() should be
// one of std::hex, std::dec or std::oct
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}

LaplaceDenom::LaplaceDenom(std::shared_ptr<Vector> eps_occ,
                           std::shared_ptr<Vector> eps_vir,
                           double delta,
                           double omega,
                           int rank) :
    eps_occ_(eps_occ),
    eps_vir_(eps_vir),
    delta_(delta),
    omega_(omega),
    rank_(rank)
{
}
LaplaceDenom::~LaplaceDenom()
{
}
void LaplaceDenom::compute(const std::string& occ_name, const std::string& vir_name,double power)
{
    if (omega_ != 0.0) {
        throw PSIEXCEPTION("LaplaceDenom: Rob implemented this wrong");
    }

    if (omega_ < 0.0) {
        throw PSIEXCEPTION("LaplaceDenom: omega cannot be negative");
    }

    int nocc = eps_occ_->dimpi()[0];
    int nvir = eps_vir_->dimpi()[0];

    double E_LOMO = eps_occ_->get(0, 0);
    double E_HOMO = eps_occ_->get(0, nocc - 1);
    double E_LUMO = eps_vir_->get(0, 0);
    double E_HUMO = eps_vir_->get(0, nvir - 1);

    double* e_o = eps_occ_->pointer();
    double* e_v = eps_vir_->pointer();

    for (int i = 0; i < nocc; i++) {
        E_LOMO = (E_LOMO < e_o[i] ? E_LOMO : e_o[i]);
        E_HOMO = (E_HOMO > e_o[i] ? E_HOMO : e_o[i]);
    }

    for (int i = 0; i < nvir; i++) {
        E_LUMO = (E_LUMO < e_v[i] ? E_LUMO : e_v[i]);
        E_HUMO = (E_HUMO > e_v[i] ? E_HUMO : e_v[i]);
    }

    double A = rank_ * (E_LUMO - E_HOMO);
    double B = rank_ * (E_HUMO - E_LOMO);
    double R = B / A;

    if (R <= 0.0) {
        throw PSIEXCEPTION("Intrinsic R is 0 or negative.");
    }

    // Pick appropriate quadrature file and read contents
    std::string PSIDATADIR = Process::environment("PSIDATADIR");
    std::string err_table_filename = PSIDATADIR + "/quadratures/1_x/error.bin";
    std::string R_filename = PSIDATADIR + "/quadratures/1_x/R_avail.bin";

    std::ifstream err_table_file(err_table_filename.c_str(), std::ios::in | std::ios::binary);
    std::ifstream R_avail_file(R_filename.c_str(), std::ios::in | std::ios::binary);

    if (!err_table_file)
        throw PSIEXCEPTION("LaplaceQuadrature: Cannot locate error property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/error.bin)");
    if (!R_avail_file)
        throw PSIEXCEPTION("LaplaceQuadrature: Cannot locate R property file for quadrature rules (should be PSIDATADIR/quadratures/1_x/R_avail.bin)");

    int nk = 53;
    int nR = 99;

    // Read in the R available
    double* R_availp = new double[nR];
    R_avail_file.read((char*) R_availp, nR*sizeof(double));

    SharedMatrix err_table(new Matrix("Error Table (nR x nk)", nR, nk));
    double** err_tablep = err_table->pointer();
    err_table_file.read((char*) err_tablep[0], nR*nk*sizeof(double));

    R_avail_file.close();
    err_table_file.close();

    //for (int r2 = 0; r2 < nR; r2++)
    //    outfile->Printf( "  R[%4d] = %20.14E\n", r2+1, R_availp[r2]);
    //err_table->print();

    int indR;
    for (indR = 0; indR < nR; indR++) {
        if (R < R_availp[indR])
            break;
    }
    if (indR == nR) {
        // TODO: Relax this
        throw PSIEXCEPTION("Laplace Quadrature requested for (E_HUMO - E_LOMO)/(E_LUMO-E_HOMO) > 7.0 * 10^12, quadratures are not designed for this range.");
    }

    double accuracy;
    int k, r;
    bool found = false;
    for (k = 0; k < nk; k++) {
        for (r = indR; r < nR; r++) {
            double err = err_tablep[r][k];
            if (err != 0.0 && err < delta_) {
                accuracy = err;
                found = true;
                break;
            }
        }
        if (found)
            break;
    }


    if (!found) {
        throw PSIEXCEPTION("Laplace Quadrature rule could not be found with specified accuracy for this system");
    }

    npoints_ = k + 1;

    // A bit hacky, but OK
    int exponent = (int) floor(log(R_availp[r])/log(10.0));
    int mantissa = (int) round(R_availp[r]/pow(10.0,exponent));
    if (mantissa == 10) {
        exponent++;
        mantissa = 1;
    }

    std::stringstream st;
    st << std::setfill('0');
    st << "1_xk" <<  std::setw(2) << npoints_;
    st << "_" << mantissa;
    st << "E" << exponent;

    std::string quadfile = PSIDATADIR + "/quadratures/1_x/" + st.str().c_str();

    outfile->Printf( "\n  ==> Laplace Denominator <==\n\n");
    outfile->Printf( "  This system has an intrinsic R = (E_HUMO - E_LOMO)/(E_LUMO - E_HOMO) of %7.4E.\n", R);
    outfile->Printf( "  A %d point minimax quadrature with R of %1.0E will be used for the denominator.\n", npoints_, R_availp[r]);
    outfile->Printf( "  The worst-case Chebyshev norm for this quadrature rule is %7.4E.\n", accuracy);
    outfile->Printf( "  Quadrature rule read from file %s.\n\n", quadfile.c_str());

    // The quadrature is defined as \omega_v exp(-\alpha_v x) = 1/x
    double* alpha = new double[npoints_];
    double* omega = new double[npoints_];

    std::vector<std::string> lines;
    std::string text;
    std::ifstream infile(quadfile.c_str());
    if (!infile)
        throw PSIEXCEPTION("LaplaceDenom: Unable to open quadrature rule file: " + quadfile);
    while (infile.good()) {
        getline(infile, text);
        lines.push_back(text);
    }

    std::regex numberline("^\\s*(" NUMBER ").*");
    std::smatch what;

    // We'll be rigorous, the files are extremely well defined
    int lineno = 0;
    for (int index = 0; index < npoints_; index++) {
        std::string line  = lines[lineno++];
        if (!std::regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenom: Unable to read grid file line: \n" + line);
        if (!from_string<double>(omega[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenom: Unable to convert grid file line: \n" + line);
    }
    for (int index = 0; index < npoints_; index++) {
        std::string line  = lines[lineno++];
        if (!std::regex_match(line, what, numberline))
            throw PSIEXCEPTION("LaplaceDenom: Unable to read grid file line: \n" + line);
        if (!from_string<double>(alpha[index], what[1], std::dec))
            throw PSIEXCEPTION("LaplaceDenom: Unable to convert grid file line: \n" + line);
    }

    //for (int k = 0; k < npoints_; k++)
    //    printf("  %24.16E, %24.16E\n", omega[k], alpha[k]);

    // Cast weights back to problem size
    for (int k = 0; k < npoints_; k++) {
        alpha[k] /= A;
        omega[k] /= A;
    }

    tau_occ_ = CoreTensor::build(occ_name, "w", npoints_, "o", nocc);
    tau_vir_ = CoreTensor::build(vir_name, "w", npoints_, "v", nvir);

    double* dop = tau_occ_->pointer();
    double* dvp = tau_vir_->pointer();

    for (int k = 0; k < npoints_; k++) {
        for (int i = 0; i < nocc; i++) {
            dop[k * nocc + i] = exp(alpha[k]*e_o[i]);
        }
        for (int a = 0; a < nvir; a++) {
            dvp[k * nvir + a] = exp(-alpha[k]*e_v[a]);
        }
        C_DSCAL(nocc,pow(omega[k] * exp(-alpha[k] * omega_), 1.0 / (2.0 * rank_)), dop + k * nocc, 1);
        C_DSCAL(nvir,pow(omega[k] * exp(-alpha[k] * omega_), 1.0 / (2.0 * rank_)), dvp + k * nvir, 1);
    }

    for (int k = 0; k < npoints_; k++) {
        for (int i = 0; i < nocc; i++) {
            dop[k * nocc + i] = pow(dop[k * nocc + i],power);
        }
        for (int a = 0; a < nvir; a++) {
            dvp[k * nvir + a] = pow(dvp[k * nvir + a],power);
        }
    }

    delete[] alpha;
    delete[] omega;
    delete[] R_availp;
}

}
