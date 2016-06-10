/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <boost/python.hpp>
#include <libefp_solver/efp_solver.h>
#include <liboptions/liboptions.h>

using namespace boost::python;
using namespace psi;
using namespace psi::efp;

void export_efp()
{
    // because there is no default constructor for libefp, need flag
    // "no_init" and the constructor definition, def(init<Options&>())
    class_<EFP, boost::shared_ptr<EFP> >("EFP", "Class interfacing with libefp", no_init).
        def(init<Options&>()).
        def("compute", &EFP::compute, "Computes libefp energies and, if active, torque").
        def("set_qm_atoms", &EFP::set_qm_atoms, "Provides libefp with QM fragment information").
        def("nfragments", &EFP::get_frag_count, "Returns the number of EFP fragments in the molecule").
        def("print_out", &EFP::print_out, "Prints options settings and EFP and QM geometries");
}