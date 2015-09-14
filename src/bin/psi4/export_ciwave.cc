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

#include <boost/python.hpp>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <../bin/detci/ciwave.h>

using namespace boost::python;
using namespace psi;
using namespace psi::detci;

void export_ciwave()
{
    class_<CIWavefunction, boost::shared_ptr<CIWavefunction> >("CIWavefunction", "docstring", no_init)
        .def(init<boost::shared_ptr<Wavefunction> >())
        .def("get_active_tpdm", &CIWavefunction::get_active_tpdm, "docstring")
        .def("compute_mcscf", &CIWavefunction::compute_mcscf, "docstring");
}

