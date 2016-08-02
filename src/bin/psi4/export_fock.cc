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
#include <libmints/mints.h>
#include <libfock/jk.h>
#include <libfock/soscf.h>
#include <libdiis/diismanager.h>

using namespace boost::python;
using namespace psi;

// Just a little patch until we can figure out options python-side.
boost::shared_ptr<JK> py_build_JK(boost::shared_ptr<BasisSet> basis){
    return JK::build_JK(basis, Process::environment.options);
}

void export_fock()
{
    class_<JK, boost::shared_ptr<JK>, boost::noncopyable>("JK", "docstring", no_init)
//            .def(init<boost::shared_ptr<BasisSet> >())
            // .def("build_JK", &JK::build_JK, "docstring")
            .def("build_JK", py_build_JK, "docstring")
            .staticmethod("build_JK")
            .def("initialize", &JK::initialize)
            .def("compute", &JK::compute)
            .def("finalize", &JK::finalize)
            .def("C_left", &JK::C_left, return_internal_reference<>())
            .def("C_right", &JK::C_right, return_internal_reference<>())
            .def("J", &JK::J, return_internal_reference<>())
            .def("K", &JK::K, return_internal_reference<>())
            .def("wK", &JK::wK, return_internal_reference<>())
            .def("D", &JK::D, return_internal_reference<>())
            .def("print_header", &JK::print_header, "docstring");

    class_<SOMCSCF, boost::shared_ptr<SOMCSCF>, boost::noncopyable>("SOMCSCF", "docstring", no_init)
            // .def(init<boost::shared_ptr<JK>, SharedMatrix, SharedMatrix >())
            .def("Ck", &SOMCSCF::Ck)
            .def("rhf_energy", &SOMCSCF::rhf_energy)
            .def("update", &SOMCSCF::update)
            .def("approx_solve", &SOMCSCF::approx_solve)
            .def("solve", &SOMCSCF::solve)
            .def("H_approx_diag", &SOMCSCF::H_approx_diag)
            .def("compute_Hk", &SOMCSCF::Hk)
            .def("compute_Q", &SOMCSCF::compute_Q)
            .def("compute_Qk", &SOMCSCF::compute_Qk)
            .def("compute_AFock", &SOMCSCF::compute_AFock)
            .def("current_total_energy", &SOMCSCF::current_total_energy)
            .def("current_docc_energy", &SOMCSCF::current_docc_energy)
            .def("current_ci_energy", &SOMCSCF::current_ci_energy)
            .def("current_AFock", &SOMCSCF::current_AFock)
            .def("current_IFock", &SOMCSCF::current_IFock)
            .def("zero_redundant", &SOMCSCF::zero_redundant)
            .def("gradient", &SOMCSCF::gradient)
            .def("gradient_rms", &SOMCSCF::gradient_rms);

    class_<DFSOMCSCF, boost::shared_ptr<DFSOMCSCF>, bases<SOMCSCF> >("DFSOMCSCF", "docstring", no_init);
    class_<DiskSOMCSCF, boost::shared_ptr<DiskSOMCSCF>, bases<SOMCSCF> >("DiskSOMCSCF", "docstring", no_init);

    typedef void (DIISManager::*set_wrapper_sm)(SharedMatrix);
    typedef bool (DIISManager::*add_entry_wrapper_sm)(SharedMatrix, SharedMatrix);
    typedef bool (DIISManager::*extrapolate_wrapper_sm)(SharedMatrix);

    class_<DIISManager, boost::shared_ptr<DIISManager>, boost::noncopyable>("DIISManager", "docstring", no_init)
            .def(init<int, const std::string >())
            // .def(init<int, const std::string, int, int>())
            // .def("remove_entry", &DIISManager::remove_entry)
            .def("reset_subspace", &DIISManager::reset_subspace)
            .def("delete_diis_file", &DIISManager::delete_diis_file)
            .def("subspace_size", &DIISManager::subspace_size)
            .def("set_error_vector_size", set_wrapper_sm(&DIISManager::set_error_vector_size))
            .def("set_vector_size", set_wrapper_sm(&DIISManager::set_vector_size))
            .def("add_entry", add_entry_wrapper_sm(&DIISManager::add_entry))
            .def("extrapolate", extrapolate_wrapper_sm(&DIISManager::extrapolate));

}
