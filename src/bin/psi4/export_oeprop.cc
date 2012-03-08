#include <boost/python.hpp>
#include <libmints/oeprop.h>

using namespace boost::python;
using namespace psi;

void export_oeprop()
{
    class_<OEProp, boost::shared_ptr<OEProp> >("OEProp", "docstring").
        def("add", &OEProp::oepy_add, "docstring").
        def("compute", &OEProp::oepy_compute, "docstring").
        def("set_title", &OEProp::oepy_set_title, "docstring");

    class_<GridProp, boost::shared_ptr<GridProp> >("GridProp", "docstring").
        def("add", &GridProp::gridpy_add, "docstring").
        def("set_filename", &GridProp::set_filename, "docstring").
        def("add_alpha_mo", &GridProp::add_alpha_mo, "docstring").
        def("add_beta_mo", &GridProp::add_beta_mo, "docstring").
        def("add_basis_fun", &GridProp::add_basis_fun, "docstring").
        def("build_grid_overages", &GridProp::build_grid_overages, "docstring").
        def("set_n", &GridProp::set_n, "docstring").
        def("set_o", &GridProp::set_o, "docstring").
        def("set_l", &GridProp::set_l, "docstring").
        def("get_n", &GridProp::get_n, "docstring").
        def("get_o", &GridProp::get_o, "docstring").
        def("get_l", &GridProp::get_l, "docstring").
        def("set_caxis", &GridProp::set_caxis, "docstring").
        def("set_format", &GridProp::set_format, "docstring").
        def("compute", &GridProp::gridpy_compute, "docstring");
}

