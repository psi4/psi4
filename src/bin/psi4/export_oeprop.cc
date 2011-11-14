#include <boost/python.hpp>
#include <libmints/oeprop.h>

using namespace boost::python;
using namespace psi;

void export_oeprop()
{
    class_<OEProp, boost::shared_ptr<OEProp> >("OEProp").
        def("add", &OEProp::oepy_add).
        def("compute", &OEProp::oepy_compute).
        def("set_title", &OEProp::oepy_set_title);

    class_<GridProp, boost::shared_ptr<GridProp> >("GridProp").
        def("add", &GridProp::gridpy_add).
        def("set_filename", &GridProp::set_filename).
        def("add_alpha_mo", &GridProp::add_alpha_mo).
        def("add_beta_mo", &GridProp::add_beta_mo).
        def("add_basis_fun", &GridProp::add_basis_fun).
        def("build_grid_overages", &GridProp::build_grid_overages).
        def("set_n", &GridProp::set_n).
        def("set_o", &GridProp::set_o).
        def("set_l", &GridProp::set_l).
        def("get_n", &GridProp::get_n).
        def("get_o", &GridProp::get_o).
        def("get_l", &GridProp::get_l).
        def("set_caxis", &GridProp::set_caxis).
        def("set_format", &GridProp::set_format).
        def("compute", &GridProp::gridpy_compute);
}

