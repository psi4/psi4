#include <boost/python.hpp>
#include <libpsio/psio.hpp>

using namespace boost;
using namespace boost::python;
using namespace psi;

void export_psio()
{
    class_<PSIO, boost::shared_ptr<PSIO> >( "IO", "docstring" ).
        def( "state", &PSIO::state, "docstring" ).
        def( "open", &PSIO::open, "docstring" ).
        def( "close", &PSIO::close, "docstring" ).
        def( "rehash", &PSIO::rehash, "docstring" ).
        def( "open_check", &PSIO::open_check, "docstring" ).
        def( "tocclean", &PSIO::tocclean, "docstring" ).
        def( "tocprint", &PSIO::tocprint, "docstring" ).
        def( "tocwrite", &PSIO::tocwrite, "docstring" ).
        def( "shared_object", &PSIO::shared_object).
        def( "set_pid", &PSIO::set_pid, "docstring" ).
        staticmethod("shared_object").
        def( "get_default_namespace", &PSIO::get_default_namespace, "docstring").
        staticmethod("get_default_namespace").
        def( "set_default_namespace", &PSIO::set_default_namespace, "docstring").
        staticmethod("set_default_namespace").
        def( "change_file_namespace", &PSIO::change_file_namespace, "docstring").
        staticmethod("change_file_namespace");

    class_<PSIOManager, boost::shared_ptr<PSIOManager> >( "IOManager", "docstring" ).
        def( "shared_object", &PSIOManager::shared_object, "docstring" ).
        staticmethod("shared_object").
        def( "print_out", &PSIOManager::print_out, "docstring" ).
        def( "psiclean", &PSIOManager::psiclean, "docstring" ).
        def( "crashclean", &PSIOManager::crashclean, "docstring" ).
        def( "mark_file_for_retention", &PSIOManager::mark_file_for_retention, "docstring" ).
        def( "write_scratch_file", &PSIOManager::write_scratch_file, "docstring").
        def( "set_default_path", &PSIOManager::set_default_path, "docstring" ).
        def( "set_specific_path", &PSIOManager::set_specific_path, "docstring" ).
        def( "get_file_path", &PSIOManager::get_file_path, "docstring" ).
        def( "set_specific_retention", &PSIOManager::set_specific_retention, "docstring" ).
        def( "get_default_path", &PSIOManager::get_default_path, "docstring" );
}
