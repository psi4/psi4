#include <boost/python.hpp>
#include <libpsio/psio.hpp>

using namespace boost;
using namespace boost::python;
using namespace psi;

void export_psio()
{
    class_<PSIO, boost::shared_ptr<PSIO> >( "IO" ).
        def( "state", &PSIO::state ).
        def( "open", &PSIO::open ).
        def( "close", &PSIO::close ).
        def( "rehash", &PSIO::rehash ).
        def( "open_check", &PSIO::open_check ).
        def( "tocclean", &PSIO::tocclean ).
        def( "tocprint", &PSIO::tocprint ).
        def( "tocwrite", &PSIO::tocwrite ).
        def( "shared_object", &PSIO::shared_object).
        staticmethod("shared_object").
        def( "get_default_namespace", &PSIO::get_default_namespace).
        staticmethod("get_default_namespace").
        def( "set_default_namespace", &PSIO::set_default_namespace).
        staticmethod("set_default_namespace").
        def( "change_file_namespace", &PSIO::change_file_namespace).
        staticmethod("change_file_namespace");

    class_<PSIOManager, boost::shared_ptr<PSIOManager> >( "IOManager" ).
        def( "shared_object", &PSIOManager::shared_object ).
        staticmethod("shared_object").
        def( "print_out", &PSIOManager::print_out ).
        def( "psiclean", &PSIOManager::psiclean ).
        def( "crashclean", &PSIOManager::crashclean ).
        def( "mark_file_for_retention", &PSIOManager::mark_file_for_retention ).
        def( "write_scratch_file", &PSIOManager::write_scratch_file).
        def( "set_default_path", &PSIOManager::set_default_path ).
        def( "set_specific_path", &PSIOManager::set_specific_path ).
        def( "get_file_path", &PSIOManager::get_file_path ).
        def( "set_specific_retention", &PSIOManager::set_specific_retention ).
        def( "get_default_path", &PSIOManager::get_default_path );
}
