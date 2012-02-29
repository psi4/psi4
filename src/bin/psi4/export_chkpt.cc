#include <boost/python.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>

using namespace boost;
using namespace boost::python;
using namespace psi;

void export_chkpt()
{
    class_<Chkpt, boost::shared_ptr<Chkpt> >( "Checkpoint", "docstring", init<PSIO*, int>() ).
        add_property( "enuc", &Chkpt::rd_enuc, &Chkpt::wt_enuc, "docstring").
        add_property( "label", &Chkpt::rd_label, &Chkpt::wt_label, "docstring").
        add_property( "escf", &Chkpt::rd_escf, &Chkpt::wt_escf, "docstring").
        add_property( "eref", &Chkpt::rd_eref, &Chkpt::wt_eref, "docstring").
        add_property( "ecorr", &Chkpt::rd_ecorr, &Chkpt::wt_ecorr, "docstring").
        add_property( "efzc", &Chkpt::rd_efzc, &Chkpt::wt_efzc, "docstring").
        add_property( "etot", &Chkpt::rd_etot, &Chkpt::wt_etot, "docstring").
        add_property( "disp", &Chkpt::rd_disp, &Chkpt::wt_disp, "docstring").
        add_property( "eccsd", &Chkpt::rd_eccsd, &Chkpt::wt_eccsd, "docstring").
        add_property( "e_t", &Chkpt::rd_e_t, &Chkpt::wt_e_t, "docstring").
        add_property( "emp2", &Chkpt::rd_emp2, &Chkpt::wt_emp2, "docstring").
        def( "shared_object", &Chkpt::shared_object, "docstring").
        staticmethod("shared_object");
}
