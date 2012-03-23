#include <boost/shared_ptr.hpp>
#include <psi4-dec.h>
#include "functional.h"
#include "xfunctional.h"
#include "wpbex_functional.h"

namespace psi {

boost::shared_ptr<Functional> Functional::build_base(const std::string& alias)
{
    Functional* fun;

    if (alias == "S_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::GGA_None;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else {
        throw PSIEXCEPTION("Functional::build_base: Unrecognized base Functional.");
    }

    return boost::shared_ptr<Functional>(fun);
}

}
