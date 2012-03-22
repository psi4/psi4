#include <boost/shared_ptr.hpp>
#include <psi4-dec.h>
#include "functional.h"
#include "xfunctional.h"
#include "wpbex_functional.h"

namespace psi {

boost::shared_ptr<Functional> Functional::build(const std::string& alias)
{
    boost::shared_ptr<Functional> fun;

    if (alias == "S_X") {
        fun = Functional::build_base("S_X"); 
    }

    return fun;
}
boost::shared_ptr<Functional> Functional::build_base(const std::string& alias)
{
    Functional* fun;

    if (alias == "S_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::GGA_None;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;

        x->set_gga(false);
        x->set_meta(false);
        x->set_alpha(1.0);
        x->set_omega(1.0);
        x->set_name("S_X");
        x->set_description("    Slater Exchange\n");
        x->set_citation("    J.C. Slater, Phys. Rev., 81(3):385-390, 1951\n");

        fun = static_cast<Functional*>(x);
    }

    return boost::shared_ptr<Functional>(fun);
}

}
