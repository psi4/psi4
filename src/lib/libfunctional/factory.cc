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
    } else if (alias == "B88_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B88;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "PBE_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PBE;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "PW91_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PW91;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "B97_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B97;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "M_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PBE;
        x->meta_type_ = XFunctional::Becke;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wS_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::GGA_None;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::LSDA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wB88_X_LSDA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B88;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::LSDA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPBE_X_LSDA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PBE;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::LSDA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPW91_X_LSDA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PW91;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::LSDA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wB97_X_LSDA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B97;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::LSDA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wM_X_LSDA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PBE;
        x->meta_type_ = XFunctional::Becke;
        x->sr_type_   = XFunctional::LSDA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wB88_X_GGA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B88;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::GGA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPBE_X_GGA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PBE;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::GGA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPW91_X_GGA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PW91;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::GGA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wB97_X_GGA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B97;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::GGA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wM_X_GGA") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PBE;
        x->meta_type_ = XFunctional::Becke;
        x->sr_type_   = XFunctional::GGA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPBE_X") {
        fun = new wPBEXFunctional(); 
    } else {
        throw PSIEXCEPTION("Functional::build_base: Unrecognized base Functional.");
    }

    return boost::shared_ptr<Functional>(fun);
}

}
