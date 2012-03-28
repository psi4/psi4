#include <boost/shared_ptr.hpp>
#include <psi4-dec.h>
#include "functional.h"
#include "xfunctional.h"
#include "wpbex_functional.h"
#include "LYP_Cfunctional.h"
#include "FT97B_Xfunctional.h"
#include "PZ81_Cfunctional.h"
#include "P86_Cfunctional.h"
#include "VWN_Cfunctional.h"
#include "PW91_Cfunctional.h"
#include "PW92_Cfunctional.h"
#include "PBE_Cfunctional.h"
#include "FT97_Cfunctional.h"
#include "B972_Cfunctional.h"
#include "B974_Cfunctional.h"

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
    } else if (alias == "wB97_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B97;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::LSDA;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPBE_X") {
        fun = new wPBEXFunctional(); 
    } else if (alias == "LYP_C") {
        fun = new LYP_CFunctional();
    } else if (alias == "FT97B_X") {
        fun = new FT97B_XFunctional();
    } else if (alias == "PZ81_C") {
        fun = new PZ81_CFunctional();
    } else if (alias == "P86_C") {
        fun = new P86_CFunctional();
    } else if (alias == "VWN_C") {
        fun = new VWN_CFunctional();
    } else if (alias == "PW91_C") {
        fun = new PW91_CFunctional();
    } else if (alias == "PW92_C") {
        fun = new PW92_CFunctional();
    } else if (alias == "PBE_C") {
        fun = new PBE_CFunctional();
    } else if (alias == "FT97_C") {
        fun = new FT97_CFunctional();
    } else if (alias == "B972_C") {
        fun = new B972_CFunctional();
    } else if (alias == "B974_C") {
        fun = new B974_CFunctional();
    } else {
        throw PSIEXCEPTION("Functional::build_base: Unrecognized base Functional.");
    }

    return boost::shared_ptr<Functional>(fun);
}

}
