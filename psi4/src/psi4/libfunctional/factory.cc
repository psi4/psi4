/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
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

 #include "psi4/pragma.h"
 PRAGMA_WARNING_PUSH
 PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
 #include <memory>
 PRAGMA_WARNING_POP
#include "psi4/psi4-dec.h"
#include "functional.h"
#include "xfunctional.h"
#include "cfunctional.h"
#include "wpbex_functional.h"
#include "wpbec_functional.h"
#include "LYP_Cfunctional.h"
#include "FT97B_Xfunctional.h"
#include "PZ81_Cfunctional.h"
#include "P86_Cfunctional.h"
#include "PW91_Cfunctional.h"
#include "PBE_Cfunctional.h"
#include "FT97_Cfunctional.h"
#include "VWN3_Cfunctional.h"
#include "VWN5_Cfunctional.h"

namespace psi {

std::shared_ptr<Functional> Functional::build_base(const std::string& alias)
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
    } else if (alias == "B86B_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::B86B;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "PW86_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PW86;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "PBE_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::PBE;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "RPBE_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::RPBE;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "SOGGA_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::SOGGA;
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
        x->gga_ = true;
        x->meta_ = true;
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
        x->gga_ = true;
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPBE_X") {
        wPBEXFunctional* x = new wPBEXFunctional();
        x->set_B88(false);
        fun = static_cast<Functional*>(x);
    } else if (alias == "wB88_X") {
        wPBEXFunctional* x = new wPBEXFunctional();
        x->set_B88(true);
        fun = static_cast<Functional*>(x);
    } else if (alias == "PW92_C") {
        CFunctional* x = new CFunctional();
        x->lsda_type_  = CFunctional::PW92;
        x->gga_type_  = CFunctional::GGA_None;
        x->meta_type_ = CFunctional::Meta_None;
        fun = static_cast<Functional*>(x);
    } else if (alias == "B_C") {
        CFunctional* x = new CFunctional();
        x->lsda_type_  = CFunctional::PW92;
        x->gga_type_  = CFunctional::B97;
        x->meta_type_ = CFunctional::Meta_None;
        x->gga_ = true;
        fun = static_cast<Functional*>(x);
    } else if (alias == "M_C") {
        CFunctional* x = new CFunctional();
        x->lsda_type_  = CFunctional::PW92;
        x->gga_type_  = CFunctional::B97;
        x->meta_type_ = CFunctional::B95;
        x->gga_ = true;
        x->meta_ = true;
        fun = static_cast<Functional*>(x);
    } else if (alias == "LYP_C") {
        fun = new LYP_CFunctional();
    } else if (alias == "FT97B_X") {
        fun = new FT97B_XFunctional();
    } else if (alias == "PZ81_C") {
        fun = new PZ81_CFunctional();
    } else if (alias == "P86_C") {
        fun = new P86_CFunctional();
    } else if (alias == "PW91_C") {
        fun = new PW91_CFunctional();
    } else if (alias == "PBE_C") {
        fun = new PBE_CFunctional();
    } else if (alias == "FT97_C") {
        fun = new FT97_CFunctional();
    } else if (alias == "VWN3_C") {
        fun = new VWN3_CFunctional();
    } else if (alias == "VWN5_C") {
        fun = new VWN5_CFunctional();
    } else if (alias == "PW92A_C") {
        wPBECFunctional* x = new wPBECFunctional();
        x->set_wPBEC_type(wPBECFunctional::pw92c_type);
        fun = static_cast<Functional*>(x);
    } else if (alias == "PBEA_C") {
        wPBECFunctional* x = new wPBECFunctional();
        x->set_wPBEC_type(wPBECFunctional::pbec_type);
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPW92_C") {
        wPBECFunctional* x = new wPBECFunctional();
        x->set_wPBEC_type(wPBECFunctional::pw92c_sr_type);
        fun = static_cast<Functional*>(x);
    } else if (alias == "wPBE_C") {
        wPBECFunctional* x = new wPBECFunctional();
        x->set_wPBEC_type(wPBECFunctional::pbec_sr_type);
        fun = static_cast<Functional*>(x);
    } else if (alias == "HF_X") {
        XFunctional* x = new XFunctional();
        x->gga_type_  = XFunctional::GGA_None;
        x->meta_type_ = XFunctional::Meta_None;
        x->sr_type_   = XFunctional::SR_None;
        fun = static_cast<Functional*>(x);
     }
   else {
        throw PSIEXCEPTION("Functional::build_base: Unrecognized base Functional.");
    }

    return std::shared_ptr<Functional>(fun);
}

}
