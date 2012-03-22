#include <boost/shared_ptr.hpp>
#include <psi4-dec.h>
#include "superfunctional.h"
#include "functional.h"

namespace psi {

boost::shared_ptr<SuperFunctional> SuperFunctional::build(const std::string& alias, int max_points, int deriv)
{
    SuperFunctional* super = new SuperFunctional();
    super->set_max_points(max_points);
    super->set_deriv(deriv);

    std::vector<boost::shared_ptr<Functional> >& x_funs = super->x_functionals();
    std::vector<boost::shared_ptr<Functional> >& c_funs = super->c_functionals();

    if (alias == "S_X") {
        super->set_name("S_X");
        super->set_description("    Slater Exchange Functional\n");
        super->set_citation("    J.C. Slater, Phys. Rev., 81(3):385-390, 1951\n");

        x_funs.push_back(Functional::build("S_X"));

        super->set_x_omega(0.0); 
        super->set_c_omega(0.0); 
        super->set_x_alpha(0.0); 
        super->set_c_alpha(0.0); 
    } else {
        throw PSIEXCEPTION("Superfunctional: Unrecognized superfunctional alias");
    }

    return boost::shared_ptr<SuperFunctional>(super);
}

}
