#include "functional.h"

namespace psi {

Functional::Functional()
{
    common_init();
}
Functional::~Functional()
{
}
void Functional::common_init()
{
    lrc_ = false;
    gga_ = false;
    meta_ = false;
    name_ = "";
    description_ = "";
    citation_ = "";
    alpha_ = 1.0;
    omega_ = 0.0;
    
    lsda_cutoff_ = 1.0E-20;
    meta_cutoff_ = 1.0E-20;
}
void Functional::set_parameter(const std::string& key, double val) 
{
    parameters_[key] = val;
}
void Functional::print(FILE* out, int level) const
{
    if (level < 1) return;

    fprintf(out, "   => %s Functional <=\n\n", name_.c_str());

    fprintf(out, "%s", description_.c_str());
    fprintf(out, "\n");
    
    fprintf(out, "%s", citation_.c_str());
    fprintf(out, "\n");
    
    fprintf(out, "    GGA   = %14s\n", (gga_ ? "TRUE" : "FALSE"));
    fprintf(out, "    Meta  = %14s\n", (meta_ ? "TRUE" : "FALSE"));
    fprintf(out, "    LRC   = %14s\n", (lrc_ ? "TRUE" : "FALSE"));
    fprintf(out, "    Alpha = %14.6E\n", alpha_);
    fprintf(out, "    Omega = %14.6E\n", omega_);
    fprintf(out, "\n");
    
    if (level > 2) {
        fprintf(out, "    > Parameters <\n\n");
        for (std::map<std::string, double>::const_iterator it = parameters_.begin();  
            it != parameters_.end(); ++it) {
            fprintf(out,"    %11s = %24.16E\n", (*it).first.c_str(), (*it).second);
        }    
        fprintf(out, "\n");
    }
} 

}
