/*
 * EFP solver
 */ 

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <vector>
#include <utility> 
             
#include <libmints/mints.h>
#include <psifiles.h>
#include <physconst.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libiwl/iwl.hpp>
#include <libqt/qt.h>
#include <liboptions/python.h>
#include <psifiles.h>

#include "efp_solver.h"
//#include "../../../libefp/install/include/efp.h"

#include <psi4-dec.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace boost;
using namespace std;
using namespace psi;

// TODO: change allocated memory to shared pointers and ditch the deletes

namespace psi { namespace efp {
EFP::EFP(Options& options): options_(options) 
{
    common_init();
}

EFP::~EFP(){
       for (int i = 0; frag_name[i]; i++)
               delete[] frag_name[i];
       delete[] frag_name;

       efp_shutdown(efp);
}

static int string_compare(const void *a, const void *b)
{
       const char *s1 = *(const char **)a;
       const char *s2 = *(const char **)b;

       return strcasecmp(s1, s2);
}

/* this is from libefp/efpmd/parse.c */
static char **make_potential_file_list(const char **frag_name,
                                      const char *fraglib_path,
                                      const char *userlib_path)
{
       /* This function constructs the list of library potential data files.
        * For each unique fragment if fragment name contains an _l suffix
        * append fraglib_path prefix and remove _l suffix. Otherwise append
        * userlib_path prefix. Add .efp extension in both cases. */

       int n_frag = 0;

       while (frag_name[n_frag])
               n_frag++;

       const char **unique = new const char*[n_frag];

       for (int i = 0; i < n_frag; i++)
               unique[i] = frag_name[i];

       qsort(unique, n_frag, sizeof(char *), string_compare);

       int n_unique = 1;

       for (int i = 1; i < n_frag; i++)
               if (strcasecmp(unique[i - 1], unique[i])) {
                       unique[n_unique] = unique[i];
                       n_unique++;
               }

       char **list = new char*[n_unique + 1];

       for (int i = 0; i < n_unique; i++) {
               const char *name = unique[i];
               size_t len = strlen(name);
               char *path;

               if (len > 2 && name[len - 2] == '_' && name[len - 1] == 'l') {
                       path = new char[strlen(fraglib_path) + len + 4];
                       strcat(strcpy(path, fraglib_path), "/");
                       strcat(strncat(path, name, len - 2), ".efp");
               }
               else {
                       path = new char[strlen(userlib_path) + len + 6];
                       strcat(strcpy(path, userlib_path), "/");
                       strcat(strcat(path, name), ".efp");
               }

               list[i] = path;
       }

       list[n_unique] = NULL;
       delete[] unique;
       return list;
}

void EFP::common_init() {
       enum efp_result res;
       int print = options_.get_int("PRINT");

       struct efp_opts opts;
       memset(&opts, 0, sizeof(struct efp_opts));

       elst_enabled = options_.get_bool("EFP_ELST");
       pol_enabled  = options_.get_bool("EFP_POL");
       disp_enabled = options_.get_bool("EFP_DISP");
       exch_enabled = options_.get_bool("EFP_EXCH");

       std::string dertype = options_.get_str("DERTYPE");
       do_grad = false;
       if (dertype == "FIRST")
           do_grad = true;

       if (elst_enabled)
               opts.terms |= EFP_TERM_ELEC;
       if (pol_enabled)
               opts.terms |= EFP_TERM_POL;
       if (disp_enabled)
               opts.terms |= EFP_TERM_DISP;
       if (exch_enabled)
               opts.terms |= EFP_TERM_XR;

       std::string elst_damping = options_.get_str("EFP_ELST_DAMPING");
       std::string disp_damping = options_.get_str("EFP_DISP_DAMPING");

       if (elst_damping == "SCREEN") {
               opts.elec_damp = EFP_ELEC_DAMP_SCREEN;
       } else if (elst_damping == "OVERLAP") {
               opts.elec_damp = EFP_ELEC_DAMP_OVERLAP;
       } else if (elst_damping == "OFF") {
               opts.elec_damp = EFP_ELEC_DAMP_OFF;
       }

       if (disp_damping == "TT") {
           opts.disp_damp = EFP_DISP_DAMP_TT;
       } else if (disp_damping == "OVERLAP") {
           opts.disp_damp = EFP_DISP_DAMP_OVERLAP;
       } else if (disp_damping == "OFF") {
           opts.disp_damp = EFP_DISP_DAMP_OFF;
       }

       nfrag = options_["FRAGS"].size();
       molecule = Process::environment.molecule();

       frag_name = new char*[nfrag + 1]; /* NULL-terminated */

       for(int i=0; i<nfrag; i++) {
               std::string name = options_["FRAGS"][i].to_string();
               frag_name[i] = new char[name.length() + 1];
               strcpy(frag_name[i], name.c_str());
               for (char *p = frag_name[i]; *p; p++)
                       *p = tolower(*p);
       }
       frag_name[nfrag] = NULL;

       std::string psi_data_dir = Process::environment("PSIDATADIR");
       std::string frag_lib_path = psi_data_dir + "/fraglib";

       char **potential_file_list = make_potential_file_list((const char **)frag_name,
               frag_lib_path.c_str(), frag_lib_path.c_str());

       if ((res = efp_init(&efp, &opts, NULL, (const char **)potential_file_list, (const char **)frag_name))) {
               fprintf(outfile, efp_result_to_string(res));
                throw PsiException("efp",__FILE__,__LINE__);
       }

       for (int i = 0; potential_file_list[i]; i++)
                       delete[] potential_file_list[i];
       delete[] potential_file_list;

       fprintf(outfile, "  ==> Calculation Information <==\n\n");

       fprintf(outfile, "  Electrostatics damping: %12s\n", elst_damping.c_str());

       if (disp_enabled)
               fprintf(outfile, "  Dispersion damping:     %12s\n", disp_damping.c_str());

       fprintf(outfile, "\n");
}

void EFP::SetGeometry(){
       enum efp_result res;
       fprintf(outfile, "\n\n");
       fprintf(outfile, efp_banner());
       fprintf(outfile, "\n\n");

       fprintf(outfile, "  ==> Geometry <==\n\n");

       double *coords = NULL;

       molecule->print();
       if (molecule->nfragments() != nfrag)
               throw InputException("Molecule doesn't have FRAGS number of fragments.", "FRAGS", nfrag, __FILE__, __LINE__);

       // array of coordinates, 9 numbers for each fragment - first three atoms
       coords = new double[9 * nfrag];
        double * pcoords = coords;

       for (int i = 0; i < nfrag; i++) {
               std::vector<int> realsA;
               realsA.push_back(i);
               std::vector<int> ghostsA;
               for (int j = 0; j < nfrag; j++) {
                               if (i != j)
                                               ghostsA.push_back(j);
               }
               boost::shared_ptr<Molecule> monomerA = molecule->extract_subsets(realsA, ghostsA);
               monomerA->print();
               monomerA->print_in_bohr();

               int natomA = 0;
               for (int n=0; n<monomerA->natom(); n++)
                       if (monomerA->Z(n))
                               natomA++;

               if (natomA != 3)
                       throw InputException("Fragment doesn't have three coordinate triples.", "natomA", natomA, __FILE__, __LINE__);
       
               SharedMatrix xyz = SharedMatrix (new Matrix("Fragment Cartesian Coordinates(x,y,z)", monomerA->natom(), 3));
               double** xyzp = xyz->pointer();
       
               for (int j = 0; j < monomerA->natom(); j++) {
                       if (monomerA->Z(j)) {
                               *pcoords++ = xyzp[j][0] = monomerA->x(j);
                               *pcoords++ = xyzp[j][1] = monomerA->y(j);
                               *pcoords++ = xyzp[j][2] = monomerA->z(j);
                       }
               }

        fprintf(outfile, "%s\n", frag_name[i]);
       xyz->print();
       }

       if ((res = efp_set_coordinates_2(efp, coords))) {
               fprintf(outfile, efp_result_to_string(res));
                throw PsiException("efp",__FILE__,__LINE__);
       }
} // end of SetGeometry()

// compute efp energy components and/or gradient
void EFP::Compute() {
       enum efp_result res;
        double *grad = NULL;

       /* Main EFP computation routine */
       if ((res = efp_compute(efp, do_grad ? 1 : 0))) {
               fprintf(outfile, efp_result_to_string(res));
                throw PsiException("efp",__FILE__,__LINE__);
       }

       struct efp_energy energy;

       if ((res = efp_get_energy(efp, &energy))) {
               fprintf(outfile, efp_result_to_string(res));
                throw PsiException("efp",__FILE__,__LINE__);
       }

       if (do_grad) {
                       grad = new double[6 * nfrag];
                       if ((res = efp_get_gradient(efp, 6 * nfrag, grad))) {
                               fprintf(outfile, efp_result_to_string(res));
                                throw PsiException("efp",__FILE__,__LINE__);
                       }

                       fprintf(outfile, "  ==> EFP Gradient <==\n\n");

                       double *pgrad = grad;
                       for (int i = 0; i < nfrag; i++) {
                                       for (int j = 0; j < 6; j++) {
                                                       fprintf(outfile, "%14.6lf", *pgrad++);
                                       }
                                       fprintf(outfile, "\n");
                       }
                       fprintf(outfile, "\n");

                       SharedMatrix smgrad(new Matrix("EFP Gradient", nfrag, 6)); //
                       double ** psmgrad = smgrad->pointer();
                       pgrad = grad;
                       for (int i = 0; i < nfrag; i++) {
                                       for (int jj = 0; jj < 6; jj++) {
                                               psmgrad[i][jj] = *pgrad++;
                                       }
                       }

                       psi::Process::environment.set_gradient(smgrad);
                       smgrad->print();
       }

    fprintf(outfile, "  ==> Energetics <==\n\n");

    fprintf(outfile, "  Electrostatics Energy = %24.16f [H] %s\n", energy.electrostatic, elst_enabled ? "*" : "");
    fprintf(outfile, "  Polarization Energy =   %24.16f [H] %s\n", energy.polarization, pol_enabled ? "*" : "");
    fprintf(outfile, "  Dispersion Energy =     %24.16f [H] %s\n", energy.dispersion, disp_enabled ? "*" : "");
    fprintf(outfile, "  Exchange Energy =       %24.16f [H] %s\n", energy.exchange_repulsion, exch_enabled ? "*" : "");
    fprintf(outfile, "  Total Energy =          %24.16f [H] %s\n", energy.total, "*");

    Process::environment.globals["EFP ELST ENERGY"] = energy.electrostatic;
    Process::environment.globals["EFP POL ENERGY"] = energy.polarization;
    Process::environment.globals["EFP DISP ENERGY"] = energy.dispersion;
    Process::environment.globals["EFP EXCH ENERGY"] = energy.exchange_repulsion;
    Process::environment.globals["CURRENT ENERGY"] = energy.total;


}

}} // End namespaces

