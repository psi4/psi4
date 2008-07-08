
#include "optking.h"
#include "intcos.h"
#include "salc.h"
#include "except.h"
#include "geom_opt.h"

extern "C" { char *gprgid() { char *prgid = "OPT08"; return(prgid); } }

namespace psi { namespace optking {
  void intro(void);       // print header
  void get_optinfo(void); // get parameters
}}

int main(int argc, char **argv) {
  using namespace psi::optking;
  int parsed = 1;
  FILE *fp_intco;

  psi_start(&infile,&outfile,&psi_file_prefix,argc-parsed,argv+parsed,0);

  ffile_noexit(&fp_intco, "intco.dat", 2);
  if (fp_intco) {
    ip_append(fp_intco, outfile) ;
    fclose(fp_intco);
  }

  string progid(":");
  progid += gprgid();
  ip_cwk_add(&progid[0]);
  psio_init();
  psio_ipv1_config();

  intro();
  get_optinfo();

  Geom_opt geom_opt;

  geom_opt.build_B();

/*
 Intcos intcos;
 intcos.add_intcos_from_input();
 intcos.compute();
 intcos.compute_s();
 intcos.print(1);
 intcos.print_s();
 intcos.add_intcos_from_input();
 intcos.print(1);
 intcos.build_B();
 intcos.print_B();

 if (ip_exist("SYMM",0)) {
   try { Salc_set salc_set("SYMM",0,optinfo.natom,intcos.size()); }
   catch (bad_intco_io & bi) {
     bi.mesg();
     fprintf(outfile,"Sorry: optking.cc can't read SYMM\n");
   }
 }

 Salc_set salc_set("SYMM",0,optinfo.natom,intcos.size());
 salc_set.print();
 salc_set.build_matrix();
 salc_set.print_matrix();

 double **B = build_B(salc_set,intcos);
 print_mat(B,salc_set.size(),3*optinfo.natom,outfile);
 double **G = build_G(B, salc_set.size(), 3*optinfo.natom, 0) {
*/

 psi_stop(infile,outfile,psi_file_prefix);

  return(0);
}

namespace psi { namespace optking {

/***  INTRO   prints into ***/
void intro() {
  fprintf(outfile, "\n\t----------------------------------\n");
  fprintf(outfile, "\t  OPT08: for we're not sure yet   \n");
  fprintf(outfile, "\t       - R. A. King               \n");
  fprintf(outfile, "\t----------------------------------\n");
}

}}

