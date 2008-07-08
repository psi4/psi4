#include <masses.h>

namespace psi { namespace optking {

extern double **build_B(const Salc_set &, const Intcos &);
extern double **build_G(double **B, int nrows, int ncols, bool use_masses,
    double *masses = NULL);

class Geom_opt {
  int step;
  double **geoms;
  double **q;
  double **fq;
  int B_present;
  double **B;
  double *masses;

 public:
  Intcos * simples;
  Salc_set * salcs;

  Geom_opt();
  ~Geom_opt();
  void build_B(void);
  void print_B(void);
 
  void load_status();
  int opt_step();
  int newgeom();
  void save_status();
};

Geom_opt::Geom_opt(void)
{
  B_present = 0;
  printf("geom_opt created");

  simples = new Intcos();
  simples->add_intcos_from_input();
}


int Geom_opt::opt_step(void)
{
  int nints = salcs->size();
  build_B();

 // double **G = build_G(B, nints, 3*salcs.get_natom(), 1, simples.masses); 
 // mat_print(G,nints,3*salcs.get_natom(),outfile);
  return 1;
}

void Geom_opt::build_B(void) {
  simples->build_Bsimp();
  if (B_present) free_block(B);
  int ncols = 3*salcs->get_natom();
  int nrows = salcs->size();
  int nlinks = simples->size();
  B = block_matrix(nrows,ncols);
  mmult(salcs->get_matrix(),0,simples->get_Bsimp(),0,B,0,nrows,nlinks,ncols,0);
  return;
}

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

Geom_opt::~Geom_opt(void) {
  if (B_present) free_block(B);
}

}}

