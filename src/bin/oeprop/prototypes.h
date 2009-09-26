/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/

#ifndef _psi_bin_oeprop_prototypes_h_
#define _psi_bin_oeprop_prototypes_h_

namespace psi { namespace oeprop {

void start_io(int argc, char *argv[]);
void stop_io();
//void punt(const char *);
void parsing();
void grid_unitvec();
void compute_density();
void read_density();
void get_nmo();
void compute_overlap();
void compute_oeprops();
void compute_grid();
void read_basset_info();
void read_zvec();
void MI_OSrecurs(double, double, double, double, double, double, double, int, int, int);
void AI_OSrecurs(double, double, double, double, double, double, double, double, double, double, int, int);
void populate();
void compute_mp_ref_xyz();
void move2ref();
void init_xyz();
double ***init_box(int, int, int);
void free_box(double ***, int, int);
void print_intro();
void print_tasks();
void print_params();
void print_pop_header();
void print_mp();
void print_lm();
void print_esp();
void print_grid();
void print_misc();
void compute_connectivity();
void get_opdm_lbl(void);
void grid_parse();

}} // namespace psi::oeprop

#endif // header guard
