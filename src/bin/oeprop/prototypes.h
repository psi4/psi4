/*! \file
    \ingroup OEPROP
    \brief Enter brief description of file here 
*/

#ifndef _psi_bin_oeprop_prototypes_h_
#define _psi_bin_oeprop_prototypes_h_

namespace psi { namespace oeprop {

void parsing(Options &);
void compute_density(void);
void read_density(Options &);
void get_nmo(void);
void compute_overlap(void);
void compute_oeprops(void);
void compute_grid(void);
void read_basset_info(void);
void read_zvec(void);
void MI_OSrecurs(double, double, double, double, double, double, double, int, int, int);
void AI_OSrecurs(double, double, double, double, double, double, double, double, double, double, int, int);
void populate(void);
void compute_mp_ref_xyz(void);
void move2ref(void);
void init_xyz(void);
double ***init_box(int, int, int);
void free_box(double ***, int, int);
void print_intro(void);
void print_tasks(void);
void print_params(void);
void print_pop_header(void);
void print_mp(void);
void print_lm(void);
void print_esp(void);
void print_grid(void);
void print_misc(void);
void compute_connectivity(void);
void get_opdm_lbl(Options &);
void grid_parse(Options &);

}} // namespace psi::oeprop

#endif // header guard
