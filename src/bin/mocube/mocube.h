/*! \file
    \ingroup MOCUBE
    \brief Enter brief description of file here 
*/

#ifndef _psi_bin_mocube_mocube_h_
#define _psi_bin_mocube_modube_h_

namespace psi { namespace mocube {

struct cube_struct {
  char title[80];
  char subtitle[80];
  int natom;
  int *mos_to_plot;
  int nmo_to_plot;
  double grid_start[3]; /* x, y and z */
  double grid_end[3]; /* x, y and z */
  double border;
  int ngrid[3];         /* x, y and z */
  double step_size[3];   /* x, y and z */
  double *zvals;
  double **geom;
  double ****grid;
} ;

struct Params {
  char *label;
  int nao;
  int natom;
  int nirreps;
  int nirreps_present;
  int nprim;
  int nso;
  int nmo;
  int nclsd;
  int *clsdpi;
  int *openpi;
  int *orbspi;
  int homo;
  int lumo;
} ;

#ifdef EXTERN
#undef EXTERN
#define EXTERN extern
#else
#define EXTERN
#endif

EXTERN struct Params params;
EXTERN struct cube_struct cube;

}} // namespace psi::mocube

#endif // header guard
