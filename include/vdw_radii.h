/*------------------------------------------
  standard van der Waals radii for atoms
  (in angstrom)

  EV, March 30, 2001
 ------------------------------------------*/

#ifndef _psi_include_vdwradii_h_
#define _psi_include_vdwradii_h_

#define LAST_VDW_RADII_INDEX 9

#ifdef __cplusplus
extern "C" {
#endif

double atomic_vdw_radii[] = { 2.0,  /* default element or ghost */
                    1.2,  /* hydrogen */
		    1.4,  /* helium */
		    1.82, /* lithium */
		    1.8,  /* berrilium (there's no info in literature) */
		    1.6,  /* boron (there's no info in literature) */
		    1.70, /* carbon */
		    1.55, /* nitrogen */
		    1.52, /* oxygen */
		    1.47  /* fluorine */
};

#ifdef __cplusplus
}
#endif

#endif /* header guard */
