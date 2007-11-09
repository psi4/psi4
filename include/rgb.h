/*------------------------------------------
  Default RGB colors for atoms

  EV, March 30, 2001
 ------------------------------------------*/

#ifndef _psi_include_rgb_h_
#define _psi_include_rgb_h_

#define LAST_RGB_INDEX 9

#ifdef __cplusplus
extern "C" {
#endif

double atomic_rgb[][3] = { {0.40, 0.40, 0.40},  /* default element or ghost */
                    {1.00, 1.00, 1.00},  /* hydrogen */
		    {0.80, 0.80, 0.50},  /* helium */
		    {0.30, 0.80, 0.30},  /* lithium */
		    {0.65, 0.80, 0.25},  /* berrilium */
		    {0.50, 0.80, 0.15},  /* boron */
		    {0.25, 0.25, 0.25},  /* carbon */
		    {0.00, 0.00, 1.00},  /* nitrogen */
		    {1.00, 0.00, 0.00},  /* oxygen */
		    {0.75, 0.40, 0.15}  /* fluorine */
};

#ifdef __cplusplus
}
#endif

#endif /* header guard */
