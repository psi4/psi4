#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "xregion.h"

void Setcmap()
/*
  Make the color map ... 
*/
{
  int i;
  XColor color, spare;
  double cscale = 1.0 / ((double) (MAX_COL-1));
  double hue, saturation, value;
  double redvar,red=0xba/255.0, green, blue=0xd2/255.0;
  redvar = red;

  colormap = DefaultColormap(display, screen);

  /* Linear interpolation on green */

  for (i=0; i<MAX_COL; i++) 
  {
    if (i == 0) 
    {
      /* Assign white as the first color 
      color.red   = 65535;
      color.blue  = 65535;
      color.green = 65535; */
      color.pixel = CANVAS_COLOR;
      XQueryColor(display, colormap, &color);
    }
    else 
    {
      if(i>=(MAX_COL-5))
      {
        redvar = .4*red + .6*red*(MAX_COL-i)/5.; 
      }
      color.red   = (short) (redvar   * 65535.0);
      color.blue  = (short) (blue  * 65535.0);
      green = cscale * (MAX_COL - i);
      color.green = (short) (green * 65535.0);
    }

    if (XAllocColor(display, colormap, &color) == 0)
    {
      Error("couldn't assign color",i);
    }

    cmap[i] = color.pixel;
 
/* now set the "accessed" color for regions that were accessed */
    color.red   = 65535;
    color.green = 65535;
    color.blue  = 65535 * (220.0/255.0);
    if (XAllocColor(display, colormap, &color) == 0)
    {
      Error("couldn't assign accessed color",99);
    }
    cmap[MAX_COL] = color.pixel; 

/*
    (void) printf("Colour %d red=%x, green=%x, blue=%x, pixel=%x\n",
		  i, color.red, color.green, color.blue, color.pixel);
*/
  }

}
