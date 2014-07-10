#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*  This module contains the functions for displaying the access color map
    that appears on the right side of the display window 
*/

#include "xregion.h"

void DrawColorMap()
{
/* Actually fill in the colors */
int i, black = 1;
unsigned width = 20, height = 350 / MAX_COL;
int x, y, length,index;
XColor color;

  for (i=0,y=0,x=0;i<MAX_COL;i++)
  {
    index = (animation==False && i==0)? MAX_COL : i; 
    XSetForeground(display, gc_map, (unsigned long) cmap[index]);
    XFillRectangle(display,window_map,gc_map,x,y,width,height);
    y += height;
    /* XSetForeground(display, gc_map, (unsigned long) cmap[index]);*/
  }
}


void PrintColorMapText()
{
/* Print Legend numbers */
int i;
unsigned width = 20, height = 350 / MAX_COL;
int x, y, length;
char string[9];
double factor, intval;
int formatE;

  XSetForeground(display, gc_map, DEFAULT_FG);

  x = width+4;
  y = 0;
  factor  = maxval/(MAX_COL);
  formatE = (maxval >= 10. || maxval < .001) ? 1 : 0; 
  for (i=0;i<MAX_COL;i++)
  {
    y += height;
    if(animation)
    {
      sprintf(string,"%2d",i);
    }
    else
    {
      intval = factor*i;
      if(formatE)
      {
        sprintf(string,"%5.1e",intval); 
      }
      else
      {
        sprintf(string,"%6.4f",intval);
      }
    }
    length = strlen(string);
    XDrawString(display,window_map,gc_map,x,y-2,string, length);
  }
}

