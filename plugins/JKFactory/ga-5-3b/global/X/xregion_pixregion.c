#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "xregion.h"

void UpdatePixRegion(int ilo, int ihi, int jlo, int jhi, int increment, 
                    double time)
{
  register int i, j, k, l, index;
  register u_char *from, *to, *tempk, *tempkl, value, *pflag;
  register double *pintegr, *pltime, corr;
  
  for (i = ilo; i < ihi + 1; i++)
  {
    for (j = jlo; j < jhi + 1; j++) 
    {
      to = pict + (i * pict_width * scale + j) * scale; 
      from = grid + (i * pict_width + j);
      pflag = flag + (i * pict_width + j);
      pintegr = integr + (i * pict_width + j);
      pltime = ltime  + (i * pict_width + j);

      /* increment == 0 means animation is done and displaying integrals */

      if (animation)
      {
        corr = *from > 1 ? *from - 1 : 0;
        corr = corr  * (time - *pltime);

        if(corr < 0.0)
        {
          fprintf(stderr, "error: time =%f ltime =%f height =%d \n",
                  time,*pltime,*from);
        }

        *pintegr += corr;

        *pltime   = time;

        *from = *from + increment;
        index = *from;
        index = GA_MIN(index,MAX_COL-1);

        if(increment)
        {
          *pflag = 1;
        }

        /*  calculate max value of integrals */
        if(*pintegr > maxval)
        {
          maxval = *pintegr;
        }
      }
      else /* done with animation display integral */
      {
        index = (int) (((*pintegr) / maxval) * MAX_COL);
        index = GA_MIN(index, MAX_COL - 1);
        if(!index && *pflag)
        {
          index = MAX_COL;  /* sets the "accessed" color */
        }
      }
          
      value = cmap[index];

      for (k = 0, tempk = to; k < scale; k++, tempk += pict_width * scale)
      {
	for (l = 0, tempkl = tempk; l < scale; l++, tempkl++)
        {
	  *tempkl = value;
        }
      }

    }  /* end for j */
  } /* end for i */
}


void DisplayPixRegion(int ilo, int ihi, int jlo, int jhi)
{
  int count, x, y, height, width;

  y = ilo * scale;
  x = jlo * scale;
  height = (ihi - ilo + 1) * scale;
  width = (jhi - jlo + 1) * scale;
  
  XPutImage(display, window, gc, image, x, y, x, y, width, height);
  XFlush(display);

/*
  for(count = 0; count < rows; count++)
  {
    XDrawLine(display, XtWindow(canvas_widget), gc, 0, 
             (overlay_row[count] - 1 - top_edge) * scale,
              pict_width * scale - 1, 
             (overlay_row[count] - 1 - top_edge) * scale);
  }
  for(count = 0; count < cols; count++)
  {
    XDrawLine(display, XtWindow(canvas_widget), gc, 
             (overlay_col[count] - 1 - left_edge) * scale, 0, 
             (overlay_col[count] - 1 - left_edge) * scale, 
              pict_height * scale - 1);
  }
*/

  for(count = 0; count < rows; count++)
  {
    XDrawLine(display, XtWindow(canvas_widget), gc, 0,
             (overlay_row[count]  - top_edge + .5) * scale,
              pict_width * scale - 1,
             (overlay_row[count] - top_edge + .5) * scale);
  }
  for(count = 0; count < cols; count++)
  {
    XDrawLine(display, XtWindow(canvas_widget), gc,
             (overlay_col[count]  - left_edge + .5) * scale, 0,
             (overlay_col[count]  - left_edge + .5) * scale,
              pict_height * scale - 1);
  }


}

