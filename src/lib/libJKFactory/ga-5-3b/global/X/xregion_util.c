#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "xregion.h"

void Error(char *message, int err_num)
{
  (void) fflush(stdout);
  (void) fprintf(stderr,"\n\nError was called.\n");
  (void) fprintf(stderr,message);
  (void) fprintf(stderr," %d (%#x).\n",err_num,err_num);
  exit(1);
}

int in_display_region(int ilo, int ihi, int jlo, int jhi)
{
  return( !((ilo > bottom_edge) || 
            (ihi < top_edge) || 
            (jlo > right_edge) || 
            (jhi < left_edge)));

}

void set_config()
{
  int found = True;
  Pixel pix_val;
  Display *disp;
  XColor xcolor;
  XColor spare;
  Colormap defcmap;
  char *bufptr, name[40], buf[80];
  char path[128];
  FILE *fd;

  disp = XtDisplay(overview_shell);
  defcmap = DefaultColormapOfScreen(XtScreen(overview_shell));

  if (XAllocNamedColor(disp, defcmap, "black", &xcolor, &spare) != 0)
  {
    DEFAULT_FG = xcolor.pixel; 
  }
  else
  {
    fprintf(stderr, "unable to allocate color 'black'\n");
  }

  if (XAllocNamedColor(disp, defcmap, "whitesmoke", &xcolor, &spare) != 0)
  {
    DEFAULT_BG = xcolor.pixel; 
  }
  else
  {
    fprintf(stderr, "unable to allocate color 'whitesmoke'\n");
  }

  if (XAllocNamedColor(disp, defcmap, "green", &xcolor, &spare) != 0)
  {
    SELECT_COLOR = xcolor.pixel; 
    RBAND_COLOR = xcolor.pixel; 
  }
  else
  {
    fprintf(stderr, "unable to allocate color 'green'\n");
  }

  if (XAllocNamedColor(disp, defcmap, "black", &xcolor, &spare) != 0)
  {
    CANVAS_COLOR = xcolor.pixel; 
  }
  else
  {
    fprintf(stderr, "unable to allocate color 'black'\n");
  }

  if (XAllocNamedColor(disp, defcmap, "red", &xcolor, &spare) != 0)
  {
    GRID_COLOR = xcolor.pixel; 
  }
  else
  {
    fprintf(stderr, "unable to allocate color 'red'\n");
  }


  sprintf(path, "./xregion.config");

  if ((fd = fopen(path,"r")) == NULL)
  {
    fprintf(stderr,"config file '%s' not found, using default colors\n", path);
    found = False;
  }
 
  while (found && fgets(buf, sizeof(buf), fd) != NULL)
  {
    if (buf[0] == '#') continue; /* a comment */
    
    bufptr = strtok(buf, " ");
    if (bufptr)
    {
      strcpy(name, bufptr);
    } 

    bufptr = strtok(NULL, "\n");
    if (bufptr)
    {
      if (XAllocNamedColor(disp, defcmap, bufptr, &xcolor, &spare) == 0)
      {
        strcat(name, " Not Allocated"); 
      }
      else
      {
        pix_val = xcolor.pixel;
      } 
    } 
   
    if (strstr(name, "Not"))
    {
      fprintf(stderr, "%s: check color name '%s'\n", name, bufptr);
    }
    else if (strstr(name, "foreground"))
    {
      DEFAULT_FG = pix_val; 
    }
    else if (strstr(name, "background")) 
    {
      DEFAULT_BG = pix_val; 
    } 
    else if (strstr(name, "selection")) 
    {
      SELECT_COLOR = pix_val; 
    } 
    else if (strstr(name, "rubberband")) 
    {
      RBAND_COLOR = pix_val; 
    }
    else if (strstr(name, "grid")) 
    {
      GRID_COLOR = pix_val; 
    }
    else if (strstr(name, "canvas")) 
    {
      CANVAS_COLOR = pix_val; 
    }
    else
    {
       fprintf(stderr, "Unknown entry: '%s'\n", name);
    }

  }

  fclose(fd);
}
