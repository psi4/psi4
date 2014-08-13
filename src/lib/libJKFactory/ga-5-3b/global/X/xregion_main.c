#if HAVE_CONFIG_H
#   include "config.h"
#endif

#define FIRST_TIME
#include "xregion.h"
#if HAVE_MATH_H
#   include <math.h>
#endif

int main(int argc, char **argv)
{
  int i, argn;
  char filename[128];

  /* initialize global variables */
  first_time = True;  /* Used to set scroll bar on first expose */
  interval_max = 2000;
  interval = 500;         /* 0.5s between exposures by default */
  slowdown_max = 100.;
  slowdown_min = .1;
  slowdown = INITSLOW;        /* slowdown factor for animation */
  oldslowdown = INITSLOW;
  cur_time = 0; /* current time */
  maxval = 0.;       /* max value of integral, zero is default */
  working = False; 
  animation = True;
  cur_event = 0;

  /* First read the argument list */ 
  for(i = 1, argn = 0; i < argc; i++)
  {
    if (argv[i][0] == '-') 
    {
      break;
    }
    else 
    {
      argn = i; 
    }
  }
  argn ++;

  if (argn < 2 )
  {
    printf("Usage:\n");
    printf("xregion <filename> \n");
    exit(1);
  }

  sscanf(argv[1],"%s", filename);

  ReadEventFile(filename);
  
  scale = 1;
  overview_scale = 499.0 / GA_MAX((grid_x), (grid_y));
  overview_width = ceil( (grid_x) * (500.0 / GA_MAX((grid_x), (grid_y))) );
  overview_height = ceil( (grid_y) * (500.0 / GA_MAX((grid_x), (grid_y))) );

  printf("overview_scale %lf\n", overview_scale);
  printf("overview_width %d\n", overview_width);
  printf("overview_height %d\n", overview_height);

  /* Realize everything */

  xregion_app = create_overview(argc, argv);
  XtRealizeWidget(overview_shell);

  /* Enter the event loop */

  XtAppMainLoop(xregion_app);

}

