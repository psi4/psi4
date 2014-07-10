#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "xregion.h"

void ReadEventFile(char *filename)
{
  FILE *fin;
  long int i, k, act_events = 0;
  char errmsg[128];
 
  fin = fopen(filename,"r");

  if (!fin)
  {
    Error("Input File Not Found",-2);
  }

  if (!fscanf(fin, "%d%d%d", &grid_y, &grid_x, &num_events))
  {
    sprintf(errmsg, "Unable to read data from file %s", filename);
    Error(errmsg, -2);
  }

  if (!fscanf(fin, "%d%d", &rows, &cols))
  {
    sprintf(errmsg, "Unable to read data from file %s", filename);
    Error(errmsg, -2);
  }

  if (!(overlay_row = (int*) malloc(rows * sizeof(int))))
  {
    Error("couldn't allocate memory",-1);
  }

  if (!(overlay_col = (int*) malloc(cols * sizeof(int))))
  {
    Error("couldn't allocate memory",-1);
  }
  
  for(i = 0; i < rows; i++)
  {
    if (!fscanf(fin, "%d", (overlay_row + i)))
    {
      sprintf(errmsg, "Unable to read data from file %s", filename);
      Error(errmsg, -2);
    }
  }

  for(i = 0; i < cols; i++)
  {
    if (!fscanf(fin, "%d", (overlay_col + i)))
    {
      sprintf(errmsg, "Unable to read data from file %s", filename);
      Error(errmsg, -2);
    }
  }

  if (!(record = (int*) malloc(RECLEN * num_events * sizeof(int))))
  {
    Error("couldn't allocate memory",-1);
  }

  if (!(ev_times = (unsigned long int *) malloc(num_events * sizeof(unsigned long))))
  {
    Error("couldn't allocate memory",-2);
  }

  for(i = 0; i < num_events; i++)
  {
    for(k = 0; k < RECLEN; k++)
    {
      fscanf(fin, "%d", (i * RECLEN + k) + record);
    }
    if(fscanf(fin, "%lu", ev_times + i))
    {
      act_events++;
    }
    
    /* Adjust from Fortran to C base addressing */
    
    for (k = 2; k <= 5; k++)
    {
      (*((i * RECLEN + k) + record))--;
    }

    if (feof(fin))
    {
      break;
    }
  }

  num_events = act_events;
  printf("File %s has been read. %d events are displayed\n",filename,num_events);

  fclose(fin);
}
