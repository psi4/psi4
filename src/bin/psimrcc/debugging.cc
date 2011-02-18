/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <liboptions/liboptions.h>
#include "debugging.h"

namespace psi{ namespace psimrcc{

        Debugging::Debugging(Options &options):
                options_(options)
{
  level = new bool[11];
  for(int n=0;n<=10;n++)
    level[n]=false;

  int maxn = options_.get_int("DEBUG");
  if(maxn>10)
    maxn = 10;
  for(int n=0;n<=maxn;n++)
    level[n]=true;
}
}} /* End Namespaces */
