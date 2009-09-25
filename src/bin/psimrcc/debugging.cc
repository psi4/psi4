/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include <liboptions/liboptions.hpp>
#include "debugging.h"

namespace psi{ namespace psimrcc{

Debugging::Debugging()
{
  level = new bool[11];
  for(int n=0;n<=10;n++)
    level[n]=false;

  int maxn = options_get_int("DEBUG");
  if(maxn>10)
    maxn = 10;
  for(int n=0;n<=maxn;n++)
    level[n]=true;
}
}} /* End Namespaces */
