#ifndef _psi_src_bin_psimrcc_debugging_h
#define _psi_src_bin_psimrcc_debugging_h
/***************************************************************************
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***********
****************************************************************/
#include <liboptions/liboptions.h>

namespace psi{ namespace psimrcc{

// #define NODEBUG
#ifdef NODEBUG
#define DEBUGGING(level,statements)
#else
#define DEBUGGING(level,statements) \
  if(debugging->is_level(level)){   \
    statements                      \
  }
#endif 

#ifdef NODEBUG
#define START_TIMER(level,title)
#else
#define START_TIMER(level,title)             \
  Timer timer;                               \
  if(debugging->is_level(level)){            \
    fprintf(outfile,"\n  %-48s ...",title);  \
    fflush(outfile);                         \
  }
#endif

#ifdef NODEBUG
#define END_TIMER(level)
#else
#define END_TIMER(level)                                    \
  if(debugging->is_level(level)){                           \
    fprintf(outfile," done. Timing %10.4f s",timer.get());  \
    fflush(outfile);                                        \
  }
#endif

/**
	@author Francesco Evangelista <frank@ccc.uga.edu>
*/
class Debugging{
public:
  Debugging(Options &options);
  ~Debugging() {delete[] level;}

  bool is_level(int n) {return(n <= 10 ? level[n] : false);}
private:
  Options &options_;
  bool* level;
};

extern Debugging* debugging;

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_debugging_h
