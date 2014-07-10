#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/stats.c,v 1.4 1995-02-24 02:17:54 d3h325 Exp $ */

#include <stdio.h>

#include "sndrcvP.h"
#include "sndrcv.h"

void STATS_()
/*
  Print out communication statistics for this node
*/
{
  long me = NODEID_();
  long nproc = NNODES_();
  long i, msg_s, kb_s, s_s, r_s, msg_r, kb_r, s_r, r_r;


  (void) printf("Communication statistics for node %ld\n",me);
  (void) printf("-------------------------------------\n\n");
  (void) printf("\
                 sending                     receiving\n\
        -------------------------   -------------------------\n\
 node   #msgs.   #kb   secs  kb/s   #msgs.   #kb   secs  kb/s\n\
 ----   ------  -----  ----  ----   ------  -----  ----  ----\n");

  for (i=0; i<nproc; i++)
    if (i != me) {
      msg_s = SR_proc_info[i].n_snd;
      kb_s = SR_proc_info[i].nb_snd/1024.0;
      s_s = SR_proc_info[i].t_snd;
      r_s = (s_s > 0) ? kb_s / s_s : 0;
      msg_r = SR_proc_info[i].n_rcv;
      kb_r = SR_proc_info[i].nb_rcv/1024.0;
      s_r = SR_proc_info[i].t_rcv;
      r_r = (s_r > 0) ? kb_r / s_r : 0;
      
      (void) printf("%5ld%9ld%7ld%6ld%6ld%9ld%7ld%6ld%6ld\n", i,
		    msg_s, kb_s, s_s, r_s, 
		    msg_r, kb_r, s_r, r_r); 
    }
}

