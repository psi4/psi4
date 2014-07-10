/** @file
  This header file provides definitions for c for the names of the
  c message passing routines accessible from FORTRAN. It need not
  be included directly in user c code, assuming that sndrcv.h has already.

  It is needed as the FORTRAN naming convention varies between machines
  and it is the FORTRAN interface that is portable, not the c interface.
  However by coding with the macro defnition names c portability is
  ensured.
*/
#ifndef SRFTOC_H_
#define SRFTOC_H_

#define BRDCST_    armci_tcgmsg_brdcst
#define DGOP_      armci_tcgmsg_dgop
#define DRAND48_   armci_tcgmsg_drand48
#define IGOP_      armci_tcgmsg_igop
#define LLOG_      armci_tcgmsg_llog
#define MDTOB_     armci_tcgmsg_mdtob
#define MDTOI_     armci_tcgmsg_mdtoi
#define MITOB_     armci_tcgmsg_mitob
#define MITOD_     armci_tcgmsg_mitod
#define MTIME_     armci_tcgmsg_mtime
#define NICEFTN_   armci_tcgmsg_niceftn
#define NNODES_    armci_tcgmsg_nnodes
#define NODEID_    armci_tcgmsg_nodeid
#define NXTVAL_    armci_tcgmsg_nxtval
#define PARERR_    armci_tcgmsg_parerr
#define PBEGINF_   armci_tcgmsg_pbeginf
#define PBGINF_    armci_tcgmsg_pbginf
#define PEND_      armci_tcgmsg_pend
#define PFCOPY_    armci_tcgmsg_pfcopy
#define PFILECOPY_ armci_tcgmsg_pfilecopy
#define PROBE_     armci_tcgmsg_probe
#define RCV_       armci_tcgmsg_rcv
#define SETDBG_    armci_tcgmsg_setdbg
#define SND_       armci_tcgmsg_snd
#define SRAND48_   armci_tcgmsg_srand48
#define STATS_     armci_tcgmsg_stats
#define SYNCH_     armci_tcgmsg_synch
#define TCGREADY_  armci_tcgmsg_tcgready
#define TCGTIME_   armci_tcgmsg_tcgtime
#define WAITCOM_   armci_tcgmsg_waitcom

#endif /* SRFTOC_H_ */
