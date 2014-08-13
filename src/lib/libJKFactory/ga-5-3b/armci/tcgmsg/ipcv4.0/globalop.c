#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/globalop.c,v 1.8 2004-04-01 02:04:56 manoj Exp $ */
#include <stdlib.h>
#ifdef SEQUENT
#include <strings.h>
#else
#include <string.h>
#endif
#include "sndrcv.h"
#include "msgtypesc.h"

#define TCG_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define TCG_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define TCG_ABS(a) (((a) >= 0) ? (a) : (-(a)))

extern void free();

#ifndef IPSC
#include "sndrcvP.h"

#define GOP_BUF_SIZE 81920

static void idoop(n, op, x, work)
     long n;
     char *op;
     long *x, *work;
{
  if (strncmp(op,"+",1) == 0)
    while(n--)
      *x++ += *work++;
  else if (strncmp(op,"*",1) == 0)
    while(n--)
      *x++ *= *work++;
  else if (strncmp(op,"max",3) == 0)
    while(n--) {
      *x = TCG_MAX(*x, *work);
      x++; work++;
    }
  else if (strncmp(op,"min",3) == 0)
    while(n--) {
      *x = TCG_MIN(*x, *work);
      x++; work++;
    }
  else if (strncmp(op,"absmax",6) == 0)
    while(n--) {
      register long x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
      *x = TCG_MAX(x1, x2);
      x++; work++;
    }
  else if (strncmp(op,"absmin",6) == 0)
    while(n--) {
      register long x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
      *x = TCG_MIN(x1, x2);
      x++; work++;
    }
  else if (strncmp(op,"or",2) == 0) 
    while(n--) {
      *x |= *work;
      x++; work++;
    }
  else
    Error("idoop: unknown operation requested", (long) n);
}

static void ddoop(n, op, x, work)
     long n;
     char *op;
     double *x, *work;
{
  if (strncmp(op,"+",1) == 0)
    while(n--)
      *x++ += *work++;
  else if (strncmp(op,"*",1) == 0)
    while(n--)
      *x++ *= *work++;
  else if (strncmp(op,"max",3) == 0)
    while(n--) {
      *x = TCG_MAX(*x, *work);
      x++; work++;
    }
  else if (strncmp(op,"min",3) == 0)
    while(n--) {
      *x = TCG_MIN(*x, *work);
      x++; work++;
    }
  else if (strncmp(op,"absmax",6) == 0)
    while(n--) {
      register double x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
      *x = TCG_MAX(x1, x2);
      x++; work++;
    }
  else if (strncmp(op,"absmin",6) == 0)
    while(n--) {
      register double x1 = TCG_ABS(*x), x2 = TCG_ABS(*work);
      *x = TCG_MIN(x1, x2);
      x++; work++;
    }
  else
    Error("ddoop: unknown operation requested", (long) n);
}

/*ARGSUSED*/
void DGOP_(ptype, x, pn, op, len)
     double *x;
     long *ptype, *pn;
     char *op;
     int len;
/*
  Global summation optimized for networks of clusters of processes.

  This routine is directly callable from C only.  There is a
  wrapper that makes fortran work (see bottom of this file).
*/
{
  long me = NODEID_();
  long master = SR_clus_info[SR_clus_id].masterid;
  long nslave = SR_clus_info[SR_clus_id].nslave;
  long slaveid = me - master;
  long synch = 1;
  long type = (*ptype & MSGDBL) ? *ptype : *ptype + MSGDBL;
  long nleft = *pn;
  long buflen = TCG_MIN(nleft,GOP_BUF_SIZE); /* Try to get even sized buffers */
  long nbuf   = (nleft-1) / buflen + 1;
  long zero = 0;
  double *tmp = x;
  double *work;
  long nb, ndo, lenmes, from, up, left, right;

  buflen = (nleft-1) / nbuf + 1;
  if (!(work = (double *) malloc((unsigned) (buflen*sizeof(double)))))
     Error("DGOP: failed to malloc workspace", nleft);

  /* This loop for pipelining and to avoid caller
     having to provide workspace */

  while (nleft) {
    ndo = TCG_MIN(nleft, buflen);
    nb  = ndo * sizeof(double);

    /* Do summation amoung slaves in a cluster */

    up    = master + (slaveid-1)/2;
    left  = master + 2*slaveid + 1;
    right = master + 2*slaveid + 2;

    if (left < (master+nslave)) {
      RCV_(&type, (char *) work, &nb, &lenmes, &left, &from, &synch);
      ddoop(ndo, op, x, work);
    }
    if (right < (master+nslave)) {
      RCV_(&type, (char *) work, &nb, &lenmes, &right, &from, &synch);
      ddoop(ndo, op, x, work);
    }
    if (me != master)
      SND_(&type, (char *) x, &nb, &up, &synch);

    /* Do summation amoung masters */

    if (me == master) {
      up    = (SR_clus_id-1)/2;
      left  = 2*SR_clus_id + 1;
      right = 2*SR_clus_id + 2;
      up = SR_clus_info[up].masterid;
      left = (left < SR_n_clus) ? SR_clus_info[left].masterid : -1;
      right = (right < SR_n_clus) ? SR_clus_info[right].masterid : -1;

      if (left > 0) {
        RCV_(&type, (char *) work, &nb, &lenmes, &left, &from, &synch);
        ddoop(ndo, op, x, work);
      }
      if (right > 0) {
        RCV_(&type, (char *) work, &nb, &lenmes, &right, &from, &synch);
        ddoop(ndo, op, x, work);
      }
      if (me != 0)
        SND_(&type, (char *) x, &nb, &up, &synch);
    }
    nleft -= ndo;
    x     += ndo;
    type  += 13;   /* Temporary hack for hippi switch */
  }
  free((char *) work);

  /* Zero has the results ... broadcast them back */
  nb = *pn * sizeof(double);
  BRDCST_(&type, (char *) tmp, &nb, &zero);
}

void IGOP_(ptype, x, pn, op, len)
     long *x;
     long *ptype, *pn;
     char *op;
     int len;
/*
  Global summation optimized for networks of clusters of processes.

  This routine is directly callable from C only.  There is a
  wrapper that makes fortran work (see the bottom of this file).
*/
{
  long me = NODEID_();
  long master = SR_clus_info[SR_clus_id].masterid;
  long nslave = SR_clus_info[SR_clus_id].nslave;
  long slaveid = me - master;
  long synch = 1;
  long type = (*ptype & MSGINT) ? *ptype : *ptype + MSGINT;
  long nleft = *pn;
  long zero = 0;
  long *tmp = x;
  long *work;
  long nb, ndo, lenmes, from, up, left, right;

  if (!(work = (long *) 
	malloc((unsigned) (TCG_MIN(nleft,GOP_BUF_SIZE)*sizeof(long)))))
     Error("IGOP: failed to malloc workspace", nleft);

  /* This loop for pipelining and to avoid caller
     having to provide workspace */

  while (nleft) {
    ndo = TCG_MIN(nleft, GOP_BUF_SIZE);
    nb  = ndo * sizeof(long);
     /* Do summation amoung slaves in a cluster */

    up    = master + (slaveid-1)/2;
    left  = master + 2*slaveid + 1;
    right = master + 2*slaveid + 2;

    if (left < (master+nslave)) {
      RCV_(&type, (char *) work, &nb, &lenmes, &left, &from, &synch);
      idoop(ndo, op, x, work);
    }
    if (right < (master+nslave)) {
      RCV_(&type, (char *) work, &nb, &lenmes, &right, &from, &synch);
     idoop(ndo, op, x, work);
    }
    if (me != master)
      SND_(&type, (char *) x, &nb, &up, &synch);

    /* Do summation amoung masters */

    if (me == master) {
      up    = (SR_clus_id-1)/2;
      left  = 2*SR_clus_id + 1;
      right = 2*SR_clus_id + 2;
      up = SR_clus_info[up].masterid;
      left = (left < SR_n_clus) ? SR_clus_info[left].masterid : -1;
      right = (right < SR_n_clus) ? SR_clus_info[right].masterid : -1;

      if (left > 0) {
        RCV_(&type, (char *) work, &nb, &lenmes, &left, &from, &synch);
        idoop(ndo, op, x, work);
      }
      if (right > 0) {
        RCV_(&type, (char *) work, &nb, &lenmes, &right, &from, &synch);
        idoop(ndo, op, x, work);
      }
      if (me != 0)
        SND_(&type, (char *) x, &nb, &up, &synch);
    }
    nleft -= ndo;
    x     += ndo;
    type  += 13;   /* Temporary hack for hippi switch */
  }
  (void) free((char *) work);

  /* Zero has the results ... broadcast them back */
  nb = *pn * sizeof(long);
  BRDCST_(&type, (char *) tmp, &nb, &zero);
}

#endif

/* Wrapper for fortran interface ... UGH ... note that
   string comparisons above do NOT rely on NULL termination
   of the operation character string */

#ifdef CRAY
#include <fortran.h>
#endif
#ifdef ARDENT
struct char_desc {
  char *string;
  int len;
};
#endif

/*ARGSUSED*/
#if defined(CRAY) || defined(CRAY)
#ifdef ARDENT
void dgop_(ptype, x, pn, arg)
     long *ptype, *pn;
     double *x;
     struct char_desc *arg;
{
  char *op = arg->string;
  int len_op = arg->len;
#endif
#if defined(CRAY)
void dgop_(ptype, x, pn, arg)
     long *ptype, *pn;
     double *x;
     _fcd arg;
{
  char *op = _fcdtocp(arg);
  int len_op = _fcdlen(arg);
#endif
  DGOP_(ptype, x, pn, op);
}
#endif
/* This crap to handle FORTRAN character strings */

/*ARGSUSED*/
#if defined(CRAY) || defined(CRAY)
#ifdef ARDENT
void igop_(ptype, x, pn, arg)
     long *ptype, *pn;
     long *x;
     struct char_desc *arg;
{
  char *op = arg->string;
  int len_op = arg->len;
#endif
#if defined(CRAY)
void igop_(wrap_ptype, x, wrap_pn, arg)
     long *wrap_ptype, *wrap_pn;
     long *x;
     _fcd arg;
{
  long ptype, pn;
  ptype = (long) *ptype;
  char *op = _fcdtocp(arg);
  int len_op = _fcdlen(arg);
#endif
  IGOP_(ptype, x, pn, op);
}
#endif
