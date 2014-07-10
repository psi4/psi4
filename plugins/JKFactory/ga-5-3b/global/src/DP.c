#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: DP.c,v 1.14 2003-02-18 00:24:32 manoj Exp $ */
#include "global.h"
#include "globalp.h"
#include "macommon.h"
#include "typesf2c.h"
#include "ga-papi.h"
#include "ga-wapi.h"

/*\ check if I own the patch
\*/
static logical own_patch(g_a, ilo, ihi, jlo, jhi)
     Integer g_a, ilo, ihi, jlo, jhi;
{
   Integer ilop, ihip, jlop, jhip, me=pnga_nodeid();
   Integer lo[2],hi[2];

   pnga_distribution(g_a, me, lo, hi);
   ilop = lo[0];
   jlop = lo[1];
   ihip = hi[0];
   jhip = hi[1];
   if(ihip != ihi || ilop != ilo || jhip != jhi || jlop != jlo) return(FALSE);
   else return(TRUE);
}

static logical patch_intersect(ilo, ihi, jlo, jhi, ilop, ihip, jlop, jhip)
     Integer ilo, ihi, jlo, jhi;
     Integer *ilop, *ihip, *jlop, *jhip;
{
     /* check consistency of patch coordinates */
     if( ihi < ilo || jhi < jlo)     return FALSE; /* inconsistent */
     if( *ihip < *ilop || *jhip < *jlop) return FALSE; /* inconsistent */

     /* find the intersection and update (ilop: ihip, jlop: jhip) */
     if( ihi < *ilop || *ihip < ilo) return FALSE; /* don't intersect */
     if( jhi < *jlop || *jhip < jlo) return FALSE; /* don't intersect */
     *ilop = GA_MAX(ilo,*ilop);
     *ihip = GA_MIN(ihi,*ihip);
     *jlop = GA_MAX(jlo,*jlop);
     *jhip = GA_MIN(jhi,*jhip);

     return TRUE;
}


/*\ COPY A PATCH 
 *
 *  . identical shapes 
 *  . copy by column order - Fortran convention
\*/
#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_copy_patch_dp = pnga_copy_patch_dp
#endif
void pnga_copy_patch_dp(t_a, g_a, ailo, aihi, ajlo, ajhi,
                   g_b, bilo, bihi, bjlo, bjhi)
     Integer g_a, ailo, aihi, ajlo, ajhi;
     Integer g_b, bilo, bihi, bjlo, bjhi;
     char *t_a;
{
Integer atype, btype, adim1, adim2, bdim1, bdim2;
Integer ilos, ihis, jlos, jhis;
Integer ilod, ihid, jlod, jhid, corr, nelem;
Integer me= pnga_nodeid(), ld, i,j;
Integer lo[2], hi[2];
Integer ldT;
char transp;
DoublePrecision *dbl_ptrA=NULL, *dbl_ptrB=NULL;
Integer ndim, dims[2];

   pnga_check_handle(g_a, "pnga_copy_patch_dp");
   pnga_check_handle(g_b, "pnga_copy_patch_dp");

   /* if(g_a == g_b) pnga_error("pnga_copy_patch_dp: arrays have to different ", 0L); */

   pnga_inquire(g_a, &atype, &ndim, dims);
   adim1 = dims[0];
   adim2 = dims[1];
   pnga_inquire(g_b, &btype, &ndim, dims);
   bdim1 = dims[0];
   bdim2 = dims[1];

   if(atype != btype || (atype != C_DBL ))
      pnga_error("pnga_copy_patch_dp: wrong types ", 0L);

   /* check if patch indices and dims match */
   if (ailo <= 0 || aihi > adim1 || ajlo <= 0 || ajhi > adim2)
       pnga_error(" pnga_copy_patch_dp: g_a indices out of range ", 0L);
   if (bilo <= 0 || bihi > bdim1 || bjlo <= 0 || bjhi > bdim2)
       pnga_error(" pnga_copy_patch_dp: g_b indices out of range ", 0L);

   /* check if numbers of elements in two patches match each other */
   if (((bihi - bilo + 1)  != (aihi - ailo + 1)) || 
      ( (bjhi - bjlo + 1)  != (ajhi - ajlo + 1)) )
       pnga_error(" pnga_copy_patch_dp: shapes two of patches do not match ", 0L);

    /* is transpose operation required ? */
   transp = (*t_a == 'n' || *t_a =='N')? 'n' : 't';

   /* now find out cordinates of a patch of g_a that I own */
   pnga_distribution(g_a, me, lo, hi);
   ilos = lo[0];
   jlos = lo[1];
   ihis = hi[0];
   jhis = hi[1];

   if(patch_intersect(ailo, aihi, ajlo, ajhi, &ilos, &ihis, &jlos, &jhis)){
      pnga_access_ptr(g_a, lo, hi, &dbl_ptrA, &ld);
      
      nelem = (ihis-ilos+1)*(jhis-jlos+1);
      
      if ( transp == 'n' ) {
	  corr  = bilo - ailo;
	  ilod  = ilos + corr; 
	  ihid  = ihis + corr;
	  corr  = bjlo - ajlo;
	  jlod  = jlos + corr; 
	  jhid  = jhis + corr;
      } else {
	/* If this is a transpose copy, we need local scratch space */
	dbl_ptrB = (DoublePrecision*) ga_malloc(nelem,MT_F_DBL,"copypatch_dp");

	  /* Copy from the source into this local array, transposed */
	  ldT = jhis-jlos+1;
	  
	  for(j=0; j< jhis-jlos+1; j++)
	      for(i=0; i< ihis-ilos+1; i++)
		  *(dbl_ptrB + i*ldT + j) = *(dbl_ptrA + j*ld + i);

	  /* Now we can reset index to point to the transposed stuff */
      pnga_release(g_a, lo, hi);
	  dbl_ptrA = dbl_ptrB;
	  ld = ldT;

	  /* And finally, figure out what the destination indices are */
	  corr  = bilo - ajlo;
	  ilod  = jlos + corr; 
	  ihid  = jhis + corr;
	  corr  = bjlo - ailo;
	  jlod  = ilos + corr; 
	  jhid  = ihis + corr;
      }
	  
      /* Put it where it belongs */
      lo[0] = ilod;
      lo[1] = jlod;
      hi[0] = ihid;
      hi[1] = jhid;
      pnga_put(g_b, lo, hi, dbl_ptrA, &ld);

      /* Get rid of local memory if we used it */
      if( transp == 't') ga_free(dbl_ptrB);
  }
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_ddot_patch_dp = pnga_ddot_patch_dp
#endif
DoublePrecision pnga_ddot_patch_dp(g_a, t_a, ailo, aihi, ajlo, ajhi,
                                  g_b, t_b, bilo, bihi, bjlo, bjhi)
     Integer g_a, ailo, aihi, ajlo, ajhi;    /* patch of g_a */
     Integer g_b, bilo, bihi, bjlo, bjhi;    /* patch of g_b */
     char    *t_a, *t_b;                          /* transpose operators */
{
Integer atype, btype, adim1, adim2, bdim1, bdim2;
Integer iloA, ihiA, jloA, jhiA, ldA;
Integer iloB, ihiB, jloB, jhiB, ldB;
Integer alo[2], ahi[2];
Integer blo[2], bhi[2];
Integer g_A = g_a;
Integer me= pnga_nodeid(), i, j, temp_created=0;
Integer corr, nelem;
char    transp, transp_a, transp_b;
DoublePrecision  sum = 0.;
DoublePrecision *dbl_ptrA;
DoublePrecision *dbl_ptrB;
Integer ndim, dims[2];

   pnga_check_handle(g_a, "pnga_ddot_patch_dp");
   pnga_check_handle(g_b, "pnga_ddot_patch_dp");

   pnga_inquire(g_a, &atype, &ndim, dims);
   adim1 = dims[0];
   adim2 = dims[1];
   pnga_inquire(g_b, &btype, &ndim, dims);
   bdim1 = dims[0];
   bdim2 = dims[1];

   if(atype != btype || (atype != C_DBL ))
      pnga_error("pnga_ddot_patch_dp: wrong types ", 0L);

  /* check if patch indices and g_a dims match */
   if (ailo <= 0 || aihi > adim1 || ajlo <= 0 || ajhi > adim2)
      pnga_error(" pnga_ddot_patch_dp: g_a indices out of range ", 0L);

   /* check if patch indices and g_b dims match */
   if (bilo <= 0 || bihi > bdim1 || bjlo <= 0 || bjhi > bdim2)
       pnga_error(" pnga_ddot_patch_dp: g_b indices out of range ", 0L);


   /* is transpose operation required ? */
   /* -- only if for one array transpose operation requested*/
   transp_a = (*t_a == 'n' || *t_a =='N')? 'n' : 't';
   transp_b = (*t_b == 'n' || *t_b =='N')? 'n' : 't';
   transp   = (transp_a == transp_b)? 'n' : 't';
   if(transp == 't')
          pnga_error(" pnga_ddot_patch_dp: transpose operators don't match: ", me);


   /* find out coordinates of patches of g_A and g_B that I own */
   pnga_distribution(g_A, me, alo, ahi);
   iloA = alo[0];
   jloA = alo[1];
   ihiA = ahi[0];
   jhiA = ahi[1];

   if (patch_intersect(ailo, aihi, ajlo, ajhi, &iloA, &ihiA, &jloA, &jhiA)){

       pnga_access_ptr(g_A, alo, ahi, &dbl_ptrA, &ldA);
       nelem = (ihiA-iloA+1)*(jhiA-jloA+1);

       corr  = bilo - ailo;
       iloB  = iloA + corr;
       ihiB  = ihiA + corr;
       corr  = bjlo - ajlo;
       jloB  = jloA + corr;
       jhiB  = jhiA + corr;
       blo[0] = iloB; blo[1] = jloB;
       bhi[0] = ihiB; bhi[1] = jhiB;

      if(own_patch(g_b, iloB, ihiB, jloB, jhiB)){
         /* all the data is local */
         pnga_access_ptr(g_b, blo, bhi, &dbl_ptrB, &ldB);
      }else{
         /* data is remote -- get it to temp storage*/
         temp_created =1;
	 dbl_ptrB = (DoublePrecision*)ga_malloc(nelem, MT_F_DBL, "ddot_dp_b");

         ldB   = ihiB-iloB+1; 
         pnga_get(g_b, blo, bhi, dbl_ptrB, &ldB);
      }

      sum = 0.;
      for(j=0; j< jhiA-jloA+1; j++)
          for(i=0; i< ihiA-iloA+1; i++)
             sum += *(dbl_ptrA + j*ldA + i) * 
                    *(dbl_ptrB + j*ldB + i);
      pnga_release(g_A, alo, ahi);

      if(temp_created) ga_free(dbl_ptrB);
      else pnga_release(g_b, blo, bhi);
   }
   return sum;
}

