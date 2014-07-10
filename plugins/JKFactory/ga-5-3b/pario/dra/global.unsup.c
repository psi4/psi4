#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * This file contains some data parallel GA operations that are not part of 
 * the official GA distribution.
 * 
 * Date:   05.10.95
 * Author: Jarek Nieplocha
 */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "global.h"
#include "macommon.h"

#define ga_idot_ F77_FUNC_(ga_idot,GA_IDOT)

#define DRA_TYPE_GSM 32760 - 6

/**
 * Integer version of ga_ddot 
 */
Integer ga_idot_(Integer *g_a, Integer *g_b)
{
    Integer  atype, adim1, adim2, btype, bdim1, bdim2, ald, bld;
    Integer  ailo,aihi, ajlo, ajhi, bilo, bihi, bjlo, bjhi;
    register Integer i,j;
    Integer  me, sum;
    Integer  index_a, index_b;

    pnga_sync();

    me = ga_nodeid_();

    ga_check_handle(g_a, "ga_idot");
    ga_check_handle(g_b, "ga_idot");

    ga_inquire_(g_a,  &atype, &adim1, &adim2);
    ga_inquire_(g_b,  &btype, &bdim1, &bdim2);

    if(atype != btype || atype != MT_F_INT)
        ga_error("ga_idot: types not correct", 0L);

    if (adim1!=bdim1 || adim2 != bdim2)
        ga_error("ga_idot: arrays not conformant", 0L);

    if (DBL_MB == (DoublePrecision*)0 || INT_MB == (Integer*)0)
        ga_error("ga_idot: null pointer for base array",0L);

    ga_distribution_(g_a, &me, &ailo, &aihi, &ajlo, &ajhi);
    ga_distribution_(g_b, &me, &bilo, &bihi, &bjlo, &bjhi);

    if (ailo!=bilo || aihi != bihi || ajlo!=bjlo || ajhi != bjhi){
        /*
           fprintf(stderr,"\nme =%d: %d-%d %d-%d vs %d-%d %d-%d dim:%dx%d\n",me,
           ailo,aihi, ajlo, ajhi, bilo, bihi, bjlo, bjhi,adim1,adim2);
           */
        ga_error("ga_idot: distributions not identical",0L);
    }

    sum = 0.;
    if (  aihi>0 && ajhi>0 ){
        ga_access_(g_a, &ailo, &aihi, &ajlo, &ajhi,  &index_a, &ald);
        if(g_a == g_b){
            index_b = index_a; bld =ald;
        }else
            ga_access_(g_b, &bilo, &bihi, &bjlo, &bjhi,  &index_b, &bld);

        index_a --;  /* Fortran to C correction of starting address */
        index_b --;  /* Fortran to C correction of starting address */


        /* compute "local" contribution to the dot product */
        for(j=0; j<ajhi-ajlo+1; j++)
            for(i=0; i<aihi-ailo+1; i++)
                sum += INT_MB[index_a +j*ald + i]  *
                    INT_MB[index_b +j*bld + i];

        /* release access to the data */
        ga_release_(g_a, &ailo, &aihi, &ajlo, &ajhi);
        ga_release_(g_b, &bilo, &bihi, &bjlo, &bjhi);
    }

    ga_igop((Integer)DRA_TYPE_GSM, &sum, (Integer)1, "+");

    pnga_sync();

    return (sum);
}
