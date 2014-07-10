#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * operations on sections of 2D arrays
 */

#include "global.h"       /* used only to define datatypes */
#include "drap.h"

#define PARIO_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define PARIO_MIN(a,b) (((a) <= (b)) ? (a) : (b))


/**
 * check if two patches are conforming (dimensions are divisible)
 */
logical dai_patches_conforming(
        Integer* ailo, Integer* aihi, Integer* ajlo, Integer* ajhi, 
        Integer* bilo, Integer* bihi, Integer* bjlo, Integer* bjhi)
{
    Integer mismatch;
    Integer adim1, bdim1, adim2, bdim2;
    adim1 = *aihi - *ailo +1;
    adim2 = *ajhi - *ajlo +1;
    bdim1 = *bihi - *bilo +1;
    bdim2 = *bjhi - *bjlo +1;
    mismatch  = (adim1<bdim1) ? bdim1%adim1 : adim1%bdim1;
    mismatch += (adim2<bdim2) ? bdim2%adim2 : adim2%bdim2;
    if(mismatch)return(FALSE);
    else return(TRUE);
}


/**
 * check if patches are identical
 */
logical dai_comp_patch(
        Integer* ilo, Integer* ihi, Integer* jlo, Integer* jhi, 
        Integer* ilop, Integer* ihip, Integer* jlop, Integer* jhip)
{
    if(*ihip != *ihi || *ilop != *ilo || *jhip != *jhi || *jlop != *jlo)
        return(FALSE);
    else return(TRUE);
}


/**
 * check if patches have a nontrivial intersection
 * if yes, find the intersection and update [ilop:ihip, jlop:jhip]
 */
logical dai_patch_intersect(
        Integer ilo, Integer ihi, Integer jlo, Integer jhi,
        Integer* ilop, Integer* ihip, Integer* jlop, Integer* jhip)
{
    /* check consistency of patch coordinates */
    if( ihi < ilo || jhi < jlo)     return FALSE; /* inconsistent */
    if( *ihip < *ilop || *jhip < *jlop) return FALSE; /* inconsistent */

    /* find the intersection and update (ilop: ihip, jlop: jhip) */
    if( ihi < *ilop || *ihip < ilo) return FALSE; /* don't intersect */
    if( jhi < *jlop || *jhip < jlo) return FALSE; /* don't intersect */
    *ilop = PARIO_MAX(ilo,*ilop);
    *ihip = PARIO_MIN(ihi,*ihip);
    *jlop = PARIO_MAX(jlo,*jlop);
    *jhip = PARIO_MIN(jhi,*jhip);

    return(TRUE);
}


/**
 * check if sections have a nontrivial intersection
 * if yes, find the intersection and update [ilop:ihip, jlop:jhip]
 * section format
 */
logical dai_section_intersect(section_t sref, section_t *sadj)
{
    Integer ndim = sref.ndim;
    Integer i;
    logical isconsistent = TRUE;
    /* check that patches have same dimension */
    if (ndim != sadj->ndim) isconsistent = FALSE;
    /* check consistency of patch coordinates */
    if (isconsistent) {
        for (i=0; i<ndim; i++) {
            if (sref.hi[i] < sref.lo[i]) isconsistent = FALSE;
            if (sadj->hi[i] < sadj->lo[i]) isconsistent = FALSE;
        }
    }
    /* check to see if there is an intersection */
    if (isconsistent) {
        for (i=0; i<ndim; i++) {
            if (sref.hi[i] < sadj->lo[i]) isconsistent = FALSE;
            if (sadj->hi[i] < sref.lo[i]) isconsistent = FALSE;
        }
    }
    /* if there is an intersection then return it in sadj */
    if (isconsistent) {
        for (i=0; i<ndim; i++) {
            sadj->lo[i] = PARIO_MAX(sref.lo[i],sadj->lo[i]);
            sadj->hi[i] = PARIO_MIN(sref.hi[i],sadj->hi[i]);
        }
    }
    return (isconsistent);
}
