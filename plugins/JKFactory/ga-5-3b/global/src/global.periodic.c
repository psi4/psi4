#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*              
 * module: global.periodic.c
 * author: Author: Jialin Ju, PNNL
 * description: periodic put/get/acc
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with United States
 * Government support under Contract Number DE-AC06-76RLO-1830 awarded by
 * the United States Department of Energy.  The United States Government
 * retains a paid-up non-exclusive, irrevocable worldwide license to
 * reproduce, prepare derivative works, perform publicly and display
 * publicly by or for the US Government, including the right to
 * distribute to other US Government contractors.
 */


#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#include "globalp.h"
#include "ga-papi.h"
#include "ga-wapi.h"

#define RANGE_NUM 3
#define RANGE_BOUND 6
#define IS_REGULAR_PATCH 100

typedef struct {
    Integer lo;
    Integer hi;
    int isvalid;
} Range_bound;

typedef struct {
    Range_bound low;
    Range_bound mid;
    Range_bound hig;
} Range;


static int ngai_peri_get_range_(Integer ndim, Integer *dims, Integer *lo_orig, Integer *hi_orig,
                         Integer range[][RANGE_BOUND], Integer range_num[],
                         Integer offset[][RANGE_BOUND/2], Integer op_code)
/*
  ndim                   - dimension of array [input]
  dims[]                 - dimensions of each array component [input]
  lo_orig[]              - lower bounds of request
  hi_orig[]              - upper bounds of request
  range[][]              - adjusted ranges for requests [output]
  range_num[]            - number of ranges per axis 
  offset[]               - offset between original request and shifted request [output]
  op_code                - not used
*/
{
    Integer lo[MAXDIM], hi[MAXDIM];
    Integer i, rmndr, id;
    Range range_raw[MAXDIM];
    Integer is_regular_patch, no_shift;
 
    /* check the patch is valid or not */
    for (i=0; i<ndim; i++) {
      lo[i] = lo_orig[i];
      hi[i] = hi_orig[i];
      if ((hi[i] - lo[i] + 1) < 1 || (hi[i]-lo[i]+1) > 2*dims[i]) return 0;
    }
    /* move request so that it at least partially overlaps the array
     * (if necessary)
     */
    no_shift = 1;
    for(i=0; i<ndim; i++) {
      if (lo[i] > dims[i]) {
        rmndr = lo[i]%dims[i];
        id = (lo[i]-rmndr)/dims[i]; 
        lo[i] = lo[i] - id*dims[i];
        hi[i] = hi[i] - id*dims[i];
        if (id != 0) {
          no_shift = 0;
        }
      }
      if (hi[i] < 1) {
        rmndr = -hi[i]%dims[i];
        id = (hi[i]+rmndr)/dims[i]+1; 
        lo[i] = lo[i] + id*dims[i];
        hi[i] = hi[i] + id*dims[i];
        if (id != 0) {
          no_shift = 0;
        }
      }
    }
    /* break the lo and hi into the corresponding ranges
     *
     * lo (if < 1)         1                dims[i]      hi (if > dims[i])
     * |-------------------|------------------|------------------|
     * Example:
     *    lo[0] = -3 hi[0] = 5 dims[0] = 10 then
     *    range_raw[0].low.lo = -3  range_raw[0].low.hi = 0
     *    range_raw[0].low.isvalid = 1
     *    range_raw[0].mid.lo = 1   range_raw[0].mid.hi = 5
     *    range_raw[0].mid.isvalid = 1
     *    range_raw[0].hig.isvalid = 0
     */
    for(i=0; i<ndim; i++) {
        /* initialize the dim range */
        range_raw[i].low.lo = 0; range_raw[i].low.hi = 0;
        range_raw[i].low.isvalid = 0;
        range_raw[i].mid.lo = 1; range_raw[i].mid.hi = dims[i];
        range_raw[i].mid.isvalid = 0;
        range_raw[i].hig.lo = dims[i]+1; range_raw[i].hig.hi = 0;
        range_raw[i].hig.isvalid = 0;
        
        if(lo[i] < 1) {
            range_raw[i].low.lo = lo[i];
            range_raw[i].low.isvalid = 1;
        } else if (lo[i] >= 1) {
            range_raw[i].mid.lo = lo[i];
            range_raw[i].mid.isvalid = 1;
        }

        if (lo[i] < 1 && hi[i] > dims[i]) {
            range_raw[i].mid.isvalid = 1;
        }


        if(hi[i] > dims[i]) {
            range_raw[i].hig.hi = hi[i];
            range_raw[i].hig.isvalid = 1;
        } else if (hi[i] <= dims[i]) {
            range_raw[i].mid.hi = hi[i];
            range_raw[i].mid.isvalid = 1;
        }

    }

    /* check if this is a regular patch, not periodic operation needed */
    is_regular_patch = 1;
    for(i=0; i<ndim; i++) 
        if(range_raw[i].low.isvalid || range_raw[i].hig.isvalid)
            is_regular_patch = 0;

    if(is_regular_patch && no_shift) return IS_REGULAR_PATCH;
    
    /* adjust the range so that they are in the range of 1 to dims[i] */
    for(i=0; i<ndim; i++) {
        range_num[i] = 0;
        if(range_raw[i].low.isvalid) {
            range[i][range_num[i]*2] = range_raw[i].low.lo + dims[i];
            range[i][range_num[i]*2+1] = range_raw[i].low.hi + dims[i];
            offset[i][range_num[i]] = range_raw[i].low.lo - lo[i];
            range_num[i]++;
        }
        if(range_raw[i].mid.isvalid) {
            range[i][range_num[i]*2] = range_raw[i].mid.lo;
            range[i][range_num[i]*2+1] = range_raw[i].mid.hi;
            offset[i][range_num[i]] = range_raw[i].mid.lo - lo[i];
            range_num[i]++;
        }
        if(range_raw[i].hig.isvalid) {
            range[i][range_num[i]*2] = range_raw[i].hig.lo - dims[i];
            range[i][range_num[i]*2+1] = range_raw[i].hig.hi - dims[i];
            offset[i][range_num[i]] = range_raw[i].hig.lo - lo[i];
            range_num[i]++;
        }  
    }
    
    return 1;
}

#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak wnga_periodic = pnga_periodic
#endif
void pnga_periodic(Integer g_a, Integer *lo, Integer *hi, void *buf,
                    Integer *ld, void *alpha, Integer op_code)
{
    int i, j, counter[MAXDIM], done;
    Integer type, ndim, dims[MAXDIM];
    Integer range[MAXDIM][RANGE_BOUND], range_num[MAXDIM];
    Integer offset[MAXDIM][RANGE_BOUND/2], my_offset, temp_offset;
    Integer lop[MAXDIM], hip[MAXDIM];
    int get_range;
    
    pnga_inquire(g_a, &type, &ndim, dims);

    get_range = ngai_peri_get_range_(ndim, dims, lo, hi, range, range_num,
                                     offset, op_code);
    
    if(!get_range) pnga_error("g_a indices are invalid ", 0L);

    /* If this is a regular patch, not periodic operation needed */
    if(get_range == IS_REGULAR_PATCH) {
        switch(op_code) {
          case PERIODIC_GET:
              pnga_get(g_a, lo, hi, buf, ld);
              break;
          case PERIODIC_PUT:
              pnga_put(g_a, lo, hi, buf, ld);
              break;
          case PERIODIC_ACC:    
              pnga_acc(g_a, lo, hi, buf, ld, alpha);
              break;
          default:
              pnga_error("This operation is invalid ", 0L);
        }
        return;
    }

    /* periodic operation continues */
    for(i=0; i<ndim; i++) counter[i] = 0;

    done = 0;
    do {
        my_offset = 0;
        
        for(i=0; i<ndim; i++) {
            lop[i] = range[i][counter[i]*2];
            hip[i] = range[i][counter[i]*2+1];

            /* calculate the offset */
            if(i == 0) my_offset += offset[i][counter[i]];
            else {
                temp_offset = offset[i][counter[i]];
                for(j=0; j<i; j++)
                    /* temp_offset *= dims[j]; */
                    temp_offset *= ld[j];
                my_offset += temp_offset;
            }
        }
        
        /* deal with this patch */
        switch(op_code) {
          case PERIODIC_GET:
              pnga_get(g_a, lop, hip,
                       (char *)buf+my_offset*GAsizeofM(type), ld);
              break;
          case PERIODIC_PUT:
              pnga_put(g_a, lop, hip,
                       (char *)buf+my_offset*GAsizeofM(type), ld);
              break;
          case PERIODIC_ACC:
              pnga_acc(g_a, lop, hip,
                       (char *)buf+my_offset*GAsizeofM(type), ld, alpha);
              break;
          default:
              pnga_error("This operation is invalid ", 0L);
        }      

        counter[0]++;
        for(i=0; i<ndim; i++) {
            if(counter[i] == range_num[i]) {
                counter[i] = 0;
                if(i == (ndim-1)) { done = 1; break; }
                else counter[i+1]++;
            }
        }
    } while(!done);
}

