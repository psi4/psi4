/*! \file
    \ingroup TRANSQT
    \brief Enter brief description of file here 
*/
#ifndef _psi_bin_transqt_global_h_
#define _psi_bin_transqt_global_h_

namespace psi { namespace transqt {

extern int *ioff;
extern struct MOInfo moinfo;
extern struct Params params;

/* Globals needed for the post-backtransform sort */
extern int nbuckets;      /* number of sorting buckets   */
extern int *shell;        /* AO -> shell                 */
extern int *shell_size;   /* AOs in shell                */
extern int *bucket_map;   /* shell-pair -> sort bucket   */
extern int *bucket_offset;/* bucket -> quartet offset    */
extern int *bucket_quarts;/* no. of quartets in a bucket */
extern int *bucket_firstpq; /* First pq in bucket          */
extern int *bucket_lastpq;  /* Last pq in bucket           */

/*-------------------
  Global definitions
 -------------------*/
#define MAKE_GGGG 0            /* Make integrals of type (pq|rs) */
#define MAKE_OGOG 1            /* Make integrals of type (ip|jq) */
#define MAKE_OVOV 2            /* Make integrals of type (ia|jb) */

#define ERI 0                  /* ERIs - two-electron ints which have 8-fold 
                                  permutational symmetry:
				  (pq|rs) = (pq|sr) = (qp|rs) = (qp|sr) = 
                                  (rs|pq) = (sr|pq) = (rs|qp) = (sr|qp)     */
#define R12 0                  /* integrals of r12 - have the same 
                                  permutational symmetry as ERIs */
#define R12T1 1                /* integrals of [r12,T1] operator - have the 
                                  following symmetry:
                                  (pq|rs) = - (qp|rs) = (pq|sr) = - (qp|sr) */

#define MODE_NORMAL        0     /* normal ERI transformation   */
#define MODE_MP2R12AERI    1     /* transform ERI's for R12     */
#define MODE_MP2R12AR12    2     /* transform ints of r12       */
#define MODE_MP2R12AR12T1  3     /* transform ints of [r12,T1]  */

}} // end namespace psi::transqt
#endif // header guard
