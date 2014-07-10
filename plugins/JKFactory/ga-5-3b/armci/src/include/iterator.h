/** @file iterator.h
 *  @author Sriram Krishnamoorthy
 *  @brief Stride iterator.
 *  An iterator for the stride descriptor to reuse common traversal
 *  functionality. More functionality related to the strided
 *  descriptor reusable across files will be extracted here as well. 
 */
#ifndef _STRIDE_INFO_H_
#define _STRIDE_INFO_H_

#include "armci.h" /*for ARMCI_MAX_STRIDE_LEVEL and dassert*/

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

typedef struct {
  char *base_ptr;
  int stride_levels;
  int stride_arr[ARMCI_MAX_STRIDE_LEVEL];
  int seg_count[ARMCI_MAX_STRIDE_LEVEL+1];  

  int size, pos, itr[ARMCI_MAX_STRIDE_LEVEL];
} stride_info_t;

void armci_stride_info_init(stride_info_t *sinfo,
		      void *base_ptr,
		      int stride_levels,
		      const int *stride_arr,
		      const int *seg_count);
  
void  armci_stride_info_destroy(stride_info_t *sinfo);
int   armci_stride_info_count(stride_info_t *sinfo);
int   armci_stride_info_pos(stride_info_t *sinfo);
void  armci_stride_info_next(stride_info_t *sinfo);
void *armci_stride_info_seg_ptr(stride_info_t *sinfo);
int   armci_stride_info_seg_size(stride_info_t *sinfo);
int   armci_stride_info_seg_off(stride_info_t *sinfo);
int   armci_stride_info_size(stride_info_t *sinfo);
int   armci_stride_info_has_more(stride_info_t *sinfo);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif /*_STRIDE_INFO_H_*/
