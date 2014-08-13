/** @file iterator.h
 *  @author Sriram Krishnamoorthy
 *  @brief Stride iterator.
 *  An iterator for the stride descriptor to reuse common traversal
 *  functionality. More functionality related to the strided
 *  descriptor reusable across files will be extracted here as well. 
 */
#ifndef _STRIDE_ITR_H_
#define _STRIDE_ITR_H_

#if defined(__cplusplus) || defined(c_plusplus)
extern "C" {
#endif

typedef void *stride_itr_t;

stride_itr_t armci_stride_itr_init(void *base_ptr,
				   int stride_levels,
				   const int *stride_arr,
				   const int *seg_count);

void  armci_stride_itr_destroy(stride_itr_t *psitr);
int   armci_stride_itr_count(stride_itr_t sitr);
int   armci_stride_itr_pos(stride_itr_t sitr);
void  armci_stride_itr_next(stride_itr_t sitr);
void *armci_stride_itr_seg_ptr(stride_itr_t sitr);
int   armci_stride_itr_seg_size(stride_itr_t sitr);
int   armci_stride_itr_seg_off(stride_itr_t sitr);
int   armci_stride_itr_size(stride_itr_t sitr);
int   armci_stride_itr_has_more(stride_itr_t sitr);

#if defined(__cplusplus) || defined(c_plusplus)
}
#endif

#endif /*_STRIDE_ITR_H_*/
