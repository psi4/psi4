#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file iterator.h
 *  @author Sriram Krishnamoorthy
 *  @brief Stride iterator.
 *  An iterator for the stride descriptor to reuse common traversal
 *  functionality. More functionality related to the strided
 *  descriptor reusable across files will be extracted here as well. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "iterator.h"
#include "armci.h"


typedef struct {
  void *base_ptr;
  int stride_levels;
  int stride_arr[ARMCI_MAX_STRIDE_LEVEL];
  int seg_count[ARMCI_MAX_STRIDE_LEVEL+1];  

  int size, pos, itr[ARMCI_MAX_STRIDE_LEVEL];
} stride_info_t;


/**Create a stride iterator.
 * @param base_ptr IN Starting pointer for stride descriptor
 * @param stride_levels IN #stride levels
 * @param stride_arr IN the strides (arr of size[stride_levels])
 * @param seg_count IN #segments in each stride
 * level([stride_levels+1])  
 * @return Handle to stride iterator created
 */
stride_itr_t armci_stride_itr_init(void *base_ptr,
				   int stride_levels,
				   const int *stride_arr,
				   const int *seg_count) {
  int i;
  stride_info_t *sinfo = malloc(sizeof(stride_info_t));
  assert(sinfo!=NULL);
  assert(stride_levels>=0);
  assert(stride_levels<=ARMCI_MAX_STRIDE_LEVEL);
  for(i=0; i<stride_levels; i++) {
    if(i==0) 
      assert(stride_arr[0] >= seg_count[0]);
    else 
      assert(stride_arr[i] >= stride_arr[i-1]*seg_count[i]);
  }

  sinfo->base_ptr= base_ptr;
  sinfo->stride_levels = stride_levels;
  for(i=0; i<stride_levels; i++) {
    sinfo->stride_arr[i] = stride_arr[i];
  }
  for(i=0; i<stride_levels+1; i++) {
    sinfo->seg_count[i] = seg_count[i];
  }
  sinfo->size=1;
  for(i=1; i<stride_levels+1; i++) {
    sinfo->size *= sinfo->seg_count[i];
  }
  assert(sinfo->size>0);
  sinfo->pos=0;
  for(i=0; i<stride_levels+1; i++) {
    sinfo->itr[i] = 0;
  }
  return sinfo;
}


/**Destroy a stride iterator.
 * @param psitr IN/OUT Pointer to stride iterator
 * @return void
 */
void armci_stride_itr_destroy(stride_itr_t *psitr) {
  free(*psitr);
  *psitr = NULL; /*for a safe crash*/
}


/**Size of the stride iterator. Defined as total #contiguous
 * segments in the stride iterator.
 * @param sitr IN Handle to stride iterator
 * @return Size of the stride iterator
 */
int armci_stride_itr_size(stride_itr_t sitr) {
  stride_info_t *sinfo = (stride_info_t *)sitr;
  assert(sinfo!=NULL);
  return sinfo->size;
}


/**Position of the stride iterator. Between 0 and (size-1),
 * inclusive. Position is the index of the contiguous segment
 * currently traversed by the iterator.
 * @param sitr IN Handle to stride descriptor
 * @return Position of the iterator
 */
int armci_stride_itr_pos(stride_itr_t sitr) {
  stride_info_t *sinfo = (stride_info_t *)sitr;
  assert(sinfo!=NULL);
  return sinfo->pos;
}


/**Move the iterator to the next position. Assumes position<=size.
 * @param sitr IN Handle to stride descriptor
 * @return void
 */
void armci_stride_itr_next(stride_itr_t sitr) {
  int i;
  stride_info_t *sinfo = (stride_info_t *)sitr;
  assert(sinfo!=NULL);
  assert(sinfo->pos <sinfo->size);
  sinfo->pos += 1;
  if(sinfo->stride_levels>0) {
    sinfo->itr[0] += 1;
    for(i=0; i<sinfo->stride_levels-1 && sinfo->itr[i]==sinfo->seg_count[i+1]; i++) {
      sinfo->itr[i] = 0;
      sinfo->itr[i+1] += 1;
    }
    assert(sinfo->itr[i] <= sinfo->seg_count[i+1]);
  }
}

/**Get pointer to the contiguous segment currently being
 * traversed. This is the pointer to the user buffer.
 * @param sitr IN Handle to stride descriptor
 * @return pointer to current contiguous segment
 */
void *armci_stride_itr_seg_ptr(stride_itr_t sitr) {
  stride_info_t *sinfo = (stride_info_t *)sitr;
  assert(sinfo!=NULL);
  return sinfo->base_ptr + armci_stride_itr_seg_off(sitr);
}

/**Get the size of the current segment.
 * @param sitr IN Handle to stride descriptor
 * @return Size of the current segment
 */
int armci_stride_itr_seg_size(stride_itr_t sitr) {
  stride_info_t *sinfo = (stride_info_t *)sitr;
  assert(sinfo!=NULL);
  return sinfo->seg_count[0];
}


/**Get the offset of the current segment with respect to the start of
 * the first segment (a.k.a src_ptr)
 * @param sitr IN Handle to stride descriptor
 * @return Offset of the current segment
 */
int armci_stride_itr_seg_off(stride_itr_t sitr) {
  int i;
  int off;
  stride_info_t *sinfo = (stride_info_t *)sitr;
  assert(sinfo!=NULL);
  
  off=0;
  for(i=0; i<sinfo->stride_levels; i++) {
    off += sinfo->itr[i] * sinfo->stride_arr[i];
  }
  return off;
}


/**Check if there are more segments to iterate over. (a.k.a
 * position<size).
 * @param sitr IN Handle to stride descriptor
 * @return Zero if current position is past the size of the
 * iterator. Non-zero otherwise. 
 */
int armci_stride_itr_has_more(stride_itr_t sitr) {
  stride_info_t *sinfo = (stride_info_t *)sitr;
  assert(sinfo!=NULL);
  return sinfo->pos<sinfo->size;
}


void armci_write_strided(
        void *ptr, int stride_levels, int stride_arr[], int count[], char *buf)
{
  const int seg_size = count[0];
  int off=0;
  assert(count[0]>0);
  stride_itr_t sitr=armci_stride_itr_init(ptr,stride_levels,stride_arr,count);
  while(armci_stride_itr_has_more(sitr)) {
    char *sptr = armci_stride_itr_seg_ptr(sitr);
    memcpy(&buf[off],sptr,seg_size);
    off += seg_size;
    armci_stride_itr_next(sitr);
  }
  armci_stride_itr_destroy(&sitr);
}


void armci_read_strided(
        void *ptr, int stride_levels, int stride_arr[], int count[], char *buf)
{
  const int seg_size = count[0];
  int off=0;
  assert(count[0]>0);
  stride_itr_t sitr=armci_stride_itr_init(ptr,stride_levels,stride_arr,count);
  while(armci_stride_itr_has_more(sitr)) {
    char *dptr = armci_stride_itr_seg_ptr(sitr);
    memcpy(dptr,&buf[off],seg_size);
    off += seg_size;
    armci_stride_itr_next(sitr);
  }
  armci_stride_itr_destroy(&sitr);
}

