#include "utils.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

template <typename val_t, typename idx_t>
static void
qs(val_t* item, idx_t* index, idx_t left, idx_t right)
{
  register idx_t i,j;
  val_t x;
  idx_t y;
 
  i=left; j=right;
  x=item[index[(left+right)/2]];
 
  do {
    while(item[index[i]]<x && i<right) i++;
    while(x<item[index[j]] && j>left) j--;
 
    if (i<=j) {
      if (item[index[i]] != item[index[j]]) {
        y=index[i];
        index[i]=index[j];
        index[j]=y;
        }
      i++; j--;
      }
    } while(i<=j);
       
  if (left<j) qs(item,index,left,j);
  if (i<right) qs(item,index,i,right);
}

template <typename val_t, typename idx_t>
static void
quicksort(val_t* item, idx_t* index, idx_t n)
{
  idx_t i;
  if (n<=0) return;
  for (i=0; i<n; i++) {
    index[i] = i;
    }
  qs<val_t, idx_t>(item,index,0,n-1);
}

static void
quicksort(int* item, int* index, int n)
{
    quicksort<int,int>(item, index, n);
}

static void
quicksort(usi* item, usi* index, usi n)
{
    quicksort<usi,usi>(item, index, n);
}
  
}

