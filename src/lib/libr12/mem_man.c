/*! \file
    \ingroup R12
    \brief Enter brief description of file here 
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"mem_man.h"
#define MAXALLOC 10000

int free_block[MAXALLOC];   /* Marks of free-occupied blocks */
int block_length[MAXALLOC]; /* Keeps sizes of all blocks */
int last_free;              /* Number of the last free block */
int max_mem;                /* Amount of memory used by classes up to this point */
int mem_top;                /* Total memory allocated up to this point */

/* initialize memory stack */
void init_mem(int memory)
{
  last_free = 1;
  memset(free_block,0,MAXALLOC*sizeof(int));    /* All blocks are free */
  memset(block_length,0,MAXALLOC*sizeof(int));  /* All of them are zero in length ... */
  block_length[0] = memory;                     /* except for the first one */
  max_mem = 0;
  mem_top = memory;
}

/* add memory */
void add_mem(int memory)
{
  int i;
  int addto = 0;
  int himem = 0;

  mem_top += memory;
  /* Find the last free block and add all this memory to it */
  for(i=0; i<last_free; i++){
    if(free_block[i] > himem){
      himem = free_block[i];
      addto = i;
      }
    }
#if 0
  printf("adding %d to memory\n\tnew top = %d\n", memory,
    free_block[addto]+block_length[addto]+memory);
#endif
  block_length[addto] = block_length[addto]+memory;
}

/* Find a block of the requested size */
int get_mem(int size)
{
  int i, j;

/* try to find one that fits exactly */
  for(i=last_free-1; i>=0; i--){
    if(block_length[i] == size){
      j = free_block[i];
      if(free_block[i]+block_length[i]==mem_top) add_mem(500);
      use(i, size);
      return j;
      }
    }

/* ok, try to find a bigger one that will work */
  for(i=last_free-1; i>=0; i--){
    if(block_length[i] > size){
      j = free_block[i];
      use(i, size);
      return j;
      }
    }

  
/* last resort, expand memory */
  add_mem(1000);
  return get_mem(size);

}


void use(int n, int s)
{
  int i;
#if 0
  printf("issuing %d doubles starting at %d out of free block %d\n",
        s, free_block[n], n);
  printf("last_free = %d\n", last_free);
#endif
  if(s==block_length[n]){
    for(i=n; i<last_free; i++){
      free_block[i] = free_block[i+1];
      block_length[i] = block_length[i+1];
      }
    last_free--;
    }
  else if(s<block_length[n]){
    free_block[n] = free_block[n]+s;
    block_length[n] = block_length[n]-s;
/* keep track of how much memory we use overall */
    if(free_block[n]>max_mem) max_mem = free_block[n];
    }
  else exit(1);

}

void free_mem(int n, int size)
{
  int i, j;

  free_block[last_free] = n;
  block_length[last_free] = size;
  last_free++;

  consolidate();

}

void consolidate()
{
  int i, j;
  int right_bound_i;
  int done = 1;

  do {
    done = 1;
    for(i=0; i<last_free; i++){
      right_bound_i = free_block[i]+block_length[i];
      for(j=0; j<last_free; j++){
        if(free_block[j]==right_bound_i){
          block_length[i]+=block_length[j];
          use(j, block_length[j]);
          done = 0;
          }
        }
      }
    } while (!done);
#if 0
  printf("after consolidation, free_blocks looks like:\n");
  for(i=1; i<=last_free; i++){
    printf("%d[%d]\t", free_block[i-1], block_length[i-1]);
    if(!(i%5)) printf("\n");
    }
  printf("\n");
#endif

}

int get_total_memory()
{
  return max_mem;
}

