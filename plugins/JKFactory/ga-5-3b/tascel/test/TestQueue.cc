#if HAVE_CONFIG_H
# include "config.h"
#endif

#include <cstdio>

#include <ga.h>
#include <mpi.h>

#if SPLIT_QUEUE
# include "UniformTaskCollectionSplit.h"
#elif SHARED_QUEUE
# include "UniformTaskCollectionShared.h"
#else
# error one of SPLIT_QUEUE or SHARED_QUEUE must be defined
#endif

using namespace std;
using namespace tascel;


#define MAX_TASKS 100
static int me;


struct TskDscr {
  int v;
  TskDscr(int _v): v(_v) {}
}; 

void fn(UniformTaskCollection *coll, void *desc, int dscr_len, void *pldata,
        int pldata_len, std::vector<void*> data_bufs) {
  TskDscr *dsc = (TskDscr*)desc;
  printf("%d: val=%d\n", me, dsc->v);
  if(dsc->v>0) {
    TskDscr dsc2(dsc->v-1);
    coll->addTask(&dsc2, sizeof(dsc2));
  }
}


void the_test() {
  TslFuncRegTbl frt;
  TslFunc tf = frt.add(fn);
  TaskCollProps props;
  props.functions(tf,frt).taskSize(sizeof(TskDscr)).maxTasks(MAX_TASKS);
#if SPLIT_QUEUE
  if(me==0) {
    printf("SplitQueue version\n");
  }
  UniformTaskCollectionSplit utc(props);
#elif SHARED_QUEUE
  if(me==0) {
    printf("SharedQueue version\n");
  }
  UniformTaskCollectionShared utc(props);
#endif
  if(me==0) {
    TskDscr dsc(10);
    utc.addTask(&dsc, sizeof(dsc));
  }
  utc.process();  
}


int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  GA_Initialize();
  the_test();
  GA_Terminate();
  MPI_Finalize();
  return 0;
}
