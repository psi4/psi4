/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libqt/qt.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {

void dpd_file4_cache_init(void)
{
  dpd_main.file4_cache = NULL;
  dpd_main.file4_cache_most_recent = 0;
  dpd_main.file4_cache_least_recent = 1;
  dpd_main.file4_cache_lru_del = 0;
  dpd_main.file4_cache_low_del = 0;
}

void dpd_file4_cache_close(void)
{
  int dpdnum;
  struct dpd_file4_cache_entry *this_entry, *next_entry;
  dpdfile4 Outfile;

  this_entry = dpd_main.file4_cache;

  /* save the current dpd_default */  
  dpdnum = dpd_default;

  while(this_entry != NULL) {

    dpd_set_default(this_entry->dpdnum);

    /* Clean out each file4_cache entry */
    dpd_file4_init(&Outfile, this_entry->filenum, this_entry->irrep,
		   this_entry->pqnum, this_entry->rsnum, this_entry->label);

    next_entry = this_entry->next;
      
    dpd_file4_cache_del(&Outfile);
    dpd_file4_close(&Outfile);

    this_entry = next_entry;
  }

  /* return the dpd_default to its original value */
  dpd_set_default(dpdnum);
}

struct dpd_file4_cache_entry 
*dpd_file4_cache_scan(int filenum, int irrep, int pqnum, int rsnum, const char *label, int dpdnum)
{
  struct dpd_file4_cache_entry *this_entry;

#ifdef DPD_TIMER
  timer_on("file4_cache");
#endif

  this_entry = dpd_main.file4_cache;

  while(this_entry != NULL) {
    if(this_entry->filenum == filenum      &&
       this_entry->irrep == irrep          &&
       this_entry->pqnum == pqnum          &&
       this_entry->rsnum == rsnum          &&
       this_entry->dpdnum == dpdnum && 
       !strcmp(this_entry->label,label)) {
#ifdef DPD_TIMER
      timer_off("file4_cache");
#endif
      /* increment the access timers */
      dpd_main.file4_cache_most_recent++;
      this_entry->access = dpd_main.file4_cache_most_recent;
         
      /* increment the usage counter */
      this_entry->usage++;
	      
      return(this_entry);
    }
       
    this_entry = this_entry->next;
  }

#ifdef DPD_TIMER
  timer_off("file4_cache");
#endif
  return(this_entry);
}

struct dpd_file4_cache_entry *dpd_file4_cache_last(void)
{
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_main.file4_cache;

  while(this_entry !=NULL) {
      if(this_entry->next == NULL) return(this_entry);
      this_entry = this_entry->next;
    }

  return(NULL);
}

int dpd_file4_cache_add(dpdfile4 *File, unsigned int priority)
{
  int h, dpdnum;
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_file4_cache_scan(File->filenum, File->my_irrep,
				   File->params->pqnum, File->params->rsnum,
                                   File->label, File->dpdnum);

  if((this_entry != NULL && !(File->incore)) ||
     (this_entry == NULL && (File->incore))) {
      /* Either the file4 appears in the cache but incore is not set,
	 or incore is set and the file4 isn't in the cache */
      dpd_error("File4 cache add error!", stderr);
    }
  else if(this_entry != NULL && File->incore) {
      /* We already have this one in cache, but change its priority level */
      this_entry->priority = priority; 
      return 0;
    }
  else if(this_entry == NULL && !(File->incore)) { /* New cache entry */

      this_entry = (struct dpd_file4_cache_entry *) 
                       malloc(sizeof(struct dpd_file4_cache_entry));

      /* save the current dpd_default value */
      dpdnum = dpd_default;
      dpd_set_default(File->dpdnum);

      /* Read all data into core */
      this_entry->size = 0;
      for(h=0; h < File->params->nirreps; h++) {
          this_entry->size +=
              File->params->rowtot[h] * File->params->coltot[h^(File->my_irrep)];
          dpd_file4_mat_irrep_init(File, h);
          dpd_file4_mat_irrep_rd(File, h);
        }

      this_entry->dpdnum = File->dpdnum;
      this_entry->filenum = File->filenum;
      this_entry->irrep = File->my_irrep;
      this_entry->pqnum = File->params->pqnum;
      this_entry->rsnum = File->params->rsnum;
      strcpy(this_entry->label,File->label);
      this_entry->next = NULL;
      this_entry->last = dpd_file4_cache_last();

      this_entry->lock = 0;
      
      if(this_entry->last != NULL) this_entry->last->next = this_entry;
      else dpd_main.file4_cache = this_entry;

      /* increment the access timers */
      dpd_main.file4_cache_most_recent++;
      this_entry->access = dpd_main.file4_cache_most_recent;

      /* initialize the usage counter */
      this_entry->usage = 1;

      /* Set the clean flag */
      this_entry->clean = 1;

      /* Set the priority level */
      this_entry->priority = priority;
	      
      this_entry->matrix = File->matrix;

      File->incore = 1;

      /* Adjust the global cache size value */
      dpd_main.memcache += this_entry->size;

      /* return dpd_value to its original value */
      dpd_set_default(dpdnum);

      return 0;
  }

  return 0;
}

int dpd_file4_cache_del(dpdfile4 *File)
{
  int h, dpdnum;
  struct dpd_file4_cache_entry *this_entry, *next_entry, *last_entry;

  this_entry = dpd_file4_cache_scan(File->filenum, File->my_irrep,
				    File->params->pqnum, File->params->rsnum,
				    File->label, File->dpdnum);

  if((this_entry == NULL && File->incore) ||
     (this_entry != NULL && !(File->incore)) ||
     (this_entry == NULL && !(File->incore))) {
    dpd_error("File4 cache delete error!", stderr);
  }
  else {

    /* Save the current dpd_default */
    dpdnum = dpd_default;
    dpd_set_default(File->dpdnum);

    /* Unlock the entry first */
    dpd_file4_cache_unlock(File);

    File->incore = 0;

    /* Write all the data to disk and free the memory */
    for(h=0; h < File->params->nirreps; h++) {
      if(!(this_entry->clean)) dpd_file4_mat_irrep_wrt(File, h);
      dpd_file4_mat_irrep_close(File, h);
    }

    next_entry = this_entry->next;
    last_entry = this_entry->last;

    /* Adjust the global cache size value */
    dpd_main.memcache -= this_entry->size;

    /* Are we deleting the top of the tree? */
    if(this_entry == dpd_main.file4_cache) 
      dpd_main.file4_cache = next_entry;

    free(this_entry);

    /* Reassign pointers for adjacent entries in the list */
    if(next_entry != NULL) next_entry->last = last_entry;
    if(last_entry != NULL) last_entry->next = next_entry;

    /* Return the dpd_default to original value */
    dpd_set_default(dpdnum);

  }

  return 0;
}

void dpd_file4_cache_print_screen(void)
{
  int total_size=0;
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_main.file4_cache;

  fprintf(stdout, "\n\tDPD File4 Cache Listing:\n\n");
  fprintf(stdout,
	  "Cache Label            DPD File symm  pq  rs  use acc clean    pri lock size(kB)\n");
  fprintf(stdout,
	  "--------------------------------------------------------------------------------\n");
  while(this_entry != NULL) {
    fprintf(stdout,
	    "%-22s  %1d   %3d   %1d   %2d  %2d  %3d %3d    %1d  %6d   %1d  %8.1f\n",
	    this_entry->label, this_entry->dpdnum, this_entry->filenum, this_entry->irrep,
	    this_entry->pqnum, this_entry->rsnum, this_entry->usage, this_entry->access,
	    this_entry->clean, this_entry->priority,this_entry->lock,
	    (this_entry->size)*sizeof(double)/1e3);
    total_size += this_entry->size;
    this_entry = this_entry->next;
  }
  fprintf(stdout,
	  "--------------------------------------------------------------------------------\n");
  fprintf(stdout, "Total cached: %9.1f kB; MRU = %6d; LRU = %6d\n",
	  (total_size*sizeof(double))/1e3,dpd_main.file4_cache_most_recent,
	  dpd_main.file4_cache_least_recent);
  fprintf(stdout, "#LRU deletions = %6d; #Low-priority deletions = %6d\n", 
          dpd_main.file4_cache_lru_del,dpd_main.file4_cache_low_del);
  fprintf(stdout, "Core max size:  %9.1f kB\n", (dpd_main.memory)*sizeof(double)/1e3);
  fprintf(stdout, "Core used:      %9.1f kB\n", (dpd_main.memused)*sizeof(double)/1e3);
  fprintf(stdout, "Core available: %9.1f kB\n", dpd_memfree()*sizeof(double)/1e3);
  fprintf(stdout, "Core cached:    %9.1f kB\n", (dpd_main.memcache)*sizeof(double)/1e3);
  fprintf(stdout, "Locked cached:  %9.1f kB\n", (dpd_main.memlocked)*sizeof(double)/1e3);
  fprintf(stdout, "Most recent entry  = %d\n", dpd_main.file4_cache_most_recent);
  fprintf(stdout, "Least recent entry = %d\n", dpd_main.file4_cache_least_recent);
}

void dpd_file4_cache_print(FILE *out)
{
  int total_size=0;
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_main.file4_cache;

  fprintf(out, "\n\tDPD File4 Cache Listing:\n\n");
  fprintf(out,
	  "Cache Label            DPD File symm  pq  rs  use acc clean    pri lock size(kB)\n");
  fprintf(out,
	  "--------------------------------------------------------------------------------\n");
  while(this_entry != NULL) {
    fprintf(outfile,
	    "%-22s  %1d   %3d   %1d   %2d  %2d  %3d %3d    %1d  %6d   %1d  %8.1f\n",
	    this_entry->label, this_entry->dpdnum, this_entry->filenum, this_entry->irrep,
	    this_entry->pqnum, this_entry->rsnum, this_entry->usage, this_entry->access,
	    this_entry->clean, this_entry->priority, this_entry->lock,
	    (this_entry->size)*sizeof(double)/1e3);
    total_size += this_entry->size;
    this_entry = this_entry->next;
  }
  fprintf(outfile,
	  "--------------------------------------------------------------------------------\n");
  fprintf(out, "Total cached: %8.1f kB; MRU = %6d; LRU = %6d\n",
	  (total_size*sizeof(double))/1e3,dpd_main.file4_cache_most_recent,
	  dpd_main.file4_cache_least_recent);
  fprintf(out, "#LRU deletions = %6d; #Low-priority deletions = %6d\n", 
          dpd_main.file4_cache_lru_del,dpd_main.file4_cache_low_del);
  fprintf(out, "Core max size:  %9.1f kB\n", (dpd_main.memory)*sizeof(double)/1e3);
  fprintf(out, "Core used:      %9.1f kB\n", (dpd_main.memused)*sizeof(double)/1e3);
  fprintf(out, "Core available: %9.1f kB\n", dpd_memfree()*sizeof(double)/1e3);
  fprintf(out, "Core cached:    %9.1f kB\n", (dpd_main.memcache)*sizeof(double)/1e3);
  fprintf(out, "Locked cached:  %9.1f kB\n", (dpd_main.memlocked)*sizeof(double)/1e3);
  fprintf(out, "Most recent entry  = %d\n", dpd_main.file4_cache_most_recent);
  fprintf(out, "Least recent entry = %d\n", dpd_main.file4_cache_least_recent);
}

struct dpd_file4_cache_entry *dpd_file4_cache_find_lru(void)
{
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_main.file4_cache;

  if(this_entry == NULL) return(NULL);

  /* find the first unlocked entry */
  while(this_entry != NULL) {
    if(this_entry->lock) this_entry = this_entry->next;
    else break; /* Is this right? */
  }

  while(dpd_main.file4_cache_least_recent <= dpd_main.file4_cache_most_recent) {
    while(this_entry !=NULL) {
      if(this_entry->access <= dpd_main.file4_cache_least_recent && !this_entry->lock) 
	return(this_entry);
      this_entry = this_entry->next;
    }
    dpd_main.file4_cache_least_recent++;
    this_entry = dpd_main.file4_cache;
  }

  /*
  dpd_file4_cache_print(stderr);
  fprintf(stderr, "Possibly out of memory!\n");
  dpd_error("Error locating file4_cache LRU!", stderr);
  */
  return(NULL);
}

int dpd_file4_cache_del_lru(void)
{
  int dpdnum;
  dpdfile4 File;
  struct dpd_file4_cache_entry *this_entry;

#ifdef DPD_TIMER
  timer_on("cache_lru");
#endif

  this_entry = dpd_file4_cache_find_lru();

  if(this_entry == NULL) {
#ifdef DPD_TIMER
      timer_off("cache_lru");
#endif
      return 1; /* there is no cache or all entries are locked */
    }
  else { /* we found the LRU so delete it */
#ifdef DPD_DEBUG
      printf("Deleteing LRU: %-22s %3d %2d %2d %6d %1d %6d %8.1f\n", 
             this_entry->label,
	     this_entry->filenum, this_entry->pqnum, this_entry->rsnum,
	     this_entry->usage,this_entry->clean,this_entry->priority,
             (this_entry->size*sizeof(double))/1e3);
#endif

      /* increment the global LRU deletion counter */
      dpd_main.file4_cache_lru_del++;

      /* Save the current dpd_default */
      dpdnum = dpd_default;
      dpd_set_default(this_entry->dpdnum);
      
      dpd_file4_init(&File, this_entry->filenum, this_entry->irrep,
		     this_entry->pqnum, this_entry->rsnum, this_entry->label);

      dpd_file4_cache_del(&File);
      dpd_file4_close(&File);

      /* Return the default DPD to its original value */
      dpd_set_default(dpdnum);

#ifdef DPD_TIMER
  timer_off("cache_lru");
#endif

      return 0;
    }
}

void dpd_file4_cache_dirty(dpdfile4 *File)
{
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_file4_cache_scan(File->filenum, File->my_irrep,
                                    File->params->pqnum, File->params->rsnum,
                                    File->label, File->dpdnum);

  if((this_entry == NULL && File->incore) ||
     (this_entry != NULL && !File->incore) ||
     (this_entry == NULL && !File->incore))
     dpd_error("Error setting file4_cache dirty flag!", stderr);
  else {
     this_entry->clean = 0;
    }
}

int dpd_file4_cache_get_priority(dpdfile4 *File)
{
  struct dpd_file4_cache_entry *this_entry;
  
  this_entry = dpd_main.file4_cache_priority;

  while(this_entry != NULL) {
      if(this_entry->filenum == File->filenum     &&
         this_entry->irrep == File->my_irrep      &&
         this_entry->pqnum == File->params->pqnum &&
         this_entry->rsnum == File->params->rsnum &&
         !strcmp(this_entry->label,File->label)) 
              return(this_entry->priority);

      this_entry = this_entry->next;
    }

  return(0);
}

struct dpd_file4_cache_entry *dpd_file4_cache_find_low(void)
{
  struct dpd_file4_cache_entry *this_entry, *low_entry;

  this_entry = dpd_main.file4_cache;

  if(this_entry == NULL) return(NULL);

  /* find the first unlocked entry */
  while(this_entry != NULL) {
    if(this_entry->lock) this_entry = this_entry->next;
    else break; /* Is this right? */
  }

  /* Now search for the lowest priority entry */
  low_entry = this_entry;
  while(this_entry != NULL && low_entry != NULL) {
    if((this_entry->priority < low_entry->priority) && !this_entry->lock)
      low_entry = this_entry;
    this_entry = this_entry->next;
  }

  return low_entry;
}

int dpd_file4_cache_del_low(void)
{
  int dpdnum;
  dpdfile4 File;
  struct dpd_file4_cache_entry *this_entry;

#ifdef DPD_TIMER
  timer_on("cache_low");
#endif

  this_entry = dpd_file4_cache_find_low();

  if(this_entry == NULL) {
#ifdef DPD_TIMER
      timer_off("cache_low");
#endif
      return 1; /* there is no cache or everything is locked */
    }
  else { /* we found the LOW so delete it */
#ifdef DPD_DEBUG
      printf("Delete LOW: %-22s %3d %2d %2d %6d %1d %6d %8.1f\n",
             this_entry->label,
             this_entry->filenum, this_entry->pqnum, this_entry->rsnum,
             this_entry->usage,this_entry->clean,this_entry->priority,
             (this_entry->size*sizeof(double))/1e3);
#endif

      /* increment the global LOW deletion counter */
      dpd_main.file4_cache_low_del++;

      /* save the current dpd default value */
      dpdnum = dpd_default;

      dpd_set_default(this_entry->dpdnum);

      dpd_file4_init(&File, this_entry->filenum, this_entry->irrep,
                     this_entry->pqnum, this_entry->rsnum, this_entry->label);
      dpd_file4_cache_del(&File);
      dpd_file4_close(&File);

      /* return the default dpd to its original value */
      dpd_set_default(dpdnum);

#ifdef DPD_TIMER
  timer_off("cache_low");
#endif

      return 0;
    }
}

void dpd_file4_cache_lock(dpdfile4 *File)
{
  int h;
  struct dpd_file4_cache_entry *this_entry;

  this_entry = dpd_file4_cache_scan(File->filenum, File->my_irrep,
                                    File->params->pqnum, File->params->rsnum,
                                    File->label, File->dpdnum);

  if(this_entry != NULL && !this_entry->lock) {

    /* Increment the locked cache memory counter */
    for(h=0; h < File->params->nirreps; h++) {
      dpd_main.memlocked += File->params->rowtot[h] * File->params->coltot[h^(File->my_irrep)];
    }

    this_entry->lock = 1;
  }
}

void dpd_file4_cache_unlock(dpdfile4 *File)
{ 
  int h;
  struct dpd_file4_cache_entry *this_entry;
  
  this_entry = dpd_file4_cache_scan(File->filenum, File->my_irrep,
                                    File->params->pqnum, File->params->rsnum,
                                    File->label, File->dpdnum);
  
  if(this_entry != NULL && this_entry->lock) {

    this_entry->lock = 0;

    /* Decrement the locked cache memory counter */
    for(h=0; h < File->params->nirreps; h++) {
      dpd_main.memlocked -= File->params->rowtot[h] * File->params->coltot[h^(File->my_irrep)];
    }

  }
}


} // namespace psi
