/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {

void dpd_file2_cache_init(void)
{
  dpd_main.file2_cache = NULL;
}

void dpd_file2_cache_close(void)
{
  int dpdnum;
  struct dpd_file2_cache_entry *this_entry, *next_entry;
  dpdfile2 Outfile;

  this_entry = dpd_main.file2_cache;

  dpdnum = dpd_default;  

  while(this_entry != NULL) {

    dpd_set_default(this_entry->dpdnum);

    /* Clean out each file2_cache entry */
    dpd_file2_init(&Outfile, this_entry->filenum, this_entry->irrep,
		   this_entry->pnum, this_entry->qnum, this_entry->label);

    next_entry = this_entry->next;
      
    dpd_file2_cache_del(&Outfile);
    dpd_file2_close(&Outfile);

    this_entry = next_entry;
  }

  dpd_set_default(dpdnum);
}

struct dpd_file2_cache_entry 
*dpd_file2_cache_scan(int filenum, int irrep, int pnum, int qnum, const char *label, int dpdnum)
{
  struct dpd_file2_cache_entry *this_entry;

  this_entry = dpd_main.file2_cache;

  while(this_entry != NULL) {
    if(this_entry->filenum == filenum       &&
       this_entry->irrep == irrep           &&
       this_entry->pnum == pnum             &&
       this_entry->qnum == qnum             &&
       this_entry->dpdnum == dpdnum         &&
       !strcmp(this_entry->label,label)) return(this_entry);
       
    this_entry = this_entry->next;
  }

  return(this_entry);
}

struct dpd_file2_cache_entry *dpd_file2_cache_last(void)
{
  struct dpd_file2_cache_entry *this_entry;

  this_entry = dpd_main.file2_cache;

  while(this_entry !=NULL) {
      if(this_entry->next == NULL) return(this_entry);
      this_entry = this_entry->next;
    }

  return(NULL);
}

int dpd_file2_cache_add(dpdfile2 *File)
{
  int h, dpdnum;
  struct dpd_file2_cache_entry *this_entry;

  if(File->incore) return 0; /* Already have this one in cache */

  this_entry = dpd_file2_cache_scan(File->filenum, File->my_irrep,
				    File->params->pnum, File->params->qnum,
				    File->label, File->dpdnum);

  if(this_entry == NULL) { /* New cache entry */
    this_entry = (struct dpd_file2_cache_entry *) 
      malloc(sizeof(struct dpd_file2_cache_entry));

    dpdnum = dpd_default;
    dpd_set_default(File->dpdnum);

    this_entry->dpdnum = File->dpdnum;
    this_entry->filenum = File->filenum;
    this_entry->irrep = File->my_irrep;
    this_entry->pnum = File->params->pnum;
    this_entry->qnum = File->params->qnum;
    strcpy(this_entry->label,File->label);
    this_entry->next = NULL;
    this_entry->last = dpd_file2_cache_last();
      
    if(this_entry->last != NULL) this_entry->last->next = this_entry;
    else dpd_main.file2_cache = this_entry;

    this_entry->size = 0;
    for(h=0; h < File->params->nirreps; h++)
      this_entry->size += 
	File->params->rowtot[h] * File->params->coltot[h^File->my_irrep];

    /* Read all data into core */
    dpd_file2_mat_init(File);
    dpd_file2_mat_rd(File);

    this_entry->clean = 1;

    this_entry->matrix = File->matrix;

    File->incore = 1;

    dpd_set_default(dpdnum);

    return 0;
  }

  /* The Buffer appears in the cache, but incore is not set */
  dpd_error("File2 cache add error!", stderr);
  
  return 0;
}

int dpd_file2_cache_del(dpdfile2 *File)
{
  int dpdnum;
  struct dpd_file2_cache_entry *this_entry, *next_entry, *last_entry;

  /* The input buffer isn't in the cache! */
  if(!File->incore) dpd_error("File2 cache delete error!", stderr);

  this_entry = dpd_file2_cache_scan(File->filenum, File->my_irrep,
				    File->params->pnum, File->params->qnum,
				    File->label, File->dpdnum);


  if(this_entry == NULL) dpd_error("File2 cache delete error!", stderr);
  else {
    File->incore = 0;

    dpdnum = dpd_default;
    dpd_set_default(File->dpdnum);
      
    /* Write all the data to disk and free the memory */
    if(!(this_entry->clean)) dpd_file2_mat_wrt(File);
    dpd_file2_mat_close(File);

    next_entry = this_entry->next;
    last_entry = this_entry->last;

    /* Are we deleting the top of the tree? */
    if(this_entry == dpd_main.file2_cache) 
      dpd_main.file2_cache = next_entry;
	  
    free(this_entry);

    /* Reassign pointers for adjacent entries in the list */
    if(next_entry != NULL) next_entry->last = last_entry;
    if(last_entry != NULL) last_entry->next = next_entry;

    dpd_set_default(dpdnum);
  }

  return 0;
}

void dpd_file2_cache_print(FILE *out)
{
  int total_size=0;
  struct dpd_file2_cache_entry *this_entry;

  this_entry = dpd_main.file2_cache;

  fprintf(out, "\n\tDPD File2 Cache Listing:\n\n");
  fprintf(out,
    "Cache Label                     File symm  p  q  size(kB)\n");
  fprintf(out,
    "---------------------------------------------------------\n");
  while(this_entry != NULL) {
      fprintf(out,
      "%-32s %3d    %1d  %1d  %1d  %8.1f\n",
	      this_entry->label, this_entry->filenum, this_entry->irrep,
	      this_entry->pnum, this_entry->qnum,
	      (this_entry->size)*sizeof(double)/1e3);
      total_size += this_entry->size;
      this_entry = this_entry->next;
    }
  fprintf(out,
    "---------------------------------------------------------\n");
  fprintf(out, "Total cached: %8.1f kB\n", total_size*sizeof(double)/1e3);
}

void dpd_file2_cache_dirty(dpdfile2 *File)
{
  struct dpd_file2_cache_entry *this_entry;

  this_entry = dpd_file2_cache_scan(File->filenum, File->my_irrep,
                                    File->params->pnum, File->params->qnum,
                                    File->label, File->dpdnum);

  if((this_entry == NULL && File->incore) ||
     (this_entry != NULL && !File->incore) ||
     (this_entry == NULL && !File->incore))
     dpd_error("Error setting file4_cache dirty flag!", stderr);
  else {
     this_entry->clean = 0;
    }
}

} // namespace psi
