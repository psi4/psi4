/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libqt/qt.h"
#include "dpd.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

void DPD::file4_cache_init(void)
{
    dpd_main.file4_cache = NULL;
    dpd_main.file4_cache_most_recent = 0;
    dpd_main.file4_cache_least_recent = 1;
    dpd_main.file4_cache_lru_del = 0;
    dpd_main.file4_cache_low_del = 0;
}

void DPD::file4_cache_close(void)
{
    int dpdnum;
    dpd_file4_cache_entry *this_entry, *next_entry;
    dpdfile4 Outfile;

    this_entry = dpd_main.file4_cache;

    /* save the current dpd_default */
    dpdnum = dpd_default;

    while(this_entry != NULL) {

        dpd_set_default(this_entry->dpdnum);

        /* Clean out each file4_cache entry */
        file4_init(&Outfile, this_entry->filenum, this_entry->irrep,
                   this_entry->pqnum, this_entry->rsnum, this_entry->label);

        next_entry = this_entry->next;

        file4_cache_del(&Outfile);
        file4_close(&Outfile);

        this_entry = next_entry;
    }

    /* return the dpd_default to its original value */
    dpd_set_default(dpdnum);
}

dpd_file4_cache_entry*
DPD::file4_cache_scan(int filenum, int irrep, int pqnum, int rsnum, const char *label, int dpdnum)
{
    dpd_file4_cache_entry *this_entry;

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

dpd_file4_cache_entry*
DPD::file4_cache_last(void)
{
    dpd_file4_cache_entry *this_entry;

    this_entry = dpd_main.file4_cache;

    while(this_entry !=NULL) {
        if(this_entry->next == NULL) return(this_entry);
        this_entry = this_entry->next;
    }

    return(NULL);
}

int DPD::file4_cache_add(dpdfile4 *File, unsigned int priority)
{
    int h, dpdnum;
    dpd_file4_cache_entry *this_entry;

    this_entry = file4_cache_scan(File->filenum, File->my_irrep,
                                  File->params->pqnum, File->params->rsnum,
                                  File->label, File->dpdnum);

    if((this_entry != NULL && !(File->incore)) ||
            (this_entry == NULL && (File->incore))) {
        /* Either the file4 appears in the cache but incore is not set,
     or incore is set and the file4 isn't in the cache */
        dpd_error("File4 cache add error!", "outfile");
    }
    else if(this_entry != NULL && File->incore) {
        /* We already have this one in cache, but change its priority level */
        this_entry->priority = priority;
        return 0;
    }
    else if(this_entry == NULL && !(File->incore)) { /* New cache entry */

        this_entry = (dpd_file4_cache_entry *)
                malloc(sizeof(dpd_file4_cache_entry));

        /* save the current dpd_default value */
        dpdnum = dpd_default;
        dpd_set_default(File->dpdnum);

        /* Read all data into core */
        this_entry->size = 0;
        for(h=0; h < File->params->nirreps; h++) {
            this_entry->size +=
                    File->params->rowtot[h] * File->params->coltot[h^(File->my_irrep)];
            file4_mat_irrep_init(File, h);
            file4_mat_irrep_rd(File, h);
        }

        this_entry->dpdnum = File->dpdnum;
        this_entry->filenum = File->filenum;
        this_entry->irrep = File->my_irrep;
        this_entry->pqnum = File->params->pqnum;
        this_entry->rsnum = File->params->rsnum;
        strcpy(this_entry->label,File->label);
        this_entry->next = NULL;
        this_entry->last = file4_cache_last();

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

int DPD::file4_cache_del(dpdfile4 *File)
{
    int h, dpdnum;
    dpd_file4_cache_entry *this_entry, *next_entry, *last_entry;

    this_entry = file4_cache_scan(File->filenum, File->my_irrep,
                                  File->params->pqnum, File->params->rsnum,
                                  File->label, File->dpdnum);

    if((this_entry == NULL && File->incore) ||
            (this_entry != NULL && !(File->incore)) ||
            (this_entry == NULL && !(File->incore))) {
        dpd_error("File4 cache delete error!", "outfile");
    }
    else {

        /* Save the current dpd_default */
        dpdnum = dpd_default;
        dpd_set_default(File->dpdnum);

        /* Unlock the entry first */
        file4_cache_unlock(File);

        File->incore = 0;

        /* Write all the data to disk and free the memory */
        for(h=0; h < File->params->nirreps; h++) {
            if(!(this_entry->clean)) file4_mat_irrep_wrt(File, h);
            file4_mat_irrep_close(File, h);
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

void DPD::file4_cache_print_screen(void)
{
    int total_size=0;
    dpd_file4_cache_entry *this_entry;

    this_entry = dpd_main.file4_cache;

    outfile->Printf("\n\tDPD File4 Cache Listing:\n\n");
    outfile->Printf(
            "Cache Label            DPD File symm  pq  rs  use acc clean    pri lock size(kB)\n");
    outfile->Printf(
            "--------------------------------------------------------------------------------\n");
    while(this_entry != NULL) {
        outfile->Printf(
                "%-22s  %1d   %3d   %1d   %2d  %2d  %3d %3d    %1d  %6d   %1d  %8.1f\n",
                this_entry->label, this_entry->dpdnum, this_entry->filenum, this_entry->irrep,
                this_entry->pqnum, this_entry->rsnum, this_entry->usage, this_entry->access,
                this_entry->clean, this_entry->priority,this_entry->lock,
                (this_entry->size)*sizeof(double)/1e3);
        total_size += this_entry->size;
        this_entry = this_entry->next;
    }
    outfile->Printf(
            "--------------------------------------------------------------------------------\n");
    outfile->Printf( "Total cached: %9.1f kB; MRU = %6d; LRU = %6d\n",
            (total_size*sizeof(double))/1e3,dpd_main.file4_cache_most_recent,
            dpd_main.file4_cache_least_recent);
    outfile->Printf( "#LRU deletions = %6d; #Low-priority deletions = %6d\n",
            dpd_main.file4_cache_lru_del,dpd_main.file4_cache_low_del);
    outfile->Printf( "Core max size:  %9.1f kB\n", (dpd_main.memory)*sizeof(double)/1e3);
    outfile->Printf( "Core used:      %9.1f kB\n", (dpd_main.memused)*sizeof(double)/1e3);
    outfile->Printf( "Core available: %9.1f kB\n", dpd_memfree()*sizeof(double)/1e3);
    outfile->Printf( "Core cached:    %9.1f kB\n", (dpd_main.memcache)*sizeof(double)/1e3);
    outfile->Printf( "Locked cached:  %9.1f kB\n", (dpd_main.memlocked)*sizeof(double)/1e3);
    outfile->Printf( "Most recent entry  = %d\n", dpd_main.file4_cache_most_recent);
    outfile->Printf( "Least recent entry = %d\n", dpd_main.file4_cache_least_recent);
}

void DPD::file4_cache_print(std::string out)
{
    int total_size=0;
    std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
             std::shared_ptr<OutFile>(new OutFile(out)));
    dpd_file4_cache_entry *this_entry;

    this_entry = dpd_main.file4_cache;

    printer->Printf( "\n\tDPD File4 Cache Listing:\n\n");
    printer->Printf(
            "Cache Label            DPD File symm  pq  rs  use acc clean    pri lock size(kB)\n");
    printer->Printf(
            "--------------------------------------------------------------------------------\n");
    while(this_entry != NULL) {
        printer->Printf(
                "%-22s  %1d   %3d   %1d   %2d  %2d  %3d %3d    %1d  %6d   %1d  %8.1f\n",
                this_entry->label, this_entry->dpdnum, this_entry->filenum, this_entry->irrep,
                this_entry->pqnum, this_entry->rsnum, this_entry->usage, this_entry->access,
                this_entry->clean, this_entry->priority, this_entry->lock,
                (this_entry->size)*sizeof(double)/1e3);
        total_size += this_entry->size;
        this_entry = this_entry->next;
    }
    printer->Printf(
            "--------------------------------------------------------------------------------\n");
    printer->Printf( "Total cached: %8.1f kB; MRU = %6d; LRU = %6d\n",
            (total_size*sizeof(double))/1e3,dpd_main.file4_cache_most_recent,
            dpd_main.file4_cache_least_recent);
    printer->Printf( "#LRU deletions = %6d; #Low-priority deletions = %6d\n",
            dpd_main.file4_cache_lru_del,dpd_main.file4_cache_low_del);
    printer->Printf( "Core max size:  %9.1f kB\n", (dpd_main.memory)*sizeof(double)/1e3);
    printer->Printf( "Core used:      %9.1f kB\n", (dpd_main.memused)*sizeof(double)/1e3);
    printer->Printf( "Core available: %9.1f kB\n", dpd_memfree()*sizeof(double)/1e3);
    printer->Printf( "Core cached:    %9.1f kB\n", (dpd_main.memcache)*sizeof(double)/1e3);
    printer->Printf( "Locked cached:  %9.1f kB\n", (dpd_main.memlocked)*sizeof(double)/1e3);
    printer->Printf( "Most recent entry  = %d\n", dpd_main.file4_cache_most_recent);
    printer->Printf( "Least recent entry = %d\n", dpd_main.file4_cache_least_recent);
}

dpd_file4_cache_entry*
DPD::file4_cache_find_lru(void)
{
    dpd_file4_cache_entry *this_entry;

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
  dpd_file4_cache_print("outfile");
  outfile->Printf( "Possibly out of memory!\n");
  dpd_error("Error locating file4_cache LRU!", "outfile");
  */
    return(NULL);
}

int DPD::file4_cache_del_lru(void)
{
    int dpdnum;
    dpdfile4 File;
    dpd_file4_cache_entry *this_entry;

#ifdef DPD_TIMER
    timer_on("cache_lru");
#endif

    this_entry = file4_cache_find_lru();

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

        file4_init(&File, this_entry->filenum, this_entry->irrep,
                   this_entry->pqnum, this_entry->rsnum, this_entry->label);

        file4_cache_del(&File);
        file4_close(&File);

        /* Return the default DPD to its original value */
        dpd_set_default(dpdnum);

#ifdef DPD_TIMER
        timer_off("cache_lru");
#endif

        return 0;
    }
}

void DPD::file4_cache_dirty(dpdfile4 *File)
{
    dpd_file4_cache_entry *this_entry;

    this_entry = file4_cache_scan(File->filenum, File->my_irrep,
                                  File->params->pqnum, File->params->rsnum,
                                  File->label, File->dpdnum);

    if((this_entry == NULL && File->incore) ||
            (this_entry != NULL && !File->incore) ||
            (this_entry == NULL && !File->incore))
        dpd_error("Error setting file4_cache dirty flag!", "outfile");
    else {
        this_entry->clean = 0;
    }
}

int DPD::file4_cache_get_priority(dpdfile4 *File)
{
    dpd_file4_cache_entry *this_entry;

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

dpd_file4_cache_entry*
dpd_file4_cache_find_low(void)
{
    dpd_file4_cache_entry *this_entry, *low_entry;

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

int DPD::file4_cache_del_low(void)
{
    int dpdnum;
    dpdfile4 File;
    dpd_file4_cache_entry *this_entry;

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

        file4_init(&File, this_entry->filenum, this_entry->irrep,
                   this_entry->pqnum, this_entry->rsnum, this_entry->label);
        file4_cache_del(&File);
        file4_close(&File);

        /* return the default dpd to its original value */
        dpd_set_default(dpdnum);

#ifdef DPD_TIMER
        timer_off("cache_low");
#endif

        return 0;
    }
}

void DPD::file4_cache_lock(dpdfile4 *File)
{
    int h;
    dpd_file4_cache_entry *this_entry;

    this_entry = file4_cache_scan(File->filenum, File->my_irrep,
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

void DPD::file4_cache_unlock(dpdfile4 *File)
{
    int h;
    dpd_file4_cache_entry *this_entry;

    this_entry = file4_cache_scan(File->filenum, File->my_irrep,
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


}
