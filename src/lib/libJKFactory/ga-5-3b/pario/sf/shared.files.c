#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * $Id: shared.files.c,v 1.14 2002-11-09 06:05:40 sohirata Exp $
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with United States
 * Government support under Contract Number DE-AC06-76RLO-1830 awarded by
 * the United States Department of Energy.  The United States Government
 * retains a paid-up non-exclusive, irrevocable worldwide license to
 * reproduce, prepare derivative works, perform publicly and display
 * publicly by or for the US Government, including the right to
 * distribute to other US Government contractors.
 */
#if HAVE_STRING_H
#   include <string.h>
#endif
#include "ga.h"
#include "elio.h"
#include "sf.h"
#include "sff2c.h"

#define _max_shared_files 100
#define SF_OFFSET 3000
#define SF_FAIL (Integer)1

typedef struct{
    Integer handle;
    Integer actv;
    SFsize_t soft_size; 
    SFsize_t hard_size; 
    Fd_t fd;
    char fname[200];
} SF_t;

SF_t SF[_max_shared_files];

#include "coms.h"

/******************** internal macros ********************************/

#define sfi_check_handleM(_handle, msg)\
{\
    if(_handle+SF_OFFSET >=_max_shared_files||_handle+SF_OFFSET<0) { \
        fprintf(stderr,"%s, %ld --",msg,(long) _max_shared_files); \
        ERROR("invalid SF handle",(int)_handle); \
    } \
    if( SF[(_handle+SF_OFFSET)].actv == 0) { \
        fprintf(stderr,"%s:",msg); \
        ERROR("disk array not active",(int)_handle); \
    } \
}

/*********************************************************************/

Integer sfi_get_handle()
{
    Integer sf_handle =-1, candidate = 0;

    do{
        if(!SF[candidate].actv){
            sf_handle=candidate;
            SF[candidate].actv =1;
        }
        candidate++;
    }while(candidate < _max_shared_files && sf_handle == -1);
    return(sf_handle);
}


void sfi_release_handle(Integer *handle)
{
    SF[*handle+SF_OFFSET].actv =0;
    *handle = 0;
}


Integer sfi_create(char *fname, SFsize_t *size_hard_limit,
        SFsize_t *size_soft_limit, SFsize_t *req_size, Integer *handle)
{
    Integer hndl;
    SYNC();

    /*** Get next free SF handle ***/
    if( (hndl = sfi_get_handle()) == -1)
        ERROR("sf_create: too many shared files ",(int)_max_shared_files);
    *handle = hndl - SF_OFFSET;

    /*generate file name(s) */
    sprintf(SF[hndl].fname,"%s.%d",fname, (int)hndl);

    if (ME() == 0) SF[hndl].fd = elio_open(SF[hndl].fname,ELIO_RW, ELIO_SHARED);
    SYNC();
    if (ME() != 0) SF[hndl].fd = elio_open(SF[hndl].fname,ELIO_RW, ELIO_SHARED);

    if(SF[hndl].fd==NULL) ERROR("sf_create: could not open file",0);
    if(SF[hndl].fd->fd==-1) ERROR("sf_create: descriptor -1",0);

    SF[hndl].soft_size = *size_soft_limit;
    SF[hndl].hard_size = *size_hard_limit;
    SYNC();
    return(ELIO_OK);
}


/**
 * takes an additional input suffix which is appended to the filename instead
 * of the handle. This will give users complete control over the filename as it
 * is written to and read from disk. 
 */
Integer sfi_create_suffix(char *fname, SFsize_t *size_hard_limit,
        SFsize_t *size_soft_limit, SFsize_t *req_size,
        Integer *handle, Integer *suffix)
{
    Integer hndl;

    SYNC();

    /*** Get next free SF handle ***/
    if( (hndl = sfi_get_handle()) == -1) 
    {
        ERROR("sf_create_suffix: too many shared files ",
                (int)_max_shared_files);
    }

    *handle = hndl - SF_OFFSET;

    /*generate file name(s) */
    sprintf(SF[hndl].fname,"%s.%d",fname, (int) *suffix);

    if (ME() == 0) SF[hndl].fd = elio_open(SF[hndl].fname,ELIO_RW, ELIO_SHARED);
    SYNC();
    if (ME() != 0) SF[hndl].fd = elio_open(SF[hndl].fname,ELIO_RW, ELIO_SHARED);

    if(SF[hndl].fd==NULL) ERROR("sf_create_suffix: could not open file",0);
    if(SF[hndl].fd->fd==-1) ERROR("sf_create_suffix: descriptor -1",0);

    SF[hndl].soft_size = *size_soft_limit;
    SF[hndl].hard_size = *size_hard_limit;
    SYNC();

    return(ELIO_OK);
}

/*****************************************************************************/

Integer FATR sf_destroy_(Integer *s_a /**< [in] SF handle */)
{
    Integer handle = *s_a+SF_OFFSET;

    SYNC();

    sfi_check_handleM(*s_a,"sf_delete");

    elio_close(SF[handle].fd); /* fix from Peter Knowles */

    SYNC(); /* this sync is unnecessary under Unix */

    if(ME() == 0)elio_delete(SF[handle].fname);
    sfi_release_handle(s_a);

    SYNC();

    return(ELIO_OK);
}


/**
 * close rw file and open as read only
 */
Integer FATR sf_rwtor_(Integer *s_a /**< [in] SF handle */)
{
    Integer handle = *s_a+SF_OFFSET;

    elio_close(SF[handle].fd);
    SF[handle].fd = elio_open(SF[handle].fname,ELIO_R, ELIO_SHARED);

    if(SF[handle].fd==NULL) ERROR("sf_open: could not open file",0);
    if(SF[handle].fd->fd==-1) ERROR("sf_open: descriptor -1",0);

    return(ELIO_OK);
}


/**
 * open
 */
Integer FATR sf_open_(Integer *s_a /**< [in] SF handle */)
{
    Integer handle = *s_a+SF_OFFSET;

    SF[handle].fd = elio_open(SF[handle].fname,ELIO_RW, ELIO_SHARED);

    if(SF[handle].fd==NULL) ERROR("sf_open: could not open file",0);
    if(SF[handle].fd->fd==-1) ERROR("sf_open: descriptor -1",0);

    return(ELIO_OK);
}


/**
 * close
 */
Integer FATR sf_close_(Integer *s_a)
{
    Integer handle = *s_a+SF_OFFSET;

    elio_close(SF[handle].fd);
    return(ELIO_OK);
}

/**
 * fsync
 */
Integer FATR sf_fsync_(Integer *s_a)
{
    Integer handle = *s_a+SF_OFFSET;

    elio_fsync(SF[handle].fd);
    return(ELIO_OK);
}


/**
 * asynchronous write to shared file
 */
Integer FATR sf_write_(Integer *s_a, SFsize_t *offset, SFsize_t *bytes,
        char *buffer, Integer *req_id)
{
    Integer handle = *s_a+SF_OFFSET;
    int status;
    io_request_t id;

    sfi_check_handleM(*s_a,"sf_write");
    status = elio_awrite(SF[handle].fd, (Off_t)*offset, buffer, 
            (Size_t)*bytes, &id);
    *req_id = (Integer)id;
    return((Integer)status);
    /*
       status = elio_write(SF[handle].fd, (Off_t)*offset, buffer,(Size_t)*bytes);
     *req_id = (Integer)ELIO_DONE;
     if(status != (int)*bytes)
     return((Integer)SF_FAIL);
     else
     return((Integer)ELIO_OK);
     */
}


/**
 * asynchronous read from shared file
 */
Integer FATR sf_read_(Integer *s_a, SFsize_t *offset, SFsize_t *bytes,
        char *buffer, Integer *req_id)
{
    Integer handle = *s_a+SF_OFFSET;
    int status;
    io_request_t id;

    sfi_check_handleM(*s_a,"sf_read");
    status = elio_aread(SF[handle].fd, (Off_t)*offset, buffer, 
            (Size_t)*bytes, &id);
    *req_id = (Integer)id;
    return((Integer)status);
}


/**
 * block calling process until I/O operation completes
 */
Integer FATR sf_wait_(Integer *req_id)
{
    int status;
    io_request_t id = (io_request_t) *req_id;
    status = elio_wait(&id);
    *req_id = (Integer)id;
    return((Integer)status);
}


/**
 * block calling process until all I/O operations associated
 *  with id in the list complete
 */
Integer FATR sf_waitall_(Integer *list, Integer *num)
{
    int i;
    int status, fail=0;

    for(i=0;i<*num;i++){
        io_request_t id = (io_request_t) list[i];
        status = elio_wait(&id);
        if(status != ELIO_OK) fail++;
        list[i] = (Integer)id;
    }
    if (fail)return((Integer)SF_FAIL);
    else     return((Integer)ELIO_OK);
}


/**
 * retrieve message associated with error code
 */
void sfi_errmsg(int code, char *msg)
{
    if(code==SF_FAIL)
        (void) strcpy(msg, "SF operation failed");
    else
        elio_errmsg(code, msg);
}
