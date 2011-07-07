/*
 $Id: crayio.c 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
   Simulate the cray 64 bit word addressable I/O routines and
   allow for large buffering in core.


   CALL WOPEN(UNIT, NAME, BLOCKS, STATS, IERR)
                    ----

   CALL WCLOSE(UNIT, IERR)

   CALL GETWA(UNIT, RESULT, ADDR, COUNT, IERR)

   CALL PUTWA(UNIT, SOURCE, ADDR, COUNT, IERR)

   Currently the I/O is syncronous and unbuffered */

#include <stdio.h>

#include <sys/time.h>
#include <sys/types.h>
#include <sys/file.h>
#include <fcntl.h>

typedef int integer;

#ifndef L_XTND
#define L_XTND SEEK_END
#endif

#define max_file 99

extern char *malloc();
extern off_t lseek();
extern char *strncpy();
/*extern char *sprintf();*/

static int first_call = 1;   /* need to do stuff on the first call */

static struct w_file {
  int fds;                   /* file descriptor */
  off_t length;              /* file length in bytes */
  off_t position;            /* current file position in bytes a la lseek */
  char *path;                /* file name */
  int stats;                 /* boolean flag to collect statistics */
  double words_write;        /* total no. of words written */
  double words_read;         /* total no. of words read */
  double time_write;         /* total wall time writing */
  double time_read;          /* total wall time reading */
  int n_read;                /* no. of read requests */
  int n_write;               /* no. of write reqeusts */
  int seek_read;             /* no. of seeks on read */
  int seek_write;            /* no. of seeks on write */
} file_array[max_file];

void walltm_(double *ai)
/* return ai with the wall clock time in seconds as a double.
   it might be accurate to about 0.01s at best */
{
  struct timeval tp;
  struct timezone tzp;

  (void) gettimeofday(&tp,&tzp);
  *ai = (double) tp.tv_sec + ((double) tp.tv_usec) * 1.0e-6;
}

static int CheckUnit(integer unit)
{
  if ( (unit < 0) || (unit >= max_file) )
    return -1;

  if ( file_array[unit].fds == -1 )
    return -1;

  return 0;
}

static int CheckAddr(integer addr)
{
  if (addr <= 0)
    return -4;
  else
    return 0;
}

static int CheckCount(integer count)
{
  if (count < 0)
    return -4;
  else
    return 0;
}

void InitFileStats(struct w_file *file)
{
  file->stats = 1;
  file->words_write = 0.0e0;
  file->words_read = 0.0e0;
  file->time_read = 0.0e0;
  file->time_write = 0.0e0;
  file->n_read = 0;
  file->n_write = 0;
  file->seek_read = 0;
  file->seek_write = 0;
}

void PrintFileStats(struct w_file *unit, struct w_file *file)
{
  double ave_read=0.0e0, ave_write=0.0e0;
  double rate_read=0.0e0, rate_write=0.0e0;

  (void) fflush(stdout);

  if (file->n_read) {
    ave_read = file->words_read / (double) file->n_read;
    if (file->time_read > 0.0e0)
      rate_read = file->words_read / (1000000.0e0 * file->time_read);
  }

  if (file->n_write) {
    ave_write = file->words_write / (double) file->n_write;
    if (file->time_write > 0.0e0)
      rate_write = file->words_write / (1000000.0e0 * file->time_write);
  }

  (void) fflush(stdout);
  (void) fflush(stderr);
  (void) fprintf(stdout,"CRAYIO: Statistics for unit %d, file '%s', length=%d bytes.\n",
		 unit, file->path, file->length);
  (void) fprintf(stdout,"CRAYIO: oper :  #req.  :  #seek  :   #words  :");
  (void) fprintf(stdout," #w/#req : time(s) :  MW/s \n");
  (void) fprintf(stdout,"CRAYIO: read : %7d : %7d : %9d : %7d : %7.1f : %6.3f\n",
		 file->n_read, file->seek_read, (long) file->words_read,
		 (long) ave_read, file->time_read, rate_read);
  (void) fprintf(stdout,"CRAYIO:write : %7d : %7d : %9d : %7d : %7.1f : %6.3f\n",
		 file->n_write, file->seek_write, (long) file->words_write,
		 (long) ave_write, file->time_write, rate_write);
  (void) fflush(stdout);
}

void InitFileData(struct w_file *file)
{
  file->fds = -1;
  file->length = -1;
  file->path = (char *) NULL;
  file->position = (off_t) -1;
}

void FirstCall()
     /* Initialization on first call to anything */
{
  int i;

  for (i=0; i<max_file; ++i) {

    InitFileData(&file_array[i]);

    InitFileStats(&file_array[i]);
  }

  first_call = 0;
}

void wclose_(integer* unit, integer* ierr)
{
  struct w_file *file;

  if (first_call)
    FirstCall();

  if (*ierr = CheckUnit(*unit))
    return;

  file = file_array + *unit;

  *ierr = close(file->fds);

  if (file->stats)
    PrintFileStats(*unit, file);

  InitFileData(file);

  InitFileStats(file);
}

/* ARGSUSED */
void wopen_(integer* unit, char* name, integer* blocks, integer* stats, integer* ierr, int lennam)
{
  struct w_file *file;

  *ierr = (integer) 0;

  if (first_call)
    FirstCall();

  if ( (*unit < 0) || (*unit >= max_file) ) {
    *ierr = -1;
    return;
  }

  file = file_array + *unit;

  file->stats = *stats;

  if (lennam > 0) {
    file->path = malloc((unsigned) (lennam + 1));
    (void) strncpy(file->path,name,lennam);
    file->path[lennam] = 0;
    }
  else {
    file->path = malloc((unsigned) 8);
    (void) sprintf(file->path,"fort.%.2d",*unit);
  }

  if (( file->fds = open(file->path, (int)(O_RDWR|O_CREAT), (int) 0660))
      == -1) {
    *ierr = -6;
    return;
  }

  file->length = lseek(file->fds, (off_t) 0, (int) L_XTND);
  file->position = lseek(file->fds, (off_t) 0, (int) L_SET);

}

void getwa_(integer* unit, double* result, integer* addr, integer* count, integer* ierr)
{
  long nbytes;
  off_t where;
  double start, end;
  struct w_file *file;

  if (first_call)
    FirstCall();

  if (*ierr = CheckUnit(*unit))
    return;

  if (*ierr = CheckAddr(*addr))
    return;

  if (*ierr = CheckCount(*count))
    return;

  file = file_array + *unit;

  nbytes = *count * 8;
  where = (*addr - 1) * 8;
  if ( (where+nbytes) > file->length ) {
    *ierr = -5;
    return;
  }

  if (file->stats)
    walltm_(&start);

  if (where != file->position) {
    ++(file->seek_read);
    if ( (file->position = lseek(file->fds, where, L_SET)) == (off_t) -1) {
      *ierr = -4;
      return;
    }
  }

  if ((long) read(file->fds, (char *) result, (int) nbytes) != nbytes) {
    *ierr = -6;
    return;
  }

  file->position += nbytes;

  if (file->stats) {
    walltm_(&end);
    ++(file->n_read);
    file->words_read += (double) *count;
    file->time_read +=  end - start;
  }

  *ierr = 0;
}

void putwa_(integer* unit, double* source, integer* addr, integer* count, integer* ierr)
{
  long nbytes;
  off_t where;
  double start, end;
  struct w_file *file;

  if (first_call)
    FirstCall();

  if ( *ierr = CheckUnit(*unit))
    return;

  if (*ierr = CheckAddr(*addr))
    return;

  if (*ierr = CheckCount(*count))
    return;

  file = file_array + *unit;

  nbytes = *count * 8;
  where = (*addr - 1) * 8;

  if (file->stats)
    walltm_(&start);

  if (where != file->position) {
    ++(file->seek_write);
    if ( (file->position = lseek(file->fds, where, L_SET)) == (off_t) -1) {
      *ierr = -4;
      return;
    }
  }

  if ((long) write(file->fds, (char *) source, (int) nbytes) != nbytes) {
    *ierr = -6;
    return;
  }

  where += nbytes;
  file->position += nbytes;
  if (file->length < where)
    file->length = where;


  if (file->stats) {
    walltm_(&end);
    ++(file->n_write);
    file->words_write += (double) *count;
    file->time_write +=  end - start;
  }

  *ierr = 0;
}

