#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#define MAX_TOKEN 2048

/**
 * Read next token from file. Tokens are separated by the character
 * delimeter.
 * Tokens are returned as NUL terminated character strings which
 * almost certainly should always be freed with free() after use.
 *
 * There is an internally imposed maximum token size of MAX_TOKEN bytes.
 *
 * All errors are handled by returning a NULL pointer.
 */
char *ReadToken(FILE *file, char delimiter)
{
    char *buf = malloc((unsigned) MAX_TOKEN);
    char *temp;
    int input, used=0;

    if (!buf) {
        return (char *) NULL;
    }

    temp = buf;

    while ( (input = getc(file)) != EOF ) {

        used++;

        if (input == delimiter) {
            *temp = '\0';
            break;
        }
        else {
            *temp++ = (char) input;
            if (used == MAX_TOKEN) {
                used = 0; break;
            }
        }
    }

    /* duplicate string to minimize problems if string is not
       freed in calling program */

    if (used) {
        temp = strdup(buf);
    }
    else {
        temp = (char *) NULL;
    }

    (void) free(buf);

    return temp;
}


int FindToken(char *token, FILE *file, char delimiter)
{
    char *input;

    while (input = ReadToken(file, delimiter)) {
        if (strcmp(input, token) == 0) {
            return 1;
        }
        else {
            (void) free(input);
        }
    }

    return 0;
}


/**
 * Return the starting time and duration of the events file.
 */
void GetTimeSpan(FILE *file, DoublePrecision *start, DoublePrecision *duration)
{
    char *input;
    DoublePrecision end;

    end = *start = 0.0;

    if (FindToken("TIME", file, ':')) {
        end = *start = atof(input = ReadToken(file, ':'));
        (void) free(input);
    }

    while (FindToken("TIME", file, ':')) {
        end = atof(input = ReadToken(file, ':'));
        (void) free(input);
    }

    *duration = end - *start;

    (void) fseek(file, 0L, 0);
}
  

int main(int argc, char **argv)
{
    char filename[11];
    FILE *file, *plot;
    char *token;
    char *event;
    DoublePrecision time, start, duration=0.0, otim, span, margin, comms, useful;
    int newstate, state, i, nproc, lo, hi;

    if (argc == 1) {
        lo = 0;
        hi = 127;
    }
    else if (argc == 3) {
        lo = atoi(argv[1]);
        hi = atoi(argv[2]);
    }
    else {
        (void) fprintf(stderr, "usage: %s [lo hi]\n", argv[0]);
        (void) fprintf(stderr, "... with no arguments parse all event files\n");
        (void) fprintf(stderr, "... or with lo & hi only files in this range\n");
        (void) fprintf(stderr, "... e.g.  parse 16 31\n");
        return 1;
    }

    /* open the file that will have the plot stuff in it */
    /* change of heart here ... just write to stdout */

    /* 
       if (!(plot = fopen("plot", "w"))) {
       perror("failed to open plot file");
       return 1;
       }
       */

    plot = stdout;

    /* Determine how many processes there are and maximum time span */

    nproc = 0;
    for (i=lo; i<=hi; i++) {

        (void) sprintf(filename, "events.%03d", nproc);

        if ( !(file = fopen(filename, "r")) ) {
            break;
        }

        GetTimeSpan(file, &start, &span);
        (void) fclose(file);

        if (span > duration) {
            duration = span;
        }

        nproc++;
    }      

    margin = duration * 0.1;

    (void) fprintf(plot, "s %d %d %d %d\n",0,0,
                   (int) ((margin*2.0+duration)*100.0), 5*nproc);
    /* (void) fprintf(stderr, "nproc=%d, duration=%4.2f\n", nproc, duration); */

    /* Now go thru the files and actually parse the contents */

    for (i=lo; i<=hi; i++) {
        (void) sprintf(filename, "events.%03d", i);

        if ( !(file = fopen(filename, "r")) ) {
            break;
        }

        GetTimeSpan(file, &start, &span);

        comms = 0.0;
        state = 5*(i-lo);
        otim = 0.0 + margin;

        (void) fprintf(plot, "t %d %d %d\n", 0, state, i);

        while ( token = ReadToken(file, ':') ) {
            if (strcmp(token, "BEGIN") == 0)  {
                newstate = 5*(i-lo) + 1;
            }
            else if (strcmp(token, "END") == 0) {
                newstate = 5*(i-lo);
            }
            else {
                continue;
            }

            /* Have a BEGIN or END ... only process Snd/Rcv at moment */

            event = ReadToken(file, ':');
            if ((strcmp(event, "Snd") == 0)
                    || (strcmp(event,"Rcv") == 0)
                    || (strcmp(event,"Waitcom") == 0)) {

                if (FindToken("TIME", file, ':')) {
                    time = atof(ReadToken(file, ':')) - start + margin;

                    (void) fprintf(plot, "l %d %d %d %d\n",
                                   (int) (otim*100.0), state,
                                   (int) (time*100.0), state);
                    (void) fprintf(plot, "l %d %d %d %d\n",
                                   (int) (time*100.0), state,
                                   (int) (time*100.0), newstate);

                    /* Accumulate the time spent in communication */

                    if (newstate == (5*(i-lo)))  {
                        comms = comms + time - otim;
                    }

                    otim = time;
                    state = newstate;
                }
            }
            else if (strcmp(event, "Process") == 0) {
                if (FindToken("TIME", file, ':')) {
                    time = atof(ReadToken(file, ':')) - start + margin;

                    (void) fprintf(plot, "l %d %d %d %d\n", (int) (otim*100.0), state,
                                   (int) (time*100.0), state);
                    otim = time;
                }
            }
        }

        /* Assume that non-communication time is useful */

        useful = 100.0 * (span - comms) / duration;

        (void) fprintf(plot, "t %d %d %4.1f%%\n",
                       (int) (100.0*duration+150.0*margin), 5*(i-lo), useful);

        (void) fflush(plot);
        (void) fclose(file);
    }
    return 0;
}
