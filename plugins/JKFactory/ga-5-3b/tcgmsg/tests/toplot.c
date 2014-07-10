#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_LIBPLOT

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_CTYPE_H
#   include <ctype.h>
#endif

extern void openpl();
extern void erase();
extern void label();
extern void line();
extern void circle();
extern void arc();
extern void move();
extern void cont();
extern void point();
extern void linemod();
extern void space();
extern void closepl();

void Error(char *string, int integer)
{
    (void) fflush(stdout);
    (void) fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    (void) fprintf(stderr,"%s %d (0x%x)\n",string, integer, integer);
    (void) fprintf(stderr,"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    (void) fflush(stderr);
    (void) exit(1);
}

static int xmin = 0;
static int ymin = 0;

#define SCALE(X, XMIN, XSCALE) \
2048 + (int) (XSCALE * (DoublePrecision) ( (X) - (XMIN) ) )


/**
 * Crude ascii interface into UNIX binary plot file format.
 *
 * All co-ordinates are integer values.
 *
 * s xmin ymin xmax ymax  ... defines user co-ordinates
 *
 * f 1 ... solid lines
 * f 2 ... dotted
 * f 3 ... shortdashed
 * f 4 ... longdashed
 *
 * t x y text_string_with_no_spaces_in_it ... draw text at coords x,y
 *
 * l x1 y1 x2 y2 ... draw line from (x1,y1)->(x2,y2)
 */
int main(int argc, char **argv)
{
    int test, i1,i2,i3,i4,nline=0,x1,x2,y1,y2;
    char string[1024];
    DoublePrecision scalex, scaley;

    openpl();
    space(0, 0, 32767, 32767);    /* open as large as possible */

    while ( (test = getchar()) != EOF ) {
        if ( !isspace(test) ) {
            nline++;
            switch (test) {
                case 's':
                    if (scanf("%d %d %d %d",&i1,&i2,&i3,&i4) != 4)
                        Error("plot: scanning space arguments, line=",nline);
                    else {
                        xmin = i1;
                        ymin = i2;
                        scalex = (32767.0-4096.0) / (DoublePrecision) (i3 - i1);
                        scaley = (32767.0-4096.0) / (DoublePrecision) (i4 - i2);
                    }
                    break;

                case 'l':
                    if (scanf("%d %d %d %d",&i1,&i2,&i3,&i4) != 4)
                        Error("plot: scanning line arguments, line=",nline);
                    x1 = SCALE(i1, xmin, scalex);
                    x2 = SCALE(i3, xmin, scalex);
                    y1 = SCALE(i2, ymin, scaley);
                    y2 = SCALE(i4, ymin, scaley);
                    line(x1, y1, x2, y2);
                    break;

                case 'f':
                    if (scanf("%d",&i1) != 1)
                        Error("plot: scanning linemode arguments",-1);
                    switch (i1) {
                        case 1:
                            linemod("solid");
                            break;
                        case 2:
                            linemod("dotted");
                            break;
                        case 3:
                            linemod("shortdashed");
                            break;
                        case 4:
                            linemod("longdashed");
                            break;
                        default:
                            Error("plot: unknown linemode",i1);
                    }
                    break;

                case 't':
                    if (scanf("%d %d %s",&i1, &i2, string) != 3)
                        Error("plot: scanning text arguments", -1);

                    x1 = SCALE(i1, xmin, scalex);
                    y1 = SCALE(i2, ymin, scaley);
                    move(x1, y1);
                    label(string);
                    break;

                default:
                    Error("plot: unknown directive, line=",nline);
                    break;
            }
        }
    }

    closepl();

    (void) fflush(stdout);

    return 0;
}

#else /* HAVE_LIBPLOT */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

int main(int argc, char **argv)
{
    printf("no plot utility on this system\n");
    return 0;
}

#endif /* HAVE_LIBPLOT */
