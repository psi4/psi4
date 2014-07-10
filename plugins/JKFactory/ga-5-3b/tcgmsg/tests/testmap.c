#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "typesf2c.h"

int ncols;
int nrows;

/* Given the col and row no. return the actual process no. */
#define MAP(Row,Col) (ncols*(Row) + (Col))

/* Given the node return the row no. */
#define ROW(Node) ((Node) / ncols)

/* Given the node return the column no. */
#define COL(Node) ((Node) - ncols*((Node)/ncols))

int main(int argc, char **argv)
{
    int node, type, len, me;

    (void) printf("Input nrows, ncols ");
    (void) scanf("%d %d",&nrows, &ncols);

    node = 0;
    type = 1;
    len = 4;

    for (me=0; me<(nrows*ncols); me++) {
        (void) printf(" me=%d row=%d col=%d map=%d\n",me,
                      ROW(me),COL(me),MAP(ROW(me),COL(me)));
    }
    return 0;
}
