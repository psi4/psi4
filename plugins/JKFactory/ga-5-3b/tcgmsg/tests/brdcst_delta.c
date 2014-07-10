#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "sndrcv.h"

/* Given the col and row no. return the actual process no. */
#define MAP(Row,Col) (ncols*(Row) + (Col))

/* Given the node return the row no. */
#define ROW(Node) ((Node) / ncols)

/* Given the node return the column no. */
#define COL(Node) ((Node) - ncols*((Node)/ncols))

static Integer ncols;
static Integer nrows;

/**
 * broadcast buffer to all other processes from process originator
 * ... all processes call this routine specifying the same
 * orginating process
 *
 * Modified for grid nrows * ncols.
 *
 * Always send message to process 0. Then proceed from there.
 *
 * Algorithm ... recursively bisect the grid horizontally
 *                                  and then vertically
 *
 * Need to modify to include pipelining.
 *
 * 4 x 6 grid is numbered (nrows=4, ncols=6)
 *
 * 0  1  2  3   4  5
 * 6  7  8  9  10 11
 * 12 13 14 15 16 17
 * 18 19 20 21 22 23
 */
void brdcst_delta_(Integer *type, char *buf, Integer *lenbuf, Integer *originator)
{
    Integer me = NODEID_();
    Integer mycol = COL(me);
    Integer myrow = ROW(me);
    Integer sync = 1;
    Integer from, id, left, middle, right, lenmes;

    /* First try implementation ... always send data to process 0 */

    if (*originator != 0) {
        if (me == 0) {
            (void) printf("a %d receiving from %d %d\n",me,*originator);
            (void) fflush(stdout);
            RCV_(type, buf, lenbuf, &lenmes, originator, &from, &sync);
        }
        else if (me == *originator) {
            id = 0;
            (void) printf("a %d sending to %d type %d\n",me,*originator,*type);
            (void) fflush(stdout);
            SND_(type, buf, lenbuf, &id, &sync);
        }
    }

    /* Now broadcast from process 0 */

    /* Bisect aInteger top horizonal edge of mesh */

    if (myrow == 0) {
        (void) printf("%d myrow == 0\n",NODEID_());
        left = 0;
        right = ncols-1;
        while (left != right) {
            middle = (left + right + 1) / 2;
            if (mycol == left) {
                id = MAP((Integer) 0,middle);
                (void) printf("b %d sending to %d type %d\n",me,id,*type);
                (void) fflush(stdout);
                SND_(type, buf, lenbuf, &id, &sync);
            }
            else if (mycol == middle) {
                id = MAP((Integer) 0,left);
                (void) printf("b %d receiving from %d\n",me,id);
                (void) fflush(stdout);
                RCV_(type, buf, lenbuf, &lenmes, &id, &from, &sync);
            }
            if (mycol < middle) {
                right = middle-1;
            } else {
                left = middle;
            }
        }
    }

    /* Bisect down vertical columns of mesh */

    left = 0;
    right = nrows-1;
    while (left != right) {
        middle = (left + right + 1) / 2;
        if (myrow == left) {
            id = MAP(middle,mycol);
            (void) printf("c %d sending to %d type %d\n",me,id,*type);
            (void) fflush(stdout);
            SND_(type, buf, lenbuf, &id, &sync);
        }
        else if (myrow == middle) {
            id = MAP(left,mycol);
            (void) printf("c %d receiving from %d %d\n",me,id);
            (void) fflush(stdout);
            RCV_(type, buf, lenbuf, &lenmes, &id, &from, &sync);
        }
        if (myrow < middle) {
            right = middle-1;
        } else {
            left = middle;
        }
    }
    (void) fflush(stdout);
}

int main(int argc, char **argv)
{
    Integer row, col, node, data, type, len, me;

    pbegin(argc, argv);
    LLOG_();

    if (NODEID_() == 0) {
        (void) printf("Input nrows, ncols ");
        (void) scanf("%d %d",&nrows, &ncols);
    }

    node = 0;
    type = 1;
    len = 4;

    BRDCST_(&type, &nrows, &len, &node);
    BRDCST_(&type, &ncols, &len, &node);

    me = NODEID_();
    (void) printf(" me=%d row=%d col=%d map=%d\n",me,
                  ROW(me),COL(me),MAP(ROW(me),COL(me)));

    /*  SETDBG_(&type); */

    for (node=0; node<NNODES_(); node++) {
        type = node*10 + 1;
        if (NODEID_() == node) {
            data = node*1000;
        } else {
            data = -1;
        }
        brdcst_delta_(&type, &data, &len, &node);
        (void) fflush(stdout);
        (void) fflush(stderr);
        (void) SYNCH_(&type);
        if (data != (node*1000)) {
            Error("Invalid data", data);
        }
    }

    PEND_();

    return 0;
}
