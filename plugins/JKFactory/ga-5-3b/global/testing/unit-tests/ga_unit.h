#ifndef GA_UNIT_H
#define GA_UNIT_H

#include "mp3.h"

#define GA_PRINT_MSG() printf("Test Completed\n")

#define GA_COMPLETE_MSG() printf("Test Completed\n")

#define GA_ERROR_MSG() printf("GA ERROR\n")

#define GA_ERROR_MSG2() printf("GA ERROR\n")

/* the C long long type is rarely implemented for most GA operations */
#define TEST_LONGLONG 0

#if TEST_LONGLONG
#   define NUM_TYPES 7
#else
#   define NUM_TYPES 6
#endif
int TYPES[NUM_TYPES] = {
    C_INT,
    C_LONG,
#if TEST_LONGLONG
    C_LONGLONG,
#endif
    C_FLOAT,
    C_DBL,
    C_SCPL,
    C_DCPL,
};

char* TYPE_NAMES[NUM_TYPES] = {
    "C_INT",
    "C_LONG",
#if TEST_LONGLONG
    "C_LONGLONG",
#endif
    "C_FLOAT",
    "C_DBL",
    "C_SCPL",
    "C_DCPL",
};

enum dist_type {
    DIST_REGULAR=0,
    DIST_CHUNK,
    DIST_IRREGULAR,
    DIST_BLOCK_CYCLIC,
    DIST_SCALAPACK,
    NUM_DISTS,
};

int DIST_TYPES[NUM_DISTS] = {
    DIST_REGULAR,
    DIST_CHUNK,
    DIST_IRREGULAR,
    DIST_BLOCK_CYCLIC,
    DIST_SCALAPACK,
};

char* DIST_NAMES[NUM_DISTS] = {
    "DIST_REGULAR",
    "DIST_CHUNK",
    "DIST_IRREGULAR",
    "DIST_BLOCK_CYCLIC",
    "DIST_SCALAPACK",
};

#define NUM_SHAPES 4
static int SHAPES_ZERO[] = {10};
#define    SHAPES_ZERO_NDIM 1
#define    SHAPES_ZERO_NAME "10"
static int SHAPES_ONE[] = {2,3};
#define    SHAPES_ONE_NDIM 2
#define    SHAPES_ONE_NAME "2x3"
static int SHAPES_TWO[] = {2,3,4};
#define    SHAPES_TWO_NDIM 3
#define    SHAPES_TWO_NAME "2x3x4"
static int SHAPES_THREE[] = {2,3,4,5};
#define    SHAPES_THREE_NDIM 4
#define    SHAPES_THREE_NAME "2x3x4x5"

static int* SHAPES[] = {
    SHAPES_ZERO,
    SHAPES_ONE,
    SHAPES_TWO,
    SHAPES_THREE,
    NULL
};
static char* SHAPE_NAMES[] = {
    SHAPES_ZERO_NAME,
    SHAPES_ONE_NAME,
    SHAPES_TWO_NAME,
    SHAPES_THREE_NAME,
    NULL
};
static int SHAPES_NDIM[] = {
    SHAPES_ZERO_NDIM,
    SHAPES_ONE_NDIM,
    SHAPES_TWO_NDIM,
    SHAPES_THREE_NDIM,
};

//#define OP_TYPES 6
//char operators[OP_TYPES] = {'+', '*', 'max', 'min', 'absmax', 'absmin'};

#define TEST_SETUP    GA_Initialize_args(&argc, &argv)
#define TEST_TEARDOWN GA_Terminate(); MP_FINALIZE()

static void aprint(char *name, int *array, int size)
{
    int i;
    printf("%s={", name);
    if (size > 0) {
        printf("%d", array[0]);
    }
    for (i=1; i<size; ++i) {
        printf(",%d", array[i]);
    }
    printf("}\n");
}

static int create_regular(int type, int ndim, int *shape) {
    return NGA_Create(type, ndim, shape, "name", NULL);
}

/* chunk is based on a random value between 0 and the dimension's limit */
static int create_chunk(int type, int ndim, int *shape) {
#if 1
    int i;
    int chunk[GA_MAX_DIM];

    /* the chunks must be the same on each process! */
    if (0 == GA_Nodeid()) {
        for (i=0; i<ndim; ++i) {
            chunk[i] = rand() % shape[i];
        }
    }
    GA_Brdcst(chunk, ndim*sizeof(int), 0);
    if (0 == GA_Nodeid()) {
        printf("\tcreating chunked array\n");
        printf("\t");
        aprint("shape", shape, ndim);
        printf("\t");
        aprint("chunk", chunk, ndim);
    }
    return NGA_Create(type, ndim, shape, "name", chunk);
#else
    return NGA_Create(type, ndim, shape, "name", NULL);
#endif
}

static int create_irregular(int type, int ndim, int *shape) {
    return NGA_Create(type, ndim, shape, "name", NULL);
}

static int create_block_cyclic(int type, int ndim, int *shape) {
    return NGA_Create(type, ndim, shape, "name", NULL);
}

static int create_scalapack(int type, int ndim, int *shape) {
    return NGA_Create(type, ndim, shape, "name", NULL);
}

typedef int(*creator)(int,int,int*);

static creator create_function[] = {
    create_regular,
    create_chunk,
    create_irregular,
    create_block_cyclic,
    create_scalapack,
};

#endif
