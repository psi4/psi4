#include "ga.h"
#include "mock.h"
#include "ga_unit.h"

static int test(int shape_idx, int type_idx, int dist_idx)
{
    int type = TYPES[type_idx];
    int *dims = SHAPES[shape_idx];
    int ndim = SHAPES_NDIM[shape_idx];
    mock_ga_t *mock_a, *result_a;
    int g_a;
    int buffer[100];
    int lo[GA_MAX_DIM], hi[GA_MAX_DIM], ld[GA_MAX_DIM], shape[GA_MAX_DIM];
    int result=0, error_index=-1, error_proc=-1;

    /* create the local array and result array */
    mock_a = Mock_Create(type, ndim, dims, "mock", NULL);
    result_a = Mock_Create(type, ndim, dims, "mock", NULL);

    /* create the global array */
    g_a = create_function[dist_idx](type, ndim, dims);

    /* create meaningful data for local array */
    mock_data(mock_a, g_a);

    /* init global array with same data as local array */
    mock_to_global(mock_a, g_a);

    /* call the local routine */
    Mock_Abs_value(mock_a);

    /* call the global routine */
    GA_Abs_value(g_a);

    /* get the results from the global array */
    global_to_mock(g_a, result_a);

    /* compare the results */
    result = neq_mock(mock_a, result_a, &error_index);
    if (0 != result) {
        error_proc = GA_Nodeid();
    }
    /* make sure all procs get same result so they can die gracefully */
    GA_Igop(&result, 1, "+");
    /* if error occured, find the highest failing node ID */
    GA_Igop(&error_proc, 1, "max");
    /* clear the error index for all but the highest failing node ID */
    if (error_proc != GA_Nodeid()) {
        error_index = 0;
    }
    /* make sure all procs get the error index on the highest failing node ID */
    GA_Igop(&error_index, 1, "+");
    if (0 != result) {
        if (error_proc == GA_Nodeid()) {
            printf("ERROR: local result failed to compare to global result\n");
            printf("\terror_proc=%d\n", error_proc);
            printf("\terror_index=%d\n", error_index);
            printf("***LOCAL RESULT***\n");
            Mock_Print(mock_a);
            printf("***GLOBAL RESULT***\n");
            Mock_Print(result_a);
            printf("\tprinting array distribution\n");
        }
        GA_Sync();
        GA_Print(g_a);
        GA_Print_distribution(g_a);
        return 1;
    }

    /* clean up */
    Mock_Destroy(mock_a);
    Mock_Destroy(result_a);
    GA_Destroy(g_a);

    return 0;
}

int main(int argc, char **argv)
{
    TEST_SETUP;

    int shape_idx=0, type_idx=0, dist_idx=0;
    int return_code=0;

    for (shape_idx=0; shape_idx < NUM_SHAPES; ++shape_idx) {
        for (type_idx=0; type_idx < NUM_TYPES; ++type_idx) {
            for (dist_idx=0; dist_idx < NUM_DISTS; ++dist_idx) {
                if (0 == GA_Nodeid()) {
                    printf("%s\t%s\t%s\n",
                            SHAPE_NAMES[shape_idx],
                            TYPE_NAMES[type_idx],
                            DIST_NAMES[dist_idx]
                            );
                }
                GA_Sync();
                return_code = test(shape_idx, type_idx, dist_idx);
                if (0 != return_code) {
                    break;
                }
            }
            if (0 != return_code) {
                break;
            }
        }
        if (0 != return_code) {
            break;
        }
    }

    TEST_TEARDOWN;
    return return_code;
}
