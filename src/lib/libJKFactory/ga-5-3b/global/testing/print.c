#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include <ga.h>
#include <macdecls.h>
#include <mp3.h>

int main(int argc, char **argv)
{
    int size_dst = 15;
    int g_a = 0;
    int I_NEG_ONE = -1;
    long L_NEG_ONE = -1;
    long long LL_NEG_ONE = -1;
    int FIVE = 5;
    int TEN = 10;
    int lo;
    int hi;
    int *ptr;
    int i;

    MP_INIT(argc,argv);

    GA_INIT(argc,argv);

    for (i=0; i<3; ++i) {
        if (0 == i) {
            g_a = NGA_Create(C_INT, 1, &size_dst, "dst", NULL);
            GA_Fill(g_a, &I_NEG_ONE);
        } else if (1 == i) {
            g_a = NGA_Create(C_LONG, 1, &size_dst, "dst", NULL);
            GA_Fill(g_a, &L_NEG_ONE);
        } else if (2 == i) {
            g_a = NGA_Create(C_LONGLONG, 1, &size_dst, "dst", NULL);
            GA_Fill(g_a, &LL_NEG_ONE);
        }
        GA_Sync();
        GA_Print(g_a);
        NGA_Print_patch(g_a, &FIVE, &TEN, 0);
        NGA_Print_patch(g_a, &FIVE, &TEN, 1);
        NGA_Distribution(g_a, GA_Nodeid(), &lo, &hi);
        NGA_Access(g_a, &lo, &hi, &ptr, NULL);
        printf("[%d] (%d)=%d\n", GA_Nodeid(), lo, *ptr);
        NGA_Release(g_a, &lo, &hi);
    }

    GA_Terminate();
    MP_FINALIZE();
    exit(EXIT_SUCCESS);
}
