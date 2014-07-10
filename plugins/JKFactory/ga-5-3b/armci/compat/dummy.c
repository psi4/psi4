#if HAVE_CONFIG_H
#   include "config.h"
#endif

int dummy_func_for_nonempty_libcompat_armci(int dummy)
{
    int dont_optimize_away = 4;
    return dont_optimize_away + dummy;
}
