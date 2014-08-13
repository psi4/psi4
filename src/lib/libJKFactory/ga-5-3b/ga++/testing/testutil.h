#define add_range_       F77_FUNC_(add_range,ADD_RANGE)
#define compare_patches_ F77_FUNC_(compare_patches,COMPARE_PATCHES)
#define copy_range_      F77_FUNC_(copy_range,COPY_RANGE)
#define dot_range_       F77_FUNC_(dot_range,DOT_RANGE)
#define init_array_      F77_FUNC_(init_array,INIT_ARRAY)
#define print_range_     F77_FUNC_(print_range,PRINT_RANGE)
#define scale_patch_     F77_FUNC_(scale_patch,SCALE_PATCH)
#define util_mdtob_      F77_FUNC_(util_mdtob,UTIL_MDTOB)
#define util_mitob_      F77_FUNC_(util_mitob,UTIL_MITOB)
#define util_timer_      F77_FUNC_(util_timer,UTIL_TIMER)

#ifdef __cplusplus
extern "C" {
#endif

extern void get_range( int ndim, int dims[], int lo[], int hi[]);
extern void new_range(int ndim, int dims[], int lo[], int hi[],
                      int new_lo[], int new_hi[]);
extern void print_range(char *pre,int ndim, int lo[], int hi[], char* post);
extern void print_subscript(char *pre,int ndim, int subscript[], char* post);
extern void print_distribution(int g_a);

#ifdef __cplusplus
}
#endif
