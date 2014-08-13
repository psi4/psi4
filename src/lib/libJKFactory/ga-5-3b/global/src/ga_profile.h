/* $Id: ga_profile.h,v 1.3 2005-07-21 08:14:30 manoj Exp $ */

#define GA_PROFILE_PUT 1
#define GA_PROFILE_GET 2
#define GA_PROFILE_ACC 3

extern void ga_profile_init();
extern void ga_profile_terminate();
extern void ga_profile_start(int g_a, long bytes, int ndim, Integer *lo, 
			     Integer *hi, int comm_type);
extern void ga_profile_stop();
