extern struct
{
    __int64 myt0[10];
    __int64 myt1[10];
    double mypass[10];
} mydata_;

#if 1
inline void mytimer0_ (int *n)
{
    int idx;
    int i;
    idx = *n;
    __int64 t0;
    t0 = __rdtsc ();
    mydata_.myt0[idx] = t0;
}

inline void mytimer1_ (int *n)
{
    int idx;
    idx = *n;
    __int64 t1;
    t1 = __rdtsc ();
    mydata_.myt1[idx] = t1;
    mydata_.mypass[idx] = 
        mydata_.mypass[idx] +
        mydata_.myt1[idx] -
        mydata_.myt0[idx];
}

inline void getmytimer (int idx)
{
    printf ("timer %d takes %le secs\n", idx, mydata_.mypass[idx]/(3.33*1e9));
}
#else
inline void mytimer0_ (int *n)
{
}

inline void mytimer1_ (int *n)
{
}

inline void getmytimer (int idx)
{
}
#endif
