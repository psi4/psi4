#ifndef _macommon_h
#define _macommon_h

#define MA_FALSE 0
#define MA_TRUE  1

#define MA_DEFAULT_SPACE (-1)

#define MA_NAMESIZE    32

#define MT_BASE        1000

#define MT_C_CHAR     (MT_BASE + 0)
#define MT_C_INT      (MT_BASE + 1)
#define MT_C_LONGINT  (MT_BASE + 2)
#define MT_C_FLOAT    (MT_BASE + 3)
#define MT_C_DBL      (MT_BASE + 4)
#define MT_C_LDBL     (MT_BASE + 5)
#define MT_C_SCPL     (MT_BASE + 6)
#define MT_C_DCPL     (MT_BASE + 7)
#define MT_C_LDCPL    (MT_BASE + 8)
                     
#define MT_F_BYTE     (MT_BASE + 9)
#define MT_F_INT      (MT_BASE + 10)
#define MT_F_LOG      (MT_BASE + 11)
#define MT_F_REAL     (MT_BASE + 12)
#define MT_F_DBL      (MT_BASE + 13)
#define MT_F_SCPL     (MT_BASE + 14)
#define MT_F_DCPL     (MT_BASE + 15)

#define MT_C_LONGLONG (MT_BASE + 16)

#define MT_FIRST      MT_C_CHAR
#define MT_LAST       MT_C_LONGLONG
#define MT_NUMTYPES   (MT_LAST - MT_FIRST + 1)

#endif /* _macommon_h */
