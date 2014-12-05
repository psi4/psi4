#ifndef EFPMD_MSG_H
#define EFPMD_MSG_H

#include <stdarg.h>
#include <stdio.h>

void msg(const char *, ...);
void fmsg(FILE *, const char *, ...);
void vmsg(const char *, va_list);
void vfmsg(FILE *, const char *, va_list);

#endif /* EFPMD_MSG_H */
