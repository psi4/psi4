/*-
 * Copyright (c) 2012-2013 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#ifndef EFPMD_CFG_H
#define EFPMD_CFG_H

#include <stdbool.h>
#include <stdio.h>

struct cfg;

struct cfg *cfg_create(void);
void cfg_add_int(struct cfg *, const char *, int);
void cfg_add_double(struct cfg *, const char *, double);
void cfg_add_bool(struct cfg *, const char *, bool);
void cfg_add_string(struct cfg *, const char *, const char *);
void cfg_add_enum(struct cfg *, const char *, int, const char *, const int *);
void cfg_set_int(struct cfg *, const char *, int);
void cfg_set_double(struct cfg *, const char *, double);
void cfg_set_bool(struct cfg *, const char *, bool);
void cfg_set_string(struct cfg *, const char *, const char *);
void cfg_set_enum(struct cfg *, const char *, int);
int cfg_get_int(const struct cfg *, const char *);
double cfg_get_double(const struct cfg *, const char *);
bool cfg_get_bool(const struct cfg *, const char *);
const char *cfg_get_string(const struct cfg *, const char *);
int cfg_get_enum(const struct cfg *, const char *);
bool cfg_parse_line(struct cfg *, const char *);
void cfg_print(const struct cfg *, FILE *);
const char *cfg_get_last_error(const struct cfg *);
void cfg_free(struct cfg *);

#endif /* EFPMD_CFG_H */
