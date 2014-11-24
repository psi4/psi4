/*-
 * Copyright (c) 2012-2014 Ilya Kaliman
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

#ifndef LIBEFP_STREAM_H
#define LIBEFP_STREAM_H

#include <stddef.h>

struct stream;

struct stream *efp_stream_open(const char *);
void efp_stream_set_split_char(struct stream *, char);
char efp_stream_get_split_char(struct stream *);
const char *efp_stream_get_ptr(struct stream *);
void efp_stream_next_line(struct stream *);
void efp_stream_reset_line(struct stream *);
char efp_stream_get_char(struct stream *);
char efp_stream_current_char(struct stream *);
int efp_stream_parse_int(struct stream *, int *);
int efp_stream_parse_double(struct stream *, double *);
int efp_stream_advance(struct stream *, size_t);
void efp_stream_skip_space(struct stream *);
void efp_stream_skip_nonspace(struct stream *);
int efp_stream_eol(struct stream *);
int efp_stream_eof(struct stream *);
void efp_stream_close(struct stream *);

#endif /* LIBEFP_STREAM_H */
