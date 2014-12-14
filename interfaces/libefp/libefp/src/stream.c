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

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

#include "stream.h"

struct stream {
	char *buffer;
	char *ptr;
	FILE *in;
	char split;
};

static char *read_line(FILE *in, char split)
{
	size_t size = 128;
	size_t i = 0;
	char *buffer = (char *)malloc(size);

	if (!buffer)
		return NULL;

	for(;;) {
		int ch = getc(in);

		if (split != '\0' && ch == split) {
			ch = getc(in);

			if (ch == '\n') {
				continue;
			}
			else {
				ungetc(ch, in);
				ch = split;
			}
		}

		switch(ch) {
		case EOF:
			if (i == 0) {
				free(buffer);
				return NULL;
			}

			/* FALL THROUGH */

		case '\n':
			if (i == size)
				buffer = (char *)realloc(buffer, size + 1);

			buffer[i] = '\0';
			return buffer;

		default:
			buffer[i++] = (char)ch;

			if (i == size) {
				size *= 2;
				buffer = (char *)realloc(buffer, size);
			}
		}
	}

	assert(0);
}

struct stream *efp_stream_open(const char *path)
{
	struct stream *stream;

	assert(path);

	if ((stream = (struct stream *)calloc(1, sizeof(struct stream))) == NULL)
		return NULL;

	stream->in = fopen(path, "r");

	if (!stream->in) {
		free(stream);
		return NULL;
	}

	return stream;
}

void efp_stream_set_split_char(struct stream *stream, char c)
{
	assert(stream);

	stream->split = c;
}

char efp_stream_get_split_char(struct stream *stream)
{
	assert(stream);

	return stream->split;
}

const char *efp_stream_get_ptr(struct stream *stream)
{
	assert(stream);

	return stream->ptr;
}

void efp_stream_next_line(struct stream *stream)
{
	assert(stream);

	if (stream->buffer)
		free(stream->buffer);

	stream->buffer = read_line(stream->in, stream->split);
	stream->ptr = stream->buffer;
}

void efp_stream_reset_line(struct stream *stream)
{
	assert(stream);

	stream->ptr = stream->buffer;
}

char efp_stream_get_char(struct stream *stream)
{
	assert(stream);

	return stream->ptr && *stream->ptr ? *stream->ptr++ : '\0';
}

char efp_stream_current_char(struct stream *stream)
{
	assert(stream);

	return stream->ptr ? *stream->ptr : '\0';
}

int efp_stream_parse_int(struct stream *stream, int *out)
{
	assert(stream);

	char *endptr;
	int x = (int)strtol(stream->ptr, &endptr, 10);

	if (endptr == stream->ptr)
		return 0;

	if (out)
		*out = x;

	stream->ptr = endptr;
	return 1;
}

int efp_stream_parse_double(struct stream *stream, double *out)
{
	assert(stream);

	char *endptr;
	double x = strtod(stream->ptr, &endptr);

	if (endptr == stream->ptr)
		return 0;

	if (out)
		*out = x;

	stream->ptr = endptr;
	return 1;
}

int efp_stream_advance(struct stream *stream, size_t cnt)
{
	assert(stream);

	while (cnt--)
		if (efp_stream_get_char(stream) == '\0')
			return 0;

	return 1;
}

void efp_stream_skip_space(struct stream *stream)
{
	assert(stream);

	if (!stream->ptr)
		return;

	while (*stream->ptr && isspace(*stream->ptr))
		stream->ptr++;
}

void efp_stream_skip_nonspace(struct stream *stream)
{
	assert(stream);

	if (!stream->ptr)
		return;

	while (*stream->ptr && !isspace(*stream->ptr))
		stream->ptr++;
}

int efp_stream_eol(struct stream *stream)
{
	assert(stream);

	return stream->ptr == NULL || *stream->ptr == '\0';
}

int efp_stream_eof(struct stream *stream)
{
	assert(stream);

	return feof(stream->in);
}

void efp_stream_close(struct stream *stream)
{
	if (!stream)
		return;

	if (stream->buffer)
		free(stream->buffer);

	if (stream->in)
		fclose(stream->in);

	free(stream);
}
