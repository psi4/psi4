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
#include <stdint.h>
#include <stdlib.h>

#include "bvec.h"

struct bvec {
	size_t len;
	uint64_t *data;
};

#define EBITS 64

struct bvec *efp_bvec_create(size_t len)
{
	struct bvec *bvec = (struct bvec *)malloc(sizeof(struct bvec));

	if (!bvec)
		return NULL;

	size_t cnt = len / EBITS + (len % EBITS == 0 ? 0 : 1);

	bvec->len = len;
	bvec->data = (uint64_t *)calloc(cnt, sizeof(uint64_t));

	if (!bvec->data) {
		free(bvec);
		return NULL;
	}

	return bvec;
}

void efp_bvec_set(struct bvec *bvec, size_t n)
{
	assert(bvec);
	assert(n < bvec->len);

	int shift = n % EBITS;

	bvec->data[n / EBITS] |= (1ULL << shift);
}

void efp_bvec_unset(struct bvec *bvec, size_t n)
{
	assert(bvec);
	assert(n < bvec->len);

	int shift = n % EBITS;

	bvec->data[n / EBITS] &= ~(1ULL << shift);
}

void efp_bvec_set_value(struct bvec *bvec, size_t n, bool value)
{
	value ? efp_bvec_set(bvec, n) : efp_bvec_unset(bvec, n);
}

bool efp_bvec_is_set(struct bvec *bvec, size_t n)
{
	assert(bvec);
	assert(n < bvec->len);

	int shift = n % EBITS;

	return bvec->data[n / EBITS] & (1ULL << shift);
}

void efp_bvec_free(struct bvec *bvec)
{
	if (bvec) {
		free(bvec->data);
		free(bvec);
	}
}
