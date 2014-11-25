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
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "cfg.h"

enum node_type {
	NODE_TYPE_INT,
	NODE_TYPE_DOUBLE,
	NODE_TYPE_BOOL,
	NODE_TYPE_STRING,
	NODE_TYPE_ENUM
};

struct node {
	char *name;
	enum node_type type;
	struct node *next;
	union xdata {
		int xint;
		double xdouble;
		bool xbool;
		char *xstring;
		const struct esv *xenum;
	} xdata;
	struct esv {
		char *str;
		int val;
		struct esv *next;
	} *esv;
};

struct cfg {
	char error[256];
	struct node *nodes;
};

static char *xstrdup(const char *str)
{
	char *dup = malloc(strlen(str) + 1);
	return strcpy(dup, str);
}

static char *xstrndup(const char *str, size_t len)
{
	char *dup;

	dup = malloc(len + 1);
	dup[len] = '\0';

	return strncpy(dup, str, len);
}

static const struct esv *esv_from_val(const struct esv *esv, int val)
{
	while (esv) {
		if (esv->val == val)
			return esv;

		esv = esv->next;
	}

	return NULL;
}

static bool xstrtobool(const char *str, const char **endptr)
{
	if (strncmp(str, "true", strlen("true")) == 0) {
		if (endptr)
			*endptr = str + strlen("true");

		return true;
	}

	if (strncmp(str, "false", strlen("false")) == 0) {
		if (endptr)
			*endptr = str + strlen("false");

		return false;
	}

	if (endptr)
		*endptr = str;

	return false;
}

static char *xstrtostr(const char *str, const char **endptr)
{
	size_t cnt = 0;
	const char *ptr;

	if (endptr)
		*endptr = str;

	if (!*str)
		return NULL;

	if (*str == '"') {
		ptr = ++str;

		while (*ptr && *ptr != '"')
			ptr++, cnt++;

		if (!*ptr)
			return NULL;

		ptr++;
	}
	else {
		cnt = strlen(str);

		while (isspace(str[cnt - 1]))
			cnt--;

		ptr = str + cnt;
	}

	if (endptr)
		*endptr = ptr;

	return xstrndup(str, cnt);
}

static const struct esv *xstrtoenum(const char *str, const char **endptr,
					const struct esv *esv)
{
	size_t len = 0;

	while (esv) {
		len = strlen(esv->str);

		if (strncmp(str, esv->str, len) == 0)
			if (str[len] == '\0' || isspace(str[len]))
				break;

		len = 0;
		esv = esv->next;
	}

	if (endptr)
		*endptr = str + len;

	return esv;
}

static bool error(struct cfg *cfg, const char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vsnprintf(cfg->error, sizeof(cfg->error), fmt, args);
	va_end(args);

	return false;
}

static struct node *node_find(struct node *nodes, const char *name)
{
	while (nodes) {
		if (strcmp(nodes->name, name) == 0)
			return nodes;

		nodes = nodes->next;
	}

	return NULL;
}

static void print_node(const struct node *node, FILE *out)
{
	if (!node)
		return;

	/* use recursion to print in correct (reverse) order */
	print_node(node->next, out);

	fprintf(out, "%s ", node->name);

	switch (node->type) {
	case NODE_TYPE_INT:
		fprintf(out, "%d\n", node->xdata.xint);
		break;
	case NODE_TYPE_DOUBLE:
		fprintf(out, "%g\n", node->xdata.xdouble);
		break;
	case NODE_TYPE_BOOL:
		fprintf(out, "%s\n", node->xdata.xbool ? "true" : "false");
		break;
	case NODE_TYPE_STRING:
		fprintf(out, "%s\n", node->xdata.xstring);
		break;
	case NODE_TYPE_ENUM:
		fprintf(out, "%s\n", node->xdata.xenum->str);
		break;
	}
}

static void add_node(struct cfg *cfg, const char *name, enum node_type type)
{
	struct node *node;

	assert(node_find(cfg->nodes, name) == NULL);

	node = calloc(1, sizeof(struct node));
	node->name = xstrdup(name);
	node->type = type;
	node->next = cfg->nodes;
	cfg->nodes = node;
}

struct cfg *cfg_create(void)
{
	return calloc(1, sizeof(struct cfg));
}

void cfg_add_int(struct cfg *cfg, const char *name, int def)
{
	assert(cfg);
	assert(name);

	add_node(cfg, name, NODE_TYPE_INT);
	cfg->nodes->xdata.xint = def;
}

void cfg_add_double(struct cfg *cfg, const char *name, double def)
{
	assert(cfg);
	assert(name);

	add_node(cfg, name, NODE_TYPE_DOUBLE);
	cfg->nodes->xdata.xdouble = def;
}

void cfg_add_bool(struct cfg *cfg, const char *name, bool def)
{
	assert(cfg);
	assert(name);

	add_node(cfg, name, NODE_TYPE_BOOL);
	cfg->nodes->xdata.xbool = def;
}

void cfg_add_string(struct cfg *cfg, const char *name, const char *def)
{
	assert(cfg);
	assert(name);
	assert(def);

	add_node(cfg, name, NODE_TYPE_STRING);
	cfg->nodes->xdata.xstring = xstrdup(def);
}

void cfg_add_enum(struct cfg *cfg, const char *name, int def,
			const char *keys, const int *vals)
{
	assert(cfg);
	assert(name);
	assert(keys);
	assert(vals);

	add_node(cfg, name, NODE_TYPE_ENUM);
	struct node *node = cfg->nodes;

	while (*keys) {
		struct esv *esv = malloc(sizeof(struct esv));
		const char *nxt = strchr(keys, '\n');
		assert(nxt);

		esv->str = xstrndup(keys, (size_t)(nxt - keys));
		esv->val = *vals++;
		esv->next = node->esv;
		node->esv = esv;
		keys = nxt + 1;
	}

	node->xdata.xenum = esv_from_val(node->esv, def);
	assert(node->xdata.xenum);
}

void cfg_set_int(struct cfg *cfg, const char *name, int val)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_INT);

	node->xdata.xint = val;
}

void cfg_set_double(struct cfg *cfg, const char *name, double val)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_DOUBLE);

	node->xdata.xdouble = val;
}

void cfg_set_bool(struct cfg *cfg, const char *name, bool val)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_BOOL);

	node->xdata.xbool = val;
}

void cfg_set_string(struct cfg *cfg, const char *name, const char *val)
{
	assert(cfg);
	assert(name);
	assert(val);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_STRING);

	free(node->xdata.xstring);
	node->xdata.xstring = xstrdup(val);
}

void cfg_set_enum(struct cfg *cfg, const char *name, int val)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_ENUM);

	node->xdata.xenum = esv_from_val(node->esv, val);
	assert(node->xdata.xenum);
}

int cfg_get_int(const struct cfg *cfg, const char *name)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_INT);

	return node->xdata.xint;
}

double cfg_get_double(const struct cfg *cfg, const char *name)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_DOUBLE);

	return node->xdata.xdouble;
}

bool cfg_get_bool(const struct cfg *cfg, const char *name)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_BOOL);

	return node->xdata.xbool;
}

const char *cfg_get_string(const struct cfg *cfg, const char *name)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_STRING);

	return node->xdata.xstring;
}

int cfg_get_enum(const struct cfg *cfg, const char *name)
{
	assert(cfg);
	assert(name);

	struct node *node = node_find(cfg->nodes, name);
	assert(node);
	assert(node->type == NODE_TYPE_ENUM);

	return node->xdata.xenum->val;
}

bool cfg_parse_line(struct cfg *cfg, const char *line)
{
	assert(cfg);
	assert(line);

	bool res = true;
	union xdata xdata;
	const char *endptr;
	struct node *node = cfg->nodes;

	while (*line && isspace(*line))
		line++;

	if (!*line)
		return true;

	while (node) {
		size_t len = strlen(node->name);

		if (strncmp(node->name, line, len) == 0)
			if (line[len] == '\0' || isspace(line[len]))
				break;

		node = node->next;
	}

	if (!node) {
		endptr = line;

		while (*endptr && !isspace(*endptr))
			endptr++;

		return error(cfg, "%.*s: unknown option", endptr - line, line);
	}

	line += strlen(node->name);

	while (*line && isspace(*line))
		line++;

	switch (node->type) {
	case NODE_TYPE_INT:
		xdata.xint = (int)strtol(line, (char **)&endptr, 10);
		break;
	case NODE_TYPE_DOUBLE:
		xdata.xdouble = strtod(line, (char **)&endptr);
		break;
	case NODE_TYPE_BOOL:
		xdata.xbool = xstrtobool(line, &endptr);
		break;
	case NODE_TYPE_STRING:
		xdata.xstring = xstrtostr(line, &endptr);
		break;
	case NODE_TYPE_ENUM:
		xdata.xenum = xstrtoenum(line, &endptr, node->esv);
		break;
	}

	if (endptr == line) {
		res = error(cfg, "%s: incorrect format", node->name);
		goto fail;
	}

	while (*endptr) {
		if (!isspace(*endptr)) {
			res = error(cfg, "%s: unexpected end of line", node->name);
			goto fail;
		}
		endptr++;
	}

fail:
	if (node->type == NODE_TYPE_STRING) {
		if (res)
			free(node->xdata.xstring);
		else
			free(xdata.xstring);
	}

	if (res)
		node->xdata = xdata;

	return res;
}

void cfg_print(const struct cfg *cfg, FILE *out)
{
	assert(cfg);
	assert(out);

	print_node(cfg->nodes, out);
}

const char *cfg_get_last_error(const struct cfg *cfg)
{
	assert(cfg);
	return cfg->error;
}

void cfg_free(struct cfg *cfg)
{
	if (cfg) {
		while (cfg->nodes) {
			struct node *node = cfg->nodes;
			cfg->nodes = node->next;

			if (node->type == NODE_TYPE_STRING) {
				free(node->xdata.xstring);
			}
			else if (node->type == NODE_TYPE_ENUM) {
				while (node->esv) {
					struct esv *esv = node->esv;
					node->esv = esv->next;

					free(esv->str);
					free(esv);
				}
			}

			free(node->name);
			free(node);
		}

		free(cfg);
	}
}
