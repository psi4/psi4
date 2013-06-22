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

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define PI 3.14159265358979323846

struct pt {
	double x, y, z;
	struct pt *next;
};

static void fail_bad_format(void)
{
	fputs("bad file format", stderr);
	exit(128);
}

static double get_dist(const struct pt *p1, const struct pt *p2,
			double xbox, double ybox, double zbox)
{
	double dx = p2->x - p1->x;
	double dy = p2->y - p1->y;
	double dz = p2->z - p1->z;

	if (xbox > 1.0e-8 && ybox > 1.0e-8 && zbox > 1.0e-8) {
		dx -= xbox * round(dx / xbox);
		dy -= ybox * round(dy / ybox);
		dz -= zbox * round(dz / zbox);
	}

	return sqrt(dx * dx + dy * dy + dz * dz);
}

static void usage(const char *progname)
{
	fprintf(stderr, "Compute Fragment Center of Mass Radial Distribution Function\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n");
	fprintf(stderr, "\t%s [-r <max_dist>] [-n <n_bins>]\n", progname);
	fprintf(stderr, "\t%s [-h]\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Defaults:\n");
	fprintf(stderr, "\tmax_dist 10.0\n");
	fprintf(stderr, "\t  n_bins 100\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Input is read from stdin\n");
	fprintf(stderr, "Output is written to stdout\n");
}

int main(int argc, char **argv)
{
	const char *progname = argv[0];
	int ch, n_bins = 100;
	double max_dist = 10.0;

	while ((ch = getopt(argc, argv, "r:n:h")) != -1) {
		switch (ch) {
		case 'n':
			n_bins = atoi(optarg);
			break;
		case 'r':
			max_dist = strtod(optarg, NULL);
			break;
		default:
			usage(progname);
			exit(127);
		}
	}

	char buffer[512];
	struct pt *list = NULL;
	size_t *hist = calloc(n_bins, sizeof(size_t));
	double xbox, ybox, zbox;

	while (fgets(buffer, sizeof(buffer), stdin)) {
		if (strstr(buffer, "RESTART")) {
			struct pt *pt = list;

			if (fgets(buffer, sizeof(buffer), stdin) == NULL)
				fail_bad_format();

			while (fgets(buffer, sizeof(buffer), stdin)) {
				if (strstr(buffer, "fragment") == NULL)
					break;

				struct pt *next;

				if (!pt) {
					pt = malloc(sizeof(struct pt));
					pt->next = list;
					list = pt;
					next = NULL;
				}
				else {
					next = pt->next;
				}

				if (fgets(buffer, sizeof(buffer), stdin) == NULL)
					fail_bad_format();

				sscanf(buffer, "%lf %lf %lf", &pt->x, &pt->y, &pt->z);

				for (int i = 0; i < 3; i++)
					if (fgets(buffer, sizeof(buffer), stdin) == NULL)
						fail_bad_format();

				pt = next;
			}

			if (pt)
				fail_bad_format();

			while (fgets(buffer, sizeof(buffer), stdin)) {
				const char *ptr = strstr(buffer, "PERIODIC BOX SIZE");

				if (ptr) {
					ptr += strlen("PERIODIC BOX SIZE");
					sscanf(ptr, "%lf %lf %lf", &xbox, &ybox, &zbox);
					break;
				}
				else if (strstr(buffer, "STATE AFTER")) {
					xbox = ybox = zbox = 0.0;
					break;
				}
			}

			for (struct pt *p1 = list; p1; p1 = p1->next) {
				for (struct pt *p2 = p1->next; p2; p2 = p2->next) {
					double dist = get_dist(p1, p2, xbox, ybox, zbox);

					if (dist < max_dist)
						hist[(int)(dist / max_dist * n_bins)] += 2;
				}
			}
		}
	}

	for (int i = 1; i < n_bins; i++) {
		double r = max_dist * i / n_bins;
		double volume = 4.0 * PI * r * r * (max_dist / n_bins);
		double rdf = (double)hist[i] / volume;

		printf("%10.3lf %10.3lf\n", r, rdf);
	}

	while (list) {
		struct pt *pt = list;

		list = list->next;
		free(pt);
	}

	free(hist);
	return 0;
}
