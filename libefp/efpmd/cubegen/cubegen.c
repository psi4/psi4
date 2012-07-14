/* Generate a cube of fragments */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PI 3.14159265358979323846

static int random_int(int max)
{
	return rand() % max;
}

static double random_real(double max)
{
	return (double)rand() / RAND_MAX * max;
}

static void print_frag(const char *name, double space, int a, int b, int c)
{
	double xyzabc[6];

	xyzabc[0] = space * a;
	xyzabc[1] = space * b;
	xyzabc[2] = space * c;
	xyzabc[3] = random_real(2.0 * PI);
	xyzabc[4] = random_real(PI);
	xyzabc[5] = random_real(2.0 * PI);

	printf("fragment %s\n", name);

	for (int i = 0; i < 6; i++)
		printf(" %8.3lf", xyzabc[i]);

	printf("\n");
}

static int select_frag(int n_frag, const int *ratio)
{
	int num = random_int(ratio[n_frag - 1]);

	for (int i = 0; i < n_frag; i++)
		if (ratio[i] > num)
			return i;

	assert(0);
}

static int parse_frag(char *name_arg, char *ratio_arg,
		      int *n_frag_out, char ***name_out, int **ratio_out)
{
	int n_frag = 1;

	for (char *ptr = name_arg; *ptr; ptr++)
		if (*ptr == ':') {
			*ptr = '\0';
			n_frag++;
		}

	char **name = malloc(n_frag * sizeof(char *));
	int *ratio = malloc(n_frag * sizeof(int));

	for (int i = 0; i < n_frag; i++) {
		name[i] = name_arg;
		name_arg += strlen(name_arg) + 1;
		ratio[i] = strtol(ratio_arg, &ratio_arg, 10);

		if (*ratio_arg != ':' && i < n_frag - 1) {
			free(name);
			free(ratio);
			return 1;
		}
		if (i > 0)
			ratio[i] += ratio[i - 1];

		ratio_arg++;
	}

	*n_frag_out = n_frag;
	*name_out = name;
	*ratio_out = ratio;
	return 0;
}

static void print_usage(void)
{
	puts("usage: cubegen <n1:n2:...> <r1:r2:...> <space> <nx> <ny> <nz>");
	puts("");
	puts("  <n1:n2:...>  a colon separated list of fragment names");
	puts("  <r1:r2:...>  a colon separated list of fragment ratios");
	puts("      <space>  distance between the molecules");
	puts("         <nx>  number of molecules in x direction");
	puts("         <ny>  number of molecules in y direction");
	puts("         <nz>  number of molecules in z direction");
	puts("");
	puts("example: cubegen c2h5oh_l:h2o_l 40:60 3.0 20 20 20 > vodka.in");
}

int main(int argc, char **argv)
{
	if (argc < 7) {
		print_usage();
		return EXIT_FAILURE;
	}

	srand(time(NULL));

	char **name;
	int n_frag, *ratio;

	if (parse_frag(argv[1], argv[2], &n_frag, &name, &ratio)) {
		print_usage();
		return EXIT_FAILURE;
	}

	double space = strtod(argv[3], NULL);
	int na = atoi(argv[4]);
	int nb = atoi(argv[5]);
	int nc = atoi(argv[6]);

	for (int a = 0; a < na; a++)
		for (int b = 0; b < nb; b++)
			for (int c = 0; c < nc; c++) {
				int idx = select_frag(n_frag, ratio);
				print_frag(name[idx], space, a, b, c);
			}

	free(name);
	free(ratio);
	return EXIT_SUCCESS;
}
