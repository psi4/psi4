/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <libciomr/libciomr.h>
#include "dpd.h"
#include "dpd.gbl"

namespace psi {

typedef struct {
    int *left_orbspi;
    int *right_orbspi;
    int *left_orbsym;
    int *right_orbsym;
    int *left_orboff;
    int *right_orboff;
    int permlr;
    int ler;
} dpdpair;

int dpd_init(int dpd_num, int nirreps, long int memory, int cachetype,
             int *cachefiles, int **cachelist, 
             struct dpd_file4_cache_entry *priority, int num_subspaces, ...)
{
  int h,h0,h1,cnt,***dp,l_irrep,r_irrep,p,q;
  int i,j,k,l,*count,offset1,offset2;
  int num_pairs, **orbspi, **orbsym, *numorbs, **orboff;
  int ***pairidx, **pairtot, ****pairorb;
  int **orbidx2, ***orbs2;
  int *tmparray;
  dpdpair *pairs;
  va_list ap;
  dpd_data *this_dpd;
  int nump, nrows, Gp, offset;

  this_dpd = &(dpd_list[dpd_num]);

  this_dpd->nirreps = nirreps;
  this_dpd->num_subspaces = num_subspaces;
  
  dpd_main.memory = memory/sizeof(double);  /* Available memory in doubles */
  dpd_main.memused = 0; /* At first... */
  dpd_main.memcache = 0; /* At first... */
  dpd_main.memlocked = 0; /* At first... */

  dpd_main.cachetype = cachetype;
  dpd_main.cachelist = cachelist;
  dpd_main.cachefiles = cachefiles;
  dpd_main.file4_cache_priority = priority;

  /* Construct binary direct product array */
  dp = (int ***) malloc(nirreps * sizeof(int **));
  for(h=0; h < nirreps; h++) {
    dp[h] = init_int_matrix(nirreps,2);
    cnt=0;
    for(h0=0; h0 < nirreps; h0++) {
      for(h1=0; h1 < nirreps; h1++) {
	if((h0^h1)==h) {
	  dp[h][cnt][0] = h0;
	  dp[h][cnt++][1] = h1;
	}
      }
    }
  }

  /* Grab the irrep population and orbital symmetry arrays from the arg list */
  va_start(ap, num_subspaces);
  orbspi = (int **) malloc(sizeof(int *) * num_subspaces);
  orbsym = (int **) malloc(sizeof(int *) * num_subspaces);
  numorbs = (int *) malloc(num_subspaces * sizeof(int));
  for(i=0; i < num_subspaces; i++) {
    orbspi[i] = (int *) malloc(sizeof(int) * nirreps);
    tmparray = va_arg(ap, int *);
    for(j=0; j < nirreps; j++) orbspi[i][j] = tmparray[j];

    /* Compute the number of orbitals in this subspace */
    numorbs[i] = 0;
    for(h=0; h < nirreps; h++)
      numorbs[i] += orbspi[i][h];

    orbsym[i] = (int *) malloc(sizeof(int) * numorbs[i]);
    tmparray = va_arg(ap, int *);
    for(j=0; j < numorbs[i]; j++) orbsym[i][j] = tmparray[j];
  }
  va_end(ap);
  this_dpd->orbspi = orbspi;
  this_dpd->orbsym = orbsym;
  this_dpd->numorbs = numorbs;

  /* Compute the orbital offset arrays */
  orboff = (int **) malloc(num_subspaces * sizeof(int *));
  for(i=0; i < num_subspaces; i++) {
    orboff[i] = init_int_array(nirreps);
    for(j=1; j < nirreps; j++)
      orboff[i][j] = orboff[i][j-1] + orbspi[i][j-1];
  }
  this_dpd->orboff = orboff;

  /* Compute the number of bra or ket index combinations */
  num_pairs = (num_subspaces * (num_subspaces - 1)) + (5 * num_subspaces);
  this_dpd->num_pairs = num_pairs;

  /* Set up the pair structs for later use */
  pairs = (dpdpair *) malloc(num_pairs * sizeof(dpdpair));

  /* Build the row/column dimension arrays */
  pairtot = init_int_matrix(num_pairs, nirreps);
  this_dpd->pairtot = pairtot;
  /* Loop over the groups of five "diagonal" pairs */
  for(i=0; i < num_subspaces; i++) {

    pairs[5*i].left_orbspi = orbspi[i];
    pairs[5*i].left_orbsym = orbsym[i];
    pairs[5*i].left_orboff = orboff[i];
    pairs[5*i].right_orbspi = orbspi[i];
    pairs[5*i].right_orbsym = orbsym[i];
    pairs[5*i].right_orboff = orboff[i];
    pairs[5*i].permlr = 0;
    pairs[5*i].ler = 0;

    pairs[5*i+1].left_orbspi = orbspi[i];
    pairs[5*i+1].left_orbsym = orbsym[i];
    pairs[5*i+1].left_orboff = orboff[i];
    pairs[5*i+1].right_orbspi = orbspi[i];
    pairs[5*i+1].right_orbsym = orbsym[i];
    pairs[5*i+1].right_orboff = orboff[i];
    pairs[5*i+1].permlr = 1;
    pairs[5*i+1].ler = 0;

    pairs[5*i+2].left_orbspi = orbspi[i];
    pairs[5*i+2].left_orbsym = orbsym[i];
    pairs[5*i+2].left_orboff = orboff[i];
    pairs[5*i+2].right_orbspi = orbspi[i];
    pairs[5*i+2].right_orbsym = orbsym[i];
    pairs[5*i+2].right_orboff = orboff[i];
    pairs[5*i+2].permlr = -1;
    pairs[5*i+2].ler = 0;

    pairs[5*i+3].left_orbspi = orbspi[i];
    pairs[5*i+3].left_orbsym = orbsym[i];
    pairs[5*i+3].left_orboff = orboff[i];
    pairs[5*i+3].right_orbspi = orbspi[i];
    pairs[5*i+3].right_orbsym = orbsym[i];
    pairs[5*i+3].right_orboff = orboff[i];
    pairs[5*i+3].permlr = 1;
    pairs[5*i+3].ler = 1;

    pairs[5*i+4].left_orbspi = orbspi[i];
    pairs[5*i+4].left_orbsym = orbsym[i];
    pairs[5*i+4].left_orboff = orboff[i];
    pairs[5*i+4].right_orbspi = orbspi[i];
    pairs[5*i+4].right_orbsym = orbsym[i];
    pairs[5*i+4].right_orboff = orboff[i];
    pairs[5*i+4].permlr = -1;
    pairs[5*i+4].ler = 1;
  		  
    for(j=0; j < nirreps; j++) 
      for(k=0; k < nirreps; k++) {
	l_irrep = dp[j][k][0]; r_irrep = dp[j][k][1];

	/* orbspi,orbspi */
	pairtot[5*i][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];

	if(l_irrep > r_irrep) {
	  /* orbspi < orbspi, +1 */
	  pairtot[5*i+1][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
	  /* orbspi < orbspi, -1 */
	  pairtot[5*i+2][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
	  /* orbspi <= orbspi, +1 */
	  pairtot[5*i+3][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
	  /* orbspi <= orbspi, -1 */
	  pairtot[5*i+4][j] += orbspi[i][l_irrep] * orbspi[i][r_irrep];
	}
	else if(l_irrep == r_irrep) {
	  /* orbspi < orbspi, +1 */
	  pairtot[5*i+1][j] +=
	    (orbspi[i][l_irrep] * (orbspi[i][l_irrep]-1))/2;
	  /* orbspi < orbspi, -1 */
	  pairtot[5*i+2][j] +=
	    (orbspi[i][l_irrep] * (orbspi[i][l_irrep]-1))/2;
	  /* orbspi <= orbspi, +1 */
	  pairtot[5*i+3][j] +=
	    (orbspi[i][l_irrep] * (orbspi[i][l_irrep]+1))/2;
	  /* orbspi <= orbspi, -1 */
	  pairtot[5*i+4][j] +=
	    (orbspi[i][l_irrep] * (orbspi[i][l_irrep]+1))/2;
	}
      }
  }
	  
  /* Loop over the remaining "off diagonal" pairs */
  for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++) 
    for(j=i+1; j < num_subspaces; j++,cnt+=2) {

      pairs[cnt].left_orbspi = orbspi[i];
      pairs[cnt].left_orbsym = orbsym[i];
      pairs[cnt].left_orboff = orboff[i];
      pairs[cnt].right_orbspi = orbspi[j];
      pairs[cnt].right_orbsym = orbsym[j];
      pairs[cnt].right_orboff = orboff[j];
      pairs[cnt].permlr = 0;
      pairs[cnt].ler = 0;

      pairs[cnt+1].left_orbspi = orbspi[j];
      pairs[cnt+1].left_orbsym = orbsym[j];
      pairs[cnt+1].left_orboff = orboff[j];
      pairs[cnt+1].right_orbspi = orbspi[i];
      pairs[cnt+1].right_orbsym = orbsym[i];
      pairs[cnt+1].right_orboff = orboff[i];
      pairs[cnt+1].permlr = 0;
      pairs[cnt+1].ler = 0;
  	  
      for(k=0; k < nirreps; k++)
	for(l=0; l < nirreps; l++) {
	      
	  l_irrep = dp[k][l][0]; r_irrep = dp[k][l][1];

	  /* orbspi[i],orbspi[j] */
	  pairtot[cnt][k] += orbspi[i][l_irrep] * orbspi[j][r_irrep];
	  /* orbspi[j],orbspi[i] */
	  pairtot[cnt+1][k] += orbspi[j][l_irrep] * orbspi[i][r_irrep];

	}
    }

  /* Temporary check until I'm sure I'm doing this right */
  if(num_pairs != cnt) { printf("Error in dpd_init()!\n"); exit(PSI_RETURN_FAILURE); }

  /* Build the row/column index lookup arrays */
  pairidx = (int ***) malloc(num_pairs * sizeof(int **));
  pairorb = (int ****) malloc(num_pairs * sizeof(int ***));
  this_dpd->pairidx = pairidx;
  this_dpd->pairorb = pairorb;
  count = init_int_array(nirreps);
  /* Loop over the groups of five "diagonal" pairs */
  for(i=0; i < num_subspaces; i++) {

    for(l=0; l < 5; l++) {
      pairidx[5*i+l] = init_int_matrix(numorbs[i],numorbs[i]);
      for(j=0; j < numorbs[i]; j++)
	for(k=0; k < numorbs[i]; k++) 
	  pairidx[5*i+l][j][k] = -1;
      
      pairorb[5*i+l] = (int ***) malloc(nirreps * sizeof(int **));
      for(j=0; j < nirreps; j++) {
	pairorb[5*i+l][j] =
	  pairtot[5*i+l][j] ? init_int_matrix(pairtot[5*i+l][j],2) : NULL;
	for(k=0; k < pairtot[5*i+l][j]; k++) {
	  pairorb[5*i+l][j][k][0] = -1;
	  pairorb[5*i+l][j][k][1] = -1;
	}
      }
    }

    zero_int_array(count,nirreps);

    /* orbspi[i],orbspi[i] */
    for(j=0; j < nirreps; j++)
      for(k=0; k < nirreps; k++) {
	h0 = dp[j][k][0]; h1 = dp[j][k][1];
	offset1 = orboff[i][h0];  offset2 = orboff[i][h1];
	for(p=0; p < orbspi[i][h0]; p++)
	  for(q=0; q < orbspi[i][h1]; q++) {
	    pairorb[5*i][j][count[j]][0] = p+offset1;
	    pairorb[5*i][j][count[j]][1] = q+offset2;
	    pairidx[5*i][p+offset1][q+offset2] = count[j]++;
	  }
      }

    zero_int_array(count, nirreps);

    /* orbspi[i] < orbspi[i], +1 */
    for(j=0; j < nirreps; j++)
      for(k=0; k < nirreps; k++) {
	h0 = dp[j][k][0]; h1 = dp[j][k][1];
	offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	if(h0 == h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q < p; q++) {
	      pairorb[5*i+1][j][count[j]][0] = p+offset1;
	      pairorb[5*i+1][j][count[j]][1] = q+offset2;
	      pairidx[5*i+1][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+1][q+offset2][p+offset1] = count[j]++;
	    }
	}
	else if(h0 > h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q < orbspi[i][h1]; q++) {
	      pairorb[5*i+1][j][count[j]][0] = p+offset1;
	      pairorb[5*i+1][j][count[j]][1] = q+offset2;
	      pairidx[5*i+1][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+1][q+offset2][p+offset1] = count[j]++;
	    }
	}
      }

    zero_int_array(count, nirreps);

    /* orbspi[i] < orbspi[i], -1 */
    for(j=0; j < nirreps; j++)
      for(k=0; k < nirreps; k++) {
	h0 = dp[j][k][0]; h1 = dp[j][k][1];
	offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	if(h0 == h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q < p; q++) {
	      pairorb[5*i+2][j][count[j]][0] = p+offset1;
	      pairorb[5*i+2][j][count[j]][1] = q+offset2;
	      pairidx[5*i+2][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+2][q+offset2][p+offset1] = count[j]++;
	    }
	}
	else if(h0 > h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q < orbspi[i][h1]; q++) {
	      pairorb[5*i+2][j][count[j]][0] = p+offset1;
	      pairorb[5*i+2][j][count[j]][1] = q+offset2;
	      pairidx[5*i+2][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+2][q+offset2][p+offset1] = count[j]++;
	    }
	}
      }

    zero_int_array(count, nirreps);

    /* orbspi[i] <= orbspi[i], +1 */
    for(j=0; j < nirreps; j++)
      for(k=0; k < nirreps; k++) {
	h0 = dp[j][k][0]; h1 = dp[j][k][1];
	offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	if(h0 == h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q <= p; q++) {
	      pairorb[5*i+3][j][count[j]][0] = p+offset1;
	      pairorb[5*i+3][j][count[j]][1] = q+offset2;
	      pairidx[5*i+3][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+3][q+offset2][p+offset1] = count[j]++;
	    }
	}
	else if(h0 > h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q < orbspi[i][h1]; q++) {
	      pairorb[5*i+3][j][count[j]][0] = p+offset1;
	      pairorb[5*i+3][j][count[j]][1] = q+offset2;
	      pairidx[5*i+3][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+3][q+offset2][p+offset1] = count[j]++;
	    }
	}
      }

    zero_int_array(count, nirreps);

    /* orbspi[i] <= orbspi[i], -1 */
    for(j=0; j < nirreps; j++)
      for(k=0; k < nirreps; k++) {
	h0 = dp[j][k][0]; h1 = dp[j][k][1];
	offset1 = orboff[i][h0]; offset2 = orboff[i][h1];
	if(h0 == h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q <= p; q++) {
	      pairorb[5*i+4][j][count[j]][0] = p+offset1;
	      pairorb[5*i+4][j][count[j]][1] = q+offset2;
	      pairidx[5*i+4][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+4][q+offset2][p+offset1] = count[j]++;
	    }
	}
	else if(h0 > h1) {
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q < orbspi[i][h1]; q++) {
	      pairorb[5*i+4][j][count[j]][0] = p+offset1;
	      pairorb[5*i+4][j][count[j]][1] = q+offset2;
	      pairidx[5*i+4][p+offset1][q+offset2] = count[j];
	      pairidx[5*i+4][q+offset2][p+offset1] = count[j]++;
	    }
	}
      }

  }

  /* Loop over the remaining "off diagonal" pairs */
  for(i=0,cnt=5*num_subspaces; i < num_subspaces; i++) {
    for(j=i+1; j < num_subspaces; j++,cnt+=2) {

      pairidx[cnt] = init_int_matrix(numorbs[i],numorbs[j]);
      for(k=0; k < numorbs[i]; k++)
	for(l=0; l < numorbs[j]; l++)
	  pairidx[cnt][k][l] = -1;
      pairidx[cnt+1] = init_int_matrix(numorbs[j],numorbs[i]);
      for(k=0; k < numorbs[j]; k++)
	for(l=0; l < numorbs[i]; l++)
	  pairidx[cnt+1][k][l] = -1;

      pairorb[cnt] = (int ***) malloc(nirreps * sizeof(int **));
      pairorb[cnt+1] = (int ***) malloc(nirreps * sizeof(int **));
      for(k=0; k < nirreps; k++) {
	pairorb[cnt][k] =
	  pairtot[cnt][k] ? init_int_matrix(pairtot[cnt][k],2) : NULL;
	pairorb[cnt+1][k] =
	  pairtot[cnt+1][k] ? init_int_matrix(pairtot[cnt+1][k],2) : NULL;
	for(l=0; l < pairtot[cnt][k]; l++) {
	  pairorb[cnt][k][l][0] = -1;
	  pairorb[cnt][k][l][1] = -1;
	}
	for(l=0; l < pairtot[cnt+1][k]; l++) {
	  pairorb[cnt+1][k][l][0] = -1;
	  pairorb[cnt+1][k][l][1] = -1;
	}

      }

      zero_int_array(count, nirreps);
	  
      for(k=0; k < nirreps; k++)
	for(l=0; l < nirreps; l++) {
	  h0 = dp[k][l][0]; h1 = dp[k][l][1];
	  offset1 = orboff[i][h0];  offset2 = orboff[j][h1];
	  for(p=0; p < orbspi[i][h0]; p++)
	    for(q=0; q < orbspi[j][h1]; q++) {
	      pairorb[cnt][k][count[k]][0] = p+offset1;
	      pairorb[cnt][k][count[k]][1] = q+offset2;
	      pairidx[cnt][p+offset1][q+offset2] = count[k]++;
	    }
	}

      zero_int_array(count, nirreps);

      for(k=0; k < nirreps; k++)
	for(l=0; l < nirreps; l++) {
	  h0 = dp[k][l][0]; h1 = dp[k][l][1];
	  offset1 = orboff[j][h0];  offset2 = orboff[i][h1];
	  for(p=0; p < orbspi[j][h0]; p++)
	    for(q=0; q < orbspi[i][h1]; q++) {
	      pairorb[cnt+1][k][count[k]][0] = p+offset1;
	      pairorb[cnt+1][k][count[k]][1] = q+offset2;
	      pairidx[cnt+1][p+offset1][q+offset2] = count[k]++;
	    }
	}
    }
  }

  /* Temporary check until I'm sure I'm doing this right */
  if(num_pairs != cnt) { printf("Error in dpd_init()!\n"); exit(PSI_RETURN_FAILURE); }

  /* Now generate the global list of DPD parameters */
  this_dpd->params4 =
    (dpdparams4 **) malloc(num_pairs*sizeof(dpdparams4 *));
  for(i=0; i < num_pairs; i++)
    this_dpd->params4[i] =
      (dpdparams4 *) malloc(num_pairs*sizeof(dpdparams4));

  for(i=0; i < num_pairs; i++) {
    for(j=0; j < num_pairs; j++) {
      this_dpd->params4[i][j].nirreps = nirreps;

      this_dpd->params4[i][j].pqnum = i;
      this_dpd->params4[i][j].rsnum = j;

      this_dpd->params4[i][j].rowtot = pairtot[i];
      this_dpd->params4[i][j].coltot = pairtot[j];

      this_dpd->params4[i][j].rowidx = pairidx[i];
      this_dpd->params4[i][j].colidx = pairidx[j];

      this_dpd->params4[i][j].roworb = pairorb[i];
      this_dpd->params4[i][j].colorb = pairorb[j];

      this_dpd->params4[i][j].ppi = pairs[i].left_orbspi;
      this_dpd->params4[i][j].qpi = pairs[i].right_orbspi;

      this_dpd->params4[i][j].rpi = pairs[j].left_orbspi;
      this_dpd->params4[i][j].spi = pairs[j].right_orbspi;

      this_dpd->params4[i][j].psym = pairs[i].left_orbsym;
      this_dpd->params4[i][j].qsym = pairs[i].right_orbsym;

      this_dpd->params4[i][j].rsym = pairs[j].left_orbsym;
      this_dpd->params4[i][j].ssym = pairs[j].right_orbsym;

      this_dpd->params4[i][j].poff = pairs[i].left_orboff;
      this_dpd->params4[i][j].qoff = pairs[i].right_orboff;

      this_dpd->params4[i][j].roff = pairs[j].left_orboff;
      this_dpd->params4[i][j].soff = pairs[j].right_orboff;

      this_dpd->params4[i][j].perm_pq = pairs[i].permlr;
      this_dpd->params4[i][j].perm_rs = pairs[j].permlr;
      this_dpd->params4[i][j].peq = pairs[i].ler;
      this_dpd->params4[i][j].res = pairs[j].ler;

    }
  }

  /* generate the start13 lookup array */
  for(i=0; i < num_pairs; i++) {
    for(j=0; j < num_pairs; j++) {

      for(h=0,nump=0; h < nirreps; h++) nump += this_dpd->params4[i][j].ppi[h];

      this_dpd->params4[i][j].start13 = init_int_matrix(nirreps, nump);

      for(h=0; h < nirreps; h++) { /* h = Gamma_pq */
	for(p=0; p < nump; p++) this_dpd->params4[i][j].start13[h][p] = -1; /* error checking */
	nrows = 0;
	for(Gp=0; Gp < nirreps; Gp++) { /* Gamma_p */
	  for(p=0; p < this_dpd->params4[i][j].ppi[Gp]; p++) {
	    offset = this_dpd->params4[i][j].poff[Gp];
	    if(this_dpd->params4[i][j].qpi[Gp^h])
	      this_dpd->params4[i][j].start13[h][offset + p] = nrows;
	    nrows += this_dpd->params4[i][j].qpi[Gp^h];
	  } /* p */
	} /* Gp */
      } /* h */

    } /* j */
  } /* i */

  /* Now generate the global list of one-electron DPD parameters */
  orbidx2 = (int **) malloc(num_subspaces*sizeof(int *));
  for(i=0; i < num_subspaces; i++) {
    orbidx2[i] = init_int_array(numorbs[i]);
    for(j=0; j < numorbs[i]; j++)
      orbidx2[i][j] = -1;
  }
  this_dpd->orbidx2 = orbidx2;

  orbs2 = (int ***) malloc(num_subspaces*sizeof(int **));
  for(i=0;i < num_subspaces; i++) {
    orbs2[i] = (int **) malloc(nirreps*sizeof(int *));
    for(j=0; j < nirreps; j++) {
      orbs2[i][j] = orbspi[i][j] ? init_int_array(orbspi[i][j]) : NULL;
      for(k=0; k < orbspi[i][j]; k++)
	orbs2[i][j][k] = -1;
    }
  }
  this_dpd->orbs2 = orbs2;

  for(i=0; i < num_subspaces; i++) {
    zero_int_array(count, nirreps);
    for(j=0; j < nirreps; j++) {
      offset1 = orboff[i][j];
      for(p=0; p < orbspi[i][j]; p++) {
	orbs2[i][j][count[j]] = p+offset1;
	orbidx2[i][p+offset1] = count[j]++;
      }
    }
  }

  this_dpd->params2 =
    (dpdparams2 **) malloc(num_subspaces*sizeof(dpdparams2 *));
  for(i=0; i < num_subspaces; i++)
    this_dpd->params2[i] =
      (dpdparams2 *) malloc(num_subspaces*sizeof(dpdparams2));

  for(i=0,cnt=0; i < num_subspaces; i++) {
    for(j=0; j < num_subspaces; j++,cnt++) {
      this_dpd->params2[i][j].nirreps = nirreps;

      this_dpd->params2[i][j].pnum = i;
      this_dpd->params2[i][j].qnum = j;

      this_dpd->params2[i][j].rowtot = orbspi[i];
      this_dpd->params2[i][j].coltot = orbspi[j];

      this_dpd->params2[i][j].rowidx = orbidx2[i];
      this_dpd->params2[i][j].colidx = orbidx2[j];

      this_dpd->params2[i][j].roworb = orbs2[i];
      this_dpd->params2[i][j].colorb = orbs2[j];

      this_dpd->params2[i][j].ppi = orbspi[i];
      this_dpd->params2[i][j].qpi = orbspi[j];

      this_dpd->params2[i][j].poff = orboff[i];
      this_dpd->params2[i][j].qoff = orboff[j];

      this_dpd->params2[i][j].psym = orbsym[i];
      this_dpd->params2[i][j].qsym = orbsym[j];
    }
  }

  free(count);
  
  free(pairs);

  for(h=0; h < nirreps; h++) free_int_matrix(dp[h]);
  free(dp);

  /* Set the default DPD set to the current one */
  dpd_default = dpd_num;

  /* Init the Cache Linked Lists */
  dpd_file2_cache_init();
  dpd_file4_cache_init();

  return 0;
}


} // namespace psi
