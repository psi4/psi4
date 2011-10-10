/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <ccfiles.h>
#include <psifiles.h>
#include <exception.h>

namespace psi {
extern FILE *outfile;
namespace ccenergy {

void cache_abcd_rhf(int **cachelist);
void cache_iabc_rhf(int **cachelist);
void cache_ijab_rhf(int **cachelist);
void cache_iajb_rhf(int **cachelist);
void cache_ijka_rhf(int **cachelist);
void cache_ijkl_rhf(int **cachelist);

void cache_abcd_uhf(int **cachelist);
void cache_iabc_uhf(int **cachelist);
void cache_ijab_uhf(int **cachelist);
void cache_iajb_uhf(int **cachelist);
void cache_ijka_uhf(int **cachelist);
void cache_ijkl_uhf(int **cachelist);

int **cacheprep_uhf(int level, int *cachefiles)
{
  int **cachelist;

  /* The listing of CC files whose entries may be cached */
  cachefiles[CC_AINTS] = 1;
  cachefiles[CC_CINTS] = 1;
  cachefiles[CC_DINTS] = 1;
  cachefiles[CC_EINTS] = 1;
  cachefiles[CC_DENOM] = 1;
  cachefiles[CC_TAMPS] = 1;
  cachefiles[CC_LAMPS] = 1;
  cachefiles[CC_HBAR] = 1;

  /* The listing of DPD patterns which may be cached */
  cachelist = init_int_matrix(32,32);

  if(level == 0) return cachelist;
  else if(level == 1) {

      /*** Cache oooo and ooov ***/
      cache_ijkl_uhf(cachelist);
      cache_ijka_uhf(cachelist);

      return cachelist;
    }
  else if(level == 2) {

      /*** Cache oooo, ooov, oovv, and ovov ***/
      cache_ijkl_uhf(cachelist);
      cache_ijka_uhf(cachelist);
      cache_ijab_uhf(cachelist);
      cache_iajb_uhf(cachelist);

      return cachelist;
    }
  else if(level == 3) {

      /*** Cache, oooo, oov, oovv, ovov, and ovvv ***/

      cache_ijkl_uhf(cachelist);
      cache_ijka_uhf(cachelist);
      cache_ijab_uhf(cachelist);
      cache_iajb_uhf(cachelist);
      cache_iabc_uhf(cachelist);

      return cachelist;
    }
  else if(level == 4) {

      /*** Cache everything ***/
      cache_ijkl_uhf(cachelist);
      cache_ijka_uhf(cachelist);
      cache_ijab_uhf(cachelist);
      cache_iajb_uhf(cachelist);
      cache_iabc_uhf(cachelist);
      cache_abcd_uhf(cachelist);

      return cachelist;
    }
  else {
    printf("Error: invalid cache level!\n");
    throw InputException("Invalid cache level!", "CACHELEV", level, __FILE__, __LINE__);
  }
}

int **cacheprep_rhf(int level, int *cachefiles)
{
  int **cachelist;

  /* The listing of CC files whose entries may be cached */
  cachefiles[CC_AINTS] = 1;
  cachefiles[CC_CINTS] = 1;
  cachefiles[CC_DINTS] = 1;
  cachefiles[CC_EINTS] = 1;
  cachefiles[CC_DENOM] = 1;
  cachefiles[CC_TAMPS] = 1;
  cachefiles[CC_LAMPS] = 1;
  cachefiles[CC_HBAR] = 1;

  /* The listing of DPD patterns which may be cached */
  cachelist = init_int_matrix(12,12);

  if(level == 0) return cachelist;
  else if(level == 1) {

      /*** Cache oooo and ooov ***/
      cache_ijkl_rhf(cachelist);
      cache_ijka_rhf(cachelist);

      return cachelist;
    }
  else if(level == 2) {

      /*** Cache oooo, ooov, oovv, and ovov ***/
      cache_ijkl_rhf(cachelist);
      cache_ijka_rhf(cachelist);
      cache_ijab_rhf(cachelist);
      cache_iajb_rhf(cachelist);

      return cachelist;
    }
  else if(level == 3) {

      /*** Cache, oooo, oov, oovv, ovov, and ovvv ***/

      cache_ijkl_rhf(cachelist);
      cache_ijka_rhf(cachelist);
      cache_ijab_rhf(cachelist);
      cache_iajb_rhf(cachelist);
      cache_iabc_rhf(cachelist);

      return cachelist;
    }
  else if(level == 4) {

      /*** Cache everything ***/
      cache_ijkl_rhf(cachelist);
      cache_ijka_rhf(cachelist);
      cache_ijab_rhf(cachelist);
      cache_iajb_rhf(cachelist);
      cache_iabc_rhf(cachelist);
      cache_abcd_rhf(cachelist);

      return cachelist;
    }
  else {
    printf("Error: invalid cache level!\n");
    throw InputException("Invalid cache level!", "CACHELEV", level, __FILE__, __LINE__);
  }
}

void cache_abcd_uhf(int **cachelist)
{
  /* <ab|cd> */
  cachelist[5][5] = 1;
  cachelist[5][6] = 1;
  cachelist[5][7] = 1;
  cachelist[5][8] = 1;
  cachelist[5][9] = 1;
  cachelist[6][5] = 1;
  cachelist[6][6] = 1;
  cachelist[6][7] = 1;
  cachelist[6][8] = 1;
  cachelist[6][9] = 1;
  cachelist[7][5] = 1;
  cachelist[7][6] = 1;
  cachelist[7][7] = 1;
  cachelist[7][8] = 1;
  cachelist[7][9] = 1;
  cachelist[8][5] = 1;
  cachelist[8][6] = 1;
  cachelist[8][7] = 1;
  cachelist[8][8] = 1;
  cachelist[8][9] = 1;
  cachelist[9][5] = 1;
  cachelist[9][6] = 1;
  cachelist[9][7] = 1;
  cachelist[9][8] = 1;
  cachelist[9][9] = 1;
  /* <AB|CD> */
  cachelist[15][15] = 1;
  cachelist[15][16] = 1;
  cachelist[15][17] = 1;
  cachelist[15][18] = 1;
  cachelist[15][19] = 1;
  cachelist[16][15] = 1;
  cachelist[16][16] = 1;
  cachelist[16][17] = 1;
  cachelist[16][18] = 1;
  cachelist[16][19] = 1;
  cachelist[17][15] = 1;
  cachelist[17][16] = 1;
  cachelist[17][17] = 1;
  cachelist[17][18] = 1;
  cachelist[17][19] = 1;
  cachelist[18][15] = 1;
  cachelist[18][16] = 1;
  cachelist[18][17] = 1;
  cachelist[18][18] = 1;
  cachelist[18][19] = 1;
  cachelist[19][15] = 1;
  cachelist[19][16] = 1;
  cachelist[19][17] = 1;
  cachelist[19][18] = 1;
  cachelist[19][19] = 1;
  /* <Ab|Cd> */
  cachelist[28][28] = 1;
  cachelist[29][29] = 1;
  cachelist[28][29] = 1;
  cachelist[29][28] = 1;
}

void cache_abcd_rhf(int **cachelist)
{
  /* <ab|cd> */
  cachelist[5][5] = 1;
  cachelist[5][6] = 1;
  cachelist[5][7] = 1;
  cachelist[5][8] = 1;
  cachelist[5][9] = 1;
  cachelist[6][5] = 1;
  cachelist[6][6] = 1;
  cachelist[6][7] = 1;
  cachelist[6][8] = 1;
  cachelist[6][9] = 1;
  cachelist[7][5] = 1;
  cachelist[7][6] = 1;
  cachelist[7][7] = 1;
  cachelist[7][8] = 1;
  cachelist[7][9] = 1;
  cachelist[8][5] = 1;
  cachelist[8][6] = 1;
  cachelist[8][7] = 1;
  cachelist[8][8] = 1;
  cachelist[8][9] = 1;
  cachelist[9][5] = 1;
  cachelist[9][6] = 1;
  cachelist[9][7] = 1;
  cachelist[9][8] = 1;
  cachelist[9][9] = 1;
}

void cache_iabc_rhf(int **cachelist)
{
  /* <ia|bc> */
  cachelist[10][5] = 1;
  cachelist[10][6] = 1;
  cachelist[10][7] = 1;
  cachelist[10][8] = 1;
  cachelist[10][9] = 1;
  cachelist[11][5] = 1;
  cachelist[11][6] = 1;
  cachelist[11][7] = 1;
  cachelist[11][8] = 1;
  cachelist[11][9] = 1;
  /* <ab|ci> */
  cachelist[5][10] = 1;
  cachelist[5][11] = 1;
  cachelist[6][10] = 1;
  cachelist[6][11] = 1;
  cachelist[7][10] = 1;
  cachelist[7][11] = 1;
  cachelist[8][10] = 1;
  cachelist[8][11] = 1;
  cachelist[9][10] = 1;
  cachelist[9][11] = 1;
}

void cache_iabc_uhf(int **cachelist)
{
  /* <IA|BC> */
  cachelist[20][5] = 1;
  cachelist[20][6] = 1;
  cachelist[20][7] = 1;
  cachelist[20][8] = 1;
  cachelist[20][9] = 1;
  cachelist[21][5] = 1;
  cachelist[21][6] = 1;
  cachelist[21][7] = 1;
  cachelist[21][8] = 1;
  cachelist[21][9] = 1;
  /* <AB|CI> */
  cachelist[5][20] = 1;
  cachelist[5][21] = 1;
  cachelist[6][20] = 1;
  cachelist[6][21] = 1;
  cachelist[7][20] = 1;
  cachelist[7][21] = 1;
  cachelist[8][20] = 1;
  cachelist[8][21] = 1;
  cachelist[9][20] = 1;
  cachelist[9][21] = 1;

  /* <ia|bc> */
  cachelist[30][15] = 1;
  cachelist[30][16] = 1;
  cachelist[30][17] = 1;
  cachelist[30][18] = 1;
  cachelist[30][19] = 1;
  cachelist[31][15] = 1;
  cachelist[31][16] = 1;
  cachelist[31][17] = 1;
  cachelist[31][18] = 1;
  cachelist[31][19] = 1;
  /* <ab|ci> */
  cachelist[15][30] = 1;
  cachelist[15][31] = 1;
  cachelist[16][30] = 1;
  cachelist[16][31] = 1;
  cachelist[17][30] = 1;
  cachelist[17][31] = 1;
  cachelist[18][30] = 1;
  cachelist[18][31] = 1;
  cachelist[19][30] = 1;
  cachelist[19][31] = 1;

  /* <Ia|Bc> */
  cachelist[24][28] = 1;
  cachelist[24][29] = 1;
  cachelist[25][28] = 1;
  cachelist[25][29] = 1;

  /* <Ab|Ci> */
  cachelist[28][24] = 1;
  cachelist[28][25] = 1;
  cachelist[29][24] = 1;
  cachelist[29][25] = 1;
}

void cache_ijab_rhf(int **cachelist)
{
  /* <ij|ab> */
  cachelist[0][5] = 1;
  cachelist[0][6] = 1;
  cachelist[0][7] = 1;
  cachelist[0][8] = 1;
  cachelist[0][9] = 1;
  cachelist[1][5] = 1;
  cachelist[1][6] = 1;
  cachelist[1][7] = 1;
  cachelist[1][8] = 1;
  cachelist[1][9] = 1;
  cachelist[2][5] = 1;
  cachelist[2][6] = 1;
  cachelist[2][7] = 1;
  cachelist[2][8] = 1;
  cachelist[2][9] = 1;
  cachelist[3][5] = 1;
  cachelist[3][6] = 1;
  cachelist[3][7] = 1;
  cachelist[3][8] = 1;
  cachelist[3][9] = 1;
  cachelist[4][5] = 1;
  cachelist[4][6] = 1;
  cachelist[4][7] = 1;
  cachelist[4][8] = 1;
  cachelist[4][9] = 1;
  /* <ab|ij> */
  cachelist[5][0] = 1;
  cachelist[5][1] = 1;
  cachelist[5][2] = 1;
  cachelist[5][3] = 1;
  cachelist[5][4] = 1;
  cachelist[6][0] = 1;
  cachelist[6][1] = 1;
  cachelist[6][2] = 1;
  cachelist[6][3] = 1;
  cachelist[6][4] = 1;
  cachelist[7][0] = 1;
  cachelist[7][1] = 1;
  cachelist[7][2] = 1;
  cachelist[7][3] = 1;
  cachelist[7][4] = 1;
  cachelist[8][0] = 1;
  cachelist[8][1] = 1;
  cachelist[8][2] = 1;
  cachelist[8][3] = 1;
  cachelist[8][4] = 1;
  cachelist[9][0] = 1;
  cachelist[9][1] = 1;
  cachelist[9][2] = 1;
  cachelist[9][3] = 1;
  cachelist[9][4] = 1;
}

void cache_ijab_uhf(int **cachelist)
{
  /* <IJ|AB> */
  cachelist[0][5] = 1;
  cachelist[0][6] = 1;
  cachelist[0][7] = 1;
  cachelist[0][8] = 1;
  cachelist[0][9] = 1;
  cachelist[1][5] = 1;
  cachelist[1][6] = 1;
  cachelist[1][7] = 1;
  cachelist[1][8] = 1;
  cachelist[1][9] = 1;
  cachelist[2][5] = 1;
  cachelist[2][6] = 1;
  cachelist[2][7] = 1;
  cachelist[2][8] = 1;
  cachelist[2][9] = 1;
  cachelist[3][5] = 1;
  cachelist[3][6] = 1;
  cachelist[3][7] = 1;
  cachelist[3][8] = 1;
  cachelist[3][9] = 1;
  cachelist[4][5] = 1;
  cachelist[4][6] = 1;
  cachelist[4][7] = 1;
  cachelist[4][8] = 1;
  cachelist[4][9] = 1;
  /* <AB|IJ> */
  cachelist[5][0] = 1;
  cachelist[5][1] = 1;
  cachelist[5][2] = 1;
  cachelist[5][3] = 1;
  cachelist[5][4] = 1;
  cachelist[6][0] = 1;
  cachelist[6][1] = 1;
  cachelist[6][2] = 1;
  cachelist[6][3] = 1;
  cachelist[6][4] = 1;
  cachelist[7][0] = 1;
  cachelist[7][1] = 1;
  cachelist[7][2] = 1;
  cachelist[7][3] = 1;
  cachelist[7][4] = 1;
  cachelist[8][0] = 1;
  cachelist[8][1] = 1;
  cachelist[8][2] = 1;
  cachelist[8][3] = 1;
  cachelist[8][4] = 1;
  cachelist[9][0] = 1;
  cachelist[9][1] = 1;
  cachelist[9][2] = 1;
  cachelist[9][3] = 1;
  cachelist[9][4] = 1;

  /* <ij|ab> */
  cachelist[10][15] = 1;
  cachelist[10][16] = 1;
  cachelist[10][17] = 1;
  cachelist[10][18] = 1;
  cachelist[10][19] = 1;
  cachelist[11][15] = 1;
  cachelist[11][16] = 1;
  cachelist[11][17] = 1;
  cachelist[11][18] = 1;
  cachelist[11][19] = 1;
  cachelist[12][15] = 1;
  cachelist[12][16] = 1;
  cachelist[12][17] = 1;
  cachelist[12][18] = 1;
  cachelist[12][19] = 1;
  cachelist[13][15] = 1;
  cachelist[13][16] = 1;
  cachelist[13][17] = 1;
  cachelist[13][18] = 1;
  cachelist[13][19] = 1;
  cachelist[14][15] = 1;
  cachelist[14][16] = 1;
  cachelist[14][17] = 1;
  cachelist[14][18] = 1;
  cachelist[14][19] = 1;
  /* <ab|ij> */
  cachelist[15][10] = 1;
  cachelist[15][11] = 1;
  cachelist[15][12] = 1;
  cachelist[15][13] = 1;
  cachelist[15][14] = 1;
  cachelist[16][10] = 1;
  cachelist[16][11] = 1;
  cachelist[16][12] = 1;
  cachelist[16][13] = 1;
  cachelist[16][14] = 1;
  cachelist[17][10] = 1;
  cachelist[17][11] = 1;
  cachelist[17][12] = 1;
  cachelist[17][13] = 1;
  cachelist[17][14] = 1;
  cachelist[18][10] = 1;
  cachelist[18][11] = 1;
  cachelist[18][12] = 1;
  cachelist[18][13] = 1;
  cachelist[18][14] = 1;
  cachelist[19][10] = 1;
  cachelist[19][11] = 1;
  cachelist[19][12] = 1;
  cachelist[19][13] = 1;
  cachelist[19][14] = 1;

  /* <Ij|Ab> */
  cachelist[22][28] = 1;
  cachelist[23][28] = 1;
  cachelist[22][29] = 1;
  cachelist[23][29] = 1;
  /* <Ab|Ij> */
  cachelist[28][22] = 1;
  cachelist[28][23] = 1;
  cachelist[29][22] = 1;
  cachelist[29][23] = 1;
}

void cache_iajb_rhf(int **cachelist)
{
  /* <ia|jb> */
  cachelist[10][10] = 1;
  cachelist[10][11] = 1;
  cachelist[11][10] = 1;
  cachelist[11][11] = 1;
}

void cache_iajb_uhf(int **cachelist)
{
  /* <IA|JB> */
  cachelist[20][20] = 1;
  cachelist[20][21] = 1;
  cachelist[21][20] = 1;
  cachelist[21][21] = 1;
  /* <ia|jb> */
  cachelist[30][30] = 1;
  cachelist[30][31] = 1;
  cachelist[31][30] = 1;
  cachelist[31][31] = 1;
  /* <Ia|Jb> */
  cachelist[24][24] = 1;
  cachelist[24][25] = 1;
  cachelist[25][24] = 1;
  cachelist[25][25] = 1;
}

void cache_ijka_rhf(int **cachelist)
{
  /* <ij|ka> */
  cachelist[0][10] = 1;
  cachelist[0][11] = 1;
  cachelist[1][10] = 1;
  cachelist[1][11] = 1;
  cachelist[2][10] = 1;
  cachelist[2][11] = 1;
  cachelist[3][10] = 1;
  cachelist[3][11] = 1;
  cachelist[4][10] = 1;
  cachelist[4][11] = 1;
  /* <ia|jk> */
  cachelist[10][0] = 1;
  cachelist[10][1] = 1;
  cachelist[10][2] = 1;
  cachelist[10][3] = 1;
  cachelist[10][4] = 1;
  cachelist[11][0] = 1;
  cachelist[11][1] = 1;
  cachelist[11][2] = 1;
  cachelist[11][3] = 1;
  cachelist[11][4] = 1;
}

void cache_ijka_uhf(int **cachelist)
{
  /* <IJ|KA> */
  cachelist[0][20] = 1;
  cachelist[0][21] = 1;
  cachelist[1][20] = 1;
  cachelist[1][21] = 1;
  cachelist[2][20] = 1;
  cachelist[2][21] = 1;
  cachelist[3][20] = 1;
  cachelist[3][21] = 1;
  cachelist[4][20] = 1;
  cachelist[4][21] = 1;
  /* <IA|JK> */
  cachelist[20][0] = 1;
  cachelist[20][1] = 1;
  cachelist[20][2] = 1;
  cachelist[20][3] = 1;
  cachelist[20][4] = 1;
  cachelist[21][0] = 1;
  cachelist[21][1] = 1;
  cachelist[21][2] = 1;
  cachelist[21][3] = 1;
  cachelist[21][4] = 1;

  /* <ij|ka> */
  cachelist[10][30] = 1;
  cachelist[10][31] = 1;
  cachelist[11][30] = 1;
  cachelist[11][31] = 1;
  cachelist[12][30] = 1;
  cachelist[12][31] = 1;
  cachelist[13][30] = 1;
  cachelist[13][31] = 1;
  cachelist[14][30] = 1;
  cachelist[14][31] = 1;
  /* <ia|jk> */
  cachelist[30][10] = 1;
  cachelist[30][11] = 1;
  cachelist[30][12] = 1;
  cachelist[30][13] = 1;
  cachelist[30][14] = 1;
  cachelist[31][10] = 1;
  cachelist[31][11] = 1;
  cachelist[31][12] = 1;
  cachelist[31][13] = 1;
  cachelist[31][14] = 1;

  /* <Ij|Ka> */
  cachelist[22][24] = 1;
  cachelist[22][25] = 1;
  cachelist[23][24] = 1;
  cachelist[23][25] = 1;
  /* <Ka|Ij> */
  cachelist[24][22] = 1;
  cachelist[25][22] = 1;
  cachelist[24][23] = 1;
  cachelist[25][23] = 1;
}

void cache_ijkl_rhf(int **cachelist)
{
  /* <ij|kl> */
  cachelist[0][0] = 1;
  cachelist[0][1] = 1;
  cachelist[0][2] = 1;
  cachelist[0][3] = 1;
  cachelist[0][4] = 1;
  cachelist[1][0] = 1;
  cachelist[1][1] = 1;
  cachelist[1][2] = 1;
  cachelist[1][3] = 1;
  cachelist[1][4] = 1;
  cachelist[2][0] = 1;
  cachelist[2][1] = 1;
  cachelist[2][2] = 1;
  cachelist[2][3] = 1;
  cachelist[2][4] = 1;
  cachelist[3][0] = 1;
  cachelist[3][1] = 1;
  cachelist[3][2] = 1;
  cachelist[3][3] = 1;
  cachelist[3][4] = 1;
  cachelist[4][0] = 1;
  cachelist[4][1] = 1;
  cachelist[4][2] = 1;
  cachelist[4][3] = 1;
  cachelist[4][4] = 1;
}

void cache_ijkl_uhf(int **cachelist)
{
  /* <IJ|KL> */
  cachelist[0][0] = 1;
  cachelist[0][1] = 1;
  cachelist[0][2] = 1;
  cachelist[0][3] = 1;
  cachelist[0][4] = 1;
  cachelist[1][0] = 1;
  cachelist[1][1] = 1;
  cachelist[1][2] = 1;
  cachelist[1][3] = 1;
  cachelist[1][4] = 1;
  cachelist[2][0] = 1;
  cachelist[2][1] = 1;
  cachelist[2][2] = 1;
  cachelist[2][3] = 1;
  cachelist[2][4] = 1;
  cachelist[3][0] = 1;
  cachelist[3][1] = 1;
  cachelist[3][2] = 1;
  cachelist[3][3] = 1;
  cachelist[3][4] = 1;
  cachelist[4][0] = 1;
  cachelist[4][1] = 1;
  cachelist[4][2] = 1;
  cachelist[4][3] = 1;
  cachelist[4][4] = 1;
  /* <ij|kl> */
  cachelist[10][10] = 1;
  cachelist[10][11] = 1;
  cachelist[10][12] = 1;
  cachelist[10][13] = 1;
  cachelist[10][14] = 1;
  cachelist[11][10] = 1;
  cachelist[11][11] = 1;
  cachelist[11][12] = 1;
  cachelist[11][13] = 1;
  cachelist[11][14] = 1;
  cachelist[12][10] = 1;
  cachelist[12][11] = 1;
  cachelist[12][12] = 1;
  cachelist[12][13] = 1;
  cachelist[12][14] = 1;
  cachelist[13][10] = 1;
  cachelist[13][11] = 1;
  cachelist[13][12] = 1;
  cachelist[13][13] = 1;
  cachelist[13][14] = 1;
  cachelist[14][10] = 1;
  cachelist[14][11] = 1;
  cachelist[14][12] = 1;
  cachelist[14][13] = 1;
  cachelist[14][14] = 1;
  /* <Ij|Kl> */
  cachelist[22][22] = 1;
  cachelist[22][23] = 1;
  cachelist[23][22] = 1;
  cachelist[23][23] = 1;
}

void cachedone_uhf(int **cachelist)
{
  free_int_matrix(cachelist);
}

void cachedone_rhf(int **cachelist)
{
  free_int_matrix(cachelist);
}
}} // namespace psi::ccenergy
