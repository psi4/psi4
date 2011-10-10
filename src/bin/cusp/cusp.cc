/*! \defgroup CUSP Add a description of the group CUSP */

/*! \file
    \ingroup CUSP
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <psifiles.h>
#include <libiwl/iwl.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <ccfiles.h>
#include <physconst.h>

extern "C" {
  FILE *infile, *outfile;
  char *psi_file_prefix;
}

namespace psi { namespace cusp {

void init_io(int argc, char *argv[]);
void exit_io(void);
void compute_delta(double **delta, double x, double y, double z);
int **cacheprep(int level, int *cachefiles);
void local(void);


int cusp(int argc, char *argv[])
{
  int n, i, j, h, offset;
  int nmo, nso, natom, nirreps;
  int *clsdpi, *openpi, *orbspi;
  int nactive, row, col;
  int p,q,r,s,P,Q,R,S;
  int *qt_occ, *qt_vir;
  double **delta, **geom, **D;
  double **delta1, **delta2;
  double x, y, z, dens, average, theta;
  unsigned long int next;
  dpdbuf4 G, A;
  dpdfile2 g;
  long int memory;
  int *occpi, *virtpi, *occ_sym, *vir_sym, *occ_off, *vir_off;
  int *socc;
  int **cachelist, *cachefiles;
  FILE *data;
  FILE *data1;
  FILE *data2;
  double two_energy;
  double **scf_local;

  init_io(argc, argv);

  chkpt_init(PSIO_OPEN_OLD);
  nmo = chkpt_rd_nmo();
  nso = chkpt_rd_nso();
  geom = chkpt_rd_geom();
  natom = chkpt_rd_natom();
  nirreps = chkpt_rd_nirreps();
  clsdpi = chkpt_rd_clsdpi();
  openpi = chkpt_rd_openpi();
  orbspi = chkpt_rd_orbspi();
  chkpt_close();

  /*  local();
      throw PsiException("Elvis has left the building", __FILE__, __LINE__) */

  /*** Build the QT-ordered SCF density ***/
  /*
  D = block_matrix(nmo,nmo);
  for(h=offset=0; h < nirreps; h++) {

    for(i=0; i < clsdpi[h]; i++)
      D[i+offset][i+offset] = 2.0;

    for(i=clsdpi[h]; i < clsdpi[h] + openpi[h]; i++)
      D[i+offset][i+offset] = 1.0;

    offset += clsdpi[h] + openpi[h];
  }
  free(clsdpi);
  free(openpi);
  free(orbspi);
  */

  /*
    fprintf(outfile, "\nSCF Density (MO):\n");
    fprintf(outfile,   "-----------------\n");
    print_mat(D,nmo,nmo,outfile);
  */

  /*** Compute the SCF density at each nucleus ***/
  /*
  delta = block_matrix(nmo,nmo);
  for(n=0; n < natom; n++) {
    x = geom[n][0];
    y = geom[n][1];
    z = geom[n][2];

    compute_delta(delta, x, y, z);

    dens = 0.0;
    for(i=0; i < nmo; i++)
      for(j=0; j < nmo; j++)
        dens += delta[i][j] * D[i][j];

    fprintf(outfile,"\n\tdens(%5.3f,%5.3f,%5.3f) = %20.10f\n",x,y,z,dens);
  }
  free_block(delta);
  */

  psio_read_entry(CC_INFO, "No. of Active Orbitals", (char *) &(nactive),
                  sizeof(int));
  qt_occ = init_int_array(nactive);
  qt_vir = init_int_array(nactive);
  psio_read_entry(CC_INFO, "CC->QT Active Occ Order",
                  (char *) qt_occ, sizeof(int)*nactive);
  psio_read_entry(CC_INFO, "CC->QT Active Virt Order",
                  (char *) qt_vir, sizeof(int)*nactive);

  occpi = init_int_array(nirreps);
  virtpi = init_int_array(nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orbs Per Irrep",
                  (char *) occpi, sizeof(int)*nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orbs Per Irrep",
                  (char *) virtpi, sizeof(int)*nirreps);

  occ_sym = init_int_array(nactive);
  vir_sym = init_int_array(nactive);
  psio_read_entry(CC_INFO, "Active Occ Orb Symmetry",
                  (char *) occ_sym, sizeof(int)*nactive);
  psio_read_entry(CC_INFO, "Active Virt Orb Symmetry",
                  (char *) vir_sym, sizeof(int)*nactive);

  occ_off = init_int_array(nirreps);
  vir_off = init_int_array(nirreps);
  psio_read_entry(CC_INFO, "Active Occ Orb Offsets",
                  (char *) occ_off, sizeof(int)*nirreps);
  psio_read_entry(CC_INFO, "Active Virt Orb Offsets",
                  (char *) vir_off, sizeof(int)*nirreps);

  socc = init_int_array(nactive);
  psio_read_entry(CC_INFO, "Active Socc Orbital Boolean",
                  (char *) socc, sizeof(int)*nactive);


  /*** Get the correlated density ***/
  /*
  D = block_matrix(nactive,nactive);
  rfile(PSIF_MO_OPDM);
  next = 0;
  for(i=0; i < nactive; i++)
    wreadw(PSIF_MO_OPDM, (char *) D[i], sizeof(double)*nactive, next, &next);
  rclose(PSIF_MO_OPDM, 3);
  */

  /*
    fprintf(outfile, "\nCCSD Density (MO):\n");
    fprintf(outfile,   "------------------\n");
    print_mat(D,nmo,nmo,outfile);
  */

  /*** Compute the CCSD density at each nucleus ***/
  /*
  delta = block_matrix(nactive,nactive);
  for(n=0; n < natom; n++) {
    x = geom[n][0];
    y = geom[n][1];
    z = geom[n][2];

    compute_delta(delta, x, y, z);

    dens = 0.0;
    for(i=0; i < nactive; i++)
      for(j=0; j < nactive; j++)
        dens += delta[i][j] * D[i][j];

    fprintf(outfile,"\n\tdens(%5.3f,%5.3f,%5.3f) = %20.10f\n",x,y,z,dens);
  }
  free_block(delta);
  free_block(geom);
  free_block(D);
  */

  /* Initialize the DPD library */
  cachefiles = init_int_array(PSIO_MAXUNIT);
  cachelist = cacheprep(2, cachefiles);

  fndcor(&(memory),infile,outfile);

  dpd_init(0, nirreps, memory, 0, cachefiles, cachelist, NULL,
           2, occpi, occ_sym, virtpi, vir_sym);

  delta1 = block_matrix(nmo,nmo);
  delta2 = block_matrix(nmo,nmo);

  ffile(&data, "data", 0);

  compute_delta(delta1, 0.0,0.0,1.0);
  x = 0.0;
  for(theta=-_pi; theta <= _pi; theta += _pi/64.0) {

    dens = 0.0;
    two_energy = 0.0;

    y = sin(theta);
    z = cos(theta);

    compute_delta(delta2, x, y, z);

    dpd_buf4_init(&A, CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 0, 0, 0, 0, "GIjKl");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

      dpd_buf4_mat_irrep_init(&A, h);
      dpd_buf4_mat_irrep_rd(&A, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
        p = G.params->roworb[h][row][0]; P = qt_occ[p];
        q = G.params->roworb[h][row][1]; Q = qt_occ[q];
        for(col=0; col < G.params->coltot[h]; col++) {
          r = G.params->colorb[h][col][0]; R = qt_occ[r];
          s = G.params->colorb[h][col][1]; S = qt_occ[s];

          dens += G.matrix[h][row][col] * delta1[P][R] * delta2[Q][S];

          /* Reference contributions to two-electron density */
          if(P==R && Q==S) {
            if(!socc[P] && !socc[Q]) {
              /*	      dens += 2.0 * delta1[P][R] * delta2[Q][S]; */
              two_energy += 2.0 * A.matrix[h][row][col];
            }
            else if(socc[P] && !socc[Q]) {
              /*	      dens += 2.0 * delta1[P][R] * delta2[Q][S]; */
              two_energy += 2.0 * A.matrix[h][row][col];
            }
            else if(socc[P] && socc[Q]) {
              /*	      dens += 0.5 * delta1[P][R] * delta2[Q][S]; */
              two_energy += 0.5 * A.matrix[h][row][col];
            }
          }
          if (P==S && Q==R) {
            if(!socc[P] && !socc[Q]) {
              /*	      dens -= delta1[P][R] * delta2[Q][S]; */
              two_energy -= A.matrix[h][row][col];
            }
            else if(socc[P] && !socc[Q]) {
              /*	      dens -= delta1[P][R] * delta2[Q][S]; */
              two_energy -= A.matrix[h][row][col];
            }
            else if(socc[P] && socc[Q]) {
              /*	      dens -= 0.5 * delta1[P][R] * delta2[Q][S]; */
              two_energy -= 0.5 * A.matrix[h][row][col];
            }
          }

        }
      }
      dpd_buf4_mat_irrep_close(&A, h);
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&A);
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 5, 5, 5, 5, 0, "GAbCd");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
        p = G.params->roworb[h][row][0]; P = qt_vir[p];
        q = G.params->roworb[h][row][1]; Q = qt_vir[q];
        for(col=0; col < G.params->coltot[h]; col++) {
          r = G.params->colorb[h][col][0]; R = qt_vir[r];
          s = G.params->colorb[h][col][1]; S = qt_vir[s];

          dens += G.matrix[h][row][col] * delta1[P][R] * delta2[Q][S];

        }
      }
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 10, 0, 10, 0, "GIjKa");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
        p = G.params->roworb[h][row][0]; P = qt_occ[p];
        q = G.params->roworb[h][row][1]; Q = qt_occ[q];
        for(col=0; col < G.params->coltot[h]; col++) {
          r = G.params->colorb[h][col][0]; R = qt_occ[r];
          s = G.params->colorb[h][col][1]; S = qt_vir[s];

          average = 0.5 * (delta1[P][R] * delta2[Q][S] +
                           delta1[Q][S] * delta2[P][R]);

          dens += 2.0* G.matrix[h][row][col] * average;

        }
      }
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 11, 5, 11, 5, 0, "GCiAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);
      for(row=0; row < G.params->rowtot[h]; row++) {
        p = G.params->roworb[h][row][0]; P = qt_vir[p];
        q = G.params->roworb[h][row][1]; Q = qt_occ[q];
        for(col=0; col < G.params->coltot[h]; col++) {
          r = G.params->colorb[h][col][0]; R = qt_vir[r];
          s = G.params->colorb[h][col][1]; S = qt_vir[s];

          average = 0.5 * (delta1[P][R] * delta2[Q][S] +
                           delta1[Q][S] * delta2[P][R]);

          dens += 2.0* G.matrix[h][row][col] * average;
        }
      }
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 0, 5, 0, 5, 0, "GIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
        p = G.params->roworb[h][row][0]; P = qt_occ[p];
        q = G.params->roworb[h][row][1]; Q = qt_occ[q];
        for(col=0; col < G.params->coltot[h]; col++) {
          r = G.params->colorb[h][col][0]; R = qt_vir[r];
          s = G.params->colorb[h][col][1]; S = qt_vir[s];

          average = 0.5 * (delta1[P][R] * delta2[Q][S] +
                           delta1[Q][S] * delta2[P][R]);

          dens += G.matrix[h][row][col] * delta1[P][R] * delta2[Q][S];
        }
      }
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    dpd_buf4_init(&G, CC_GAMMA, 0, 10, 10, 10, 10, 0, "GIBJA");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&G, h);
      dpd_buf4_mat_irrep_rd(&G, h);

      for(row=0; row < G.params->rowtot[h]; row++) {
        p = G.params->roworb[h][row][0]; P = qt_occ[p];
        q = G.params->roworb[h][row][1]; Q = qt_vir[q];
        for(col=0; col < G.params->coltot[h]; col++) {
          r = G.params->colorb[h][col][0]; R = qt_occ[r];
          s = G.params->colorb[h][col][1]; S = qt_vir[s];

          average = 0.5 * (delta1[P][R] * delta2[Q][S] +
                           delta1[Q][S] * delta2[P][R]);

          dens += G.matrix[h][row][col] * average;
        }
      }
      dpd_buf4_mat_irrep_close(&G, h);
    }
    dpd_buf4_close(&G);

    fprintf(data, "%4.2f %20.10f\n", theta, dens);
  }

  printf("SCF Two-Electron Energy = %20.15f\n", two_energy);

  dpd_close(0);
  exit_io();
  exit(PSI_RETURN_SUCCESS);
}

void init_io(int argc, char *argv[])
{
  int i;

  tstart();

  /* Open all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_open(i,1);
}

void exit_io(void)
{
  int i;
  /* Close all dpd data files here */
  for(i=CC_MIN; i <= CC_MAX; i++) psio_close(i,1);

  tstop();
}

}} // namespace psi::cusp
