/*! \file
  \ingroup CUSP
  \brief Enter brief description of file here
  */
/*
 ** delta(): Compute the MO-basis delta-function for a given point.
 ** TDC, June 2001
 */

#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libchkpt/chkpt.h>
#include <libqt/qt.h>
#include <psifiles.h>


namespace psi { namespace cusp {

    int nmo, nao;
    double **scf, **u;

    void compute_phi(double *phi, double x, double y, double z);
    void setup_delta(void);

    void compute_delta(double **delta, double x, double y, double z)
    {
        int i, j;
        double *phi_ao, *phi_so, *phi_mo;

        setup_delta();

        phi_ao = init_array(nao);  /* AO function values */
        phi_so = init_array(nmo);  /* SO function values */
        phi_mo = init_array(nmo);  /* MO function values */

        compute_phi(phi_ao, x, y, z);

        /*  for(i=0; i < nao; i++) printf("%d %20.10f\n", i, phi_ao[i]); */

        /* Transform the basis function values to the MO basis */
        C_DGEMV('n', nmo, nao, 1.0, &(u[0][0]), nao, &(phi_ao[0]), 1,
                0.0, &(phi_so[0]), 1);

        C_DGEMV('t', nmo, nmo, 1.0, &(scf[0][0]), nmo, &(phi_so[0]), 1,
                0.0, &(phi_mo[0]), 1);

        /* for(i=0; i < nmo; i++) printf("%d %20.10f\n", i, phi_mo[i]); */


        /* Build the MO-basis delta function */
        for(i=0; i < nmo; i++)
            for(j=0; j < nmo; j++)
                delta[i][j] = phi_mo[i] * phi_mo[j];

        free(phi_ao);
        free(phi_so);
        free(phi_mo);
    }

    void setup_delta(void)
    {
        static int done=0;
        int i, I, j, errcod;
        int nirreps, nfzc, nfzv;
        int *order, *clsdpi, *openpi, *orbspi, *fruocc, *frdocc;
        double **scf_pitzer;

        if(done) return;

        chkpt_init(PSIO_OPEN_OLD);
        nmo = chkpt_rd_nmo();
        nao = chkpt_rd_nao();
        nirreps = chkpt_rd_nirreps();
        clsdpi = chkpt_rd_clsdpi();
        openpi = chkpt_rd_openpi();
        orbspi = chkpt_rd_orbspi();
        scf_pitzer = chkpt_rd_scf();
        u = chkpt_rd_usotao();
        chkpt_close();

        frdocc = options.get_int_array("FROZEN_DOCC");
        fruocc = options.get_int_array("FROZEN_UOCC");

        nfzc = nfzv = 0;
        for(i=0; i < nirreps; i++) {
            nfzc += frdocc[i];
            nfzv += fruocc[i];
        }

        /*
           if(nfzc || nfzv) {
           printf("Frozen orbitals not yet coded!\n");
           exit(PSI_RETURN_FAILURE);
           }
           */

        /*** Get the Pitzer -> QT reordering array ***/
        order = init_int_array(nmo);
        reorder_qt(clsdpi, openpi, frdocc, fruocc, order, orbspi, nirreps);

        /*** Arrange the SCF eigenvectors into QT ordering ***/
        scf = block_matrix(nmo, nmo);
        for(i=0; i < nmo; i++) {
            I = order[i];  /* Pitzer --> QT */
            for(j=0; j < nmo; j++) scf[j][I] = scf_pitzer[j][i];
        }

        free(order);
        free(clsdpi);
        free(openpi);
        free(orbspi);
        free(fruocc);
        free(frdocc);
        free_block(scf_pitzer);

        done = 1;

        return;
    }




}} // namespace psi::cusp
