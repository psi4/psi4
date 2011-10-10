/*!
  \file
  \ingroup CHKPT
*/

#include <cstdio>
#include <cstdlib>
#include <psifiles.h>
#include <boost/shared_ptr.hpp>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.h>
#include <libchkpt/chkpt.hpp>

using namespace psi;

double **Chkpt::rd_scf(void)
{
        double **scf;
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("MO coefficients");

        if (psio->tocscan(PSIF_CHKPT, keyword) != NULL) {

                nmo = rd_nmo();
                nso = rd_nso();

                scf = matrix<double>(nso,nmo);
                psio->read_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                        nso*nmo*sizeof(double));
        }
        else {
                scf = rd_alpha_scf();
        }

        free(keyword);
        return scf;
}

double **Chkpt::rd_alpha_scf(void)
{
        double **scf;
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("Alpha MO coefficients");

        if (psio->tocscan(PSIF_CHKPT, keyword) != NULL) {
                nmo = rd_nmo();
                nso = rd_nso();

                scf = matrix<double>(nso,nmo);
                psio->read_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                        nso*nmo*sizeof(double));
        }
        else
                scf = NULL;

        free(keyword);
        return scf;
}

double **Chkpt::rd_beta_scf(void)
{
        double **scf;
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("Beta MO coefficients");

        if (psio->tocscan(PSIF_CHKPT, keyword) != NULL) {
                nmo = rd_nmo();
                nso = rd_nso();

                scf = matrix<double>(nso,nmo);
                psio->read_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                        nso*nmo*sizeof(double));
        }
        else
                scf = NULL;

        free(keyword);
        return scf;
}

void Chkpt::wt_scf(double **scf)
{
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("MO coefficients");

        nmo = rd_nmo();
        nso = rd_nso();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                nso*nmo*sizeof(double));
        free(keyword);
}

void Chkpt::wt_alpha_scf(double **scf)
{
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("Alpha MO coefficients");

        nmo = rd_nmo();
        nso = rd_nso();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                nso*nmo*sizeof(double));

        free(keyword);
}

void Chkpt::wt_beta_scf(double **scf)
{
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("Beta MO coefficients");

        nmo = rd_nmo();
        nso = rd_nso();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                nso*nmo*sizeof(double));

        free(keyword);
}

double **Chkpt::rd_scf_irrep(int irrep)
{
        int i, j, row, col;
        int *sopi, *mopi;
        double **scf, **scf_full;

        sopi = rd_sopi();
        mopi = rd_orbspi();

        if (!sopi[irrep] || !mopi[irrep]) {
                free(sopi);
                free(mopi);
                return NULL;
    }

        scf_full = rd_scf();
        if (scf_full == NULL) {
                free(sopi);
                free(mopi);
                return NULL;
        }

        scf = matrix<double>(sopi[irrep],mopi[irrep]);

/* compute row and column offsets */
        for(i=0,row=0,col=0; i < irrep; i++) {
                row += sopi[i];
                col += mopi[i];
        }

        for(i=0; i < sopi[irrep]; i++)
                for(j=0; j < mopi[irrep]; j++)
                scf[i][j] = scf_full[i+row][j+col];

        free(scf_full);
        free(sopi);
        free(mopi);

        return scf;
}

double **Chkpt::rd_alpha_scf_irrep(int irrep)
{
        int i, j, row, col;
        int nirreps, nso, nmo;
        int *sopi, *mopi;
        double **scf, **scf_full;

        nirreps = rd_nirreps();
        sopi = rd_sopi();
        mopi = rd_orbspi();
        nso = rd_nso();
        nmo = rd_nmo();

        scf = matrix<double>(sopi[irrep],mopi[irrep]);
        scf_full = rd_alpha_scf();
        if (scf_full == NULL) {
                free(scf);
                free(sopi);
                free(mopi);
                return NULL;
        }

/* compute row and column offsets */
        for(i=0,row=0,col=0; i < irrep; i++) {
                row += sopi[i];
                col += mopi[i];
        }

        for(i=0; i < sopi[irrep]; i++)
                for(j=0; j < mopi[irrep]; j++)
                scf[i][j] = scf_full[i+row][j+col];

        free(scf_full);
        free(sopi);
        free(mopi);

        return scf;
}

double **Chkpt::rd_beta_scf_irrep(int irrep)
{
        int i, j, row, col;
        int nirreps, nso, nmo;
        int *sopi, *mopi;
        double **scf, **scf_full;

        nirreps = rd_nirreps();
        sopi = rd_sopi();
        mopi = rd_orbspi();
        nso = rd_nso();
        nmo = rd_nmo();

        scf = matrix<double>(sopi[irrep],mopi[irrep]);
        scf_full = rd_beta_scf();
        if (scf_full == NULL) {
                free(scf);
                free(sopi);
                free(mopi);
                return NULL;
        }

/* compute row and column offsets */
        for(i=0,row=0,col=0; i < irrep; i++) {
                row += sopi[i];
                col += mopi[i];
        }

        for(i=0; i < sopi[irrep]; i++)
                for(j=0; j < mopi[irrep]; j++)
                scf[i][j] = scf_full[i+row][j+col];

        free(scf_full);
        free(sopi);
        free(mopi);

        return scf;
}

void Chkpt::wt_scf_irrep(double **scf, int irrep)
{
        int i, j, row, col;
        int nirreps, nso, nmo;
        int *sopi, *mopi;
        double **scf_full;

        nirreps = rd_nirreps();
        sopi = rd_sopi();
        mopi = rd_orbspi();
        nso = rd_nso();
        nmo = rd_nmo();

        scf_full = rd_scf();

/* compute row and column offsets */
        for(i=0,row=0,col=0; i < irrep; i++) {
                row += sopi[i];
                col += mopi[i];
        }

        for(i=0; i < sopi[irrep]; i++)
                for(j=0; j < mopi[irrep]; j++)
                scf_full[i+row][j+col] = scf[i][j];

        wt_scf(scf_full);
        free(scf_full);
        free(sopi);
        free(mopi);
}

void Chkpt::wt_alpha_scf_irrep(double **scf, int irrep)
{
        int i, j, row, col;
        int nirreps, nso, nmo;
        int *sopi, *mopi;
        double **scf_full;

        nirreps = rd_nirreps();
        sopi = rd_sopi();
        mopi = rd_orbspi();
        nso = rd_nso();
        nmo = rd_nmo();

        scf_full = rd_alpha_scf();

/* compute row and column offsets */
        for(i=0,row=0,col=0; i < irrep; i++) {
                row += sopi[i];
                col += mopi[i];
        }

        for(i=0; i < sopi[irrep]; i++)
                for(j=0; j < mopi[irrep]; j++)
                scf_full[i+row][j+col] = scf[i][j];

        wt_alpha_scf(scf_full);
        free(scf_full);
        free(sopi);
        free(mopi);
}

void Chkpt::wt_beta_scf_irrep(double **scf, int irrep)
{
        int i, j, row, col;
        int nirreps, nso, nmo;
        int *sopi, *mopi;
        double **scf_full;

        nirreps = rd_nirreps();
        sopi = rd_sopi();
        mopi = rd_orbspi();
        nso = rd_nso();
        nmo = rd_nmo();

        scf_full = rd_beta_scf();

/* compute row and column offsets */
        for(i=0,row=0,col=0; i < irrep; i++) {
                row += sopi[i];
                col += mopi[i];
        }

        for(i=0; i < sopi[irrep]; i++)
                for(j=0; j < mopi[irrep]; j++)
                scf_full[i+row][j+col] = scf[i][j];

        wt_beta_scf(scf_full);
        free(scf_full);
        free(sopi);
        free(mopi);
}

double **Chkpt::set_mo_phases(double **coeff, int nrows, int ncols)
{
        int col, row;

        for (col=0; col<ncols; col++) {
                if (coeff[0][col] >= 0.0) continue;
                for (row=0; row<nrows; row++) {
                        coeff[row][col] = -coeff[row][col];
                }
        }

        wt_phase_check(1);
        return coeff;
}

double **Chkpt::rd_local_scf(void)
{
        double **scf;
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("Local MO coefficients");

        if (psio->tocscan(PSIF_CHKPT, keyword) != NULL) {
                nmo = rd_nmo();
                nso = rd_nso();

                scf = matrix<double>(nso,nmo);
                psio->read_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                        nso*nmo*sizeof(double));
        }
        else
                scf = NULL;

        return scf;
}

void Chkpt::wt_local_scf(double **scf)
{
        int nmo, nso;
        char *keyword;
        keyword = build_keyword("Local MO coefficients");

        nmo = rd_nmo();
        nso = rd_nso();

        psio->write_entry(PSIF_CHKPT, keyword, (char *) scf[0],
                nso*nmo*sizeof(double));

        free(keyword);
}

extern "C" {
/*!
** chkpt_rd_scf():  Reads in the full SCF eigenvector matrix for RHF/ROHF.
**
**   takes no arguments.
**
**   returns: double **scf = This rectangular matrix has dimensions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector would
**     come out looking something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries
**    represent (double)0.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_scf(void)
        {
                return _default_chkpt_lib_->rd_scf();
        }

/*!
** chkpt_rd_alpha_scf(): Reads in the full alpha SCF eigenvector matrix for UHF
**
**   takes no arguments.
**
**   returns: double **scf =  This rectangular matrix has dimensions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector would
**     come out looking something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries
**    represent (double)0.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_alpha_scf(void)
        {
                return _default_chkpt_lib_->rd_alpha_scf();
        }

/*!
** chkpt_rd_beta_scf():  Reads in the full beta SCF eigenvector matrix for UHF.
**
**   takes no arguments.
**
**   returns: double **scf = This rectangular matrix has dimensions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector would
**     come out looking something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries
**    represent (double)0.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_beta_scf(void)
        {
                return _default_chkpt_lib_->rd_beta_scf();
        }

/*!
** chkpt_wt_scf():  Writes the full SCF eigenvector matrix for RHF/ROHF.
**
** \param scf = This rectangular matrix has dimensions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector would
**     come out looking something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries
**    represent (double)0.
**
** returns: none
**
** NOTE: The input scf matrix must occupy a contiguous block of nmo x
** nso memory.  Use matrix<double>() to allocate space for
** the matrix.
**
** \ingroup CHKPT
*/
        void chkpt_wt_scf(double **scf)
        {
                _default_chkpt_lib_->wt_scf(scf);
        }

/*!
** chkpt_wt_alpha_scf():  Writes the full alpha SCF eigenvector matrix for UHF.
**
** \param scf = This rectangular matrix has dimensions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector would
**     come out looking something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries
**    represent (double)0.
**
** returns: none
**
** NOTE: The input scf matrix must occupy a contiguous block of nmo x
** nso memory.  Use matrix<double>() to allocate space for
** the matrix.
**
** \ingroup CHKPT
*/
        void chkpt_wt_alpha_scf(double **scf)
        {
                _default_chkpt_lib_->wt_alpha_scf(scf);
        }

/*!
** chkpt_wt_beta_scf():  Writes the full beta SCF eigenvector matrix for UHF.
**
** \param scf = This rectangular matrix has dimensions nso
**     by nmo (see: rd_nmo()). For STO water, scf_vector would
**     come out looking something like the following:
**
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         *** *** *** *** 0.0 0.0 0.0
**         0.0 0.0 0.0 0.0 *** 0.0 0.0
**         0.0 0.0 0.0 0.0 0.0 *** ***
**         0.0 0.0 0.0 0.0 0.0 *** ***
**
**    where the *** represent the non-zero values, and the 0.0 entries
**    represent (double)0.
**
** returns: none
**
** NOTE: The input scf matrix must occupy a contiguous block of nmo x
** nso memory.  Use matrix<double>() to allocate space for
** the matrix.
**
** \ingroup CHKPT
*/
        void chkpt_wt_beta_scf(double **scf)
        {
                _default_chkpt_lib_->wt_beta_scf(scf);
        }

/*!
** chkpt_rd_scf_irrep(): Reads a single irrep of the SCF eigenvectors for
** RHF/ROHF.
**
** \param irrep = The desired irreducible representation.
**
** returns: double **scf   A rectangualr sopi[irrep] by orbspi[irrep] matrix.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_scf_irrep(int irrep)
        {
                return _default_chkpt_lib_->rd_scf_irrep(irrep);
        }

/*!
** chkpt_rd_alpha_scf_irrep(): Reads a single irrep of the alpha SCF
** eigenvectors for UHF.
**
** \param irrep = The desired irreducible representation.
**
** returns: double **scf = A rectangualr sopi[irrep] by orbspi[irrep] matrix.
** \ingroup CHKPT
*/
        double **chkpt_rd_alpha_scf_irrep(int irrep)
        {
                return _default_chkpt_lib_->rd_alpha_scf_irrep(irrep);
        }

/*!
** chkpt_rd_beta_scf_irrep(): Reads a single irrep of the beta SCF
** eigenvectors for UHF.
**
** \param irrep = The desired irreducible representation.
**
** returns: double **scf = A rectangualr sopi[irrep] by orbspi[irrep] matrix.
**
** \ingroup CHKPT
*/
        double **chkpt_rd_beta_scf_irrep(int irrep)
        {
                return _default_chkpt_lib_->rd_beta_scf_irrep(irrep);
        }

/*!
** chkpt_wt_scf_irrep(): Writes a single irrep of the SCF eigenvectors for
** RHF/ROHF.
**
** \param irrep = The desired irreducible representation.
**
** returns: double **scf = A rectangualr sopi[irrep] by orbspi[irrep] matrix.
**
** \ingroup CHKPT
*/
        void chkpt_wt_scf_irrep(double **scf, int irrep)
        {
                _default_chkpt_lib_->wt_scf_irrep(scf, irrep);
        }

/*!
** chkpt_wt_alpha_scf_irrep(): Writes a single irrep of the alpha SCF
**                             eigenvectors for RHF/ROHF.
**
** \param irrep = The desired irreducible representation.
**
** returns: double **scf = A rectangualr sopi[irrep] by orbspi[irrep] matrix.
**
** \ingroup CHKPT
*/
        void chkpt_wt_alpha_scf_irrep(double **scf, int irrep)
        {
                _default_chkpt_lib_->wt_alpha_scf_irrep(scf, irrep);
        }

/*!
** chkpt_wt_beta_scf_irrep(): Writes a single irrep of the beta SCF
**                             eigenvectors for RHF/ROHF.
**
** \param irrep = The desired irreducible representation.
**
** returns: double **scf = A rectangualr sopi[irrep] by orbspi[irrep] matrix.
**
** \ingroup CHKPT
*/
        void chkpt_wt_beta_scf_irrep(double **scf, int irrep)
        {
                _default_chkpt_lib_->wt_beta_scf_irrep(scf, irrep);
        }

/*!
** chkpt_set_mo_phases(): Set the phase of the MO's according to the standard
** that the first element in every column is always positive.  This may
** help keep the phase consistent for more complicated problems like
** natural orbitals.
**
** \param coeff = MO coefficient matrix
** \param nrows = number of rows in MO coefficient matrix
** \param ncols = number of columns in MO coefficient matrix
**
** Note: since it only looks at the first element in each column, it will
** not work for matrices where that element can be zero by symmetry.  So,
** this function is only helpful when called for an irrep block at a time.
**
** David Sherrill, July 2002
**
** returns: none
*/
        double **chkpt_set_mo_phases(double **coeff, int nrows, int ncols)
        {
                return _default_chkpt_lib_->set_mo_phases(coeff, nrows, ncols);
        }

/*!
** chkpt_rd_local_scf(): Reads the localized mo coefficents
**
** returns: double **scf = A rectangular nso by nmo matrix.
**
** Nick J. Russ, December 2003
*/
        double **chkpt_rd_local_scf(void)
        {
                return _default_chkpt_lib_->rd_local_scf();
        }

/*!
** chkpt_wt_local_scf(): Writes the localized mo coefficents
**
** returns: none
**
** Nick J. Russ, December 2003
*/
        void chkpt_wt_local_scf(double **scf)
        {
                _default_chkpt_lib_->wt_local_scf(scf);
        }
}
