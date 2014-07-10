#include <stdio.h>
#include <stdlib.h>

#include "one_electron.h"


static void matrix_block_write (double *matrix, int startrow, int startcol, int ldm,
                                double *block, int nrows, int ncols)
{
    int i;
    int j;
    int k;
    int l;

    for (k = 0; k < nrows; k++)
    {
        for (l = 0; l < ncols; l++)
        {
            i = startrow + k;
            j = startcol + l;
            matrix[i * ldm + j] = block[k + nrows * l];
        }
    }
}


void compute_S (PFock_t pfock, BasisSet_t basis,
                int startshellrow, int endshellrow,
                int startshellcol, int endshellcol,
                double *S)
{
    OED_t oed;
    int A;
    int B;
    int row_id_1;
    int row_id_2;
    int col_id_1;
    int col_id_2;
    int start_row_id;
    int start_col_id;
    int end_col_id;
    int ldS;
    int nrows;
    int ncols;
    int startrow;
    int startcol;
    double *integrals;
    int nints;  
    
    CInt_createOED (basis, &oed);

    start_row_id = pfock->f_startind[startshellrow];
    start_col_id = pfock->f_startind[startshellcol];
    end_col_id = pfock->f_startind[endshellcol + 1] - 1;
    ldS = end_col_id - start_col_id + 1;
    for (A = startshellrow; A <= endshellrow; A++)
    {
        row_id_1 = pfock->f_startind[A];
        row_id_2 = pfock->f_startind[A + 1] - 1;
        startrow = row_id_1 - start_row_id;
        nrows = row_id_2 - row_id_1 + 1;
        for (B = startshellcol; B <= endshellcol; B++)
        {
            col_id_1 = pfock->f_startind[B];
            col_id_2 = pfock->f_startind[B + 1] - 1;
            startcol = col_id_1 - start_col_id;
            ncols = col_id_2 - col_id_1 + 1;
            CInt_computePairOvl (basis, oed, A, B, &integrals, &nints);
            if (nints != 0)
            {
                matrix_block_write (S, startrow, startcol, ldS,
                                    integrals, nrows, ncols);
            }
        }
    }

    CInt_destroyOED (oed);
}


void compute_H (PFock_t pfock, BasisSet_t basis,
                int startshellrow, int endshellrow,
                int startshellcol, int endshellcol,
                double *H)
{
    OED_t oed;
    int A;
    int B;
    int row_id_1;
    int row_id_2;
    int col_id_1;
    int col_id_2;
    int start_row_id;
    int start_col_id;
    int end_col_id;
    int ldH;
    int nrows;
    int ncols;
    int startrow;
    int startcol;
    double *integrals;
    int nints;
    
    CInt_createOED (basis, &oed);    

    start_row_id = pfock->f_startind[startshellrow];
    start_col_id = pfock->f_startind[startshellcol];
    end_col_id = pfock->f_startind[endshellcol + 1] - 1;
    ldH = end_col_id - start_col_id + 1;
    for (A = startshellrow; A <= endshellrow; A++)
    {
        row_id_1 = pfock->f_startind[A];
        row_id_2 = pfock->f_startind[A + 1] - 1;
        startrow = row_id_1 - start_row_id;
        nrows = row_id_2 - row_id_1 + 1;
        for (B = startshellcol; B <= endshellcol; B++)
        {
            col_id_1 = pfock->f_startind[B];
            col_id_2 = pfock->f_startind[B + 1] - 1;
            startcol = col_id_1 - start_col_id;
            ncols = col_id_2 - col_id_1 + 1;
            CInt_computePairCoreH (basis, oed, A, B, &integrals, &nints);
            if (nints != 0)
            {
                matrix_block_write (H, startrow, startcol, ldH,
                                    integrals, nrows, ncols);
            }
        }
    }

    CInt_destroyOED (oed);
}
