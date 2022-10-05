#include "psi4/libmints/matrix.h"

#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

SharedMatrix block_to_matrix(double ** block) {
    auto mat = std::make_shared<Matrix>("Matrix", moinfo.Ca->colspi(), moinfo.Ca->colspi());
    int mo_offset = 0;

    for (int h = 0; h < moinfo.Ca->nirrep(); h++) {
        int nmo = moinfo.orbspi[h];
        int nfv = moinfo.fruocc[h];
        int nmor = nmo - nfv;
        if (!nmo || !nmor) continue;

        // Loop over QT, convert to Pitzer
        double **matp = mat->pointer(h);
        for (int i = 0; i < nmo; i++) { // nmor in original code, density matrix
            for (int j = 0; j < nmo; j++) { // nmor in original code, density matrix
                int I = moinfo.pitzer2qt[i + mo_offset];
                int J = moinfo.pitzer2qt[j + mo_offset];
                matp[i][j] = block[I][J];
            }
        }
        mo_offset += nmo;
    }

    return mat;
} 

}
}
