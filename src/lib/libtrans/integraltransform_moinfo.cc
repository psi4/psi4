#include "integraltransform.h"
#include <libchkpt/chkpt.hpp>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <algorithm>
#include "mospace.h"
#include "spaceinfo.h"

namespace psi{ namespace libtrans{

void
IntegralTransform::process_spaces(std::vector<shared_ptr<MOSpace> > spaces)
{
    shared_ptr<PSIO> psio(new PSIO); psiopp_ipv1_config(psio);
    shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

    std::vector<shared_ptr<MOSpace> >::iterator space;

    _nirreps = chkpt->rd_nirreps();
    _nmo     = chkpt->rd_nmo();
    _nso     = chkpt->rd_nso();
    _nao     = chkpt->rd_nao();
    _labels  = chkpt->rd_irr_labs();
    _enuc    = chkpt->rd_enuc();
    _escf    = chkpt->rd_escf();
    _sopi    = chkpt->rd_sopi();
    _mopi    = chkpt->rd_orbspi();
    _clsdpi  = chkpt->rd_clsdpi();
    _openpi  = chkpt->rd_openpi();
    _frzcpi  = chkpt->rd_frzcpi();
    _frzvpi  = chkpt->rd_frzvpi();
    _nTriSo  = _nso * (_nso + 1) / 2;

    _sosym = init_int_array(_nso);
    int count = 0;
    for(int h = 0; h < _nirreps; ++h){
        for(int i = 0; i < _sopi[h]; ++i, ++count){
            _sosym[count] = h;
        }
    }

    _nfzc = _nfzv = 0;
    for(int h = 0; h < _nirreps; ++h){
        if(_frozenOrbitals == VirOnly || _frozenOrbitals == None){
            _frzcpi[h] = 0;
        }
        if(_frozenOrbitals == OccOnly || _frozenOrbitals == None){
            _frzvpi[h] = 0;
        }
        _nfzc += _frzcpi[h];
        _nfzv += _frzvpi[h];
    }

    // Read the eigenvectors from the checkpoint file
    if(_transformationType == Restricted){
        // Set up for a restricted transformation
        _fullCb = _fullCa = chkpt->rd_scf();
        _Ca = new double**[_nirreps];
        for(int h = 0; h < _nirreps; ++h){
            _Ca[h] = chkpt->rd_scf_irrep(h);
        }
        _Cb = _Ca;
    }else{
        // Set up for an unrestricted transformation
        _fullCa = chkpt->rd_alpha_scf();
        _fullCb = chkpt->rd_beta_scf();
        _Ca = new double**[_nirreps];
        _Cb = new double**[_nirreps];
        for(int h = 0; h < _nirreps; ++h){
            _Ca[h] = chkpt->rd_alpha_scf_irrep(h);
            _Cb[h] = chkpt->rd_beta_scf_irrep(h);
        }
    }

    // Start by adding the AO orbital space - this is always needed
    _spacesUsed.push_back(MOSpace::nil);
    _spaceArrays.push_back(_sopi);
    _spaceArrays.push_back(_sosym);
    int spaceNumber = 0;
    _aSpaceNum[MOSpace::nil->label()] = spaceNumber;
    _bSpaceNum[MOSpace::nil->label()] = spaceNumber++;
    
    for(space = spaces.begin(); space != spaces.end(); ++space){
        shared_ptr<MOSpace> moSpace = *space;
        // If this space has already been added, move on
        if(_aSpaceNum.count(moSpace->label())) continue;
        double ***Ca = new double**[_nirreps];
        double ***Cb = Ca;
        int *aOrbsPI = new int[_nirreps];
        int *bOrbsPI = aOrbsPI;
        int *aIndex, *bIndex;
        int *aOrbSym, *bOrbSym;
        if(_transformationType != Restricted){
            Cb      = new double**[_nirreps];
            bOrbsPI = new int[_nirreps];
        }
        if(moSpace->label() == MOSPACE_OCC){
            // This is the occupied space
            int numAOcc = 0, numBOcc = 0, aOccCount = 0, bOccCount = 0;
            for(int h = 0; h < _nirreps; ++h){
                aOrbsPI[h] = _clsdpi[h] + _openpi[h] - _frzcpi[h];
                numAOcc += aOrbsPI[h];
                Ca[h] = block_matrix(_sopi[h], aOrbsPI[h]);
                // Copy over the occupied eigenvectors for the occupied orbitals in this irrep
                if(_sopi[h] * aOrbsPI[h])
                    for(int mu = 0; mu < _sopi[h]; ++mu){
                        for(int i = 0; i < aOrbsPI[h]; ++i){
                            Ca[h][mu][i] = _Ca[h][mu][i + _frzcpi[h]];
                        }
                    }
                if(_transformationType != Restricted){
                    bOrbsPI[h] = _clsdpi[h] - _frzcpi[h];
                    numBOcc += bOrbsPI[h];
                    Cb[h] = block_matrix(_sopi[h], bOrbsPI[h]);
                    // Copy over the occupied eigenvectors for the occupied orbitals in this irrep
                    if(_sopi[h] * bOrbsPI[h])
                        for(int mu = 0; mu < _sopi[h]; ++mu){
                            for(int i = 0; i < bOrbsPI[h]; ++i){
                                Cb[h][mu][i] = _Cb[h][mu][i + _frzcpi[h]];
                            }
                        }
                }
            }
            bOrbSym = aOrbSym = new int[numAOcc];
            bIndex  = aIndex  = new int[numAOcc];
            if(_transformationType != Restricted){
                bOrbSym = new int[numBOcc];
                bIndex  = new int[numBOcc];
            }
            // Build the reindexing arrays
            int aPitzerCount = 0, bPitzerCount = 0, aOrbCount = 0, bOrbCount = 0;
            int pitzerOffset = 0;
            if(_moOrdering == PitzerOrder){
                for(int h = 0; h < _nirreps; ++h){
                    aPitzerCount = bPitzerCount = pitzerOffset + _frzcpi[h];
                    for(int n = 0; n < aOrbsPI[h]; ++n)
                        aIndex[aOrbCount++] = aPitzerCount++;
                    if(_transformationType != Restricted)
                        for(int n = 0; n < bOrbsPI[h]; ++n)
                            bIndex[bOrbCount++] = bPitzerCount++;
                    pitzerOffset += _mopi[h];
                }
            }
            // Compute the orbital symmetries
            for(int h = 0; h < _nirreps; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aOccCount++] = h;
                if(_transformationType != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bOccCount++] = h;
            }
        }else if(moSpace->label() == MOSPACE_ALL){
            // This is the full MO space
            int numActMO = 0;
            for(int h = 0; h < _nirreps; ++h){
                bOrbsPI[h] = aOrbsPI[h] = _mopi[h] - _frzcpi[h] - _frzvpi[h];
                numActMO += aOrbsPI[h];
                Ca[h] = block_matrix(_sopi[h], aOrbsPI[h]);
                // Copy over the eigenvectors for the active orbitals in this irrep
                if(_sopi[h] * aOrbsPI[h])
                    for(int mu = 0; mu < _sopi[h]; ++mu){
                        for(int i = 0; i < aOrbsPI[h]; ++i){
                            Ca[h][mu][i] = _Ca[h][mu][i + _frzcpi[h]];
                        }
                    }
                if(_transformationType != Restricted){
                    Cb[h] = block_matrix(_sopi[h], bOrbsPI[h]);
                    // Copy over the eigenvectors for the active orbitals in this irrep
                    if(_sopi[h] * bOrbsPI[h])
                        for(int mu = 0; mu < _sopi[h]; ++mu){
                            for(int i = 0; i < bOrbsPI[h]; ++i){
                                Cb[h][mu][i] = _Cb[h][mu][i + _frzcpi[h]];
                            }
                        }
                }
            }
            bOrbSym = aOrbSym = new int[numActMO];
            bIndex  = aIndex  = new int[numActMO];
            // Build the reindexing arrays and orbital symmetries
            int actMOCount = 0, pitzerCount = 0, pitzerOffset = 0;
            if(_moOrdering == PitzerOrder){
                for(int h = 0; h < _nirreps; ++h){
                    pitzerCount = pitzerOffset + _frzcpi[h];
                    for(int n = 0; n < aOrbsPI[h]; ++n){
                        aIndex[actMOCount] = pitzerCount++;
                        aOrbSym[actMOCount++] = h;
                    }
                    pitzerOffset += _mopi[h];
                }
            }
        }else if(moSpace->label() == MOSPACE_VIR){
            // This is the virtual space
            int numAVir = 0, numBVir = 0, aVirCount = 0, bVirCount = 0;
            for(int h = 0; h < _nirreps; ++h){
                if(_transformationType == Restricted){
                    bOrbsPI[h] = aOrbsPI[h] = _mopi[h] - _clsdpi[h] - _frzvpi[h];
                    Ca[h] = block_matrix(_sopi[h], aOrbsPI[h]);
                    // Copy over the eigenvectors for the active orbitals in this irrep
                    if(_sopi[h] * aOrbsPI[h])
                        for(int mu = 0; mu < _sopi[h]; ++mu){
                            for(int i = 0; i < aOrbsPI[h]; ++i){
                                Ca[h][mu][i] = _Ca[h][mu][i + _clsdpi[h]];
                            }
                        }
                }else{
                    aOrbsPI[h] = _mopi[h] - _clsdpi[h] - _frzvpi[h] - _openpi[h];
                    bOrbsPI[h] = _mopi[h] - _clsdpi[h] - _frzvpi[h];
                    if(_sopi[h] * aOrbsPI[h]){
                        Ca[h] = block_matrix(_sopi[h], aOrbsPI[h]);
                        // Copy over the eigenvectors for the active orbitals in this irrep
                        if(_sopi[h] * aOrbsPI[h])
                            for(int mu = 0; mu < _sopi[h]; ++mu){
                                for(int i = 0; i < aOrbsPI[h]; ++i){
                                    Ca[h][mu][i] = _Ca[h][mu][i +  + _clsdpi[h] + _openpi[h]];
                                }
                            }
                    }
                    if(_sopi[h] * bOrbsPI[h]){
                        Cb[h] = block_matrix(_sopi[h], bOrbsPI[h]);
                        // Copy over the eigenvectors for the active orbitals in this irrep
                        if(_sopi[h] * bOrbsPI[h])
                            for(int mu = 0; mu < _sopi[h]; ++mu){
                                for(int i = 0; i < bOrbsPI[h]; ++i){
                                    Cb[h][mu][i] = _Cb[h][mu][i + _clsdpi[h]];
                                }
                            }
                    }
                }
                numAVir += aOrbsPI[h];
                numBVir += bOrbsPI[h];
            }
            bOrbSym = aOrbSym = new int[numAVir];
            bIndex  = aIndex  = new int[numAVir];
            if(_transformationType != Restricted){
                bOrbSym = new int[numBVir];
                bIndex  = new int[numBVir];
            }
            // Build the reindexing arrays
            int aPitzerCount = 0, bPitzerCount = 0, aOrbCount = 0, bOrbCount = 0;
            int pitzerOffset = 0;
            if(_moOrdering == PitzerOrder){
                for(int h = 0; h < _nirreps; ++h){
                    if(_transformationType == Restricted){
                        aPitzerCount = bPitzerCount = pitzerOffset + _clsdpi[h];
                    }else{
                        aPitzerCount = pitzerOffset + _clsdpi[h];
                        bPitzerCount = pitzerOffset + _clsdpi[h] + _openpi[h];
                    }
                    for(int n = 0; n < aOrbsPI[h]; ++n)
                        aIndex[aOrbCount++] = aPitzerCount++;
                    if(_transformationType != Restricted)
                        for(int n = 0; n < bOrbsPI[h]; ++n)
                            bIndex[bOrbCount++] = bPitzerCount++;
                    pitzerOffset += _mopi[h];
                }
            }
            // Compute the orbital symmetries
            for(int h = 0; h < _nirreps; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aVirCount++] = h;
                if(_transformationType != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bVirCount++] = h;
            }
        }




        if(_print > 5){
            fprintf(outfile, "The full Alpha SCF Eigenvectors\n");
            print_mat(_fullCa, _nso, _nmo, outfile);
            fprintf(outfile, "The full Beta SCF Eigenvectors\n");
            print_mat(_fullCb, _nso, _nmo, outfile);
            int nAOrbs = 0, nBOrbs = 0;
            fprintf(outfile, "Adding arrays for space %c:-\n",moSpace->label());
            fprintf(outfile, "\talpha orbsPI = ");
            for(int h = 0; h < _nirreps; ++h){
                fprintf(outfile, "%d ", aOrbsPI[h]);
                nAOrbs += aOrbsPI[h];
            }
            fprintf(outfile, "\n\tbeta orbsPI  = ");
            for(int h = 0; h < _nirreps; ++h){
                fprintf(outfile, "%d ", bOrbsPI[h]);
                nBOrbs += bOrbsPI[h];
            }
            fprintf(outfile, "\n\talpha orbSym = ");
            for(int i = 0; i < nAOrbs; ++i) fprintf(outfile, "%d ", aOrbSym[i]);
            fprintf(outfile, "\n\tbeta orbSym  = ");
            for(int i = 0; i < nBOrbs; ++i) fprintf(outfile, "%d ", bOrbSym[i]);
            fprintf(outfile, "\n\talpha Indexing Array = ");
            for(int i = 0; i < nAOrbs; ++i) fprintf(outfile, "%d ", aIndex[i]);
            fprintf(outfile, "\n\tbeta Indexing Array  = ");
            for(int i = 0; i < nBOrbs; ++i) fprintf(outfile, "%d ", bIndex[i]);
            for(int h = 0; h < _nirreps; ++h){
                fprintf(outfile, "\nAlpha orbitals for irrep %d\n", h);
                print_mat(Ca[h], _sopi[h], aOrbsPI[h], outfile);
            }
            for(int h = 0; h < _nirreps; ++h){
                fprintf(outfile, "\nBeta orbitals for irrep %d\n", h);
                print_mat(Cb[h], _sopi[h], bOrbsPI[h], outfile);
            }
            fprintf(outfile, "\n\n");
        }
        _spacesUsed.push_back(moSpace);
        _aSpaceNum[moSpace->label()] = spaceNumber;
        // This is about to be overridden if it's not a restricted transformation
        _bSpaceNum[moSpace->label()] = spaceNumber++;
        _spaceArrays.push_back(aOrbsPI);
        _spaceArrays.push_back(aOrbSym);
        if(_transformationType != Restricted){
            _bSpaceNum[moSpace->label()] = spaceNumber++;
            _spaceArrays.push_back(bOrbsPI);
            _spaceArrays.push_back(bOrbSym);
        }
        _aMOCoefficients[moSpace->label()] = Ca;
        _bMOCoefficients[moSpace->label()] = Cb;
        _aOrbsPI[moSpace->label()] = aOrbsPI;
        _bOrbsPI[moSpace->label()] = bOrbsPI;
        _aIndices[moSpace->label()]  = aIndex;
        _bIndices[moSpace->label()]  = bIndex;
    }
   
    // Test the DPD indexing lookup...
    if(_print > 5){
        std::vector<shared_ptr<MOSpace> >::iterator it1;
        std::vector<shared_ptr<MOSpace> >::iterator it2;
        fprintf(outfile, "DPD IDs: there are %d unique spaces\n", _spaceArrays.size()/2);
        for(it1 = _spacesUsed.begin(); it1 != _spacesUsed.end(); ++it1){
            for(it2 = _spacesUsed.begin(); it2 != _spacesUsed.end(); ++it2){
                shared_ptr<MOSpace> a = *it1;
                shared_ptr<MOSpace> b = *it2;
                fprintf(outfile, "Alpha: Space 1: %c   Space 2: %c   Unpacked ID: %3d  Packed ID: %3d\n",
                    a->label(), b->label(), DPD_ID(a, b, Alpha, 0), DPD_ID(a, b, Alpha, 1));
                fprintf(outfile, "Beta:  Space 1: %c   Space 2: %c   Unpacked ID: %3d  Packed ID: %3d\n",
                    a->label(), b->label(), DPD_ID(a, b, Beta, 0), DPD_ID(a, b, Beta, 1));
            }
        }
    }
}

}} // End namespaces

