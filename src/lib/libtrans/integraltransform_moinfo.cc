#include "integraltransform.h"
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libmints/molecule.h>
#include <libmints/matrix.h>
#include <libmints/view.h>
#include <libmints/wavefunction.h>
#include <libqt/qt.h>
#include <sstream>
#include "mospace.h"

using namespace boost;
using namespace psi;

void IntegralTransform::common_moinfo_initialize()
{
    _nTriSo  = _nso * (_nso + 1) / 2;
    _nTriMo  = _nmo * (_nmo + 1) / 2;
    _sosym   = init_int_array(_nso);
    _zeros   = init_int_array(_nirreps);

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
}

/**
 * Sets up the DPD information for the transformation
 * by querying the MO spaces passed into the constructor
 */
void
IntegralTransform::process_spaces()
{
    std::vector<shared_ptr<MOSpace> >::const_iterator space;

//    for(int h = 0; h < _nirreps; ++h){
//        fprintf(outfile, "docc = %d socc = %d frzcpi = %d frvirt = %d, mopi = %d, sopi = %d\n",
//                _clsdpi[h], _openpi[h], _frzcpi[h], _frzvpi[h], _mopi[h], _sopi[h]);fflush(outfile);
//    }
    // Build the Pitzer -> Qt lookup, if needed
    int *aQT, *bQT;
    if(_moOrdering == QTOrder){
        aQT = init_int_array(_nmo);
        if(_transformationType == Restricted){
            reorder_qt(_clsdpi, _openpi, _frzcpi, _frzvpi, aQT, _mopi, _nirreps);
        }else{
            bQT = init_int_array(_nmo);
            reorder_qt_uhf(_clsdpi, _openpi, _frzcpi, _frzvpi, aQT, bQT, _mopi, _nirreps);
        }
    }

    // Start by adding the AO orbital space - this is always needed
    _spacesUsed.push_back(MOSPACE_NIL);
    _spaceArrays.push_back(_sopi);
    _spaceArrays.push_back(_sosym);
    _aOrbsPI[MOSPACE_NIL] = _zeros;
    _bOrbsPI[MOSPACE_NIL] = _zeros;

    for(space = _uniqueSpaces.begin(); space != _uniqueSpaces.end(); ++space){
        shared_ptr<MOSpace> moSpace = *space;
        int *aOrbsPI = new int[_nirreps];
        int *bOrbsPI = aOrbsPI;
        int *aIndex, *bIndex;
        int *aOrbSym, *bOrbSym;
        if(_transformationType != Restricted){
            bOrbsPI = new int[_nirreps];
        }
        if(moSpace->label() == MOSPACE_OCC){
            // This is the occupied space
            int numAOcc = 0, numBOcc = 0, aOccCount = 0, bOccCount = 0;
            for(int h = 0; h < _nirreps; ++h){
                aOrbsPI[h] = _clsdpi[h] + _openpi[h] - _frzcpi[h];
                numAOcc += aOrbsPI[h];
                if(_transformationType != Restricted){
                    bOrbsPI[h] = _clsdpi[h] - _frzcpi[h];
                    numBOcc += bOrbsPI[h];
                }
            }
            bOrbSym = aOrbSym = new int[numAOcc];
            bIndex  = aIndex  = new int[numAOcc];
            if(_transformationType != Restricted){
                bOrbSym = new int[numBOcc];
                bIndex  = new int[numBOcc];
            }
            // Build the reindexing arrays for Pitzer ordering
            int aPitzerCount = 0, bPitzerCount = 0, aOrbCount = 0, bOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < _nirreps; ++h){
                aPitzerCount = bPitzerCount = pitzerOffset + _frzcpi[h];
                for(int n = 0; n < aOrbsPI[h]; ++n)
                    aIndex[aOrbCount++] = aPitzerCount++;
                if(_transformationType != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)
                        bIndex[bOrbCount++] = bPitzerCount++;
                pitzerOffset += _mopi[h];
            }
            if(_moOrdering == QTOrder){
                for(int n = 0; n < numAOcc; ++n) aIndex[n] = aQT[aIndex[n]];
                if(_transformationType != Restricted)
                    for(int n = 0; n < numBOcc; ++n) bIndex[n] = bQT[bIndex[n]];
            };
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
            }
            bOrbSym = aOrbSym = new int[numActMO];
            bIndex  = aIndex  = new int[numActMO];
            // Build the reindexing arrays and orbital symmetries
            int actMOCount = 0, pitzerCount = 0, pitzerOffset = 0;
            for(int h = 0; h < _nirreps; ++h){
                pitzerCount = pitzerOffset + _frzcpi[h];
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[actMOCount] = pitzerCount++;
                    aOrbSym[actMOCount++] = h;
                }
                pitzerOffset += _mopi[h];
            }
            if(_moOrdering == QTOrder)
                for(int n = 0; n < numActMO; ++n) aIndex[n] = aQT[aIndex[n]];
        }else if(moSpace->label() == MOSPACE_VIR){
            // This is the virtual space
            int numAVir = 0, numBVir = 0, aVirCount = 0, bVirCount = 0;
            for(int h = 0; h < _nirreps; ++h){
                if(_transformationType == Restricted){
                    bOrbsPI[h] = aOrbsPI[h] = _mopi[h] - _clsdpi[h] - _frzvpi[h];
                }else{
                    aOrbsPI[h] = _mopi[h] - _clsdpi[h] - _frzvpi[h] - _openpi[h];
                    bOrbsPI[h] = _mopi[h] - _clsdpi[h] - _frzvpi[h];
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
            if(_moOrdering == QTOrder){
                for(int n = 0; n < numAVir; ++n) aIndex[n] = aQT[aIndex[n]];
                if(_transformationType != Restricted)
                    for(int n = 0; n < numBVir; ++n) bIndex[n] = bQT[bIndex[n]];
            };

            // Compute the orbital symmetries
            for(int h = 0; h < _nirreps; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aVirCount++] = h;
                if(_transformationType != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bVirCount++] = h;
            }
        }else{
            // This must be a custom MOSpace that the user provided
            aOrbsPI = const_cast<int*>(moSpace->aOrbsPI());
            bOrbsPI = const_cast<int*>(moSpace->bOrbsPI());
            aOrbSym = const_cast<int*>(moSpace->aOrbSym());
            bOrbSym = const_cast<int*>(moSpace->bOrbSym());
            aIndex  = const_cast<int*>(moSpace->aIndex());
            bIndex  = const_cast<int*>(moSpace->bIndex());
            if(_useIWL && aIndex == NULL){
                std::string error("You must provide an indexing array for space ");
                error += moSpace->label();
                error += " or disable IWL output by changing OutputType.";
                throw SanityCheckError(error, __FILE__, __LINE__);
            }
        }

        if(_print > 5){
            int nAOrbs = 0, nBOrbs = 0;
            fprintf(outfile, "Adding arrays for space %c:-\n",moSpace->label());
            fprintf(outfile, "\n\talpha orbsPI = ");
            for(int h = 0; h < _nirreps; nAOrbs += aOrbsPI[h], ++h)
                fprintf(outfile, "%d ", aOrbsPI[h]);
            fprintf(outfile, "\n\tbeta orbsPI = ");
            for(int h = 0; h < _nirreps; nBOrbs += bOrbsPI[h], ++h)
                fprintf(outfile, "%d ", bOrbsPI[h]);
            fprintf(outfile, "\n\talpha orbSym = ");
            for(int i = 0; i < nAOrbs; ++i) fprintf(outfile, "%d ", aOrbSym[i]);
            fprintf(outfile, "\n\tbeta orbSym  = ");
            for(int i = 0; i < nBOrbs; ++i) fprintf(outfile, "%d ", bOrbSym[i]);
            fprintf(outfile, "\n\talpha Indexing Array = ");
            for(int i = 0; i < nAOrbs; ++i) fprintf(outfile, "%d ", aIndex[i]);
            fprintf(outfile, "\n\tbeta Indexing Array  = ");
            for(int i = 0; i < nBOrbs; ++i) fprintf(outfile, "%d ", bIndex[i]);
            fprintf(outfile, "\n\n");
        }

        _spacesUsed.push_back(toupper(moSpace->label()));
        _spaceArrays.push_back(aOrbsPI);
        _spaceArrays.push_back(aOrbSym);
        if(_transformationType != Restricted){
            _spacesUsed.push_back(tolower(moSpace->label()));
            _spaceArrays.push_back(bOrbsPI);
            _spaceArrays.push_back(bOrbSym);
        }
        _aOrbsPI[moSpace->label()]  = aOrbsPI;
        _bOrbsPI[moSpace->label()]  = bOrbsPI;
        _aIndices[moSpace->label()] = aIndex;
        _bIndices[moSpace->label()] = bIndex;
    }// End loop over spaces

    /* Populate the DPD indexing map.  The string class is used instead of a char*
     * because I can't be bothered to roll my own char* container with comparison
     * operators, which are needed for maps.  The memory overhead of this approach
     * is negligible and it's leak-free*/
    int pairCount = 0;
    for(int a = 0; a < _spacesUsed.size(); ++a){
        std::stringstream stream;
        char s = _spacesUsed[a];
        stream << "[" << s << "," << s << "]";
        _dpdLookup[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">" << s << "]+";
        _dpdLookup[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">" << s << "]-";
        _dpdLookup[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">=" << s << "]+";
        _dpdLookup[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">=" << s << "]-";
        _dpdLookup[stream.str()] = pairCount++;
    }

    for(int a = 0; a < _spacesUsed.size(); ++a){
        for(int b = a+1; b < _spacesUsed.size(); ++b){
            std::stringstream stream;
            stream << "[" << _spacesUsed[a] << "," << _spacesUsed[b] << "]";
            _dpdLookup[stream.str()] = pairCount++;
            stream.str("");
            stream << "[" << _spacesUsed[b] << "," << _spacesUsed[a] << "]";
            _dpdLookup[stream.str()] = pairCount++;
        }
    }

    if(_print > 5) print_dpd_lookup();

    if(_moOrdering == QTOrder){
        free(aQT);
        if(_transformationType != Restricted) free(bQT);
    }
}


/**
 * Re-reads the MO coeffiecients from the checkpoint file and sets up
 * the computation.  For default spaces, it is sufficient to call this.  For
 * custom MOSpaces, the MOSpace objects' MO coefficients should be updated first.
 * All one-electron integrals are generate in the new basis, including the Fock matrix.
 */
void
IntegralTransform::update_orbitals()
{
  if(_transformationType == SemiCanonical){
      throw FeatureNotImplemented("Libtrans", " update of semicanonical orbitals",
              __FILE__, __LINE__);
  }
  process_eigenvectors();
  generate_oei();
}


/**
 * Sets up the eigenvectors for the transformation by querying the MO spaces
 * passed into the constructor.  This is done seperately to the DPD setup
 * to give us a chance to semicanonicalize the orbitals if necessary.
 */
void
IntegralTransform::process_eigenvectors()
{
    std::vector<shared_ptr<MOSpace> >::const_iterator space;

    if(_print > 4){
        _Ca->print();
        _Cb->print();
    }

    // N.B. The frozen orbitals have been zeroed, if appropriate
    Dimension focc = _frzcpi;
    Dimension aocc = _clsdpi + _openpi - _frzcpi;
    Dimension bocc = _clsdpi - _frzcpi;
    Dimension avir = _mopi - _clsdpi - _openpi - _frzvpi;
    Dimension bvir = _mopi - _clsdpi - _frzvpi;
    Dimension aall = _mopi - _frzcpi - _frzvpi;
    Dimension ball = _mopi - _frzcpi - _frzvpi;
    Dimension zero = Dimension(_nirreps);

    for(space = _uniqueSpaces.begin(); space != _uniqueSpaces.end(); ++space){
        shared_ptr<MOSpace> moSpace = *space;
        SharedMatrix Ca;
        SharedMatrix Cb = Ca;
        if(moSpace->label() == MOSPACE_OCC){
            // This is the occupied space
            View Vaocc(_Ca, _sopi, aocc, zero, focc);
            Ca = Vaocc();
            Ca->set_name("Alpha occupied orbitals");
            if(_transformationType != Restricted){
                View Vbocc(_Cb, _sopi, bocc, zero, focc);
                Cb = Vbocc();
                Cb->set_name("Beta occupied orbitals");
            }
        }else if(moSpace->label() == MOSPACE_ALL){
            // This is the full space, sans frozen orbitals
            View Vaall(_Ca, _sopi, aall, zero, focc);
            Ca = Vaall();
            Ca->set_name("All alpha orbitals");
            if(_transformationType != Restricted){
                View Vball(_Cb, _sopi, ball, zero, focc);
                Cb = Vball();
                Cb->set_name("All beta orbitals");
            }
        }else if(moSpace->label() == MOSPACE_VIR){
            // This is the virtual space
            if(_transformationType == Restricted){
                View Vavir(_Ca, _sopi, avir, zero, _clsdpi);
                Ca = Vavir();
                Ca->set_name("Alpha virtual orbitals");
            }else{
                View Vavir(_Ca, _sopi, avir, zero, _clsdpi + _openpi);
                Ca = Vavir();
                Ca->set_name("Alpha virtual orbitals");
                View Vbvir(_Cb, _sopi, bvir, zero, _clsdpi);
                Cb = Vbvir();
                Cb->set_name("Beta virtual orbitals");
            }
        }else if(moSpace->label() == MOSPACE_NIL){
            // Do nothing!
        }else{
            //TODO This is a custom space work on this!
        }

        _aMOCoefficients[moSpace->label()] = Ca;
        _bMOCoefficients[moSpace->label()] = Cb;

        if(_print > 5){
            fprintf(outfile, "Orbitals for space %c:-\n",moSpace->label());
            Ca->print();
            Cb->print();
            fprintf(outfile, "\n\n");
        }
    }// End loop over spaces
}


/**
 * Prints the dpd index for each pair type used in this object
 */
void
IntegralTransform::print_dpd_lookup()
{
    std::map<std::string, int>::iterator iter;
    fprintf(outfile, "The DPD mappings used in this transformation:-\n");
    for(iter = _dpdLookup.begin(); iter != _dpdLookup.end(); ++iter){
        fprintf(outfile, "Pair %-10s ID = %d\n", iter->first.c_str(), iter->second);
    }
}
