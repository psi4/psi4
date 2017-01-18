/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "integraltransform.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/view.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libqt/qt.h"
#include <sstream>
#include "mospace.h"

using namespace psi;

void IntegralTransform::common_initialize()
{
    aaIntName_ = "";
    abIntName_ = "";
    bbIntName_ = "";

    keepHtInts_ = 0;

    nTriSo_  = nso_ * (nso_ + 1) / 2;
    nTriMo_  = nmo_ * (nmo_ + 1) / 2;
    sosym_   = init_int_array(nso_);
    mosym_   = init_int_array(nmo_);
    zeros_   = init_int_array(nirreps_);

    write_dpd_so_tpdm_ = false;

    int count = 0;
    for(int h = 0; h < nirreps_; ++h){
        for(int i = 0; i < sopi_[h]; ++i, ++count){
            sosym_[count] = h;
        }
    }

    count = 0;
    for(int h = 0; h < nirreps_; ++h){
        for(int i = 0; i < mopi_[h]; ++i, ++count){
            mosym_[count] = h;
        }
    }

    nfzc_ = nfzv_ = 0;
    for(int h = 0; h < nirreps_; ++h){
        if(frozenOrbitals_ == VirOnly || frozenOrbitals_ == None){
            frzcpi_[h] = 0;
        }
        if(frozenOrbitals_ == OccOnly || frozenOrbitals_ == None){
            frzvpi_[h] = 0;
        }
        nfzc_ += frzcpi_[h];
        nfzv_ += frzvpi_[h];
    }
}

/**
 * Sets up the DPD information for the transformation
 * by querying the MO spaces passed into the constructor
 */
void
IntegralTransform::process_spaces()
{
    std::vector<std::shared_ptr<MOSpace> >::const_iterator space;

    //    for(int h = 0; h < _nirreps; ++h){
    //        outfile->Printf( "docc = %d socc = %d frzcpi = %d frvirt = %d, mopi = %d, sopi = %d\n",
    //                _clsdpi[h], _openpi[h], _frzcpi[h], _frzvpi[h], _mopi[h], _sopi[h]);
    //    }

    bool qt_order = (moOrdering_ == QTOrder);  // If false, we assume Pitzer below

    for(space = uniqueSpaces_.begin(); space != uniqueSpaces_.end(); ++space){

        std::shared_ptr<MOSpace> moSpace = *space;
        int *aOrbsPI = new int[nirreps_];
        int *aIndex;
        int *aOrbSym;

        char label = moSpace->label();
        if(labelsUsed_.count(label)){
            std::string error("Space ");
            error += label;
            error += " is already in use.  Choose a unique name for the custom MOSpace.";
            throw SanityCheckError(error, __FILE__, __LINE__);
        }
        ++labelsUsed_[label];

        if(moSpace->label() == MOSPACE_FZC){
            // This is the frozen occupied space
            int numAOcc = 0, aOccCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                aOrbsPI[h] = frzcpi_[h];
                numAOcc += aOrbsPI[h];
            }
            aOrbSym = new int[numAOcc];
            aIndex  = new int[numAOcc];
            // Build the reindexing arrays for Pitzer ordering
            int aPitzerCount = 0, aOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                aPitzerCount = pitzerOffset;
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aOccCount++] = h;
            }
        }else if(moSpace->label() == MOSPACE_OCC){
            // This is the occupied space
            int numAOcc = 0, aOccCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                aOrbsPI[h] = nalphapi_[h] - frzcpi_[h];
                numAOcc += aOrbsPI[h];
            }
            aOrbSym = new int[numAOcc];
            aIndex  = new int[numAOcc];
            // Build the reindexing arrays for Pitzer ordering
            int aPitzerCount = 0, aOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                aPitzerCount = pitzerOffset + frzcpi_[h];
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aOccCount++] = h;
            }
        }else if(moSpace->label() == MOSPACE_ALL){
            // This is the full MO space
            int numActMO = 0;
            for(int h = 0; h < nirreps_; ++h){
                aOrbsPI[h] = mopi_[h] - frzcpi_[h] - frzvpi_[h];
                numActMO += aOrbsPI[h];
            }
            aOrbSym = new int[numActMO];
            aIndex  = new int[numActMO];
            // Build the reindexing arrays and orbital symmetries
            int aOrbCount = 0, pitzerCount = 0, pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                pitzerCount = pitzerOffset + frzcpi_[h];
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount] = (qt_order ? aQT_[pitzerCount] : pitzerCount);
                    aOrbSym[aOrbCount] = h;
                    aOrbCount++;
                    pitzerCount++;
                }
                pitzerOffset += mopi_[h];
            }
        }else if(moSpace->label() == MOSPACE_VIR){
            // This is the virtual space
            int numAVir = 0, aVirCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                if(transformationType_ == Restricted){
                    aOrbsPI[h] = mopi_[h] - nbetapi_[h] - frzvpi_[h]; // TDC -- keep soccs
                }else{
                    aOrbsPI[h] = mopi_[h] - nalphapi_[h] - frzvpi_[h];
                }
                numAVir += aOrbsPI[h];
            }
            aOrbSym = new int[numAVir];
            aIndex  = new int[numAVir];
            // Build the reindexing arrays
            int aPitzerCount = 0, aOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    if(transformationType_ == Restricted){
                        aPitzerCount = pitzerOffset + nalphapi_[h];
                    }else{
                        aPitzerCount = pitzerOffset + nalphapi_[h];
                    }
                }
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aVirCount++] = h;
            }
        }else if(moSpace->label() == MOSPACE_FZV){
            // This is the frozen virtual space
            int numAVir = 0, aVirCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                aOrbsPI[h] = frzvpi_[h];
                numAVir += aOrbsPI[h];
            }
            aOrbSym = new int[numAVir];
            aIndex  = new int[numAVir];
            // Build the reindexing arrays
            int aPitzerCount = 0, aOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                aPitzerCount = pitzerOffset + mopi_[h] - frzvpi_[h];
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h)
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aVirCount++] = h;
        }else if(moSpace->label() == MOSPACE_DUM){
            // This is the dummy single-function-per-irrep space
            aOrbSym = new int[1];
            aIndex  = new int[1];
            for(int h = 0; h < nirreps_; ++h)
                aOrbsPI[h] = 0;
            aOrbsPI[0] = 1;
            aOrbSym[0] = 0;
            aIndex[0] = 0;
        }else{
            // This must be a custom MOSpace that the user provided
            if(moSpace->placeholder()){
                // This is an AO space, such as an auxilliary basis set
                const std::vector<int> &aorbspi  = moSpace->aOrbs();
                int numOrbs = 0;
                for(int h = 0; h < nirreps_; ++h){
                    aOrbsPI[h] = aorbspi[h];
                    numOrbs += aorbspi[h];
                }
                aOrbSym = new int[numOrbs];
                int orbCount = 0;
                for(int h = 0; h < nirreps_; ++h)
                    for(int orb = 0; orb < aOrbsPI[h]; ++orb)
                        aOrbSym[orbCount++] = h;
            }else{
                const std::vector<int> aorbs  = moSpace->aOrbs();
                const std::vector<int> aindex = moSpace->aIndex();

                // Figure out how many orbitals per irrep, and group all orbitals by irrep
                int nAOrbs = aorbs.size();
                aOrbSym = new int[nAOrbs];
                ::memset(aOrbsPI, '\0', nirreps_*sizeof(int));
                aIndex = (aindex.empty() ? 0 : new int[nAOrbs]);
                for(int h = 0, count = 0; h < nirreps_; ++h){
                    for(int n = 0; n < nAOrbs; ++n){
                        int orb = aorbs[n];
                        if(mosym_[orb] == h){
                            aOrbsPI[h]++;
                            aOrbSym[count] = h;
                            if(aIndex) aIndex[count] = aindex[n];
                            count++;
                        }
                    }
                }
                // Check that the indexing array was provided, if needed
                if(useIWL_ && aindex.empty()){
                    std::string error("You must provide an indexing array for space ");
                    error += moSpace->label();
                    error += " or disable IWL output by changing OutputType.";
                    throw SanityCheckError(error, __FILE__, __LINE__);
                }
            }
        }

        if(print_ > 5){
            int nAOrbs = 0;
            outfile->Printf( "Adding arrays for space %c:-\n",moSpace->label());
            outfile->Printf( "\n\talpha orbsPI = ");
            for(int h = 0; h < nirreps_; nAOrbs += aOrbsPI[h], ++h)
                outfile->Printf( "%d ", aOrbsPI[h]);
            outfile->Printf( "\n\talpha orbSym = ");
            for(int i = 0; i < nAOrbs; ++i) outfile->Printf( "%d ", aOrbSym[i]);
            outfile->Printf( "\n\talpha Indexing Array = ");
            for(int i = 0; i < nAOrbs; ++i) outfile->Printf( "%d ", aIndex[i]);
            outfile->Printf( "\n\n");
        }

        spacesUsed_.push_back(toupper(moSpace->label()));
        spaceArray_.push_back(aOrbsPI);
        spaceArray_.push_back(aOrbSym);
        aOrbsPI_[moSpace->label()]  = aOrbsPI;
        aIndices_[moSpace->label()] = aIndex;
//        if(transformationType_ == Restricted){
//            spacesUsed_.push_back(tolower(moSpace->label()));
//            spaceArray_.push_back(bOrbsPI);
//            spaceArray_.push_back(bOrbSym);
//            bOrbsPI_[moSpace->label()]  = bOrbsPI;
//            bIndices_[moSpace->label()] = bIndex;
//        }

    } // End loop over alpha spaces

    // And now the beta spaces, if needed
    if(transformationType_ != Restricted){
        for(space = uniqueSpaces_.begin(); space != uniqueSpaces_.end(); ++space){
            std::shared_ptr<MOSpace> moSpace = *space;
            int *bOrbsPI = new int[nirreps_];;
            int *bIndex;
            int *bOrbSym;

            if(moSpace->label() == MOSPACE_FZC){
                // This is the frozen occupied space
                int numBOcc = 0, bOccCount = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bOrbsPI[h] = frzcpi_[h];
                    numBOcc += bOrbsPI[h];
                }
                bOrbSym = new int[numBOcc];
                bIndex  = new int[numBOcc];
                // Build the reindexing arrays for Pitzer ordering
                int bPitzerCount = 0, bOrbCount = 0;
                int pitzerOffset = 0;
                for(int h = 0; h < nirreps_; ++h){
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[bPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                    pitzerOffset += mopi_[h];
                }
                // Compute the orbital symmetries
                for(int h = 0; h < nirreps_; ++h)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bOccCount++] = h;
            }else if(moSpace->label() == MOSPACE_OCC){
                // This is the occupied space
                int numBOcc = 0, bOccCount = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bOrbsPI[h] = nbetapi_[h] - frzcpi_[h];
                    numBOcc += bOrbsPI[h];
                }
                bOrbSym = new int[numBOcc];
                bIndex  = new int[numBOcc];
                // Build the reindexing arrays for Pitzer ordering
                int bPitzerCount = 0, bOrbCount = 0;
                int pitzerOffset = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bPitzerCount = pitzerOffset + frzcpi_[h];
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[bPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                    pitzerOffset += mopi_[h];
                }
                // Compute the orbital symmetries
                for(int h = 0; h < nirreps_; ++h)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bOccCount++] = h;
            }else if(moSpace->label() == MOSPACE_ALL){
                // This is the full MO space
                int numActMO = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bOrbsPI[h] = mopi_[h] - frzcpi_[h] - frzvpi_[h];
                    numActMO += bOrbsPI[h];
                }
                bOrbSym = new int[numActMO];
                bIndex  = new int[numActMO];
                // Build the reindexing arrays and orbital symmetries
                int bOrbCount = 0, pitzerCount = 0, pitzerOffset = 0;
                for(int h = 0; h < nirreps_; ++h){
                    pitzerCount = pitzerOffset + frzcpi_[h];
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount] = (qt_order ? bQT_[pitzerCount] : pitzerCount);
                        bOrbSym[bOrbCount++] = h;
                        pitzerCount++;
                    }
                    pitzerOffset += mopi_[h];
                }
            }else if(moSpace->label() == MOSPACE_VIR){
                // This is the virtual space
                int numBVir = 0,bVirCount = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bOrbsPI[h] = mopi_[h] - nbetapi_[h] - frzvpi_[h];
                    numBVir += bOrbsPI[h];
                }
                bOrbSym = new int[numBVir];
                bIndex  = new int[numBVir];
                // Build the reindexing arrays
                int bPitzerCount = 0, bOrbCount = 0;
                int pitzerOffset = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bPitzerCount = pitzerOffset + nbetapi_[h];
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[bPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                    pitzerOffset += mopi_[h];
                }
                // Compute the orbital symmetries
                for(int h = 0; h < nirreps_; ++h)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bVirCount++] = h;
            }else if(moSpace->label() == MOSPACE_FZV){
                // This is the frozen virtual space
                int numBVir = 0, bVirCount = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bOrbsPI[h] = frzvpi_[h];
                    numBVir += bOrbsPI[h];
                }
                bOrbSym = new int[numBVir];
                bIndex  = new int[numBVir];
                // Build the reindexing arrays
                int bPitzerCount = 0, bOrbCount = 0;
                int pitzerOffset = 0;
                for(int h = 0; h < nirreps_; ++h){
                    bPitzerCount = pitzerOffset + mopi_[h] - frzvpi_[h];
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[bPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                    pitzerOffset += mopi_[h];
                }
                // Compute the orbital symmetries
                for(int h = 0; h < nirreps_; ++h)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bVirCount++] = h;
            }else if(moSpace->label() == MOSPACE_DUM){
                // This is the dummy single-function space
                bOrbSym = new int[1];
                bIndex  = new int[1];
                for(int h = 0; h < nirreps_; ++h)
                    bOrbsPI[h] = 0;
                bOrbsPI[0] = 1;
                bOrbSym[0] = 0;
                bIndex[0] = 0;
            }else{
                if(moSpace->placeholder()){
                    // This is an AO space, such as an auxilliary basis set
                    const std::vector<int> borbspi  = moSpace->aOrbs();
                    int numOrbs = 0;
                    for(int h = 0; h < nirreps_; ++h){
                        bOrbsPI[h] = borbspi[h];
                        numOrbs += borbspi[h];
                    }
                    bOrbSym = new int[numOrbs];
                    int orbCount = 0;
                    for(int h = 0; h < nirreps_; ++h)
                        for(int orb = 0; orb < bOrbsPI[h]; ++orb)
                            bOrbSym[orbCount++] = h;

                }else{
                    // This must be a custom MOSpace that the user provided
                    const std::vector<int> borbs  = moSpace->bOrbs();
                    const std::vector<int> bindex = moSpace->bIndex();

                    // Figure out how many orbitals per irrep, and group all orbitals by irrep
                    int nBOrbs = borbs.size();
                    bOrbSym = new int[nBOrbs];
                    bIndex = (bindex.empty() ? 0 : new int[nBOrbs]);
                    ::memset(bOrbsPI, '\0', nirreps_*sizeof(int));
                    for(int h = 0, count = 0; h < nirreps_; ++h){
                        for(int n = 0; n < nBOrbs; ++n){
                            int orb = borbs[n];
                            if(mosym_[orb] == h){
                                bOrbsPI[h]++;
                                bOrbSym[count] = h;
                                if(bIndex) bIndex[count] = bindex[n];
                                count++;
                            }
                        }
                    }
                    // Check that the indexing array was provided, if needed
                    if(useIWL_ && bindex.empty()){
                        std::string error("You must provide a beta indexing array for space ");
                        error += moSpace->label();
                        error += " or disable IWL output by changing OutputType.";
                        throw SanityCheckError(error, __FILE__, __LINE__);
                    }
                }
            }
            if(print_ > 5){
                int nAOrbs = 0, nBOrbs = 0;
                outfile->Printf( "Adding arrays for space %c:-\n",moSpace->label());
                outfile->Printf( "\n\tbeta orbsPI = ");
                for(int h = 0; h < nirreps_; nBOrbs += bOrbsPI[h], ++h)
                    outfile->Printf( "%d ", bOrbsPI[h]);
                outfile->Printf( "\n\tbeta orbSym  = ");
                for(int i = 0; i < nBOrbs; ++i) outfile->Printf( "%d ", bOrbSym[i]);
                outfile->Printf( "\n\tbeta Indexing Array  = ");
                for(int i = 0; i < nBOrbs; ++i) outfile->Printf( "%d ", bIndex[i]);
                outfile->Printf( "\n\n");
            }

            spacesUsed_.push_back(tolower(moSpace->label()));
            spaceArray_.push_back(bOrbsPI);
            spaceArray_.push_back(bOrbSym);
            bOrbsPI_[moSpace->label()]  = bOrbsPI;
            bIndices_[moSpace->label()] = bIndex;
        } // End loop over spaces
    }

    // End by adding the AO orbital space - this is always needed
    spacesUsed_.push_back(MOSPACE_NIL);
    spaceArray_.push_back(sopi_);
    spaceArray_.push_back(sosym_);
    aOrbsPI_[MOSPACE_NIL] = sopi_;
    bOrbsPI_[MOSPACE_NIL] = sopi_;

    /* Populate the DPD indexing map.  The string class is used instead of a char*
     * because I can't be bothered to roll my own char* container with comparison
     * operators, which are needed for maps.  The memory overhead of this approach
     * is negligible and it's leak-free*/
    int pairCount = 0;
    for(int a = 0; a < spacesUsed_.size(); ++a){
        std::stringstream stream;
        char s = spacesUsed_[a];
        stream << "[" << s << "," << s << "]";
        dpdLookup_[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">" << s << "]+";
        dpdLookup_[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">" << s << "]-";
        dpdLookup_[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">=" << s << "]+";
        dpdLookup_[stream.str()] = pairCount++;
        stream.str("");
        stream << "[" << s << ">=" << s << "]-";
        dpdLookup_[stream.str()] = pairCount++;
    }

    for(int a = 0; a < spacesUsed_.size(); ++a){
        for(int b = a+1; b < spacesUsed_.size(); ++b){
            std::stringstream stream;
            stream << "[" << spacesUsed_[a] << "," << spacesUsed_[b] << "]";
            dpdLookup_[stream.str()] = pairCount++;
            stream.str("");
            stream << "[" << spacesUsed_[b] << "," << spacesUsed_[a] << "]";
            dpdLookup_[stream.str()] = pairCount++;
        }
    }

    if(print_ > 5) print_dpd_lookup();
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
    if(transformationType_ == SemiCanonical){
        throw FeatureNotImplemented("Libtrans", " update of semicanonical orbitals",
                                    __FILE__, __LINE__);
    }
    process_eigenvectors();
    generate_oei();
}
/**
 * Sets the orbital matrix, but touches nothing else. This is used for a MCSCF wavefunction
 * and is a bit of a hack, use at your own risk.
**/
void
IntegralTransform::set_orbitals(SharedMatrix C)
{
    Ca_ = C->clone();
    Cb_ = Ca_;
    process_eigenvectors();
}

/**
 * Sets up the eigenvectors for the transformation by querying the MO spaces
 * passed into the constructor.  This is done seperately to the DPD setup
 * to give us a chance to semicanonicalize the orbitals if necessary.
 */
void
IntegralTransform::process_eigenvectors()
{
    std::vector<std::shared_ptr<MOSpace> >::const_iterator space;

    if(print_ > 4){
        Ca_->print();
        Cb_->print();
    }

    // N.B. The frozen orbitals have been zeroed, if appropriate
    Dimension focc = frzcpi_;
    Dimension aocc = nalphapi_ - frzcpi_;
    Dimension bocc = nbetapi_ - frzcpi_;
    Dimension avir = mopi_ - nalphapi_ - frzvpi_;
    Dimension bvir = mopi_ - nbetapi_ - frzvpi_;
    Dimension aall = mopi_ - frzcpi_ - frzvpi_;
    Dimension fvir = frzvpi_;
    Dimension ball = mopi_ - frzcpi_ - frzvpi_;
    Dimension zero = Dimension(nirreps_);

    for(space = uniqueSpaces_.begin(); space != uniqueSpaces_.end(); ++space){
        std::shared_ptr<MOSpace> moSpace = *space;
        SharedMatrix Ca, Cb;
        if(moSpace->label() == MOSPACE_FZC){
            // This is the frozen occupied space
            View Vafzc(Ca_, sopi_, focc);
            Ca = Vafzc();
            Ca->set_name("Alpha frozen occupied orbitals");
            if(transformationType_ != Restricted){
                View Vbfzc(Cb_, sopi_, focc);
                Cb = Vbfzc();
                Cb->set_name("Beta frozen occupied orbitals");
            }
        }else if(moSpace->label() == MOSPACE_OCC){
            // This is the occupied space
            View Vaocc(Ca_, sopi_, aocc, zero, focc);
            Ca = Vaocc();
            Ca->set_name("Alpha occupied orbitals");
            if(transformationType_ != Restricted){
                View Vbocc(Cb_, sopi_, bocc, zero, focc);
                Cb = Vbocc();
                Cb->set_name("Beta occupied orbitals");
            }
        }else if(moSpace->label() == MOSPACE_ALL){
            // This is the full space, sans frozen orbitals
            View Vaall(Ca_, sopi_, aall, zero, focc);
            Ca = Vaall();
            Ca->set_name("All alpha orbitals");
            if(transformationType_ != Restricted){
                View Vball(Cb_, sopi_, ball, zero, focc);
                Cb = Vball();
                Cb->set_name("All beta orbitals");
            }
        }else if(moSpace->label() == MOSPACE_VIR){
            // This is the virtual space
            if(transformationType_ == Restricted){
                // Take the true virtual orbitals, and then append the SOCC orbitals
                View Vavir(Ca_, sopi_, avir, zero, nalphapi_);
                View Vasoc(Ca_, sopi_, openpi_, zero, clsdpi_);
                std::vector<SharedMatrix> virandsoc;
                virandsoc.push_back(Vavir());
                virandsoc.push_back(Vasoc());
                Ca = Matrix::horzcat(virandsoc);
                Ca->set_name("Alpha virtual orbitals");
            }else{
                View Vavir(Ca_, sopi_, avir, zero, nalphapi_);
                Ca = Vavir();
                Ca->set_name("Alpha virtual orbitals");
                View Vbvir(Cb_, sopi_, bvir, zero, nbetapi_);
                Cb = Vbvir();
                Cb->set_name("Beta virtual orbitals");
            }
        }else if(moSpace->label() == MOSPACE_FZV){
            // This is the frozen virtual space
            View Vafzv(Ca_, sopi_, fvir, zero, mopi_ - frzvpi_);
            Ca = Vafzv();
            Ca->set_name("Alpha frozen virtual orbitals");
            if(transformationType_ != Restricted){
                View Vbfzv(Cb_, sopi_, fvir, zero,  mopi_ - frzvpi_);
                Cb = Vbfzv();
                Cb->set_name("Beta frozen virtual orbitals");
            }
        }else if(moSpace->label() == MOSPACE_NIL){
            // Do nothing!
        }else{
            char label = moSpace->label();
            if(aOrbsPI_.count(label) == 0){
                std::string err("Space ");
                err += label;
                err += " has not been initialized properly.";
                throw PSIEXCEPTION(err);
            }

            const std::vector<int> & aorbs = moSpace->aOrbs();
            const std::vector<int> & borbs = moSpace->bOrbs();
            int nAOrbs = aorbs.size();
            int nBOrbs = borbs.size();
            std::string name("Alpha orbitals for space ");
            name += label;
            Ca = SharedMatrix(new Matrix(name, nirreps_, sopi_, aOrbsPI_[label]));
            int mo_offsets[8];
            mo_offsets[0] = 0;
            for(int h = 1; h < nirreps_; ++h)
                mo_offsets[h] = mo_offsets[h-1] + mopi_[h-1];

            for(int h = 0; h < nirreps_; ++h){
                int count = 0;
                for(int n = 0; n < nAOrbs; ++n){
                    int orb = aorbs[n];
                    if(mosym_[orb] == h){
                        for(int so = 0; so < sopi_[h]; ++so){
                            Ca->set(h, so, count, Ca_->get(h, so, orb - mo_offsets[h]));
                        }
                        count++;
                    }
                }
            }
            if(transformationType_ != Restricted){
                name = "Beta orbitals for space " + std::string(1, label);
                Cb = SharedMatrix(new Matrix(name, nirreps_, sopi_, bOrbsPI_[label]));
                for(int h = 0; h < nirreps_; ++h){
                    int count = 0;
                    for(int n = 0; n < nBOrbs; ++n){
                        int orb = borbs[n];
                        if(mosym_[orb] == h){
                            for(int so = 0; so < sopi_[h]; ++so)
                                Cb->set(h, so, count, Cb_->get(h, so, orb - mo_offsets[h]));
                            count++;
                        }
                    }
                }
            }
        }

        if(transformationType_ == Restricted) Cb = Ca;

        aMOCoefficients_[moSpace->label()] = Ca;
        bMOCoefficients_[moSpace->label()] = Cb;

        if(print_ > 5){
            outfile->Printf( "Orbitals for space %c:-\n",moSpace->label());
            Ca->print();
            if (transformationType_ != Restricted)
                Cb->print();
            outfile->Printf( "\n\n");
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
    outfile->Printf( "The DPD mappings used in this transformation:-\n");
    for(iter = dpdLookup_.begin(); iter != dpdLookup_.end(); ++iter){
        outfile->Printf( "Pair %-10s ID = %d\n", iter->first.c_str(), iter->second);
    }
}
