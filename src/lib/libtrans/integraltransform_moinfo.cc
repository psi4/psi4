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
    std::vector<shared_ptr<MOSpace> >::const_iterator space;

    //    for(int h = 0; h < _nirreps; ++h){
    //        fprintf(outfile, "docc = %d socc = %d frzcpi = %d frvirt = %d, mopi = %d, sopi = %d\n",
    //                _clsdpi[h], _openpi[h], _frzcpi[h], _frzvpi[h], _mopi[h], _sopi[h]);fflush(outfile);
    //    }


    bool qt_order = (moOrdering_ == QTOrder);  // If false, we assume Pitzer below

    for(space = uniqueSpaces_.begin(); space != uniqueSpaces_.end(); ++space){
        shared_ptr<MOSpace> moSpace = *space;
        int *aOrbsPI = new int[nirreps_];
        int *bOrbsPI = aOrbsPI;
        int *aIndex, *bIndex;
        int *aOrbSym, *bOrbSym;
        if(transformationType_ != Restricted){
            bOrbsPI = new int[nirreps_];
        }


        if(moSpace->label() == MOSPACE_FZC){
            // This is the frozen occupied space
            int numAOcc = 0, numBOcc = 0, aOccCount = 0, bOccCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                aOrbsPI[h] = frzcpi_[h];
                numAOcc += aOrbsPI[h];
                if(transformationType_ != Restricted){
                    bOrbsPI[h] = frzcpi_[h];
                    numBOcc += bOrbsPI[h];
                }
            }
            bOrbSym = aOrbSym = new int[numAOcc];
            bIndex  = aIndex  = new int[numAOcc];
            if(transformationType_ != Restricted){
                bOrbSym = new int[numBOcc];
                bIndex  = new int[numBOcc];
            }
            // Build the reindexing arrays for Pitzer ordering
            int aPitzerCount = 0, bPitzerCount = 0, aOrbCount = 0, bOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                aPitzerCount = bPitzerCount = pitzerOffset;
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[aPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aOccCount++] = h;
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bOccCount++] = h;
            }
        }else if(moSpace->label() == MOSPACE_OCC){
            // This is the occupied space
            int numAOcc = 0, numBOcc = 0, aOccCount = 0, bOccCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                aOrbsPI[h] = clsdpi_[h] + openpi_[h] - frzcpi_[h];
                numAOcc += aOrbsPI[h];
                if(transformationType_ != Restricted){
                    bOrbsPI[h] = clsdpi_[h] - frzcpi_[h];
                    numBOcc += bOrbsPI[h];
                }
            }
            bOrbSym = aOrbSym = new int[numAOcc];
            bIndex  = aIndex  = new int[numAOcc];
            if(transformationType_ != Restricted){
                bOrbSym = new int[numBOcc];
                bIndex  = new int[numBOcc];
            }
            // Build the reindexing arrays for Pitzer ordering
            int aPitzerCount = 0, bPitzerCount = 0, aOrbCount = 0, bOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                aPitzerCount = bPitzerCount = pitzerOffset + frzcpi_[h];
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[aPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aOccCount++] = h;
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bOccCount++] = h;
            }
        }else if(moSpace->label() == MOSPACE_ALL){
            // This is the full MO space
            int numActMO = 0;
            for(int h = 0; h < nirreps_; ++h){
                bOrbsPI[h] = aOrbsPI[h] = mopi_[h] - frzcpi_[h] - frzvpi_[h];
                numActMO += aOrbsPI[h];
            }
            bOrbSym = new int[numActMO];
            aOrbSym = new int[numActMO];
            aIndex  = new int[numActMO];
            bIndex  = new int[numActMO];
            // Build the reindexing arrays and orbital symmetries
            int aOrbCount = 0, bOrbCount = 0, pitzerCount = 0, pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                pitzerCount = pitzerOffset + frzcpi_[h];
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    bIndex[aOrbCount] = aIndex[aOrbCount] = (qt_order ? aQT_[pitzerCount] : pitzerCount);
                    bOrbSym[aOrbCount] = aOrbSym[aOrbCount] = h;
                    aOrbCount++;
                    if(transformationType_ != Restricted){
                        bIndex[bOrbCount] = (qt_order ? bQT_[pitzerCount] : pitzerCount);
                        bOrbSym[bOrbCount++] = h;
                    }
                    pitzerCount++;
                }
                pitzerOffset += mopi_[h];
            }
        }else if(moSpace->label() == MOSPACE_VIR){
            // This is the virtual space
            int numAVir = 0, numBVir = 0, aVirCount = 0, bVirCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                if(transformationType_ == Restricted){
                    bOrbsPI[h] = aOrbsPI[h] = mopi_[h] - clsdpi_[h] - frzvpi_[h];
                }else{
                    aOrbsPI[h] = mopi_[h] - clsdpi_[h] - frzvpi_[h] - openpi_[h];
                    bOrbsPI[h] = mopi_[h] - clsdpi_[h] - frzvpi_[h];
                }
                numAVir += aOrbsPI[h];
                numBVir += bOrbsPI[h];
            }
            bOrbSym = aOrbSym = new int[numAVir];
            bIndex  = aIndex  = new int[numAVir];
            if(transformationType_ != Restricted){
                bOrbSym = new int[numBVir];
                bIndex  = new int[numBVir];
            }
            // Build the reindexing arrays
            int aPitzerCount = 0, bPitzerCount = 0, aOrbCount = 0, bOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                if(transformationType_ == Restricted){
                    aPitzerCount = bPitzerCount = pitzerOffset + clsdpi_[h];
                }else{
                    aPitzerCount = pitzerOffset + clsdpi_[h];
                    bPitzerCount = pitzerOffset + clsdpi_[h] + openpi_[h];
                }
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[aPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aVirCount++] = h;
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bVirCount++] = h;
            }
        }else if(moSpace->label() == MOSPACE_FZV){
            // This is the frozen virtual space
            int numAVir = 0, numBVir = 0, aVirCount = 0, bVirCount = 0;
            for(int h = 0; h < nirreps_; ++h){
                aOrbsPI[h] = frzvpi_[h];
                numAVir += aOrbsPI[h];
                if(transformationType_ != Restricted){
                    bOrbsPI[h] = frzvpi_[h];
                    numBVir += bOrbsPI[h];
                }
            }
            bOrbSym = aOrbSym = new int[numAVir];
            bIndex  = aIndex  = new int[numAVir];
            if(transformationType_ != Restricted){
                bOrbSym = new int[numBVir];
                bIndex  = new int[numBVir];
            }
            // Build the reindexing arrays
            int aPitzerCount = 0, bPitzerCount = 0, aOrbCount = 0, bOrbCount = 0;
            int pitzerOffset = 0;
            for(int h = 0; h < nirreps_; ++h){
                aPitzerCount = bPitzerCount = pitzerOffset + mopi_[h] - frzvpi_[h];
                for(int n = 0; n < aOrbsPI[h]; ++n){
                    aIndex[aOrbCount++] = (qt_order ? aQT_[aPitzerCount] : aPitzerCount);
                    aPitzerCount++;
                }
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n){
                        bIndex[bOrbCount++] = (qt_order ? bQT_[aPitzerCount] : bPitzerCount);
                        bPitzerCount++;
                    }
                pitzerOffset += mopi_[h];
            }
            // Compute the orbital symmetries
            for(int h = 0; h < nirreps_; ++h){
                for(int n = 0; n < aOrbsPI[h]; ++n)  aOrbSym[aVirCount++] = h;
                if(transformationType_ != Restricted)
                    for(int n = 0; n < bOrbsPI[h]; ++n)  bOrbSym[bVirCount++] = h;
            }
        }else{
            // This must be a custom MOSpace that the user provided
            const std::vector<int> aorbs  = moSpace->aOrbs();
            const std::vector<int> borbs  = moSpace->bOrbs();
            const std::vector<int> aindex = moSpace->aIndex();
            const std::vector<int> bindex = moSpace->bIndex();

            // Figure out how many orbitals per irrep, and group all orbitals by irrep
            int nAOrbs = aorbs.size();
            int nBOrbs = borbs.size();
            bOrbSym = aOrbSym = new int[nAOrbs];
            ::memset(aOrbsPI, '\0', nirreps_*sizeof(int));
            bIndex = aIndex = (aindex.empty() ? 0 : new int[nAOrbs]);
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
            if(transformationType_ != Restricted){
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
            }
            // Check that the indexing array was provided, if needed
            if(useIWL_ && aindex.empty()){
                std::string error("You must provide an indexing array for space ");
                error += moSpace->label();
                error += " or disable IWL output by changing OutputType.";
                throw SanityCheckError(error, __FILE__, __LINE__);
            }
        }

        if(print_ > 5){
            int nAOrbs = 0, nBOrbs = 0;
            fprintf(outfile, "Adding arrays for space %c:-\n",moSpace->label());
            fprintf(outfile, "\n\talpha orbsPI = ");
            for(int h = 0; h < nirreps_; nAOrbs += aOrbsPI[h], ++h)
                fprintf(outfile, "%d ", aOrbsPI[h]);
            fprintf(outfile, "\n\tbeta orbsPI = ");
            for(int h = 0; h < nirreps_; nBOrbs += bOrbsPI[h], ++h)
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

        spacesUsed_.push_back(toupper(moSpace->label()));
        spaceArray_.push_back(aOrbsPI);
        spaceArray_.push_back(aOrbSym);
        if(transformationType_ != Restricted){
            spacesUsed_.push_back(tolower(moSpace->label()));
            spaceArray_.push_back(bOrbsPI);
            spaceArray_.push_back(bOrbSym);
        }
        aOrbsPI_[moSpace->label()]  = aOrbsPI;
        bOrbsPI_[moSpace->label()]  = bOrbsPI;
        aIndices_[moSpace->label()] = aIndex;
        bIndices_[moSpace->label()] = bIndex;
    }// End loop over spaces

    // End by adding the AO orbital space - this is always needed
    spacesUsed_.push_back(MOSPACE_NIL);
    spaceArray_.push_back(sopi_);
    spaceArray_.push_back(sosym_);
    aOrbsPI_[MOSPACE_NIL] = zeros_;
    bOrbsPI_[MOSPACE_NIL] = zeros_;

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
 * Sets up the eigenvectors for the transformation by querying the MO spaces
 * passed into the constructor.  This is done seperately to the DPD setup
 * to give us a chance to semicanonicalize the orbitals if necessary.
 */
void
IntegralTransform::process_eigenvectors()
{
    std::vector<shared_ptr<MOSpace> >::const_iterator space;

    if(print_ > 4){
        Ca_->print();
        Cb_->print();
    }

    // N.B. The frozen orbitals have been zeroed, if appropriate
    Dimension focc = frzcpi_;
    Dimension aocc = clsdpi_ + openpi_ - frzcpi_;
    Dimension bocc = clsdpi_ - frzcpi_;
    Dimension avir = mopi_ - clsdpi_ - openpi_ - frzvpi_;
    Dimension bvir = mopi_ - clsdpi_ - frzvpi_;
    Dimension aall = mopi_ - frzcpi_ - frzvpi_;
    Dimension fvir = frzvpi_;
    Dimension ball = mopi_ - frzcpi_ - frzvpi_;
    Dimension zero = Dimension(nirreps_);

    for(space = uniqueSpaces_.begin(); space != uniqueSpaces_.end(); ++space){
        shared_ptr<MOSpace> moSpace = *space;
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
                // This is a slightly strange one, but we actually take the beta
                // orbitals here, so that the singly occupied orbitals are included.
                // Makes no difference for closed shell cases, of course.
                View Vavir(Ca_, sopi_, bvir, zero, clsdpi_);
                Ca = Vavir();
                Ca->set_name("Alpha virtual orbitals");
            }else{
                View Vavir(Ca_, sopi_, avir, zero, clsdpi_ + openpi_);
                Ca = Vavir();
                Ca->set_name("Alpha virtual orbitals");
                View Vbvir(Cb_, sopi_, bvir, zero, clsdpi_);
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
                name = "Beta orbitals for space " + label;
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
            fprintf(outfile, "Orbitals for space %c:-\n",moSpace->label());
            Ca->print();
            if (transformationType_ != Restricted)
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
    for(iter = dpdLookup_.begin(); iter != dpdLookup_.end(); ++iter){
        fprintf(outfile, "Pair %-10s ID = %d\n", iter->first.c_str(), iter->second);
    }
}
