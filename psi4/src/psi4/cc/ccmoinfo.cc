/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "ccmoinfo.h"

#include <vector>

#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/dimension.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

// FIXME Figure out a way to not duplicate reference info in 3 places
#include "ccparams.h"

namespace psi {
namespace cc {

std::vector<int> get_pitzer2qt(std::vector<Dimension> &spaces) {
    int nirreps = spaces[0].n();

    Dimension total(nirreps);
    for (int h = 0; h < nirreps; h++)
        for (int i = 0; i < spaces.size(); i++) total[h] += spaces[i][h];
    int nmo = total.sum();

    std::vector<int> order(nmo);
    order.assign(nmo, 0);

    Dimension offset(nirreps);
    offset[0] = 0;
    for (int h = 1; h < nirreps; h++) offset[h] = offset[h - 1] + total[h - 1];

    int count = 0;

    for (int j = 0; j < spaces.size(); j++)
        for (int h = 0; h < nirreps; h++) {
            int this_offset = offset[h];
            for (int k = 0; k < j; k++) this_offset += spaces[k][h];
            for (int i = 0; i < spaces[j][h]; i++) order[this_offset + i] = count++;
        }

    return order;
}

CCMOInfo::CCMOInfo(SharedWavefunction wfn, Reference ref) {
    nirreps = wfn->nirrep();
    nmo = wfn->nmo();
    nso = wfn->nso();
    nao = wfn->basisset()->nao();

    labels = wfn->molecule()->irrep_labels();

    enuc = wfn->molecule()->nuclear_repulsion_energy(wfn->get_dipole_field_strength());
    escf = 0.0;
    if (wfn->reference_wavefunction()) {
        escf = wfn->reference_wavefunction()->reference_energy();
    } else {
        escf = wfn->reference_energy();
    }
    // TODO we probably need to get the PCM polarization energy too?

    orbspi = wfn->nmopi();
    sopi = wfn->nsopi();
    frdocc = wfn->frzcpi();
    fruocc = wfn->frzvpi();
    clsdpi = wfn->doccpi() - frdocc;
    openpi = wfn->soccpi();
    uoccpi = orbspi - clsdpi - openpi - frdocc - fruocc;
    occpi = clsdpi + openpi;
    aoccpi = clsdpi + openpi;
    boccpi = clsdpi;
    virtpi = uoccpi + openpi;
    avirtpi = uoccpi;
    bvirtpi = uoccpi + openpi;
    nvirt = virtpi.sum();

    // Build Pitzer->QT and QT->Pitzer reordering arrays
    std::vector<Dimension> subspaces;
    if (ref == Reference::UHF) {
        subspaces.push_back(frdocc);
        subspaces.push_back(aoccpi);
        subspaces.push_back(avirtpi);
        subspaces.push_back(fruocc);
        pitzer2qt_a = get_pitzer2qt(subspaces);
        qt2pitzer_a.reserve(nmo);
        for (int i = 0; i < nmo; i++) qt2pitzer_a[pitzer2qt_a[i]] = i;
        subspaces.clear();
        subspaces.push_back(frdocc);
        subspaces.push_back(boccpi);
        subspaces.push_back(bvirtpi);
        subspaces.push_back(fruocc);
        pitzer2qt_b = get_pitzer2qt(subspaces);
        qt2pitzer_b.reserve(nmo);
        for (int i = 0; i < nmo; i++) qt2pitzer_b[pitzer2qt_b[i]] = i;
    } else {
        subspaces.push_back(frdocc);
        subspaces.push_back(clsdpi);
        subspaces.push_back(openpi);
        subspaces.push_back(uoccpi);
        subspaces.push_back(fruocc);
        pitzer2qt = get_pitzer2qt(subspaces);
        qt2pitzer.reserve(nmo);
        for (int i = 0; i < nmo; i++) qt2pitzer[pitzer2qt[i]] = i;
    }

    int nfzc = frdocc.sum();
    int nclsd = clsdpi.sum();
    int nopen = openpi.sum();
    int nuocc = uoccpi.sum();
    nactive = nclsd + nopen + nuocc;

    // Build QT->CC and CC->QT reordering arrays
    if (ref == Reference::UHF) {  // UHF/semicanonical
        aocc_off.push_back(0);
        avir_off.push_back(0);
        bocc_off.push_back(0);
        bvir_off.push_back(0);
        int aocount = aoccpi[0];
        int avcount = avirtpi[0];
        int bocount = boccpi[0];
        int bvcount = bvirtpi[0];
        for (int h = 1; h < nirreps; h++) {
            aocc_off.push_back(aocount);
            aocount += aoccpi[h];
            avir_off.push_back(avcount);
            avcount += avirtpi[h];
            bocc_off.push_back(bocount);
            bocount += boccpi[h];
            bvir_off.push_back(bvcount);
            bvcount += bvirtpi[h];
        }

        cc_aocc.assign(nactive, -1);
        cc_bocc.assign(nactive, -1);
        cc_avir.assign(nactive, -1);
        cc_bvir.assign(nactive, -1);
        qt_aocc.assign(nactive, -1);
        qt_bocc.assign(nactive, -1);
        qt_avir.assign(nactive, -1);
        qt_bvir.assign(nactive, -1);
        aocc_sym.assign(nactive, -1);
        bocc_sym.assign(nactive, -1);
        avir_sym.assign(nactive, -1);
        bvir_sym.assign(nactive, -1);
        for (int h = 0, count = 0, offset = 0; h < nirreps; h++) {
            if (h) offset += clsdpi[h - 1] + openpi[h - 1];
            for (int i = 0; i < clsdpi[h] + openpi[h]; i++, count++) {
                cc_aocc[offset + i] = count;
                qt_aocc[count] = nfzc + offset + i;
                aocc_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, offset = 0; h < nirreps; h++) {
            if (h) offset += clsdpi[h - 1];
            for (int i = 0; i < clsdpi[h]; i++, count++) {
                cc_bocc[offset + i] = count;
                qt_bocc[count] = nfzc + offset + i;
                bocc_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, offset = nclsd + nopen; h < nirreps; h++) {
            if (h) offset += uoccpi[h - 1];
            for (int i = 0; i < uoccpi[h]; i++, count++) {
                cc_avir[offset + i] = count;
                qt_avir[count] = nfzc + offset + i;
                avir_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, offset = nclsd; h < nirreps; h++) {
            if (h) offset += uoccpi[h - 1] + openpi[h - 1];
            for (int i = 0; i < uoccpi[h] + openpi[h]; i++, count++) {
                cc_bvir[offset + i] = count;
                qt_bvir[count] = nfzc + offset + i;
                bvir_sym[count] = h;
            }
        }
    } else {  // RHF/ROHF
        occ_off.push_back(0);
        vir_off.push_back(0);
        int ocount = occpi[0];
        int vcount = virtpi[0];
        for (int h = 1; h < nirreps; h++) {
            occ_off.push_back(ocount);
            ocount += occpi[h];
            vir_off.push_back(vcount);
            vcount += virtpi[h];
        }

        cc_occ.assign(nactive, -1);
        cc_vir.assign(nactive, -1);
        qt_occ.assign(nactive, -1);
        qt_vir.assign(nactive, -1);
        occ_sym.assign(nactive, -1);
        vir_sym.assign(nactive, -1);
        for (int h = 0, count = 0, cl_offset = 0, op_offset = nclsd; h < nirreps; h++) {
            if (h) cl_offset += clsdpi[h - 1];
            for (int i = 0; i < clsdpi[h]; i++, count++) {
                cc_occ[cl_offset + i] = count;
                qt_occ[count] = nfzc + cl_offset + i;
                occ_sym[count] = h;
            }
            if (h) op_offset += openpi[h - 1];
            for (int i = 0; i < openpi[h]; i++, count++) {
                cc_occ[op_offset + i] = count;
                qt_occ[count] = nfzc + op_offset + i;
                occ_sym[count] = h;
            }
        }
        for (int h = 0, count = 0, vr_offset = nclsd + nopen, op_offset = nclsd; h < nirreps; h++) {
            if (h) vr_offset += uoccpi[h - 1];
            for (int i = 0; i < uoccpi[h]; i++, count++) {
                cc_vir[vr_offset + i] = count;
                qt_vir[count] = nfzc + vr_offset + i;
                vir_sym[count] = h;
            }
            if (h) op_offset += openpi[h - 1];
            for (int i = 0; i < openpi[h]; i++, count++) {
                cc_vir[op_offset + i] = count;
                qt_vir[count] = nfzc + op_offset + i;
                vir_sym[count] = h;
            }
        }
    }

    // Build sosym array (for AO-basis BT2)
    sosym.reserve(nso);
    for (int h = 0, q = 0; h < nirreps; h++)
        for (int p = 0; p < sopi[h]; p++) sosym[q++] = h;
}

void print_ccmoinfo(const CCMOInfo &info) {
    outfile->Printf("\n\tWfn Parameters:\n");
    outfile->Printf("\t--------------------\n");
    outfile->Printf("\tNumber of irreps     = %d\n", info.nirreps);
    outfile->Printf("\tNumber of MOs        = %d\n", info.nmo);
    outfile->Printf("\tNumber of active MOs = %d\n", info.nactive);
    outfile->Printf("\n");

    outfile->Printf("\tIRREP\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
    outfile->Printf("\t-----\t-----\t------\t------\t------\t------\t------\n");
    for (int i = 0; i < info.nirreps; i++) {
        outfile->Printf("\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n", info.labels[i].c_str(),
                        info.orbspi[i], info.frdocc[i], info.clsdpi[i], info.openpi[i], info.uoccpi[i], info.fruocc[i]);
    }
    outfile->Printf("    Nuclear Rep. energy (wfn)     = %20.15f\n", info.enuc);
    outfile->Printf("    SCF energy          (wfn)     = %20.15f\n", info.escf);
}

}  // namespace cc
}  // namespace psi
