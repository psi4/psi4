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

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/view.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libdpd/dpd.h"


using std::vector;
namespace psi{ namespace cctransort {

int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);

void memcheck(int reference);

vector<int> pitzer2qt(vector<Dimension> &spaces);

void sort_tei_rhf(std::shared_ptr<PSIO> psio, int print);
void sort_tei_uhf(std::shared_ptr<PSIO> psio, int print);

void c_sort(int reference);
void d_sort(int reference);
void e_sort(int reference);
void f_sort(int reference);
void b_spinad(std::shared_ptr<PSIO>);
void a_spinad();
void d_spinad();
void e_spinad();

void fock_rhf(std::shared_ptr<Wavefunction> ref, Dimension &occpi, Dimension &openpi,
              Dimension &virpi, Dimension &frzcpi, int print);
void fock_uhf(std::shared_ptr<Wavefunction> ref, Dimension &aoccpi, Dimension &boccpi,
              Dimension &avirpi, Dimension &bvirpi, Dimension &frzcpi, int print);

double scf_check(int reference, Dimension &openpi);

void denom_rhf(Dimension &openpi);
void denom_uhf();

PsiReturnType cctransort(SharedWavefunction ref, Options& options)
{
  tstart();

  std::shared_ptr<PSIO> psio(_default_psio_lib_);
  if(!ref) throw PSIEXCEPTION("SCF has not been run yet!");

  int print = options.get_int("PRINT");

  bool semicanonical = false;
  int reference = 0;
  if(options.get_str("REFERENCE") =="RHF") reference = 0;
  else if(options.get_str("REFERENCE") =="ROHF" &&
          (options.get_str("WFN")=="MP2" || options.get_str("WFN")=="CCSD_T" || options.get_str("WFN")=="CCSD_AT" ||
           options.get_str("WFN")=="CC3" || options.get_str("WFN")=="EOM_CC3" ||
           options.get_str("WFN")=="CC2" || options.get_str("WFN")=="EOM_CC2" ||
           options.get_str("WFN")=="BCCD" || options.get_str("WFN")=="BCCD_T")) {
    reference = 2;
    semicanonical = true;
  }
  else if(options.get_str("REFERENCE") == "ROHF") reference = 1;
  else if(options.get_str("REFERENCE") == "UHF") reference = 2;
  else {
    outfile->Printf("Invalid value of input keyword REFERENCE: %s\n", options.get_str("REFERENCE").c_str());
    throw PsiException("ccsort failure", __FILE__, __LINE__);
  }

  // Allow user to force semicanonical
  if(options["SEMICANONICAL"].has_changed()) {
   semicanonical = options.get_bool("SEMICANONICAL");
   if (semicanonical) reference = 2;
  }

  int nirreps = ref->nirrep();
  int nmo = ref->nmo();
  char **labels = ref->molecule()->irrep_labels();
  double enuc = ref->molecule()->nuclear_repulsion_energy();
  double escf;
  if(ref->reference_wavefunction())
      escf = ref->reference_wavefunction()->reference_energy();
  else
      escf = ref->reference_energy();

  Dimension nmopi = ref->nmopi();
  Dimension nsopi = ref->nsopi();
  Dimension frzcpi = ref->frzcpi();
  Dimension frzvpi = ref->frzvpi();
  Dimension openpi = ref->soccpi();
  Dimension clsdpi = ref->doccpi() - frzcpi;
  Dimension uoccpi = nmopi - clsdpi - openpi - frzcpi - frzvpi;

  Dimension aoccpi = clsdpi + openpi;
  Dimension boccpi = clsdpi;
  Dimension avirpi = uoccpi;
  Dimension bvirpi = uoccpi + openpi;
  Dimension occpi = clsdpi + openpi;
  Dimension virpi = uoccpi + openpi;

  // Build Pitzer->QT and QT->Pitzer reordering arrays
  vector<int> p2qt, p2qt_a, p2qt_b;
  vector<int> qt2p(nmo), qt2p_a(nmo), qt2p_b(nmo);
  vector<Dimension> subspaces;
  if(reference == 2) {
    subspaces.push_back(frzcpi);
    subspaces.push_back(aoccpi);
    subspaces.push_back(avirpi);
    subspaces.push_back(frzvpi);
    p2qt_a = pitzer2qt(subspaces);
    for(int i=0; i < nmo; i++) qt2p_a[p2qt_a[i]] = i;
    subspaces.clear();
    subspaces.push_back(frzcpi);
    subspaces.push_back(boccpi);
    subspaces.push_back(bvirpi);
    subspaces.push_back(frzvpi);
    p2qt_b = pitzer2qt(subspaces);
    for(int i=0; i < nmo; i++) qt2p_b[p2qt_b[i]] = i;
  }
  else {
    subspaces.push_back(frzcpi);
    subspaces.push_back(clsdpi);
    subspaces.push_back(openpi);
    subspaces.push_back(uoccpi);
    subspaces.push_back(frzvpi);
    p2qt = pitzer2qt(subspaces);
    for(int i=0; i < nmo; i++) qt2p[p2qt[i]] = i;
  }

  int nfzc = frzcpi.sum();
  int nfzv = frzvpi.sum();
  int nclsd = clsdpi.sum();
  int nopen = openpi.sum();
  int nuocc = uoccpi.sum();
  int nactive = nclsd + nopen + nuocc;

  psio->open(PSIF_CC_INFO, PSIO_OPEN_OLD);

  psio->write_entry(PSIF_CC_INFO, "Reference Wavefunction", (char *) &(reference), sizeof(int));
  psio->write_entry(PSIF_CC_INFO, "Frozen Core Orbs Per Irrep", (char *) (int *) frzcpi, sizeof(int)*nirreps);
  psio->write_entry(PSIF_CC_INFO, "Frozen Virt Orbs Per Irrep", (char *) (int *) frzvpi, sizeof(int)*nirreps);
  psio->write_entry(PSIF_CC_INFO, "No. of Active Orbitals", (char *) &(nactive), sizeof(int));

  // Build QT->CC and CC->QT reordering arrays
  vector<int> cc_aocc, cc_bocc, cc_avir, cc_bvir;
  vector<int> qt_aocc, qt_bocc, qt_avir, qt_bvir;
  vector<int> aocc_sym, bocc_sym, avir_sym, bvir_sym;
  vector<int> aocc_off, bocc_off, avir_off, bvir_off;
  vector<int> cc_occ, cc_vir;
  vector<int> qt_occ, qt_vir;
  vector<int> occ_sym, vir_sym;
  vector<int> occ_off, vir_off;
  if(reference == 2) { // UHF/semicanonical

    aocc_off.push_back(0);
    avir_off.push_back(0);
    bocc_off.push_back(0);
    bvir_off.push_back(0);
    int aocount = aoccpi[0];
    int avcount = avirpi[0];
    int bocount = boccpi[0];
    int bvcount = bvirpi[0];
    for(int h=1; h < nirreps; h++) {
      aocc_off.push_back(aocount); aocount += aoccpi[h];
      avir_off.push_back(avcount); avcount += avirpi[h];
      bocc_off.push_back(bocount); bocount += boccpi[h];
      bvir_off.push_back(bvcount); bvcount += bvirpi[h];
    }

    cc_aocc.assign(nactive, -1); cc_bocc.assign(nactive, -1);
    cc_avir.assign(nactive, -1); cc_bvir.assign(nactive, -1);
    qt_aocc.assign(nactive, -1); qt_bocc.assign(nactive, -1);
    qt_avir.assign(nactive, -1); qt_bvir.assign(nactive, -1);
    aocc_sym.assign(nactive, -1); bocc_sym.assign(nactive, -1);
    avir_sym.assign(nactive, -1); bvir_sym.assign(nactive, -1);
    for(int h=0, count=0, offset=0; h < nirreps; h++) {
      if(h) offset += clsdpi[h-1] + openpi[h-1];
      for(int i=0; i < clsdpi[h] + openpi[h]; i++, count++) {
        cc_aocc[offset+i] = count;
        qt_aocc[count] = nfzc + offset + i;
        aocc_sym[count] = h;
      }
    }
    for(int h=0, count=0, offset=0; h < nirreps; h++) {
      if(h) offset += clsdpi[h-1];
      for(int i=0; i < clsdpi[h]; i++, count++) {
        cc_bocc[offset+i] = count;
        qt_bocc[count] = nfzc + offset + i;
        bocc_sym[count] = h;
      }
    }
    for(int h=0, count=0, offset=nclsd+nopen; h < nirreps; h++) {
      if(h) offset += uoccpi[h-1];
      for(int i=0; i < uoccpi[h]; i++, count++) {
        cc_avir[offset+i] = count;
        qt_avir[count] = nfzc + offset + i;
        avir_sym[count] = h;
      }
    }
    for(int h=0, count=0, offset=nclsd; h < nirreps; h++) {
      if(h) offset += uoccpi[h-1] + openpi[h-1];
      for(int i=0; i < uoccpi[h] + openpi[h]; i++, count++) {
        cc_bvir[offset+i] = count;
        qt_bvir[count] = nfzc + offset + i;
        bvir_sym[count] = h;
      }
    }

    psio->write_entry(PSIF_CC_INFO, "Active Alpha Occ Orbs Per Irrep", (char *) (int *) aoccpi, sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Occ Orbs Per Irrep", (char *) (int *) boccpi, sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Virt Orbs Per Irrep", (char *) (int *) avirpi, sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Virt Orbs Per Irrep", (char *) (int *) bvirpi, sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Offsets", (char *) aocc_off.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Occ Orb Offsets", (char *) bocc_off.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Offsets", (char *) avir_off.data(),  sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Virt Orb Offsets", (char *) bvir_off.data(),  sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Occ Orb Symmetry", (char *) aocc_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Occ Orb Symmetry", (char *) bocc_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Alpha Virt Orb Symmetry", (char *) avir_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Beta Virt Orb Symmetry", (char *) bvir_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Alpha Active Occ Order", (char *) cc_aocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Beta Active Occ Order", (char *) cc_bocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Alpha Active Virt Order", (char *) cc_avir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Beta Active Virt Order", (char *) cc_bvir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Alpha Active Occ Order", (char *) qt_aocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Beta Active Occ Order", (char *) qt_bocc.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Alpha Active Virt Order", (char *) qt_avir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Beta Active Virt Order", (char *) qt_bvir.data(), sizeof(int)*nactive);
  }
  else { // RHF/ROHF
    occ_off.push_back(0);
    vir_off.push_back(0);
    int ocount = occpi[0];
    int vcount = virpi[0];
    for(int h=1; h < nirreps; h++) {
      occ_off.push_back(ocount); ocount += occpi[h];
      vir_off.push_back(vcount); vcount += virpi[h];
    }

    cc_occ.assign(nactive, -1); cc_vir.assign(nactive, -1);
    qt_occ.assign(nactive, -1); qt_vir.assign(nactive, -1);
    occ_sym.assign(nactive, -1); vir_sym.assign(nactive, -1);
    for(int h=0, count=0, cl_offset=0, op_offset=nclsd; h < nirreps; h++) {
      if(h) cl_offset += clsdpi[h-1];
      for(int i=0; i < clsdpi[h]; i++, count++) {
        cc_occ[cl_offset+i] = count;
        qt_occ[count] = nfzc + cl_offset+i;
        occ_sym[count] = h;
      }
      if(h) op_offset += openpi[h-1];
      for(int i=0; i < openpi[h]; i++, count++) {
        cc_occ[op_offset+i] = count;
        qt_occ[count] = nfzc + op_offset+i;
        occ_sym[count] = h;
      }
    }
    for(int h=0, count=0, vr_offset=nclsd+nopen, op_offset=nclsd; h < nirreps; h++) {
      if(h) vr_offset += uoccpi[h-1];
      for(int i=0; i < uoccpi[h]; i++, count++) {
        cc_vir[vr_offset+i] = count;
        qt_vir[count] = nfzc + vr_offset+i;
        vir_sym[count] = h;
      }
      if(h) op_offset += openpi[h-1];
      for(int i=0; i < openpi[h]; i++, count++) {
        cc_vir[op_offset+i] = count;
        qt_vir[count] = nfzc + op_offset+i;
        vir_sym[count] = h;
      }
    }

    psio->write_entry(PSIF_CC_INFO, "Active Occ Orbs Per Irrep", (char *) (int *) occpi, sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Virt Orbs Per Irrep", (char *) (int *) virpi, sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Occ Orb Offsets", (char *) occ_off.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Virt Orb Offsets", (char *) vir_off.data(), sizeof(int)*nirreps);
    psio->write_entry(PSIF_CC_INFO, "Active Occ Orb Symmetry", (char *) occ_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "Active Virt Orb Symmetry", (char *) vir_sym.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Active Occ Order", (char *) cc_occ.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "QT->CC Active Virt Order", (char *) cc_vir.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Active Occ Order", (char *) qt_occ.data(), sizeof(int)*nactive);
    psio->write_entry(PSIF_CC_INFO, "CC->QT Active Virt Order", (char *) qt_vir.data(), sizeof(int)*nactive);
  }

  psio->close(PSIF_CC_INFO, 1);

  outfile->Printf("\n\tWfn Parameters:\n");
  outfile->Printf("\t--------------------\n");
  outfile->Printf("\tWavefunction         = %s\n", options.get_str("WFN").c_str());
  outfile->Printf("\tNumber of irreps     = %d\n", nirreps);
  outfile->Printf("\tNumber of MOs        = %d\n", nmo);
  outfile->Printf("\tNumber of active MOs = %d\n", nactive);
  outfile->Printf("\tAO-Basis             = %s\n", options.get_str("AO_BASIS").c_str());
  outfile->Printf("\tSemicanonical        = %s\n", semicanonical ? "true" : "false");
  if(semicanonical)
    outfile->Printf("\tReference            = ROHF changed to UHF for semicanonical orbitals\n");
  else
    outfile->Printf("\tReference            = %s\n", reference == 2 ? "UHF" : (reference == 1 ? "ROHF" : "RHF"));
  outfile->Printf("\tPrint Level          = %d\n", print);
  outfile->Printf("\n");

  outfile->Printf("\tIRREP\t# MOs\t# FZDC\t# DOCC\t# SOCC\t# VIRT\t# FZVR\n");
  outfile->Printf(
            "\t-----\t-----\t------\t------\t------\t------\t------\n");
  for(int i=0; i < nirreps; i++) {
    outfile->Printf("\t %s\t   %d\t    %d\t    %d\t    %d\t    %d\t    %d\n",
                    labels[i],nmopi[i],frzcpi[i],clsdpi[i],openpi[i],uoccpi[i],frzvpi[i]);
  }

  // Transformation

  outfile->Printf("\tTransforming integrals...\n");
  std::vector<std::shared_ptr<MOSpace> > transspaces;
  transspaces.push_back(MOSpace::occ);
  transspaces.push_back(MOSpace::vir);

  IntegralTransform *ints;
  if(options.get_str("REFERENCE") == "RHF")
    ints = new IntegralTransform(ref, transspaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
  else if(options.get_str("REFERENCE") == "ROHF") {
    if(semicanonical)
      // Importantly the transform is handled python-side so we technically have unrestricted orbitals at this point
      ints = new IntegralTransform(ref, transspaces, IntegralTransform::Unrestricted, IntegralTransform::DPDOnly);
    else
      ints = new IntegralTransform(ref, transspaces, IntegralTransform::Restricted, IntegralTransform::DPDOnly);
  }
  else if(options.get_str("REFERENCE") == "UHF")
    ints = new IntegralTransform(ref, transspaces, IntegralTransform::Unrestricted, IntegralTransform::DPDOnly);
  else
    throw PSIEXCEPTION("Invalid choice of reference wave function.");

  dpd_set_default(ints->get_dpd_id());
  ints->set_keep_dpd_so_ints(true);
  if(!options.get_bool("DELETE_TEI") || options.get_str("AO_BASIS") == "DISK") {
    outfile->Printf("\tIWL integrals will be retained.\n");
    ints->set_keep_iwl_so_ints(true);
  }
  else {
    outfile->Printf("\tIWL integrals will be deleted.\n");
    ints->set_keep_iwl_so_ints(false);
  }

  // On the second and later passes of Brueckner, the presort is already done
  // TDC: Always re-compute the presorted integrals until the frozen-core operator is included
  bool presort_predone = false;
  /*
  if(psio->tocentry_exists(PSIF_SO_PRESORT, "SO Ints (nn|nn)")) {
    outfile->Printf("\tPresorted integrals already available.\n");
    ints->set_tei_already_presorted(true);
    presort_predone = true;
  }
  else {
    outfile->Printf("\tPresorted integrals will be generated.\n");
    ints->set_tei_already_presorted(false);
  }
  */

  outfile->Printf("\t(OO|OO)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
  outfile->Printf("\t(OO|OV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
  outfile->Printf("\t(OO|VV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::occ, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);

  outfile->Printf("\t(OV|OO)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
  outfile->Printf("\t(OV|OV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
  outfile->Printf("\t(OV|VV)...\n");
  ints->transform_tei(MOSpace::occ, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);

  if(options.get_bool("DELETE_TEI")) ints->set_keep_dpd_so_ints(false);

  outfile->Printf("\t(VV|OO)...\n");
  ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::occ, IntegralTransform::MakeAndKeep);
  outfile->Printf("\t(VV|OV)...\n");
  ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::occ, MOSpace::vir, IntegralTransform::ReadAndKeep);
  outfile->Printf("\t(VV|VV)...\n");
  ints->transform_tei(MOSpace::vir, MOSpace::vir, MOSpace::vir, MOSpace::vir, IntegralTransform::ReadAndNuke);

  double efzc;
  psio->open(PSIF_CC_INFO, PSIO_OPEN_OLD);
  if(presort_predone) psio->read_entry(PSIF_CC_INFO, "Frozen-Core Energy", (char *) &(efzc), sizeof(double));
  else {
    efzc = ints->get_frozen_core_energy();
    psio->write_entry(PSIF_CC_INFO, "Frozen-Core Energy", (char *) &(efzc), sizeof(double));
  }
  psio->close(PSIF_CC_INFO, 1);
  outfile->Printf(  "\tFrozen core energy     =  %20.14f\n", efzc);

  if(nfzc && (fabs(efzc) < 1e-7)) {
    outfile->Printf( "\tCCSORT Error: Orbitals are frozen in input,\n");
    outfile->Printf( "\tbut frozen core energy is small!\n");
    outfile->Printf( "\tCalculation will be aborted...\n");
    exit(PSI_RETURN_FAILURE);
  }
  else if(!nfzc && fabs(efzc)) {
    outfile->Printf( "\tCCSORT Warning: No orbitals are frozen,\n");
    outfile->Printf( "\tbut the frozen-core energy in wfn is non-zero.\n");
    outfile->Printf( "\tCalculation will continue with zero efzc...\n");
    efzc = 0.0;
  }

  delete ints;

  // Set up DPD object
  std::vector<DPDMOSpace> spaces;
  int *cachefiles = init_int_array(PSIO_MAXUNIT);
  int **cachelist;
  if(reference == 2) { // UHF/semicanonical
    cachelist = cacheprep_uhf(options.get_int("CACHELEVEL"), cachefiles);

    DPDMOSpace aocc('O', "IJKLMN", aoccpi);
    DPDMOSpace avir('V', "ABCDEF", avirpi);
    DPDMOSpace bocc('o', "ijklmn", boccpi);
    DPDMOSpace bvir('v', "abcdef", bvirpi);
    spaces.push_back(aocc);
    spaces.push_back(avir);
    spaces.push_back(bocc);
    spaces.push_back(bvir);

    if(dpd_list[0]) throw PSIEXCEPTION("Attempting to initilize new DPD instance before the old one was freed.");
    dpd_list[0] = new DPD(0, nirreps, Process::environment.get_memory(), 0, cachefiles, cachelist, NULL, 4, spaces);
    dpd_default = 0;
    global_dpd_ = dpd_list[0];
  }
  else { // RHF/ROHF
    cachelist = cacheprep_rhf(options.get_int("CACHELEVEL"), cachefiles);

    DPDMOSpace occ('o', "ijklmn", occpi);
    DPDMOSpace vir('v', "abcdef", virpi);
    spaces.push_back(occ);
    spaces.push_back(vir);

    if(dpd_list[0]) throw PSIEXCEPTION("Attempting to initilize new DPD instance before the old one was freed.");
    dpd_list[0] = new DPD(0, nirreps, Process::environment.get_memory(), 0, cachefiles, cachelist, NULL, 2, spaces);
    dpd_default = 0;
    global_dpd_ = dpd_list[0];
  }

  memcheck(reference);

  // Sort two-electron integrals into six main categories
  psio->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
  if(reference == 2) sort_tei_uhf(psio, print);
  else sort_tei_rhf(psio, print);
  psio->close(PSIF_LIBTRANS_DPD, 0); // delete file

  for(int i =PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio->open(i,1);

  // Sort one-electron integals into three main categories

  if(reference == 2) { // UHF/Semicanonical
    dpdfile2 H;
    int ntri_all = nmo * (nmo + 1)/2;
    double *tmp_oei = init_array(ntri_all);

    iwl_rdone(PSIF_OEI, PSIF_MO_A_FZC, tmp_oei, ntri_all, 0, 0, "outfile");

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 0, 0, "h(I,J)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int i=0; i < aoccpi[h]; i++) {
        for(int j=0; j < aoccpi[h]; j++) {
          int I = qt2p_a[qt_aocc[aocc_off[h] + i]];
          int J = qt2p_a[qt_aocc[aocc_off[h] + j]];
          H.matrix[h][i][j] = tmp_oei[INDEX(I,J)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 1, 1, "h(A,B)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int a=0; a < avirpi[h]; a++) {
        for(int b=0; b < avirpi[h]; b++) {
          int A = qt2p_a[qt_avir[avir_off[h] + a]];
          int B = qt2p_a[qt_avir[avir_off[h] + b]];
          H.matrix[h][a][b] = tmp_oei[INDEX(A,B)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 0, 1, "h(I,A)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int i=0; i < aoccpi[h]; i++) {
        for(int a=0; a < avirpi[h]; a++) {
          int I = qt2p_a[qt_aocc[aocc_off[h] + i]];
          int A = qt2p_a[qt_avir[avir_off[h] + a]];
          H.matrix[h][i][a] = tmp_oei[INDEX(I,A)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    iwl_rdone(PSIF_OEI, PSIF_MO_B_FZC, tmp_oei, ntri_all, 0, 0, "outfile");

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 2, 2, "h(i,j)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int i=0; i < boccpi[h]; i++) {
        for(int j=0; j < boccpi[h]; j++) {
          int I = qt2p_b[qt_bocc[bocc_off[h] + i]];
          int J = qt2p_b[qt_bocc[bocc_off[h] + j]];
          H.matrix[h][i][j] = tmp_oei[INDEX(I,J)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 3, 3, "h(a,b)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int a=0; a < bvirpi[h]; a++) {
        for(int b=0; b < bvirpi[h]; b++) {
          int A = qt2p_b[qt_bvir[bvir_off[h] + a]];
          int B = qt2p_b[qt_bvir[bvir_off[h] + b]];
          H.matrix[h][a][b] = tmp_oei[INDEX(A,B)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 2, 3, "h(i,a)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int i=0; i < boccpi[h]; i++) {
        for(int a=0; a < bvirpi[h]; a++) {
          int I = qt2p_b[qt_bocc[bocc_off[h] + i]];
          int A = qt2p_b[qt_bvir[bvir_off[h] + a]];
          H.matrix[h][i][a] = tmp_oei[INDEX(I,A)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);
  }
  else { // RHF/ROHF
    dpdfile2 H;
    int ntri_all = nmo * (nmo + 1)/2;
    double *tmp_oei = init_array(ntri_all);

    iwl_rdone(PSIF_OEI, PSIF_MO_FZC, tmp_oei, ntri_all, 0, 0, "outfile");

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 0, 0, "h(i,j)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int i=0; i < occpi[h]; i++) {
        for(int j=0; j < occpi[h]; j++) {
          int I = qt2p[qt_occ[occ_off[h] + i]];
          int J = qt2p[qt_occ[occ_off[h] + j]];
          H.matrix[h][i][j] = tmp_oei[INDEX(I,J)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 1, 1, "h(a,b)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int a=0; a < virpi[h]; a++) {
        for(int b=0; b < virpi[h]; b++) {
          int A = qt2p[qt_vir[vir_off[h] + a]];
          int B = qt2p[qt_vir[vir_off[h] + b]];
          H.matrix[h][a][b] = tmp_oei[INDEX(A,B)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);

    global_dpd_->file2_init(&H, PSIF_CC_OEI, 0, 0, 1, "h(i,a)");
    global_dpd_->file2_mat_init(&H);
    for(int h=0; h < nirreps; h++) {
      for(int i=0; i < occpi[h]; i++) {
        for(int a=0; a < virpi[h]; a++) {
          int I = qt2p[qt_occ[occ_off[h] + i]];
          int A = qt2p[qt_vir[vir_off[h] + a]];
          H.matrix[h][i][a] = tmp_oei[INDEX(I,A)];
        }
      }
    }
    global_dpd_->file2_mat_wrt(&H);
    global_dpd_->file2_mat_close(&H);
    global_dpd_->file2_close(&H);
  }

  // Generate additional orderings of basic integrals
  c_sort(reference);
  d_sort(reference);
  e_sort(reference);
  f_sort(reference);
  if(reference == 0) {
    b_spinad(psio);
    a_spinad();
    d_spinad();
    e_spinad();
  }

  // Organize Fock matrices
  if(reference == 2) fock_uhf(ref, aoccpi, boccpi, avirpi, bvirpi, frzcpi, print);
  else fock_rhf(ref, occpi, openpi, virpi, frzcpi, print);

  if(reference == 2) denom_uhf();
  else denom_rhf(openpi);

  outfile->Printf("\tNuclear Rep. energy          =  %20.14f\n", enuc);
  outfile->Printf("\tSCF energy                   =  %20.14f\n", escf);
  double eref = scf_check(reference, openpi) + enuc + efzc;
  outfile->Printf("\tReference energy             =  %20.14f\n", eref);
  psio->write_entry(PSIF_CC_INFO, "Reference Energy", (char *) &(eref), sizeof(double));

  SharedMatrix Ca_vir, Cb_vir, Ca_occ;
  Dimension zero(nirreps);
  psio_address next;
  if(reference == 2) {
    View VCa_vir(ref->Ca(), nsopi, avirpi, zero, aoccpi+frzcpi);
    Ca_vir = VCa_vir();
    Ca_vir->set_name("Alpha virtual orbitals");

    next = PSIO_ZERO;
    for(int h=0; h < nirreps; h++)
      if(avirpi[h])
        psio->write(PSIF_CC_INFO, "UHF Active Alpha Virtual Orbs", (char *) Ca_vir->pointer(h)[0],
                    nsopi[h]*avirpi[h]*sizeof(double), next, &next);

    View VCb_vir(ref->Cb(), nsopi, bvirpi, zero, boccpi+frzcpi);
    Cb_vir = VCb_vir();
    Cb_vir->set_name("Beta virtual orbitals");

    next = PSIO_ZERO;
    for(int h=0; h < nirreps; h++)
      if(bvirpi[h])
        psio->write(PSIF_CC_INFO, "UHF Active Beta Virtual Orbs", (char *) Cb_vir->pointer(h)[0],
                    nsopi[h]*bvirpi[h]*sizeof(double), next, &next);
  }
  else {
    View VCa_occ(ref->Ca(), nsopi, occpi, zero, frzcpi);
    Ca_occ = VCa_occ();
    Ca_occ->set_name("Occupied orbitals");

    next = PSIO_ZERO;
    for(int h=0; h < nirreps; h++)
      if(occpi[h])
        psio->write(PSIF_CC_INFO, "RHF/ROHF Active Occupied Orbitals", (char *) Ca_occ->pointer(h)[0],
                    nsopi[h]*occpi[h]*sizeof(double), next, &next);

    View Vavir(ref->Ca(), nsopi, uoccpi, zero, frzcpi+clsdpi+openpi);
    View Vasoc(ref->Ca(), nsopi, openpi, zero, frzcpi+clsdpi);
    std::vector<SharedMatrix> virandsoc;
    virandsoc.push_back(Vavir());
    virandsoc.push_back(Vasoc());
    Ca_vir = Matrix::horzcat(virandsoc);
    Ca_vir->set_name("Virtual orbitals");

    next = PSIO_ZERO;
    for(int h=0; h < nirreps; h++)
      if(virpi[h])
        psio->write(PSIF_CC_INFO, "RHF/ROHF Active Virtual Orbitals", (char *) Ca_vir->pointer(h)[0],
                    nsopi[h]*virpi[h]*sizeof(double), next, &next);
  }

  dpd_close(0);
  if(reference == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  for(int i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio->close(i,1);
  for(int i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio->close(i,0); /* delete CC_TMP files */
  for(int i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio->close(i,1);

  tstop();

  return Success;
}

}} // End namespaces
