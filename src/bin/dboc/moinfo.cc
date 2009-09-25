/*! \file
    \ingroup DBOC
    \brief Enter brief description of file here 
*/

#include <stdexcept>
#include <libchkpt/chkpt.h>
#include <libipv1/ip_lib.h>
#include <libqt/qt.h>
#include "moinfo.h"

namespace psi { namespace dboc {

extern MOInfo_t MOInfo;

void read_moinfo()
{
  chkpt_init(PSIO_OPEN_OLD);
  int nirreps = chkpt_rd_nirreps();
  //  if (nirreps != 1)
  //    done("DBOC computations currently possible only in C1 symmetry");

  MOInfo.num_so = chkpt_rd_nso();
  MOInfo.num_mo = chkpt_rd_nmo();
  int* orbspi = chkpt_rd_orbspi();
  int* clsdpi = chkpt_rd_clsdpi();
  MOInfo.ndocc = 0;
  for(int irrep=0; irrep<nirreps; irrep++)
    MOInfo.ndocc += clsdpi[irrep];
  int* openpi = chkpt_rd_openpi();
  MOInfo.nsocc = 0;
  for(int irrep=0; irrep<nirreps; irrep++)
    MOInfo.nsocc += openpi[irrep];
  MOInfo.nalpha = MOInfo.ndocc + MOInfo.nsocc;
  MOInfo.nbeta = MOInfo.ndocc;

  chkpt_close();

  delete[] clsdpi;
  delete[] openpi;
  delete[] orbspi;
}

void mo_maps(short int** qt2p_map, short int** act_qt2p_map)
{
  chkpt_init(PSIO_OPEN_OLD);
  int nirreps = chkpt_rd_nirreps();
  int* orbspi = chkpt_rd_orbspi();
  int* clsdpi = chkpt_rd_clsdpi();
  int* openpi = chkpt_rd_openpi();

  // Parse frozen_docc and frozen_uocc arrays from input file
  int* nfzcpi = get_frzcpi();
  int* nfzvpi = get_frzvpi();
  MOInfo.nfzc = 0;
  MOInfo.nfzv = 0;
  for(int irr=0; irr<nirreps; irr++) {
    MOInfo.nfzc += nfzcpi[irr];
    MOInfo.nfzv += nfzvpi[irr];
  }
  chkpt_close();

  if (MOInfo.nfzc >= MOInfo.num_mo)
    throw std::runtime_error("ERROR: mo_maps -- nfzc > nmo");
  if (MOInfo.nfzv >= MOInfo.num_mo)
    throw std::runtime_error("ERROR: mo_maps -- nfzv > nmo");
  if ( (MOInfo.nfzc + MOInfo.nfzv) >= MOInfo.num_mo)
    throw std::runtime_error("ERROR: mo_maps -- nfzv + nfzc > nmo");

  MOInfo.nact = MOInfo.num_mo - MOInfo.nfzc - MOInfo.nfzv;

  short int* aq2pmap = new short int[MOInfo.nact];
  short int* q2pmap = new short int[MOInfo.num_mo];

  //
  // Compute the map
  //
  int QTS_count = 0;
  int act_QTS_count = 0;

  // Add frozen doccs only to QTS_count
  int pitzer_offset = 0;
  for(int irr=0; irr<nirreps; irr++) {
    int nfzc = nfzcpi[irr];
    for(int mo=0; mo<nfzc; mo++)
      q2pmap[QTS_count++] = pitzer_offset++;
    pitzer_offset += orbspi[irr] - nfzcpi[irr];
  }
  // Add active doccs
  pitzer_offset = 0;
  for(int irr=0; irr<nirreps; irr++) {
    int nact = clsdpi[irr] - nfzcpi[irr];
    pitzer_offset += nfzcpi[irr];
    for(int mo=0; mo<nact; mo++) {
      q2pmap[QTS_count++] = pitzer_offset;
      aq2pmap[act_QTS_count++] = pitzer_offset++;
    }
    pitzer_offset += orbspi[irr] - clsdpi[irr];
  }
  // Add soccs
  pitzer_offset = 0;
  for(int irr=0; irr<nirreps; irr++) {
    int nact = openpi[irr];
    pitzer_offset += clsdpi[irr];
    for(int mo=0; mo<nact; mo++) {
      q2pmap[QTS_count++] = pitzer_offset;
      aq2pmap[act_QTS_count++] = pitzer_offset++;
    }
    pitzer_offset += orbspi[irr] - openpi[irr] - clsdpi[irr];
  }
  // Add active uoccs
  pitzer_offset = 0;
  for(int irr=0; irr<nirreps; irr++) {
    int nact = orbspi[irr] - clsdpi[irr] - openpi[irr] - nfzvpi[irr];
    pitzer_offset += clsdpi[irr] + openpi[irr];
    for(int mo=0; mo<nact; mo++) {
      q2pmap[QTS_count++] = pitzer_offset;
      aq2pmap[act_QTS_count++] = pitzer_offset++;
    }
    pitzer_offset += nfzvpi[irr];
  }
  // Add active uoccs to QTS_count
  pitzer_offset = 0;
  for(int irr=0; irr<nirreps; irr++) {
    int nfzv = nfzvpi[irr];
    pitzer_offset += orbspi[irr] - nfzv;
    for(int mo=0; mo<nfzv; mo++) {
      q2pmap[QTS_count++] = pitzer_offset++;
    }
  }

  delete[] clsdpi;
  delete[] openpi;
  delete[] orbspi;

  *qt2p_map = q2pmap;
  *act_qt2p_map = aq2pmap;
}

}} // namespace psi::dboc
