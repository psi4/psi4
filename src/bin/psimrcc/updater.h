#ifndef _psi_src_bin_psimrcc_updater_h_
#define _psi_src_bin_psimrcc_updater_h_

/**
 *  @file updater.h
 *  @ingroup (PSIMRCC)
*/

namespace psi{ namespace psimrcc{

class Hamiltonian;

/**
 *  @class Updater
 *  @brief Containts the procedure for updating the amplitudes
*/
class Updater{
public:
  Updater();
  virtual ~Updater();
  virtual void update(int cycle,Hamiltonian* heff) = 0;
  void zero_internal_amps();
  void zero_t1_internal_amps();
  void zero_internal_delta_amps();
};

class MkUpdater : public Updater{
public:
  MkUpdater();
  virtual ~MkUpdater();
  virtual void update(int cycle,Hamiltonian* heff);
};

class BWUpdater : public Updater{
public:
  BWUpdater();
  virtual ~BWUpdater();
  virtual void update(int cycle,Hamiltonian* heff);
};

}} /* End Namespaces */

#endif // _psi_src_bin_psimrcc_updater_h_
