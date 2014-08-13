#ifndef __AccessMode_h__
#define __AccessMode_h__

namespace tascel {
  /**
   * Data access modes.
   *
   * The enumeration lists mode of access for data collections. This
   * information is used by the runtime to schedule communication.
   */
  enum AccessMode {
    MODE_RONLY, /**< Read-only accessmode*/
    MODE_RDWR , /**< Read at start and write at end*/
    MODE_ACC    /**< Accumulate result*/
  }; /* AccessMode*/

}; /*tascel*/

#endif /*__AccessMode_h__*/

