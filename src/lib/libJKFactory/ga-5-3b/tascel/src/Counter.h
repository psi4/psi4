#ifndef __tascel_Counter_h__
#define __tascel_Counter_h__

namespace tascel {
  class Counter {
  private:
    int num;   /**< Current value of counter */
  public:
    
    /**
     * Constructor
     */
  Counter(): num(0) {}
    
    /**
     * Increment counter by argument
     *
     * @param[in] n Value to increment counter by
     */
    inline void inc(int n=1) { 
      num += n;
    }

    /**
     * Returns current count
     */
    inline int count() const {
      return num;
    }
  }; /* Counter */
}; /* tascel */

#endif /* __tascel_Counter_h__ */

