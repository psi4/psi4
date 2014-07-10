/** @file
 * Define EVENT and KEY values used when calling evlog.
 */
#ifndef EVLOG_H_
#define EVLOG_H_

extern void evlog(int farg_key, ...);

/* Values of keys in key value pairs */

#define EVKEY_LAST_ARG  0    /**> Terminates list ... takes no value */

#define EVKEY_BEGIN     1    /**> Push (char *) value onto state stack */
#define EVKEY_END       2    /**> Pop  (char *) value off state stack  */
#define EVKEY_EVENT     3    /**> Record (char *) value, no stack change */

#define EVKEY_MSG_LEN   4    /**> Value is (int) mesage length SND/RCV only */
#define EVKEY_MSG_TO    5    /**> Value is (int) to process id SND/RCV only */
#define EVKEY_MSG_FROM  6    /**> Value is (int) from process  SND/RCV only */
#define EVKEY_MSG_TYPE  7    /**> Value is (int) message type  SND/RCV only */
#define EVKEY_MSG_SYNC  8    /**> Value is (int) message sync  SND/RCV only */

#define EVKEY_STR_INT   9    /**> User data value pair (char *), (int) */
#define EVKEY_STR_DBL  10    /**> User data value pair (char *), (double) */
#define EVKEY_STR      11    /**> User data value (char *) */

#define EVKEY_ENABLE   12    /**> Enable logging  ... takes no value */
#define EVKEY_DISABLE  13    /**> Disable logging ... takes no value */

#define EVKEY_DUMP     14    /**> Dump out the current buffer to disk */

#define EVKEY_FILENAME 15    /**> Set the name of the events file */

#define EVENT_SND      "Snd" /**> Predefined strings for internal events */
#define EVENT_RCV      "Rcv"   
#define EVENT_PROCESS  "Process"

#endif /* EVLOG_H_ */
