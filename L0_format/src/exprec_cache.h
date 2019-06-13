#ifndef __EXPREC_CACHE_METHODS__
#define __EXPREC_CACHE_METHODS__ 1

#include <iocsdpc.h>
#include <libconfig.h>

typedef struct Exprec_Cache_Method_Type Exprec_Cache_Method_Type;

/* Cache implements a FIFO queue of records of the form (file, file_index) */
struct Exprec_Cache_Method_Type
{
   int (*cache_erec)(Exprec_Cache_Method_Type *, const char *, size_t);
   /**< Append a record to the end of the queue */

   int (*cache_num_recs)(Exprec_Cache_Method_Type *, size_t *);
   /**< Return the numbers of records in the queue */

   int (*cache_open)(Exprec_Cache_Method_Type *);
   /**< Prepare to process the queue */

   int (*cache_close)(Exprec_Cache_Method_Type *);
   /**< Finish processing the queue */

   IOCSDPC_Exprec_Type *(*cache_erec_get)(Exprec_Cache_Method_Type *);
   /**< Retrieve the record at the front of the queue (the oldest one) */

   int (*cache_erec_path)(Exprec_Cache_Method_Type *, char *, size_t);
   /**< Copy front (oldest) record's path into a buffer.
    *   @return len=-1 on error, len<buflen on success,
    *           len>=buflen indicates path truncated
    */

   int (*cache_erec_bad)(Exprec_Cache_Method_Type *);
   /**< Regard the front record as "bad", and advance to the next one */

   int (*cache_erec_done)(Exprec_Cache_Method_Type *);
   /**< Regard the front record as "done", and advance to the next one */

   void (*cache_delete)(Exprec_Cache_Method_Type *);
   /**< Delete this Exprec_Cache_Method_Type object */

#ifdef EXPREC_CACHE_METHOD_PRIVATE_DATA
   EXPREC_CACHE_METHOD_PRIVATE_DATA
#endif
};

extern Exprec_Cache_Method_Type *open_erec_cache_disk (config_t *cfg);
extern Exprec_Cache_Method_Type *open_erec_cache_mem (config_t *cfg);

#endif
