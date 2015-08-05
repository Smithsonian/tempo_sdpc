/** @file
 *  @brief Private utility functions
 */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>

#include "netcdf.h"
#include "tell.h"
#include "tio.h"
#include "tio_template.h"
#include "_tio.h"

int _pTIO_define_processing_level (int grp, int level)
{
   int status, enum_typeid;
   static _pEnum_Type enum_table[] =
     {
        {"level-0",  TIO_PROC_LEVEL_0},
        {"level-1a", TIO_PROC_LEVEL_1A},
        {"level-1b", TIO_PROC_LEVEL_1B},
        {"level-2",  TIO_PROC_LEVEL_2},
        {"level-3",  TIO_PROC_LEVEL_3},
        _pENUM_TABLE_END
     };

   if (-1 == _pTIO_define_enum (grp, "processing_level_enum", enum_table, &enum_typeid))
     return -1;
   status = nc_put_att (grp, NC_GLOBAL,
                        "processing_level", enum_typeid, 1, &level);
   if (NC_NOERR != status)
     {
        Tell_verror (TELL_IO_WRITE_ERROR,
                     "%s: defining processing_level attribute (%s)",
                     __func__, nc_strerror(status));
        return -1;
     }

   return 0;
}
