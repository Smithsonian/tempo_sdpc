#ifndef __TEMPO_FILEDB_INCLUDE_H__
#define __TEMPO_FILEDB_INCLUDE_H__ 1

#include <stdlib.h>
#include <time.h>

#ifndef MALLOC
# define MALLOC malloc
#endif

#ifndef FREE
# define FREE free
#endif

typedef struct Filedb_Type Filedb_Type;
typedef struct Filedb_Entry_Type Filedb_Entry_Type;

struct Filedb_Type
{
   Filedb_Entry_Type *lst;
   const char *basename_pattern;
   char *root_dir;
   char *lookup_table;
   int (*parse_timestamp)(const char *, struct tm *);
};

extern int read_config_common (Filedb_Type *fdb, config_t *cfg, const char *name);

/* extern int config_ephemeris (Filedb_Type *fdb, config_t *cfg, const char *name); */
extern int config_snow (Filedb_Type *fdb, config_t *cfg, const char *name);
extern int config_tempo (Filedb_Type *fdb, config_t *cfg, const char *name);
extern int config_met (Filedb_Type *fdb, config_t *cfg, const char *name);

#endif
