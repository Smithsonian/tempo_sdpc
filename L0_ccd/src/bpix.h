#ifndef __BPIX_INCLUDE__
#define __BPIX_INCLUDE__ 1

typedef unsigned short Badpix_Bitmap_Type;
#define BADPIX_BITMAP_TIO_TYPE TIO_USHORT

typedef struct
{
   Badpix_Bitmap_Type *bits;
   int num_rows;
   int num_cols;
}
Badpix_Map_Type;

extern void bpix_free (Badpix_Map_Type *);
extern Badpix_Map_Type *bpix_new (int num_rows, int num_cols);
extern Badpix_Map_Type *bpix_read (const char *file);

extern int bpix_logical_or (Badpix_Map_Type *a,
                            const Badpix_Bitmap_Type *bits);

typedef struct Badpix_Map_Occur_Type Badpix_Map_Occur_Type;

extern void bpix_occur_close (Badpix_Map_Occur_Type *ot);
extern Badpix_Map_Occur_Type *bpix_occur_open (int num_rows, int num_cols,
                                               Badpix_Bitmap_Type mask);

extern int bpix_occur_incr (Badpix_Map_Occur_Type *ot,
                            const Badpix_Bitmap_Type *bits);

extern int bpix_occur_set (const Badpix_Map_Occur_Type *ot,
                           int num_threshold, Badpix_Bitmap_Type mask,
                           Badpix_Bitmap_Type *bits);
#endif
