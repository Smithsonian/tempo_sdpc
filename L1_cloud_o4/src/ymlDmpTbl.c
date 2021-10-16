/*------------------------------------------------------------------------------*\
** ymlDmpTbl.c
**
** Dump table with parsed YAML data.
**
**   LDx = level of detail
**   nry = number of elements in tbl[]
**   tbl = table of YAML data
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlDmpTbl(int LDx, int nry, yTB tbl[])
{
    /*--------------------------------------------------------------------------*/

    char* clb[] = {					/* column labels	*/
			"tbl",
			"lvl",
			"blk",
			"usd",
			"key",
			"val",
                  };

    size_t nrc = sizeof(clb)/sizeof(clb[0]);		/* number of columns	*/

    /*--------------------------------------------------------------------------*/

    size_t wid[nrc];					/* column widths	*/

    int ix1; for( ix1 = 0; ix1 < nrc; ix1++ ) { wid[ix1] = strlen(clb[ix1]); }

    /*--------------------------------------------------------------------------*/

    int ix2; for( ix2 = 0; ix2 < nry; ix2++ )		/* override some widths	*/
    {
        char* key = tbl[ix2].key;			/* convenience variable	*/
        char* val = tbl[ix2].val;			/* convenience variable	*/

        if ( key[0] == '\0' ) { break; }		/* done; end of table	*/

        size_t szK = strlen(key);			/* get string length	*/
        size_t szV = strlen(val);			/* get string length	*/

        if ( wid[4] < szK ) { wid[4] = szK; }		/* update width for key	*/
        if ( wid[5] < szV ) { wid[5] = szV; }		/* update width for val	*/
    }

    /*--------------------------------------------------------------------------*/
    /* create column formats							*/

    char b1[C_STR_LEN];					/* str buf for format 1	*/
    char b2[C_STR_LEN];					/* str buf for format 2	*/
    char b3[C_STR_LEN];					/* str buf for dashes	*/
    char b4[20];					/* short tmp str buffer	*/

    int  s1 = sizeof(b1);				/* size of b1		*/
    int  m1 = s1;					/* max remaining in b1	*/
    int  t1 = 0;					/* total written to b1	*/
    int  p1 = 0;					/* number chars printed	*/

    int  s2 = sizeof(b2);				/* size of b2		*/
    int  m2 = s2;					/* max remaining in b2	*/
    int  t2 = 0;					/* total written to b2	*/
    int  p2 = 0;					/* number chars printed	*/

    int  s3 = sizeof(b3);				/* size of b3		*/
    int  m3 = s3;					/* max remaining in b3	*/
    int  t3 = 0;					/* total written to b3	*/
    int  p3 = 0;					/* number chars printed	*/

    int  s4 = sizeof(b4);				/* size of b4		*/
    int  m4 = s4;					/* max remaining in b4	*/
    int  t4 = 0;					/* total written to b4	*/
    int  p4 = 0;					/* number chars printed	*/

    /*--------------------------------------------------------------------------*/

    char oks[] = "OK";					/* string fragment	*/

    p1 = snprintf(b1+t1, m1, oks);			/* add string fragment	*/
    p2 = snprintf(b2+t2, m2, oks);			/* add string fragment	*/
    p3 = snprintf(b3+t3, m3, oks);			/* add string fragment	*/

    if ( p1  <  0 ) { p1 = m1; }			/* pre-glibc-2.0.6	*/
    if ( p2  <  0 ) { p2 = m2; }			/* pre-glibc-2.0.6	*/
    if ( p3  <  0 ) { p3 = m3; }			/* pre-glibc-2.0.6	*/

    if ( p1 >= m1 ) { p1 = m1 - 1; } t1 += p1;		/* update total chars	*/
    if ( p2 >= m2 ) { p2 = m2 - 1; } t2 += p2;		/* update total chars	*/
    if ( p3 >= m3 ) { p3 = m3 - 1; } t3 += p3;		/* update total chars	*/

    b1[t1] = '\0';					/* terminate string	*/
    b2[t2] = '\0';					/* terminate string	*/
    b3[t3] = '\0';					/* terminate string	*/

    m1 = s1 - p1;					/* new remaining in b1	*/
    m2 = s2 - p2;					/* new remaining in b2	*/
    m3 = s3 - p3;					/* new remaining in b3	*/

    char f4[] = "  %%%ds";				/* fmt fragment: "%7s"	*/
    char f2[] = "  %%%d%c";				/* fmt fragment: "%3d"	*/

    int ix3; for( ix3 = 0; ix3 < nrc; ix3++ )		/* add widths to format	*/
    {
        int wd = wid[ix3];				/* one column ix3 width	*/

        p4 = snprintf(b4, m4, f4, wd);			/* add width number	*/

        if ( p4  <  0 ) { p4 = m4; }			/* pre-glibc-2.0.6	*/
        if ( p4 >= m4 ) { p4 = m4 - 1; } t4 = p4;	/* update total chars	*/

        b4[t4] = '\0';					/* terminate string	*/

        /*----------------------------------------------------------------------*/

        int ty = 'd';					/* field type (d or s)	*/

        if ( ix3 == 4 ) { wd = -wid[ix3]; ty = 's'; }	/* left justified: key	*/
        if ( ix3 == 5 ) { wd = -wid[ix3]; ty = 's'; }	/* left justified: val	*/

        p1 = snprintf(b1+t1, m1, b4, clb[ix3]);		/* add column label	*/

        if ( p1  <  0 ) { p1 = m1; }			/* pre-glibc-2.0.6	*/
        if ( p1 >= m1 ) { p1 = m1 - 1; } t1 += p1;	/* update total chars	*/

        b1[t1] = '\0';		m1 = s1 - t1;		/* terminate; new limit	*/

        /*----------------------------------------------------------------------*/

        p2 = snprintf(b2+t2, m2, f2, wd, ty);		/* add column type	*/

        if ( p2  <  0 ) { p2 = m2; }			/* pre-glibc-2.0.6	*/
        if ( p2 >= m2 ) { p2 = m2 - 1; } t2 += p2;	/* update total chars	*/

        b2[t2] = '\0';		m2 = s2 - t2;		/* terminate; new limit	*/

        /*----------------------------------------------------------------------*/

        p3 = 2;						/* add this many chars	*/

        if ( p3 >= m3 ) { p3 = m3 - 1; }		/* prevents overflow	*/

        memset(b3+t3, ' ', p3); t3 += p3;		/* add spaces		*/

        b3[t3] = '\0';		m3 = s3 - t3;		/* terminate; new limit	*/

        /*----------------------------------------------------------------------*/

        p3 = wid[ix3];					/* add this many chars	*/

        if ( p3 >= m3 ) { p3 = m3 - 1; }		/* prevents overflow	*/

        memset(b3+t3, '-', p3); t3 += p3;		/* add dashes		*/

        b3[t3] = '\0';		m3 = s3 - t3;		/* terminate; new limit	*/

        /*----------------------------------------------------------------------*/
    }

    Info(LDx, b1);					/* show column labels	*/
    Info(LDx, b3);					/* show line of dashes	*/

    int ix4; for( ix4 = 0; ix4 < nry; ix4++ )		/* loop rows of tbl[]	*/
    {
        int   lvl = tbl[ix4].lvl;			/* convenience variable	*/
        int   blk = tbl[ix4].blk;			/* convenience variable	*/
        int   usd = tbl[ix4].usd;			/* convenience variable	*/
        char* key = tbl[ix4].key;			/* convenience variable	*/
        char* val = tbl[ix4].val;			/* convenience variable	*/

        Info(LDx, b2, ix4, lvl, blk, usd, key, val);	/* show one row of data	*/

        if ( key[0] == '\0' ) { break; }		/* done; end of table	*/
    }

    /*--------------------------------------------------------------------------*/

    return(0);
}

/*------------------------------------------------------------------------------*/
