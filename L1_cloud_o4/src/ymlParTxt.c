/*------------------------------------------------------------------------------*\
** ymlParTxt.c
**
** Parses a YAML input stream, stores parameter names and values into tbl[].
**
**   LDx = level of detail
**   toe = 'T' (parse by TOKEN) or 'E' (parse by EVENT)
**   raw = pointer to memory buffer with raw YAML text to parse
**   nrp = pointer to number of elements in tbl[]
**   ptr = pointer to table for YAML data
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlParTxt(int LDx, int toe, char* raw, int* nrp, yTB** ptr)
{
    int LDb = L7;				/* level of detail		*/

    /*--------------------------------------------------------------------------*/
    /* show version of the YAML library being used				*/

    const char* yst = yaml_get_version_string();

    Info(LDx, "OK LibYAML version: <%s>", yst);

    /*--------------------------------------------------------------------------*/
    /* initialize the YAML parser						*/

    if ( ymlChkTxt(LDx, raw) != 0 ) return(-7);		/* find junk, if any	*/

    yaml_parser_t parser;				/* parser internal info	*/

    if ( ! yaml_parser_initialize(&parser) )		/* there are errors	*/
    {
       ymlErrors(LDx, parser);				/* display error detail	*/
       return(-8);
    }

    yaml_parser_set_input_string(&parser, raw, strlen(raw));

    Info(LDx, "OK init YAML parser");

    /*--------------------------------------------------------------------------*/

    ymlIniTbl(LDx, nrp, ptr);				/* tbl[] initialization	*/

    int  nry = *nrp;				/* number of elements in tbl[]	*/
    yTB* tbl = *ptr;				/* tbl[] is what ptr points to	*/

    /*--------------------------------------------------------------------------*/

    Info(LDx, "OK starting to parse YAML");

    /*	There are two methods that can be used for parsing YAML files.  One
	method is based on TOKENs while the other is based on EVENTs.  I used
	the token-based method at first, but it got convoluted as I needed to
	parse ever-more complex YAML files.  The event-based method seems
	cleaner.  Still, at this time, both do the same job.			*/

    if ( toe != 'T' && toe != 'E' )
    {
       Error("parsing mode parameter <%c> must be either <T> or <E>");
       return(-15);
    }

    if ( toe == 'T' )
    {
       int er1 = ParUseTok(LDx, LDb, nry, tbl, parser); ReturnIfFatal(er1);
    }
    else
    {
       int er2 = ParUseEvt(LDx, LDb, nry, tbl, parser); ReturnIfFatal(er2);
    }

    Info(LDx, "OK finished parsing YAML");

    /*--------------------------------------------------------------------------*/

    yaml_parser_delete(&parser);			/* close up shop	*/

    /*--------------------------------------------------------------------------*/
    /* analyze if tbl[] was allocated properly					*/

    int idx; for( idx = 0; idx < nry; idx++ )		/* loop rows of tbl[]	*/
    {
        if ( tbl[idx].key[0] == '\0' ) { break; }	/* done; end of table	*/
    }

    if ( idx < nry ) { idx++; }				/* because zero-based	*/

    Info(LDx, "%d entries are used in tbl[%d]", idx, nry);

    /*--------------------------------------------------------------------------*/

    return(0);
}

/*------------------------------------------------------------------------------*/

int ShowToken(int LDx,yaml_token_type_t ytp,char* ynm) { Show_TorE(LDx,ytp,ynm); }
int ShowEvent(int LDx,yaml_event_type_t ytp,char* ynm) { Show_TorE(LDx,ytp,ynm); }

int Show_TorE(int LDx,yaml_event_type_t ytp,char* ynm)
{
    Info(LDx, "OK %02d = %s", ytp, ynm);
}

/*------------------------------------------------------------------------------*/

int ShowFlags(int LDx, int lvl, int blk, int nxK, int nxV, int nxA)
{
    Info(LDx, "OK    lvl=%d   blk=%d   nxK=%d   nxV=%d   nxA=%d",
                     lvl,     blk,     nxK,     nxV,     nxA    );
}

/*------------------------------------------------------------------------------*/

int ResetPars(int LDx, int* lvl, int* blk, int* nxV, int* nxK, int* nxA)
{
    *lvl = -1;						/* YAML level number	*/
    *blk = -1;						/* YAML block number	*/
    *nxK =  0;						/* set when key is next	*/
    *nxV =  0;						/* set when val is next	*/
    *nxA =  0;						/* set when val array	*/
}

/*------------------------------------------------------------------------------*/

int ParUseTok(int LDx, int LDy, int nry, yTB tbl[], yaml_parser_t parser)
{
    int err = 0;					/* code from function	*/

    int lvl;						/* YAML level number	*/
    int blk;						/* YAML block number	*/
    int nxK;						/* set when key is next	*/
    int nxV;						/* set when val is next	*/
    int nxA;						/* set when val array	*/

    char* itm;						/* convenience variable	*/
    char* anc;						/* convenience variable	*/

    ResetPars(LDx, &lvl, &blk, &nxK, &nxV, &nxA);	/* reset flags		*/

    ShowFlags(LDx, lvl, blk, nxK, nxV, nxA);		/* show initial values	*/

    int max = YML_MAX_STK; char* stk[max];		/* parts accumulation	*/

    int cnt; for( cnt = 0; cnt < max; cnt++ ) { stk[cnt] = '\0'; }	/* init	*/

    yaml_token_t token;					/* YAML data structure	*/

    while ( token.type != YAML_STREAM_END_TOKEN )	/* loop YAML stream	*/
    {
       if ( ! yaml_parser_scan(&parser, &token) )	/* get next token	*/
       {
          Error("failed to parse YAML text");		/* there are errors	*/

          ymlErrors(LDx, parser);			/* display error detail	*/

          yaml_token_delete(&token);			/* close up shop	*/
          yaml_parser_delete(&parser);			/* close up shop	*/

          return(-9);
       }

       switch ( token.type )					/* act on token	*/
       {
       case YAML_NO_TOKEN:
            ShowToken(LDx, token.type, "YAML_NO_TOKEN");	/* do nothing	*/
            break;
       case YAML_VERSION_DIRECTIVE_TOKEN:		/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_VERSION_DIRECTIVE_TOKEN");
            Info(LDy, "OK    major     = %d", token.data.version_directive.major);
            Info(LDy, "OK    minor     = %d", token.data.version_directive.minor);
            break;
       case YAML_TAG_DIRECTIVE_TOKEN:			/*   not implemented	*/
            ShowToken(LDx, token.type, "YAML_TAG_DIRECTIVE_TOKEN");
            Info(LDy, "OK    handle    = %s", token.data.tag_directive.handle);
            Info(LDy, "OK    prefix    = %s", token.data.tag_directive.prefix);
            break;
       case YAML_STREAM_START_TOKEN:			/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_STREAM_START_TOKEN");
            Info(LDy, "OK    encoding  = %d", token.data.stream_start.encoding);
            Info(LDy, "OK    index     = %d", token.start_mark.index);
            Info(LDy, "OK    line      = %d", token.start_mark.line);
            Info(LDy, "OK    column    = %d", token.start_mark.column);
            break;
       case YAML_DOCUMENT_START_TOKEN:
            ShowToken(LDx, token.type, "YAML_DOCUMENT_START_TOKEN");
            Info(LDy, "OK    index     = %d", token.start_mark.index);
            Info(LDy, "OK    line      = %d", token.start_mark.line);
            Info(LDy, "OK    column    = %d", token.start_mark.column);
            ResetPars(LDx, &lvl, &blk, &nxK, &nxV, &nxA);	/* reset flags	*/
            break;
       case YAML_BLOCK_MAPPING_START_TOKEN:		/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_BLOCK_MAPPING_START_TOKEN");
            if ( nxV )			/* value expected; got new key instead	*/
            {				/* save off previous (old) key (no val)	*/
               err = ymlAddKey(LDx, nry, tbl, lvl, stk, max);	/* concatenate	*/
               err = ymlAddVal(LDx, nry, tbl, lvl, blk, "");	/* save NO val	*/
            }
            nxK = 1;				/* next scalar is  a  key	*/
            nxV = 0;				/* next scalar is not val	*/
            nxA = 0;				/* next scalar not part of arr	*/
            lvl++;				/* one level up			*/
            break;
       case YAML_BLOCK_SEQUENCE_START_TOKEN:
            ShowToken(LDx, token.type, "YAML_BLOCK_SEQUENCE_START_TOKEN");
            nxA = 1;				/* val array started		*/
            nxV = 1;				/* next scalar is  a  val	*/
            break;
       case YAML_FLOW_MAPPING_START_TOKEN:		/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_FLOW_MAPPING_START_TOKEN");
            break;
       case YAML_FLOW_SEQUENCE_START_TOKEN:
            ShowToken(LDx, token.type, "YAML_FLOW_SEQUENCE_START_TOKEN");
            nxA = 1;				/* val array started		*/
            nxV = 1;				/* next scalar is  a  val	*/
            break;
       case YAML_BLOCK_ENTRY_TOKEN:			/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_BLOCK_ENTRY_TOKEN");
            break;
       case YAML_FLOW_ENTRY_TOKEN:			/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_FLOW_ENTRY_TOKEN");
            break;
       case YAML_KEY_TOKEN:
            ShowToken(LDx, token.type, "YAML_KEY_TOKEN");
            nxK = 1;				/* next scalar is  a  key	*/
            nxV = 0;				/* next scalar is not val	*/
            nxA = 0;				/* next scalar not part of arr	*/
            break;
       case YAML_VALUE_TOKEN:				/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_VALUE_TOKEN");
            break;
       case YAML_ANCHOR_TOKEN:
            ShowToken(LDx, token.type, "YAML_ANCHOR_TOKEN");
            Info(LDy, "OK    value     = %s", token.data.anchor.value);
            nxK = 0;				/* next scalar is not key	*/
            nxV = 1;				/* next scalar is  a  val	*/
            nxA = 0;				/* next scalar not part of arr	*/
            break;
       case YAML_TAG_TOKEN:				/*   not implemented	*/
            ShowToken(LDx, token.type, "YAML_TAG_TOKEN");
            Info(LDy, "OK    handle     = %s", token.data.tag.handle);
            Info(LDy, "OK    suffix     = %s", token.data.tag.suffix);
            break;
       case YAML_ALIAS_TOKEN:
            ShowToken(LDx, token.type, "YAML_ALIAS_TOKEN");
            Info(LDy, "OK    value     = %s", token.data.alias.value);
            break;
       case YAML_SCALAR_TOKEN:
            ShowToken(LDx, token.type, "YAML_SCALAR_TOKEN");
            Info(LDy, "OK    value     = %s", token.data.scalar.value);
            Info(LDy, "OK    length    = %d", token.data.scalar.length);
            Info(LDy, "OK    style     = %d", token.data.scalar.style);
            itm = token.data.scalar.value;			/* item to save	*/
            if ( nxK )
            {
               blk++;
               err = ymlStkPsh(LDx, itm, lvl, stk, max);      ReturnIfFatal(err);
               nxK = 0;						/* one key done	*/
               nxV = 1;						/* val is next	*/
            }
            else
            {
               err = ymlAddKey(LDx, nry, tbl, lvl, stk, max); ReturnIfFatal(err);
               err = ymlAddVal(LDx, nry, tbl, lvl, blk, itm); ReturnIfFatal(err);
               nxV = 0;						/* one val done	*/
               if ( ! nxA ) { nxK = 1; }			/* key is next	*/
            }
            break;
       case YAML_FLOW_SEQUENCE_END_TOKEN:
            ShowToken(LDx, token.type, "YAML_FLOW_SEQUENCE_END_TOKEN");
            nxV = 0;				/* next scalar is not val	*/
            nxA = 0;				/* val array ended		*/
            break;
       case YAML_FLOW_MAPPING_END_TOKEN:		/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_FLOW_MAPPING_END_TOKEN");
            break;
       case YAML_BLOCK_END_TOKEN:
            ShowToken(LDx, token.type, "YAML_BLOCK_END_TOKEN");
            if ( lvl >= 0 )
            {
               ymlStkPsh(LDx, '\0', lvl, stk, max);		/* remove old	*/
            }
            lvl--;						/* level ends	*/
            if ( lvl < -1 ) { lvl++; }	/* ### kluge ### */
            nxK = 1; nxV = 0; nxA = 0;				/* key is next	*/
            break;
       case YAML_DOCUMENT_END_TOKEN:			/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_DOCUMENT_END_TOKEN");
            Info(LDy, "OK    index     = %d", token.end_mark.index);
            Info(LDy, "OK    line      = %d", token.end_mark.line);
            Info(LDy, "OK    column    = %d", token.end_mark.column);
            break;
       case YAML_STREAM_END_TOKEN:			/* ignoring this token	*/
            ShowToken(LDx, token.type, "YAML_STREAM_END_TOKEN");
            Info(LDy, "OK    index     = %d", token.end_mark.index);
            Info(LDy, "OK    line      = %d", token.end_mark.line);
            Info(LDy, "OK    column    = %d", token.end_mark.column);
            break;
       default:
            Error("YAML TOKEN not recognized; this should never happen.");
            err = -16;
            break;
       }

       ShowFlags(LDx, lvl, blk, nxK, nxV, nxA);		/* show final values	*/

       ReturnIfFatal(err);
    }

    if ( lvl != -1 )
    {
       Error("parsing YAML by TOKEN did not go well.");
       Info(0, "    Value of lvl should be -1, but it is %d instead.", lvl);
       err = -18;
    }

    yaml_token_delete(&token);				/* close up shop	*/

    return(err);
}

/*------------------------------------------------------------------------------*/

int ParUseEvt(int LDx, int LDy, int nry, yTB tbl[], yaml_parser_t parser)
{
    int err = 0;					/* code from function	*/

    int lvl;						/* YAML level number	*/
    int blk;						/* YAML block number	*/
    int nxK;						/* set when key is next	*/
    int nxV;						/* set when val is next	*/
    int nxA;						/* set when arr is next	*/

    char* itm;						/* item value to save	*/

    ResetPars(LDx, &lvl, &blk, &nxK, &nxV, &nxA);	/* reset flags		*/

    ShowFlags(LDx, lvl, blk, nxK, nxV, nxA);		/* show initial values	*/

    int max = YML_MAX_STK; char* stk[max];		/* parts accumulation	*/

    int cnt; for( cnt = 0; cnt < max; cnt++ ) { stk[cnt] = '\0'; }	/* init	*/

    yaml_event_t event;					/* YAML data structure	*/

    while ( event.type != YAML_STREAM_END_EVENT )	/* loop YAML stream	*/
    {
       if ( ! yaml_parser_parse(&parser, &event) )	/* get next event	*/
       {
          Error("failed to parse YAML text");		/* there are errors	*/

          ymlErrors(LDx, parser);			/* display error detail	*/

          yaml_event_delete(&event);			/* close up shop	*/
          yaml_parser_delete(&parser);			/* close up shop	*/

          return(-9);
       }

       switch ( event.type )					/* act on event	*/
       {
       case YAML_NO_EVENT:
            ShowEvent(LDx, event.type, "YAML_NO_EVENT");	/* do nothing	*/
            break;
       case YAML_STREAM_START_EVENT:			/* ignoring this token	*/
            ShowEvent(LDx, event.type, "YAML_STREAM_START_EVENT");
            Info(LDy, "OK    encoding  = %d", event.data.stream_start.encoding);
            Info(LDy, "OK    index     = %d", event.start_mark.index);
            Info(LDy, "OK    line      = %d", event.start_mark.line);
            Info(LDy, "OK    column    = %d", event.start_mark.column);
            break;
       case YAML_DOCUMENT_START_EVENT:
            ShowEvent(LDx, event.type, "YAML_DOCUMENT_START_EVENT");
            Info(LDy, "OK    implicit  = %d", event.data.document_start.implicit);
            Info(LDy, "OK    index     = %d", event.start_mark.index);
            Info(LDy, "OK    line      = %d", event.start_mark.line);
            Info(LDy, "OK    column    = %d", event.start_mark.column);
            ResetPars(LDx, &lvl, &blk, &nxK, &nxV, &nxA);	/* reset flags	*/
            break;
       case YAML_MAPPING_START_EVENT:
            ShowEvent(LDx, event.type, "YAML_MAPPING_START_EVENT");
            Info(LDy, "OK    anchor    = %s", event.data.mapping_start.anchor);
            Info(LDy, "OK    implicit  = %d", event.data.mapping_start.implicit);
            Info(LDy, "OK    style     = %d", event.data.mapping_start.style);
            Info(LDy, "OK    index     = %d", event.end_mark.index);
            Info(LDy, "OK    line      = %d", event.end_mark.line);
            Info(LDy, "OK    column    = %d", event.end_mark.column);
            itm = event.data.mapping_start.anchor;		/* item to save	*/
            if ( nxV )			/* value expected; got new key instead	*/
            {				/* save off previous (old) key (no val)	*/
               err = ymlAddKey(LDx, nry, tbl, lvl, stk, max);	/* concatenate	*/
					/* save the anchor or a blank string	*/
               if ( ! itm ) { itm = ""; }			/* blank value	*/
               err = ymlAddVal(LDx, nry, tbl, lvl, blk, itm);	/* save NO val	*/
            }
            nxK = 1;				/* next scalar is  a  key	*/
            nxV = 0;				/* next scalar is not val	*/
            nxA = 0;				/* next scalar not part of arr	*/
            lvl++;				/* one level up			*/
            break;
       case YAML_SEQUENCE_START_EVENT:
            ShowEvent(LDx, event.type, "YAML_SEQUENCE_START_EVENT");
            Info(LDy, "OK    anchor    = %s", event.data.sequence_start.anchor);
            Info(LDy, "OK    implicit  = %d", event.data.sequence_start.implicit);
            Info(LDy, "OK    style     = %d", event.data.sequence_start.style);
            Info(LDy, "OK    index     = %d", event.start_mark.index);
            Info(LDy, "OK    line      = %d", event.start_mark.line);
            Info(LDy, "OK    column    = %d", event.start_mark.column);
            nxK = 0; nxV = 1; nxA = 1;				/* val arr start*/
            break;
       case YAML_ALIAS_EVENT:
            ShowEvent(LDx, event.type, "YAML_ALIAS_EVENT");
            Info(LDy, "OK    anchor    = %s", event.data.alias.anchor);
            itm = event.data.alias.anchor;			/* item to save	*/
            if ( nxK )
            {
               blk++;
               err = ymlStkPsh(LDx, itm, lvl, stk, max);      ReturnIfFatal(err);
               nxK = 0;						/* one key done	*/
               nxV = 1;						/* val is next	*/
            }
            else
            {
               err = ymlAddKey(LDx, nry, tbl, lvl, stk, max); ReturnIfFatal(err);
               err = ymlAddVal(LDx, nry, tbl, lvl, blk, itm); ReturnIfFatal(err);
               nxV = 0;						/* one val done	*/
               if ( ! nxA ) { nxK = 1; }			/* key is next	*/
            }
            break;
       case YAML_SCALAR_EVENT:
            ShowEvent(LDx, event.type, "YAML_SCALAR_EVENT");
            Info(LDy, "OK    anchor    = %s", event.data.scalar.anchor);
            Info(LDy, "OK    value     = %s", event.data.scalar.value);
            Info(LDy, "OK    length    = %d", event.data.scalar.length);
            Info(LDy, "OK    style     = %d", event.data.scalar.style);
            itm = event.data.alias.anchor;			/* item to save	*/
            if ( ! itm ) { itm = event.data.scalar.value; }	/* item to save	*/
            if ( nxK )
            {
               blk++;
               err = ymlStkPsh(LDx, itm, lvl, stk, max);      ReturnIfFatal(err);
               nxK = 0;						/* one key done	*/
               nxV = 1;						/* val is next	*/
            }
            else
            {
               err = ymlAddKey(LDx, nry, tbl, lvl, stk, max); ReturnIfFatal(err);
               err = ymlAddVal(LDx, nry, tbl, lvl, blk, itm); ReturnIfFatal(err);
               nxV = 0;						/* one val done	*/
               if ( ! nxA ) { nxK = 1; }			/* key is next	*/
            }
            break;
       case YAML_SEQUENCE_END_EVENT:
            ShowEvent(LDx, event.type, "YAML_SEQUENCE_END_EVENT");
            Info(LDy, "OK    index     = %d", event.end_mark.index);
            Info(LDy, "OK    line      = %d", event.end_mark.line);
            Info(LDy, "OK    column    = %d", event.end_mark.column);
            nxK = 1; nxV = 0; nxA = 0;				/* val arr end	*/
            break;
       case YAML_MAPPING_END_EVENT:
            ShowEvent(LDx, event.type, "YAML_MAPPING_END_EVENT");
            Info(LDy, "OK    index     = %d", event.end_mark.index);
            Info(LDy, "OK    line      = %d", event.end_mark.line);
            Info(LDy, "OK    column    = %d", event.end_mark.column);
            if ( lvl >= 0 )
            {
               ymlStkPsh(LDx, '\0', lvl, stk, max);		/* remove old	*/
            }
            lvl--;						/* level ends	*/
            nxK = 1; nxV = 0; nxA = 0;				/* key is next	*/
            break;
       case YAML_DOCUMENT_END_EVENT:			/* ignoring this token	*/
            ShowEvent(LDx, event.type, "YAML_DOCUMENT_END_EVENT");
            Info(LDy, "OK    implicit  = %d", event.data.document_end.implicit);
            Info(LDy, "OK    index     = %d", event.end_mark.index);
            Info(LDy, "OK    line      = %d", event.end_mark.line);
            Info(LDy, "OK    column    = %d", event.end_mark.column);
            break;
       case YAML_STREAM_END_EVENT:			/* ignoring this token	*/
            ShowEvent(LDx, event.type, "YAML_STREAM_END_EVENT");
            Info(LDy, "OK    index     = %d", event.end_mark.index);
            Info(LDy, "OK    line      = %d", event.end_mark.line);
            Info(LDy, "OK    column    = %d", event.end_mark.column);
            break;
       default:
            Error("YAML EVENT not recognized; this should never happen.");
            err = -17;
            break;
       }

       ShowFlags(LDx, lvl, blk, nxK, nxV, nxA);		/* show final values	*/

       ReturnIfFatal(err);
    }

    if ( lvl != -1 )
    {
       Error("parsing YAML by EVENT did not go well.");
       Info(0, "    Value of lvl should be -1, but it is %d instead.", lvl);
       err = -18;
    }

    yaml_event_delete(&event);				/* close up shop	*/

    return(err);
}

/*------------------------------------------------------------------------------*/
