/*------------------------------------------------------------------------------*\
** ymlErrors.c
**
** Shows detail for error returned by YAML parser.
**
**   LDx = level of detail
**   parser = data structure with all the info from the YAML library
\*------------------------------------------------------------------------------*/

#include "ymlPkgInc.h"

/*------------------------------------------------------------------------------*/

int ymlErrors(int LDx, yaml_parser_t parser)
{
    switch ( parser.error )
    {
      case YAML_MEMORY_ERROR: Error("not enough memory for parsing");
           break;

      case YAML_READER_ERROR:
           if (parser.problem_value != -1)
           {
              Error("reader error: %s: #%X at %d",
                     parser.problem, parser.problem_value,
                     parser.problem_offset);
           }
           else
           {
              Error("reader error: %s at %d",
                     parser.problem, parser.problem_offset);
           }
           break;

      case YAML_SCANNER_ERROR:
           if (parser.context)
           {
              Error("scanner error: %s at line %d, column %d",
                     parser.context, parser.context_mark.line+1,
                     parser.context_mark.column+1);
              Error("               %s at line %d, column %d",
                     parser.problem, parser.problem_mark.line+1,
                     parser.problem_mark.column+1);
           }
           else
           {
              Error("scanner error: %s at line %d, column %d",
                     parser.problem, parser.problem_mark.line+1,
                     parser.problem_mark.column+1);
           }
           break;

      case YAML_PARSER_ERROR:
           if (parser.context)
           {
              Error("parser error: %s at line %d, column %d",
                     parser.context, parser.context_mark.line+1,
                     parser.context_mark.column+1);
              Error("              %s at line %d, column %d",
                     parser.problem, parser.problem_mark.line+1,
                     parser.problem_mark.column+1);
           }
           else
           {
              Error("parser error: %s at line %d, column %d",
                     parser.problem, parser.problem_mark.line+1,
                     parser.problem_mark.column+1);
           }
           break;

      default: Error("internal error that should never happen");
           break;
   }

   return(0);
}

/*------------------------------------------------------------------------------*/
