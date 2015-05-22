#!/usr/bin/env slsh

private variable Script_Version_String = "0.2.0";

require ("cmdopt");
require ("setfuns");

private variable Dot = struct
{
   caller_callee_entries,
   cluster_definitions,
   re_caller_exclude_list,
   re_callee_exclude_list,
};
private variable Dot_Mode = 0;

private define dot_init ()
{
   Dot_Mode = 1;
   Dot.cluster_definitions = {};
   Dot.caller_callee_entries = {};

   variable a = {};
   Dot.re_callee_exclude_list = a;
   list_append (a, "^err_");
   list_append (a, "^error_check");
   list_append (a, "^ezspline");
   list_append (a, "^interpolat");
   list_append (a, "^l1bread");
   list_append (a, "^array");
   list_append (a, "^round");
   list_append (a, "^convert");
   list_append (a, "^h5");
   list_append (a, "^he5");
   list_append (a, "^dpolft");

   a = {};
   Dot.re_caller_exclude_list = a;
   list_append (a, "^elsunc");
}

private define dot_add_caller_callee (caller, callee)
{
   variable re;
   foreach re (Dot.re_caller_exclude_list)
     {
	if (string_match (caller, re))
	  return 0;		       %  do not decend
     }
   foreach re (Dot.re_callee_exclude_list)
     {
	if (string_match (callee, re))
	  return 0;
     }
   list_append (Dot.caller_callee_entries, "$caller -> $callee"$);
   return 1;
}

private define write_dot_output (fp)
{
   variable cc_entries = list_to_array (Dot.caller_callee_entries);
   cc_entries = cc_entries[unique(cc_entries)];

   variable fontsize=16;
   variable header =
`        rankdir="LR"; splines="line"; clusterrank=global;
        #ratio=fill; nodesep=0.25; ranksep=0.2;
    edge [fontname="Helvetica",fontsize="$fontsize",labelfontname="Helvetica",labelfontsize="$fontsize"];
    node [fontname="Helvetica",fontsize="$fontsize",shape=record];
`$;
   variable indent = "    ";
   variable c_c   = strjoin (indent + cc_entries, ";\n");
   variable cldef = strjoin (indent + list_to_array (Dot.cluster_definitions), "\n\n");
   cldef="";
   () = fprintf (fp, "digraph \"G\"\n{\n%s%s;\n\n%s\n}\n", header, c_c, cldef);
}


% The general structure of a fortran file is:
%     module
%       use foo
%         .
%     subroutine
%       use bar
%       call bam
%     end subroutine
%     end module
%     program name
%       .
%     end program

% There may be multiple modules per file.  Here we assume just one.
private define extract_fields (line)
{
   line = strreplace (line, "!", " !");
   return strtok (strtrans (line, "^0-9A-Za-z_!'\"", " "));
}

private define get_file_deps (file)
{
   variable fp = fopen (file, "r");
   if (fp == NULL)
     {
	() = fprintf (stderr, "**WARNING: Could not open %s, skipped\n", file);
	return NULL;
     }
   variable lines = fgetslines (fp);
   lines = strcompress (lines, " \t\n");
   lines = strtrim (lines);
   lines = strlow (lines[where (strncmp (lines, "!", 1))]);

   variable ok = Char_Type[length(lines)];
   variable use_flag = 1, end_flag = 2, module_flag = 3,
     subroutine_flag = 4, call_flag = 5, interface_flag = 6,
     call_func_flag = 7,
     recursive_flag = 0x80, function_flag = 0x40;

   ok[wherenot (strncmp (lines, "use ", 4))] = use_flag;
   ok[wherenot (strncmp (lines, "end ", 4))] = end_flag;
   ok[wherenot (strncmp (lines, "module ", 7))] = module_flag;
   ok[wherenot (strncmp (lines, "subroutine ", 11))] = subroutine_flag;
   ok[wherenot (strncmp (lines, "recursive subroutine ", 21))] = subroutine_flag|recursive_flag;
   ok[wherenot (strncmp (lines, "call ", 5))] = call_flag;
   ok[wherenot (strncmp (lines, "program ", 8))] = subroutine_flag;
   ok[wherenot (strncmp (lines, "interface", 9))] = interface_flag;

   ok[where ((ok == 0) and is_substr (lines, "function "))] = function_flag;
   ok[where ((ok == 0)
	     and (0 == strncmp (lines, "if", 2))
	     and is_substr (lines, "call "))]
     = call_flag;

   variable fun_call_re = "^[a-zA-Z0-9_]+ ?= ?\\([a-zA-Z0-9_]+\\) ?(";
   variable fun_call_names = String_Type[length(ok)];
   variable i;
   _for i (0, length(ok)-1, 1)
     {
	if (ok[i]) continue;
	variable matches = string_matches (lines[i], fun_call_re);
	if (matches == NULL) continue;
	ok[i] = call_func_flag;
	fun_call_names[i] = matches[1];
     }

   i = where (ok);
   lines = lines[i];
   fun_call_names = fun_call_names[i];
   %print (lines);
   ok = ok[i];

   variable num_ok = length(ok);
   if (num_ok == 0)
     return NULL;

   variable in_module = 0, in_subroutine = 0, in_function = 0,
     curr_module = NULL, curr_subroutine = NULL, in_interface = 0;

   variable s = struct
     {
	module_name = NULL,
	use_list = {},
	subroutines = {},
	functions = {},
	call_lists = {},
	fun_call_lists = {},
     };

   variable call_list = {};
   variable fun_call_list = {};
   variable j, field;
   i = 0;
   while (i < num_ok)
     {
	variable fields = extract_fields (lines[i]),
	  ok_i = ok[i], nfields = length(fields);

	i++;

	if (in_interface)
	  {
	     if (ok_i != end_flag)
	       continue;
	     if ((nfields >= 2) && (fields[1] == "interface"))
	       in_interface = 0;
	     continue;
	  }

	if (ok_i & recursive_flag)
	  {
	     fields = fields[where (fields != "recursive")];
	     nfields = length(fields);
	     ok_i = ok_i & ~recursive_flag;
	  }

	if (ok_i == end_flag)
	  {
	     if (nfields < 2)
	       continue;

	     if ((fields[1] == "subroutine") || (fields[1] == "program"))
	       {
		  if (in_subroutine == 0)
		    {
		       () = fprintf (stderr, "%s: Unexpected END SUBROUTINE encountered\n", file);
		       continue;
		    }

		  call_list = call_list[unique(call_list)];
		  list_append (s.call_lists, call_list);
		  call_list = {};

		  fun_call_list = fun_call_list[unique(fun_call_list)];
		  list_append (s.fun_call_lists, fun_call_list);
		  fun_call_list = {};

		  in_subroutine = 0;
		  continue;
	       }

	     if (fields[1] == "function")
	       {
		  if (in_function == 0)
		    {
		       () = fprintf (stderr, "%s: Unexpected END FUNCTION encountered\n", file);
		       continue;
		    }
		  call_list = call_list[unique(call_list)];
		  list_append (s.call_lists, call_list);
		  call_list = {};

		  fun_call_list = fun_call_list[unique(fun_call_list)];
		  list_append (s.fun_call_lists, fun_call_list);
		  fun_call_list = {};

		  in_function = 0;
		  continue;
	       }

	     if (fields[1] == "module")
	       {
		  if (in_module == 0)
		    {
		       () = fprintf (stderr, "%s: unexpected END MODULE encountered\n", file);
		    }
		  in_module = 0;
		  continue;
	       }

	     if (in_interface && (fields[1] == "interface"))
	       {
		  in_interface = 0;
		  continue;
	       }

	     continue;
	  }

	if (in_module == -1)
	  continue;

	if (ok_i == call_flag)
	  {
	     if (fields[0] == "if")
	       {
		  j = wherefirst (fields == "call");
		  if (j == NULL)
		    continue;
		  fields = fields[[j:]];
		  nfields = length (fields);
	       }

	     if (nfields < 2)
	       {
		  () = fprintf (stderr, "%s: Continued CALL statement not processed\n", file);
		  continue;
	       }
	     if (in_subroutine + in_function == 0)
	       {
		  () = fprintf (stderr, "%s: %S is called from something not a subroutine/function\n",
				file, fields[1]);
		  continue;
	       }
	     list_append (call_list, fields[1]);
	     continue;
	  }

	if (ok_i == call_func_flag)
	  {
	     list_append (fun_call_list, fields[1]);
	     continue;
	  }

	if (ok_i == use_flag)
	  {
	     if (nfields == 1)
	       {
		  () = fprintf (stderr, "***WARNING: %s: Continued USE not handled\n", file);
		  continue;
	       }
	     list_append (s.use_list, fields[1]);
	     continue;
	  }

	if (ok_i == subroutine_flag)
	  {
	     if (in_subroutine)
	       () = fprintf (stderr, "%s: END SUBROUTINE statement missing for %S\n", file, curr_subroutine);

	     if (nfields < 2)
	       {
		  () = fprintf (stderr, "%s: continued subroutine statement not handled\n", file);
		  curr_subroutine = NULL;
	       }
	     else curr_subroutine = fields[1];
	     in_subroutine = 1;
	     list_append (s.subroutines, curr_subroutine);
	     call_list = {};
	     fun_call_list = {};
	     continue;
	  }

	if (ok_i == function_flag)
	  {
	     % This is tricky, since "function" may be preceeded by the
	     % return type, also "function" may be part of a string or comment.
	     j = 0;
	     ok_i = 0;
	     foreach field (fields)
	       {
		  j++;
		  if (field == "function")
		    {
		       if (j < nfields)
			 {
			    ok_i = function_flag;
			    curr_subroutine = fields[j];
			    break;
			 }
		       continue;
		    }
		  if (field != strtrans (field, "'!\"", ""))
		    break;
	       }
	     if (ok_i)
	       {
		  list_append (s.functions, curr_subroutine);
		  list_append (s.subroutines, curr_subroutine);
		  in_function = 1;
		  call_list = {};
		  fun_call_list = {};
	       }
	     continue;
	  }

	if (ok_i == module_flag)
	  {
	     if (nfields == 1)
	       {
		  () = fprintf (stderr, "***WARNING: %s: continuation of MODULE not handled\n", file);
		  continue;
	       }
	     if (s.module_name != NULL)
	       {
		  () = fprintf (stderr, "***WARNING: %s: multiple modules in a single file are not supported\n");
		  in_module = -1;
		  continue;
	       }

	     if (s.module_name == NULL)
	       s.module_name = fields[1];
	     in_module = 1;
	     continue;
	  }

	if (ok_i == interface_flag)
	  in_interface = 1;
     }

   if (length (s.use_list))
     s.use_list = s.use_list[unique(s.use_list)];
#iffalse
   _for i (0, length(s.subroutines)-1, 1)
     {
	() = fprintf (stdout, "SUBROUTINE: %S\n", s.subroutines[i]);
	print (s.fun_call_lists[i]);
     }
#endif
   return s;
}

private define expand_deps_1 (files, ifile, ifile_deps_map, is_processed, new_list, chain);
private define expand_deps_1 (files, ifile, ifile_deps_map, is_processed, new_list, chain)
{
   list_append (chain, ifile);

   is_processed[ifile] = -1;

   variable list = ifile_deps_map[ifile];
   variable list2, dep;
   foreach dep (list)
     {
	variable ok = is_processed[dep];
	if (ok == -1)
	  {
	     () = fprintf (stderr, "Dependency loop detected: %S\n", files[dep]);
	     () = fprintf (stderr, " -> %S\n", strjoin (files[list_to_array(chain)], " -> "));
	     is_processed[dep] = 1;
	     continue;
	  }

	list_append (new_list, dep);
	if (ok == 1)
	  list2 = ifile_deps_map[dep];
	else
	  {
	     list2 = {};
	     expand_deps_1 (files, dep, ifile_deps_map, is_processed, list2, chain);
	  }
	list_join (new_list, list2);
     }

   is_processed[ifile] = 1;
   () = list_pop (chain, -1);

   if (length (new_list))
     new_list = new_list[unique (new_list)];
   ifile_deps_map[ifile] = new_list;
}


private define expand_deps (files, ifile_deps_map)
{
   variable file, dep, list;
   variable i, num = length(ifile_deps_map);
   variable is_processed = Char_Type[num];

   _for i (0, num-1, 1)
     {
	if (is_processed[i])
	  continue;

	variable new_list = {};
	expand_deps_1 (files, i, ifile_deps_map, is_processed, new_list, {});
	%print (file_deps_map[file]);
     }
}

private variable Call_Tree_Fp = stdout;
private define recurs_generate_tree ();
private define recurs_generate_tree (files, ifile, ifile_deps_map, indent, depth, max_depth,
				     rev_map, depth_map)
{
   variable ch = " +-";
   depth++;
   if (depth > max_depth) return;

   foreach (depth_map)
     {
	variable f = ();
	list_append (rev_map[ifile], f);
     }

   variable deps = ifile_deps_map[ifile];
   if (deps == NULL) return;

   variable file = files[ifile];
   variable i = 0, ndeps = length (deps);
   variable indent1 = indent + " | ";

   list_append (depth_map, ifile);

   while (i < ndeps)
     {
	variable idep = deps[i];
	variable dep = files[idep];
	i++;
	variable count = 1;
	while ((i < ndeps) && (idep == deps[i]))
	  {
	     count++;
	     i++;
	  }
	if (i == ndeps)
	  {
	     indent1 = indent + "   ";
	     ch = " `-";
	  }

	if (Dot_Mode)
	  {
	     if (0 == dot_add_caller_callee (file, dep))
	       continue;	       %  do not descend
	  }
	else
	  {
	     %() = fputs(indent + ch, Call_Tree_Fp);
	     if (count > 1)
	       () = fprintf (Call_Tree_Fp, "%s%s%s[%d] <--- %s\n", indent, ch, dep, count, file);
	     else
	       () = fprintf (Call_Tree_Fp, "%s%s%s <--- %s\n", indent, ch, dep, file);
	  }

	recurs_generate_tree (files, idep, ifile_deps_map, indent1, depth, max_depth,
			      rev_map, depth_map);
     }
   list_delete (depth_map, -1);
   %() = fputs ("\n", Call_Tree_Fp);
}

private define generate_tree (files, ifile, ifile_deps_map, indent, depth, max_depth, rev_map)
{
   variable file = files[ifile];
   if (Dot_Mode == 0)
     {
	() = fputs (file, Call_Tree_Fp);
	() = fputs ("\n", Call_Tree_Fp);
     }
   variable depth_map = {};

   recurs_generate_tree (files, ifile, ifile_deps_map, indent, depth, max_depth, rev_map, depth_map);
}


private define exit_version ()
{
   () = fprintf (stdout, "Version: %S\n", Script_Version_String);
   exit (0);
}

private define exit_usage ()
{
   () = fprintf (stderr, "Usage: %s [opts] file.f90...\n", __argv[0]);
   () = fprintf (stderr, "\
Options:\n\
  -d|--detect-loops     Check for circular dependencies\n\
  -h|--help             This usage message\n\
  -r|--rules            Generate makefile rules\n\
  --tree[=maxdepth]     Generate a file-dependency tree\n\
  -v|--version          Show version number\n\
"
		);
   exit (1);
}

define slsh_main ()
{
   variable do_rules = 0;
   variable do_loop_check = 0;
   variable max_depth = -1, max_ct_depth = -1;
   variable do_tree = 0;
   variable do_call_tree = 0;
   variable do_dot = 0;
   variable main_name = NULL;
   variable c = cmdopt_new ();
   c.add ("h|help", &exit_usage);
   c.add ("v|version", &exit_version);
   c.add ("call-tree", &max_ct_depth; type="int", optional=INT_MAX);
   c.add ("l|loop-detect", &do_loop_check);
   c.add ("r|rules", &do_rules);
   c.add ("d|dot", &do_dot; type="int", optional=INT_MAX);
   c.add ("m|main", &main_name; type="str");
   c.add ("tree", &max_depth; type="int", optional=INT_MAX);

   variable i = c.process (__argv, 1);
   if (i == __argc)
     exit_usage();

   if (max_depth != -1) do_tree = 1;
   if (max_ct_depth != -1) do_call_tree = 1;

   if (do_dot)
     {
	dot_init ();
	max_ct_depth = do_dot;
     }

   variable files = __argv[[i:]];

   variable numfiles = length (files);

   variable module_ifile_map = Assoc_Type[Int_Type, -1];
   variable ifile_module_map = String_Type[numfiles];
   variable ifile_use_modules_map = List_Type[numfiles];

   variable subroutine_int_map = Assoc_Type[Int_Type, -1];
   variable subroutine_call_list_map = Assoc_Type[List_Type];
   variable subroutine_file_map = Assoc_Type[String_Type];
   variable subroutines = {};

   variable function_call_list_map = Assoc_Type[List_Type];

   variable file, deps, module, dep, subroutine;

   variable j, num_subroutines = 0;

   _for i (0, numfiles-1, 1)
     {
	file = files[i];
	variable s = get_file_deps (file);
	deps = s.use_list;
	module = s.module_name;

        % define clusters for 'dot' format call graph
        if (do_dot && (length(s.subroutines) > 0))
          {
             list_append (Dot.cluster_definitions,
                          sprintf("subgraph cluster_%s {label=\"%s\"; %s; color=\"blue\"; }"$,
                                  path_basename_sans_extname(file),
                                  path_basename(file),
                                  strjoin (list_to_array(s.subroutines), ";")));
          }

	if (module != NULL)
	  {
	     module_ifile_map[module] = i;
	     ifile_module_map[i] = module;
	  }
	ifile_use_modules_map [i] = deps;

	j = 0;
	foreach subroutine (s.subroutines)
	  {
	     if (subroutine_int_map[subroutine] != -1)
	       {
		  () = fprintf (stderr, "Subroutine/Function %S in %S is duplicated in %S\n",
				subroutine, file, subroutine_file_map[subroutine]);
		  continue;
	       }
	     subroutine_file_map[subroutine] = file;

	     subroutine_int_map[subroutine] = num_subroutines;
	     list_append (subroutines, subroutine);
	     num_subroutines++;
	     try
	       {
		  subroutine_call_list_map[subroutine] = s.call_lists[j];
		  function_call_list_map[subroutine] = s.fun_call_lists[j];
	       }
	     catch AnyError: { print (s); throw;}
	     j++;
	  }
     }

   variable ifile_deps_map = List_Type[numfiles];

   _for i (0, numfiles-1, 1)
     {
	variable list = {};
	foreach module (ifile_use_modules_map[i])
	  {
	     j = module_ifile_map[module];
	     if (j != -1)
	       list_append (list, j);
	  }
	ifile_deps_map[i] = list;
     }

   % Now add the subroutines that were called, but have no corresponding
   % definition.
   variable ilist, name;
   foreach subroutine, list (subroutine_call_list_map)
     {
	foreach name (list)
	  {
	     j = subroutine_int_map[name];
	     if (j == -1)
	       {
		  subroutine_int_map[name] = num_subroutines;
		  num_subroutines++;
		  list_append (subroutines, name);
	       }
	  }
     }

   % Do the same for functions, but only those that are defined.
   foreach subroutine, list (function_call_list_map)
     {
	foreach name (list)
	  {
	     j = subroutine_int_map[name];
	     if (j != -1)
	       {
		  list_append (subroutine_call_list_map[subroutine], name);
	       }
	  }
     }


   variable isub_call_map = List_Type[num_subroutines];

   foreach subroutine, list (subroutine_call_list_map)
     {
	ilist = {};
	foreach name (list)
	  {
	     j = subroutine_int_map[name];
	     %if (j == -1) continue;
	     list_append (ilist, j);
	  }
	isub_call_map[subroutine_int_map[subroutine]] = ilist;
     }

   if (do_loop_check)
     expand_deps (files, ifile_deps_map);

   variable odir = "$(ODIR)";
   variable sdir = "$(SDIR)";
   variable compile_rule = "$(FC) -c $(FCFLAGS) $(IFLAGS) $sfile -o $ofile $(MODOUTFLG) $(ODIR)";

   if (do_rules) _for i (0, numfiles-1, 1)
     {
	file = files[i];
	variable ofile = path_sans_extname (file) + ".o";
	ofile = path_concat (odir, ofile);
	variable sfile = path_concat(sdir, file);
	() = fprintf (stdout, "%s: %s", ofile, sfile);
	deps = files[list_to_array(ifile_deps_map[i], Int_Type)];
	foreach file (deps)
	  {
	     file = path_sans_extname (file) + ".o";
	     () = fprintf (stdout, " \\\n  %s", path_concat (odir, file));
	  }
	() = fputs ("\n", stdout);
	() = fputs ("\t" + _$(compile_rule) + "\n", stdout);
     }

   variable reverse_deps_map;

   if (do_call_tree || do_dot)
     {
	reverse_deps_map = List_Type[num_subroutines];
	_for i (0, num_subroutines-1, 1)
	  reverse_deps_map [i] = {};

	_for i (0, num_subroutines-1, 1)
	  {
	     if ((main_name != NULL) && (subroutines[i] != main_name))
	       continue;

	     generate_tree (subroutines, i, isub_call_map, "", 0, max_ct_depth, reverse_deps_map);
	     if (do_dot == 0)
	       () = fputs ("\n", stdout);
	  }
	if (do_dot == 0) _for i (0, num_subroutines-1, 1)
	  {
	     variable rev_list = reverse_deps_map[i];
	     () = fprintf (stdout, "%s :", subroutines[i]);
	     if (length (rev_list))
	       {
		  rev_list = rev_list[unique (rev_list)];
		  foreach j (rev_list)
		    {
		       () = fprintf (stdout, " %s", subroutines[j]);
		    }
	       }
	     () = fputs ("\n\n", stdout);
	  }

	if (do_dot)
	  write_dot_output (stdout);
     }

   if (do_tree)
     {
	reverse_deps_map = List_Type[numfiles];
	_for i (0, numfiles-1, 1)
	  reverse_deps_map [i] = {};

	_for i (0, numfiles-1, 1)
	  {
	     generate_tree (files, i, ifile_deps_map, "", 0, max_depth, reverse_deps_map);
	     () = fputs ("\n\n", stdout);
	  }
     }
}

