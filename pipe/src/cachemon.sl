#! /usr/bin/env slsh

require ("cmdopt");
require ("process");
require ("sysconf");
require ("rand");
require ("pcre");

$1 = path_dirname (__FILE__);
prepend_to_slang_load_path (path_concat ($1, "../share/slsh/local-packages"));
prepend_to_slang_load_path ($1);
set_import_module_path (path_concat ($1, "../lib/slang/v2/modules") + ":" + get_import_module_path ());
require ("pipeutil");

private variable Num_Cpus = sysconf("_SC_NPROCESSORS_ONLN");
if (Num_Cpus < 1)
{
   throw ApplicationError, "Invalid number of cpus: _SC_NPROCESSORS_ONLN = $Num_Cpus"$;
}

private variable Service_Name = NULL;

private variable Monitor_State;
private variable Sigterm_Received;
private variable Verbose = 0;

private variable Host_Name = uname().nodename;
private variable My_Pid = getpid();

private variable Num_Running = 0;

private variable LOG_INFO = "INFO";
private variable LOG_ERR = "ERROR";

private define write_log (severity, level, msg)
{
   if (Verbose >= level)
     () = fprintf (stderr, "cachemon[%d] %s %s\n", My_Pid, severity, msg);
}

private define file_age (file)
{
   variable s = stat_file (file);
   if (s == NULL) return _NaN;
   return double(_time - s.st_mtime);
}

private define newest_first (files)
{
   variable ages = array_map (Double_Type, &file_age, files);
   ages = ages[wherenot (isnan(ages))];
   return files[array_sort(ages)];
}

private define oldest_first (files)
{
   return reverse (newest_first(files));
}

private define alpha_order (files)
{
   return files[array_sort(files)];
}

private define parse_utc_timestamp_in_filename (filename)
{
   % e.g. TEMPO_RAD_L0_V01_20190614T000743Z_S1013G06.nc
   %      TEMPO_DRK_L0_V01_20190613T013740Z.nc
   variable tok = strtok (path_basename (filename), "_");
   return tok[4];
}
private define oldest_utc_timestamp_first (files)
{
   variable timestamps = array_map (String_Type, &parse_utc_timestamp_in_filename, files);
   return files[array_sort (timestamps)];
}

define paths_matching_pattern (patc, paths)
{
   variable p, m = {};
   foreach p (paths)
     {
        variable tok = pcre_matches (patc, p);
        if ((tok != NULL) && (0 == strcmp (tok[0], path_basename (p))))
          list_append (m, p);
     }

   if (0 == length(m))
     return NULL;

   return list_to_array(m);
}

private define dir_monitor (obj, order)
{
   variable dirlist = glob (path_concat (obj.incoming_dir, obj.glob));
   if (dirlist == NULL)
     {
        obj.dirlist = NULL;
        return 0;
     }

   if (obj.patc != NULL)
     {
        dirlist = paths_matching_pattern (obj.patc, dirlist);
        if (dirlist == NULL)
          {
             obj.dirlist = NULL;
             return 0;
          }
     }
   obj.dirlist = dirlist;

   % Because other processes may have changed the directory state
   % since we last checked, we'll loop over the list until either
   % we successfully process a file or we run out of files.
   % In the latter case, the assumption is that all the files
   % were handled by other processes before we got to them.

   variable file_list;

   switch (order)
     {
      case "newest":
	file_list = newest_first (dirlist);
     }
     {
      case "oldest":
	file_list = oldest_first (dirlist);
     }
     {
      case "oldest_utc_timestamp":
	file_list = oldest_utc_timestamp_first (dirlist);
     }
     {
      case "alpha":
	file_list = alpha_order (dirlist);
     }
     {
	% default:
	throw ApplicationError, "Unsupported sort order: $order"$;
     }

   if (obj.process == NULL || length(file_list) == 0)
     return 0;

   variable file;
   foreach file (file_list)
     {
        if (Sigterm_Received)
          break;
        if (1 == obj.process (file))
          break;
     }

   return 0;
}

private define sleep_loop (dt)
{
   % In this context, sleep(dt) is somewhat unreliable, because
   % SIGCHLD keeps interrupting it.  This loop is a hack
   % to work around that.  If implemented in C, we'd need
   % something similar anyway.
   variable trem = dt;
   variable tend = _time() + trem;
   forever
     {
        sleep (trem);
        variable tnow = _time();
        if ((tnow >= tend) || (Sigterm_Received != 0))
          break;
        trem = tend - tnow;
     }
}

private define do_wait (obj)
{
   sleep_loop (obj.delay);
}

private variable Disable_File_Path = "$SDPC_PIPE_DIR/ctrl/disable"$;

private define disable_read_config_tokens (disable_conf_file)
{
   variable lines = NULL;
   variable fp = fopen (disable_conf_file, "r");
   if (fp != NULL)
     {
        lines = fgetslines (fp ; trim=3);
     }
   if (lines == NULL || length(lines) == 0)
     {
        write_log (LOG_INFO, 1, sprintf ("invalid config file: %s", disable_conf_file));
        return NULL;
     }

   lines = strtrim (lines, "\n");

   % Parse the non-empty config file.
   % Skip comment lines prefixed with '#'.
   % Valid config lines look like:
   %   <service-specifier>:<config-string>
   %
   % EXAMPLES:
   % 1) disable.conf file to facilitate a quick SDPC shutdown:
   %    #
   %    # Disable the level1a service upon arrival of the next G01 radiance
   %    # file, and disable all other cachemon services immediately:
   %    #
   %    default:disable
   %    level1a:match:default
   %
   % 2) disable.conf a file to facilitate recovery from a quick shutdown:
   %    #
   %    # Allow all services to run normally, but disable the level1a
   %    # service upon arrival of the next day's first radiance file.
   %    #
   %    default:continue
   %    level1a:time:20260307T024500Z:default
   %

   variable config_tokens = Assoc_Type[];
   variable s, tok;
   foreach s (lines)
       {
          if (strlen(s) == 0 or s[0] == '#') continue;
          tok = strtok (s, ":");
          if (length(tok) < 2) continue;
          config_tokens[tok[0]] = tok[[1:]];
       }

   % If no service-specific tokens are found, use the default:

   % Service_Name == NULL should never happen
   if (0 != assoc_key_exists (config_tokens, Service_Name))
     {
        tok = config_tokens[Service_Name];
     }
   else if (0 != assoc_key_exists (config_tokens, "default"))
     {
        tok = config_tokens["default"];
     }
   else
     {
        write_log (LOG_INFO, 1, sprintf ("invalid config file: %s", disable_conf_file));
        return NULL;
     }

   return tok;
}

private define disable_ctrl_file_check (obj)
{
   if (obj.disable_ctrl_file_seen)
     return 0;

   variable disable_file_path = Disable_File_Path;
   if (0 != access (disable_file_path, F_OK))
     return 0;

   % If the disable file exists, look in the same directory
   % for a "disable.conf" file containing service-specific
   % configuration parameters.
   variable disable_conf_file = disable_file_path + ".conf";
   variable st = stat_file (disable_conf_file);
   if ((st == NULL) || (st.st_size == 0))
     {
        Monitor_State = 0;
        write_log (LOG_INFO, 1, sprintf ("disable: file exists: %s (no config parameters)", disable_file_path));
     }
   else
     {
        variable tokens = disable_read_config_tokens (disable_conf_file);
        % On read error, initialization fails
        if (NULL == tokens)
          return 0;
        obj.disable_config_tokens = tokens;

        % This cachemon instance might be disabled unconditionally:
        if (tokens[0] == "disable")
          {
             Monitor_State = 0;
             write_log (LOG_INFO, 1, sprintf ("disable on default config string"));
          }
     }

   % initialization succeeded
   obj.disable_ctrl_file_seen = 1;

   return 1;
}

private define disable_by_filename_pattern (obj, basename)
{
   if (obj.disable_config_tokens == NULL)
     return 0;

   variable tok = obj.disable_config_tokens;
   variable num_tokens = length(tok);
   variable config_string = strjoin (tok, ":");

   % Valid token strings are:
   % match:<pattern>
   % time:<when>:<pattern>
   %
   % where <pattern> = default | <regex>
   %          <when> = time stamp string YYYYMMDDThhmmssZ

   % config file string examples:
   % match
   % match:default
   % match:TEMPO_RAD_L\d{1}_[^G]*G01
   % match:TEMPO_RAD_L\d{1}_V\d{2}_\d{8}T\d{6}Z_S\d{3}G01
   % time:20260215T150000Z
   % time:20260215T150000Z:default
   % time:20260215T150000Z:TEMPO_RAD_L1_V\d{2}_(\d{8}T\d{6}Z)_S\d{3}G\d{02}

   variable m, regex, tstamp;
   switch (tok[0])
     {
      case "continue":
        % do nothing
     }
     {
        % disable upon appearance of a matching filename
        % (e.g. radiance scan granule number, either L0 or L1)
      case "match":
          regex = (num_tokens > 1) ? tok[1] : "default";
          if (regex == "default")
            {
               regex = "TEMPO_RAD_L\d{1}_[^G]*G01"R;
            }
          m = pcre_matches (regex, basename);
          if (m != NULL)
             {
                write_log (LOG_INFO, 1, sprintf ("disable on filename match: %s", basename));
                return 1;
             }
     }
     {
        % disable upon filename timestamp comparison
      case "time":
          if (num_tokens < 2)
            {
               write_log (LOG_INFO, 1, sprintf ("invalid config string: %s", config_string));
               return 0;
            }
          tstamp = tok[1];
          m = pcre_matches ("\d{8}T\d{6}Z"R, tstamp);  % minimal validation
          if (m == NULL)
            {
               write_log (LOG_INFO, 1, sprintf ("invalid config string: %s", config_string));
               return 0;
            }
          regex = (num_tokens > 2) ? tok[2] : "default";
          if (regex == "default")
             {
                regex = "TEMPO_RAD_L\d{1}_V\d{2}_(\d{8}T\d{6}Z)_S\d{3}G\d{2}"R;
             }
          m = pcre_matches (regex, basename);
          if ((m == NULL) || (length(m) < 2))
             return 0;
          if (tstamp < m[1])
             {
               write_log (LOG_INFO, 1, sprintf ("disable on filename timestamp comparison: %s < %s", tstamp, m[1]));
               return 1;
             }
     }
     {
        % default
        write_log (LOG_INFO, 1, sprintf ("invalid config string: %s", config_string));
     }

   return 0;
}

private define new_dirmon (incoming_dir)
{
   variable host_pid = sprintf ("%s_%d", Host_Name, getpid());

   variable obj = struct
     {
        incoming_dir = incoming_dir,
        glob = qualifier ("glob", "*"),
        patc = qualifier ("pcre_pat", NULL),
        substr_match_patc = qualifier ("substr_match_pat", NULL),
        substr_match_dir = qualifier ("substr_match_dir", NULL),
        monitor = &dir_monitor,
        process = qualifier ("task", NULL),
        disable_ctrl_file_seen = 0,
        disable_config_tokens = NULL,
        disable_ctrl_file_check = &disable_ctrl_file_check,
        disable_by_filename_pattern = &disable_by_filename_pattern,
        wait = &do_wait,
        delay = 1.0,
        num_processed = 0,
        dirlist = NULL,
        client_data = qualifier ("client_data", NULL),
        host_pid = host_pid
     };

   if (obj.patc != NULL)
     {
        obj.patc = pcre_compile (obj.patc);
     }

   if (obj.substr_match_patc != NULL)
     {
        obj.substr_match_patc = pcre_compile (obj.substr_match_patc);
     }

   return struct_combine (obj, __qualifiers);
}

private define exit_msg (w)
{
   if (w.exited)
     return sprintf ("exited normally: exit_status %d", w.exit_status);

   if (w.signal)
     return sprintf ("terminated by signal %d %s", w.signal,
                     (w.coredump ? "(coredump)" : ""));
   else if (w.stopped)
     return sprintf ("stopped by signal %d", w.stopped);
   else if (w.continued)
     return "continued";
}

private define sigchld_handler (sig);
private define sigchld_handler (sig)
{
   forever
     {
        variable w = waitpid (-1, WNOHANG|WUNTRACED);
        if ((w == NULL) || (w.pid == 0))
          break;
	write_log (LOG_INFO, 1, sprintf ("[%d] %s", w.pid, exit_msg(w)));
        Num_Running--;
     }
   signal (SIGCHLD, &sigchld_handler);
}
private define catch_sigchild ()
{
   sigprocmask (SIG_BLOCK, SIGCHLD);
   signal (SIGCHLD, &sigchld_handler);
   sigprocmask (SIG_UNBLOCK, SIGCHLD);
}

private define sigterm_handler (sig);
private define sigterm_handler (sig)
{
   Sigterm_Received++;
   signal (SIGTERM, &sigterm_handler);
}
private define catch_sigterm ()
{
   sigprocmask (SIG_BLOCK, SIGTERM);
   signal (SIGTERM, &sigterm_handler);
   sigprocmask (SIG_UNBLOCK, SIGTERM);
}

private define sigusr1_handler (sig);
private define sigusr1_handler (sig)
{
   Monitor_State = 1;
   signal (SIGUSR1, &sigusr1_handler);
}
private define catch_sigusr1 ()
{
   sigprocmask (SIG_BLOCK, SIGUSR1);
   signal (SIGUSR1, &sigusr1_handler);
   sigprocmask (SIG_UNBLOCK, SIGUSR1);
}

private define sigusr2_handler (sig);
private define sigusr2_handler (sig)
{
   Monitor_State = 0;
   signal (SIGUSR2, &sigusr2_handler);
}
private define catch_sigusr2 ()
{
   sigprocmask (SIG_BLOCK, SIGUSR2);
   signal (SIGUSR2, &sigusr2_handler);
   sigprocmask (SIG_UNBLOCK, SIGUSR2);
}

private variable Monitor_Last_State;
private define monitoring_enabled ()
{
   if (Monitor_State != Monitor_Last_State)
     {
        variable msg = sprintf ("New monitor state = %s", Monitor_State ? "ENABLE" : "DISABLE");
        write_log (LOG_INFO, 0, msg);
        Monitor_Last_State = Monitor_State;
     }

   return Monitor_State;
}
private define init_monitor_enable_switch ()
{
   Monitor_State = 1; % non-zero means ENABLE
   Monitor_Last_State = Monitor_State;
   catch_sigusr1();
   catch_sigusr2();
}

private define wait_for_processes_to_exit ()
{
   if (Num_Running == 0)
     return;

   write_log (LOG_INFO, 0, "waiting for $Num_Running processes to exit"$);

   while (Num_Running > 0)
     {
        sleep(1);
     }

   write_log (LOG_INFO, 0, "done");
}

private variable Exec = NULL;
private define set_executable (p, argv)
{
   if (p.exec_name == NULL)
     return;

   variable s = new_process (["which", p.exec_name];
                             fp2=1, stdout="/dev/null").wait();
   if (s.exit_status != 0)
     {
        throw ApplicationError,
          sprintf("*** Error: executable '%s' not found", p.exec_name);
     }

   Exec = struct
     {
        name=p.exec_name,
        wait_arg = p.wait_arg,
        argv=argv
     };

   if (p.exec_root_dir == NULL)
     return;

   if (0 != mkdir_p (p.exec_root_dir))
     {
        variable msg = sprintf ("*** Error: creating directory: %s",
                                p.exec_root_dir);
        throw ApplicationError, msg;
     }
}

private define run_executable (obj, file, run_dir)
{
   if (Exec == NULL)
     return 0;

   if (file == NULL)
     throw ApplicationError;

   variable argv = [Exec.name, file];
   if (Exec.argv != NULL && length(Exec.argv) > 0)
     {
        argv = [argv, Exec.argv];
     }

   % Note that if $run_dir is NFS mounted, then redirecting stdout
   % to a file in $run_dir is probably a bad idea.  In that case,
   % it may not then be possible for the new process to clean up after
   % itself by deleting $run_dir before exiting -- because the OS
   % may create a .nfs file in $run_dir, and that file can't be
   % deleted until _after_ the process itself exits.
   % Since stdout is usually empty anyway, we'll usually just
   % redirect output to /dev/null, but provide a config file
   % alternative for debugging.

   variable p, dir_str = "";
   if (run_dir != NULL)
     {
        p = new_process (argv; dir=run_dir);
        dir_str = sprintf (" in %s", run_dir);
     }
   else
     {
        p = new_process (argv);
     }

   variable s;
   if (Exec.wait_arg != WNOHANG)
     s = p.wait();
   else s = p.wait(WNOHANG);
   if (s == NULL)
     throw OSError, "waitpid failed: " + errno_string ();

   variable argv_rest = "";
   if (length(argv) > 1)
     {
        argv_rest = strjoin (argv[[1:]], " ");
     }
   variable msg = sprintf ("[%d] started%s: %s %s",
                           p.pid, dir_str, path_basename(argv[0]),
                           argv_rest);
   write_log (LOG_INFO, 1, msg);

   if (s.pid == 0)
     {
        Num_Running += 1;
        return 0;
     }

   if (s.exited)
     {
        write_log (LOG_INFO, 1, sprintf ("child exited with status %d", s.exit_status));
     }
   else if (s.signal)
     {
	write_log (LOG_INFO, 1, sprintf ("child terminated by signal %d %s",
				      s.signal, (s.coredump ? "(coredump)" : "")));
     }

   return -1;
}

private define have_substr_match (patc, dir, basename)
{
   % extract the substring we'll search for:
   variable num = pcre_exec (patc, basename);
   if (num == 0)
     return 0;

   variable seek_substr = pcre_nth_substr (patc, basename, 1);
   if (seek_substr == NULL)
     return 0;

   variable file_list = listdir (dir);
   if (length (file_list) == 0)
     return 0;

   variable file;
   foreach file (file_list)
     {
        if (is_substr (file, seek_substr))
          return 1;
     }

   return 0;
}

private define claim_with_rename (obj, path)
{
   if (path == NULL)
     return -1;

   variable cl = obj.client_data;

   if (cl.use_cpu_limiter != 0)
     {
        % If the machine is fully loaded, try again later
        if (Num_Running > Num_Cpus)
          return 0;
     }

   variable bname = path_basename (path);
   variable dname = path_dirname (path);
   variable claimed_bname = ".${bname}"$;
   variable claimed_path = path_concat (dname, claimed_bname);

   if (0 != obj.disable_by_filename_pattern (bname))
     {
        Monitor_State = 0;
        return 1;  % return non-zero to stop file processing
     }

   if (cl.substr_match_dir != NULL)
     {
        % Run executable only when a file with a matching substring exists in substr_match_dir
        if (have_substr_match (obj.substr_match_patc, obj.substr_match_dir, bname) == 0)
          return 0;
     }

   if (rename (path, claimed_path) != 0)
     return 0;

   if (run_executable (obj, claimed_path, cl.exec_root_dir) < 0)
     return -1;

   obj.num_processed++;

   return 1;
}

private define usage ()
{
   variable msg =
`Usage: cachemon.sl [opts] <config-file> [-- <exec-args>]
Options:
    -r|--rename      Claim files using only rename
`;
   () = fprintf (stderr, msg);
   exit (0);
}

private define cmdopt_error (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

% Because this variable is used in the configuration file,
% it must not be declared static or private.
variable _P = struct
{
   order = "alpha",
   incoming_dir = NULL,
   file_glob = "*",
   file_pcre_pat = NULL,
   substr_match_pat = NULL,
   substr_match_dir = NULL,
   wait_sec = 1.0,
   wait_arg = WNOHANG,
   use_cpu_limiter = 1,     % boolean
   exec_root_dir = NULL,
   exec_name = NULL,
   verbose = Verbose
};

private define load_config_file (file)
{
   if (path_is_absolute(file))
     () = evalfile (file);
   else
     () = evalfile (path_concat (".", file));
   return _P;
}

define slsh_main()
{
   variable exec_args = NULL;

   variable opts = cmdopt_new (&cmdopt_error);
   variable i = opts.process (__argv,1);

   if (__argc == i || i < 0)
     usage();

   variable config_file = __argv[i];
   if (NULL == stat_file (config_file))
     {
        () = fprintf (stderr, "*** Error: cannot read config file: $config_file"$);
        exit(1);
     }

   Service_Name = path_basename_sans_extname (config_file);

   if (__argc > i+1 && __argv[i+1] == "--")
     {
        exec_args = __argv[[i+2:]];
     }

   variable p = load_config_file (config_file);

   Verbose = p.verbose;

   set_executable (p, exec_args);

   variable dir = new_dirmon (p.incoming_dir; glob = p.file_glob,
                              pcre_pat = p.file_pcre_pat,
                              substr_match_pat = p.substr_match_pat,
                              substr_match_dir = p.substr_match_dir,
                              task = &claim_with_rename,
                              delay = p.wait_sec, client_data = p);

   Sigterm_Received = 0;
   catch_sigterm();
   catch_sigchild();

   init_monitor_enable_switch();

   write_log (LOG_INFO, 0, "started");

   while (Sigterm_Received == 0)
     {
        () = dir.disable_ctrl_file_check();

        if (monitoring_enabled())
          {
             if (-1 == dir.monitor (p.order))
               break;
          }
        dir.wait();
     }

   if (Sigterm_Received)
     write_log (LOG_INFO, 0, "received SIGTERM");

   wait_for_processes_to_exit ();
}

