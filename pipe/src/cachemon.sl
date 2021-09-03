#! /usr/bin/env slsh

require ("cmdopt");
require ("process");
require ("sysconf");
require ("rand");
%require ("pcre");

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

private variable Sigterm_Received;
private variable Verbose = 0;

private variable Host_Name = uname().nodename;
private variable My_Pid = getpid();

private variable Num_Running = 0;

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

private define dir_monitor (obj, order)
{
   variable dirlist = glob (path_concat (obj.incoming_dir, obj.glob));
   if (dirlist == NULL)
     {
        obj.dirlist = NULL;
        return 0;
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

private define new_dirmon (incoming_dir)
{
   variable host_pid = sprintf ("%s_%d", Host_Name, getpid());

   variable obj = struct
     {
        incoming_dir = incoming_dir,
        glob = qualifier ("glob", "*"),
        monitor = &dir_monitor,
        process = qualifier ("task", NULL),
        wait = &do_wait,
        delay = 1.0,
        num_processed = 0,
        dirlist = NULL,
        client_data = qualifier ("client_data", NULL),
        host_pid = host_pid
     };

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

private variable LOG_INFO = "INFO";
private variable LOG_ERR = "ERROR";

private define write_log (severity, level, msg)
{
   if (Verbose >= level)
     () = fprintf (stderr, "cachemon[%d] %s %s\n", My_Pid, severity, msg);
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

   if (__argc > i+1 && __argv[i+1] == "--")
     {
        exec_args = __argv[[i+2:]];
     }

   variable p = load_config_file (config_file);

   Verbose = p.verbose;

   set_executable (p, exec_args);

   variable dir = new_dirmon (p.incoming_dir; glob = p.file_glob,
                              task = &claim_with_rename,
                              delay = p.wait_sec, client_data = p);

   Sigterm_Received = 0;
   catch_sigterm();
   catch_sigchild();

   write_log (LOG_INFO, 0, "started");

   while (Sigterm_Received == 0)
     {
        if (-1 == dir.monitor (p.order))
          break;
        dir.wait();
     }

   if (Sigterm_Received)
     write_log (LOG_INFO, 0, "received SIGTERM");

   wait_for_processes_to_exit ();
}

