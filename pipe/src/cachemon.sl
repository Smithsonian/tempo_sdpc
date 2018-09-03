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
require ("daemon");

private variable Host_Name = uname().nodename;
private variable Num_Processors = sysconf("_SC_NPROCESSORS_ONLN");

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

private define dir_monitor (obj, newest)
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

   if (newest != 0)
     file_list = newest_first (dirlist);
   else
     file_list = oldest_first (dirlist);

   if (obj.process == NULL || length(file_list) == 0)
     return 0;

   variable file;
   foreach file (file_list)
     {
        if (1 == obj.process (file))
          break;
     }

   return 0;
}

private define sleep_loop (dt)
{
   % In this context, sleep(dt) is somewhat unreliable, because
   % (I think) SIGCHLD keeps interrupting it.  This loop is a hack
   % to work around that.  If implemented in C, we'd need
   % something similar anyway.
   variable trem = dt;
   variable tend = _time() + trem;
   forever
     {
        sleep (trem);
        variable tnow = _time();
        if (tnow >= tend)
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

private define sigchld_handler (sig);
private define sigchld_handler (sig)
{
   forever
     {
        variable w = waitpid (-1, WNOHANG|WUNTRACED);
        if ((w == NULL) || (w.pid == 0))
          break;
        daemon_log (LOG_INFO, sprintf ("[%d] %s", w.pid, exit_msg(w)));
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

private variable Sighup_Received;
private define sighup_handler (sig);
private define sighup_handler (sig)
{
   Sighup_Received++;
   signal (SIGHUP, &sighup_handler);
}
private define catch_sighup ()
{
   sigprocmask (SIG_BLOCK, SIGHUP);
   signal (SIGHUP, &sighup_handler);
   sigprocmask (SIG_UNBLOCK, SIGHUP);
}

private define wait_for_processes_to_exit ()
{
   if (Num_Running == 0)
     return;

   daemon_log (LOG_INFO, "waiting for $Num_Running processes to exit"$);

   while (Num_Running > 0)
     {
        sleep(1);
     }

   daemon_log (LOG_INFO, "done");
}

private define have_idle_processors (num_parallel)
{
   % /proc/loadavg gives a measure of the total system load averaged
   % over the last 1, 5, and 15 minutes
   variable load_avg_col;
   variable num_read = readascii ("/proc/loadavg", &load_avg_col;
                                  format="%f", nrows=1);
   variable sys_load = (num_read == 1) ? load_avg_col[0] : 0.0;

   % The system load average has an associated lag time and will essentially
   % never be controlled by the current process.  On the other hand,
   % the load that the current process has contributed is known exactly.
   variable self_load = Num_Running * num_parallel;

   % Try using a weighted average of these numbers to constrain the load:
   variable weighted_load = 0.5 * (sys_load + self_load);

   return weighted_load < Num_Processors;
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

private define ensure_logfile_exists (logfile)
{
   if (NULL != stat_file (logfile))
     return;

   variable fp = fopen (logfile, "w");
   if (fp == NULL)
     throw IOError, "opening logfile $logfile"$;

   if (0 != fclose (fp))
     {
        variable msg = sprintf ("closing logfile %s (%s)",
                                logfile, errno_string(errno));
        throw IOError, msg;
     }
}

private define run_executable (obj, file, run_dir, logfile)
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

   ensure_logfile_exists (logfile);

   variable p, dir_str = "";
   if (run_dir != NULL)
     {
        p = new_process (argv; dir=run_dir, stdout=">>${logfile}"$, dup2=1);
        dir_str = sprintf (" in %s", run_dir);
     }
   else
     {
        p = new_process (argv; stdout=">>${logfile}"$, dup2=1);
     }

   variable s = p.wait(WNOHANG);
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
   daemon_log (LOG_INFO, msg);

   if (s.pid == 0)
     {
        Num_Running += 1;
        return 0;
     }

   if (s.exited)
     {
        daemon_log (LOG_INFO,
                    sprintf ("child exited with status %d", s.exit_status));
     }
   else if (s.signal)
     {
        daemon_log (LOG_INFO,
                    sprintf ("child terminated by signal %d %s",
                             s.signal, (s.coredump ? "(coredump)" : "")));
     }

   return -1;
}

private define claim_with_rename (obj, path)
{
   if (path == NULL)
     return -1;

   variable cl = obj.client_data;

   % If we're busy, let somebody else handle it.
   if (0 == have_idle_processors (cl.exec_num_parallel))
     return 0;

   variable bname = path_basename (path);
   variable dname = path_dirname (path);
   variable claimed_bname = ".${bname}"$;
   variable claimed_path = path_concat (dname, claimed_bname);
   if (rename (path, claimed_path) != 0)
     return 0;

   if (run_executable (obj, claimed_path, cl.exec_root_dir, cl.exec_logfile) < 0)
     return -1;

   obj.num_processed++;

   return 1;
}

private define claim_with_subdir_rename (obj, path)
{
   if (path == NULL)
     return -1;

   variable cl = obj.client_data;

   % If we're busy, let somebody else handle it.
   if (0 == have_idle_processors (cl.exec_num_parallel))
     return 0;

   variable
     bname = path_basename (path),
     dname = path_dirname (path);

   % First, claim the file by renaming.
   % If we fail, somebody else got their first, but that's ok.

   variable claimed_bname = sprintf (".%s@%s", obj.host_pid, bname);
   variable claimed_path = path_concat (dname, claimed_bname);
   if (rename (path, claimed_path) != 0)
     return 0;

   % Now that we've claimed the file, we're free to operate
   % without competition as long as we create nothing that
   % matches the current search pattern. It's the user's
   % responsibility to choose search patterns and naming
   % patterns that don't conflict.

   % Create a subdirectory and rename into it, reverting to
   % the original basename.
   variable e;
   try (e)
     {
        variable exec_subdir = _$(cl.exec_subdir);
        if (cl.exec_subdir_is_regex)
          {
             variable m = pcre_matches (exec_subdir, bname);
             if (m == NULL)
               throw RunTimeError, "regex mismatch: $bname"$;
             exec_subdir = m[-1];
          }
        variable rename_dir = path_concat (obj.incoming_dir, exec_subdir);
        if (0 != mkdir_p (rename_dir))
          throw IOError, "creating rename directory: $rename_dir"$;
     }
   catch IOError, RunTimeError:
     {
        daemon_log (LOG_ERR, sprintf ("%s: %s", e.descr, e.message));
        () = rename (claimed_path, path);
        return -1;
     }

   try (e)
     {
        variable rename_path = path_concat (rename_dir, bname);
        if (rename (claimed_path, rename_path) != 0)
          throw IOError, "renaming to $rename_path"$;
     }
   catch IOError:
     {
        daemon_log (LOG_ERR, sprintf ("%s: %s", e.descr, e.message));
        () = rename (claimed_path, path);
        () = rmdir (rename_dir);
        return -1;
     }

   % Now that the claimed file is in the desired subdirectory
   % under its original basename, we can run the executable
   % if one has been specified.

   variable run_dir;
   if (cl.exec_root_dir == NULL)
     run_dir = rename_dir;
   else
     {
        try (e)
          {
             run_dir = path_concat (cl.exec_root_dir, exec_subdir);
             if (0 != mkdir_p (run_dir))
               throw IOError, "creating run directory $run_dir"$;
          }
        catch IOError:
          {
             daemon_log (LOG_ERR, sprintf ("%s: %s", e.descr, e.message));
             () = rename (rename_path, path);
             () = rmdir (rename_dir);
             return -1;
          }
     }

   % From here on, the executable is responsible for
   % rename_dir and all its contents (and run_dir as well).
   if (run_executable (obj, rename_path, run_dir, cl.exec_logfile) < 0)
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
    -d|--daemon      Run as a daemon
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
   newest = 0,
   incoming_dir = NULL,
   file_glob = "*",
   wait_sec = 1.0,
   logfile_name = "/tmp/cachemon.log",
   exec_subdir = Host_Name,
   exec_subdir_is_regex = 0,
   exec_root_dir = NULL,
   exec_name = NULL,
   exec_num_parallel = 1,
   exec_logfile = "/dev/null"
};

private define load_config_file (file)
{
   if (path_is_absolute(file))
     () = evalfile (file);
   else
     () = evalfile (path_concat (".", file));
   return _P;
}

private variable Temp_Pid_File = NULL;
private define delete_pidfile ()
{
   if (Temp_Pid_File != NULL)
     () = remove (Temp_Pid_File);
}

private define make_pidfile (pgm_name)
{
   variable pid_dir = "/var/tmp/${USER}/${pgm_name}"$;
   if (0 != mkdir_p (pid_dir))
     {
        throw ApplicationError, "*** Error creating $pid_dir"$;
     }
   variable pid_file = path_concat (pid_dir, string(getpid()));
   variable fp = fopen (pid_file, "w");
   if ((fp == NULL) || (0 != fclose (fp)))
     {
        variable msg = sprintf ("*** Error creating %s (%s)",
                                pid_file, errno_string(errno));
        throw IOError, msg;
     }
   Temp_Pid_File = pid_file;
   atexit (&delete_pidfile);
}

define slsh_main()
{
   variable run_as_daemon = 0;
   variable claim_via_rename_only = 0;
   variable exec_args = NULL;

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("r|rename", &claim_via_rename_only; inc);
   opts.add ("d|daemon", &run_as_daemon; inc);
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

   set_executable (p, exec_args);

   variable dirmon_task = claim_via_rename_only ?
     &claim_with_rename : &claim_with_subdir_rename;

   variable dir = new_dirmon (p.incoming_dir; glob = p.file_glob,
                              task = dirmon_task,
                              delay = p.wait_sec, client_data = p);

   if (run_as_daemon)
     {
        variable log_dir = path_dirname (p.logfile_name);
        if (0 != mkdir_p (log_dir))
          {
             throw ApplicationError,
               "*** Error: creating daemon log directory: $log_dir"$;
          }
        daemonize ("cachemon", p.logfile_name);
        make_pidfile ("cachemon");
     }

   Sighup_Received = 0;
   catch_sighup();
   catch_sigchild();

   while (Sighup_Received == 0)
     {
        if (-1 == dir.monitor (p.newest))
          break;
        dir.wait();
     }

   if (Sighup_Received)
     daemon_log (LOG_INFO, "received SIGHUP");

   wait_for_processes_to_exit ();
}

