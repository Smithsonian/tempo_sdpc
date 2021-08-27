#! /usr/bin/env slsh

require ("cmdopt");

$1 = path_dirname (__FILE__);
prepend_to_slang_load_path (path_concat ($1, "../share/slsh/local-packages"));
prepend_to_slang_load_path ($1);
require ("pipeutil");

private define usage ()
{
   variable msg =
`Usage: o3p_partition.sl [opts] <host_index-num_hosts>
Options:
    -m|--mksub DIR     Create processing subdirectories in DIR
`;
   () = fprintf (stderr, msg);
   exit (0);
}

private define cmdopt_error (msg)
{
   () = fprintf (stderr, "%s\n", msg);
   usage ();
}

% Partition num objects into blocks which may differ in size
% by no more than one.
private define define_nearly_equal_blocks (this_block, num_blocks, num)
{
   variable block_size = num / num_blocks;
   variable residual = num - num_blocks * block_size;
   variable beg_step = this_block * block_size;

   if (residual > 0)
     {
        if (this_block < residual)
          {
             block_size += 1;
             beg_step += this_block;
          }
        else beg_step += residual;
     }

   % end_step is one past the end
   variable end_step = beg_step + block_size;
   if (end_step > num)
     {
        end_step = num;
     }

   return struct
     {
        beg_step=beg_step,
        end_step=end_step
     };
}

private variable Num_Cores_Per_Host = 20;

private define __init_num_cores_using_ntasks ()
{
   variable env = getenv ("SLURM_NTASKS");
   if (env == NULL)
     return;

   try
     {
        variable n = eval(env);
        if (_typeof(n) == Integer_Type)
          {
             Num_Cores_Per_Host = n;
          }
     }
   catch AnyError:
     {
        return;
     }
}

__init_num_cores_using_ntasks();

define host_partition (this_host)
{
   variable lst = {};

   variable num_cores_per_host = qualifier ("num_cores_per_host", Num_Cores_Per_Host);
   variable num_hosts = qualifier ("num_hosts", 3);
   variable bin_factor = qualifier ("bin_factor", 4);

   variable max_xtrack = 2047;

   % Subdivide on the binned coordinate
   variable
     nmin = 0,
     nmax = max_xtrack/bin_factor;

   variable num = nmax - nmin + 1;

   % divide pixels among hosts
   variable h = define_nearly_equal_blocks (this_host, num_hosts, num);
   h.beg_step = (h.beg_step + nmin) * bin_factor;
   h.end_step = (h.end_step + nmin) * bin_factor;

   % on each host, divide bins among cores
   variable hnum = h.end_step - h.beg_step + 1;
   variable num_bins = hnum / bin_factor;
   variable i, b, beg, end;

   _for i (0, num_cores_per_host-1, 1)
     {
        b = define_nearly_equal_blocks (i, num_cores_per_host, num_bins);
        beg = h.beg_step + b.beg_step * bin_factor;
        end = h.beg_step + b.end_step * bin_factor;

        % By construction, end-beg is a multiple of bin_factor
        variable p = struct
          {
             beg = beg,
             end = end,
             blkid = this_host * num_cores_per_host + i,
             core = i,
             num_bins = end-beg
          };
        list_append (lst, p);
     }

   return lst;
}

private define print_entry (p)
{
   variable status;
   status = fprintf (stdout,  "%d,%d,%d,%d,%d\n",
                     p.beg, p.end-1,
                     p.blkid,
                     p.core,
                     p.num_bins);
   if (status < 0)
     throw IOError;
}

private define create_dir (p, target_dir)
{
   variable dirname = path_concat (target_dir, sprintf ("block_%03d", p.blkid));
   variable pathname = path_concat (dirname, "block.txt");

   if (mkdir_p (dirname) != 0)
     throw IOError, "creating directory $dirname"$;

   variable fp = fopen (pathname, "w");
   if (fp == NULL)
     throw IOError, "opening $pathname for writing"$;

   if (fprintf (fp, "%d %d\n", p.beg, p.end-1) < 0)
     throw IOError, "writing to $pathname"$;

   if (fclose (fp) < 0)
     throw IOError, "closing $pathname"$;
}

private define write_array_bounds_file (ary, target_dir, filename)
{
   variable dirname, pathname;

   if (target_dir != NULL)
     {
        dirname = target_dir;
        pathname = path_concat (dirname, path_basename (filename));
     }
   else
     {
        dirname = path_dirname (filename);
        pathname = filename;
     }

   if (mkdir_p (dirname) != 0)
     throw IOError, "creating directory $dirname"$;

   variable fp = fopen (pathname, "w");
   if (fp == NULL)
     throw IOError, "opening $pathname for writing";

   if (fprintf (fp, "%d-%d\n", ary[0].blkid, ary[-1].blkid) < 0)
     throw IOError, "writing to $pathname";

   if (fclose (fp) < 0)
     throw IOError, "closing $pathname";
}

define slsh_main ()
{
   variable this_host, num_hosts;
   variable target_dir = NULL;
   variable array_bounds_file = NULL;

   variable opts = cmdopt_new (&cmdopt_error);
   opts.add ("m|mksub", &target_dir; type="string");
   opts.add ("b|bounds", &array_bounds_file; type="string");
   variable i = opts.process (__argv,1);

   if (__argc == 1 || i < 1)
     usage();

   if (2 != sscanf (__argv[i], "%d-%d", &this_host, &num_hosts))
     usage();

   if ((this_host >= num_hosts)
       || (this_host < 0))
     usage();

   variable lst = host_partition (this_host; num_hosts=num_hosts);
   variable ary = list_to_array (lst);

   if (target_dir == NULL)
     {
        array_map (Void_Type, &print_entry, ary);
     }
   else
     {
        array_map (Void_Type, &create_dir, ary, target_dir);
     }

   if (array_bounds_file != NULL)
     {
        write_array_bounds_file (ary, target_dir, array_bounds_file);
     }
}
