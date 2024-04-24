#! /usr/bin/env slsh

private variable Num_Plan_Days = 21;
private variable Window_Start_Days = [0.0, 4.65139, 7.0, 9.2139, 14.0];
private variable Window_Duration = 3000L;
private variable Window_Margin = 300L;

private variable Epoch = 315964800L;  % time_t for 1980-01-06T00:00:00Z (GPS epoch)
private variable Sec_Per_Day = 86400L;

define print_maneuver_windows (fp, t0)
{
   variable window_start_time = t0 + int(Sec_Per_Day * Window_Start_Days);
   variable window_end_time = window_start_time + Window_Duration;
   variable window_duration = window_end_time - window_start_time;
   variable start_time = window_start_time + Window_Margin;
   variable end_time = window_end_time - Window_Margin;
   variable centroid_time = (start_time + end_time)/2;
   variable dummy_remainder = "MAINT,20,1.2,3.4,5.6,7,Bogus example entry";

   variable n = length(window_start_time);
   foreach ([0:n-1])
     {
	variable i = ();
	() = fprintf (fp, "%ld,%ld,%ld,%ld,%ld,%ld,%s\n",
		      window_start_time[i],
		      window_end_time[i],
		      window_duration[i],
		      start_time[i],
		      end_time[i],
		      centroid_time[i],
		      dummy_remainder);
     }
}

define print_maneuver_table (fp, t_beg, t_end, plan_id)
{
   variable col_heads = "window_start_time,window_end_time,window_duration,start_time,end_time,centroid_time,maneuver_type,thruster_set,dv_t,dv_n,dv_r,num_pulses,comments";
   variable dummy_remainder = "0,0.0,0,0,0,NULL,0,0.0,0.0,0.0,0,NULL";

   () = fprintf (fp, "%s\n", col_heads);
   () = fprintf (fp, ":table_begin_time=%ld,%s\n", t_beg, dummy_remainder);
   () = fprintf (fp, ":table_end_time=%ld,%s\n", t_end, dummy_remainder);
   () = fprintf (fp, ":plan_id=%s,%s\n", plan_id, dummy_remainder);
   () = fprintf (fp, ":tempo_epoch=%ld,%s\n", Epoch, dummy_remainder);

   print_maneuver_windows (fp, t_beg);
}

define slsh_main()
{
   if (__argc == 1)
     {
        vmessage ("Usage:  %s YYYY-MM-DD", path_basename(__argv[0]));
        exit (0);
     }

   variable day_string = __argv[1];
   variable year, mon, day;
   if (3 != sscanf (day_string, "%d-%d-%d", &year, &mon, &day))
     throw ApplicationError, "Parsing date specifier: $day_string";

   variable tm = gmtime(0);
   tm.tm_year = year - 1900;
   tm.tm_mon = mon - 1;
   tm.tm_mday = day;

   variable t_beg = timegm(tm);
   variable t_end = t_beg + Num_Plan_Days*Sec_Per_Day + Window_Duration;
   variable plan_id = "test-$day_string"$;
   variable fp = stdout;

   print_maneuver_table (fp, t_beg, t_end, plan_id);
}
