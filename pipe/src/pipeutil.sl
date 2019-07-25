
define mkdir_p (path)
{
   variable mode = qualifier ("mode", 0777);

   variable st = stat_file (path);
   if (st != NULL)
     {
        if (0 != stat_is ("dir", st.st_mode))
          return 0;
        else return -1;
     }

   ifnot (is_substr (path, "/"))
     {
        if (mkdir (path, mode) != 0)
          return (errno == EEXIST) ? 0 : -1;
     }

   variable dirs = strtok (path, "/");
   if (path[0] == '/') dirs[0] = "/" + dirs[0];

   variable i, n = length(dirs);
   variable s = "";
   _for i (0, n-1, 1)
     {
        s = path_concat (s, dirs[i]);
        if (mkdir(s, mode) != 0)
          {
             if (errno != EEXIST)
               return -1;
          }
     }
   return 0;
}
