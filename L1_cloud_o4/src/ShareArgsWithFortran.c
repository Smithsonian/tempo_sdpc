/*------------------------------------------------------------------------------*\
** ShareArgsWithFortran
**
** This routine passes command-line arguments from C to Fortran.  It causes
** the getarg() routine to work properly.  A unique method is needed for each
** different compiler, so we use conditional compilation based on the Fortran
** compiler currently used.
\*------------------------------------------------------------------------------*/

#ifdef pgf
      int    __argc_save;
      char** __argv_save;
      static int init_for_pghpf_init = 0;
#endif

/*------------------------------------------------------------------------------*/

int ShareArgsWithFortran(int argc, char** argv)
{

#ifdef pgf
     __argc_save = argc;
     __argv_save = argv;
       pghpf_init(&init_for_pghpf_init);
#endif

#ifdef ifort
       for_rtl_init_(&argc, argv);
#endif

#ifdef g77
       f_setarg(argc, argv);
#endif

#ifdef gfortran
      _gfortran_set_args(argc, argv);
#endif

    return(0);
}

/*------------------------------------------------------------------------------*/
