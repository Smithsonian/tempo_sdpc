module m_initialize

private
public initialize

contains

  subroutine initialize(rc)

    use m_vars
    use m_LUN_set
    use m_pgs_include
    implicit none
    !-------------------------------------------------------------------------
    !         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
    !-------------------------------------------------------------------------
    !BOP
    !
    ! !ROUTINE:  initialize
    ! 
    ! !DESCRIPTION: initialize reads level 1b data including
    !		solar irradiance and observed radiance
    !
    ! !CALLING SEQUENCE: 
    !
    !        call initialize
    !     
    ! !INPUT PARAMETERS: none (in modules)   
    !
    ! !OUTPUT PARAMETERS:  
    integer, intent(out)         :: rc        ! Error return code:
    !  0   all is well
    !  1   files not found
    !
    ! !SEE ALSO:  
    !
    ! !REVISION HISTORY: 
    !
    !  05Jan01   Joiner     original fortran 90
    !  14Mar02   Vasilkov   modified to read OMCTPo.pcf, modifications 
    !			marked with **********
    !
    !EOP
    !-------------------------------------------------------------------------
    !
    integer :: i               !,j, iret
    integer :: iarg=0
    integer :: argc, iargc
    character*255 ::  argv
    character*255 ::  myname, pcfpath 
    !*********************************************************************
    !include 'PGS_PC.f'
    !include 'PGS_PC_9.f'
    !include 'PGS_SMF.f'
    !include 'PGS_IO.f'
    !include 'PGS_IO_1.f'
    !*********************************************************************
    integer(kind=4), EXTERNAL :: pgs_pc_getnumberoffiles, pgs_pc_getreference
    integer(kind=4), EXTERNAL :: pgs_pc_getuniversalref, pgs_pc_getconfigdata
    integer(kind=4) :: returnstatus, pcf_int
    CHARACTER(LEN=200) :: buf
    !*********************************************************************

    myname = trim('initialize: ')

    !===============================
    ! read command line information
    ! and do argument check
    !===============================
    argc = iargc()
    if ( argc < 0 ) then
      print *, trim(myname)//' not enough inputs'
      call ret_usage()
    else
      do i = 1, 32767
        iarg = iarg + 1
        if ( iarg > argc ) go to 111
        call GetArg ( iArg, argv )
        if(index(argv,'-h ') > 0) then
          call ret_usage()
        else if(index(argv,'-p ') > 0) then
          if ( iarg+1 > argc ) call ret_usage()
          iarg = iarg + 1
          call GetArg ( iArg, argv )
          read(argv,*,err=500) iprt
        else if(index(argv,'-i ') > 0) then
          if ( iarg+1 > argc ) call ret_usage()
          iarg = iarg + 1
          call GetArg ( iArg, argv )
          read(argv,*,err=500) interp
        else if(index(argv,'-rc ') > 0) then
          if ( iarg+1 > argc ) call ret_usage()
          iarg = iarg + 1
          call GetArg ( iArg, argv )
          read(argv,*,err=500) resource_file
          read_resource_file=.true.
        else if(index(argv,'-noret ') > 0) then
          noret = .true.
        else
          if (nfiles < maxfiles) then
            nfiles = nfiles + 1
            filenames(nfiles) = argv
          else
            print *,'l1bnamelist: too many files ',i, maxfiles
          endif
        endif
      enddo
    endif
111 continue

    !if no filename specified on command line, the set nfiles=1
    !==========================================================
    if (nfiles == 0) then
      nfiles=1
      no_cl_filename=.true. ! no command line filename
    endif

    !***********************************************************************
    !read OMCLDRR.pcf
    !check if PCF exists assuming PGE resides in $(PGE_NAME)/bin directory
    !---------------------------------------------------------------------
    rc=0
    status=1
    returnstatus=1
    if (iprt > 0) print *,'initialize: checking for pcf file'
    call getenv('PGS_PC_INFO_FILE',pcfpath)
    inquire(file=pcfpath,exist=ex)
    if (iprt > 0) print *,'initialize: file status ',ex
    if (ex) then
      version = 1
      status = pgs_pc_getreference ( L1B_LUN, version, filename)
      if (iprt > 0) print *,trim(myname)//' get_pc_reference status l1b filename',&
           status, trim(filename)

      version = 1

      returnstatus = pgs_pc_getconfigdata(using_resid_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        using_resid = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting using_resid = ',using_resid
      endif

      returnstatus = pgs_pc_getconfigdata(write_resid_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        write_resid = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting write_resid = ',write_resid
      endif

      returnstatus = pgs_pc_getconfigdata(do_o3_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        do_o3 = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting write_resid = ',write_resid
      endif

      returnstatus = pgs_pc_getconfigdata(write_obs_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        write_obs = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting write_obs = ',write_obs
      endif

      returnstatus = pgs_pc_getconfigdata(using_cal_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        using_cal = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting using_cal   = ',using_cal
      endif

      returnstatus = pgs_pc_getconfigdata(no_ret_ps_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        no_ret_ps = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting no_ret_ps = ',no_ret_ps
      endif

      returnstatus = pgs_pc_getconfigdata(no_ret_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        noret = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting noret = ',noret
      endif

      returnstatus = pgs_pc_getconfigdata(transient_chk,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        transient_check = pcf_int == 1
        if (iprt >= 1) print *,'initialize: setting transient_check = ',transient_check
      endif

      returnstatus = pgs_pc_getconfigdata(wmin_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) wmin
        if (iprt >= 1) print *,'initialize: setting wmin = ',wmin
        set_wmin=.true.
      endif
      returnstatus = pgs_pc_getconfigdata(wmax_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) wmax
        if (iprt >= 1) print *,'initialize: setting wmax = ',wmax
        set_wmax=.true.
      endif
    endif

    if (cloud_clear) ny=2

    !***********************************************************************

    return
500 call ret_usage()

  end subroutine initialize

  !*****************************************************************************
  subroutine ret_usage()

    implicit none

    print *
    print *, &
         'Usage:  cloud_ret.x [-p iprt] [-i interp] [-noret]'
    print *, '[input files] '
    print *
    print *,  'where'
    print *
    print *, '-p  iprt    printout level flag (1:least amount of '
    print *, '             (printout, >1 more printouts, default=1)'
    print *, '                                             '
    print *, '-i interp   interpolation level (2 lowest, 4 medium, 5 highest)'
    print *
    print *, '-noret      do not actually perform retrieval (for testing)'    
    print *
    print *, 'input files level 1b radiance files (process several at a time)'
    print *
    print *, 'Last Revised: 28 March 2002      (J. Joiner)  '
    print *
    call exit(7)

  end subroutine ret_usage

end module m_initialize
