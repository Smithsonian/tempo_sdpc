!>Routines to read operating settings, and print usage statement if incorrect
module m_initialize

  private ret_usage
  public initialize 

contains

  !>Read operating settings from PCF file, and optionally from command line
  !-------------------------------------------------------------------------
  !
  ! !ROUTINE:  initialize
  ! 
  ! !CALLING SEQUENCE: 
  !
  !        call initialize
  !     
  ! !INPUT PARAMETERS: none (in modules)   
  !> @param errstat error return code, non-zero indicates failure
  !
  ! !OUTPUT PARAMETERS:  
  !
  ! !SEE ALSO:  
  !
  ! !REVISION HISTORY: 
  !
  !> @author  05Jan01  Joiner      original fortran 90
  !> @author  14Mar02  Vasilkov    modified to read OMCTPo.pcf
  !>  23Mar15  O'Sullivan  updating for TEMPO
  !
  !-------------------------------------------------------------------------
  subroutine initialize(errstat)

    use m_vars
    use m_LUN_set
    use m_pgs_include
    use tell_module
    implicit none
    !-------------------------------------------------------------------------
    ! !INPUT PARAMETERS: none (in modules)   
    integer, intent(inout)         :: errstat        ! Error return code:
    !  0   all is well
    !  -1  problem
    !-------------------------------------------------------------------------
    ! Local variables
    integer :: i
    integer :: iprt=0   ! verbosity level
    integer :: iarg=0
    integer :: argc, iargc
    character*255 ::  argv, logmsg
    character*255 ::  pcfpath 
    integer(kind=4), EXTERNAL :: pgs_pc_getnumberoffiles, pgs_pc_getreference
    integer(kind=4), EXTERNAL :: pgs_pc_getuniversalref, pgs_pc_getconfigdata
    integer(kind=4) :: returnstatus, pcf_int
    CHARACTER(LEN=200) :: buf
    !*********************************************************************

    if (errstat /= 0) return


    !===============================
    ! read command line information
    ! and do argument check
    !===============================
    argc = iargc()
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
      else if(index(argv,'-nc_swath ') > 0) then
        if ( iarg+1 > argc ) call ret_usage()
        iarg = iarg + 1
        call GetArg ( iArg, argv )
        read(argv,*,err=500) nc_swathname
      else if(index(argv,'-nc_only ') > 0) then
        read_he4 = .false.
      else if(index(argv,'-noret ') > 0) then
        noret = .true.
      endif
    enddo
111 continue


    ! Set logging level from iprt
    call tell_set_log_level(iprt)

    !***********************************************************************
    !read OMCLDRR.pcf
    !check if PCF exists based on PGS_PC_INFO_FILE environment variable
    !---------------------------------------------------------------------
    status=1
    returnstatus=1
    call tell_log(1,'initialize: checking for pcf file')
    call get_environment_variable('PGS_PC_INFO_FILE',pcfpath)
    inquire(file=pcfpath,exist=ex)
    write(logmsg,"(A,L1)") 'initialize: file status ',ex
    call tell_log(1,logmsg)
    if (ex) then
      version = 1
      status = pgs_pc_getreference ( L1B_LUN, version, filename)
      write(logmsg,"(A12,A24,I1)") 'initialize: ', &
           'get_pc_reference status ',status
      call tell_log(1,logmsg)
      write(logmsg,"(A13,A)") 'l1b filename ',trim(filename)
      call tell_log(1,logmsg)

      version = 1

      returnstatus = pgs_pc_getconfigdata(using_resid_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        using_resid = pcf_int == 1
        write(logmsg,"(A34,L1)") 'initialize: setting using_resid = ', &
             using_resid
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(write_resid_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        write_resid = pcf_int == 1
        write(logmsg,"(A34,L1)") 'initialize: setting write_resid = ', &
             write_resid
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(do_o3_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        do_o3 = pcf_int == 1
        write(logmsg,"(A28,L1)") 'initialize: setting do_o3 = ', &
             do_o3
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(write_obs_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        write_obs = pcf_int == 1
        write(logmsg,"(A32,L1)") 'initialize: setting write_obs = ', &
             write_obs
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(no_ret_ps_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        no_ret_ps = pcf_int == 1
        write(logmsg,"(A32,L1)") 'initialize: setting no_ret_ps = ', &
             no_ret_ps
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(no_ret_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        noret = pcf_int == 1
        write(logmsg,"(A28,L1)") 'initialize: setting noret = ', &
             noret
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(transient_chk,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        transient_check = pcf_int == 1
        write(logmsg,"(A38,L1)") 'initialize: setting transient_check = ', &
             transient_check
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(test_solar_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        test_solar = pcf_int == 1
        write(logmsg,"(A33,L1)") 'initialize: setting test_solar = ', &
             test_solar
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(add_shift_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) add_shift
        write(logmsg,"(A32,L1)") 'initialize: setting add_shift = ', &
             add_shift
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(using_spline_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) pcf_int
        using_spline = pcf_int == 1
        write(logmsg,"(A35,L1)") 'initialize: setting using_spline = ', &
             using_spline
        call tell_log(1,logmsg)
      endif

      returnstatus = pgs_pc_getconfigdata(wmin_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) wmin
        write(logmsg,"(A27,L1)") 'initialize: setting wmin = ', &
             wmin
        call tell_log(1,logmsg)
        set_wmin=.true.
      endif
      returnstatus = pgs_pc_getconfigdata(wmax_LUN,buf)
      IF(returnstatus == 0 ) THEN
        read(buf,*) wmax
        write(logmsg,"(A27,L1)") 'initialize: setting wmax = ', &
             wmin
        call tell_log(1,logmsg)
        set_wmax=.true.
      endif

    else  !ex=.false., PCF does not exist or environment variable not set

      call tell_log(1,'initialize: PCF file not found')
      errstat=-1
      call tell_error (tell_io_error, &
           "read_cld_dimensions: failed", &
           errstat)
      return

    endif


    !***********************************************************************

    return
500 call ret_usage()

  end subroutine initialize

  !****************************************************************************
  !>Print usage statement to the command line
  subroutine ret_usage()

    implicit none

    print *
    print *, &
         'Usage:  cloud_ret.x [-p iprt] [-nc_swath swathname] [-noret] [-nc_only]'
    print *
    print *,  'where'
    print *
    print *, '-p  iprt    printout level flag (1:least amount of '
    print *, '             (printout, >1 more printouts, default=1)'
    print *, '                                             '
    print *, '-nc_swath <swathname> override the default netCDF'
    print *, '                      swathname'
    print *
    print *, '-nc_only    read and write netCDF files only, no he4/5'
    print *
    print *, '-noret      do not actually perform retrieval (for testing)'    
    print *
    print *, "Last Revised: 20 August 2014      (E. O'Sullivan)  "
    print *

    stop 1

  end subroutine ret_usage

end module m_initialize
