MODULE GEMSTOOL_binset

  INTEGER , PARAMETER :: mset = 10
  INTEGER :: nbin, nwav
  TYPE :: binset
     INTEGER :: nbin,nwav
     INTEGER, DIMENSION(mset):: neof, cc
     REAL, DIMENSION (mset)  :: maxtau, dtau
  END TYPE binset
  PRIVATE 
  PUBLIC :: binset, setwin, setbin, get_bin,  bin_UV
  LOGICAL :: do_debug =.false.
CONTAINS

  SUBROUTINE setwin (do_uv, nw, waves, npca, winpix,winid)
  IMPLICIT NONE
  ! input variables
  LOGICAL, INTENT(IN) :: do_uv
  INTEGER, INTENT(IN) :: nw
  REAL (KIND=8), INTENT(IN), DIMENSION (nw) :: waves
  ! output variables
  INTEGER, INTENT(OUT) :: npca
  INTEGER, ALLOCATABLE :: winid(:), winpix(:,:)
  !local variables
  INTEGER :: mpca, w1, w2,ipca, ntmp, tmpid(13), tmpwin(13,2)
  REAL (KIND=8),DIMENSION(1:13):: pcalim
  IF (do_uv) THEN
    mpca = 7
    pcalim(1:mpca) = (/340.0,350.0,360.0,390.,420.,460.,490.,500./)
  ELSE
    mpca = 13
    pcalim(1:mpca) = (/520., 539.,550.,565.,595.,625.,650., &
                    686.75, 690.,705.,716.0,735.,752.0/)
  ENDIF
  npca = 0 ; w1 = 1
  do ipca = 1, mpca
    w2= count( waves(1:nw)<=pcalim(ipca))
    ntmp = w2 - w1 + 1
    if (ntmp == 0 ) cycle
    npca = npca + 1
    tmpid(npca ) = ipca
    tmpwin(npca,1) = w1 ; tmpwin(npca,2) = w2
    if (ntmp < 5) then
      write(*,*) 'warning in setting binning windows in pcartm !!!', ipca, ntmp
      npca = npca -1
      tmpwin(npca, 2) = w2
    endif
    w1 = w2 + 1
  enddo
  allocate(winid(npca), winpix(npca,2))
  winid = tmpid(1:npca) ; winpix(1:npca,:) = tmpwin(1:npca,:)
  if (w2 /= nw ) then
   write (*,*) 'errors in setting binning window '; STOP 1
  endif
  if (do_debug) print * , 'Finishing binset'
END SUBROUTINE setwin

  SUBROUTINE setbin (do_uv,id,tau1,tau2,sd,amf,bin)
  LOGICAL, INTENT(IN) :: do_uv
  INTEGER, INTENT(IN) :: id
  REAL (KIND=4), INTENT(IN) :: tau1,tau2, sd,amf
  TYPE(binset), INTENT(OUT) :: bin
  ! local variables
  INTEGER :: NEOF
  REAL (KIND=4) :: dtau
  TYPE(binset)  :: bin1, bin2

  IF (do_uv) THEN
    dtau=(tau2-tau1)
    IF (ID ==  1) THEN  ! O3 is dominant
      CALL bin_o3band (tau1,tau2,amf,bin)
      !CALL bin_uv(tau1, tau2, 0.5, 2, bin)
    ELSE IF (ID == 2) THEN ! 340-350
      ! 350 - 365, NO2>O4> O3, SO2,HCHO
      ! 365 - 390, NO2>O4> O3>SO2, NO HCHO
     CALL bin_UV (tau1,tau2+0.001,1.0, 2,bin) ! o3 dominant
     print * , 'bin_uv: id =2'
    ELSE IF (ID == 3) THEN !350-360
     CALL bin_UV (tau1,tau2,dtau, 4,bin) ! o4
    ! CALL bin_UV (tau1,tau2,dtau/2.0, 4,bin) ! o4, no2
    ELSE IF (ID == 4) THEN ! 360-390
     CALL bin_UV (tau1,tau2,dtau/2.0, 4,bin) ! very compelxt
    ELSE IF (ID == 5) THEN !> 390
     CALL bin_UV (tau1,tau2,dtau, 4,bin) ! No2 dominant
    ELSE IF (ID ==6) THEN !>420
     CALL bin_UV (tau1,tau2,dtau/2.0, 4,bin) ! No2 dominant
    ELSE IF (ID ==7) THEN !>460
     CALL bin_UV (tau1,tau2,dtau/3.0, 4,bin) ! O3 dominent
    ELSE IF (ID==8) THEN !>490
     CALL bin_UV (tau1,tau2,0.5, 2,bin) ! O3 dominent
    ENDIF

 ELSE ! vis
  ! cc, mininum number of samplings. if sampling is insufficient, they are
  ! included into previous group.
      !neof = 4
      !IF (sd < 1.0E-3 ) THEN
      !  CALL bin_uv(tau1, tau2, tau2-tau1, 3, bin)
      !ELSE
      IF (sd < 1.0) THEN
      nbin = 3
      bin%maxtau(1:nbin) = [0.01,0.1,100.]
      bin%dtau(1:nbin)   = 100
      bin%neof(1:nbin)   = [4, 4, 3]
      bin%cc = 10
      bin%nbin   = nbin
      ELSE
      nbin =3
      bin%maxtau(1:nbin) = [0.1,1.0,100.]
      bin%dtau(1:nbin)   = 100
      bin%neof(1:nbin)   = [4,3,3]
      bin%cc = 10
      bin%nbin   = nbin
      ENDIF
      bin%neof(1:nbin)   = [3, 3, 3]
      !ENDIF
   ENDIF
 END SUBROUTINE setbin

  SUBROUTINE get_bin (wav1,wav2,tau1,tau2,mtau,sd,amf,bin)
  REAL (KIND=4), INTENT(IN) :: wav1, wav2
  REAL (KIND=4), INTENT(IN) :: tau1,tau2, mtau,sd,amf
  TYPE(binset), INTENT(OUT) :: bin
  ! local variables
  INTEGER :: NEOF
  REAL (KIND=4) :: dtau
  TYPE(binset)  :: bin1, bin2
    IF (wav2 <= 340) THEN  ! O3 is dominant
      CALL bin_o3band (tau1,tau2,amf,bin) 
      !CALL bin_uv(tau1, tau2, 0.5, 2, bin)
    ELSE IF (wav1 >= 340 .and. wav2 <= 500) THEN 
      ! 350 - 365, NO2>O4> O3, SO2,HCHO
      ! 365 - 390, NO2>O4> O3>SO2, NO HCHO
    
      CALL bin_UV (tau1,tau2,0.5, 2,bin) ! o3 dominant
      return
      dtau=(tau2-tau1)
      IF (wav1 >= 340) CALL bin_UV (tau1,tau2,0.5, 2,bin) ! o3 dominant
      IF (wav1 >= 350) CALL bin_UV (tau1,tau2,dtau/2.0, 4,bin) ! o2 vs no2
      IF (wav1 >= 360) CALL bin_UV (tau1,tau2,dtau/2.0, 4,bin) ! very compelxt
      IF (wav1 >= 390) CALL bin_UV (tau1,tau2,dtau, 4,bin) ! No2 dominant
      IF (wav1 >= 420) CALL bin_UV (tau1,tau2,dtau/2.0, 4,bin) ! No2 dominant
      IF (wav1 >= 460) CALL bin_UV (tau1,tau2,dtau/3.0, 4,bin) ! O3 dominent
      IF (wav1 >= 490) CALL bin_UV (tau1,tau2,0.5, 2,bin) ! O3 dominent
    ELSE IF (wav1 >=500 .and. wav2 <=752) THEN 
      neof = 4
      IF (sd < 1.0E-3 ) THEN 
        CALL bin_uv(tau1, tau2, tau2-tau1, 3, bin)
      ELSE  
      IF (sd < 1.0) THEN 
      nbin = 3 ; neof =4
      bin%maxtau(1:nbin) = [0.01,0.1,100.]
      bin%dtau(1:nbin)   = 100
      bin%neof(1:nbin)   = [4, 4, 3]
      bin%cc = 10
      bin%nbin   = nbin
      ELSE
      nbin =3 ; neof =4
      bin%maxtau(1:nbin) = [0.1,1.0,100.]
      bin%dtau(1:nbin)   = 100
      bin%neof(1:nbin)   = [4,3,3]
      bin%cc = 10
      bin%nbin   = nbin
      ENDIF
      ENDIF
   ELSE 
     WRITE(*,*) 'GEMSTOOL_binset:Error in get_bin', wav1, wav2
    ENDIF
  END SUBROUTINE
 
  ! ***************** UV *****************************
  ! O3
  ! ***************************************************
  SUBROUTINE bin_o3band (tau1,tau2,amf,bin)
    REAL, INTENT(IN) :: tau1, tau2, amf
    TYPE(binset), INTENT(OUT) :: bin
    ! 0.05%
    ! @ high TOZ, accuracies could be better with more bins(0.5-1.0) or EOFs (1.0-2.5) @ high SZA
    ! @ high TOZ
    nbin = 7
    bin%maxtau(1:nbin) = [-1.6,-1.2,0.0,0.5,3.5,4.5,20.0]
    bin%dtau(1:nbin)   = [2.0,0.4,0.4,0.5,0.6,1.0,2.0]
    bin%neof(1:nbin)   = [1,1,1,1,2,2,2]
    if (amf > 70.0) THEN 
     nbin=8
     bin%maxtau(1:nbin) = [-1.5,-0.7,0.4,0.7, 2.5,3.5,4.5,10.0]
     bin%dtau(1:nbin)   = [2.0,  1.2,0.35,0.3,0.6,1.0,1.0,2.0]
     bin%neof(1:nbin)   = [1,1,     1,  1,3,3,2,2]
    endif

    if (tau1 > bin%maxtau(1)-bin%dtau(1)) then 
        bin%dtau(1) = bin%maxtau(1)-tau1
    else
        bin%dtau(1) = (bin%maxtau(1)-tau1)/2.0
    endif
    bin%cc(1:nbin)     = 1
    bin%nbin   = nbin
  END SUBROUTINE

  SUBROUTINE bin_UV (tau1,tau2,dtau,neof,bin)
    REAL, INTENT(IN) :: tau1, tau2, dtau
    INTEGER, INTENT(IN) :: neof
    TYPE(binset), INTENT(OUT) :: bin
    nbin = 1
    bin%maxtau(1:nbin) = [tau2]
    bin%dtau(1:nbin)   = [dtau]
    bin%neof(1:nbin)   = [neof]
    bin%cc(1:nbin)     = 1   !5
    bin%nbin   = nbin
    !print * , dtau, neof, tau1, tau2
  END SUBROUTINE

  SUBROUTINE bin_vis (tau1,tau2,mtau,neof,bin)
    REAL, INTENT(IN) :: tau1, tau2, mtau
    INTEGER, INTENT(IN) :: neof
    TYPE(binset), INTENT(OUT) :: bin
    !if (mtau > 1.0) THEN 
    nbin = 4
    bin%maxtau(1:nbin) = [0.01,0.1,1.,100.]
    bin%dtau(1:nbin)   = [0.01,0.045,0.45,100.]
    bin%neof(1:nbin)   = neof
    !ELSE
    nbin = 3
    bin%maxtau(1:nbin) = [0.01,0.1,100.]
    bin%dtau(1:nbin)   = 100
    bin%neof(1:nbin)   = neof
    !ENDIF
    bin%nbin   = nbin
  END SUBROUTINE
END MODULE GEMSTOOL_binset
