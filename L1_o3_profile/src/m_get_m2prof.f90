!! 1) Demonstration of reconstruction of covariance matrix from its 
!!    Eigenvectors (EOFs)
!! 2) Construct a covariance matrix in a user-defined layers from 
!!    the Eigenvectors (EOFs) of covariance matrix defined on  
!!    levels. 
!! 3) JBAK, muse.f90 is modifed to be implemented for xliu's ozone profile code 


  MODULE m_get_m2prof
  USE OMSAO_variables_module, ONLY: atmdbdir
  USE m_utilities,   ONLY : reverse
  TYPE m2_output 
    REAL(KIND=8), DIMENSION (:), ALLOCATABLE  :: o3p
    REAL(KIND=8), DIMENSION (:,:), ALLOCATABLE  :: Sa
  END TYPE m2_output
  TYPE (m2_output):: m2du
  PUBLIC :: m2du, get_m2prof
  PRIVATE
  CONTAINS
  SUBROUTINE get_m2prof (which_var, the_var,nlyrs0, plev0, neof) 

    USE HDF5
    USE ISO_C_BINDING, ONLY : C_NULL_CHAR
    USE C_NULL_PTR_m, ONLY  : poo, NULLIFY_POO
    USE M2_Clim_m
    USE Zinsert_m
    USE Pres2Zstar_m,  ONLY : zstar2hPa, hPa2zstar
    USE ControlFile_m, ONLY : KVP, nKVP, ControlTxt_Reader
    USE OMSAO_variables_module, ONLY: atmdbdir, & 
        longitude=>the_lon,latitude=>the_lat, & 
        year=>the_year, month=>the_month,day=>the_day 
    USE m_utilities,   ONLY : reverse
    IMPLICIT NONE
    !===================================================
    !input variables
    !===================================================
    ! * vertical grid variables
    INTEGER, INTENT(IN)           :: nlyrs0, neof
    REAL(KIND=8), INTENT(IN)      :: plev0(0:nlyrs0) ! Ascending/Desending OK
                                                     ! plev0(0) or plev0(nlyrs0)
                                                     ! is surface pressure
    CHARACTER (LEN=3), INTENT(IN) :: which_var ! should be TPP, TO3, LAZ
    REAL(KIND=8), INTENT(IN)      :: the_var    
    !===================================================
    !local variables
    !===================================================
    REAL(KIND=8), PARAMETER :: Z1atm = 0.d0
    ! temporal output variables
    REAL(KIND=8), ALLOCATABLE, DIMENSION (:)   :: plev, zstar, lyrOz_Col,dO3prf_dVar
    REAL(KIND=8), ALLOCATABLE, DIMENSION (:,:) :: Sap
    ! helps
    LOGICAL, PARAMETER :: do_debug=.false.
    LOGICAL :: LFAIL
    LOGICAL :: SET     
    LOGICAL :: do_sort ! decide if it is decending or ascending
    CHARACTER(LEN=256) :: ClassDir, msg
    INTEGER :: ii, nlyrs, nlyrs_max, ymd(1:3)
    INTEGER :: ilev(1:2), nlyr_pt, nlyr_1atm, z1, z2
    INTEGER, DIMENSION (nlyrs0) :: ord
    REAL(KIND=8) :: Zins(1:2), Pt, Zt
    !=========================================================
    ! module name
    !=========================================================
    CHARACTER(LEN=9), PARAMETER  :: modulename = "M2CLIM_dr"

    ClassDir = ADJUSTL(TRIM(atmdbdir))//'M2CLIMO3/DAY9to5_EOF_24X18/'
    M2AnnualClimFn = 'CLASS_Y2005-2016_YYY_LT09_17_24x18_EOF_B.h5'

    ClassDir = ADJUSTL(TRIM(atmdbdir))//'M2CLIMO3/DAY9to5_EOF_01X18/'
    M2AnnualClimFn = 'CLASS_Y2005-2016_YYY_LT09_17_01x18_EOF_B.h5'

    ymd(1:3) = (/year, month, day/)
    ! Layering treatment, put additional layers for top (78km) and bottoms (1atm)
    do_sort = .false.
    IF (plev0(1) > plev0(2)) do_sort = .true. ! it should be top-down

    nlyrs = nlyrs0
    nlyrs_max = nlyrs+2
    allocate (plev(nlyrs), zstar(nlyrs_max))

    IF (do_sort) THEN 
       ord (1:nlyrs) = (/(ii,ii = nlyrs, 1, -1)/)
       Pt = plev0(0)
       plev(1:nlyrs) = plev0(ord(1:nlyrs))
    ELSE
       pt = plev0(nlyrs)
       plev(1:nlyrs) = plev0(0:nlyrs-1)
    ENDIF
    zstar(1:nlyrs) = hPa2zstar(plev(1:nlyrs))
    Zt = hPa2zstar( Pt )
    Zins(1:2) = (/Zt, Z1atm/)
    ilev(1:2) = Zinsert( nlyrs_max, nlyrs, zstar(1:nlyrs_max), Zins(1:2) )
    nlyr_pt = ilev(1); nlyr_1atm = ilev(2)

    IF (do_debug) THEN
     WRITE(*,*) "total # of layers, nlyrs =", nlyrs, "nlyrs0", nlyrs0
     WRITE(*,'(a, f8.2, i5, f8.2)') "Pt, nlyr_pt, Zt    =", Pt, nlyr_pt, zstar(nlyr_pt)
     WRITE(*,*) "nlyr_1atm, Z1atm   =", nlyr_1atm, zstar(nlyr_1atm)
    ENDIF

    ALLOCATE (lyrOz_col(nlyrs), dO3prf_dVar(nlyrs), Sap(nlyrs, nlyrs))
    CALL M2CLIM_Stage( LFAIL, msg, M2ClimDir_k = ClassDir, & 
                       nEOFs_k = neof)! 42
                       !AnnualClimFN_k=M2AnnualClimFn, nEOFs_k = 10)! 42
    IF( LFAIL ) THEN
       WRITE(*,*) modulename//TRIM(msg)//'in M2CLIM_Stage'
       STOP 1
    ENDIF

    z1 = 2
    z2 = nlyr_pt
    IF (z2-z1 + 1 /= nlyrs0 ) THEN 
      WRITE(*,*) TRIM(ADJUSTL(modulename))//'z2-z1+1 /= nlyrs0'
    ENDIF 

    IF (do_debug) THEN 
      WRITE(*,*) modulename //'| ymd(1:3) =', ymd(1:3)
      WRITE(*,*) modulename //"| which_var, the_var, noef: ", which_var, the_var, neof
    ENDIF
    SET = M2CLIM_stVar( which_var, ymd(1:3), Longitude, Latitude,  &
                        the_var, nlyrs, zstar(1:nlyrs), &
                        lyrOz_Col(1:nlyrs), LFAIL, msg, &
                  dO3prf_dVar_k = dO3prf_dVar(1:nlyrs), &
                        Sa_k = Sap(1:nlyrs,1:nlyrs) )
     
    IF (.NOT. SET) THEN 
      ymd(2) = 0
      SET = M2CLIM_stVar( which_var, ymd(1:3), Longitude, Latitude,  &
                        the_var, nlyrs, zstar(1:nlyrs), &
                        lyrOz_Col(1:nlyrs), LFAIL, msg, &
                  dO3prf_dVar_k = dO3prf_dVar(1:nlyrs), &
                        Sa_k = Sap(1:nlyrs,1:nlyrs) )
      WRITE(*,*) modulename //': Anual mean is set'
    ENDIF
    IF (.NOT. SET) THEN 
      WRITE(*,*) modulename //': Errors in M2CLIM)stVar' ; STOP 1
    ENDIF
   
    ! copy to global  output variables (m2du)
    IF (allocated(m2du%o3p)) deallocate (m2du%o3p, m2du%Sa)
    allocate (m2du%o3p(nlyrs0), m2du%Sa(nlyrs0, nlyrs0))

    lyrOz_Col(1:nlyrs0) = lyrOz_Col(z1:z2)
    Sap(1:nlyrs0, 1:nlyrs0) = Sap(z1:z2, z1:z2)
    IF (do_sort) THEN 
      lyrOz_Col(1:nlyrs0) = lyrOz_Col(ord(1:nlyrs0))
      DO ii = 1, nlyrs0
        Sap(ii, 1:nlyrs0) = Sap(z1+ii-1, ord(z1:z2))
      ENDDO
    ENDIF
    m2du%o3p(:)=lyrOz_Col(1:nlyrs0)
    m2du%Sa(:,:) = Sap(1:nlyrs0,1:nlyrs0)
    deallocate (plev, zstar, lyrOz_col, dO3prf_dVar, Sap)
    IF (do_debug) print * , M2ClimPathFn
  END SUBROUTINE get_m2prof

END MODULE m_get_m2prof
