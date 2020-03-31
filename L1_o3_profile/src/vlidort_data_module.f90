MODULE vlidort_data_module
   USE VLIDORT_PARS, ONLY: maxstokes_Sq, maxstokes, fpk, max_atmoswfs, maxlayers

   USE GEMSTOOL_PARS_m
   USE GEMSTOOL_Input_Types_m
   USE GEMSTOOL_Input_types_m, ONLY: GEMSTOOL_Config_Inputs
   USE GEMSTOOL_Geophys_types_m, ONLY: GEMSTOOL_Geophys_Setups
   USE GEMSTOOL_L_Geophys_types_m, ONLY: GEMSTOOL_L_Geophys_Setups

   IMPLICIT NONE

   !---------------------------------------------
   !  GEMST TOOL INPUT
   !----------------------------------------------
   TYPE(GEMSTOOL_Config_Inputs) :: Inputs
   TYPE(GEMSTOOL_Geophys_Setups) :: Geophys
   TYPE(GEMSTOOL_L_Geophys_Setups) :: L_Geophys

   INTEGER :: which_win
  
   TYPE  RT_Inputs_vars
     ! n_totalatmos_wfs, maxlyaers, maxwav
     integer :: ngreek_moments_input, n_totalatmos_wfs, n_surface_wfs
     !logical, dimension (gt_maxlayers) :: layer_vary_flag
     !integer, dimension (gt_maxlayers) :: layer_vary_number
     !character(len=31), dimension (GT_maxatmoswfs) :: profilewf_names
   END TYPE RT_Inputs_vars
   TYPE( RT_Inputs_vars) :: RI

   !integer :: ngreek_moments_input, n_totalatmos_wfs, n_surface_wfs
   
   integer :: ngreek_moments_input, n_totalatmos_wfs, n_surface_wfs
   ! gc variables
   LOGICAL, DIMENSION (gt_maxlayers) :: layer_vary_flag_cc
   INTEGER, DIMENSION (gt_maxlayers) :: layer_vary_number_cc
   CHARACTER(len=31), DIMENSION (gt_maxatmoswfs) :: profilewf_names_cc
 
   ! vlidort variables  
   REAL (KIND=fpk), DIMENSION (:), ALLOCATABLE :: deltau_vert_input, &
   deltaug_vert_input, omega_total_input, fr, fa, fc
   REAL (KIND=fpk), DIMENSION (:, :, :), ALLOCATABLE :: greekmat_total_input
   REAL (KIND=fpk), DIMENSION (:,:), ALLOCATABLE :: l_deltau_vert_input, &
   l_omega_total_input, l_fr, l_fa, l_fc
   REAL (KIND=fpk), DIMENSION (:, :, :, :), ALLOCATABLE :: l_greekmat_total_input
   LOGICAL, DIMENSION (maxlayers) :: layer_vary_flag
   INTEGER, DIMENSION (maxlayers) :: layer_vary_number
   CHARACTER (LEN=31), DIMENSION (max_atmoswfs) :: profilewf_names
   
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! VLIDORT ReSults
   TYPE RT_Output
      real(kind=fpk), ALLOCATABLE :: Stokes(:,:)
      real(kind=fpk), ALLOCATABLE :: Pjac(:,:,:,:)
      real(kind=fpk), ALLOCATABLE :: Sjac(:,:,:)
   END TYPE RT_Output
   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   TYPE(RT_Output) :: RO

  PRIVATE :: LIDORT_Initialize_ConfigInputs

 Contains
 
 Subroutine LIDORT_Read_Config (Configfile, &
                                aerosols, clouds,  &
                                GEMSTOOL_INPUTS, fail, message) 

  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: Configfile
  LOGICAL, INTENT (IN) :: aerosols, clouds
  LOGICAL, INTENT (OUT) :: fail
  TYPE(GEMSTOOL_Config_Inputs), INTENT(INOUT) :: GEMSTOOL_INPUTS
  CHARACTER (LEN=*), INTENT(OUT) :: message
  ! locals
  INTEGER, PARAMETER :: vunit = 1


! Initalize
  fail=.FALSE. ; message = ''
  CALL Lidort_Initialize_ConfigInputs (GEMSTOOL_INPUTS)

! Read
  OPEN(vunit, file=trim(Configfile), status='old')

!==================================================
!1. Atmospheric Constiuents
!===================================================
! Atmosphere.cfg
   read (vunit,*) 
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Tracegases         ! turn on, or off do trace gases
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_aerosols           ! turn on, or off do aerosols calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols       ! turn on, or off do aerosols Mie calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols      ! turn on, or off do aerosols T-matrix calculation
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_User_aerosols      ! turn on, or off User-defined aerosols
   read(1,*)GEMSTOOL_INPUTS%Atmosph%do_clouds             ! turn on, or off General Cloud flag
   read(1,*)GEMSTOOL_INPUTS%Atmosph%Tshift_Ref 

   !Aerosol flags, Checking expanded to include User option. 10/20/16
   !IF (GEMSTOOL_INPUTS%Atmosph%do_aerosols /= aerosols .or. GEMSTOOL_INPUTS%Atmosph%do_clouds /= clouds) THEN
   !   message = 'Inconsistency between VLIDORT v2.7 and PROFOZ for aerosols and clouds'
   !   fail = .true. ;  return
   !ENDIF
   if (GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols.and.GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols) then
      message = 'Cannot have Mie and Tmatrix aerosols, error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif
   if (GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols.and.GEMSTOOL_INPUTS%Atmosph%do_User_aerosols) then
      message = 'Cannot have Mie and User-defined aerosols, error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif
   if (GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols.and.GEMSTOOL_INPUTS%Atmosph%do_User_aerosols) then
      message = 'Cannot have Tmat and User-defined aerosols, error in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif
   if ( GEMSTOOL_INPUTS%Atmosph%do_aerosols.and. &
              ( .not.GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols  .and. &
                .not.GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols .and. &
                .not.GEMSTOOL_INPUTS%Atmosph%do_User_aerosols ) ) then
      message = 'One of Mie/Tmatrix/User aerosol must be set when general aerosol flag set '//&
                'in GEMSTOOL_Atmosphere.cfg'
      fail = .true. ;  return
   endif

!. PCA Controls
   read(1,*)
   Read(1,*)GEMSTOOL_INPUTS%PCAControl%PCA_Strategy_index
   read(1,*)GEMSTOOL_INPUTS%PCAControl%PCA_Binning_index
   read(1,*)GEMSTOOL_INPUTS%PCAControl%PCA_nbins
   if ( GEMSTOOL_INPUTS%PCAControl%PCA_nbins.gt. GT_maxbins ) then
      message = 'PCA Control: Number of Bins exceeds Maximum allowed dimension'
      fail = .true. ; return
   endif
   read(1,*)GEMSTOOL_INPUTS%PCAControl%PCA_neofs(0:GEMSTOOL_INPUTS%PCAControl%PCA_nbins-1)
   read(1,*)GEMSTOOL_INPUTS%PCAControl%do_fast_calculation  ! Default = .true.
   read(1,*)GEMSTOOL_INPUTS%PCAControl%alb_pcainclude       ! Default = .true.,normally!
   read(1,*)GEMSTOOL_INPUTS%PCAControl%do_3M_correction       ! Default = .true.,xcept Fluxes-only
   ! not set for GEMSTOOL_INPUTS%Instruments

!. Radiative transfer control
   read(1,*)
   read(1,*)GEMSTOOL_INPUTS%RTMControl%NVlidort_nstreams     ! Number of Discrete Ordinates in VLIDORT
   read(1,*)GEMSTOOL_INPUTS%RTMControl%NVlidort_nstokes      ! number of stokes parameters in VLIDORT
   read(1,*)GEMSTOOL_INPUTS%RTMControl%do_Upwelling          ! turn on flag for Upwelling
   read(1,*)GEMSTOOL_INPUTS%RTMControl%do_PathRadiance       ! turn on flag for "Path Radiance"      (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTMControl%do_SphericalAlbedo    ! turn on flag for "Spherical Albedo"   (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTMControl%do_2wayTransmittance  ! turn on flag for "2-wayTransmittance" (TBD)
   read(1,*)GEMSTOOL_INPUTS%RTMControl%do_firstorder_option  ! turn on, using new FO code
   read(1,*)GEMSTOOL_INPUTS%RTMControl%do_regular_ps   
! LinControl
   read(1,*) 
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_GasProfile_Jacobians          ! 1
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerOpdepProfile_Jacobians     ! 2
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_AerBulk_Jacobians             ! New,June 2014
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_Surface_Jacobians             ! 4
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian               ! 5
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian            ! 6
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_H2OScaling_Jacobian           ! New July 2014
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_CH4Scaling_Jacobian           ! New Feb 2015
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_SIF_Jacobians                 ! New Oct/Nov 2016
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_normalized_wfoutput           ! 9
   read(1,*)GEMSTOOL_INPUTS%LinControl%do_hitran                        ! New Feb 2015
   if (GEMSTOOL_INPUTS%LinControl%do_GasProfile_Jacobians.or.GEMSTOOL_INPUTS%LinControl%do_AerOpdepProfile_Jacobians) then
      if (  GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian     .or. &
            GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian  .or. &
            GEMSTOOL_INPUTS%LinControl%do_H2OScaling_Jacobian .or. &
            GEMSTOOL_INPUTS%LinControl%do_CH4Scaling_Jacobian .or. &
            GEMSTOOL_INPUTS%LinControl%do_AerBulk_Jacobians )  then

      print *,  GEMSTOOL_INPUTS%LinControl%do_Tshift_Jacobian , &
            GEMSTOOL_INPUTS%LinControl%do_SurfPress_Jacobian  , &
            GEMSTOOL_INPUTS%LinControl%do_H2OScaling_Jacobian , &
            GEMSTOOL_INPUTS%LinControl%do_CH4Scaling_Jacobian , &
            GEMSTOOL_INPUTS%LinControl%do_AerBulk_Jacobians

         message = 'Cannot have Column and Profile Jacobians together - change the flags'
         fail = .true. ; return
      endif
   endif


!  Wavelengths
   ! exact calculations =>  not. do_monochromatic, calculated at coarse intervals and then interpolated onto OMI instruments
   ! PCA calculations   ==> do_monochromatic  and then interpolated onto OMI instruments
   read(1,*)
   read(1,*)GEMSTOOL_INPUTS%Lambdas%do_Monochromatic        ! Flag for standard calculation 10/26/16
   read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_start            ! Initial wavelength
   read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_finish           ! Final wavelength
   read(1,*)GEMSTOOL_INPUTS%Lambdas%lambda_wavres 

!  TIMEPOS
   read(1,*)
   read(1,*) GEMSTOOL_INPUTS%TimePos%earthradius  
!=================================
   close (vunit)
   return 
 
end subroutine LIDORT_Read_Config 


subroutine LIDORT_Initialize_ConfigInputs ( GEMSTOOL_INPUTS )

!  routine has been upgraded to include new type structures
!  Also, everything is now properly initialized

   implicit none

! GEOMTOOL inputs, structure. Intent(out) here

   TYPE(GEMSTOOL_Config_Inputs), INTENT(INOUT) :: GEMSTOOL_INPUTS

!  1. Atmospheric constituents
!  ---------------------------

   GEMSTOOL_INPUTS%Atmosph%do_Tracegases    = .false.
   GEMSTOOL_INPUTS%Atmosph%do_aerosols      = .false.
   GEMSTOOL_INPUTS%Atmosph%do_Mie_aerosols  = .false.
   GEMSTOOL_INPUTS%Atmosph%do_Tmat_aerosols = .false.
   GEMSTOOL_INPUTS%Atmosph%do_User_aerosols = .false.
   GEMSTOOL_INPUTS%Atmosph%do_clouds        = .false.

!  2. Wavelengths and Wavenumbers
   GEMSTOOL_INPUTS%Lambdas%lambda_start  = GTZERO
   GEMSTOOL_INPUTS%Lambdas%lambda_finish = GTZERO
   GEMSTOOL_INPUTS%Lambdas%lambda_wavres = GTZERO

   GEMSTOOL_INPUTS%Wavenums%wavnum_start  = GTZERO
   GEMSTOOL_INPUTS%Wavenums%wavnum_finish = GTZERO
   GEMSTOOL_INPUTS%Wavenums%wavnum_res    = GTZERO

!  3. Radiative transfer control
!  -----------------------------

   GEMSTOOL_INPUTS%RTMControl%NVlidort_nstreams     = 0
   GEMSTOOL_INPUTS%RTMControl%NVlidort_nstokes      = 0
   GEMSTOOL_INPUTS%RTMControl%do_Upwelling          = .false.
   GEMSTOOL_INPUTS%RTMControl%do_PathRadiance       = .false.
   GEMSTOOL_INPUTS%RTMControl%do_SphericalAlbedo    = .false.
   GEMSTOOL_INPUTS%RTMControl%do_2wayTransmittance  = .false.
   GEMSTOOL_INPUTS%RTMControl%do_firstorder_option  = .false.
   GEMSTOOL_INPUTS%RTMControl%do_regular_ps         = .false.
!  4. Time/Location
!  ----------------

!  Hour/Minute/Second  variables added for  specification of SIF Epoch
!    10/18/16. Rob Spurr

   GEMSTOOL_INPUTS%TimePos%latitude     = GTZERO
   GEMSTOOL_INPUTS%TimePos%longitude    = GTZERO
   GEMSTOOL_INPUTS%TimePos%earthradius  = GTZERO
   GEMSTOOL_INPUTS%TimePos%year         = 0
   GEMSTOOL_INPUTS%TimePos%month        = 0
   GEMSTOOL_INPUTS%TimePos%day_of_month = 0
   GEMSTOOL_INPUTS%TimePos%Hour         = 0
   GEMSTOOL_INPUTS%TimePos%Minute       = 0
   GEMSTOOL_INPUTS%TimePos%second       = 0

!  5. Trace gas information
   GEMSTOOL_INPUTS%Tracegas%ngases      = 0
   GEMSTOOL_INPUTS%Tracegas%which_gases = '    '
   GEMSTOOL_INPUTS%Tracegas%do_gases    = .false.
   GEMSTOOL_INPUTS%Tracegas%do_gas_wfs  = .false.
   GEMSTOOL_INPUTS%Tracegas%Gas_Profile_names = ' '

   GEMSTOOL_INPUTS%Tracegas%do_H2OScaling = .false.
   GEMSTOOL_INPUTS%Tracegas%H2OScaling    = GTZERO

   GEMSTOOL_INPUTS%Tracegas%do_CH4Scaling = .false.
   GEMSTOOL_INPUTS%Tracegas%CH4Scaling    = GTZERO
!  --------------

!  Rob, 9/5/17.

   GEMSTOOL_INPUTS%PCAControl%do_fast_calculation = .false.
   GEMSTOOL_INPUTS%PCAControl%alb_pcainclude      = .false.
   GEMSTOOL_INPUTS%PCAControl%do_3M_correction    = .false.

   GEMSTOOL_INPUTS%PCAControl%PCA_Strategy_index = 0
   GEMSTOOL_INPUTS%PCAControl%PCA_Binning_index  = 0
   GEMSTOOL_INPUTS%PCAControl%PCA_nbins = 0
   GEMSTOOL_INPUTS%PCAControl%PCA_neofs = 0

  end subroutine LIDORT_Initialize_ConfigInputs

! Jbak
  SUBROUTINE allocate_od(nz, nwf)
    INTEGER, INTENT(IN) :: nz , nwf
    ! optical depth input
    allocate (deltau_vert_input(nz), deltaug_vert_input(nz),omega_total_input(nz))    
    allocate (fr(nz), fa(nz), fc(nz))
    ! linearizaed optical depth input
    allocate (l_deltau_vert_input(nwf,nz),l_omega_total_input(nwf, nz))
    allocate (l_fr(nwf,nz),l_fa(nwf,nz),l_fc(nwf,nz))
    RETURN   
  END SUBROUTINE allocate_od

  SUBROUTINE deallocate_od
   deallocate (deltau_vert_input, omega_total_input, deltaug_vert_input, fr, fa,fc)
   deallocate (l_deltau_vert_input, l_omega_total_input,l_fr, l_fa, l_fc)
   RETURN
  END SUBROUTINE deallocate_od

END MODULE
