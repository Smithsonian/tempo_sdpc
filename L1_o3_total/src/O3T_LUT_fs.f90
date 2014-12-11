!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_LUT_fs
!
! list the HDF5 data fields (including names and dimensions) to be read 
! from the Look-Up Table. 
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
!
!!Return
! None
!
!!Revision History:
! Initial version 03/26/2002  Kai Yang/UMBC
!
!!Team-unique Header:
! This software was developed by the OMI Science Team Support
! Group for the National Aeronautics and Space Administration, Goddard
! Space Flight Center, under NASA Task 916-003-1
!
!!References and Credits
! Written by
! Kai Yang
! University of Maryland Baltimore County
! email: Kai.Yang-1@nasa.gov
!
!!Design Notes
!
!!END
!!****************************************************************************

MODULE O3T_LUT_fs
    USE UTIL_lh5_class
    IMPLICIT NONE

    TYPE (DSH5_T) :: ds_RingAdded  = DSH5_T("RingAdded", -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_dx   = DSH5_T("lyrDeltaO3", -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_sza  = DSH5_T("SolarZenith", -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_vza  = DSH5_T("SensorZenith",-1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_pres = DSH5_T("pressure",    -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_wl   = DSH5_T("wlen",        -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_lgi0 = DSH5_T("lgi0",        -1, -1, &
                               -1,4,-1,5,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_z1   = DSH5_T("z1",          -1, -1, &
                               -1,4,-1,5,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_z2   = DSH5_T("z2",          -1, -1, &
                               -1,4,-1,5,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_tr   = DSH5_T("tr",          -1, -1, &
                               -1,4,-1,5,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_knb  = DSH5_T("knr2",        -1, -1, &
                               -1,4,-1,5,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_sb   = DSH5_T("sb",          -1, -1, &
                               -1,4,-1,3,(/-1,-1,-1,-1,-1,-1,-1/) )
!   TYPE (DSH5_T) :: ds_knb  = DSH5_T("knb",         -1, -1, &
!                              -1,4,-1,3,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_f0   = DSH5_T("f0flux",      -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )

    TYPE (DSH5_T) :: ds_c0   = DSH5_T("c0",        -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_c1   = DSH5_T("c1",        -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
    TYPE (DSH5_T) :: ds_c2   = DSH5_T("c2",        -1, -1, &
                               -1,4,-1,1,(/-1,-1,-1,-1,-1,-1,-1/) )
END MODULE O3T_LUT_fs
