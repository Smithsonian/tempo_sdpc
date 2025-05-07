!!****************************************************************************
!!F90
!
!!Description:
!
!  MODULE PGS_MET_class
!  contains the defintion of PGS functions that deal with the reading and 
!  writing the ECS metadata. and the data structure definintion for a L2
!  ECS metadata.
! 
! read in tables for calculation of forward model quantities:
! dN/dX, and dN/dT.
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
! email: Kai.Yang@gsfc.nasa.gov
! 
!!Design Notes
!
!!END
!!****************************************************************************
MODULE PGS_MET_class
    IMPLICIT NONE
    INCLUDE 'PGS_MET.f'

    INTEGER (KIND=4), EXTERNAL :: PGS_MET_GetPCAttr_s, &
                                  PGS_MET_GetPCAttr_i, &
                                  PGS_MET_GetPCAttr_r, &
                                  PGS_MET_GetPCAttr_d
    INTEGER (KIND=4), EXTERNAL :: PGS_MET_init, &
                                  PGS_MET_write, &
                                  PGS_MET_sfstart
    INTEGER (KIND=4), EXTERNAL :: PGS_MET_setAttr_i, &
                                  PGS_MET_setAttr_r, &
                                  PGS_MET_setAttr_d, &
                                  PGS_MET_setAttr_s, &
                                  PGS_MET_getSetAttr_s, &
                                  PGS_MET_getSetAttr_i, &
                                  PGS_MET_SetMultiAttr_s, &
                                  PGS_MET_SetMultiAttr_i,  &
                                  PGS_MET_SetMultiAttr_r,  &
                                  PGS_MET_SetMultiAttr_d,  &
                                  PGS_MET_SFend, &
                                  PGS_MET_Remove
END MODULE PGS_MET_class
