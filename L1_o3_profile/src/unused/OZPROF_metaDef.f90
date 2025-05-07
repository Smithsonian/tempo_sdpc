!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE PROFOZ_metaDef
! 
! This module contains the definitions for some meta data items that are 
! specfic to the PROFOZ output. These items should be defined in the 
! FMSO2 MCF file with DATA_LOCATION set to PGE. Each item contains the name,
! the data type, the group (core=2 or archived=3) it belongs to, the number of
! values, and the values stored in a character array.
!
!!Revision History:
! Initial version 02/16/2005  Kai Yang/UMBC
!
!!Team-unique Header:
! This software was developed by the OMI Science Team Support
! Group for the National Aeronautics and Space Administration, Goddard
! Space Flight Center, under NASA Task 916-003-1
!
!!References and Credits
! Written by 
! Kai Yang 
! GEST/UMBC
! email: Kai.Yang.1@gsfc.nasa.gov
! 
!!Design Notes
!
!!END
!!****************************************************************************
MODULE PROFOZ_metaDef
    USE OMI_metaData_class
    IMPLICIT NONE

    TYPE (ECSMETA_ITEM_T) :: it_ParameterName  =             &
         ECSMETA_ITEM_T( "ParameterName.1             ",     &
                         "STRING                      ",     &
                          2, 1                         ,     &
                         "Ozone vertical profiles     " )

    TYPE (ECSMETA_ITEM_T) :: it_QAPercentMissingData =       &
         ECSMETA_ITEM_T( "QAPercentMissingData.1      ",     &
                         "INTEGER                     ",     &
                          2, 1                         ,     &
                         "-999                        " )

    TYPE (ECSMETA_ITEM_T) :: it_QAPercentOutofBoundsData =   &
         ECSMETA_ITEM_T( "QAPercentOutofBoundsData.1  ",     &
                         "INTEGER                     ",     &
                          2, 1                         ,     &
                         "-999                        " )

    TYPE (ECSMETA_ITEM_T) :: it_QAPercentCloudCover =        &
         ECSMETA_ITEM_T( "QAPercentCloudCover.1       ",     &
                         "INTEGER                     ",     &
                          2, 1                         ,     &
                         "-999                        " )

    TYPE (ECSMETA_ITEM_T) :: it_AutoQaFlagExpl =             &
         ECSMETA_ITEM_T( "AutomaticQualityFlagExplanation.1   "// &
                         "                                    ",  &
                         "STRING                              "// &
                         "                                    ",  &
                          2, 1                                 ,  &
                         "Set to 'Failed' if processing error "// &
                         "occurred,'Passed' otherwise         " )

    TYPE (ECSMETA_ITEM_T) :: it_AutomaticQualityFlag  = &
         ECSMETA_ITEM_T( "AutomaticQualityFlag.1      ",     &
                         "STRING                      ",     &
                          2, 1                         ,     &
                         "Passed                      " )

!! The PSAs below are not currently used by OMSO2. They may be implemented
!! at a later date.
!    TYPE (ECSMETA_ITEM_T) :: it_VolcanicActivityFlag  = &
!         ECSMETA_ITEM_T( "VolcanicActivityFlag.1      ",     &
!                         "STRING                      ",     &
!                          3, 1                         ,     &
!                         "None                        " )

!    TYPE (ECSMETA_ITEM_T) :: it_PollutionActivityFlag = &
!         ECSMETA_ITEM_T( "PollutionActivityFlag.1     ",     &
!                         "STRING                      ",     &
!                          3, 1                         ,     &
!                         "None                        " )

!    TYPE (ECSMETA_ITEM_T) :: it_ClearSceneFlag =             &
!      ECSMETA_ITEM_T( "ClearSceneFlag.1            ",     &
!                         "STRING                      ",     &
!                          3, 1                         ,     &
!                         "None                        " )

!    TYPE (ECSMETA_ITEM_T) :: it_StartVolLatitude =                     &
!         ECSMETA_ITEM_T( "StartVolLatitude.1                    ",     &
!                         "DOUBLE                                ",     &
!                          3, 20                                  ,     &
!                         "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" )

!    TYPE (ECSMETA_ITEM_T) :: it_EndVolLatitude   =                     &
!         ECSMETA_ITEM_T( "EndVolLatitude.1                      ",     &
!                         "DOUBLE                                ",     &
!                          3, 20                                  ,     &
!                         "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" )

!    TYPE (ECSMETA_ITEM_T) :: it_StartPollLatitude =                    &
!         ECSMETA_ITEM_T( "StartPollLatitude.1                   ",     &
!                         "DOUBLE                                ",     &
!                          3, 20                                  ,     &
!                         "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" )

!    TYPE (ECSMETA_ITEM_T) :: it_EndPollLatitude   =                    &
!         ECSMETA_ITEM_T( "EndPollLatitude.1                     ",     &
!                         "DOUBLE                                ",     &
!                          3, 20                                  ,     &
!                         "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0" )

!    TYPE (ECSMETA_ITEM_T) :: it_AverageCloudCover =                    &
!         ECSMETA_ITEM_T( "AverageCloudCover.1                   ",     &
!                         "DOUBLE                                ",     &
!                          3, 1                                   ,     &
!                         "0                                      " )

END MODULE PROFOZ_metaDef
