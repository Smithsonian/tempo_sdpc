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
! email: Kai.Yang-1@nasa.gov
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

    INTEGER (KIND=4), PARAMETER :: NL1Bpsa_MAX = 30
    INTEGER (KIND=4), PARAMETER :: Nelm_MAX    = 20
    INTEGER (KIND=4), PARAMETER :: NPMAX       = 10
    INTEGER (KIND=4), PARAMETER :: nZones      = 5
    INTEGER (KIND=4), PARAMETER :: NRing_MAX   = 20
    INTEGER (KIND=4), PARAMETER :: NRPT        = 6

    TYPE, PUBLIC :: L1BECSMETA_T
      !! Core MetaData String Field
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: RangeBeginningDate,   &
                                                  RangeBeginningTime,   &
                                                  RangeEndingDate,      &
                                                  RangeEndingTime,      &
                                                  EquatorCrossingTime,  &
                                                  EquatorCrossingDate,  &
                                                  ShortName,            &
                                                  AutomaticQualityFlag, &
                                                  ScienceQualityFlag

      !! value get from L1B to be compared with those in L2 mcf
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: AssociatedPlatformSN,   &
                                                  AssociatedInstrumentSN, &
                                                  AssociatedSensorSN
      INTEGER (KIND=4) :: nL1Bpsa
      INTEGER (KIND=4) :: nL1BGring

      
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(NL1Bpsa_MAX) :: &
                                                L1BpsaNames
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L), &
                     DIMENSION(NL1Bpsa_MAX,Nelm_MAX) ::  L1BpsaValues

      !! Archived MetaData String Field
      CHARACTER(LEN=PGSd_MET_MAX_STRING_SET_L) :: AlgorithmBypassList,      &
                                                  ProcessingMode,           &
                                                  OrbitData,                &
                                                  ProcessingCenter,         &
                                                  LongName,                 &
                                                  ESDTDescriptorRevision,   &
                                                  OPFSmearSwitchValue,      &
                                                  OPFMeasurementStrayFlag,  &
                                                  OPFVersion,               &
                                                  OPFValid

      ! Integers and double to be read from  core L1B
      INTEGER(KIND=4) :: orbitNumber, QAPercentMissingData
      REAL(KIND=8) :: EqCrossLon

      REAL(KIND=8), DIMENSION(NRing_MAX,NRPT)    :: GRingPointLatitude, &
                                                    GRingPointLongitude
      INTEGER(KIND=4), DIMENSION(NRing_MAX,NRPT) :: GRingPointSequenceNo
      CHARACTER(LEN=1), DIMENSION(NRing_MAX)     :: ExclusionGRingFlag

      ! Real to be read from Archived L1B
      REAL (KIND=8) :: SpacecraftMinAltitude, SpacecraftMaxAltitude

    END TYPE L1BECSMETA_T

    TYPE, PUBLIC :: L2PARAM_T
      CHARACTER (LEN = PGSd_MET_MAX_STRING_SET_L) :: Name
      INTEGER (KIND=4) :: QAPercentMissingData         , &
                          QAPercentHighQualityData     
      INTEGER (KIND=4) :: NumberOfInputSamples         , &
                          NumberOfGoodInputSamples     , &
                          NumberOfLargeSZAInputSamples , &
                          NumberOfMissingInputSamples  , &
                          NumberOfBadInputSamples      , &
                          NumberOfInputWarningSamples  , &
                          NumberOfGoodOutputSamples    , &
                          NumberOfGlintCorrectedSamples, &
                          NumberOfSkippedSamples       , &
                          NumberOfStep1InvalidSamples  , &
                          NumberOfStep2InvalidSamples
      INTEGER (KIND=4) :: NumberOfIrradianceMissing,     &
                          NumberOfIrradianceError,       &
                          NumberOfIrradianceWarning,     &
                          NumberOfRadianceMissing,       &
                          NumberOfRadianceError,         &
                          NumberOfRadianceWarning
      INTEGER (KIND=4) :: NumberOfMeasurement,           &
                          NumberOfMeasurementMissing,    &
                          NumberOfMeasurementError,      &
                          NumberOfMeasurementWarning,    &
                          NumberOfMeasurementRebinned,   &
                          NumberOfMeasurementSAA,        &
                          NumberOfMeasurementManeuver,   &
                          NumberOfInstrumentSettingsError
      INTEGER (KIND=4) :: SolarProductMissing,           &
                          BackupSolarProductUsed,        &
                          SolarProductOutOfDate
      REAL (KIND=8), DIMENSION(nZones) :: ZonalOzoneMin, ZonalOzoneMax 
      REAL (KIND=8), DIMENSION(2,nZones) :: ZonalLatRange
      INTEGER (KIND=4),DIMENSION(4,16) :: QualityFlagsCounters 
    END TYPE L2PARAM_T
END MODULE PGS_MET_class
