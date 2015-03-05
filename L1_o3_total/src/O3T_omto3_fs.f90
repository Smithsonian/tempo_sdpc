!!****************************************************************************
!!F90
!
!!Description:
!
! MODULE O3T_omto3_fs
! 
! This module contains the definitions for the HE5 output of OMTO3 PGE. 
! Each geo or data field listed in this module contians the the name, 
! data type, dimeninsion names, rank, valid range, scale, offset, title, unit,
! and UniqueFieldDefinition for the field to be written out to the output file.
! To change the output geo and data field, one only need to change the 
! definition in this file. 
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
MODULE O3T_omto3_fs
    USE HE5_class
    IMPLICIT NONE

    INTEGER (KIND=4), PARAMETER :: NG_r4 = 5, NG_i2 = 2
    INTEGER (KIND=4), PARAMETER :: NG_i1 = 1
    INTEGER (KIND=4), PARAMETER :: NG_omto3 = NG_r4+NG_i2+NG_i1+1
    INTEGER (KIND=4), PARAMETER :: ND_2D_i1 = 1, ND_2D_i2 = 12, &
                                   ND_3Dwl  = 6, ND_3Dlyr = 2
    INTEGER (KIND=4), PARAMETER :: ND_2D = ND_2D_i1 + ND_2D_i2 
    INTEGER (KIND=4), PARAMETER :: ND_omto3 = ND_2D + ND_3Dwl + ND_3Dlyr

! Geolocation Fields
    TYPE (DFHE5_T) :: gf_GroundPixelQualityFlags  =                       &
           DFHE5_T( 0.0, 65534.0D0, 1.0D0, 0.0D0,                         &
                    "GroundPixelQualityFlags      ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Ground Pixel Quality Flags   ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_UINT16, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_Latitude =                                       &
           DFHE5_T( -90.0, 90.0, 1.0D0, 0.0D0,                            &
                    "Latitude                     ",                      &
                    "nXtrack,nTimes               ",                      & 
                    "deg                          ",                      &
                    "Geodetic Latitude            ",                      &
                    "TOMS-Aura-Shared             ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_Longitude =                                      &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "Longitude                    ",                      &
                    "nXtrack,nTimes               ",                      & 
                    "deg                          ",                      &
                    "Geodetic Longitude           ",                      &
                    "TOMS-Aura-Shared             ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SolarAzimuthAngle =                              &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "SolarAzimuthAngle                              ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg(EastofNorth)                               ",    &
                    "Solar Azimuth Angle                            ",    &
                    "TOMS-Aura-Shared                               ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_ViewingAzimuthAngle =                            &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "ViewingAzimuthAngle                            ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg(EastofNorth)                               ",    &
                    "Viewing Azimuth Angle                          ",    &
                    "TOMS-Aura-Shared                               ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_RelativeAzimuthAngle =                           &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "RelativeAzimuthAngle                           ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg(EastofNorth)                               ",    &
                    "Relative Azimuth Angle (sun + 180 - view)      ",    &
                    "TOMS-OMI-Shared                                ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SolarZenithAngle =                               &
           DFHE5_T( 0.0, 180.0, 1.0D0, 0.0D0,                             &
                    "SolarZenithAngle                               ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg                                            ",    &
                    "Solar Zenith Angle                             ",    &
                    "TOMS-Aura-Shared                               ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_ViewingZenithAngle =                             &
           DFHE5_T( 0.0, 70.0, 1.0D0, 0.0D0,                              &
                    "ViewingZenithAngle                             ",    &
                    "nXtrack,nTimes                                 ",    & 
                    "deg                                            ",    &
                    "Viewing Zenith Angle                           ",    &
                    "TOMS-OMI-Shared                                ",    & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_TerrainHeight  =                                 &
           DFHE5_T( -200.0, 10000.0, 1.0D0, 0.0D0,                        & 
                    "TerrainHeight                ",                      &
                    "nXtrack,nTimes               ",                      &
                    "m                            ",                      &
                    "Terrain Height               ",                      & 
                    "TOMS-Aura-Shared             ",                      & 
                     -1, HE5T_NATIVE_INT16, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: gf_Time  =                                          &
           DFHE5_T( -5.0D09, 1.0D10, 1.0D0, 0.0D0,                        & 
                    "Time                            ",                   &
                    "nTimes                          ",                   &
                    "s                               ",                   &
                    "Time at Start of Scan (TAI93)   ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_DOUBLE, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SecondsInDay =                                   &
           DFHE5_T( 0.0, 86401.0, 1.0D0, 0.0D0,                           & 
                    "SecondsInDay                    ",                   &
                    "nTimes                          ",                   &
                    "s                               ",                   &
                    "Seconds after UTC midnight      ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SpacecraftLatitude =                             &
           DFHE5_T( -90.0, 90.0, 1.0D0, 0.0D0,                            & 
                    "SpacecraftLatitude              ",                   &
                    "nTimes                          ",                   &
                    "deg                             ",                   &
                    "Spacecraft Latitude             ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SpacecraftLongitude =                            &
           DFHE5_T( -180.0, 180.0, 1.0D0, 0.0D0,                          & 
                    "SpacecraftLongitude             ",                   &
                    "nTimes                          ",                   &
                    "deg                             ",                   &
                    "Spacecraft Longitude            ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_SpacecraftAltitude =                             &
           DFHE5_T( 4.0D05, 9.0D05, 1.0D0, 0.0D0,                         & 
                    "SpacecraftAltitude              ",                   &
                    "nTimes                          ",                   &
                    "m                               ",                   &
                    "Spacecraft Altitude             ",                   & 
                    "TOMS-Aura-Shared                ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: gf_XTrackQualityFlags  =                            &
           DFHE5_T( 0, 254, 1.0D0, 0.0D0,                                 &
                    "XTrackQualityFlags           ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Cross Track Quality Flags    ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_UINT8, 2, (/ -1, -1, -1 /))

! Data Fields
    TYPE (DFHE5_T) :: df_MeasurementQualityFlags  =                       &
           DFHE5_T( 0, 254, 1.0D0, 0.0D0,                                 & 
                    "MeasurementQualityFlags         ",                   &
                    "nTimes                          ",                   &
                    "NoUnits                         ",                   &
                    "Measurement Quality Flags       ",                   & 
                    "TOMS-OMI-Shared                 ",                   & 
                     -1, HE5T_NATIVE_UINT8, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_NumberSmallPixelColumns  =                       &
           DFHE5_T( 0, 5, 1.0D0, 0.0D0,                                   & 
                    "NumberSmallPixelColumns         ",                   &
                    "nTimes                          ",                   &
                    "NoUnits                         ",                   &
                    "Number of Small Pixel Columns   ",                   & 
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_UINT8, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_SmallPixelColumn  =                              &
           DFHE5_T( 0, 1000, 1.0D0, 0.0D0,                                & 
                    "SmallPixelColumn                ",                   &
                    "nTimes                          ",                   &
                    "NoUnits                         ",                   &
                    "Small Pixel Column              ",                   &
                    "OMI-Specific                    ",                   & 
                     -1, HE5T_NATIVE_INT16, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_InstrumentConfigurationId =                          &
           DFHE5_T( 0, 200, 1.0D0, 0.0D0,                                     & 
                    "InstrumentConfigurationId                             ", &
                    "nTimes                                                ", &
                    "NoUnits                                               ", &
                    "Instrument Configuration ID                           ", &
                    "OMI-Specific                                          ", & 
                     -1, HE5T_NATIVE_UINT8, 1, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_Wavelength  =                                    &
           DFHE5_T( 300.0, 400.0, 1.0D0, 0.0D0,                           & 
                    "Wavelength                      ",                   &
                    "nWavel                          ",                   &
                    "nm                              ",                   &
                    "Wavelength                      ",                   & 
                    "TOMS-OMI-Shared                 ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 1, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_TerrainPressure  =                               &
           DFHE5_T( 0.0, 1013.25, 1.0D0, 0.0D0,                           & 
                    "TerrainPressure              ",                      &
                    "nXtrack,nTimes               ",                      &
                    "hPa                          ",                      &
                    "Terrain Pressure             ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_AlgorithmFlags =                                 &
           DFHE5_T( 0.0, 13.0, 1.0D0, 0.0D0,                              & 
                    "AlgorithmFlags               ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Algorithm Flags              ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_UINT8, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_APrioriLayerO3 =                                 &
           DFHE5_T( 0.0, 125.0, 1.0D0, 0.0D0,                             & 
                    "APrioriLayerO3               ",                      &
                    "nLayers,nXtrack,nTimes       ",                      &
                    "DU                           ",                      &
                    "A Priori Ozone Profile       ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_CloudFraction    =                               &
           DFHE5_T( 0.0, 1.0, 1.0D0, 0.0D0,                               & 
                    "RadiativeCloudFraction                   ",          &
                    "nXtrack,nTimes                           ",          &
                    "NoUnits                                  ",          &
                    "Radiative Cloud Fraction = fc*Ic331/Im331",          & 
                    "TOMS-OMI-Shared                          ",          & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
!!                  "Radiative Cloud Fraction = fc*0.8/refl360",          & 

    TYPE (DFHE5_T) :: df_fc    =                                          &
           DFHE5_T( 0.0, 1.0, 1.0D0, 0.0D0,                               & 
                    "fc                                        ",         &
                    "nXtrack,nTimes                            ",         &
                    "NoUnits                                   ",         &
                    "Mixed LER Model (Cloud Fraction) Parameter",         & 
                    "TOMS-OMI-Shared                           ",         & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_CloudPressure  =                                 &
           DFHE5_T( 0.0, 1013.25, 1.0D0, 0.0D0,                           & 
                    "CloudPressure                ",                      &
                    "nXtrack,nTimes               ",                      &
                    "hPa                          ",                      &
                    "Effective Cloud Pressure     ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_ColumnAmountO3    =                              &
           DFHE5_T( 50.0, 700.0, 1.0D0, 0.0D0,                            & 
                    "ColumnAmountO3                ",                     &
                    "nXtrack,nTimes                ",                     &
                    "DU                            ",                     &
                    "Best Total Ozone Solution     ",                     & 
                    "TOMS-OMI-Shared               ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_dNdR =                                           &
           DFHE5_T( -200.0, 0.0, 1.0D0, 0.0D0,                            & 
                    "dN_dR                                ",              &
                    "nWavel,nXtrack,nTimes                ",              &
                    "NoUnits                              ",              &
                    "Reflectivity Sensitivity Ratio, dN/dR",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_QualityFlags   =                                 &
           DFHE5_T( 0.0, 65534.0, 1.0D0, 0.0D0,                           & 
                    "QualityFlags                 ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Quality Flags                ",                      & 
                    "TOMS-OMI-Shared              ",                      & 
                     -1, HE5T_NATIVE_UINT16, 2, (/ -1, -1, -1 /))
                     
   TYPE (DFHE5_T) :: df_RadianceBadPixelFlagAccepted =                    &
           DFHE5_T( 0.0, 65534.0, 1.0D0, 0.0D0,                           & 
                    "RadianceBadPixelFlagAccepted ",                      &
                    "nXtrack,nTimes               ",                      &
                    "NoUnits                      ",                      &
                    "Radiance Bad Pixel Flag Accepted",                   & 
                    "OMI-Specific                 ",                      & 
                     -1, HE5T_NATIVE_UINT16, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_LayerEfficiency =                                &
           DFHE5_T( 0.0, 10.0, 1.0D0, 0.0D0,                              & 
                    "LayerEfficiency                      ",              &
                    "nLayers,nXtrack,nTimes               ",              &
                    "NoUnits                              ",              &
                    "Algorithmic Layer Efficiency         ",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_NValue =                                         &
           DFHE5_T( 0.0, 600.0, 1.0D0, 0.0D0,                             & 
                    "NValue                               ",              &
                    "nWavel,nXtrack,nTimes                ",              &
                    "NoUnits                              ",              &
                    "Measured N-Value                     ",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_O3BelowCloud =                                   &
           DFHE5_T( 0.0, 100.0, 1.0D0, 0.0D0,                             & 
                    "O3BelowCloud                     ",                  &
                    "nXtrack,nTimes                   ",                  &
                    "DU                               ",                  &
                    "Ozone Below Fractional Cloud     ",                  & 
                    "TOMS-OMI-Shared                  ",                  & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_Reflectivity331 =                                 &
           DFHE5_T( -15.0, 115.0, 1.0D0, 0.0D0,                             & 
                    "Reflectivity331                                     ",&
                    "nXtrack,nTimes                                      ",&
                    "%                                                   ",&
                    "Effective Surface Reflectivity at 331 nm            ",&
                    "TOMS-OMI-Shared                                     ",& 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_Reflectivity360 =                                 &
           DFHE5_T( -15.0, 115.0, 1.0D0, 0.0D0,                             & 
                    "Reflectivity360                                     ",&
                    "nXtrack,nTimes                                      ",&
                    "%                                                   ",&
                    "Effective Surface Reflectivity at 360 nm            ",&
                    "TOMS-OMI-Shared                                     ",& 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_Residual =                                       &
           DFHE5_T( -32.0, 32.0, 1.0D0, 0.0D0,                            & 
                    "Residual                             ",              &
                    "nWavel,nXtrack,nTimes                ",              &
                    "NoUnits                              ",              &
                    "N-Value Residual                     ",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_ResidualStep1 =                                  &
           DFHE5_T( -32.0, 32.0, 1.0D0, 0.0D0,                            & 
                    "ResidualStep1                        ",              &
                    "nWavel,nXtrack,nTimes                ",              &
                    "NoUnits                              ",              &
                    "Step 1 N-Value Residual              ",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_ResidualStep2 =                                  &
           DFHE5_T( -32.0, 32.0, 1.0D0, 0.0D0,                            & 
                    "ResidualStep2                        ",              &
                    "nWavel,nXtrack,nTimes                ",              &
                    "NoUnits                              ",              &
                    "Step 2 N-Value Residual              ",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_Sensitivity =                                    &
           DFHE5_T( 0.0, 1.0, 1.0D0, 0.0D0,                               & 
                    "Sensitivity                          ",              &
                    "nWavel,nXtrack,nTimes                ",              &
                    "NoUnits                              ",              &
                    "Ozone Sensitivity Ratio, dN/dOmega   ",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_dNdT =                                           &
           DFHE5_T( 0.0, 1.0, 1.0D0, 0.0D0,                               & 
                    "dN_dT                                              ",&
                    "nWavel,nXtrack,nTimes                              ",&
                    "NoUnits                                            ",&
                    "Ozone Weighted Temperature Sensitivity Ratio, dN/dT",& 
                    "TOMS-OMI-Shared                                    ",& 
                     -1, HE5T_NATIVE_FLOAT, 3, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_SO2index =                                       &
           DFHE5_T( -300.0, 300.0, 1.0D0, 0.0D0,                          & 
                    "SO2index                             ",              &
                    "nXtrack,nTimes                       ",              &
                    "NoUnits                              ",              &
                    "SO2 Index                            ",              & 
                    "TOMS-OMI-Shared                      ",              & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))
                     
    TYPE (DFHE5_T) :: df_StepOneO3 =                                      &
           DFHE5_T( 50.0, 700.0, 1.0D0, 0.0D0,                            & 
                    "StepOneO3                     ",                     &
                    "nXtrack,nTimes                ",                     &
                    "DU                            ",                     &
                    "Step 1 Ozone Solution         ",                     & 
                    "TOMS-OMI-Shared               ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_StepTwoO3 =                                      &
           DFHE5_T( 5.0, 700.0, 1.0D0, 0.0D0,                             & 
                    "StepTwoO3                     ",                     &
                    "nXtrack,nTimes                ",                     &
                    "DU                            ",                     &
                    "Step 2 Ozone Solution         ",                     & 
                    "TOMS-OMI-Shared               ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_UVAerosolIndex =                                 &
           DFHE5_T( -30.0, 30.0, 1.0D0, 0.0D0,                            & 
                    "UVAerosolIndex                ",                     &
                    "nXtrack,nTimes                ",                     &
                    "NoUnits                       ",                     &
                    "UV Aerosol Index              ",                     & 
                    "TOMS-OMI-Shared               ",                     & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

    TYPE (DFHE5_T) :: df_CalibrationAdjustment   =                        &
           DFHE5_T( -10.0, 10.0, 1.0D0, 0.0D0,                            & 
                    "CalibrationAdjustment           ",                   &
                    "nWavel,nXtrack                  ",                   &
                    "NoUnits                         ",                   &
                    "Calibration Adjustment          ",                   & 
                    "TOMS-OMI-Shared                 ",                   & 
                     -1, HE5T_NATIVE_FLOAT, 2, (/ -1, -1, -1 /))

END MODULE O3T_omto3_fs
