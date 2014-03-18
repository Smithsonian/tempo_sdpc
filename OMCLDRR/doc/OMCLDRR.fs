Shortname:    OMCLDRR
Longname:     OMI/Aura Cloud Pressure and Fraction (Raman Scattering) 1-Orbit L2 Swath 13x24km
PSF Version:  1.9.1
Date:         27 March 2012
Author(s):    Peter Leonard (SSAI), Alexander Vassilkov (SSAI), and Mike Linda (SAIC)

PGE Version:                1.9.1
Lead Algorithm Scientist:   Joanna Joiner (NASA/GSFC)
Other Algorithm Scientists: Alexander Vassilkov (SSAI)
PGE Developer:              Alexander Vassilkov (SSAI)

Description: >

  This document specifies the product format for the OMCLDRR Level 2 PGE,
  that derives cloud pressure from rotational Raman scattering in the atmosphere. 
  The product is stored as one HDF-EOS 5 swath file for
  each granule (i.e., one orbit) of OMI Level 1B data, and
  has a size range of 4 to 15 Mb.  

Global Metadata:

 - Metadata Name:     AuthorAffiliation
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       PCF
   Description:       Actual is "NASA/GSFC".

 - Metadata Name:     AuthorName
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       PCF
   Description:       Actual is "U.S. OMI Science Team".

 - Metadata Name:     GranuleDay
   Mandatory:         T
   Data Type:         HE5T_NATIVE_INT
   Number of Values:  1
   Range or Valids:   1 to 31
   Data Source:       PGE
   Description:       The day of the month at the start of the granule.

 - Metadata Name:     GranuleMonth
   Mandatory:         T
   Data Type:         HE5T_NATIVE_INT
   Number of Values:  1
   Range or Valids:   1 to 12
   Data Source:       PGE
   Description:       The month at the start of the granule.

 - Metadata Name:     GranuleYear
   Mandatory:         T
   Data Type:         HE5T_NATIVE_INT
   Number of Values:  1
   Range or Valids:   2003 to 2099
   Data Source:       PGE
   Description:       The (four-digit) year at the start of the granule.

 - Metadata Name:     HDFEOSVersion
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Automatically set by HDF-EOS.
   Data Source:       HE
   Description:       Actual is "HDFEOS_5.1.5".

 - Metadata Name:     InputVersions
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       PGE
   Description: >

     A list of every ESDT (including version) whose product was used as
     input for the processing.

 - Metadata Name:     InstrumentName
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Valids are "HIRDLS", "MLS", "OMI" and "TES".
   Data Source:       PCF
   Description:       Actual is "OMI" (see Section 6.1 of Reference 2).

 - Metadata Name:     OrbitData
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Valids are "DEFINITIVE" and "PREDICTED".
   Data Source:       L1B
   Description: >

     Indicates whether orbit data used by the L1B processor is definitive
     or predicted (not implemented yet because of no OrbitData in synthetic input data).

 - Metadata Name:     PGEVERSION
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Range is "0.0.0" to "9.9.99".
   Data Source:       PCF
   Description:       Actual is "1.9.0" (see Appendix K of Reference 3).

 - Metadata Name:     ProcessingCenter
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       PCF
   Description:       Actual is "OMIDAPS".

 - Metadata Name:     ProcessingHost
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       PCF
   Description: >

     The output from executing the Unix "uname -a" command on the
     processing machine.

 - Metadata Name:     ProcessLevel
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Valids are "1b", "2" and "3".
   Data Source:       PCF
   Description:       Actual is "2".

 - Metadata Name:     TAI93At0zOfGranule
   Mandatory:         T
   Data Type:         HE5T_NATIVE_DOUBLE
   Number of Values:  1
   Range or Valids:   0.0d+00 to 1.0d+30
   Data Source:       PGE
   Description: >

     The TAI93 time at 0z of the granule (see Section 6.1 of Reference 2).

Swath Metadata:

 - Metadata Name:     NumTimes
   Mandatory:         T
   Data Type:         HE5T_NATIVE_INT
   Number of Values:  1
   Range or Valids:   0 to 9999
   Data Source:       L1B
   Description:       The number of "scan" lines in the swath.

 - Metadata Name:     NumTimesSmallPixel
   Mandatory:         T
   Data Type:         HE5T_NATIVE_INT
   Number of Values:  1
   Range or Valids:   0 to 9999
   Data Source:       L1B
   Description:       The number of small pixel "scan" lines in the swath.

 - Metadata Name:     EarthSunDistance 
   Mandatory:         T
   Data Type:         HE5T_NATIVE_FLOAT
   Number of Values:  1
   Range or Valids:   1.0E+11 to 2.0E+11 
   Data Source:       L1B
   Description:       The Earth to Sun distance in meters.

 - Metadata Name:     SwathName
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids:   Valid is "Cloud Product".
   Data Source:       PGE
   Description:       Actual is "Cloud Product".

 - Metadata Name:     VerticalCoordinate
   Mandatory:         T
   Data Type:         HE5T_NATIVE_CHAR
   Number of Values:  1
   Range or Valids: >

     Valids are "Pressure", "Altitude", "Potential Temperature" and
     "Total Column".

   Data Source:       PGE
   Description: >

     Actual is "Total Column" (see Section 6.2 of Reference 2).

Swath Dimensions:

 - Dimension Name:    nTimes
   Data Type:         HE5T_NATIVE_INT
   Dimension Type:    FIXED
   Number of Values:  1
   Range or Valids:   0 to 9999
   Data Source:       L1B
   Description:       The number of "scan" lines in the swath.

 - Dimension Name:    nWavel
   Data Type:         HE5T_NATIVE_INT
   Dimension Type:    FIXED
   Number of Values:  1
   Range or Valids:   1 to 500
   Data Source:       PGE
   Description:       The number of wavelengths per ground pixel.

 - Dimension Name:    nXtrack
   Data Type:         HE5T_NATIVE_INT
   Dimension Type:    FIXED
   Number of Values:  1
   Range or Valids:   1 to 60
   Data Source:       L1B
   Description:       The number of ground pixels per "scan" line.

Geolocation Fields:

 - Field Name:               GroundPixelQualityFlags
   Data Type:                HE5T_NATIVE_UINT16
   Dimensions:               nXtrack,nTimes
   Range or Valids:          Not meaningful.
   Missing Value:            65535
   Units:                    NoUnits
   Data Source:              L1B
   Title:                    "Ground Pixel Quality Flags"
   Unique Field Definition:  OMI-Specific
   Description: >

     Bits 0 to 3 together contain the land/water flags:
       0    - shallow ocean
       1    - land
       2    - shallow inland water
       3    - ocean coastline/lake shoreline
       4    - ephemeral (intermittent) water
       5    - deep inland water
       6    - continental shelf ocean
       7    - deep ocean
       8-14 - not used
       15   - error flag for land/water
     Bits 4 to 6 are flags that are set to 0 for FALSE, or 1 for TRUE:
       Bit 4 - sun glint possibility flag
       Bit 5 - solar eclipse possibility flag
       Bit 6 - geolocation error flag
     Bit 7 is reserved for future use (currently set to 0).
     Bits 8 to 14 together contain the snow/ice flags (based on NISE):
       0       - snow-free land
       1-100   - sea ice concentration (percent)
       101     - permanent ice (Greenland, Antarctica)
       102     - not used
       103     - dry snow
       104     - ocean (NISE-255)
       105-123 - reserved for future use
       124     - mixed pixels at coastline (NISE-252)
       125     - suspect ice value (NISE-253)
       126     - corners undefined (NISE-254)
       127     - error
     Bit 15 - NISE nearest neighbor filling flag.
     (See Section 6.2 of Reference 4 for more details.)

 - Field Name:               Latitude
   Data Type:                HE5T_NATIVE_FLOAT
   Dimensions:               nXtrack,nTimes
   Range or Valids:          Range is -90.0 to 90.0.
   Missing Value:            -9999.0
   Units:                    deg
   Data Source:              L1B
   Title:                    "Latitude"
   Unique Field Definition:  OMI-Specific
   Description: >

     The geodetic latitude (in deg) at the center of the ground pixel.

 - Field Name:               Longitude
   Data Type:                HE5T_NATIVE_FLOAT
   Dimensions:               nXtrack,nTimes
   Range or Valids:          Range is -180.0 to 180.0.
   Missing Value:            -9999.0
   Units:                    deg
   Data Source:              L1B
   Title:                    "Longitude"
   Unique Field Definition:  OMI-Specific
   Description: >

     The geodetic longitude (in deg) at the center of the ground pixel.

 - Field Name:               RelativeAzimuthAngle
   Data Type:                HE5T_NATIVE_FLOAT
   Dimensions:               nXtrack,nTimes
   Range or Valids:          Range is 0.0 to 180.0.
   Missing Value:            -9999.0
   Units:                    deg
   Data Source:              L1B
   Title:                    "Relative Azimuth Angle"
   Unique Field Definition:  OMI-Specific
   Description: >

     The relative (sun + 180 - view) azimuth angle (in deg) at the center
     of the ground pixel. Absolute values of the azimuth angle are given. 

 - Field Name:               SolarZenithAngle
   Data Type:                HE5T_NATIVE_FLOAT
   Dimensions:               nXtrack,nTimes
   Range or Valids:          0.0 to 180.0
   Missing Value:            -9999.0
   Units:                    deg
   Data Source:              L1B
   Title:                    "Solar Zenith Angle"
   Unique Field Definition:  OMI-Specific
   Description: >

     The solar zenith angle (in deg) at the center of the ground pixel.

 - Field Name:               TerrainHeight
   Data Type:                HE5T_NATIVE_INT
   Dimensions:               nXtrack,nTimes
   Range or Valids:          Range is -100 to 10000.
   Missing Value:            65535
   Units:                    m
   Data Source:              L1B
   Title:                    "Terrain Height"
   Unique Field Definition:  Aura-Shared
   Description: >

     The terrain height (in m) at the center of the ground pixel (from
     OMI Level 1B file).

 - Field Name:               Time
   Data Type:                HE5T_NATIVE_DOUBLE
   Dimensions:               nTimes
   Range or Valids:          0.0d+00 to 1.0d+30
   Missing Value:            -9999.0
   Units:                    s
   Data Source:              L1B
   Title:                    "Time at Start of Scan (s, TAI93)"
   Unique Field Definition:  Aura-Shared
   Description: >

     The TAI93 time (in s) at the start of the "scan".

 - Field Name:               ViewingZenithAngle
   Data Type:                HE5T_NATIVE_FLOAT
   Dimensions:               nXtrack,nTimes
   Range or Valids:          0.0 to 90.0
   Missing Value:            -9999.0
   Units:                    deg
   Data Source:              L1B
   Title:                    "Viewing Zenith Angle"
   Unique Field Definition:  OMI-Specific
   Description: >

     The viewing zenith angle (in deg) at the center of the ground pixel.

 - Field Name:               XTrackQualityFlags
   Data Type:                HE5T_NATIVE_UINT8
   Dimensions:               nXtrack,nTimes
   Range or Valids:          Not meaningful.
   Missing Value:            255
   Units:                    NoUnits
   Data Source:              L1B
   Title:                    "Cross Track Quality Flags"
   Unique Field Definition:  OMI-Specific
   Description: >

     Bits 0 to 2 are combined to a single value in the range 0-7.
         value 0 - not affected by row anomaly, pixel can be used
               1 - affected, not corrected, do not use pixel
               2 - slightly affected, not corrected, use with caution
               3 - affected, corrected not optimally, use with caution
               4 - affected, corrected optimally, pixel can be used
               5,6 - not used
               7 - error during correction, do not use pixel
       3    - reserved for future use
       4    - pixel may be affected by wavelength shift
       5    - pixel may be affected by blockage
       6    - pixel may be affected by stray sunlight
       7    - pixel may be affected by stray earthshine

Data Fields:

 - Field Name:        Chlorophyll 
   Data Type:         HE5T_NATIVE_FLOAT 
   Dimensions:        nXtrack,nTimes
   Valid Range:       0.0 to 50.0
   Missing Value:     -9999.0 
   Units:             mg/m3 
   Data Source:       PGE 
   Title:             "Chlorophyll Concentration"
   Unique Field Definition:  OMI-Specific
   Description:       >
    Retrieved chlorophyll concentration.

 - Field Name:        CloudFractionforO3
   Data Type:         HE5T_NATIVE_FLOAT
   Dimensions:        nXtrack,nTimes
   Valid Range:       0.0 to 1.0
   Missing Value:     -9999.0
   Units:             NoUnits
   Data Source:       PGE
   Title:             "Cloud Fraction for O3"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The cloud fraction associated with the ground pixel.
 
 - Field Name:        CloudMask
   Data Type:         HE5T_NATIVE_UINT16
   Dimensions:        nXtrack,nTimes
   Valid Range:       0 to 1
   Missing Value:     -9999.0
   Units:             NoUnits
   Data Source:       PGE
   Title:             "Cloud Mask"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The cloud mask (0=clear; 1=cloudy) associated with the ground pixel.
 
 - Field Name:        CloudPressureforO3
   Data Type:         HE5T_NATIVE_FLOAT
   Dimensions:        nXtrack,nTimes
   Valid Range:       0.0 to 1200.0
   Missing Value:     -9999.0
   Units:             hPa
   Data Source:       PGE
   Title:             "Cloud Pressure for O3"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The effective cloud pressure (in hPa) associated with the ground pixel.
 
 - Field Name:        Convergence_factor 
   Data Type:         HE5T_NATIVE_FLOAT 
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is 0.0 to 10.0
   Missing Value:     -9999.0 
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Convergence Factor"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The chi-squared factor characterizing the fitting procedure.

 - Field Name:        Filling-In 
   Data Type:         HE5T_NATIVE_FLOAT 
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is 0.0 to 0.10
   Missing Value:     -9999.0 
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Effective Filling-in"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The effective filling-in at 352.65 nm.

 - Field Name:        MeasurementQualityFlags
   Data Type:         HE5T_NATIVE_UINT16
   Dimensions:        nTimes
   Valid Range:       Not meaningful
   Missing Value:     65535
   Units:             NoUnits
   Data Source:       L1B
   Title:             "Measurement Quality Flags"
   Unique Field Definition:  OMI-Specific
   Description:       >

     Bits contain several output error and warning flags:
       0    - missing data
       1    - radiance measurement error flag
       2    - radiance measurement warning flag
       3    - rebinned radiance measurement flag
       4    - South Atlantic anomaly flag
       5    - irradiance measurement warning flag
       6-15 - not used yet

 - Field Name:        ProcessingQualityFlagsforO3 
   Data Type:         HE5T_NATIVE_UINT16
   Dimensions:        nXtrack,nTimes
   Valid Range:       Not meaningful
   Missing Value:     65535
   Units:             NoUnits
   Title:             "Processing Quality Flags for O3"
   Unique Field Definition:  OMI-Specific
   Description:       >

     Bits contain several processing quality flags:
       0    - failed convergence check
       1    - solar zenith angle, lat., lon., out of range (SZA > 88 deg) 
       2    - cloud pressure less than low range of table
       3    - cloud pressure greater than surface pressure
       4    - matrix inversion failed
       5    - snow/ice (if second byte of GroundPixelQualityFlags >= 50 and <= 130)
       6    - reflectivity < 0 or > 1.0
       7    - bad radiances detected
       8    - aerosol index flag
       9    - radiance PixelQuality error
       10   - radiance PixelQuality warning
       11   - irradiance PixelQuality error
       12   - irradiance PixelQuality warning
       13   - effective surface pressure retrieved because cloud fraction < 0.05
       14   - missing data
       15   - geolocation error
 
 - Field Name:        RadiativeCloudFraction
   Data Type:         HE5T_NATIVE_FLOAT
   Dimensions:        nXtrack,nTimes
   Valid Range:       0.0 to 1.0
   Missing Value:     -9999.0
   Units:             NoUnits
   Data Source:       PGE
   Title:             "Radiative Cloud Fraction"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The radiative cloud fraction is f_eff*I_cld/I_obs, where f_eff is the effective
    cloud fraction, I_cld is the cloud radiance, I_obs is the OMI-measured radiance.
 
 - Field Name:        Reflectivity
   Data Type:         HE5T_NATIVE_FLOAT
   Dimensions:        nXtrack,nTimes
   Valid Range:       0.0 to 1.0
   Missing Value:     -9999.0
   Units:             NoUnits
   Data Source:       PGE
   Title:             "Reflectivity"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The Lambert Equivalent Reflectivity associated with the ground pixel.
 
 - Field Name:        Residual_bias 
   Data Type:         HE5T_NATIVE_FLOAT 
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is -1.0 to 1.0
   Missing Value:     -9999.0 
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Residual Bias"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The residual bias characterizing the fitting procedure.

 - Field Name:        Residual_stddev 
   Data Type:         HE5T_NATIVE_FLOAT 
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is 0.0 to 1.0
   Missing Value:     -9999.0 
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Residual Standard Deviation"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The residual standard deviation characterizing the fitting procedure.

 - Field Name:        SmallPixelMean 
   Data Type:         HE5T_NATIVE_FLOAT
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is 0.0 to 10.0
   Missing Value:     -9999.0
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Small Pixel Mean" 
   Unique Field Definition:  OMI-Specific
   Description:       >
    Mean value of small pixels comprising a nominal pixel.

 - Field Name:        SmallPixelStddev 
   Data Type:         HE5T_NATIVE_FLOAT 
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is 0.0 to 10.0
   Missing Value:     -9999.0 
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Small Pixel Standard Deviation"
   Unique Field Definition:  OMI-Specific
   Description:       >
    Standard deviation of small pixel radiances in a nominal pixel.

 - Field Name:        SmallPixelnPix 
   Data Type:         HE5T_NATIVE_INT
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is 2 to 5
   Missing Value:     65535 
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Small Pixel Number of Pixels"
   Unique Field Definition:  OMI-Specific
   Description:       >
    Number of small pixels comprising a nominal pixel.

 - Field Name:        SurfaceReflectivity 
   Data Type:         HE5T_NATIVE_INT
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is 0.0 to 1.0
   Missing Value:     -9999.0
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Surface Reflectivity Climatology TOMS"
   Unique Field Definition:  OMI-Specific
   Description:       >
    Surface minimum reflectivity climatology at 360 nm from TOMS measurements.

 - Field Name:        TerrainPressure 
   Data Type:         HE5T_NATIVE_FLOAT
   Dimensions:        nXtrack,nTimes
   Valid Range:       0.0 to 1200.0
   Missing Value:     -9999.0
   Units:             hPa
   Data Source:       PGE
   Title:             "Terrain Pressure"
   Unique Field Definition:  OMI-Specific
   Description:       >
    The Terrain Pressure (in hPa) associated with the ground pixel.

 - Field Name:        WavelengthShift 
   Data Type:         HE5T_NATIVE_FLOAT 
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is -1.0 to 1.0
   Missing Value:     -9999.0 
   Units:             nm
   Data Source:       PGE 
   Title:             "Wavelength Shift"
   Unique Field Definition:  OMI-Specific
   Description:       >
    Wavelength shift in nm.

 - Field Name:        dIdR
   Data Type:         HE5T_NATIVE_FLOAT
   Dimensions:        nXtrack,nTimes
   Valid Range:       Range is -100 to 100
   Missing Value:     -9999.0
   Units:             NoUnits
   Data Source:       PGE 
   Title:             "Radiance (fractional) refl. sens."
   Unique Field Definition:  OMI-Specific
   Description:       >
    Partial derivative of normalized radiance with respect to reflectivity


Core Metadata:

 - Metadata Name:     AssociatedInstrumentShortName
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids:   Valid is "OMI".
   Data Source:       MCF
   Description:       Actual is "OMI".

 - Metadata Name:     AssociatedPlatformShortName
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids:   Valid is "Aura".
   Data Source:       MCF
   Description:       Actual is "Aura".

 - Metadata Name:     AssociatedSensorShortName
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids:   Valids are "CCD Ultra Violet" and "CCD Visible".
   Data Source:       MCF
   Description:       Actual is "CCD Ultra Violet".

 - Metadata Name:     AutomaticQualityFlag
   Mandatory:         T
   Data Type:         VA64
   Number of Values:  1
   Range or Valids:   Valids are "Passed", "Suspect" and "Failed".
   Data Source:       PGE
   Description: >

     A granule-level quality flag that applies generally to the granule
     and specifically to the parameters at the granule level.

 - Metadata Name:     AutomaticQualityFlagExplanation
   Mandatory:         T
   Data Type:         VA255
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       PGE
   Description: >

     The AutomaticQualityFlag is set to
     1) "Passed" if QAPercentHighQualityData >= 80%,
     2) "Suspect" if QAPercentHightQualityData >= 20% or if the input L1B
         file does not have its AutomaticQualityFlag set to "Passed", and
     3) "Failed" if QAPercentHighQualityData < 20%.

 - Metadata Name:     DayNightFlag
   Mandatory:         T
   Data Type:         VA5
   Number of Values:  1
   Range or Valids:   Valids are "Day", "Night" and "Both".
   Data Source:       MCF
   Description:       Actual is "Day".

 - Metadata Name:     EquatorCrossingDate
   Mandatory:         T
   Data Type:         D
   Number of Values:  1
   Range or Valids:   Range is "2003-01-01" to "2099-12-31".
   Data Source:       L1B
   Description: >

     The date of the ascending equator crossing in the granule.

 - Metadata Name:     EquatorCrossingLongitude
   Mandatory:         T
   Data Type:         LF
   Number of Values:  1
   Range or Valids:   Range is -1.80d-02 to 1.80d+02.
   Data Source:       L1B
   Description: >

     The terrestrial longitude of the ascending equator crossing in the
     granule.

 - Metadata Name:     EquatorCrossingTime
   Mandatory:         T
   Data Type:         T
   Number of Values:  1
   Range or Valids:   Range is "01:00:0.000000" to "01:59:59.999999".
   Data Source:       L1B
   Description: >

     The time of the ascending equator crossing in the granule.

 - Metadata Name:     InputPointer
   Mandatory:         T
   Data Type:         VA255
   Number of Values:  0 to 10
   Range or Valids: >

     Valid file names, each in double quotes, separated by commas, all
     surrounded by curved brackets.

   Data Source:       PGE
   Description: >

    Example is
    ("OMI-Aura_L1-OML1BRVG_2002m0630t2354-o21434_v001-2003m0327t181402.he4", 
    "OMI-Aura_L1-OML1BIRR_2002m0630t2307-o21434_v001-2003m0429t155233.he4", 
    "terr_prs.txt", "chl.txt", "oc_raman_omi.dat", "ring_gomi_p_all.dat" )

 - Metadata Name:     LocalGranuleID
   Mandatory:         T
   Data Type:         VA80
   Number of Values:  1
   Range or Valids: >

     "OMI-Aura_L2-OMCLDRR_2003m0101t0000-o00000_v001-2003m0101t000000.he5" to
     "OMI-Aura_L2-OMCLDRR_2099m1231t2359-o99999_v999-2099m1231t235959.he5"

   Data Source:       PGE
   Description: >

     Example is
     "OMI-Aura_L2-OMCLDRR_2002m0630t2354-o21434_v001-2003m0815t210014.he5"
     (see Appendix E of Reference 3).

 - Metadata Name:     OperationalQualityFlag
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids: >

     Valids are "Passed", "Failed", "Being Investigated", "Not Investigated",
     "Inferred Passed", "Inferred Failed" and "Suspect".

   Data Source:       PGE
   Description: >

     A granule-level quality flag that applies generally to the granule and
     specifically to the parameters at the granule level.

 - Metadata Name:     OperationalQualityFlagExplanation
   Mandatory:         T
   Data Type:         VA255
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       PGE
   Description: >

     The criteria for setting the OperationalQualityFlag should be stated
     here (this Metadata will not appear in the granule).

 - Metadata Name:     OperationMode
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids: >

     Valids are "Calibration", "Diagnostic", "Initialization", "Launch",
     "Normal", "Roll", "Routine", "Safe", "Solar Calibration", "Standby",
     "Survival" and "Test".

   Data Source:       PCF
   Description:       Example is "Normal" (not implemented yet).

 - Metadata Name:     OrbitNumber
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   1 to 999999
   Data Source:       L1B
   Description:       The OMI orbit number.

 - Metadata Name:     ParameterName
   Mandatory:         T
   Data Type:         VA40
   Number of Values:  1
   Range or Valids:   Valid is "Cloud_Pressure".
   Data Source:       PGE
   Description: >

     The measured science parameter expressed in the granule.

 - Metadata Name:     PGEVERSION
   Mandatory:         T
   Data Type:         VA10
   Number of Values:  1
   Range or Valids:   Range is "0.0.0" to "9.9.99".
   Data Source:       PCF
   Description:       Actual is "1.9.0" (see Appendix K of Reference 3).

 - Metadata Name:     ProductionDateTime
   Mandatory:         T
   Data Type:         DT
   Number of Values:  1
   Range or Valids: >

     "2003-01-01T00:00:00.000Z" to "2099-12-31T24:59:59.999Z"

   Data Source:       TK
   Description:       The date and time of the Level 2 processing.

 - Metadata Name:     QAPercentCloudCover
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   0 to 100
   Data Source:       PGE
   Description: >

     The percent of the data in the granule that have cloud cover with cloud
     fraction greater than 0.02.

 - Metadata Name:     QAPercentMissingData
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   0 to 100
   Data Source:       PGE
   Description: >

     The percent of the data in the granule that are missing.

 - Metadata Name:     QAPercentOutofBoundsData
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   0 to 100
   Data Source:       PGE
   Description: >

     The percent of the data in the granule that are out of bounds data.

 - Metadata Name:     RangeBeginningDate
   Mandatory:         T
   Data Type:         D
   Number of Values:  1
   Range or Valids:   Range is "2003-01-01" to "2099-12-31".
   Data Source:       L1B
   Description:       The year, month and day when the granule began.

 - Metadata Name:     RangeBeginningTime
   Mandatory:         T
   Data Type:         T
   Number of Values:  1
   Range or Valids:   Range is "00:00:00.000000" to "23:59:59.999999".
   Data Source:       L1B
   Description: >

     The hour, minute, second and fraction of a second when the granule
     began.

 - Metadata Name:     RangeEndingDate
   Mandatory:         T
   Data Type:         D
   Number of Values:  1
   Range or Valids:   2003-01-01" to "2099-12-31"
   Data Source:       L1B
   Description:       The year, month and day when the granule ended.

 - Metadata Name:     RangeEndingTime
   Mandatory:         T
   Data Type:         T
   Number of Values:  1
   Range or Valids:   Range is "00:00:00.000000" to "23:59:59.999999".
   Data Source:       L1B
   Description: >

     The hour, minute, second and fraction of a second when the granule
     ended.

 - Metadata Name:     ReprocessingActual
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids: >

     Valids are "processed 1 time", "processed 2 times", etc...

   Data Source:       PCF
   Description: >

     An indication of what reprocessing has been performed on the granule.

 - Metadata Name:     LocalVersionID
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids: >

     Valids are "RFC1321 MD5 = not yet calculated" and  "RFC1321 MD5 = [0-9,a-f]{32}"

   Data Source:       PCF
   Description: >

     MD5 fingerprint of the HDF product file.

 - Metadata Name:     ReprocessingPlanned
   Mandatory:         T
   Data Type:         VA45
   Number of Values:  1
   Range or Valids: >

     Valids are "no further update anticipated", "further update anticipated"
     and "further update anticipated using enhanced PGE".

   Data Source:       DP
   Description:       Actual is "further update is anticipated".

 - Metadata Name:     ScienceQualityFlag
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids: >

     Valids are "Passed", "Failed", "Being Investigated", "Not Investigated",
     "Inferred Passed", "Inferred Failed" and "Suspect".

   Data Source:       DP
   Description:       Actual is "Not Investigated".

 - Metadata Name:     ScienceQualityFlagExplanation
   Mandatory:         T
   Data Type:         VA255
   Number of Values:  1
   Range or Valids:   Not applicable (free format).
   Data Source:       DP
   Description: >

     An explanation of the criteria used to set the science quality flag
     should go here.

 - Metadata Name:     ShortName
   Mandatory:         T
   Data Type:         VA8
   Number of Values:  1
   Range or Valids:   Valid is "OMCLDRR".
   Data Source:       MCF
   Description:       Actual is "OMCLDRR".

 - Metadata Name:     SizeMBECSDataGranule
   Mandatory:         F
   Data Type:         LF
   Number of Values:  1
   Range or Valids:   0.00d+00 to 1.00d+04
   Data Source:       DSS
   Description: >

     The volume of data contained in the granule in Mb (this Metadata will
     not appear in the granule).

 - Metadata Name:     VersionID
   Mandatory:         T
   Data Type:         SI
   Number of Values:  1
   Range or Valids:   000 to 999
   Data Source:       MCF
   Description:       Actual is 001 for test and pre-launch.

Product Specific Attributes:

 - Metadata Name:     EndBlockNr
   Mandatory:         T
   Data Type:         SI
   Number of Values:  1 to 500
   Range or Valids:   1 to 50
   Data Source:       L1B
   Description:       The number of the NOSE end block along the track.

 - Metadata Name:     ExpeditedData
   Mandatory:         T
   Data Type:         VA10
   Number of Values:  1
   Range or Valids:   Valids are "TRUE" and "FALSE".
   Data Source:       L1B
   Description:       The indicator for expedited L0 data.

 - Metadata Name:     ExposureTimes
   Mandatory:         T
   Data Type:         F
   Number of Values:  1 to 256
   Range or Valids:   0.0 to 2000.0
   Data Source:       L1B
   Description: >

     An array containing the exposure times in seconds used for the
     measurements.

 - Metadata Name:     InstrumentConfigurationIDs
   Mandatory:         T
   Data Type:         SI
   Number of Values:  1 to 256
   Range or Valids:   0 to 255
   Data Source:       L1B
   Description: >

     An array containing the instrument configuration identifiers used
     for the measurements.

 - Metadata Name:     MasterClockPeriods
   Mandatory:         T
   Data Type:         F
   Number of Values:  1 to 128
   Range or Valids:   0.0 to 10.0
   Data Source:       L1B
   Description: >

     An array containing the master clock periods in seconds used for
     the measurements.

 - Metadata Name:     NrMeasurements
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   0 to 9999
   Data Source:       L1B
   Description: >

     The number of measurements in the granule (per output product).

 - Metadata Name:     NrSpatialZoom
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   0 to 9999
   Data Source:       L1B
   Description:       The number of measurements in spatial zoom mode.

 - Metadata Name:     NrSpectralZoom
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   0 to 9999
   Data Source:       L1B
   Description:       The number of measurements in spectral zoom mode.

 - Metadata Name:     NrZoom
   Mandatory:         T
   Data Type:         I
   Number of Values:  1
   Range or Valids:   0 to 9999
   Data Source:       L1B
   Description:       The number of measurements in zoom modes.

 - Metadata Name:     PathNr
   Mandatory:         T
   Data Type:         I
   Number of Values:  1 to 500
   Range or Valids:   1 to 466
   Data Source:       L1B
   Description:       Number of the NOSE path within the repeat cycle.

 - Metadata Name:     SolarEclipse
   Mandatory:         T
   Data Type:         VA10
   Number of Values:  1
   Range or Valids:   Valids are "TRUE" and "FALSE".
   Data Source:       L1B
   Description: >

     The indicator that during part of the measurements a solar eclipse
     occurred.

 - Metadata Name:     SouthAtlanticAnomalyCrossing
   Mandatory:         T
   Data Type:         VA10
   Number of Values:  1
   Range or Valids:   Valids are "TRUE" and "FALSE".
   Data Source:       L1B
   Description: >

     The indicator that during part of the measurements the spacecraft
     was in the SAA.

 - Metadata Name:     SpacecraftManeuverFlag
   Mandatory:         T
   Data Type:         VA10
   Number of Values:  1
   Range or Valids:   Valids are "TRUE", "FALSE" and "UNKNOWN".
   Data Source:       L1B
   Description: >

     The indicator that during part of the measurements the spacecraft
     was performing a maneuver.

 - Metadata Name:     StartBlockNr
   Mandatory:         T
   Data Type:         SI
   Number of Values:  1 to 500
   Range or Valids:   1 to 50
   Data Source:       L1B
   Description:       Number of the NOSE start block along the track.

Archived Metadata:

 - Metadata Name:     ESDTDescriptorRevision
   Mandatory:         T
   Data Type:         VA20
   Number of Values:  1
   Range or Valids:   Range is "0.0.0" to "9.9.99".
   Data Source:       MCF
   Description: >

     The version of the ESDT descriptor file as determined by ECS.

 - Metadata Name:     LongName
   Mandatory:         T
   Data Type:         VA80
   Number of Values:  1
   Range or Valids: >

     Valid is
     "OMI/Aura Cloud Pressure and Fraction (Raman Scattering) 1-Orbit L2 Swath 13x24km".

   Data Source:       MCF
   Description: >

     Actual is
     "OMI/Aura Cloud Pressure and Fraction (Raman Scattering) 1-Orbit L2 Swath 13x24km"
     (see Section 7.0 of Reference 2).

References: >

  1. "OMI Algorithm Theoretical Basis Document, Volume III, Clouds, Aerosols, and 
     Surface UV Irradiance"
     (OMI-ATBD-VOL3, ATBD-OMI-03, Version 2.0, August 2002)

  2. "HDF-EOS Aura File Format Guidelines"
     (OMI-AURA-DATA-GUIDE, Version 1.3, 5 August 2003)

  3. "OMI Science Software Delivery Guide for Version 0.9"
     (OMI-SSDG-0.9.8, Version 0.9.8, 14 August 2003)

  4. "OMI GDPS Input/Output Data Specification (IODS) Volume 2"
     (OMI-GDPS-IODS-2, SD-OMIE-7200-DS-467, 9 April 2003)

  5. "Release 6A Implementation Earth Science Data Model for the ECS Project"
     (420-TP-022-002, June 2001)
     (http://edhs1.gsfc.nasa.gov/waisdata/rel6/html/tp4202202.html and
      http://edhs1.gsfc.nasa.gov/waisdata/rel6/html/tp42022_adds.html)

  6. "OMI L2 - L4 Metadata Reference Guide"
     (3 July 2002)
     (https://omiwww.gsfc.nasa.gov/mlinda/OMImetadataRefGuide.html)
