!! F90
 ! 
 ! Description:
 ! MODULE L2_class, contains data structures for storing L2 data
 !
 !!Revision History:
 ! Revision 0.1  12/26/2002  Kai Yang/UMBC
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
!!
MODULE L2_data_structure

   USE ISO_C_BINDING, ONLY: C_LONG

   IMPLICIT NONE
   INTEGER (KIND=4), PARAMETER :: NDIM_MAX = 10
   INTEGER (KIND=4), PARAMETER :: NFLDS_MAX = 40
   INTEGER (KIND=4), PARAMETER :: MAX_STR_LEN = 512

   TYPE, PUBLIC :: L2_generic_type
      CHARACTER ( LEN = MAX_STR_LEN ) :: filename, swathname
      CHARACTER ( LEN = MAX_STR_LEN ) :: dimnames
      CHARACTER( LEN = 80 ), DIMENSION(NFLDS_MAX) :: fieldname
      !INTEGER (KIND = 4) :: sw_fid, swathID
      !INTEGER (KIND = 4) :: nDims, nFields
      !INTEGER (KIND = 4) :: iLine, eLine, nLine, nTotLine
      !INTEGER (KIND = 4), DIMENSION(NDIM_MAX) :: dimSizes
      !INTEGER (KIND = 4), DIMENSION(NFLDS_MAX) :: rank
      !INTEGER (KIND = 4), DIMENSION(NFLDS_MAX,3) :: dims
      !INTEGER (KIND = 4), DIMENSION(NFLDS_MAX) :: elmSize,  pixSize
      !INTEGER (KIND = 4), DIMENSION(NFLDS_MAX) :: lineSize, blkSize
      !INTEGER (KIND = 4) :: SumElmSize, SumLineSize
      !INTEGER (KIND = 4), DIMENSION(0:NFLDS_MAX) :: accuLineSize, accuBlkSize
      INTEGER (KIND = C_LONG) :: sw_fid, swathID
      INTEGER (KIND = C_LONG) :: nDims, nFields
      INTEGER (KIND = C_LONG) :: iLine, eLine, nLine, nTotLine
      INTEGER (KIND = C_LONG), DIMENSION(NDIM_MAX) :: dimSizes
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX) :: rank
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX,3) :: dims
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX) :: elmSize,  pixSize
      INTEGER (KIND = C_LONG), DIMENSION(NFLDS_MAX) :: lineSize, blkSize
      INTEGER (KIND = C_LONG) :: SumElmSize, SumLineSize
      INTEGER (KIND = C_LONG), DIMENSION(0:NFLDS_MAX) :: accuLineSize, accuBlkSize
      LOGICAL :: initialized
      ! geo and data fields
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: data
   END TYPE L2_generic_type
END MODULE L2_data_structure
