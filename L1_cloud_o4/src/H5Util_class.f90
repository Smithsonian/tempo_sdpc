MODULE H5Util_class
!
! DESCRIPTION:
!    a set of general purpose high-level HDF5 utilities to handle dataset 
!    creation, writing and reading. It borrows many lines source code and some 
!    ideas from Kai Yang's UTIL_lh5_class.f90 of TO3_CORE. It takes care of all
!    my HDF5 I/O needs for working with OMPS SDR and EDR type of HDF5 files, 
!    and beyond.
!..............................................................................
! USAGE:
!
!---
!--- SCENARIO #1. create and write data: 
!---
!
!    1) open file:
!       call h5fcreate_f("test.h5", H5F_ACC_TRUNC_F, fileId, hdferr)
!
!    2) create a group (optional)
!       call h5gcreate_f(fileId, "/GROUP_1", groupId, hdferr) 
!
!    3) define data type: H5SDS_T
!       long_name, units, valid_range and fill_value are optional attributes: 
!          sds%name = "real_sp"
!          sds%datatype_id = H5T_NATIVE_REAL
!          sds%rank = 2
!          sds%dims(1:sds%rank) = (/ 10, 20/)
!          sds%long_name = "a real(SP) dataset"
!          sds%units = "unitless"
!          sds%valid_range = (/ 0.0, 103.0 /)
!          sds%fill_value =  FILLVALUE_SP
!
!    4) create dataset: parentId is either fileId or groupId; L_createAttribute
!       is set to false [default] if it is plain vanilla SDS without attributes.
!       errorStatus = H5Util_createDataset(parentId, sds, L_createAttribute)
!
!    5) write data_array[] to dataset:
!       errorStatus = H5Util_writeDataset(sds, data_array, start, count)
!
!    6) close dataset:
!       errorstatus = h5Util_disposeDataset(sds)
!
!    7) close group (optional)
!       call h5gclose_f(groupId, hdferr)
!
!    8) close file
!       call h5fclose_f(fileId, hdferr)
!
!
!---
!--- SCENARIO #2. read data:
!---
!
!    1) open file:
!       call h5fcreate_f("test.h5", H5F_ACC_RDONLY_F, fileId, hdferr)
!
!    2) acquire a group_id (if parent is a group): 
!       call h5gopen_f(fileId, "GROUP_1", groupId, hdferr)
!
!    3) define data type: H5SDS_T.  All you need is to set 
!       the sds%name here. 
!       type(H5SDS_T) :: sds = &
!            H5SDS_T("real_sp",-1,-1,-1,(/0,0,0,0,0,0,0/),"","",(/0,0/),0)
!
!    4) select dataset: parentId is either fileId or groupId
!       errorStatus = H5Util_selectDataset(parentId, sds)
!
!    5) read data_array[]:
!       errorStatus = H5Util_readDataset(sds, data_array, start, count)
!
!    6) close dataset:
!       errorstatus = h5Util_disposeDataset(sds)
!
!    7) close group (optional)
!       call h5gclose_f(groupId, hdferr)
!
!    8) close file
!       call h5fclose_f(fileId, hdferr)
!...............................................................................
! REVISION HISTORY:
!    initial version:  Wed Mar 30 16:11:26 EDT 2011
!    complete history: svn log file:///omps/cm/app/JLi_Fortran_Toolkit
!
! AUTHOR: 
!    Jason Li (SSAI)
!*******************************************************************************

use HDF5
use DataTypeDef, only: SP, DP, I1B, I2B, I4B 
use Constants_class, only: SUCCESS_STATE, FAILURE_STATE, WARNING_STATE
use ErrorHandler_class, only: Display_Message

IMPLICIT NONE
PRIVATE

!...............................................................................
!
! constants:
!
!...............................................................................

INTEGER (I4B), PARAMETER :: MAXRANK = 7 ! max rank supported by Fortran 90 

INTEGER (I4B), PARAMETER :: zero = 0

!...............................................................................
!
! public methods and data types:
!
!...............................................................................

! valid_range and fill_value in memory are double precision float (REAL(DP)).
! In the file,  they are stored, rightly so, as the same type as the dataset, 
! in other words it is sds%datatype_id.

! initialization:
! type(H5SDS_T) :: sds = H5SDS_T("",-1,-1,-1,(/0,0,0,0,0,0,0/),"","",(/0,0/),0)

type, public :: H5SDS_T 
      character(len=128) :: name 
      integer(HID_T) :: dataset_id
      integer(HID_T) :: datatype_id
      integer :: rank 
      integer(HSIZE_T), dimension(MAXRANK) :: dims
      character(len=128) :: long_name
      character(len=128) :: units
      character(len=2048) :: description
!     real(DP), dimension(2) :: valid_range
!      real(DP) :: fill_value
end type H5SDS_T

type, public :: H5SDS_DIM_T 
      character(len=128) :: name 
      integer(HID_T) :: dataset_id
      integer(HID_T) :: datatype_id
      integer :: rank 
      integer(HSIZE_T), dimension(MAXRANK) :: dims
      character(len=128) :: long_name
      character(len=128) :: units
      character(len=2048) :: description
!     real(DP), dimension(2) :: valid_range
      real(DP) :: fill_value
end type H5SDS_DIM_T


! initialization:
! type(H5SDS_CHAR_T) :: sds = H5SDS_CHAR_T( "", -1, -1, -1, 0, (/ 0 /) )
type, public :: H5SDS_CHAR_T 
      character(len=128) :: name 
      integer(HID_T) :: dataset_id
      integer(HID_T) :: datatype_id
      integer :: rank 
      integer :: character_length
!     support rank=1 string array:
!     integer(HSIZE_T), dimension(MAXRANK) :: dims 
      integer(HSIZE_T), dimension(1) :: dims 
end type H5SDS_CHAR_T

public :: H5Util_createDataset
public :: H5Util_selectDataset
public :: H5Util_readDataset
public :: H5Util_writeDataset
public :: H5Util_disposeDataset

!...............................................................................
! 
! function overloading:
!
!...............................................................................

interface H5Util_createDataset
   module procedure createDataset_char
   module procedure createDataset
end interface

interface H5Util_selectDataset
   module procedure selectDataset_char
   module procedure selectDataset
end interface

interface H5Util_disposeDataset
   module procedure disposeDataset_char
   module procedure disposeDataset
end interface

interface H5Util_writeDataset
   MODULE PROCEDURE writeDataset_1DIK4
   MODULE PROCEDURE writeDataset_2DIK4
   MODULE PROCEDURE writeDataset_3DIK4
   MODULE PROCEDURE writeDataset_4DIK4
   MODULE PROCEDURE writeDataset_5DIK4
   MODULE PROCEDURE writeDataset_6DIK4
   MODULE PROCEDURE writeDataset_7DIK4

   MODULE PROCEDURE writeDataset_1DRK4
   MODULE PROCEDURE writeDataset_2DRK4
   MODULE PROCEDURE writeDataset_3DRK4
   MODULE PROCEDURE writeDataset_4DRK4
   MODULE PROCEDURE writeDataset_5DRK4
   MODULE PROCEDURE writeDataset_6DRK4
   MODULE PROCEDURE writeDataset_7DRK4

   MODULE PROCEDURE writeDataset_1DRK8
   MODULE PROCEDURE writeDataset_2DRK8
   MODULE PROCEDURE writeDataset_3DRK8
   MODULE PROCEDURE writeDataset_4DRK8
   MODULE PROCEDURE writeDataset_5DRK8
   MODULE PROCEDURE writeDataset_6DRK8
   MODULE PROCEDURE writeDataset_7DRK8

  MODULE PROCEDURE writeDataset_1DCHAR
!!$   MODULE PROCEDURE writeDataset_2DCHAR
!!$   MODULE PROCEDURE writeDataset_3DCHAR
!!$   MODULE PROCEDURE writeDataset_4DCHAR
!!$   MODULE PROCEDURE writeDataset_6DCHAR
!!$   MODULE PROCEDURE writeDataset_6DCHAR
!!$   MODULE PROCEDURE writeDataset_7DCHAR

!  MODULE PROCEDURE writeDataset_1DString
end interface

interface H5Util_readDataset
   MODULE PROCEDURE readDataset_1DIK4
   MODULE PROCEDURE readDataset_2DIK4
   MODULE PROCEDURE readDataset_3DIK4
   MODULE PROCEDURE readDataset_4DIK4
   MODULE PROCEDURE readDataset_5DIK4
   MODULE PROCEDURE readDataset_6DIK4
   MODULE PROCEDURE readDataset_7DIK4

   MODULE PROCEDURE readDataset_1DRK4
   MODULE PROCEDURE readDataset_2DRK4
   MODULE PROCEDURE readDataset_3DRK4
   MODULE PROCEDURE readDataset_4DRK4
   MODULE PROCEDURE readDataset_5DRK4
   MODULE PROCEDURE readDataset_6DRK4
   MODULE PROCEDURE readDataset_7DRK4

   MODULE PROCEDURE readDataset_1DRK8
   MODULE PROCEDURE readDataset_2DRK8
   MODULE PROCEDURE readDataset_3DRK8
   MODULE PROCEDURE readDataset_4DRK8
   MODULE PROCEDURE readDataset_5DRK8
   MODULE PROCEDURE readDataset_6DRK8
   MODULE PROCEDURE readDataset_7DRK8

   MODULE PROCEDURE readDataset_1DCHAR
end interface

CONTAINS

!
!///////////////////////// P U B L I C   M E T H O D S /////////////////////////
!

function createDataset(parent_id, sds, L_createAttributes) result(errorStatus)

use H5LT
use H5DS
implicit none

character(len=*), parameter :: routineName = "H5Util_createDataset"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_T), intent(inout) :: sds
integer :: errorStatus
logical, intent(in), optional :: L_createAttributes

integer(HID_T) :: dataspace_id, datatype_id, property_id
integer :: hdferr, exists_status

!integer(HID_T) :: dimscale_nc03, dimscale_nc08, dimscale_nc10, dimscale_nc12
!integer(HID_T) :: dimscale_nc13, dimscale_nl11, dimscale_nl13, dimscale_nl15
!integer(HID_T) :: dimscale_nb05, dimscale_nb10, dimscale_nb12, dimscale_nb15
!integer(HID_T) :: dimscale_nl20, dimscale_nl20b, dimscale_nl21
!integer(HID_T) :: dimscale_nl01, dimscale_nw01, dimscale_nc01
integer(HID_T) :: dimscale_nx01, dimscale_nt01

character(len=256) :: msg

errorStatus = SUCCESS_STATE

exists_status = h5ltfind_dataset_f(parent_id, sds%name)
if(exists_status == 1) then
   errorStatus = WARNING_STATE
   write(msg,*) trim(sds%name), " already exists, cannot create it again" 
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! verify I have enough information to create a dataset:
!
if(sds%rank <= 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide rank for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

if(sds%datatype_id < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide datatype_id for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! define dataspace:
!

call H5Screate_simple_f(sds%rank, sds%dims, dataspace_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "H5Screate_simple_f() failed for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! create a copy of data type:
!
call h5tcopy_f(sds%datatype_id, datatype_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "h5tcopy_f() failed for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! create property
!
call H5Pcreate_f(H5P_DATASET_CREATE_F, property_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "H5Pcreate_f() failed for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

call H5Pset_chunk_f(property_id, sds%rank, sds%dims, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "H5Pset_chunk_f() failed for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

call H5Pset_deflate_f(property_id, 5, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "H5Pset_deflate_f() failed for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! create dataset:
!

call H5Dcreate_f(parent_id, sds%name, &
                 datatype_id, dataspace_id, sds%dataset_id, &
                 hdferr, property_id)

! ----------------------------
! Dimensions for OMCDO2N begin
! ----------------------------
if(sds%name .eq. "nTimes") then
  dimscale_nt01=sds%dataset_id
  call H5DSset_scale_f(dimscale_nt01,hdferr,"nTimes")
end if

if(sds%name .eq. "nXtrack") then
  dimscale_nx01=sds%dataset_id
  call H5DSset_scale_f(dimscale_nx01,hdferr,"nXtrack")
end if

! -------------------
! OMIAura DATA fields
! -------------------

if(sds%name .eq. "CloudPressure") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudPressureNotClipped") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudPressureUncertainty") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "EffectiveCloudFraction") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "EffectiveCloudFractionNotClipped") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "EffectiveCloudFractionUncertainty") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudRadianceFraction440") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudRadiativeFractionNotClipped440") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudRadiativeFractionUncertainty440") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudRadianceFraction466") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudRadiativeFractionNotClipped466") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "CloudRadiativeFractionUncertainty466") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "SurfaceLER466") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "SurfaceLER440") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "SceneLER466") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "SceneLER440") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "ScenePressure") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "MeasurementQualityFlags") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
end if

if(sds%name .eq. "ProcessingQualityFlags") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "GLER440") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "GLER466") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "TerrainPressure") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "TerrainPressureStdDev") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "TerrainHeight") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "TerrainHeightStdDev") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "LandAreaFraction") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "SlantColumnAmountO2O2") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "Reflectance") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

! -------------------------------
! OMIAura GEOLOCATION_DATA fields
! -------------------------------

if(sds%name .eq. "Time") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
end if

if(sds%name .eq. "SecondsInDay") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
end if

if(sds%name .eq. "Latitude") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "Longitude") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "SolarZenithAngle") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "ViewingZenithAngle") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "RelativeAzimuthAngle") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "GroundPixelQualityFlags") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

if(sds%name .eq. "XtrackQualityFlags") then
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nt01, 1, hdferr)
  call H5DSattach_scale_f(sds%dataset_id, dimscale_nx01, 2, hdferr)
end if

!
! write dataset attributes (optional):
!

if(present(L_createAttributes)) then
  if(L_createAttributes) call writeDatasetAttributes(sds, hdferr)
endif

!
! closing up shop:
!
call H5Tclose_f(datatype_id, hdferr)
call H5Pclose_f(property_id, hdferr)
call H5Sclose_f(dataspace_id, hdferr)

99999 return

end function createDataset

!...............................................................................

function createDataset_char(parent_id, sds) result (errorStatus)

use H5LT
implicit none

character(len=*), parameter :: routineName = "createDataset_char"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_CHAR_T), intent(inout) :: sds
integer :: errorStatus

integer(HID_T) :: dataspace_id, datatype_id
character(len=256) :: msg
integer :: hdferr, exists_status
integer(SIZE_T) :: char_len

errorStatus = SUCCESS_STATE

exists_status = h5ltfind_dataset_f(parent_id, sds%name)
if(exists_status == 1) then
   errorStatus = WARNING_STATE
   write(msg,*) trim(sds%name), " already exists, cannot create it again"
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

!
! verify I have enough information to create a dataset:
!
if(sds%rank <= 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide rank for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

if(sds%datatype_id /= H5T_NATIVE_CHARACTER) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide datatype_id for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

if(sds%character_length <= 0) then
   errorStatus = FAILURE_STATE
   write(msg,*) "must provide character_length for dataset: ", trim(sds%name)
   call Display_Message(routineName, trim(msg), errorStatus)
   goto 99999
endif

call h5tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
char_len=sds%character_length
call h5tset_size_f(datatype_id, char_len, hdferr)
call H5screate_simple_f(sds%rank, sds%dims, dataspace_id, hdferr)
call h5dcreate_f(parent_id, sds%name, &
                 datatype_id, dataspace_id, sds%dataset_id, hdferr)
call h5sclose_f(dataspace_id, hdferr)
call h5tclose_f(datatype_id, hdferr)

99999 return

end function createDataset_char

!*******************************************************************************

function selectDataset(parent_id, sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "selectDataset"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_T), intent(inout) :: sds
integer :: errorStatus

integer(HID_T) :: dataspace_id
integer(HSIZE_T), dimension(MAXRANK) :: maxdims
integer :: rank, hdferr

errorStatus = SUCCESS_STATE

if(sds%dataset_id > 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "dataset " // trim(sds%name) // " has already been selected", &
                errorStatus)
   goto 99999
endif

call h5dopen_f(parent_id, sds%name, sds%dataset_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "HDF5 error acquiring " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_type_f(sds%dataset_id, sds%datatype_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_type_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_space_f(sds%dataset_id, dataspace_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_space_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5sget_simple_extent_ndims_f(dataspace_id, sds%rank, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_ndims_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5sget_simple_extent_dims_f(dataspace_id, sds%dims, maxdims, rank)
if(rank < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_dims_f error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

! read in attributes if they exist:

call readDatasetAttributes(sds, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "readDatasetAttributes error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

99999 return

end function selectDataset

!...............................................................................

function selectDataset_char(parent_id, sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "selectDataset_char"

integer(HID_T), intent(in) :: parent_id
type(H5SDS_CHAR_T), intent(inout) :: sds
integer :: errorStatus

integer(HID_T) :: dataspace_id
integer(HSIZE_T), dimension(MAXRANK) :: maxdims
integer :: rank, hdferr

errorStatus = SUCCESS_STATE

if(sds%dataset_id > 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "dataset " // trim(sds%name) // " has already been selected", &
                errorStatus)
   goto 99999
endif

call h5dopen_f(parent_id, sds%name, sds%dataset_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "HDF5 error acquiring " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_type_f(sds%dataset_id, sds%datatype_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_type_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5dget_space_f(sds%dataset_id, dataspace_id, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5dget_space_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5sget_simple_extent_ndims_f(dataspace_id, sds%rank, hdferr)
if(hdferr < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_ndims_f() error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

call h5sget_simple_extent_dims_f(dataspace_id, sds%dims, maxdims, rank)
if(rank < 0) then
   errorStatus = FAILURE_STATE
   call Display_Message(routineName, &
               "h5sget_simple_extent_dims_f error on " // trim(sds%name), &
                errorStatus)
   goto 99999
endif

99999 return

end function selectDataset_char

!*******************************************************************************

function disposeDataset(sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "H5Util_disposeDataset"

type(H5SDS_T), intent(inout) :: sds
integer :: errorStatus, hdferr

errorStatus = SUCCESS_STATE

if(sds%dataset_id <= 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "invalid dataset_id for " // trim(sds%name), errorStatus)
   return
endif

call h5dclose_f(sds%dataset_id, hdferr)
sds%dataset_id = -1
sds%datatype_id = -1
sds%rank = 0
sds%dims = 0
sds%long_name = ""
sds%units = ""
sds%description = ""
! sds%valid_range = 0.0_DP
!sds%valid_min = 0.0_DP
!sds%valid_max = 0.0_DP
!sds%fill_value = 0.0_DP

end function disposeDataset

!...............................................................................

function disposeDataset_char(sds) result(errorStatus)

implicit none

character(len=*), parameter :: routineName = "disposeDataset_char"

type(H5SDS_CHAR_T), intent(inout) :: sds
integer :: errorStatus, hdferr


errorStatus = SUCCESS_STATE

if(sds%dataset_id <= 0) then
   errorStatus = WARNING_STATE
   call Display_Message(routineName, &
               "invalid dataset_id for " // trim(sds%name), errorStatus)
   return
endif

call h5dclose_f(sds%dataset_id, hdferr)

sds%dataset_id = -1
sds%datatype_id = -1
sds%rank = 0
sds%dims = 0
sds%character_length = 0

end function disposeDataset_char

!*******************************************************************************
!
! errorStatus = writeDataset_*(ds, data_in, start, count)
!
!*******************************************************************************

FUNCTION writeDataset_1DIK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 1
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ",ds%name
   call Display_Message("writeDataset_1DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_1DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_1DIK4

!...............................................................................

FUNCTION writeDataset_2DIK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 2
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ",ds%name
   call Display_Message("writeDataset_2DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_2DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_2DIK4

!...............................................................................

FUNCTION writeDataset_3DIK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 3
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_3DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_3DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_3DIK4

!...............................................................................

FUNCTION writeDataset_4DIK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 4
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_4DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_4DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_4DIK4

!...............................................................................

FUNCTION writeDataset_5DIK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 5
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_5DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_5DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_5DIK4

!...............................................................................

FUNCTION writeDataset_6DIK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 6
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_6DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_6DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_6DIK4

!...............................................................................

FUNCTION writeDataset_7DIK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 7
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_7DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_7DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_7DIK4

!...............................................................................

FUNCTION writeDataset_1DRK4( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(SP), DIMENSION(:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 1
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_REAL) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_1DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_1DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank) 
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_1DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_1DRK4

!...............................................................................

FUNCTION writeDataset_2DRK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 2
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_REAL) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_2DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ",ds%name
   call Display_Message("writeDataset_2DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_2DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_2DRK4

!...............................................................................

FUNCTION writeDataset_3DRK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 3
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_REAL) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_3DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_3DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_3DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_3DRK4

!...............................................................................

FUNCTION writeDataset_4DRK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 4
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_REAL) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_4DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_4DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_4DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_4DRK4

!...............................................................................

FUNCTION writeDataset_5DRK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 5
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_REAL) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_5DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_5DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_5DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_5DRK4

!...............................................................................

FUNCTION writeDataset_6DRK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 6
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_REAL) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_6DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_6DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_6DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_6DRK4

!...............................................................................

FUNCTION writeDataset_7DRK4( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 7
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_REAL) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("writeDataset_7DRK4", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_7DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_7DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_7DRK4

!...............................................................................

FUNCTION writeDataset_1DRK8( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(DP), DIMENSION(:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 1
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_1DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_1DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank) 
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_1DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_1DRK8

!...............................................................................

FUNCTION writeDataset_2DRK8( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 2
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_2DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ",ds%name
   call Display_Message("writeDataset_2DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_2DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_2DRK8

!...............................................................................

FUNCTION writeDataset_3DRK8( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 3
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_3DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_3DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_3DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_3DRK8

!...............................................................................

FUNCTION writeDataset_4DRK8( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 4
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_4DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_4DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_4DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_4DRK8

!...............................................................................

FUNCTION writeDataset_5DRK8( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 5
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_5DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_5DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_5DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_5DRK8

!...............................................................................

FUNCTION writeDataset_6DRK8( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 6
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_6DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_6DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_6DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_6DRK8

!...............................................................................

FUNCTION writeDataset_7DRK8( ds, data_in, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 7
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_DOUBLE) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("writeDataset_7DRK8", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_7DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5dwrite_f( ds%dataset_id, ds%datatype_id, data_in, dims, error, &
                 memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_7DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_7DRK8

!...............................................................................

!FUNCTION writeDataset_1DString(ds, data_in) RESULT (status)
!
!implicit none
!
!TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
!CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: data_in
!INTEGER          :: error
!INTEGER (HID_T)  :: datatype_id, memtype_id
!INTEGER (I4B) :: status
!CHARACTER (LEN=256) :: msg
!INTEGER(HSIZE_T), DIMENSION(1) :: dims
!INTEGER :: hdferr
!
!! sanity check:
!
!IF( ds%rank /= 1) then
!   status = FAILURE_STATE
!   write( msg,* ) "only support rank=1 string array for now: ", ds%name
!   call Display_Message("writeDataset_1DString", msg, status)
!   RETURN
!ENDIF
!
!IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!   status = FAILURE_STATE
!   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!   call Display_Message("writeDataset_1DString", msg, status)
!   RETURN
!ENDIF
!
!IF( ds%dataset_id < zero ) THEN
!   status = FAILURE_STATE
!   write( msg,* ) "dataset must be selected before write: ", ds%name
!   call Display_Message("writeDataset_1DString", msg, status)
!   RETURN
!ENDIF
!
! write data:

!dims = size(data_in)
!
!call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
!call h5tset_size_f(memtype_id, len(data_in), hdferr)
!call h5dwrite_f(ds%dataset_id, memtype_id, data_in, dims, error)
!IF( error < zero ) THEN
!   status = FAILURE_STATE
!   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
!   call Display_Message("writeDataset_1DString", msg, status)
!   RETURN
!ENDIF
!
!END FUNCTION writeDataset_1DString

!...............................................................................

FUNCTION writeDataset_1DCHAR( ds, data_in, start, count ) RESULT (status)

! optional keywords:
! start = starting indicies, zero based
! count = number of data points

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: data_in
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
INTEGER          :: error
INTEGER (HID_T)  :: dataspace, memtype_id !, datatype_id
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 1
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
INTEGER :: hdferr
INTEGER(HSIZE_T) :: len_data_in

status = SUCCESS_STATE

IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
   call Display_Message("writeDataset_1DCHAR", msg, status)
   RETURN
ENDIF

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "dataset must be selected before write: ", ds%name
   call Display_Message("writeDataset_1DCHAR", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank) 
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_in )
   call h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
len_data_in=len(data_in)
call h5tset_size_f(memtype_id, len_data_in, hdferr)
call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dims, error, &
                 memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
   call Display_Message("writeDataset_1DCHAR", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION writeDataset_1DCHAR

!!$!...............................................................................
!!$
!!$FUNCTION writeDataset_2DCHAR( ds, data_in, start, count ) RESULT (status)
!!$
!!$! optional keywords:
!!$! start = starting indicies, zero based
!!$! count = number of data points
!!$
!!$TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
!!$CHARACTER(LEN=*), DIMENSION(:,:), INTENT(IN) :: data_in
!!$INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
!!$INTEGER          :: error
!!$INTEGER (HID_T)  :: dataspace, datatype_id, memtype_id 
!!$INTEGER (I4B) :: status
!!$CHARACTER (LEN=256) :: msg
!!$INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
!!$INTEGER :: rankm = 1
!!$INTEGER(HID_T) :: memspace 
!!$INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
!!$INTEGER :: hdferr
!!$
!!$IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!!$   call Display_Message("writeDataset_2DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$IF( ds%dataset_id < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "dataset must be selected before write: ", ds%name
!!$   call Display_Message("writeDataset_2DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$dims(1:ds%rank) = ds%dims(1:ds%rank) 
!!$IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
!!$   call h5dget_space_f(ds%dataset_id, dataspace, error)
!!$   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
!!$   dimsm = SHAPE( data_in )
!!$   call h5screate_simple_f( rankm, dimsm, memspace, error)
!!$ELSE
!!$   dataspace = H5S_ALL_F
!!$   memspace = H5S_ALL_F
!!$ENDIF
!!$
!!$call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
!!$call h5tset_size_f(memtype_id, len(data_in), hdferr)
!!$call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dims, error, &
!!$                 memspace, dataspace )
!!$IF( error < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
!!$   call Display_Message("writeDataset_2DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)
!!$
!!$END FUNCTION writeDataset_2DCHAR
!!$
!!$!...............................................................................
!!$
!!$FUNCTION writeDataset_3DCHAR( ds, data_in, start, count ) RESULT (status)
!!$
!!$! optional keywords:
!!$! start = starting indicies, zero based
!!$! count = number of data points
!!$
!!$TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
!!$CHARACTER(LEN=*), DIMENSION(:,:,:), INTENT(IN) :: data_in
!!$INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
!!$INTEGER          :: error
!!$INTEGER (HID_T)  :: dataspace, datatype_id, memtype_id 
!!$INTEGER (I4B) :: status
!!$CHARACTER (LEN=256) :: msg
!!$INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
!!$INTEGER :: rankm = 1
!!$INTEGER(HID_T) :: memspace 
!!$INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
!!$INTEGER :: hdferr
!!$
!!$IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!!$   call Display_Message("writeDataset_3DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$IF( ds%dataset_id < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "dataset must be selected before write: ", ds%name
!!$   call Display_Message("writeDataset_3DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$dims(1:ds%rank) = ds%dims(1:ds%rank) 
!!$IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
!!$   call h5dget_space_f(ds%dataset_id, dataspace, error)
!!$   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
!!$   dimsm = SHAPE( data_in )
!!$   call h5screate_simple_f( rankm, dimsm, memspace, error)
!!$ELSE
!!$   dataspace = H5S_ALL_F
!!$   memspace = H5S_ALL_F
!!$ENDIF
!!$
!!$call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
!!$call h5tset_size_f(memtype_id, len(data_in), hdferr)
!!$call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dims, error, &
!!$                 memspace, dataspace )
!!$IF( error < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
!!$   call Display_Message("writeDataset_3DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)
!!$
!!$END FUNCTION writeDataset_3DCHAR
!!$
!!$!...............................................................................
!!$
!!$FUNCTION writeDataset_4DCHAR( ds, data_in, start, count ) RESULT (status)
!!$
!!$! optional keywords:
!!$! start = starting indicies, zero based
!!$! count = number of data points
!!$
!!$TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
!!$CHARACTER(LEN=*), DIMENSION(:,:,:,:), INTENT(IN) :: data_in
!!$INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
!!$INTEGER          :: error
!!$INTEGER (HID_T)  :: dataspace, datatype_id, memtype_id 
!!$INTEGER (I4B) :: status
!!$CHARACTER (LEN=256) :: msg
!!$INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
!!$INTEGER :: rankm = 1
!!$INTEGER(HID_T) :: memspace 
!!$INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
!!$INTEGER :: hdferr
!!$
!!$IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!!$   call Display_Message("writeDataset_4DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$IF( ds%dataset_id < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "dataset must be selected before write: ", ds%name
!!$   call Display_Message("writeDataset_4DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$dims(1:ds%rank) = ds%dims(1:ds%rank) 
!!$IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
!!$   call h5dget_space_f(ds%dataset_id, dataspace, error)
!!$   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
!!$   dimsm = SHAPE( data_in )
!!$   call h5screate_simple_f( rankm, dimsm, memspace, error)
!!$ELSE
!!$   dataspace = H5S_ALL_F
!!$   memspace = H5S_ALL_F
!!$ENDIF
!!$
!!$call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
!!$call h5tset_size_f(memtype_id, len(data_in), hdferr)
!!$call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dims, error, &
!!$                 memspace, dataspace )
!!$IF( error < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
!!$   call Display_Message("writeDataset_4DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)
!!$
!!$END FUNCTION writeDataset_4DCHAR
!!$
!!$!...............................................................................
!!$
!!$FUNCTION writeDataset_5DCHAR( ds, data_in, start, count ) RESULT (status)
!!$
!!$! optional keywords:
!!$! start = starting indicies, zero based
!!$! count = number of data points
!!$
!!$TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
!!$CHARACTER(LEN=*), DIMENSION(:,:,:,:,:), INTENT(IN) :: data_in
!!$INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
!!$INTEGER          :: error
!!$INTEGER (HID_T)  :: dataspace, datatype_id, memtype_id 
!!$INTEGER (I4B) :: status
!!$CHARACTER (LEN=256) :: msg
!!$INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
!!$INTEGER :: rankm = 1
!!$INTEGER(HID_T) :: memspace 
!!$INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
!!$INTEGER :: hdferr
!!$
!!$IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!!$   call Display_Message("writeDataset_5DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$IF( ds%dataset_id < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "dataset must be selected before write: ", ds%name
!!$   call Display_Message("writeDataset_5DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$dims(1:ds%rank) = ds%dims(1:ds%rank) 
!!$IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
!!$   call h5dget_space_f(ds%dataset_id, dataspace, error)
!!$   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
!!$   dimsm = SHAPE( data_in )
!!$   call h5screate_simple_f( rankm, dimsm, memspace, error)
!!$ELSE
!!$   dataspace = H5S_ALL_F
!!$   memspace = H5S_ALL_F
!!$ENDIF
!!$
!!$call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
!!$call h5tset_size_f(memtype_id, len(data_in), hdferr)
!!$call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dims, error, &
!!$                 memspace, dataspace )
!!$IF( error < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
!!$   call Display_Message("writeDataset_5DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)
!!$
!!$END FUNCTION writeDataset_5DCHAR
!!$
!!$!...............................................................................
!!$
!!$FUNCTION writeDataset_6DCHAR( ds, data_in, start, count ) RESULT (status)
!!$
!!$! optional keywords:
!!$! start = starting indicies, zero based
!!$! count = number of data points
!!$
!!$TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
!!$CHARACTER(LEN=*), DIMENSION(:,:,:,:,:,:), INTENT(IN) :: data_in
!!$INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
!!$INTEGER          :: error
!!$INTEGER (HID_T)  :: dataspace, datatype_id, memtype_id 
!!$INTEGER (I4B) :: status
!!$CHARACTER (LEN=256) :: msg
!!$INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
!!$INTEGER :: rankm = 1
!!$INTEGER(HID_T) :: memspace 
!!$INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
!!$INTEGER :: hdferr
!!$
!!$IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!!$   call Display_Message("writeDataset_6DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$IF( ds%dataset_id < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "dataset must be selected before write: ", ds%name
!!$   call Display_Message("writeDataset_6DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$dims(1:ds%rank) = ds%dims(1:ds%rank) 
!!$IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
!!$   call h5dget_space_f(ds%dataset_id, dataspace, error)
!!$   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
!!$   dimsm = SHAPE( data_in )
!!$   call h5screate_simple_f( rankm, dimsm, memspace, error)
!!$ELSE
!!$   dataspace = H5S_ALL_F
!!$   memspace = H5S_ALL_F
!!$ENDIF
!!$
!!$call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
!!$call h5tset_size_f(memtype_id, len(data_in), hdferr)
!!$call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dims, error, &
!!$                 memspace, dataspace )
!!$IF( error < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
!!$   call Display_Message("writeDataset_6DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)
!!$
!!$END FUNCTION writeDataset_6DCHAR
!!$
!!$!...............................................................................
!!$
!!$FUNCTION writeDataset_7DCHAR( ds, data_in, start, count ) RESULT (status)
!!$
!!$! optional keywords:
!!$! start = starting indicies, zero based
!!$! count = number of data points
!!$
!!$TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds 
!!$CHARACTER(LEN=*), DIMENSION(:,:,:,:,:,:,:), INTENT(IN) :: data_in
!!$INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count 
!!$INTEGER          :: error
!!$INTEGER (HID_T)  :: dataspace, datatype_id, memtype_id 
!!$INTEGER (I4B) :: status
!!$CHARACTER (LEN=256) :: msg
!!$INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
!!$INTEGER :: rankm = 1
!!$INTEGER(HID_T) :: memspace 
!!$INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
!!$INTEGER :: hdferr
!!$
!!$IF( ds%datatype_id /= H5T_NATIVE_CHARACTER) then
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "data type should be H5T_NATIVE_CHARACTER: ", ds%name
!!$   call Display_Message("writeDataset_7DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$IF( ds%dataset_id < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "dataset must be selected before write: ", ds%name
!!$   call Display_Message("writeDataset_7DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$
!!$dims(1:ds%rank) = ds%dims(1:ds%rank) 
!!$IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
!!$   call h5dget_space_f(ds%dataset_id, dataspace, error)
!!$   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
!!$   dimsm = SHAPE( data_in )
!!$   call h5screate_simple_f( rankm, dimsm, memspace, error)
!!$ELSE
!!$   dataspace = H5S_ALL_F
!!$   memspace = H5S_ALL_F
!!$ENDIF
!!$
!!$call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
!!$call h5tset_size_f(memtype_id, len(data_in), hdferr)
!!$call h5dwrite_f( ds%dataset_id, memtype_id, data_in, dims, error, &
!!$                 memspace, dataspace )
!!$IF( error < zero ) THEN
!!$   status = FAILURE_STATE
!!$   write( msg,* ) "h5dwrite_f failed on dataset: ", ds%name
!!$   call Display_Message("writeDataset_7DCHAR", msg, status)
!!$   RETURN
!!$ENDIF
!!$IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)
!!$
!!$END FUNCTION writeDataset_7DCHAR

!*******************************************************************************
!
! errorStatus = readDataset_*(ds, data_out, start, count)
!
!*******************************************************************************

FUNCTION readDataset_1DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
INTEGER (I4B), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 1
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
integer :: typeClass

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_1DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_1DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank) 
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dims, error, &
                memspace, dataspace )
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_1DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_1DIK4

!...............................................................................

FUNCTION readDataset_2DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 2
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
integer :: typeClass

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_2DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_2DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_2DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_2DIK4

!...............................................................................

FUNCTION readDataset_3DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 3
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
integer :: typeClass

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_3DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_3DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_3DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_3DIK4

!...............................................................................

FUNCTION readDataset_4DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 4
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
integer :: typeClass

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_4DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_4DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_4DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_4DIK4

!...............................................................................

FUNCTION readDataset_5DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 5
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
integer :: typeClass

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_5DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_5DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_5DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_5DIK4

!...............................................................................

FUNCTION readDataset_6DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 6
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
integer :: typeClass

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_6DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_6DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_6DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_6DIK4

!...............................................................................

FUNCTION readDataset_7DIK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
INTEGER (I4B), DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 7
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
integer :: typeClass

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_7DIK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
IF(typeClass /=  H5T_INTEGER_F) then
   status = FAILURE_STATE
   write( msg,* ) "class type should be H5T_INTEGER_F: ", ds%name
   call Display_Message("readDataset_7DIK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, H5T_NATIVE_INTEGER, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_7DIK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_7DIK4

!...............................................................................


FUNCTION readDataset_1DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(SP), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 1
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_1DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error) 
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_1DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank) 
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_1DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_1DRK4

!...............................................................................

FUNCTION readDataset_2DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 2
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_2DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_2DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_2DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_2DRK4

!...............................................................................

FUNCTION readDataset_3DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 3
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_3DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_3DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_3DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_3DRK4

!...............................................................................

FUNCTION readDataset_4DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 4
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_4DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_4DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_4DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_4DRK4

!...............................................................................

FUNCTION readDataset_5DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 5
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_5DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_5DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_5DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_5DRK4

!...............................................................................

FUNCTION readDataset_6DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 6
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_6DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_6DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_6DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_6DRK4

!...............................................................................

FUNCTION readDataset_7DRK4( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(SP), DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 7
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_7DRK4", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 4) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_REAL: ", ds%name
   call Display_Message("readDataset_7DRK4", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_7DRK4", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_7DRK4

!...............................................................................

FUNCTION readDataset_1DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds 
REAL(DP), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 1
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(1) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_1DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_1DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank) 
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_1DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_1DRK8

!...............................................................................

FUNCTION readDataset_2DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 2
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(2) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_2DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_2DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_2DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_2DRK8

!...............................................................................

FUNCTION readDataset_3DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 3
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(3) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_3DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_3DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_3DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_3DRK8

!...............................................................................

FUNCTION readDataset_4DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 4
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(4) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_4DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_4DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_4DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_4DRK8

!...............................................................................

FUNCTION readDataset_5DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 5
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(5) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_5DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_5DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_5DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_5DRK8

!...............................................................................

FUNCTION readDataset_6DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 6
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(6) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_6DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_6DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_6DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_6DRK8

!...............................................................................

FUNCTION readDataset_7DRK8( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_T), INTENT( IN ) :: ds
REAL(DP), DIMENSION(:,:,:,:,:,:,:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER          :: error
INTEGER (HID_T)  :: dataspace 
INTEGER (I4B) :: status
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(MAXRANK) :: dims
INTEGER :: rankm = 7
INTEGER(HID_T) :: memspace 
INTEGER(HSIZE_T), DIMENSION(7) :: dimsm
integer :: typeClass
integer(HSIZE_T) :: typeSize

status = SUCCESS_STATE

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_7DRK8", msg, status)
   RETURN
ENDIF

call h5tget_class_f(ds%datatype_id, typeClass, error)
call h5tget_size_f(ds%datatype_id, typeSize, error)
IF(typeClass /=  H5T_FLOAT_F .or. typeSize /= 8) then
   status = FAILURE_STATE
   write( msg,* ) "data type should be H5T_NATIVE_DOUBLE: ", ds%name
   call Display_Message("readDataset_7DRK8", msg, status)
   RETURN
ENDIF

dims(1:ds%rank) = ds%dims(1:ds%rank)
IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm = SHAPE( data_out )
   CALL h5screate_simple_f( rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

CALL h5dread_f( ds%dataset_id, ds%datatype_id, data_out, dims, error, &
                memspace, dataspace )

IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_7DRK8", msg, status)
   RETURN
ENDIF
IF( memspace .NE. H5S_ALL_F ) CALL h5sclose_f(memspace, error)

END FUNCTION readDataset_7DRK8

!...............................................................................

FUNCTION readDataset_1DCHAR( ds, data_out, start, count ) RESULT (status)

TYPE (H5SDS_CHAR_T), INTENT( IN ) :: ds
CHARACTER(LEN=*), DIMENSION(:), INTENT(OUT) :: data_out
INTEGER (HSIZE_T), DIMENSION(:), OPTIONAL, INTENT(IN) :: start, count
INTEGER :: error
INTEGER (I4B) :: status, rankm
INTEGER(HID_T) :: memtype_id, dataspace, memspace
CHARACTER (LEN=256) :: msg
INTEGER(HSIZE_T), DIMENSION(1) :: dims, dimsm
INTEGER(HSIZE_T) :: len_data_out

status = SUCCESS_STATE

rankm = 1

IF( ds%dataset_id < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "dataset must be selected before read: ", ds%name
   call Display_Message("readDataset_1DCHAR", msg, status)
   RETURN
ENDIF

IF( PRESENT(count) .AND. PRESENT( start ) ) THEN
   call h5dget_space_f(ds%dataset_id, dataspace, error)
   call h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
   dimsm(1) = size(data_out)
   call h5screate_simple_f(rankm, dimsm, memspace, error)
ELSE
   dataspace = H5S_ALL_F
   memspace = H5S_ALL_F
ENDIF

call h5tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, error)
len_data_out=len(data_out)
call h5tset_size_f(memtype_id, len_data_out, error)
CALL h5dread_f(ds%dataset_id,memtype_id,data_out,dims,error,memspace,dataspace)
IF( error < zero ) THEN
   status = FAILURE_STATE
   WRITE( msg,* ) "h5dread_f failed on dataset: ", ds%name
   call Display_Message("readDataset_1DCHAR", msg, status)
   RETURN
ENDIF

IF( memspace .NE. H5S_ALL_F ) call h5sclose_f(memspace, error)

END FUNCTION readDataset_1DCHAR

!
!//////////////////////// P R I V A T E   M E T H O D S ////////////////////////
!

subroutine readDatasetAttributes(sds, hdferr)

implicit none
type(H5SDS_T), intent(inout) :: sds
integer, intent(out) :: hdferr

integer(HID_T) :: memtype_id, attr_id
integer(SIZE_T) :: attr_len
integer(HSIZE_T), dimension(1) :: dims 

character(len=128) :: attr_name
logical :: attr_exists

!
! long_name:
!
attr_name = "long_name"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
   attr_len = len(sds%long_name)
   call H5Tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
   call h5tset_size_f(memtype_id, attr_len, hdferr)
   call h5aread_f(attr_id, memtype_id,  sds%long_name, dims, hdferr) 
   call h5tclose_f(memtype_id, hdferr)
   call h5aclose_f(attr_id, hdferr)
endif

!
! units:
!
attr_name = "units"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
   attr_len = len(sds%units)
   call H5Tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
   call h5tset_size_f(memtype_id, attr_len, hdferr)
   call h5aread_f(attr_id, memtype_id,  sds%units, dims, hdferr) 
   call h5tclose_f(memtype_id, hdferr)
   call h5aclose_f(attr_id, hdferr)
endif

!
! description:
!
attr_name = "description"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
   attr_len = len(sds%description)
   call H5Tcopy_f(H5T_NATIVE_CHARACTER, memtype_id, hdferr)
   call h5tset_size_f(memtype_id, attr_len, hdferr)
   call h5aread_f(attr_id, memtype_id,  sds%description, dims, hdferr) 
   call h5tclose_f(memtype_id, hdferr)
   call h5aclose_f(attr_id, hdferr)
endif

!
! valid range:
!
! attr_name = "ValidRange"
! call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
! if(attr_exists) then
!    call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
!    dims(1) = 2
!    call h5aread_f(attr_id, H5T_NATIVE_DOUBLE,  sds%valid_range, dims, hdferr) 
!    call h5aclose_f(attr_id, hdferr)
! endif

!
! valid min:
!
attr_name = "valid_min"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
!   call h5aread_f(attr_id, H5T_NATIVE_DOUBLE,  sds%valid_min, dims, hdferr) 
   call h5aclose_f(attr_id, hdferr)
endif

!
! valid max:
!
attr_name = "valid_max"
call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
if(attr_exists) then
   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
   dims(1) = 1
!   call h5aread_f(attr_id, H5T_NATIVE_DOUBLE,  sds%valid_max, dims, hdferr) 
   call h5aclose_f(attr_id, hdferr)
endif

!
! fill value attribute:
!
!attr_name = "_FillValue"
!call h5aexists_f(sds%dataset_id, attr_name, attr_exists, hdferr)
!if(attr_exists) then
!   call h5aopen_f(sds%dataset_id, attr_name, attr_id, hdferr)
!   dims(1) = 1
!   call h5aread_f(attr_id, H5T_NATIVE_DOUBLE,  sds%fill_value, dims, hdferr) 
!   call h5aclose_f(attr_id, hdferr)
!endif

end subroutine readDatasetAttributes

!*******************************************************************************

subroutine writeDatasetAttributes(sds, hdferr)

! Description: 
!    create dataset attributes: long_name, units, valid range and fill value. Very
!    similar to the set of attributes for HDFEOS datasets. 

implicit none
type(H5SDS_T), intent(in) :: sds
integer, intent(out) :: hdferr

integer(HID_T) :: datatype_id, dataspace_id, attr_id
integer(SIZE_T) :: attr_len
integer(HSIZE_T), dimension(1) :: dims 
integer :: rank

rank = 1

!
! long_name:
!
dims(1) = 1
attr_len = len_trim(sds%long_name)
call H5Screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
call H5Tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
call H5Tset_strpad_f(datatype_id, H5T_STR_NULLTERM_F, hdferr)
call H5Tset_size_f(datatype_id, attr_len, hdferr)
call H5Acreate_f(sds%dataset_id, "long_name",&
                 datatype_id, dataspace_id, attr_id, hdferr)
call H5Awrite_f(attr_id, datatype_id, sds%long_name, dims, hdferr)
call H5Aclose_f(attr_id, hdferr)
call H5Tclose_f(datatype_id, hdferr)
call H5Sclose_f(dataspace_id, hdferr)

!
! units:
!
dims(1) = 1
attr_len = len_trim(sds%units)
call H5Screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
call H5Tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
call H5Tset_strpad_f(datatype_id, H5T_STR_NULLTERM_F, hdferr)
call H5Tset_size_f(datatype_id, attr_len, hdferr)
call H5Acreate_f(sds%dataset_id, "units", &
                 datatype_id, dataspace_id, attr_id, hdferr)
call H5Awrite_f(attr_id, datatype_id, sds%units, dims, hdferr)
call H5Aclose_f(attr_id, hdferr)
call H5Tclose_f(datatype_id, hdferr)
call H5Sclose_f(dataspace_id, hdferr)

!
! description:
!
dims(1) = 1
attr_len = len_trim(sds%description)
call H5Screate_f(H5S_SCALAR_F, dataspace_id, hdferr)
call H5Tcopy_f(H5T_NATIVE_CHARACTER, datatype_id, hdferr)
call H5Tset_strpad_f(datatype_id, H5T_STR_NULLTERM_F, hdferr)
call H5Tset_size_f(datatype_id, attr_len, hdferr)
call H5Acreate_f(sds%dataset_id, "description",&
                 datatype_id, dataspace_id, attr_id, hdferr)
call H5Awrite_f(attr_id, datatype_id, sds%description, dims, hdferr)
call H5Aclose_f(attr_id, hdferr)
call H5Tclose_f(datatype_id, hdferr)
call H5Sclose_f(dataspace_id, hdferr)

!
! valid range:
!
! dims(1) = 2
! call H5Screate_simple_f(rank, dims, dataspace_id, hdferr)
! call h5Tcopy_f(sds%datatype_id, datatype_id, hdferr)
! call H5Acreate_f(sds%dataset_id, "ValidRange", &
!                  datatype_id, dataspace_id, attr_id, hdferr)
! call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, sds%valid_range, dims, hdferr)
! call H5Aclose_f(attr_id, hdferr)
! call H5Tclose_f(datatype_id, hdferr)
! call H5Sclose_f(dataspace_id, hdferr)

!if(sds%fill_value < 999.9) then
!
! valid min:
!
  dims(1) = 1
  call H5Screate_simple_f(rank, dims, dataspace_id, hdferr)
  call h5Tcopy_f(sds%datatype_id, datatype_id, hdferr)
!  call H5Acreate_f(sds%dataset_id, "valid_min", &
!                   datatype_id, dataspace_id, attr_id, hdferr)
!  call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, sds%valid_min, dims, hdferr)
  call H5Aclose_f(attr_id, hdferr)
  call H5Tclose_f(datatype_id, hdferr)
  call H5Sclose_f(dataspace_id, hdferr)

!
! valid max:
!
  dims(1) = 1
  call H5Screate_simple_f(rank, dims, dataspace_id, hdferr)
  call h5Tcopy_f(sds%datatype_id, datatype_id, hdferr)
!  call H5Acreate_f(sds%dataset_id, "valid_max", &
!                   datatype_id, dataspace_id, attr_id, hdferr)
!  call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, sds%valid_max, dims, hdferr)
  call H5Aclose_f(attr_id, hdferr)
  call H5Tclose_f(datatype_id, hdferr)
  call H5Sclose_f(dataspace_id, hdferr)

!
! fill value attribute:
!
!  dims(1) = 1
!  call H5Screate_simple_f(rank, dims, dataspace_id, hdferr)
!  call h5Tcopy_f(sds%datatype_id, datatype_id, hdferr)
!  call H5Acreate_f(sds%dataset_id, "_FillValue", &
!                   datatype_id, dataspace_id, attr_id, hdferr)
!  call H5Awrite_f(attr_id, H5T_NATIVE_DOUBLE, sds%fill_value, dims, hdferr)
!  call H5Aclose_f(attr_id, hdferr)
!  call H5Tclose_f(datatype_id, hdferr)
!  call H5Sclose_f(dataspace_id, hdferr)
!end if

end subroutine writeDatasetAttributes

!*******************************************************************************

END MODULE H5Util_class

