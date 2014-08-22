       module m_write_swath_field
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_write_swath_field
! 
! !DESCRIPTION: write a swath field to an EOS-HDF swath file
!
! !CALLING SEQUENCE: 
!
!        status=put_data(swid,fieldname,dimname,data)
!     
! !INPUT PARAMETERS:   
!               integer (kind = 4) swid : swath id returned by SWattach
!               character(len=*) fieldname : name of field to write
!
! !OUTPUT PARAMETERS:  
!               integer (kind = 4) status : status of the write, 0 = good
!               generic data     : array of data to write
!
! !SEE ALSO:  
!
! !REVISION HISTORY: 
!
!  05Jun01   J. Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------
!

       interface put_data
         module procedure put_data_1dr8
         module procedure put_data_1dr4
         module procedure put_data_1di2
!        module procedure put_data_1di1
         module procedure put_data_2dr4
         module procedure put_data_2di2
         module procedure put_data_2di1
         module procedure put_data_3dr8
         module procedure put_data_3dr4
         module procedure put_data_3di2
!        module procedure put_data_3di1
       end interface

       include 'hdfeos5.inc'

!       integer (kind = 4), parameter :: HE5_HDFE_NOMERGE=0
!       integer (kind = 4), parameter :: HE5_HDFE_AUTOMERGE=1
       integer (kind = 4), parameter :: merge=HE5_HDFE_NOMERGE
!       integer (kind = 4), parameter :: HE5T_NATIVE_DOUBLE=11
!       integer (kind = 4), parameter :: HE5T_NATIVE_FLOAT=10
!       integer (kind = 4), parameter :: HE5T_NATIVE_INT=0
!       integer (kind = 4), parameter :: HE5T_NATIVE_INT8=13
!       integer (kind = 4), parameter :: HE5T_NATIVE_UINT8=14
!       integer (kind = 4), parameter :: HE5T_NATIVE_INT16=15
!       integer (kind = 4), parameter :: HE5T_NATIVE_UINT16=16
!       integer (kind = 4), parameter :: HE5T_NATIVE_CHAR=56
       real    (kind = 8), parameter :: off_set_default=0.
       real    (kind = 8), parameter :: scale_factor_default=1.
       integer (kind = 4) :: nn
       contains

       function put_data_1dr8(swid,fieldname,dimname,data, &
          missingvalue, title, units, geo, offset, iprt) &
          result(status)

       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       integer,   intent(in), optional :: offset
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       real (kind = 8), dimension(:) :: data
       real (kind = 8), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=1
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld, append
       integer (KIND=C_LONG) :: nn, n1=1, n11=11

       geofld=.false.
       if (present(iprt)) iprt2=iprt
       if (present(geo)) geofld=geo
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       edge = dims(1:dim)
       !numbertype = DFNT_FLOAT64
       numbertype = HE5T_NATIVE_DOUBLE
       append=.false.
       if (present(offset)) then
        if (offset > 0) then
         append = .true.
         start = offset
        endif
       endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname,  &
          dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname,  &
          dimname, " ", numbertype, merge)
        endif
       endif
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, 1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR, nn, units)
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n11, "Aura-Shared")
       status = he5_swwrfld (swid, fieldname, &
          start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1),data(dims(1))
       end function put_data_1dr8

       function put_data_1dr4(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)
       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       integer,   intent(in), optional :: offset
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       real (kind = 4), dimension(:) :: data
       real (kind = 4), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=1
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld, append
       integer (KIND=C_LONG) :: nn, n1=1, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       edge = dims(1:dim)
       numbertype = HE5T_NATIVE_FLOAT
       append=.false.
       if (present(offset)) then
        if (offset > 0) then
         append = .true.
         start = offset
        endif
       endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        endif
       endif
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, n1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR, nn, units)
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1),data(dims(1))
       end function put_data_1dr4




       function put_data_1di2(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)
       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       logical, intent(in), optional :: geo
       integer,   intent(in), optional :: offset
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       integer (kind = 2), dimension(:) :: data
       integer (kind = 2), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=1
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld, append
       integer (KIND=C_LONG) :: nn, n1=1, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       edge = dims(1:dim)
       numbertype = HE5T_NATIVE_UINT16
       append=.false.
       if (present(offset)) then
        if (offset > 0) then
         append = .true.
         start = offset
        endif
       endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        endif
       endif
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, n1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR, nn, units)
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1),data(dims(1))
       end function put_data_1di2


!      function put_data_1di1(swid,fieldname,dimname,data,geo) result(status)
!      implicit none
!      !include 'hdf.f90'

!      integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld
!      integer (kind = 4), intent(in) :: swid
!      logical, intent(in), optional :: geo
!      character(len=*), intent(in) :: fieldname
!      integer (kind = 1), dimension(:) :: data
!      integer (kind = 4) :: status, numbertype
!      integer, parameter :: dim=1
!      integer (kind = 4), dimension(dim) :: dims
!      integer (kind = 4) start(dim), stride(dim), edge(dim)
!      character(len=*), intent(in) :: dimname
!      logical :: geofld

!      geofld=.false.
!      if (present(geo)) geofld=geo
!      start = 0
!      stride = 1
!      dims(1)=size(data,dim=1)
!      edge = dims(1:dim)
!      numbertype = DFNT_INT8
!      if (geofld) then
!        status = he5_swdefgfld (swid, fieldname, & 
!          dimname, " ", numbertype, merge)
!      else
!        status = he5_swdefdfld (swid, fieldname, & 
!          dimname, " ", numbertype, merge)
!      endif
!      status = he5_swwrfld (swid, fieldname, &
!          start, stride, edge, data)
!      write (6, *) fieldname,status,data(1),data(dims(1))
!      end function put_data_1di1


!      function put_data_2di1(swid,fieldname,dimname,data,geo) result(status)
!      implicit none
!      !include 'hdf.f90'

!      integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld
!      integer (kind = 4), intent(in) :: swid
!      logical, intent(in), optional :: geo
!      character(len=*), intent(in) :: fieldname
!      integer (kind = 1), dimension(:,:) :: data
!      integer (kind = 4) :: status, numbertype
!      integer, parameter :: dim=2
!      integer (kind = 4), dimension(dim) :: dims
!      integer (kind = 4) start(dim), stride(dim), edge(dim)
!      character(len=*), intent(in) :: dimname
!      logical :: geofld

!      geofld=.false.
!      if (present(geo)) geofld=geo
!      start = 0
!      stride = 1
!      dims(1)=size(data,dim=1)
!      dims(2)=size(data,dim=2)
!      edge = dims(1:dim)
!      numbertype = DFNT_INT8
!      if (geofld) then
!        status = he5_swdefgfld (swid, fieldname,  &
!          dimname, " ", numbertype, merge)
!      else
!        status = he5_swdefdfld (swid, fieldname,  &
!          dimname, " ", numbertype, merge)
!      endif
!      status = he5_swwrfld (swid, fieldname, &
!          start, stride, edge, data)
!      write (6, *) fieldname,status,data(1,1),data(dims(1),dims(2))
!      end function put_data_2di1


       function put_data_2di2(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)
       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       integer (kind = 2), dimension(:,:) :: data
       integer (kind = 2), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=2
       integer, dimension(dim), intent(in), optional :: offset
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld, append
       integer (KIND=C_LONG) :: nn, n1=1, n11=11, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       dims(2)=size(data,dim=2)
       edge = dims(1:dim)
       numbertype = HE5T_NATIVE_UINT16 
       append=.false.
       if (present(offset)) then
        if (any(offset > 0)) then
         append = .true.
         start = offset
        endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       endif
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        endif
       endif
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, n1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR, nn, units)
if(index(fieldname,'TerrainHeight')==1) then
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n11, "Aura-Shared")
else
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
endif
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1,1),data(dims(1),dims(2))
       end function put_data_2di2

       function put_data_2di1(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)
       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       integer (kind = 1), dimension(:,:) :: data
       integer (kind = 1), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=2
       integer, dimension(dim), intent(in), optional :: offset
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld, append
       integer (KIND=C_LONG) :: nn, n1=1, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       dims(2)=size(data,dim=2)
       edge = dims(1:dim)
       numbertype = HE5T_NATIVE_UINT8
       append=.false.
       if (present(offset)) then
        if (any(offset > 0)) then
         append = .true.
         start = offset
        endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       endif
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname, &
           dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname, &
           dimname, " ", numbertype, merge)
        endif
       endif
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR, nn, units)
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1,1),data(dims(1),dims(2))
       end function put_data_2di1

       function put_data_2dr4(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)

       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       real (kind = 4), dimension(:,:) :: data
       real (kind = 4), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=2
       integer, dimension(dim), intent(in), optional :: offset
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld, append
       integer (KIND=C_LONG), dimension(dim) :: chunk_dim
       integer (kind=4) :: chunk_rank=2
       integer (kind=4), dimension(5) :: compparm
       integer (kind=4) :: he5_swdefchunk, he5_swdefcomch
       INTEGER (KIND=C_LONG) :: n1=1, nn, n11=11, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       dims(2)=size(data,dim=2)
       edge = dims(1:dim)
       numbertype = HE5T_NATIVE_FLOAT
       append=.false.
       if (present(offset)) then
        if (any(offset > 0)) then
         append = .true.
         start = offset
        endif
       endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
chunk_dim(1)=dims(1)
chunk_dim(2)=dims(2)
status = he5_swdefchunk(swid, chunk_rank, chunk_dim)
compparm(1) = 8
compparm(2) = 8
compparm(3) = 8
compparm(4) = 8
compparm(5) = 8
status = he5_swdefcomch(swid, HE5_HDFE_COMP_DEFLATE, compparm, chunk_rank, chunk_dim)

       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname,  &
           dimname, " ", numbertype, merge)
!         write (6,*) status, 'after define geo ',fieldname
        else
         status = he5_swdefdfld (swid, fieldname,  &
           dimname, " ", numbertype, merge)
!         write (6,*) status, 'after define dat ',fieldname
        endif
       endif
!       call pzeitbeg('wr2dfl')
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, 1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR, nn, units)
if(index(fieldname,'Latitude')==1 .or. index(fieldname,'Longitude')==1) then
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n11, "Aura-Shared")
else
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
endif
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
!       call pzeitend
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1,1),data(dims(1),dims(2))
       end function put_data_2dr4

       function put_data_3dr8(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)
       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       real (KIND=8), dimension(:,:,:) :: data
       real (kind = 4), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=3
       integer, dimension(dim), intent(in), optional :: offset
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld,append
       INTEGER (KIND=C_LONG) :: n1=1, nn, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       dims(2)=size(data,dim=2)
       dims(3)=size(data,dim=3)
       edge = dims(1:dim)
       !numbertype = DFNT_FLOAT32
       numbertype = HE5T_NATIVE_FLOAT
       append=.false.
       if (present(offset)) then
        if (any(offset > 0)) then
         append = .true.
         start = offset
        endif
       endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname,  &
           dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname,  &
           dimname, " ", numbertype, merge)
        endif
       endif
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, n1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR,nn, units)
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1,1,1), &
           data(dims(1),dims(2),dims(3))
       end function put_data_3dr8

       function put_data_3dr4(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)
       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       real (kind = 4), dimension(:,:,:) :: data
       real (kind = 4), intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=3
       integer, dimension(dim), intent(in), optional :: offset
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld,append
       integer (KIND=C_LONG) :: nn, n1=1, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       dims(2)=size(data,dim=2)
       dims(3)=size(data,dim=3)
       edge = dims(1:dim)
       numbertype = HE5T_NATIVE_FLOAT
       append=.false.
       if (present(offset)) then
        if (any(offset > 0)) then
         append = .true.
         start = offset
        endif
       endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname,  &
           dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname,  &
           dimname, " ", numbertype, merge)
        endif
       endif
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, n1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR,nn, units)
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1,1,1), &
           data(dims(1),dims(2),dims(3))
       end function put_data_3dr4

       function put_data_3di2(swid,fieldname,dimname,data, &
          missingvalue, title, units,geo,offset,iprt) &
          result(status)
       USE ISO_C_BINDING, ONLY: C_LONG
       implicit none
       integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld, he5_swwrlattr, he5_swsetfill
       integer (kind = 4), intent(in) :: swid
       logical, intent(in), optional :: geo
       character(len=*), intent(in) :: fieldname
       integer,   intent(in), optional :: iprt
       integer :: iprt2=0
       integer , dimension(:,:,:) :: data
       integer , intent(in) :: missingvalue
       integer (kind = 4) :: status, numbertype
       integer, parameter :: dim=3
       integer, dimension(dim), intent(in), optional :: offset
       integer (kind = 4), dimension(dim) :: dims
!       integer (kind = 4) start(dim), stride(dim), edge(dim)
       integer (KIND=C_LONG) start(dim), stride(dim), edge(dim)
       character(len=*), intent(in) :: dimname, title, units
       logical :: geofld, append
       integer (KIND=C_LONG) ::  nn, n1=1, n12=12

       geofld=.false.
       if (present(geo)) geofld=geo
       if (present(iprt)) iprt2=iprt
       start = 0
       stride = 1
       dims(1)=size(data,dim=1)
       dims(2)=size(data,dim=2)
       dims(3)=size(data,dim=3)
       edge = dims(1:dim)
       numbertype = HE5T_NATIVE_UINT16
       append=.false.
       if (present(offset)) then
        if (any(offset > 0)) then
         append = .true.
         start = offset
        endif
       endif
       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       if (.not. append) then
        if (geofld) then
         status = he5_swdefgfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        else
         status = he5_swdefdfld (swid, fieldname, & 
           dimname, " ", numbertype, merge)
        endif
       endif
!       status = he5_swwrlattr(swid, fieldname, "_FillValue", numbertype, n1, missingvalue)
!       status = he5_swsetfill(swid, fieldname, numbertype, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "MissingValue", numbertype, n1, missingvalue)
       status = he5_swwrlattr(swid, fieldname, "Offset", HE5T_NATIVE_DOUBLE, n1, off_set_default)
       status = he5_swwrlattr(swid, fieldname, "ScaleFactor", HE5T_NATIVE_DOUBLE, n1, scale_factor_default)
       nn = len_trim(title)
       status = he5_swwrlattr(swid, fieldname, "Title", HE5T_NATIVE_CHAR, nn, title)
       nn = len_trim(units)
       status = he5_swwrlattr(swid, fieldname, "Units", HE5T_NATIVE_CHAR, nn, units)
       status = he5_swwrlattr(swid, fieldname, "UniqueFieldDefinition", HE5T_NATIVE_CHAR, n12, "OMI-Specific")
       status = he5_swwrfld (swid, fieldname, &
           start, stride, edge, data)
       if (iprt2 >= 2) &
       write (6, *) fieldname,status,data(1,1,1), &
           data(dims(1),dims(2),dims(3))
       end function put_data_3di2

!      function put_data_3di1(swid,fieldname,dimname,data,geo) result(status)
    
!      implicit none
!      !include 'hdf.f90'

!      integer (kind = 4) he5_swwrfld, he5_swdefdfld, he5_swdefgfld
!      integer (kind = 4), intent(in) :: swid
!      logical, intent(in), optional :: geo
!      character(len=*), intent(in) :: fieldname
!      integer (kind = 1), dimension(:,:,:) :: data
!      integer (kind = 4) :: status, numbertype
!      integer, parameter :: dim=3
!      integer (kind = 4), dimension(dim) :: dims
!      integer (kind = 4) start(dim), stride(dim), edge(dim)
!      character(len=*), intent(in) :: dimname
!      logical :: geofld

!      geofld=.false.
!      if (present(geo)) geofld=geo
!      start = 0
!      stride = 1
!      dims(1)=size(data,dim=1)
!      dims(2)=size(data,dim=2)
!      dims(3)=size(data,dim=3)
!      edge = dims(1:dim)
!      numbertype = DFNT_INT8
!      if (geofld) then
!        status = he5_swdefgfld (swid, fieldname, & 
!        dimname, " ", numbertype, merge)
!      else
!        status = he5_swdefdfld (swid, fieldname, & 
!        dimname, " ", numbertype, merge)
!      endif
!      status = he5_swwrfld (swid, fieldname, &
!          start, stride, edge, data)
!      write (6, *) fieldname,status,data(1,1,1), &
!          data(dims(1),dims(2),dims(3))
!      end function put_data_3di1

       end module m_write_swath_field
