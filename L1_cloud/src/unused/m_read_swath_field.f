       module m_read_swath_field
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  m_read_swath_field
! 
! !DESCRIPTION: reads a swath field from an EOS-HDF swath file
!
! !CALLING SEQUENCE: 
!
!        status=get_data(swid,fieldname,data)
!     
! !INPUT PARAMETERS:   
!               integer (kind = 4) swid : swath id returned by SWattach
!               character(len=*) fieldname : name of field to read
!
! !OUTPUT PARAMETERS:  
!               integer (kind = 4) status : status of the read, 0 = good
!               generic data     : array to store data
!
! !SEE ALSO:  
!
! !REVISION HISTORY: 
!
!  05Jun01   J. Joiner     original fortran 90
!
!EOP
!-------------------------------------------------------------------------

       real (kind=4), dimension(:,:,:), pointer :: dummy_3Dr4
       integer (kind=1), dimension(:,:,:), pointer :: dummy_3Di1
       integer (kind=2), dimension(:,:,:), pointer :: dummy_3Di2
!
       public get_data_3dr4, get_data_3di2, get_data_3di1

       interface get_data
         module procedure get_data_1dr8
         module procedure get_data_1dr4
         module procedure get_data_1di2
         module procedure get_data_1di1
         module procedure get_data_2dr4
         module procedure get_data_2di2
         module procedure get_data_2di1
         !module procedure get_data_3dr4
         !module procedure get_data_3di2
         !module procedure get_data_3di1
       end interface

       contains

       function get_data_1dr8(swid,fieldname,data) result(status)

       implicit none
       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       real (kind = 8), dimension(:), pointer :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=1
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       allocate(data(dims(1)))
       data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, data)
!       write (6, *) fieldname,status,data(1),data(dims(1))
       end function get_data_1dr8

       function get_data_1dr4(swid,fieldname,data) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       real (kind = 4), dimension(:), pointer :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=1
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       allocate(data(dims(1)))
       data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, data)
!       write (6, *) fieldname,status,data(1),data(dims(1))
       end function get_data_1dr4

       function get_data_1di2(swid,fieldname,data) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       integer (kind = 2), dimension(:) :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=1
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
!      allocate(data(dims(1)))
!      data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, data)
!       write (6, *) fieldname,status,data(1),data(dims(1))
       end function get_data_1di2


       function get_data_1di1(swid,fieldname,data) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       integer (kind = 1), dimension(:) :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=1
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
!      allocate(data(dims(1)))
!      data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, data)
!       write (6, *) fieldname,status,data(1),data(dims(1))
       end function get_data_1di1


       function get_data_2di1(swid,fieldname,data) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       integer (kind = 1), dimension(:,:) :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=2
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       !allocate(data(dims(1),dims(2)))
       !data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, data)
!       write (6, *) fieldname,status,data(1,1),data(dims(1),dims(2))
       end function get_data_2di1


       function get_data_2di2(swid,fieldname,data) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       integer (kind = 2), dimension(:,:) :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=2
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       !allocate(data(dims(1),dims(2)))
       !data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, data)
!       write (6, *) fieldname,status,data(1,1),data(dims(1),dims(2))
       end function get_data_2di2

       function get_data_2dr4(swid,fieldname,data) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       real (kind = 4), dimension(:,:) :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=2
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       !allocate(data(dims(1),dims(2)))
       !data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, data)
!       write (6, *) fieldname,status,data(1,1),data(dims(1),dims(2))
       end function get_data_2dr4

       function get_data_3dr4(swid,fieldname) result(status)
       !function get_data_3dr4(swid,fieldname,data) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
!      real (kind = 4), dimension(:,:,:), pointer :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=3
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       allocate(dummy_3Dr4(dims(1),dims(2),dims(3)))
!      data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, dummy_3Dr4)
!       write (6, *) fieldname,status,data(1,1,1),
!     .     data(dims(1),dims(2),dims(3))
       end function get_data_3dr4

!      function get_data_3di2(swid,fieldname,data) result(status)
       function get_data_3di2(swid,fieldname) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
!      integer (kind = 2), dimension(:,:,:), pointer :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=3
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       allocate(dummy_3Di2(dims(1),dims(2),dims(3)))
!      data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, dummy_3Di2)
!       write (6, *) fieldname,status,data(1,1,1),
!     .     data(dims(1),dims(2),dims(3))
       end function get_data_3di2

       !function get_data_3di1(swid,fieldname,data) result(status)
       function get_data_3di1(swid,fieldname) result(status)
       implicit none

       integer (kind = 4) swrdfld, swfldinfo
       integer (kind = 4), intent(in) :: swid
       character(len=*), intent(in) :: fieldname
       !integer (kind = 1), dimension(:,:,:), pointer :: data
       integer (kind = 4) :: status
       integer (kind = 4) :: rank, numbertype
       integer, parameter :: dim=3
       integer (kind = 4), dimension(dim) :: dims
       integer (kind = 4) start(dim), stride(dim), edge(dim)
       character(len=64), dimension(dim) :: dimname

       status = swfldinfo(swid, fieldname, rank,
     .   dims, numbertype, dimname)
       allocate(dummy_3Di1(dims(1),dims(2),dims(3)))
       !data = 0
       start = 0
       stride = 1
       edge = dims(1:dim)
       status = swrdfld (swid, fieldname,
     *     start, stride, edge, dummy_3Di1)
!       write (6, *) fieldname,status,data(1,1,1),
!     .     data(dims(1),dims(2),dims(3))
       end function get_data_3di1

       end module m_read_swath_field
