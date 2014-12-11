MODULE HDF4_class

  INCLUDE 'hdf.f90'

  !hdf functions
   
  integer  sfstart
  external sfstart 

  integer  sffinfo
  external sffinfo

  integer  sfwdata 
  external sfwdata

  integer  sfrdata 
  external sfrdata

  integer  sfendacc
  external sfendacc

  integer  sfend
  external sfend

  integer  sfn2index
  external sfn2index

  integer  sfselect
  external sfselect

  integer  sfginfo
  external sfginfo

  integer  sfcreate
  external sfcreate
 
  integer  sfscatt
  external sfscatt

  integer  sfsdtstr
  external sfsdtstr

  integer  sfscal
  external sfscal

  integer  sfsrange
  external sfsrange

  integer  sfsfill 
  external sfsfill 

  integer  sfdimid
  external sfdimid

  integer  sfsdmname
  external sfsdmname
 
  integer  sfgainfo
  external sfgainfo

  integer  sfrcdata
  external sfrcdata

END MODULE HDF4_class
