#include <HE5_HdfEosDef.h>
#include <cfortHdf.h>

/* Given HE5 datatype, find out the corresponding HDF5 data type, and then
   use HDF5 intrinsic function H5Tget_size to figure out its size in byte,
   return this size */
int
HE5Tget_size( int fdatatype )
{ hid_t  cdatatype = FAIL;            /* return data type ID  */
  size_t size;
  cdatatype = HE5_EHconvdatatype( fdatatype ); 
  if( cdatatype == FAIL ) 
     size = -1;
  else 
     size = (int) H5Tget_size( cdatatype );
  return size;
} 

FCALLSCFUN1( INT, HE5Tget_size, HE5TGET_SIZE, he5tget_size, INT )



