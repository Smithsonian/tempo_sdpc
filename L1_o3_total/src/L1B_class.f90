!! F90
 ! 
 ! Description:
 ! MODULE L1B_class includes L1B_radgeo_class and L1B_irrad_class,
 ! which in turn contains data structures for storeing radiance 
 ! and irradiance data retrived from OMI L1B files (in HE4 format).
 ! These classes also contains the functions for the creation, deletion 
 ! and access of specific information stored in the data structures.  
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
MODULE L1B_class
   USE L1B_geoang_class
   USE L1B_radirr_class
!   USE L1B_smlpix_class
END MODULE L1B_class
