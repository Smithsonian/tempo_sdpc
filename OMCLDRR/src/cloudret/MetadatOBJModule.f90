 module MetadatOBJModule
!---------------------------------------------------------------------
! This module provides the Object Names for Metadata               
!--------------------------------------------------------------------- 

 IMPLICIT NONE
 INCLUDE 'PGS_MET.f'

! INTEGER, PARAMETER :: ninvname=13
 INTEGER, PARAMETER :: ninvname=11
 CHARACTER(LEN=PGSd_MET_GROUP_NAME_L),DIMENSION(ninvname), PARAMETER :: &  
              INVOBJ = (/                                     &
!                         "LOCALGRANULEID                   ", &
                         "AUTOMATICQUALITYFLAGEXPLANATION.1", &
                         "AUTOMATICQUALITYFLAG.1           ", &
                         "QAPERCENTMISSINGDATA.1           ", &
!                         "QAPERCENTOUTOFBOUNDSDATA.1       ", &
                         "PARAMETERNAME.1                  ", &
                         "EQUATORCROSSINGDATE.1            ", &
                         "EQUATORCROSSINGTIME.1            ", &
                         "INPUTPOINTER                     ", &
                         "RANGEENDINGDATE                  ", &
                         "RANGEENDINGTIME                  ", &
                         "RANGEBEGINNINGDATE               ", &
                         "RANGEBEGINNINGTIME               "/)

END module MetadatOBJModule
