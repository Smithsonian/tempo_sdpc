MODULE OMI_SMF_class
    IMPLICIT NONE
    INCLUDE 'PGS_SMF.f'
    INCLUDE 'PGS_OZT_52050.f'          !!this external file defines the 
                                       !!L1B PGS error codes

    INTEGER (KIND=4) :: OMI_SMF_setmsg !!this external function writes the PGS
    EXTERNAL            OMI_SMF_setmsg !!SMF error messages to the Log files
END MODULE OMI_SMF_class
