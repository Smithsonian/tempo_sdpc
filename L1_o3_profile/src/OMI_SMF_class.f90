
MODULE OMI_SMF_class
    USE OMSAO_errstat_module
    IMPLICIT NONE
    INTEGER (KIND=4) :: OMI_SMF_setmsg !!this external function writes the PGS
    EXTERNAL            OMI_SMF_setmsg !!SMF error messages to the Log files
END MODULE OMI_SMF_class
