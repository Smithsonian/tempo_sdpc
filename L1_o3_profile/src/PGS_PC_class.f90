MODULE PGS_PC_class
    IMPLICIT NONE
    include 'PGS_PC.f'
    include 'PGS_PC_9.f'

    INTEGER(KIND=4), EXTERNAL :: PGS_PC_getNumberOfFiles, &
                                 PGS_PC_getReference, &
                                 PGS_PC_getConfigData

END MODULE PGS_PC_class
