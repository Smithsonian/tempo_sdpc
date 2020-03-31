MODULE m_set_tracegas
  !INTEGER, PARAMETER :: ngas = 7, nallgas = 8
  !INTEGER            :: nfgas
  !INTEGER, DIMENSION(ngas), PARAMETER   :: gasidxs = (/no2_t1_idx,&
  !     so2_idx, so2v_idx,o2o2_idx,o2_idx, h2o_idx, bro_idx/)

 CONTAINS

 SUBROUTINE set_tracegas (ngas) 

 USE OMSAO_indices_module, ONLY: so2_idx, so2v_idx, o2o2_idx, o2_idx, h2o_idx, &
       bro_idx, no2_t1_idx, hcho_idx, o2t2_idx, h2ot2_idx
 USE OMSAO_variables_module, ONLY: npsl
 USE OMSAO_errstat_module, ONLY: www_lun
 USE ozprof_data_module, ONLY: use_effcrs,  gasname, gasidxs, &  ! input
       nfgas, fgasidxs, fgassidxs, fgaspos, fgasname, & ! input
       nhgas,hgaspos, &
       so2idx, so2vidx, o4idx, o2idx, o2t2idx,h2oidx, h2ot2idx,broidx,hchoidx, no2idx, & !gasidxs
       so2fidx, so2vfidx, o4fidx, o2fidx,o2t2fidx, h2ofidx, h2ot2fidx,brofidx, hchofidx, no2fidx,  & ![fgasidxs]
       so2sfidx,  so2vsfidx, o4sfidx, o2sfidx, h2osfidx, &
       so2crsidx, o2crsidx, o4crsidx, h2ocrsidx, & !allcol[]
       use_so2dtcrs, use_o4dtcrs, use_o2dptcrs, use_h2odptcrs, &
       do_so2shi, do_o4shi, do_o2shi, do_h2oshi, &
       do_so2tmp, do_o4tmp, do_o2tmp, do_h2otmp, &
       do_so2psl, do_o4psl, do_o2psl, do_h2opsl
        
 IMPLICIT NONE
 ! INPUT/OUTPUT variables
 INTEGER, INTENT(IN) :: ngas
 ! LOCAL variables
 INTEGER :: i, j
 LOGICAL, SAVE :: first = .true.

 WRITE(www_lun, *) '***** set tracegases ******'
 IF (first) THEN 
    DO i = 1, ngas
     IF (gasidxs(i) == so2_idx)    THEN
         so2fidx  = fgasidxs(i) ; so2idx = i 
     ENDIF
     IF (gasidxs(i) == so2v_idx)   THEN
         so2vfidx = fgasidxs(i) ; so2vidx = i
     ENDIF
     IF (gasidxs(i) == o2o2_idx)   THEN
         o4fidx   = fgasidxs(i) ; o4idx = i 
     ENDIF
     IF (gasidxs(i) == o2_idx)    THEN
         o2fidx   = fgasidxs(i) ; o2idx = i 
     ENDIF
     IF (gasidxs(i) == o2t2_idx)    THEN
         o2t2fidx   = fgasidxs(i) ; o2t2idx = i 
     ENDIF
     IF (gasidxs(i) == h2o_idx)   THEN
         h2ofidx   = fgasidxs(i) ; h2oidx = i 
     ENDIF
     IF (gasidxs(i) == h2ot2_idx)   THEN
         h2ot2fidx   = fgasidxs(i) ; h2ot2idx = i 
     ENDIF
     IF (gasidxs(i) == bro_idx) THEN
         brofidx = fgasidxs(i) ; broidx = i 
     ENDIF
     IF (gasidxs(I) == no2_t1_idx) THEN
         no2fidx = fgasidxs(i) ; no2idx = i
     ENDIF
     IF (gasidxs(i) == hcho_idx ) THEN
         hchofidx = fgasidxs(i) ; hchoidx = i
     ENDIF 
    ENDDO

    so2crsidx = 0 ; o4crsidx = 0 ; o2crsidx = 0 ; h2ocrsidx = 0
    IF (fgasidxs(so2vidx) <= 0 .AND. fgasidxs(so2idx) <= 0) use_so2dtcrs =.FALSE.
    IF (fgasidxs(o4idx) <= 0)  use_o4dtcrs   = .FALSE.
    IF (fgasidxs(o2idx) <= 0 ) use_o2dptcrs  = .FALSE.
    IF (fgasidxs(h2oidx) <= 0) use_h2odptcrs = .FALSE.

    nhgas = 0
    fgasname(1:nfgas) = gasname(fgaspos(1:nfgas))
    DO i = 1, nfgas
       j = fgaspos(i)
       IF (use_so2dtcrs .AND. (gasidxs(j) == so2_idx .OR. gasidxs(j) == so2v_idx)) so2crsidx = i
       IF (use_o4dtcrs  .AND. gasidxs(j) == o2o2_idx) o4crsidx = i
       IF (use_h2odptcrs.AND. gasidxs(j) == h2o_idx) h2ocrsidx = i
       IF (use_o2dptcrs .AND. gasidxs(j) == o2_idx) o2crsidx = i
       IF (j == h2oidx .or. j == o2idx) THEN 
          nhgas = nhgas + 1
          hgaspos (nhgas) = i
       ENDIF 
       WRITE(www_lun, *) fgasname(i), ':', i, j, gasidxs(j)
    ENDDO

    ! set for dads
    so2sfidx = 0; so2vsfidx = 0 ; o4sfidx = 0 ; o2sfidx = 0 ; h2osfidx = 0  
    do_so2shi = .false. ; do_o4shi = .false. ; do_o2shi = .false. ; do_h2oshi = .false.
    DO i = 1, nfgas
       j = fgaspos(i)
       ! find indices for shift
       IF (gasidxs(j) == so2_idx)  so2sfidx  = fgassidxs(j)
       IF (gasidxs(j) == so2v_idx) so2vsfidx = fgassidxs(j)
       IF (gasidxs(j) == o2o2_idx) o4sfidx   = fgassidxs(j)
       IF (gasidxs(j) == o2_idx)   o2sfidx   = fgassidxs(j)
       IF (gasidxs(j) == h2o_idx)  h2osfidx  = fgassidxs(j)
       IF ((gasidxs(j) == so2_idx .OR. gasidxs(j) == so2v_idx) .AND. &
           fgassidxs(j) > 0) do_so2shi = .TRUE.
       IF (gasidxs(j)  == o2o2_idx .AND.fgassidxs(j) > 0) do_o4shi = .TRUE.
       IF ((gasidxs(j) == o2_idx   .OR. gasidxs(j) == o2t2_idx) .AND. fgassidxs(j) > 0) do_o2shi = .TRUE.
       IF ((gasidxs(j) == h2o_idx  .OR. gasidxs(j) == h2ot2_idx).AND. fgassidxs(j) > 0) do_h2oshi = .TRUE.
!      WRITE(www_lun, *) i, fgasidxj(i), gasname(j)
    ENDDO 

    ! set for dadt
    do_so2tmp=.FALSE.;  do_o4tmp=.FALSE.; do_o2tmp=.FALSE. ;do_h2otmp=.FALSE.  
    
    ! set for dadp
    do_so2psl=.FALSE.;  do_o4psl=.FALSE.; do_o2psl=.FALSE. ;do_h2opsl=.FALSE.  
    IF (use_effcrs .and. npsl >0) THEN 
       !IF (fgasidxs(h2oidx) > 0) do_h2opsl = .true.
       !IF (fgasidxs(o2idx) > 0) do_o2psl = .true.
       !IF (fgasidxs(o4idx) > 0) do_o4psl = .true.
    ENDIF
    first = .FALSE.
 ENDIF 
 RETURN
 END SUBROUTINE

END MODULE
