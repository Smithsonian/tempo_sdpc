module MetadataModule
!!-----------------------------------------------------------------
!!-----------------------------------------------------------------
 IMPLICIT NONE
 
 PUBLIC :: RdWrmetadata

 CONTAINS

   FUNCTION RdWrmetadata(outfilnm) RESULT(status)
     !-------------------------------------------------------------------
     ! This function reads the metadata from L1B file
     !-------------------------------------------------------------------

     USE MetadatOBJModule
     USE m_LUN_set
     USE m_vars, ONLY: cloud_pres, eff_cld_frac, n_good_input, n_good_output, &
          highqual, badqual, cloud_pres_max, cld_frac_min, qc,   &
          n_input, n_missing, using_cal
     USE m_swathnames, ONLY: vis, visz
     USE m_pgs_include

     IMPLICIT NONE

     ! INCLUDE 'PGS_SMF.f'
     ! INCLUDE 'PGS_MET_13.f'
     ! INCLUDE 'PGS_OMI_1900.f'
     ! INCLUDE 'PGS_OMCLDRR_52251.f'

     INTEGER, PARAMETER :: INVENTORY=2
     INTEGER, PARAMETER :: ARCHIVE=3 

     INTEGER :: omi_smf_setmsg
     INTEGER :: pgs_MET_getPCAttr_I, pgs_MET_getPCAttr_d, pgs_MET_getPCAttr_s
     INTEGER :: pgs_MET_setAttr_I, pgs_MET_setAttr_d, pgs_MET_setAttr_s,          &
          pgs_MET_setmultiAttr_s
     INTEGER :: pgs_met_init,pgs_met_write, pgs_pc_getreference, pgs_met_sfstart, &
          pgs_pc_getconfigdata, & !omi_localgranuleid, &
          pgs_met_sfend, pgs_met_remove
     ! INTEGER :: pgs_pc_getuniversalref

     INTEGER :: i, status, returnstatus, version, sdid, Fil_Lun, ierr, j, ind, resid_id 
     integer, parameter :: nadd = 14, ninp = 9 

     INTEGER :: OrbitNumber, OrbitNumber_PCF, Qamissingdata , Qaboundsdata, VersionID, &
          QAPercentCloudCover, PerGoodQualData, ThreshOrbitNumber

     REAL(KIND=8) :: EqCrossLon
     ! REAL(KIND=8) :: DeEqCrLon, AsEqCrLon

     CHARACTER(LEN=PGSd_MET_GROUP_NAME_L) :: GROUPS(PGSd_MET_NUM_OF_GROUPS)
     CHARACTER(LEN=100), DIMENSION(50) :: Objvalue  
     CHARACTER(LEN=100), DIMENSION(NADD) :: AddAttrNam, AddAttrVal                                     
     CHARACTER(LEN=100), DIMENSION(ninp) :: InputPnt,supflnm 
     ! CHARACTER(LEN=100), DIMENSION(2) :: LcInputID 
     ! CHARACTER(LEN=100) :: LocalInGrID, LocalGrID, DesCrRevision, OperationMode
     CHARACTER(LEN=100) :: value, L1B_AutQualFl
     CHARACTER(LEN=200) :: buf
     CHARACTER(LEN=350) :: expl="Flag set to Passed if QAPercentHighQualityData >= 80%, "// &
          "Flag set to Suspent if QAPercentHighQualityData >= 20%, "//& 
          "or L1B AutomaticQualityFlag not set to Passed, "//         &
          "otherwise Flag set to Failed" 
     character (len=*), intent(in) :: outfilnm
     !**********************************************************************

     !
     ! start reading metadata from l1b file
     !
     status = OMI_S_SUCCESS

     DO i=1,ninvname
       version = 1
       returnstatus = pgs_met_getPCAttr_s(L1B_LUN, version , "CoreMetadata.0", &
            trim(INVOBJ(i)),Objvalue(i))
       IF(returnstatus /= 0 ) THEN
         status = OMI_E_GENERAL
         GO TO 9999               
       ENDIF
     ENDDO

     returnstatus = pgs_met_getPCAttr_i(L1B_LUN, version , "CoreMetadata.0", &
          "AUTOMATICQUALITYFLAG.1",L1B_AutQualFl)
     IF(returnstatus /= 0 ) THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

     !preliminary estimates
     !ind = count(cloud_pres > cloud_pres_max .or. cloud_pres < 0.0)
     ind = count(btest(qc(:,:),2) .or. btest(qc(:,:),3))
     QAboundsdata = nint( real(ind) / real(size(cloud_pres))*100.0)

     ind = count(eff_cld_frac > cld_frac_min)
     QAPercentCloudCover= nint( real(ind) / real(size(eff_cld_frac))*100.0)

     PerGoodQualData = nint( real(n_good_output)*100.0 / real(n_good_input))
     if(PerGoodQualData >= highqual ) then
       value = "Passed"
     else if(PerGoodQualData >= badqual) then
       value = "Suspect"
     else
       value = "Failed"
     endif
     if( trim(value) == "Passed" .and. trim(L1B_AutQualFl) /= "Passed") value = "Suspect"

     QAmissingdata = nint( real(n_missing)*100.0 /2.0/ real(n_input)) ! n_missing counts twice in m_cloud_pres_ret

     returnstatus = pgs_met_getPCAttr_i(L1B_LUN, version , "CoreMetadata.0", &
          "OrbitNumber.1",OrbitNumber)
     IF(returnstatus /= 0 ) THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

     returnstatus = pgs_pc_getconfigdata(OrbNum_LUN,buf)
     read(buf,*) OrbitNumber_PCF
     if(OrbitNumber /= OrbitNumber_PCF) then
       print *,"RdWrmetadata: OrbitNumbers do not match"
       !    ierr = OMI_SMF_setmsg(OMCLDRR_W_MET, &
       !      "OrbitNumbers do not match", "MetadaModule", 0 )
     endif

     returnstatus = pgs_pc_getconfigdata(ThreshOrbNum_LUN,buf)
     read(buf,*) ThreshOrbitNumber

     !  returnstatus = pgs_pc_getconfigdata(OperationMode_LUN,OperationMode)

     returnstatus = pgs_met_getPCAttr_d(L1B_LUN, version, "CoreMetadata.0", &
          "EQUATORCROSSINGLONGITUDE.1", EqCrossLon)

     IF(returnstatus /= 0 ) THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

     returnstatus = pgs_met_getPCAttr_i(L1B_LUN, version , "CoreMetadata.0", &
          "VERSIONID",VersionID)
     IF(returnstatus /= 0 ) THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

     !get L1B PSAs (see OMI Guidelines for Migration L1B Metadata to L2 Output)
     do i=1,nadd
       j=i
       if(i>4) j=i+1
       if(vis .or. visz) then
         if(i>11) j=i+7
       else
         if(i>11) j=i+13
       endif
       write(buf,*) j
       buf=adjustl(buf) 
       returnstatus = pgs_met_getPCAttr_i(L1B_LUN, version , "CoreMetadata.0", &
            "ADDITIONALATTRIBUTENAME."//trim(buf),AddAttrNam(i))
       if(returnstatus /= 0) AddAttrNam(i) = 'Missing'
       returnstatus = pgs_met_getPCAttr_i(L1B_LUN, version , "CoreMetadata.0", &
            "PARAMETERVALUE."//trim(buf),AddAttrVal(i))
       IF(returnstatus /= 0 ) THEN
         AddAttrVal(i) = '0'
         !   status = OMI_E_GENERAL
         !   GO TO 9999
       ENDIF
     enddo

     ! Reading Archive metadata 

     !  returnstatus = pgs_met_getPCAttr_s(L1B_LUN, version , "ArchiveMetadata.0", &
     !                               "LOCALINPUTGRANULEID",LocalInGrID)
     !  IF(returnstatus /= 0 ) THEN
     !  status = OMI_E_GENERAL
     !  GO TO 9999
     !  ENDIF  

     !  returnstatus = pgs_met_getPCAttr_d(L1B_LUN, version , "ArchiveMetadata.0",&
     !                               "DESCENDINGEQUATORCROSSINGLONGITUDE",DeEqCrLon)
     !  IF(returnstatus /= 0 ) THEN
     !  status = OMI_E_GENERAL
     !  GO TO 9999
     !  ENDIF   

     !  returnstatus = pgs_met_getPCAttr_d(L1B_LUN, version , "ArchiveMetadata.0",&
     !                               "ASCENDINGEQUATORCROSSINGLONGITUDE",AsEqCrLon)
     !  IF(returnstatus /= 0 ) THEN
     !  status = OMI_E_GENERAL
     !  GO TO 9999
     !  ENDIF

     !  returnstatus = pgs_met_getPCAttr_s(L1B_LUN, version , "ArchiveMetadata.0",&
     !                               "DESCRREVISION",DesCrRevision)
     !  IF(returnstatus /= 0 ) THEN
     !  status = OMI_E_GENERAL
     !  GO TO 9999
     !  ENDIF

     ! Set metadata for L2

     if(OrbitNumber .le. ThreshOrbitNumber) then
       resid_id=resid_id_early
     else
       resid_id=resid_id_late
     endif

     DO i=1,ninp
       IF(i==1)Fil_Lun=L1B_LUN
       IF(i==2)Fil_Lun=IRR1B_file
       IF(i==3)Fil_Lun=terr_prs_id
       IF(i==4)Fil_Lun=chl_id
       IF(i==5)Fil_Lun=oc_ram_id
       IF(i==6)Fil_Lun=ring_id
       IF(i==7)Fil_Lun=thresh_id
       IF(i==8)Fil_Lun=resid_id
       IF(i==9)Fil_Lun=refl_id
       IF(i==10) then 
         if(using_cal) then 
           Fil_Lun=cal_id
         else
           supflnm(i)=''
           go to 50
         endif
       ENDIF
       version = 1
       !!  if( i > 5) then
       !!   Fil_Lun=ring_id
       !!   version = i-5
       !!  endif

       returnstatus = PGS_PC_GetReference( Fil_Lun, version, buf )
       IF( returnstatus /= 0 ) THEN
         WRITE( buf,'(A,I0)' ) "get filename failed at LUN = ", Fil_Lun 
         !        ierr = OMI_SMF_setmsg( OMCLDRR_F_FAILURE, buf, "MetadataModule", 0 )
         GO TO 9999
       ELSE
         j = INDEX( buf, '/', BACK = .TRUE. ) + 1
         supflnm(i) = TRIM( buf( j:) )
       ENDIF
       !!  returnstatus = pgs_pc_getuniversalref(Fil_Lun,version,supflnm(i))
50     continue
     ENDDO

     InputPnt=supflnm

     returnstatus = pgs_met_init(MCF_LUN, GROUPS)
     IF(returnstatus /= 0 ) THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

     DO i=1,ninvname
       returnstatus = pgs_met_setattr_s(GROUPS(INVENTORY),trim(INVOBJ(i)),Objvalue(i))
       IF(returnstatus /= 0 ) THEN
         status = OMI_E_GENERAL
         GO TO 9999
       ENDIF
     ENDDO

     do i=1,nadd
       write(buf,*) i
       buf=adjustl(buf) 
       returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),  &
            "ADDITIONALATTRIBUTENAME."//trim(buf),AddAttrNam(i))
       returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),  &
            "PARAMETERVALUE."//trim(buf),AddAttrVal(i))
       IF(returnstatus /= 0 ) THEN
         status = OMI_E_GENERAL
         GO TO 9999
       ENDIF
     enddo

     returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),"QAPERCENTMISSINGDATA.1",Qamissingdata)
     returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),"QAPERCENTCLOUDCOVER.1",QAPercentCloudCover)
     returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),"QAPERCENTOUTOFBOUNDSDATA.1",Qaboundsdata)        
     returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),"AUTOMATICQUALITYFLAG.1",value)        
     returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),"AUTOMATICQUALITYFLAGEXPLANATION.1",expl)        

     returnstatus = pgs_met_setattr_i(GROUPS(INVENTORY),"OrbitNumber.1", OrbitNumber)
     returnstatus = pgs_met_setattr_d(GROUPS(INVENTORY),"EQUATORCROSSINGLONGITUDE.1", EqCrossLon)         

     !  returnstatus = pgs_met_setattr_s(GROUPS(INVENTORY),"InputPointer",InputPnt)
     returnstatus = pgs_MET_setmultiAttr_s(GROUPS(INVENTORY),"InputPointer",ninp,InputPnt)
     !  returnstatus = pgs_met_setattr_s(GROUPS(INVENTORY),"OperationMode.1",OperationMode)
     ! temporary ParameterName 
     returnstatus = pgs_met_setattr_s(GROUPS(INVENTORY),"ParameterName.1", "Cloud_Pressure")

     IF(returnstatus /=0)THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

     ! set Attr for Archive metadata 

     !  returnstatus = pgs_met_setattr_s(GROUPS(ARCHIVE),"LOCALINPUTGRANULEID",LocalInGrID)
     !  returnstatus = pgs_met_setattr_d(GROUPS(ARCHIVE),"DESCENDINGEQUATORCROSSINGLONGITUDE",DeEqCrLon)
     !  returnstatus = pgs_met_setattr_d(GROUPS(ARCHIVE),"ASCENDINGEQUATORCROSSINGLONGITUDE",AsEqCrLon)

     version =1

     returnstatus = pgs_met_sfstart( trim(outfilnm), HDF5_ACC_RDWR,sdid)

     IF(returnstatus /=0) THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

!!!Temporarily disabled since I have no code that includes it
     !  returnstatus = omi_localgranuleid( 2, groups)

     returnstatus = pgs_met_write(groups(INVENTORY),'CoreMetadata',sdid)  
     returnstatus = pgs_met_write(groups(ARCHIVE),'ArchiveMetadata',sdid)

     IF(returnstatus /=0) THEN
       print *,"write ArchiveMetadata failed"
       ierr = OMI_SMF_setmsg( OMCLDRR_W_MET, "write ArchiveMetadata failed", &
            "MetadataModule", 0 )
     ENDIF

     returnstatus = pgs_met_sfend(sdid)
     returnstatus = pgs_met_remove()

     IF(returnstatus /=0) THEN
       status = OMI_E_GENERAL
       GO TO 9999
     ENDIF

9999 CONTINUE 

     IF (status /= OMI_S_SUCCESS) THEN
       print *,'RdWrmetadata: failed to read or write metadata'
       status = omi_smf_setmsg(OMI_E_GENERAL, &
            'failed to read or write metadata', &
            'RdWrmetadata', 1)
     ELSE
       status = omi_smf_setmsg(OMI_S_SUCCESS, &
            'Metadata part Successfull', 'RdWrmetadata',1)
       status = OMI_S_SUCCESS
     END IF

     RETURN
   END FUNCTION RdWrmetadata
END module MetadataModule
