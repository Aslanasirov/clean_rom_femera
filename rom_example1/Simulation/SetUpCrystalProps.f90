!=====================================================================72
!
SUBROUTINE SetUpCrystalProps( )
!      
!
    implicit none
!---------------------------------------------------------------------72    
!
!    write(*,*) 'Begin SetUpCrystalProps'
!---------------------------------------------------------------------72    
!
    CALL CrystalOpenIOFiles( )
!--- Read in CoefTensors and assign them to the corresponding tensors
!
    CALL READ_COEFTENS( )
      
!---  Allocate array
    CALL ArrayAllocate( )    
!
!---  Read in data
!
    CALL CrystalModelData()
    
!--- Read in Orientation data
    CALL CrystalOrientationData()     
!      
!-- Initialize some of the arrays                  
!  
    CALL CrystalInitializeArrays( )
!
!------- parameters for inerations
!
    CALL CrystalSolverData( )
!	         
!-- Close opened files
!    
    return
END 
!=====================================================================72
!
!
!=====================================================================72
!
!
SUBROUTINE READ_COEFTENS()

    
    use NumType
    use workdir
    USE COEFTENS
    use fileIO
    ! ASLAN ADDED
    use NbData
    ! ASLAN ADDED

    implicit none
    
    INTEGER(KIND=IKind) :: ICTR, I, K1, K2, J, IPH1, I1, ILine
! aslan added
    INTEGER(KIND=IKind) ::  TotalLine, ip, nnz, NumNb, Nbtmp(Npart)
! aslan added
    INTEGER(KIND=IKind) :: lenoutdir, lenCoefFile, NPartctf, Reason
!
CHARACTER*255  filename
	REAL(KIND=RKind), ALLOCATABLE :: CoefTensVec(:)
!---------------------------------------------------------------------72
! aslan added
    !--- neighboring info input
    nnz = 0
    do ip=1, NPart 
        read(NbDataFile, *)  Nbtmp
        NumNb=minloc(Nbtmp,1)-1
        nnz=nnz+(NumNb+1)
        if (.not. allocated(Elset_nb(ip)%NBGID))  then
ALLOCATE( Elset_nb(ip)%NBGID(NumNb))
            Elset_nb(ip)%NBGID=Nbtmp(1:NumNb)
            Elset_nb(ip)%NNG=NumNb
        endif
    enddo
! aslan added

    MTDIM = 6
    if (.not. allocated(MT)) ALLOCATE ( MT( MTDIM ) )
    MT = (/ 1, 2, 3, 4, 5, 6 /)
!
!
    filename = adjustl(trim(outdir))//adjustl(trim(CoefTensFile))//'.dat'
!    
!-- echo the file path of coefTnes.dat     
    if (debug==1) write(FILE_E, '(A15, 1x, A55)') 'CoefTFile:', filename
    
  OPEN(UNIT=CoefTens_I, FILE=filename, STATUS='unknown', ACCESS='sequential')
  REWIND(CoefTens_I)    
!
READ (CoefTens_I, *) NPartctf
   Write(file_e,*) 'NPartctf=', NPartctf
!
IF (NPartctf .ne. NPart) call RunTimeError(                            &
         FILE_O, 'NPart in coeftenscompute different from module')        
!
!
    if (.not. allocated(VolFrac)) ALLOCATE( VolFrac(NPart) ) 
!-------Read the NumPartPerGrain and volume fraction from the first twoo lins of CoefTens.dat
!	
READ (CoefTens_I, *) VolFrac    

!-- Check sumamtion of VolFrac
    if (abs(sum(volFrac) - 1.d0) .gt. 1.0e-04) then
        Call RunTimeError (File_O,             &
            'Summation of volume fraction is not eqyual to 1.0')
    endif
!
!
    ! aslan changed
    ! ALLOCATE(CoefTensVec(2*NPart*36+nnz*36))
    ! TotalLine= CEILING((2*NPart*36+nnz*36)/8.d0)
    ! ASLAN SINCE COEFTENS IS COMPUTED USING INHOUSE CODE WE NEED TO SKIP G TENSOR ( Npart*6)
    ALLOCATE(CoefTensVec(2*NPart*36+nnz*36 + Npart*6))
    TotalLine= CEILING((2*NPart*36+nnz*36  + Npart*6)/8.d0)
    ! aslan added
!-------Read the CoefTens into CoefTensVec
ILine=1
    ! ASLAN COMMENTED OUT
! DO WHILE(.TRUE.)
    DO WHILE(ILine .le. (TotalLine-1))
!
READ (CoefTens_I, *, IOSTAT=Reason) CoefTensVec(((ILine-1)*8+1):(8*ILine))

IF (Reason<0) EXIT
ILine=ILine+1
ENDDO
    !----now read the last line to the reminder space of CoefTensVec
    ! ASLAN ADDED
    ! READ (CoefTens_I, *, IOSTAT=Reason) CoefTensVec((8*(ILine-1)+1):(2*NPart*36+nnz*36))
    READ (CoefTens_I, *, IOSTAT=Reason) CoefTensVec((8*(ILine-1)+1):(2*NPart*36+nnz*36 + Npart*6))
    ! ASLAN ADDED
!	
CLOSE (CoefTens_I)
!

if (.not. allocated(COEFM)) ALLOCATE( COEFM(NPart, nvec,nvec) )   
    if (.not. allocated(COEFA)) ALLOCATE( COEFA(NPart,nvec,nvec) )
    if (.not. allocated(COEFP)) ALLOCATE( COEFP(nnz,nvec,nvec) )
    !
!    
!--- Read COEFM    
!
    ICTR = 0
!    
!--- Read COEFA    
!    
    COEFA = 0.d0
    DO I = 1, NPart
        DO K1 = 1, MTDIM
            DO K2 = 1, MTDIM
                ICTR = ICTR + 1
                COEFA(I,MT(K2),MT(K1)) = CoefTensVec(ICTR)
            END DO
        END DO
    END DO
!
    ! ASLAN MIXED METHOD
    ! IF COEFTENS FILE IS CREATED USING INHOUSE CODE THEN WE NEED TO SKIP OVER THE SO CALLED G TENSOR (THERMAL?)
    !ICTR = ICTR + Npart * MTDIM


    COEFM = 0.d0
    DO I = 1, NPart
        DO K1 = 1, MTDIM
            DO K2 = 1, MTDIM
                ICTR = ICTR + 1
                COEFM(I, MT(K2),MT(K1)) = CoefTensVec(ICTR) 
            END DO
        END DO
    END DO
!--- Read COEFp    
!  
    COEFP = 0.d0
    DO I = 1, nnz
        ! aslan commented out
        ! aslan commented out
        DO K1 = 1, MTDIM
            DO K2 = 1, MTDIM
                ICTR = ICTR + 1
! Original ROM Code                    
!                    COEFP(I,J,MT(K2),MT(K1)) = PROPTS(ICTR)
!
! xiang or xiaoyu or yang -> My code, as my P^(ab) is actually P^(ba) of original ROM
! 
                    ! aslan commented out
                    ! aslan changed
                COEFP(I,MT(K2),MT(K1)) = CoefTensVec(ICTR)                     
            END DO
        END DO
            
    END DO
!   
!--echo CoefTensors for debuging
    if (debug==1) then
        write(FILE_E, '(A10, 1x, I5 )')                 &
                        'NPart from CoefTensCompute', NPartctf
        DO I=1, NPart
            write(FILE_E,*) 'A_ijkl{',I, '} is:'
            write(FILE_E,'(6(SP,ES15.6E3))') COEFA(i,:,:)
            write(FILE_E,*) 'M_ijkl{',i, '} is:'
            write(FILE_E,'(6(SP,ES15.6E3))') COEFM(i,:,:)
        ENDDO
!
        ! aslan commented out
        ! aslan added
        DO I =1, nnz
            write(file_e,*) 'P_ijkl{',I,'} is:'
            write(file_e,'(6(SP,ES15.6E3))') COEFP(I,:,:)
        ENDDO
    endif
!
    DEALLOCATE (CoefTensVec)!    
    RETURN
END SUBROUTINE  READ_COEFTENS
!
!=====================================================================72
!
SUBROUTINE ArrayAllocate()
!
    use NumType
    use CPVars
    use CPVars_n
    use SlipGeo
    use OriPar
    use PlaPar
    use fileIO
    
    implicit none
    integer iline, PartID, iPart, iSlip
!-----------------------------------------------------------------------
!
    if (maxnumslip .ne. max(PhSlip(1), PhSlip(1))) call RunTimeError(File_o,                      &
            'maxnumslip not equat to  max(PhSlip)')
    !write(*,*) 'maxnuslip= ', maxnumslip
!
!---First read in the phase ID for each part, which is needed for array allocation
!
    if (NPart .ne. sum(ParPerPh))  call RunTimeError(File_o,                      &
            'Numpart  is not equal  to sum(ParPerPh) ')     
!          
!
    Plap%PhID=0
!-- read the partid for phase 1    
    do iline=1,ParPerPh(1)
        read(Ph1PartFile, *) PartID
        if(Plap%PhID(PartID) .ne. 0) call RunTimeError(File_o,                      &
            ' Plap%PhID is not initialized to be zero! ')    
        Plap%PhID(PartID)=1 
    enddo       
    close(Ph1PartFile)
!-- read the partid for phase 1    
    do iline=1,ParPerPh(2)
        read(Ph2PartFile, *) PartID
        if(Plap%PhID(PartID) .ne. 0) call RunTimeError(File_o,                      &
            ' There is overlap between ph1.dat and ph2.dat ')          
        Plap%PhID(PartID)=2 
    enddo       
    close(Ph2PartFile)
    

    
!-----------------------------------------------------------------------

    RETURN
END 
!=====================================================================72
!
!
!=====================================================================72
!
SUBROUTINE CrystalOpenIOFiles(                                         &
                              )
!
    use NumType
    use WorkDir
    use FileIO
    implicit none
!
    character :: filename1*255, filename2*255, filename3*255, filename4*255, filename5*255,  filename6*255,  filename7*255
! aslan added
    character :: filename8*255, filename9*255, filename10*255, filename11*255
!
!---------------------------------------------------------------------72
!
!------- root name of input/output files
!
!------- open files
!
    
    filename1 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.xtali'

    open(unit=FILE_I, file=filename1, status='unknown',               &
                                                 access='sequential')
    rewind(FILE_I)

    filename2 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.xtale'
    
    open(unit=FILE_E, file=filename2, status='unknown')
    
    rewind(FILE_E)
    if (debug==0) write(File_E, *) 'debug=0, nothing written into this file!'

    filename3 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.xtalo'
    
    open(unit=FILE_O, file=filename3, status='unknown')
!
    filename4 = adjustl(trim(outdir))//adjustl(trim(textureFile))//'.txti'
    
    open(unit=TXT_I, file=filename4, status='unknown',          &
                                         access='sequential')
    rewind(TXT_I)
    !
    filename5 = adjustl(trim(outdir))//adjustl(trim(FilePre))//'.itero'
    
    open(unit=iter_O, file=filename5, status='unknown')
    
    rewind(iter_O) 
    if (iterprint==1) write(iter_O, '(3A12, A20)'), '----KStep----', '----KInc----', '----Iters----', '----errorinfo----'
    !
    filename6 = adjustl(trim(outdir))//adjustl(trim(PhasePart))//'1.dat'
    
    open(unit=Ph1PartFile, file=filename6, status='unknown')
    
    rewind(Ph1PartFile) 
!
    filename7 = adjustl(trim(outdir))//adjustl(trim(PhasePart))//'2.dat'
    
    open(unit=Ph2PartFile, file=filename7, status='unknown')
    
    rewind(Ph2PartFile)       


! ASLAN ADDED
    !
    filename8 = adjustl(trim(outdir))//'NBgrain.dat'
    
    open(unit=NbDataFile, file=filename8, status='unknown')
    
    rewind(NbDataFile) 

!--> Xiaoyu: open GNMap.dat file
    filename10 = adjustl(trim(outdir))//adjustl(trim(GNFile))//'.dat'
    open(unit=GN_I, file=filename10, status='unknown')
    rewind(GN_I)

    filename11 = adjustl(trim(outdir))//'stressstrain.dat'
    open(unit=stressfile, file=filename11, status='unknown')
    rewind(stressfile)

! ASLAN ADDED

!---echo file path
!
    if (debug==1) then 
        write(File_E, '(A15, 1x, A55)') 'xtalifile:', filename1
        write(File_E, '(A15, 1x, A55)') 'xtalefile:', filename2
        write(File_E, '(A15, 1x, A55)') 'xtalofile:', filename3
        write(File_E, '(A15, 1x, A55)') 'texturefile:', filename4
        write(File_E, '(A15, 1x, A55)') 'iterofile:', filename5        
    endif
    return
END
!
!=====================================================================72
!
!
!=====================================================================72
!
SUBROUTINE CrystalModelData(                                     &
                                            )
!         
    use numtype
    use FileIO
    USE WorkDir
    use PlaPar
    use PowerLawBoundPar
!      
    implicit none
    integer(kind=ikind) NPhin,NPartin, CrystalIDin(NPh), iPart, IPh
    real(kind=rkind)    h0in, xmin, gam0in, tausiin,      &
                        taus0in, xmsin, gamss0in, crss0in
    integer iline, PartID
!
!---------------------------------------------------------------------72
!
!
!---- open single crystal filename

    read(FILE_I, *) NPhin, NPartin
!---echo some data 
    if (debug==1)                                                      &
        write(FILE_E, '(A15, 1x, I8, 1x, A20, 1x, I8)')                &
                        'CrystalID:', NPhin,                       &
                        'NPart from xtali', NPartin    
!
    if (NPart .ne. NPartin)  call RunTimeError(File_o,                      &
            'Number of parts from test.xtali is not equal from that the module') 
!
    if (NPh .ne. NPhin)  call RunTimeError(File_o,                      &
            'Nuber of phases from test.xtali is not equal from that the module') 
!--- read PhID
!
    read(FILE_I, *)      CrystalIDin
!
! --assign the parameter for each phase                        
!
    do iPh=1, NPH  
!---read plastic parameters for each phase
        read(FILE_I, *)      h0in, xmin, gam0in, tausiin,      &
                          taus0in, xmsin, gamss0in, crss0in 
!                                     
        Plap%h0(iPh)=h0in
        Plap%xm(iPh)=xmin
        Plap%gam0(iPh)=gam0in
        Plap%tausi(iPh)=tausiin
        Plap%taus0(iPh)=taus0in
        Plap%xms(iPh)=xmsin
        Plap%gamss0(iPh)=gamss0in
        Plap%crss0(iPh)=crss0in                                                                                    
        Plap%kappa0(:, iPh)=Plap%crss0(iPh)
!        
        do ipart=1,NPart     
            !write(*,*) '    Ipart= ', iPart 
            if (Plap%PhID(iPart) .eq. iph) then
                Plap%CrystalID(iPart)=CrystalIDin(iPh)
                !write(*,*) '            Plap%CrystalID(iPart) ', CrystalIDin(iPh)
            endif
        enddo   
        if (debug==1) write(FILE_E,*) 'Plap%kappa0 of phase ', iPh,  'is assigned!'
    enddo
!  
!                         
!-- echo the plastic parameters
    if (debug==1) write(FILE_E, 1000)  Plap%CrystalID(1), NPartin,                & 
                         Plap%h0(1), Plap%xm(1), Plap%gam0(1), Plap%tausi(1),      &
                         Plap%taus0(1), Plap%xms(1), Plap%gamss0(1), Plap%crss0(1)
!
    do iph=1, NPh
        CALL BoundForArgPowLaw(Plap%xm(iPh), PLBP%argmin, PLBP%argmax)
    enddo
!---- vertices of rate independent yield surface (single crystal)
!
    do ipart=1,NPart
        if (Plap%crystalID(iPart) .eq. kFCC) then
            call SetSlipSystemFCC( iPart )
        else  ! crystalID=kBCC or kHCP 
            call RunTimeError( FILE_O,                                         &
                       'Error: Current code only supports FCC lattice!')
        endif
    enddo
!

1000  format(/'*-----   Crystal Plasticity parameters of part(1)-----*'/,        &
              7x,'  elasticity type  = ',i5/,                          &
              7x,'  number of parts  = ', i5/,                         &
              7x,'  h0               = ', e12.5/,                      & 
              7x,'  xm               = ', e12.5/,                      & 
              7x,'  gam0             = ', e12.5/,                      & 
              7x,'  tausi            = ', e12.5/,                      & 
              7x,'  taus0            = ', e12.5/,                      & 
              7x,'  xms              = ', e12.5/,                      & 
              7x,'  gamss0           = ', e12.5/,                      &               
              7x,'  crss0            = ', e12.5/ )
      return
      END
!
!=====================================================================72
!
SUBROUTINE CrystalOrientationData( )
!
    USE NumType
USE FILEIO
USE WorkDir
USE OriPar
USE SlipGeo
!
IMPLICIT none
!   
INTEGER          iPart, iSlip, i
	REAL(KIND=8)  :: pi180, EulerMax
	REAL(KIND=8)  :: sps, cps, sth, cth, sph, cph     

INTEGER    NumEulerAngle, iidr, iikc
CHARACTER  filename*255
!---------------------------------------------------------------------  
!
!
    filename = adjustl(trim(outdir))//adjustl(trim(textureFile))//'.txti'
!
    open(unit=TXT_I, file=filename, status='unknown', access='sequential')
  rewind(TXT_I)
!------- read the flag for angle convention and set some constants
!-------   iikc = 0 : angles input in Kocks convention :  (psi,the,phi)
!                 1 : angles input in Canova convention : (ph,th,om)
!                 ph = 90 + phi; th = the; om = 90 - psi
!-------   iidr = 0 : angles input in degrees
!-------          1 : angles input in radians

READ (TXT_I, *)  NumEulerAngle
IF (NumEulerAngle .NE. NPart)  CALL RunTimeError(FILE_O,             &
                       'AssignCrystalODF: NumEulerAngle != NumPart')
!
READ(TXT_I, *) iikc, iidr
! ------- Euler angles for each partitition
pi180 = 4.0 * datan(1.0d+00)/180.
IF (iidr .eq. 1) pi180=1.0         
OriP%Euler=pzero  
EulerMax=pzero     
!     
!---
    IF (debug==1) WRITE(FILE_E,*) 'Euler Angles are:'       
DO iPart=1,NPart
READ(TXT_I, *) OriP%Euler(:,iPart)
OriP%Euler(:,iPart)=OriP%Euler(:,iPart)*pi180

if (debug ==1) WRITE(FILE_E,'(3F12.5)')  OriP%Euler(:,iPart)
DO i=1,3
IF (EulerMax .LT. abs(OriP%Euler(i,iPart)))                &
                EulerMax=abs(OriP%Euler(i,iPart)) 
END DO
ENDDO

            
IF (EulerMax .GT. 10)  &
CALL RunTimeError(FILE_o,                                          &
   'Check the units of Euler angles to make sure they are in Rad') 
    
CLOSE (TXT_I)
 

!     
!! ------- build rotation matrices C0: {x}_sm = [C0] {x}_cr 
!!      
DO iPart=1,NPart
CALL AnglesToRotMatrix(OriP%Euler(:,iPart),                    &
                        OriP%gcrot0(:,:,iPart), NSD)
!---Rotate the slip system into global coordinate system
        DO ISlip=1,maxNumSlip            
SlipG%VecM0(:,ISlip, IPart)=MATMUL(OriP%gcrot0(:,:,IPart),SlipG%VecM0(:,ISlip, IPart))
SlipG%VecS0(:,ISlip, IPart)=MATMUL(OriP%gcrot0(:,:,IPart),SlipG%VecS0(:,ISlip, IPart))
!
            CALL OuterProductVec(SlipG%VecS0(:,ISlip, IPart), SlipG%VecM0(:,ISlip, IPart), &
                                                    SlipG%ZTen0(:,:,ISlip,IPart), NSD)                              
            CALL ZTenToVec(SlipG%ZTen0(:,:,ISlip,IPart), SlipG%ZVec0(:, ISlip,IPart),NSD)
        ENDDO                        
!                        
if (debug==1 .and. iPart==1) then
    WRITE(FILE_E,'(A15, 1x, I5/, 3(5x, F12.5))')     &
 'gcrot0-part-', IPart, OriP%gcrot0(:,:,iPart)
            write(FILE_E,*)   '-----Slip system in global coordinates----- ' 
            write(FILE_E,*)   'vecS0 of part 1 is:' 
            write(FILE_E,'(3F10.4)')  SlipG%vecS0(:,:,1)
            write(FILE_E,*)   'vecM0 of part 1 is:' 
            write(FILE_E,'(3F10.4)')  SlipG%vecM0(:,:,1)   
            write(FILE_E,*)   'ZVec0 of part 1 is:' 
            write(FILE_E,'(6E12.5)')  SlipG%ZVec0(:,:,1)             
        endif
ENDDO 


RETURN
END

!
!=====================================================================72
!
!=====================================================================72
!
SUBROUTINE SetSlipSystemFCC(   PartNum  )
!
    use NumType
    use SlipGeo 
    use FileIO 
    use PlaPar
!    
   IMPLICIT NONE
! 
	REAL(KIND=8) ::  sDOtm
	REAL(KIND=8) :: indexM(3,12), indexS(3,12)
	REAL(KIND=8) :: InnerProductVec
INTEGER(kind=iKind)      ::  ISlip, PartNum
      !data indexM /1.,  1., -1.,   &
      !             1.,  1., -1.,   &
      !             1.,  1., -1.,   &
      !             1., -1., -1.,   &
      !             1., -1., -1.,   &
      !             1., -1., -1.,   &
      !             1., -1.,  1.,   &
      !             1., -1.,  1.,   &
      !             1., -1.,  1.,   &
      !             1.,  1.,  1.,   &
      !             1.,  1.,  1.,   &
      !             1.,  1.,  1./
      !data indexS /0.,  1.,  1.,   &
      !             1.,  0.,  1.,   &
      !             1., -1.,  0.,   &
      !             0.,  1., -1.,   &
      !             1.,  0.,  1.,   &
      !             1.,  1.,  0.,   &
      !             0.,  1.,  1.,   &
      !             1.,  0., -1.,   &
      !             1.,  1.,  0.,   &
      !             0.,  1., -1.,   &
      !             1.,  0., -1.,   &
      !             1., -1.,  0./   
      
      
!!----Same as Marin's Order      
       data indexM /1.,  1.,  1.,   &
                    1.,  1.,  1.,   &
                    1.,  1.,  1.,   &
                                    !&
                   -1.,  1.,  1.,   &
                   -1.,  1.,  1.,   &
                   -1.,  1.,  1.,  &
                                    !&
                    -1., -1.,  1.,  &
                    -1., -1.,  1.,  &
                    -1., -1.,  1.,  &
                                    !&
                    1.,  -1., 1.,   &
                    1.,  -1., 1.,   &
                    1.,  -1., 1./
       
       data indexS /0.,  1.,  -1.,  &
                    1.,  0.,  -1.,  &
                    1., -1.,  0.,   &
                                    !&
                    0.,  1., -1.,   &
                    1.,  0.,  1.,   &
                    1.,  1.,  0.,   &
                                    !&
                    0.,  1.,  1.,   &
                    1.,  0.,  1.,   &
                    1.,  -1.,  0.,  &
                                    !&
                    0.,  1., 1.,    &
                    1.,  0., -1.,   &
                    1., 1.,  0./ 

!!----Same as that in my slides   
!      data indexM /1.,  1.,  1.,   &
!                   1.,  1.,  1.,   &
!                   1.,  1.,  1.,   &
!                                   !&
!                  -1.,  1.,  1.,   &
!                  -1.,  1.,  1.,   &
!                  -1.,  1.,  1.,  &
!                                   !&
!                   1., -1.,  1.,  &
!                   1., -1.,  1.,  &
!                   1., -1.,  1.,  &
!                                   !&
!                   1.,  1., -1.,   &
!                   1.,  1., -1.,   &
!                   1.,  1., -1./
      
!      data indexS /0.,  -1.,  1.,  &
!                   1.,  0.,  -1.,  &
!                   -1., 1.,  0.,   &
!                                   !&
!                   1.,  0., 1.,   &
!                   1.,  1.,  0.,   &
!                   0.,  -1.,  1.,   &
!                   0.,  1.,  1.,   &
!                   1.,  1.,  0.,   &
!                   1.,  0.,  -1.,  &
!                   0.,  1., 1.,    &
!                   1.,  0., 1.,   &
!                   -1., 1.,  0./ 
      
!---------------------------------------------------------------------72
!
!----echo the slip systems
    if (debug==1) then
        write(FILE_E,*)  'indexS is:' 
        write(FILE_E,'(3F10.4)')  indexS
        write(FILE_E,*)  'indexM is:' 
        write(FILE_E,'(3F10.4)')  indexM     
    endif
!
!------- slip system normals and slip directions: unit vectors
!------- for FCC, numslip=12 is by default, make changes IF necessary.

        if (debug==1) write(FILE_E,*) 'Plap%kappa0(:,IPart) is assigned!' 
        Plap%kappa0(:,Plap%PhID(PartNum))=Plap%crss0(Plap%PhID(PartNum))
        DO ISlip = 1, PhSlip(Plap%PhID(PartNum))
!		
SlipG%VecS0(:, ISlip, PartNum)=indexS(:,ISlip)/DSQRT         &
            (DOT_PRODUCT(indexS(:,ISlip), indexS(:,ISlip)))
SlipG%VecM0(:, ISlip, PartNum)=indexM(:,ISlip)/DSQRT         &
            (DOT_PRODUCT(indexM(:,ISlip), indexM(:,ISlip)))
!
!------- check normality of vecS and vecM
sDOtm = DOT_PRODUCT(SlipG%vecS0(:,ISlip,PartNum),            &
                                 SlipG%vecM0(:,ISlip,PartNum))
!			                                       
IF (sDOtm .ne. 0)  CALL RunTimeError(FILE_O,  &
             'Normality of vecS and vecM is not satisfied!')
!		    			
!
            CALL OuterProductVec(SlipG%VecS0(:,iSlip,PartNum),           &
                                 SlipG%VecM0(:,iSlip,PartNum),           &
                                 SlipG%zTen0(:,:,iSlip,PartNum), NSD)
            CALL ZTenToVec(SlipG%ZTen0(:,:,iSlip,PartNum),               &
                                        SlipG%ZVec0(:,iSlip,PartNum),NSD)
ENDDO 
       
!----echo the normilized slip system
      IF (debug==1) then
         write(FILE_E,*)   '-----Slip system in FCC----- ' 
         write(FILE_E,*)   'vecS0 of part 1 is:' 
         write(FILE_E,'(3F10.4)')  SlipG%vecS0(:,:,1)
         write(FILE_E,*)   'vecM0 of part 1 is:' 
         write(FILE_E,'(3F10.4)')  SlipG%vecM0(:,:,1)
     ENDIF
!
!
!
!              
    RETURN
END
!
!=====================================================================72
!
!=====================================================================72
!
SUBROUTINE CrystalSolverData(                                          &
         )
!
    use FILEIO
    use IterPar
    implicit none
!
!
!---------------------------------------------------------------------72
!
!------- number iterations and tolerance for state iterations
!
    read(FILE_I, *) iterP%maxIterState, iterP%tolerState
!
!------- number iterations and tolerance for newton method
!
    read(FILE_I, *) iterP%maxIterNewt, iterP%tolerNewt
!
!------- echo input data
    IF (debug==1)                                                      &
         write(FILE_E, 1000) iterP%maxIterState, iterP%maxIterNewt,    &
                          iterP%tolerState, iterP%tolerNewt
!
!------- format
!
1000  format(/'*-----   Local Convergence Control Parameters  -----*'/,&
              7x,'  max iters State  = ',i5 /,                         &
              7x,'  max iters Newton = ',i5 /,                         &
              7x,'  tolerance State  = ',e12.5 /,                      &
              7x,'  tolerance Newton = ',e12.5)
    return
END
!
!=====================================================================72
!
!
!=====================================================================72
!
      SUBROUTINE CrystalInitializeArrays(                              &
         )         
!         
      use numtype
      use CPVars
      use CPVars_n
      use oriPar
      use SlipGeo
      use PlaPar
!      
      implicit none
      integer(kind=ikind) iPart
!
!---------------------------------------------------------------------72
!
!------- initialize arrays used in XTAL constitutive integration method
!
        CPV%     gfip     =   pzero
        CPV%     gsumofgamdot     =   pzero
        CPV_n%     gfip_n     =   pzero
        CPV_n%     gsumofgamdot_n     =   pzero
!
        CPV%     gstress     =   pzero
        CPV%     gestran     =   pzero
        CPV%     gkappa      =   pzero
        CPV_n%   gstress_n   =   pzero
        CPV_n%   gestran_n   =   pzero
        CPV_n%   gkappa_n    =   pzero
        CPV%     gmu         =   pzero
        CPV_n%   gmu_n       =   pzero
!        
        CPV%     ggamdot     =   pzero
!
!
    do iPart = 1, NPart
        CPV%gkappa(:, iPart) =Plap%kappa0(:, Plap%PhID(iPart))
        CPV_n%gkappa_n(:, iPart) =Plap%kappa0(:, Plap%PhID(iPart))
    enddo
!
      return
      END
!
!=====================================================================72
!
!
!=====================================================================72
!
SUBROUTINE CrystalCloseIOFiles( )
!
    use FileIO
    implicit none
!
!---------------------------------------------------------------------72
!
!------- close files
!
    close(FILE_I)
    close(FILE_E)    
    close(FILE_O)
!
    return
END
!
!=====================================================================72
