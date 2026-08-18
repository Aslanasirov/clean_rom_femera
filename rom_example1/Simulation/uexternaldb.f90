INCLUDE 'SetUpCrystalProps.f90'
INCLUDE 'mkl_pardiso.f90'
!	
SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)

    use numtype
    use workdir
    use timing
    implicit none


    INTEGER :: LOP,LRESTART,KSTEP,KINC, length, lenjobname, lenoutdir
    REAL(KIND=rkind)::TIME(2),DTIME

!---------------Timing--------------
    CHARACTER(8) :: DATE
    CHARACTER(10) :: CURRENTTIME
    CHARACTER(5) :: ZONE
    INTEGER,DIMENSION(8) :: VALUES
!-----------------------------------
    logical :: deldexists,comp
    integer :: STAT
    ! ASLAN PRINTING
    ! ALSAN PRINTING
    IF (LOP == 0 .or. LOP == 4) THEN !BEGINNING ANALYSIS
!
        FirstIncr=.true.
!
!---    initialize the time counter of the dense solver and sparse solver        
!
        tss=0.d0
        ths=0.d0
        tsj=0.d0        
!  
        write(*,*) 'uexternaldb called, LOP=',LOP
        call system_clock(count_rate=cr)
        CALL system_clock(count_max=cm)  
        call system_clock(c1) 
        
        ! obtain the simulation director
        call getoutdir( outdir,lenoutdir )
        length = index(outdir,' ') - 1
        outdir = outdir(1:length)//'/'
        
!
        call getjobname( jobname,lenjobname )   
        CALL DATE_AND_TIME(DATE,CURRENTTIME,ZONE,VALUES)
        write(*,"(A30,A,A,I2,A,I2,A,I2,A,I3,A)")&
             'Start Analysis of job ', TRIM(ADJUSTL(jobname)), ': ',  VALUES(5),&
            ':', VALUES(6), ':', VALUES(7), ' on', VALUES(3), '^th'
        
        if (debug==1) write(*,*) 'setting up crystal props'    
        call SetUpCrystalProps(    )            
!
!   ASLAN ADDED
!--> Xiaoyu: initialize pardiso solver
        if (debug==1) write(*,*) 'calling pardiso ini'    
        call PardisoIni()
!   ASLAN ADDED
    ELSE IF (LOP == 3) THEN !END OF ANALYSIS
!
    
        call system_clock(c2)    
!    
        write(*,*) 'tss= ', tss, 'ths= ', ths, 'tsj= ', tsj
        CALL DATE_AND_TIME(DATE,CURRENTTIME,ZONE,VALUES)
        write(*,"(A30,A,A,I2,A,I2,A,I2,A,I3,A)")&
             'End Analysis of job ', TRIM(ADJUSTL(jobname)), ': ',  VALUES(5),&
            ':', VALUES(6), ':', VALUES(7), ' on', VALUES(3), '^th'  
            Write(*,'(A30, I6, A, I6, A, f12.5, A)') 'Total walltime: ', floor((c2-c1)/REAL(cr)/3600),  &
                        ' hours', floor((c2-c1)/REAL(cr)/60), ' minutes',                       &
                    ((c2-c1)/REAL(cr)-floor((c2-c1)/REAL(cr)/3600)*3600- floor((c2-c1)/REAL(cr)/60)*60), ' seconds.'
                    
    ELSE                    
        FirstIncr=.false.
         
!
    ENDIF
 
    
    RETURN

END SUBROUTINE UEXTERNALDB




!--------------
SUBROUTINE FINALIZE()

    implicit none
!

    RETURN
    
END SUBROUTINE FINALIZE

! ASLAN ADDED
!--------------
SUBROUTINE PardisoIni()
    !
    !---This subroutine is to make up an matrix that has the same structure (dimension, sparsity pattern) 
    !------as the system jacobian and conduct the renumbering and symbolic factorization for pardiso
    !
        use FileIO
        use PardisoROMVar
        use PardisoDDSDDEVar
        use NumType
        use COEFTENS
        use PlaPar
        use SlipGeo
        use mkl_pardiso

        implicit none

    !  
        integer iii, jjj, kkk, npgn, ng, b1, b2, b3, nnzAccumulate, nnzCols, iiiDiag
        integer iPart, ip1, ip2
    !
        real(kind=8) lhsPart0(NVEC, NVEC), dtime0
    
        real(kind=8) dmudotdsigma0(NVEC, NVEC), I_P(NVEC,NVEC), Term1(NVEC,NVEC)
        real(kind=8) delta, phiab(NVEC, NVEC), etaab(NVEC, NVEC)
    !
    !---- Assume value for dmudotdsigma0 and dtime0 to reproduce the sparsity patter of the jabobian
        dmudotdsigma0=pone
        dtime0=pone
        if (debug==1) write(*,*) 'starting pardisoini'    
    !
    !-- start to read the grain neigobor information and calcuate the number of non-zero conponents in the jacobian
        !--- Please note that ROM linear syatem and the linear system for solving the consistent moduli share the same sparsity pattern 
        !
        read(GN_I, *) npgn
        if (npgn .ne. NPart) call RunTimeError(FILE_O, 'Npart from GNMap.dat not correct!')
        
        allocate(Krow(NPart), Krow_ddsdde(NPart))
    !   
        nnz=0 
        do iPart=1,Npart
            read(GN_I, *) ng
            nnz=nnz+ng*NVEC*NVEC
            allocate(Krow(iPart)%col(ng), Krow_ddsdde(iPart)%col(ng))
            read(GN_I, *) Krow(iPart)%col
        enddo
        Krow_ddsdde=Krow
    !
        allocate(ja(nnz), ia(NVEC*NPart+1))
        allocate(amatrix(nnz)) 
        ia(NVEC*NPart+1)=nnz+1      
    !
    !---Starts to obtain the compressed storage format for pardiso
        ! please refer to pardiso5.0 manual http://pardiso-project.org/manual/manual.pdf for more details
    !
        kkk=0              
        nnzAccumulate = 0
        do ip1=1,NPart      !ip1 is corresponding to the index beta in the equation
            nnzCols = size(krow(ip1)%col)
            do iii=1, size(krow(ip1)%col) ! ip2 is corresponding to the index eta in the equation
                ip2=krow(ip1)%col(iii)
                delta=pzero
                if (ip2 .eq. ip1) delta=pone 
                I_P=delta*Ident4th-COEFP(nnzAccumulate+iii, :, :)
                call MultAxB(I_P,dmudotdsigma0,Term1, 6)
                lhsPart0=delta*COEFM(ip1, :,:)/dtime0+Term1
                do jjj=1,NVEC
                    b1= 36*kkk+1+6*(iii-1)+6*(jjj-1)*size(krow(ip1)%col)
                    b2= 36*kkk+6+6*(iii-1)+6*(jjj-1)*size(krow(ip1)%col)
                    Amatrix(b1: b2)=lhsPart0(jjj,:)
                    ja(b1 : b2)=(/ 6*(ip2-1)+1, 6*(ip2-1)+2, 6*(ip2-1)+3, 6*(ip2-1)+4, 6*(ip2-1)+5, 6*(ip2-1)+6 /)                           
                enddo
            enddo
            nnzAccumulate = nnzAccumulate + nnzCols
    !
            b1=6*(ip1-1)+1
            b2 = 36*kkk+1
            b3 = 6*size(krow(ip1)%col)
            ia( b1) = b2
            ia( b1+ 1) = b2+b3
            ia( b1+ 2) = b2+b3*2
            ia( b1+ 3) = b2+b3*3
            ia( b1+ 4) = b2+b3*4
            ia( b1+ 5) = b2+b3*5
    !        
            kkk=kkk+size(krow(ip1)%col)
        enddo      
    !
    !-------------------------------                                           ---------------------------------------
    !-------------------------------Initiate pardiso for the ROM linear system ---------------------------------------
    !-------------------------------                                           ---------------------------------------
    !
    !
        write(*,*) 'Start to initialize pardiso'
    !-- initiate pardiso and do the factorization
    !-->    set up pardisot control parameters           
        mtype     = 1  ! real and nonsymmetric
        solver    =  0  ! use sparse direct method
        msglvl= 0       ! print out message from pardiso
        iparm(1)=0
        iparm(3)=1
    !  .. PARDISO license check and initialize solver
!        call pardisoinit(pt, mtype, iparm)
        ALLOCATE(pt_mkl(64))

        DO iii = 1, 64
                pt_mkl(iii)%DUMMY =  0
        END DO
        do iii=1,64
                iparm(iii)=0
        enddo
        allocate(perm(NPart*NVEC))
        allocate(ddum(NVEC*NPART))
    !
    !---> only do the reording and symbolic factoriztion        
            phase=11
            
        
        write(*,*) 'pardiso reorder and symblic factorization started!'




        CALL pardiso (pt_mkl, maxfct, mnum, mtype, phase, NPart*NVEC, amatrix, ia,ja,      &
                 perm, nrhs, iparm, msglvl, ddum, ddum, errorpardiso)
        
    
        if (errorpardiso .ne. 0) then
            write(*,*) 'pardiso11 error in uexternaldb for ROM jacobian, errorcode= ', errorpardiso
            call exit
        endif 
        write(*,*) 'pardiso reorder and symblic factorization finished!'   
    !
    !-------------------------------                                                       
    !-------------------------------Initiate pardiso for solving the PlasticModuli DDSDDE  
    !
    !-- start to read the grain neigobor information and calcuate the number of non-zero conponents for the linear system for solving for calculating DDSDDE
        !--- Please note that ROM linear syatem and the linear system for solving the consistent moduli share the same sparsity pattern 
        !
        nnz_ddsdde=nnz
        allocate(ja_ddsdde(nnz_ddsdde), ia_ddsdde(NVEC*NPart+1))
        allocate(amatrix_ddsdde(nnz_ddsdde)) 
        ia_ddsdde=ia
        ja_ddsdde=ja          
    ! 
        kkk=0              
        nnzAccumulate = 0
        do ip1=1,NPart    ! ip1 is corresponding to beta in the equation  
            nnzCols = size(krow_ddsdde(ip1)%col)
            ! Find diagonal term
            do iii=1, size(krow_ddsdde(ip1)%col)
                ip2=krow_ddsdde(ip1)%col(iii)
                if (ip2 .eq. ip1) then
                    iiiDiag = iii
                end if
            end do
            do iii=1, size(krow_ddsdde(ip1)%col)     
                ip2=krow_ddsdde(ip1)%col(iii) !ip2 is corresponding to alpha in the equation
                delta=0.d0
                if (ip2 .eq.ip1) delta=pone
    !---- phi^(alpha,beta) 
                phiab=MATMUL(-COEFP(nnzAccumulate+iii, :, :),dmudotdsigma0)
                
    !---- eta^(alpha,beta)  
                I_P=Ident4th-COEFP(nnzAccumulate+iiiDiag, :, :) 
                etaab=MATMUL(I_P,dmudotdsigma0)
                etaab=COEFM(ip1, :, :)/dtime0+etaab 
    !         
    !-----Assemble A and F
                lhspart0=(1-delta)*phiab+delta*etaab
    !                
                do jjj=1,NVEC
                    b1= 36*kkk+1+6*(iii-1)+6*(jjj-1)*size(krow_ddsdde(ip1)%col)
                    b2= 36*kkk+6+6*(iii-1)+6*(jjj-1)*size(krow_ddsdde(ip1)%col)
                    amatrix_ddsdde(b1: b2)=lhsPart0(jjj,:)                          
                enddo                
            enddo
            kkk=kkk+size(krow_ddsdde(ip1)%col)        
            nnzAccumulate = nnzAccumulate + nnzCols    
        enddo
    !
        write(*,*) 'Start to initialize pardiso'
    !-- initiate pardiso and do the factorization
    !-->    set up pardisot control parameters           
        mtype_ddsdde     = 1  ! real and nonsymmetric
        solver_ddsdde    =  0  ! use sparse direct method
        msglvl_ddsdde=0        ! print out message from pardiso
        iparm_ddsdde(1)=0
        iparm_ddsdde(3)=1

        ALLOCATE(pt_mkl_ddsdde(64))

        DO iii = 1, 64
                pt_mkl_ddsdde(iii)%DUMMY =  0
        END DO
        do iii=1,64
                iparm_ddsdde(iii)=0
        enddo
        allocate(perm_ddsdde(NPart*NVEC))
        allocate(ddum_ddsdde(NVEC*NPART))
 
    !  .. PARDISO license check and initialize solver
    !
    !---> only do the reording and symbolic factoriztion        
        phase_ddsdde=11
        
        write(*,*) 'pardiso reorder and symblic factorization Started!'
        CALL pardiso (pt_mkl_ddsdde, maxfct_ddsdde, mnum_ddsdde, mtype_ddsdde,phase_ddsdde, NPart*NVEC, amatrix_ddsdde, ia_ddsdde, ja_ddsdde,      &
                 perm_ddsdde, nrhs_ddsdde, iparm_ddsdde, msglvl_ddsdde, ddum_ddsdde, ddum_ddsdde, errorpardiso_ddsdde)

        if (errorpardiso .ne. 0) then
            write(*,*) 'pardiso11 error in uexternaldb for ddsdde, errorcode= ', errorpardiso_ddsdde
            call exit
        endif 
        write(*,*) 'pardiso reorder and symblic factorization finished!'                
    
        RETURN
        
    END SUBROUTINE PardisoIni
