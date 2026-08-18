    
!=====================================================================72
!
SUBROUTINE IntegrateCrystalEqns(                                       &
                s_ij, stress, estran, kappa,                           &
                mu, mu_n, eqvalues, gamdot, rss, mudot,                &
                stress_n, estran_n, kappa_n, dtime,                    &
                iterCounterS, iterCounterN,  ierr, iterState,          & 
                tolerState, tolerNewt,                                 &
                cepmod, dstrain,  kinc,   kstep,  argmin, argmax,      &
            VecM, VecS,ZVec, PhID ,                                &
                fip, fip_n, sumofgamdot, sumofgamdot_n           )  ! ASLAN ADDED 

    USE NumType
USE COEFTENS
USE FileIO
USE Timing
      
IMPLICIT NONE


	REAL(KIND=8) ::  s_ij(NSD,NSD,NPart)
	REAL(KIND=8) ::  Stress(NVEC, NPart),Stressini(NVEC, NPart), Stressvec(NVEC*NPart)
	REAL(KIND=8) ::  estran(NVEC,NPart), kappa(maxNumSlip,NPart)   
INTEGER(kind=iKind)   iterCounterS, iterCounterN, ierr
	REAL(KIND=8) ::  argmin, argmax, dtime, dstrain(NVEC)
	REAL(KIND=8) ::  search, tolerState, tolerNewt

      

	REAL(KIND=8) ::  gamdot(maxNumSlip,NPart), mu(NVEC,NPart), eqvalues(NEQVA, NPart)
	REAL(KIND=8) ::  stress_n(NVEC,NPart), estran_n(NVEC,NPart), kappa_n(maxNumSlip,NPart), mu_n(NVEC,NPart)

	REAL(KIND=8) ::  cepmod(NVEC,NVEC, NPart)
	REAL(KIND=8) ::  rss(maxNumSlip,NPart),  mudot(NVEC,NPart)


	LOGICAL converged, NConverged
	REAL(KIND=8) ::  norm_stress0(NPart), norm_kappa0(NPart), norm_stress(NPart), norm_kappa(NPart), normtau, normk, normtau0, normk0
! ASLAN VEC1 AND VEC2 ARE DEFINED BY XIAOYU BUT SEEMS LIKE ARE NOT USED

	REAL(KIND=8) ::  delstress(NVEC,NPart), stress0(NVEC,NPart)
	REAL(KIND=8) ::  delstressvec(NVEC*NPart)
    

	REAL(KIND=8) ::  rhs(NVEC*NPart)!, lhs(NVEC*NPart,NVEC*NPart)
	REAL(KIND=8) ::  rhs_norm0(NPart), rhs_norm(NPart)
    
	REAL(KIND=8) ::  ddd
! ASLAN VECMT ARE DEFINED BY XIAOYU BUT SEEMS LIKE ARE NOT USED

    INTEGER  c10, c20, PhID(NPart)

	LOGICAL ConvergeState 
	REAL(KIND=8) ::  InnerProductVec
    
INTEGER is, ip, iterState, kinc, INFO, kstep
! ASLAN
    INTEGER ip1, iplane
! ASLAN

! aslan added
	REAL(KIND = 8) :: fip(NPart, 4)
	REAL(KIND = 8) :: fip_n(NPart, 4)
    REAL(KIND = 8) :: sumofgamdot(NPart, 4), sumofgamdot_n(NPart, 4)
    REAL(KIND = 8) :: sigmap(NPart, 4)
	REAL(KIND = 8) :: sigmap_n(NPart, 4)
! aslan added

    REAL(KIND = 8) :: VecM(NSD,maxNumSlip,NPart)
    REAL(KIND = 8) :: VecS(NSD,maxNumSlip,NPart)
    REAL(KIND = 8) :: ZVec(NVEC,maxNumSlip,NPart)
	REAL(KIND = 8)  :: start_time, end_time
CHARACTER(LEN=20) ErroInfo
    integer iflag!, info
    !real*8, allocatable :: wa(:)    
!---------------------------------------------------------------------------------
      
    if (debug==1) write(FILE_E, *) '*----Begin IntegrateCrystalEqns!----*'
 
! ------- initialze slip hardening at time t
!
    kappa=kappa_n            
! ------- initial estimate for stress       
    stress=Stress_n
!---
    if (debug==1) then
        write(FILE_E, *) 'Initial estimate of stress(:,1) is:'
        write(FILE_E, '(6e24.16)') , stress(:,1)
        write(FILE_E, *) 'Initial estimate of kappa(:,1) is:'
        write(FILE_E, '(6e24.16)') , kappa(:,1)        
    endif    
!  
    stressini=stress       
!         
!! ------- initial valuies for the norm of stress and strength
!!
! 
!---> ZX cleaned-aa
!---> ZX cleaned-bb    
!
    normtau0= norm2(stress) 
    normk0  = norm2(kappa)    
!
! ------- initialized global flag to monitor Newton/State convergence
!
    ierr = XTAL_CONVERGED
!
! ------- iterate for the material state
!
    iterState    = 0
    converged = .false.
!
!---        Start the solving procedure, staggering      
!
    do while(iterState .lt. iterCounterS  .and. .not.converged )
!       
!    
        iterState = iterState + 1
        if (debug==1) write(File_e,*) 'iterState= ', iterState         
        stressvec=reshape(stress, shape(stressvec))
!---> ZX cleaned-aa        
!   
!
!  
!---> ZX cleaned-bb       
        call system_clock(c10)      
        call  StressSolve(            &
                rhs, stress, stressini, estran, rss, mudot, gamdot,           &
                kappa,  dtime, argmin, argmax, dstrain, VecM, VecS,ZVec, iterCounterN, tolerNewt, kstep, kinc, ierr, PhID)  
!
        call system_clock(c20)
        tss=tss+ (c20-c10)/REAL(cr) 
        write(File_e,*) 'ierr= ', ierr
        if (ierr .ne. XTAL_CONVERGED) return
!
!---> ZX cleaned-aa        
!---> ZX cleaned-bb   
        normtau= norm2(stress)
!               
!------- solve for hardness
!
        call system_clock(c10) 
        call IntegrateHardening(kappa, kappa_n, eqvalues, mu,mu_n, estran,       & 
                                  estran_n, stress, rss, mudot, gamdot, &
                                   dtime, kHARD_Expl, VecM, VecS,ZVec, argmin, argmax, PhID )
        call system_clock(c20) 
        ths=ths+ (c20-c10)/REAL(cr)                                            
!
!---> ZX cleaned-aa 
!---> ZX cleaned-aa 
        normk= norm2(kappa)                 
!
!---> ZX changed-aa 
            converged= convergestate(normtau, normk, normtau0, normk0, tolerState)
!---> ZX changed-aa         
!        
        ! ASLAN ADDED
        if (debug==1) write(file_E, *) '      Converged= ', converged


        !
        if (.not.converged) then 
!---> ZX changed-aa          
!---> ZX changed-bb 
           normtau0 = normtau
           normk0= normk 
        endif
!
    enddo
!      

    



    DO ip=1, NPart 
        
        sumofgamdot(ip, 1) = dabs(gamdot(1,ip)) + dabs(gamdot(2,ip)) + dabs(gamdot(3,ip))
        sumofgamdot(ip, 2) = dabs(gamdot(4,ip)) + dabs(gamdot(5,ip)) + dabs(gamdot(6,ip))
        sumofgamdot(ip, 3) = dabs(gamdot(7,ip)) + dabs(gamdot(8,ip)) + dabs(gamdot(9,ip))
        sumofgamdot(ip, 4) = dabs(gamdot(10,ip)) + dabs(gamdot(11,ip)) + dabs(gamdot(12,ip))

        do ip1=1,4
            iplane = (ip1-1)*3 + 1
            sigmap(ip,ip1) = stress(1,ip)*VecM(1,iplane,ip)*VecM(1,iplane,ip) + stress(2,ip)*VecM(2,iplane,ip)*VecM(2,iplane,ip) +stress(3,ip)*VecM(3,iplane,ip)*VecM(3,iplane,ip) + 2.0D0*stress(4,ip)*VecM(1,iplane,ip)*VecM(2,iplane,ip) + 2.0D0*stress(5,ip)*VecM(1,iplane,ip)*VecM(3,iplane,ip) + 2.0D0*stress(6,ip)*VecM(2,iplane,ip)*VecM(3,iplane,ip)

            sigmap_n(ip,ip1) = stress_n(1,ip)*VecM(1,iplane,ip)*VecM(1,iplane,ip) + stress_n(2,ip)*VecM(2,iplane,ip)*VecM(2,iplane,ip) +stress_n(3,ip)*VecM(3,iplane,ip)*VecM(3,iplane,ip) + 2.0D0*stress_n(4,ip)*VecM(1,iplane,ip)*VecM(2,iplane,ip) + 2.0D0*stress_n(5,ip)*VecM(1,iplane,ip)*VecM(3,iplane,ip) + 2.0D0*stress_n(6,ip)*VecM(2,iplane,ip)*VecM(3,iplane,ip)

            if (sigmap(ip,ip1) .lt. 0.0D0) sigmap(ip,ip1)  = 0.0D0
            if (sigmap_n(ip,ip1) .lt. 0.0D0) sigmap_n(ip,ip1)  = 0.0D0
            
            ! ASLAN CHANGED !!! 
            

        enddo
    ENDDO 

    DO ip=1, NPart  
    call Vec6x1ToMat3x3Symm(stress(1,ip), s_ij(1,1,ip), NSD)
        estran(:,ip)=estran_n(:,ip)+MATMUL(COEFM(ip, :, :), (stress(:,ip) -stressini(:,ip)))

        DO ip1 = 1,4
            fip(ip,ip1) = fip_n(ip,ip1) + ( sumofgamdot(ip,ip1)*(1.0D0 + 0.5D0*sigmap(ip,ip1)/138.0D0 ) + &
                sumofgamdot_n(ip,ip1)*(1.0D0 + 0.5D0*sigmap_n(ip,ip1)/138.0D0 ) )/2.0D0 * dtime
        enddo

    ENDDO 

    ! ! ASLAN CHECK
    ! nnztest=0
    ! ! ASLAN CHECK
    ! DO ip=1, NPart
    ! call Vec6x1ToMat3x3Symm(stress(1,ip), s_ij(1,1,ip), NSD)
    !     !estran(:,ip)=estran_n(:,ip)+MATMUL(COEFM(ip, :, :), (stress(:,ip) -stressini(:,ip)))
    !     ! ASLAN CHECK===============================================
    !     estran(:,ip)=estran_n(:,ip)+ MATMUL(COEFA(ip,:,:), dstrain)
    !     do ip1=1,NPart
    !             nnztest=nnztest+1
    !             estran(:,ip)=estran(:,ip)+ MATMUL(COEFP(nnztest,:,:),MUDOT(:,ip1))*dtime
    !     enddo
    !     ! ASLAN CHECK===============================================

    !     DO ip1 = 1,4
    !         fip(ip,ip1) = fip_n(ip,ip1) + ( sumofgamdot(ip,ip1)*(1.0D0 + 0.5D0*sigmap(ip,ip1)/138.0D0 ) + &
    !             sumofgamdot_n(ip,ip1)*(1.0D0 + 0.5D0*sigmap_n(ip,ip1)/138.0D0 ) )/2.0D0 * dtime
    !     enddo

    ! ENDDO


! ------- keep track of state iteration and check number of iterations
!
    if (iterState .ge. iterCounterS) then
        Write(FILE_O,*) 'Stress and Strength Solve: iters > maxIters when kinc=', kinc
        ierr = XTAL_MAX_ITERS_HIT
        return
    endif
!!
!
  
!
!----------------------- - CONSISTENT TANGENT
!
    call system_clock(c10) 
    call PlasticModuli(                                   &
                    cepmod, stress, kappa,                &
                    dtime, VecM, VecS,ZVec , PhID              ) 
    call system_clock(c20)    
    tsj=tsj+ (c20-c10)/REAL(cr)                 
!    
    if (debug==1) write(FILE_E,*) '----- End IntegrateCrystalEqns.f90' 
!    
      return
!---> ZX changed-aa        
!---> ZX changed-bb 
      !
END

!
!=====================================================================72
!
!
!=====================================================================72
SUBROUTINE StressSolve(            &
         rhs, stress, stressini, estran, rss, mudot, gamdot,           &
         kappa,  dtime, argmin, argmax, dstrain, VecM, VecS,ZVec, iterCounterN, tolerNewt, kstep, kinc, ierr, PhID) 
         
    use NumType
    use COEFTENS
    use FILEIO
!   ASLAN ADDED
    use PardisoROMVar
!   ASLAN ADDED

IMPLICIT NONE

    INTEGER      :: iterCounterN, kstep, kinc, info, ierr
	REAL(KIND=8) :: stress(NVEC,NPart), stress0(NVEC,NPart), stressini(NVEC, NPart)
	REAL(KIND=8) :: kappa(maxNumSlip,NPart)

	REAL(KIND=8) :: argmin, argmax, dtime,dstrain(NVEC), tolerNewt
	REAL(KIND=8) :: rhs(6*NPart), PartResid(6,NPart), rhsnorm!
	REAL(KIND=8) :: estran(NVEC, NPart)
	REAL(KIND=8) :: norm_stress0(NPart), norm_kappa0(NPart), norm_stress(NPart), norm_kappa(NPart)
        REAL(KIND=8) :: rhs_norm0
	REAL(KIND=8) ::  delstress(NVEC,NPart), delstressvec(NVEC*NPart)


INTEGER is,ip, i,j, ip1,ip2, iterNewt, PhID(NPart)
INTEGER indx(NVEC*NPart)
           
	REAL(KIND=8) ::  InnerProductVec, SSKineticEqn, search

	REAL(KIND=8) ::  crss(maxNumSlip, NPart),rss(maxNumSlip, NPart), gamdot(maxNumSlip, NPart)
      
!---> ZX cleaned-aa 
	REAL(KIND=8) ::  delIdent4th(6,6)
!				 
	REAL(KIND=8) ::   delta,term1(6), term2(6), term3(6),mudot(6,NPart), I_P(6,6) , mult1(6), mult2(6),temp1(6)
      
INTEGER StressSolveflat
	LOGICAL converged, converged2

    ! aslan added
!--> Xiaoyu: define rhs0 as pardiso input
    REAL(KIND=8) :: rhs0(6*NPart)
    ! aslan added
    REAL(KIND = 8) :: VecM(NSD,maxNumSlip,NPart)
    REAL(KIND = 8) :: VecS(NSD,maxNumSlip,NPart)
    REAL(KIND = 8) :: ZTen(NSD,NSD,maxNumSlip,NPart),ZVec(NVEC,maxNumSlip,NPart)
!---------------------------------------------------------------------------------    
!    
    if (debug==1) write(file_e,*) '-------------Begin Stress Solve-----------------'        
         
    call ComputeResidual(                        &
            rhs, stress,                         &
            stressini,                           &
            estran, rss, mudot,                  &
            gamdot, kappa,                       &
            dtime, argmin, argmax,               &
            dstrain,                          &  
            VecM, VecS,ZVec, PhID             )   
!
!
!---
    if (debug==1)   write(FILE_E, '(A30, /, 12e12.5)') 'Initial Residual is:', rhs
            
!
!! ------- initial valuies for the norm of stress of each part
!            
    rhs_norm0=norm2(rhs)   
!  
    iterNewt=0   
    converged=.false.            
    do while (iterNewt .lt. iterCounterN  .and. .not. converged )               
!
        iterNewt=iterNewt+1
        if (debug==1) write(file_e, *) 'iterNewt= ', iterNewt        
        stress0=stress  
            
!
! ---------- solve for the crystal stresses
!---------- compute local jacobian
        call ComputeJacobian( rss, kappa,       &
                                dtime, stress ,     &
                                VecM, VecS,ZVec, PhID   )
  
!  
        if (debug==1) then
            write(file_E, *) 'rhs= '      
            write(file_E, '( 6(2x, E12.5))') rhs            
        endif

! aslan added
!--> Xiaoyu: use pardiso solver
        rhs0 = rhs
        phase = 23
        CALL pardiso (pt_mkl, maxfct, mnum, mtype, phase, NPart*NVEC, amatrix, ia,ja,      &
            perm, nrhs, iparm, msglvl, rhs0, rhs, errorpardiso)

        if (errorpardiso .ne. 0) then
            write(*,*) 'pardiso factorization error in IntegrateCruystalEqns.f90, errorcode= ', errorpardiso
            call exit
        endif 

! aslan added

! aslan comment out
! !        
! !---------- solve for the increment of stress
!         !call CPU_TIME(start_time)
!         CALL DGESV( NPart*6, 1, lhs, NPart*6, indx, rhs, NPart*6, INFO )
! aslan comment out

	    !call CPU_TIME(end_time)
	    !Write(FILE_E,'(A)', ADVANCE='no') 'DGESV'
	    !Write(FILE_E,'(6e18.8)') end_time-start_time   

        INFO=1
	    if (debug==1) write(File_E,*) 'INFO= ', INFO
        if (INFO<0) then
            ierr=XTAL_SING_JACOBIAN
        !            
            write(File_o,*)                              &
                'StressSolveDeviatoric: Jacobian is singular when kinc=', kinc
            return
        endif         
!-------Please notice that here both delstress and rhs are 6*NPart by one vector
!-------we just rearrange them to be in the matrix form
        !write(File_E,*) 'Split the rhs into individual delata stress of each part'
        search = pone
        delstressVec=rhs 
!---> ZX cleaned-aa         
        delstress=reshape(delstressVec, shape(delstress))
!---> ZX cleaned-bb         
!                        
        stress=stress0-search*delstress        
    
    !
!---> ZX cleaned-aa      
!---> ZX cleaned-bb          
!
        call ComputeResidual(                        &
                rhs, stress,                         &
                stressini,                           &
                estran, rss, mudot,                  &
                gamdot, kappa,                       &
                dtime, argmin, argmax,               &
                dstrain,                          &  
                VecM, VecS,ZVec, PhID             )   
!
    if (debug==1)   write(FILE_E, '(A30, /, 12e12.5)') 'New Residual is:', rhs                
!---------- update stresses
!

        rhsnorm=norm2(rhs)

!       Add a linear search
        ! ASLAN CHANGED
        !do while (rhs_norm .gt. rhs_norm0  )   
        ! ASLAN CHANGED
        do while (rhsnorm .gt. rhs_norm0  )   
            search = search*0.5
            if (search .lt. TINY(1.d0)) then
               call WriteWarning(File_O,                               &
                    'StressSolveDeviatoric: LS Failed, search < TINY')
               ierr = XTAL_LS_FAILED
               return
            endif                  
            stress=stress0-search*delstress
!
!
            call ComputeResidual(                        &
                    rhs, stress,                         &
                    stressini,                           &
                    estran, rss, mudot,                  &
                    gamdot, kappa,                       &
                    dtime, argmin, argmax,               &
                    dstrain,                          &  
                    VecM, VecS,ZVec, PhID            ) 
!
            rhsnorm=norm2(rhs)
            write(file_e,*) 'line search is used'              
        enddo  
       ! ASLAN CHANGED
       !rhs_norm0 = rhs_norm 
       ! ASLAN CHANGED
       rhs_norm0 = rhsnorm 
       converged=converged2(rhs, tolerNewt, NPart*6)
   enddo
!
      if (iterNewt.ge. iterCounterN) then
         call WriteWarning(FILE_O,                                     &
                      'StressSolve: Netwo iters > iterCounterN')
         ierr = XTAL_MAX_ITERS_HIT
         return
      endif 
      
      return
end      
            
    
!
!=====================================================================72        
!
!=====================================================================72
!
SUBROUTINE IntegrateHardening(  &
                           kappa, kappa_n, eqvalues, mu,mu_n, estran, estran_n,  &
                           stress,rss, mudot, gamdot ,                 &
                        dtime, kInteg_Code, VecM, VecS,ZVec, argmin, argmax, PhID) 
!                           
    use NumType
    use PlaPar
    use FileIO    

    implicit none

    real*8 estran(6,NPart),estran_n(6,NPart)
    real*8 kappa(maxNumSlip,NPart), kappa_n(maxNumSlip,NPart)
    real*8 Stress(NVEC, NPart), eqvalues(NEQVA,NPart)
    real*8 crss(maxNumSlip,NPart),rss(maxNumSlip,NPart)
    real*8 gamdot(maxNumSlip,NPart), SHEARATE(maxNumSlip, NPart)
    real*8 mu(NVEC,NPart), mudot(NVEC,NPart)

    integer is, ip,  jk, kInteg_Code, PhID(NPart)

    real*8  argmin, argmax, kappa_sat
    real*8  c, g_n, g_s, g
    real*8  dkappa, gamtot_n, mu_n(NVEC,NPart), delgam, fac
    real*8  kTHETA
    data    kTHETA /1.0d0/
    real*8  dtime
       

    real*8  InnerProductVec, SSKineticEqn
    REAL(KIND = 8) :: VecM(3,maxNumSlip,NPart)
    REAL(KIND = 8) :: VecS(3,maxNumSlip,NPart)
    REAL(KIND = 8) :: ZTen(3,3,maxNumSlip,NPart),ZVec(6,maxNumSlip,NPart)
!---------------------------------------------------------------------72
!

    crss=pzero
    do ip=1,NPart
    mudot(:, ip)=pzero
        do is = 1, PhSlip(PhID(ip))
            rss(is,ip) = DOT_PRODUCT(stress(:,ip), ZVec(:,is,ip))  
            crss(is,ip)= kappa(is,ip)
            gamdot(is,ip) = SSKineticEqn(rss(is,ip),crss(is,ip),kGAMDOT, argmin, argmax, ip, PhID)
    mudot(:,ip)=mudot(:,ip) + gamdot(is, ip)*ZVec(:, is, ip)
!
      enddo
    enddo    
!
    do ip=1,NPart
        eqvalues(kSHRATE,ip) = pzero
    enddo    
!
    do ip=1,NPart
        do is = 1, PhSlip(PhID(ip))
            eqvalues(kSHRATE,ip) = eqvalues(kSHRATE,ip) + dabs(gamdot(is,ip))
        enddo

    enddo
!------- accumulated shear strain: gamtot 
!
    do ip=1,NPart
        eqvalues(kGAMTOT,ip) = eqvalues(kGAMTOT_n,ip) + eqvalues(kSHRATE,ip)*dtime
    enddo
!

            
!    
!------- inelastic strain mu and total elastic strain for each part.
    mu=mu_n
    do ip=1,NPart
        mu(:,ip) = mu_n(:,ip) + mudot(:,ip)*dtime
    enddo    
! 
    if (debug==1) then 
        write(file_e,*) 'gamdot(:,1): '
        write(file_e,'(12(1x, e12.5))') gamdot(:,1)    
        write(file_e,*) 'eqvalues(kshrate), eqvalues(kgamtot),  mu: '
        write(file_e,'(8(1x, e12.5))') eqvalues(kshrate,1), eqvalues(kgamtot,1), mu(:,1)
    endif
!------- integration of hardening law (one hardness/slip system)
!
!
    if (debug==1) write(FILE_E, *) 'kinteg_code =', kinteg_code
    if (kInteg_Code .eq. kHARD_EXPL) then         
        do ip=1,NPart
!
!---------- explicit update
        do is = 1, PhSlip(PhID(ip))
            kappa_sat = Plap%taus0(PhID(ip)) * ((eqvalues(kSHRATE,ip) / Plap%gamss0(PhID(ip)))**Plap%xms(PhID(ip))) !saturated stress
            c = dtime*Plap%h0(PhID(ip))
            g_s = kappa_sat - Plap%tausi(PhID(ip))
            g_n = kappa_n(is,ip) - Plap%tausi(PhID(ip))
!

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*(1.0-g_n/g_s)*eqvalues(kSHRATE,ip)             
            else
               g = g_n
            endif
            kappa(is,ip) = g + Plap%tausi(PhID(ip))
         enddo
      enddo
      

    else if (kInteg_Code .eq. kHARD_MIDP) then
!
!---------- generalized mid-point rule 
        do ip=1,NPart
            do is = 1, PhSlip(PhID(ip))
!
            kappa_sat = Plap%taus0(PhID(ip)) * ((eqvalues(kSHRATE,ip) / Plap%gamss0(PhID(ip)))**Plap%xms(PhID(ip)))
!
            c = dtime*Plap%h0(PhID(ip))
            g_s = kappa_sat - Plap%tausi(PhID(ip))
            g_n = kappa_n(is,ip) - Plap%tausi(PhID(ip))

            if ( (g_n/g_s) .le. 1.0 ) then
               g = g_n + c*(  &
                      (1.0-kTHETA)*(1.0-g_n/g_s)*eqvalues(kSHRATE_n,ip) &
                     +   (kTHETA)*eqvalues(kSHRATE,ip)  &
                           )
               g = g / (1.0 + c*kTHETA*eqvalues(kSHRATE,ip)/g_s)
            else
               g = g_n
            endif
            kappa(is,ip) = g + Plap%tausi(PhID(ip))
            enddo
        enddo
        

    else if (kInteg_Code .eq. kHARD_ANAL) then
          
        do ip=1,NPart
            gamtot_n = eqvalues(kGAMTOT_n,ip)
            delgam   = eqvalues(kSHRATE,ip) * dtime
            do is= 1, PhSlip(PhID(ip))
                kappa_sat = Plap%taus0(PhID(ip)) * ((eqvalues(kSHRATE,ip) / Plap%gamss0(PhID(ip)))**Plap%xms(PhID(ip)))
                g_s = kappa_sat - Plap%tausi(PhID(ip))
                fac = dabs(Plap%h0(PhID(ip))/g_s)

                dkappa = 0.0
                do jk = 1, PhSlip(PhID(ip))
                    dkappa = dkappa +    &
                    dabs(gamdot(jk,ip))*dtime/delgam 
                enddo

                dkappa = dkappa * g_s * exp(-gamtot_n*fac) *   &
                                              (1.0 - exp(-delgam*fac))
                if (debug==1) then 
                    write(file_e,*) 'kappa_sat, g_s,  fac, dkappa: ', kappa_sat, g_s,  fac, dkappa
!
                endif                                              
                kappa(is,ip) = kappa_n(is,ip) + dkappa
                if (debug==1) write(*,*) 'kappa(', is, ',', ip, ') is ', kappa(is,ip)
            enddo
        enddo
    else
!
!------- wrong code number
        call RunTimeError(FILE_O, 'IntegrateHardening: Wrong kInteg_Code!')

    endif
!
    return
    END
!
!=====================================================================72
!
! ---==================================================================72
!
SUBROUTINE ComputeResidual(            &
         rhs, stress, stressini, estran, rss, mudot, gamdot,           &
         kappa,  dtime, argmin, argmax, dstrain, VecM, VecS,ZVec, PhID      )      
         
          
       
    use NumType
    use COEFTENS
    use FILEIO
    use PardisoDDSDDEVar

IMPLICIT NONE

	REAL(KIND=8) :: stress(NVEC,NPart),stressini(NVEC, NPart)
	REAL(KIND=8) :: kappa(maxNumSlip,NPart)

	REAL(KIND=8) :: argmin, argmax, dtime,dstrain(NVEC)
	REAL(KIND=8) :: rhs(6*NPart), PartResid(6,NPart)
	REAL(KIND=8) :: estran(NVEC, NPart)

!   ASLAN CHANGED
INTEGER is,ip, i,j, ip1,ip2, PhID(NPart), nnzAccumulate, nnzCols, iii
!   ASLAN CHANGED
	REAL(KIND=8) ::  InnerProductVec, SSKineticEqn

	REAL(KIND=8) ::  crss(maxNumSlip, NPart),rss(maxNumSlip, NPart), gamdot(maxNumSlip, NPart)
      
!---> ZX cleaned-aa 
	REAL(KIND=8) ::  delIdent4th(6,6)
!---> ZX cleaned-bb 	
	REAL(KIND=8) ::   delta,term1(6), term2(6), term3(6),mudot(6,NPart), I_P(6,6) , mult1(6), mult2(6),temp1(6)
      
INTEGER StressSolveflat


    REAL(KIND = 8) :: VecM(NSD,maxNumSlip,NPart)
    REAL(KIND = 8) :: VecS(NSD,maxNumSlip,NPart)
    REAL(KIND = 8) :: ZTen(NSD,NSD,maxNumSlip,NPart),ZVec(NVEC,maxNumSlip,NPart)
!---------------------------------------------------------------------------------
!
!
    mudot=pzero
!      
    do ip=1,NPart
          do is = 1, PhSlip(PhID(ip))
             
             rss(is,ip) = InnerProductVec(stress(1,ip), ZVec(1,is,ip), NVEC)  
             crss(is, ip) = kappa(is,ip)          
             gamdot(is,ip) = SSKineticEqn(rss(is,ip), crss(is,ip), &
                                      kGAMDOT,  argmin, argmax, ip, PhID)
!---> ZX cleaned-aa                                      
             mudot(:,ip)=mudot(:,ip)+gamdot(is,ip)*ZVec(:,is,ip)
!---> ZX cleaned-aa             
          enddo
    enddo
!                                   
!---------------------------------------------------------------------72
!
!---> ZX cleaned-aa
    PartResid=pzero  
!---> ZX cleaned-bb    
!    
  !---> ZX cleaned-aa
    ! ASLAN ADDED
    nnzAccumulate = 0
    ! ASLAN ADDED
    do ip1=1,NPart !ip1 is corresponding to the index beta in the equation
        delta=pzero       
        term1=pzero
        term2=pzero
        term3=pzero

        ! ASLAN ADDED
        nnzCols = size(krow_ddsdde(ip1)%col)
        ! ASLAN ADDED

!        nnzCols = size(krow_ddsdde(ip1)%col)
        do iii=1, size(krow_ddsdde(ip1)%col)  
            ip2 = krow_ddsdde(ip1)%col(iii)
            delta=pzero 
            if (ip2.eq.ip1) delta=pone                 
            delIdent4th=delta*Ident4th
!
            I_P=delIdent4th-COEFP(nnzAccumulate+iii, :, :)
   
            mult1=MATMUL(I_P, mudot(:,ip2))

!
            term1=term1+mult1
        enddo
        nnzAccumulate = nnzAccumulate + nnzCols


        ! ASLAN COMMENTED OUT
        ! ASLAN COMMENTED OUT
          
        mult2=(stress(:,ip1)-stressini(:,ip1))/dtime
!
        term2=MATMUL(COEFM(ip1, :, :),mult2)

    term3=MATMUL(COEFA(ip1, :, :), dstrain)/dtime
!---> ZX cleaned-bb	    
!
!---> ZX changed-aa 
        rhs((6*(ip1-1)+1):6*ip1)=term1+term2-term3
!---> ZX changed-bb         
!
      enddo
      
      
!----Write the residual in each partition into the whole 6*NPart by 1 rhs vector
! 
!---> ZX cleaned-aa 
!---> ZX cleaned-bb 
    return
END
!
! ---==================================================================72
!

!=====================================================================72
    !
SUBROUTINE ComputeJacobian( rss, kappa,       &
                               dtime,  stress,      &
                           VecM, VecS,Zvec, PhID   )  
    use NumType
    use COEFTENS
    use PlaPar
    use FILEIO
    use PardisoROMVar

    implicit none

    real*8  dtime
    real*8  rss(maxNumSlip,NPart)
    real(kind=8)  kappa(maxNumSlip, NPart)
! ASLAN CHANGED
    real(kind=8) lhsPart0(6,6)
! ASLAN CHANGED
    real(kind=8) dmudotdsigma(6,6,NPart), I_P(6,6), Term1(6,6), Term2(6,6)
    real(kind=8) delta

    ! ASLAN CHANGED
    integer is, ip, ip1,ip2, i, PhID(NPart), nnzAccumulate, nnzCols, iii
    integer jjj,kkk, b1, b2
    ! ASLAN CHANGED
    real*8  crss(maxNumSlip, NPart)
    real*8  stress(6, NPart)

!---> ZX cleaned-bb
!---> ZX cleaned-bb    

    real(kind=8)  temp2(6,6),omega(maxNumSlip, NPart)
    real*8 InnerProductVec

    ! ASLAN CHANGED
    real*8 SignOf
    ! ASLAN CHANGED

    REAL(KIND = 8) :: VecM(3,maxNumSlip,NPart)
    REAL(KIND = 8) :: VecS(3,maxNumSlip,NPart)
    REAL(KIND = 8) :: ZTen(3,3,maxNumSlip,NPart),ZVec(6,maxNumSlip,NPart)
!---------------------------------------------------------------------------------


! ---------------------------------------------------------------------72
    
    
    dmudotdsigma=pzero
    omega=pzero
!---> ZX cleaned-aa    
!---> ZX cleaned-aa    
    
    
! ------- omega and dudotdsigma
!
    do ip=1, NPart

        do is = 1, PhSlip(PhID(ip))          
            omega(is,ip) = PlaP%gam0(PhID(ip))/PlaP%xm(PhID(ip))/kappa(is,ip)* (dabs(rss(is,ip))/kappa(is,ip))**(1/PlaP%xm(PhID(ip))-1)
            call OuterProductVec(ZVec(1,is,ip),ZVec(1,is,ip),temp2,6)  
            dmudotdsigma(:,:,ip)=dmudotdsigma(:,:,ip)+omega(is,ip)*temp2           
        enddo
    enddo 
!
!
! ------- local jacobian
!
!  

    kkk=0           
    nnzAccumulate = 0      
    do ip1=1,NPart !ip1 is corresponding to the index beta in the equation
    nnzCols = size(krow(ip1)%col)
        do iii=1, size(krow(ip1)%col) ! ip2 is corresponding to the index eta in the equation
            ip2=krow(ip1)%col(iii)
            delta=pzero
            if (ip2 .eq. ip1) delta=pone 
            I_P=delta*Ident4th-COEFP(nnzAccumulate+iii, :, :)
            Term1=MATMUL(I_P, dmudotdsigma(:,:,ip2))
            lhsPart0=delta*COEFM(ip1, :,:)/dtime+Term1
            do jjj=1,NVEC
                b1= 36*kkk+1+6*(iii-1)+6*(jjj-1)*size(krow(ip1)%col)
                b2= 36*kkk+6+6*(iii-1)+6*(jjj-1)*size(krow(ip1)%col)
                Amatrix(b1: b2)=lhsPart0(jjj,:)        
            enddo
        enddo
        nnzAccumulate = nnzAccumulate + nnzCols
!       
        kkk=kkk+size(krow(ip1)%col)
    enddo  


    
!---Feed each 6*6 matrix into the whole (6*NPart) by 6*NPart) matrix
!
!---> ZX cleaned-aa 
!---> ZX cleaned-bb     
!
 return  
end
 !=====================================================================72
!
SUBROUTINE PlasticModuli(                                 &
                        cepmod, stress, kappa,                &
                        dtime, VecM, VecS,ZVec, PhID             ) 
    use NumType
    use COEFTENS
    use PlaPar
    USE FILEIO
    use PardisoDDSDDEVar
   
    implicit none
      
      
    real*8 cepmod(NVEC,NVEC, NPart)
    real*8  Stress(NVEC, NPart)
    real*8  kappa(maxNumSlip,NPart) 
    real*8  dmudotdsigma(6,6,NPart),dmudotdrss
    real*8  rss(maxNumSlip,NPart), ZZT(6,6,maxNumSlip,NPart)
    integer ip, ip1,ip2,is,i,j, PhID(NPart)
!---> ZX cleaned-aa    
    real*8  delta, I_P(6,6)   
!---> ZX cleaned-bb    
    real*8  InnerProductVec
    ! aslan added
    real*8, ALLOCATABLE :: F(:,:)
    real*8  SignOf
    ! aslan added
    real*8  phiab(6,6),etaab(6,6),fb(6,6)
    real*8  dtime
    integer indx(6*NPart)
    real*8  ddd
    integer Cepmodflag

    REAL(KIND = 8) :: VecM(3,maxNumSlip,NPart)
    REAL(KIND = 8) :: VecS(3,maxNumSlip,NPart)
    REAL(KIND = 8) :: ZTen(3,3,maxNumSlip,NPart),ZVec(6,maxNumSlip,NPart)
 	REAL(KIND = 8) :: start_time, end_time

     ! ASLAN ADDED
     !--> Xiaoyu: define variables for pardiso
    real*8, ALLOCATABLE :: F0(:,:)
    real*8 lhspart0(NVEC, NVEC)
    integer iii, jjj, kkk, b1, b2, nnzAccumulate, nnzCols, iiiDiag
    ! ASLAN ADDED
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------


      
!---- Trial by using a constant Jacobian
    
    
!----Caculate d(sigma_n+1^(beta))/d(epsilonbar_(n+1))
! ASLAN CHANGED
!--> Xiaoyu: allocate F0
    allocate (F(6*NPart,6))
    Cepmodflag=0
    dmudotdsigma=pzero
    dmudotdrss=pzero
! ASLAN CHANGED
    allocate (F0(6*NPart,6))

! ASLAN CHANGED
!---> ZX cleaned-aa     
!---> ZX cleaned-bb     

!----Caculate the rss    
    do ip=1,NPart
        do is=1,PhSlip(PhID(ip))
            rss(is,ip)=InnerProductVec(stress(:,ip),ZVec(:,is,ip),6)
            Call OuterProductVec(ZVec(:,is,ip),ZVec(:,is,ip),ZZT(:,:,is,ip),6)
            dmudotdrss=PlaP%gam0(PhID(ip))/PlaP%xm(PhID(ip))/kappa(is,ip)*(dabs(rss(is,ip))/kappa(is,ip))**(1/PlaP%xm(PhID(ip))-1)
            dmudotdsigma(:,:,ip)=dmudotdsigma(:,:,ip)+dmudotdrss*ZZT(:,:,is,ip)
        enddo
    enddo
! ASLAN ADDED
    !--> Xiaoyu: use pardiso to solve for ddsdde
    kkk=0             
    nnzAccumulate = 0     
    do ip1=1,NPart      ! ip1 is corresponding to beta in the equation  
    nnzCols = size(krow_ddsdde(ip1)%col)
    ! Find diagonal term
    do iii=1, size(krow_ddsdde(ip1)%col)
            ip2=krow_ddsdde(ip1)%col(iii)
            if (ip2 .eq. ip1) then
            iiiDiag = iii
            end if
        end do
        do iii=1, size(krow_ddsdde(ip1)%col)
            ip2=krow_ddsdde(ip1)%col(iii)   !ip2 is corresponding to alpha in the equation
            delta=pzero
            if (ip2 .eq. ip1) delta=pone 
            phiab=MATMUL(-COEFP(nnzAccumulate+iii, :, :),dmudotdsigma(:,:,ip2))
            I_P=Ident4th-COEFP(nnzAccumulate+iiiDiag, :, :)
            etaab=MATMUL(I_P,dmudotdsigma(:,:,ip1))
            etaab=COEFM(ip1, :, :)/dtime+etaab             
            lhsPart0=(1-delta)*phiab+delta*etaab
            do jjj=1,NVEC
                b1= 36*kkk+1+6*(iii-1)+6*(jjj-1)*size(krow_ddsdde(ip1)%col)
                b2= 36*kkk+6+6*(iii-1)+6*(jjj-1)*size(krow_ddsdde(ip1)%col)
                amatrix_ddsdde(b1: b2)=lhsPart0(jjj,:)                          
            enddo
            F((6*(ip1-1)+1):6*ip1,:)=COEFA(ip1, :, :)/dtime        
        enddo
        nnzAccumulate = nnzAccumulate + nnzCols
!        
        kkk=kkk+size(krow_ddsdde(ip1)%col)
    enddo
    F0=F
    phase_ddsdde = 23
! ASLAN ADDED
    CALL pardiso (pt_mkl_ddsdde, maxfct_ddsdde, mnum_ddsdde, mtype_ddsdde, &
                phase_ddsdde, NPart*NVEC, amatrix_ddsdde, ia_ddsdde,   &
                ja_ddsdde, perm_ddsdde, nrhs_ddsdde, iparm_ddsdde,     &
                msglvl_ddsdde, F0, F, errorpardiso_ddsdde)

    ! ASLAN ADDED
! ASLAN COMMENTED    

!---- Solve for cepmod as a whole and cepmod in each part is stored in cepmodm  
!
!---- Split cepmodm into   cepmod(:,:,ip)   
    do ip=1,NPart
        cepmod(:,:,ip)=F((6*(ip-1)+1):6*ip,:)
    enddo
    ! ASLAN ADDED
    deallocate (F)
    ! ASLAN ADDED
!   
    return
END
!
!=====================================================================72
!
!=====================================================================72
!
logical FUNCTION Converged2(                                           &
        res, toler, n                                                  &
        )
    use NumType
    implicit none

    real*8  toler
    real*8  res(n)

    integer i,n
!
!---------------------------------------------------------------------72
!     
!------- check convergence on residual
!
    Converged2 = ( dabs(res(1)) .lt. toler )
    do i = 2, n
        Converged2 = ( (dabs(res(i)) .lt. toler) .and. Converged2)
    enddo

    return
END
!     
!=====================================================================72
!
!=====================================================================72
!
      logical FUNCTION ConvergeState(                                  &
         norm_a, norm_b, norm_a0, norm_b0, toler                       &
         )

      implicit none

      real*8  norm_a, norm_b, norm_a0, norm_b0, toler
!
!---------------------------------------------------------------------72
!
      ConvergeState = (dabs(norm_a - norm_a0) .lt. toler*norm_a0) .and. & 
                      (dabs(norm_b - norm_b0) .lt. toler*norm_b0)

      ! ASLAN ADDED
      if (norm_a0 .lt. 1e-12) then
        ConvergeState = (dabs(norm_a - norm_a0) .lt. toler) .and. &
                      (dabs(norm_b - norm_b0) .lt. toler*norm_b0)
      endif
      ! ASLAN ADDED


      if (.not.ConvergeState) then
          norm_a0 = norm_a
          norm_b0 = norm_b
      endif

      return
      END
!
!=====================================================================72
!
!=====================================================================72
!
      SUBROUTINE WriteWarning(                                         &
         io, message                                                   &
         )

      implicit none
      
      character*(*) message
      integer io
!
!---------------------------------------------------------------------72
!
      write(io, 1000) message

1000  format(/,'***WARNING Message: '/, 3x, a)

      return
      END
!
!=====================================================================72
