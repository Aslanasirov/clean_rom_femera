!
INCLUDE 'Modules.f90' 
INCLUDE 'uexternaldb.f90'
INCLUDE 'Utility.f90'
INCLUDE 'ROMSolve.f90'   
INCLUDE 'ROMUtility.f90'
!

!=====================================================================72
      
!
!  UMAT subroutine for ROM in Crystal Plasticity, this version is based 
!  on small deformation assumption such that the crystal orientation wi
!  -ll not change during the deformation history.
!
!=====================================================================72
!
    SUBROUTINE UMAT ( stress,  statev,  ddsdde,  sse,     spd,       &
                       scd,     rpl,     ddsddt,  drplde,  drpldt,   &
                       stran,  dstran, time,    dtime,   temp,       &
                        dtemp,   predef,  dpred,   cmname,  ndi,     &
                        nshr,    ntens,   nstatv,  props,   nprops,  &
                        coords,  drot,    pnewdt,  celent,  dfgrd0,  &
                        dfgrd1,  noel,    npt,     layer,   kspt,    &
                        kstep,   kinc )
!


    

    USE NumType
USE DataType
USE CPVars
use CPVars_n
use FileIO
    use PowerLawBoundPar
    implicit none                        
!---------------------------------------------------------------------72                        
!---- Variables passed into the UMAT sub
!
    character*80, intent(in)  :: cmname
    integer(kind=ikind), intent(in)  :: ntens, nstatv, nprops,         &
                                        ndi, nshr, noel, npt,          &
                                        kspt, kstep, kinc, layer
!                                        
    real(kind=rkind) :: sse, spd, scd, rpl, drpldt,                    &
                        dtime, temp, dtemp, pnewdt, celent
!
!---- Dimension arrays passed into the UMAT sub
!-----
!
    real(kind=rkind)        &
    stress(ntens),          &! Cauchy stress (vector form)
    statev(nstatv),         &! State variables
    ddsdde(ntens,ntens),    &! Tangent Stiffness Matrix
    ddsddt(ntens),          &! Change in stress per change in temperature
    drplde(ntens),          &! Change in heat generation per change in strain
    stran(ntens),           &! Strain tensor (vector form)
    dstran(ntens),          &! Strain increment tensor (vector form)
    time(2),                &! Time Step and Total Time
    predef(1),              &! Predefined state vars dependent on field variables
    dpred(1),               &! Change in predefined state variables
    props(nprops),          &! Material properties
    coords(3),              &! Coordinates of Gauss pt. being evaluated
    drot(3,3),              &! Incremental rotation matrix
    dfgrd0(3,3),            &! Deformation gradient at t_n
    dfgrd1(3,3)             ! Deformation gradient at t_(n+1)
!
    integer numel_aba, numqpt_aba    
! aslan added
    integer dex
    integer ipart
! aslan added

!
    type (xtalVars) CPV0
    type (xtalVars_n) CPV_n0
    type(POWERLOW_BOUNDPar) PLBP0
!---------------------------------------------------------------------72
!	
!

 
    CPV0=CPV
    CPV_n0= CPV_n
    PLBP0=PLBP 

!-----------------------------------------------------------
!--The first two Props record the total element number and number of integration points in each element.
        numel_aba  = nint (props(1))
        numqpt_aba = nint (props(2))
!
!--- Set up initial state at the first increment
!--- Since current code is small deformation with fix texture, 
!--- RotateSlipGeometry need only to be called onece,
!    other wise, it need to be called every time after 
!    the RecoverStateVars and before the ROM is solved
!   
     
!    print *, 'nstatv ', nstatv
    IF ( time(2)==0.d0  )  THEN 
         CALL SetUpStateVars (nstatv, statev)
   
    ENDIF
!
    !write(*,*) 'flag02'
!--- Solve ROM system
!
    CALL ROMSolve(                                      &  
                    stress,  ddsdde,stran,              & 
                    dstran, statev, nstatv,             &
                    time, dtime, kinc, kstep,           &    
                    pnewdt, PLBP0%argmin, PLBP0%argmax, &
                    CPV0%geqvalues,                     &                       
                    CPV_n0%gstress_n,                   &
                    CPV_n0%gestran_n,                   &
                    CPV_n0%gkappa_n,                    &
                    CPV0%gmu,                           &
                    CPV_n0%gmu_n                        &    
                                            )   
    if (debug==1) write(file_e,'(A30, 6(1x, e12.5))') 'Sigma(1) in SDV at the end of umat is:', statev(1:6)
    if (debug==1) write(FILE_E, '(A16, /, 12e24.16)')  'STATEV(13:24)=', STATEV(13:24)                                           
    if (debug==1) write(file_e, *) 'End of UMAT!'   

        ! this part assigns step and increments at which printing output is
        ! desired
        ! ALSAN ADDED
        print_bool=1
        if (KSTEP==1 .and. KINC==4) then
                print_bool=1
        else if (KSTEP==2 .and. KINC==4) THEN
                print_bool=1
        else if (KSTEP==3 .and. KINC== 8) THEN ! 
                print_bool=1
        else if (KSTEP==5 .and. KINC==4) then !
                print_bool=1
        else if (KSTEP==7 .and. KINC==4) then !
                print_bool=1
        else if (KSTEP==8 .and. KINC==41) then !
                print_bool=1
        else if (KSTEP==9 .and. KINC==4) THEN !
                print_bool=1
        endif
        
        if (print_bool==1) THEN
        ! ASLAN ADDED FOR PRINTING
                do ipart=1,Npart
                        dex = (ipart-1)*71
                        write(stressfile,'(I10, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.9, F15.9, F15.9, F15.9)') kstep, statev(dex+1), statev(dex+2), statev(dex+3), statev(dex+4),statev(dex+5),statev(dex+6),                 statev(dex+7), statev(dex+8), statev(dex+9), statev(dex+10),statev(dex+11),statev(dex+12),          statev(dex+40), statev(dex+41), statev(dex+42), statev(dex+43),statev(dex+44),statev(dex+45),       statev(dex+64), statev(dex+65), statev(dex+66), statev(dex+67)

                end do
        ! ASLAN ADDED FOR PRINTING 

        ENDIF

        ! ASLAN ADDED
RETURN 
END
!=====================================================================72
!
!
!=====================================================================72
!
! SET UP THE INITIAL STATE VARIABLES 
SUBROUTINE SetUpStateVars(  nstatv, statev                             &    
                            ) 

    USE NumType
    USE PlaPar
    USE OriPar 
    USE FileIO
! aslan added
    use WorkDir
    use NbData

          
        IMPLICIT NONE
!
!---------------------------------------------------------------------72
    integer nstatv
    real*8  statev(nstatv)
!
    integer dex, ip, varsPerPart1, varsPerPart2, id
!
!---------------------------------------------------------------------72
!
!---- initialize state variable vector per ip
!
    varsPerPart1 =  2*NVEC        & ! stress, estran                1-12
                  + PhSlip(1)      & ! kappa                         13-24 
                  + NEQVA/2      & ! eqps, eqstr, gam_star, gamtot 25-27              
                  + PhSlip(1)      & ! gamdot                        28-39              
                  + NVEC         & ! mu                            40-45
                  + PhSlip(1)      & ! rss                           46-57
                  + NVEC           &! mudot                         58-63
                  + 4           &! fip                         64-67
                  + 4           ! sumgamdot                   68-71

!    
    varsPerPart2 =  2*NVEC        & ! stress, estran                1-12
                  + PhSlip(2)      & ! kappa                         13-24 
                  + NEQVA/2      & ! eqps, eqstr, gam_star, gamtot 25-27              
                  + PhSlip(2)      & ! gamdot                        28-39              
                  + NVEC         & ! mu                            40-45
                  + PhSlip(2)      & ! rss                           46-57
                  + NVEC           &! mudot                         58-63
                  + 4           &! fip                         64-67	          
                  + 4           ! sumgamdot                   68-71
    if (debug == 1) write(FILE_E,*) '*----Initialize StateVars----*'
    if (debug == 1) write(FILE_E,*) 'varsPerPart1 is:', varsPerPart1, 'varsPerPart2 is:', varsPerPart2
    if (debug == 1) write(FILE_E,'(A30, /, 12(1x, e12.5))') 'Plap%kappa0(:,1) is: ', Plap%kappa0(:,1)

    if (nstatv .ne. ParPerPh(1)*varsPerPart1+ParPerPh(2)*varsPerPart2)  then
        call RunTimeError(FILE_O, 'nstatv .ne. NPart*varsPerPart')
    endif
!
    dex=0
    do ip = 1, NPart
!
        dex = dex + 1           ! stress (xtal)
        statev(dex:(dex+NVEC-1))=pzero
!
        if (debug==1) then         
            write(FILE_E,*) '*----Initial Stress in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))
        endif
!
        dex = dex + NVEC                      ! estrain
        statev(dex:(dex+NVEC-1))=pzero
!        
        if (debug==1) then          
            write(FILE_E,*) '*----Initial estrain in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))  
         endif     
!
        dex = dex + NVEC                      ! kappa         
        statev(dex:(dex+PhSlip(Plap%PhID(ip))-1))=Plap%kappa0(:,Plap%PhID(ip))
!        
        if (debug==1) then  
            write(FILE_E,*) 'PhSlip(Plap%PhID(ip))=', PhSlip(Plap%PhID(ip))         
            write(FILE_E,*) '*----Initial strength in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+PhSlip(Plap%PhID(ip))-1))
        endif
!
        dex = dex + PhSlip(Plap%PhID(ip)) ! VonMises, Shearate, gamtot
        statev(dex:(dex+NEQVA/2-1))=pzero
!        
        if (debug==1) then  
            write(FILE_E,*) '*----Initial VonMises, Shearate, gamtot in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+NEQVA/2-1))
        endif         
        
        dex = dex + NEQVA/2                                                           
        statev(dex:(dex+PhSlip(Plap%PhID(ip))-1))=pzero     ! gamdot
!        
        if (debug==1) then  
            write(FILE_E,*) '*----Initial gamdot in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+PhSlip(Plap%PhID(ip))-1))
        endif
!         
        dex = dex + PhSlip(Plap%PhID(ip))                 ! mu              
        statev(dex:(dex+NVEC-1))=pzero
!        
        if (debug==1) then  
            write(FILE_E,*) '*----Initial mu in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))
        endif        
!
         dex = dex + NVEC                   ! rss             
    statev(dex:(dex+PhSlip(Plap%PhID(ip))-1))=pzero
!        
        if (debug==1) then  
            write(FILE_E,*) '*----Initial rss in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+PhSlip(Plap%PhID(ip))-1))
        endif    
!
        dex = dex + PhSlip(Plap%PhID(ip))                  ! mudot              
    statev(dex:(dex+NVEC-1))=pzero
!        
        if (debug==1) then  
            write(FILE_E,*) '*----Initial mudot in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+NVEC-1))
        ENDIF
        
        dex = dex + NVEC                  ! fip              
    statev(dex:(dex+4-1))=pzero

        dex = dex + 4                  ! sumgamdot              
    statev(dex:(dex+4-1))=pzero
        dex = dex + 4 -1
ENDDO

RETURN
END
!
!========================================================
!
