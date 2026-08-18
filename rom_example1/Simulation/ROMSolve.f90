
INCLUDE 'IntegrateCrystalEqns.f90'

!=====================================================================72
!
SUBROUTINE RomSolve(                 &
         stress,  ddsdde,strain,     & 
         dstrain, statev, nstatv,    &
         time, dtime, kinc, kstep,   &    
         pnewdt, argmin, argmax,     &
         geqvalues,                  &          
         gstress_n,                  &
         gestran_n,                  &
         gkappa_n,                   &
         gmu,                        &
         gmu_n                       & 
                                     )  
    USE      NumType
USE      FILEIO
USE      SlipGeo
USE      COEFTENS
USE      IterPar
USE      PlaPar
! ASLAN ADDED
use  NbData
    ! ASLAN ADDED
IMPLICIT NONE

!---------------------------------------------------------------------------------       
!
! aslan added
INTEGER      ::  iNG
INTEGER      ::  nstatv ,kinc, kstep
	REAL(KIND=8) ::  argmin, argmax, dtime, pnewdt, time(2)
	REAL(KIND=8) ::  stress(NVEC), ddsdde(NVEC, NVEC), PartStress(NVEC, NPart), &
                 strain(NVEC), dstrain(NVEC), statev(nstatv)

!---         
	REAL(KIND=8) ::  gstress_n  (NVEC, NPart)
	REAL(KIND=8) ::  gstress    (NVEC, NPart)
	REAL(KIND=8) ::  stress_n   (NVEC, NPart)
	REAL(KIND=8) ::  estran     (NVEC, NPart)
	REAL(KIND=8) ::  gestran    (NVEC, NPart)
	REAL(KIND=8) ::  estran_n   (NVEC, NPart)
	REAL(KIND=8) ::  gestran_n  (NVEC, NPart)
	REAL(KIND=8) ::  geqvalues  (NEQVA, NPart)
	REAL(KIND=8) ::  eqvalues   (NEQVA, NPart)
	REAL(KIND=8) ::  gmu        (NVEC, NPart)
	REAL(KIND=8) ::  gmu_n      (NVEC, NPart)
	REAL(KIND=8) ::  mu         (NVEC, NPart)
	REAL(KIND=8) ::  mu_n       (NVEC, NPart)
	REAL(KIND=8) ::  gkappa_n   (maxNumSlip, NPart)
	REAL(KIND=8) ::  kappa      (maxNumSlip, NPart)
	REAL(KIND=8) ::  gkappa     (maxNumSlip, NPart)
	REAL(KIND=8) ::  kappa_n    (maxNumSlip, NPart)
	REAL(KIND=8) ::  ggamdot    (maxNumSlip, NPart)
	REAL(KIND=8) ::  gamdot     (maxNumSlip, NPart)
    REAL(KIND=8) ::  grss       (maxNumSlip,NPart)
	REAL(KIND=8) ::  rss        (maxNumSlip,NPart)
	REAL(KIND=8) ::  gmudot     (NVEC, NPart)
	REAL(KIND=8) ::  mudot      (NVEC, NPart)

! ASLAN ADDED ----------------
! ASLAN ADDED
	REAL(KIND = 8) :: gfip(NPart, 4)
    REAL(KIND = 8) :: gsumofgamdot(NPart, 4)
! ASLAN ADDED
	REAL(KIND = 8) :: gfip_n(NPart, 4)
	REAL(KIND = 8) :: gsumofgamdot_n(NPart, 4)  
! ASLAN ADDED  
! ASLAN ADDED ----------------
!---      
	REAL(KIND=8) ::  savg_ij(NSD, NSD, NPart) 
	REAL(KIND=8) ::  cavg_ijkl(NVEC, NVEC, NPart)
	REAL(KIND=8) ::  c_ijkl(NVEC, NVEC, NPart)
	REAL(KIND=8) ::  s_ij(NSD, NSD, NPart) 
!	
INTEGER statusFlag  
!---
	REAL(KIND = 8) :: VecM(NSD,maxNumSlip,NPart)
	REAL(KIND = 8) :: VecS(NSD,maxNumSlip,NPart)
	REAL(KIND = 8) :: ZVec(NVEC,maxNumSlip,NPart)
! 
INTEGER ip,is, IPart1, IPart2, ierr, iterState
CHARACTER(LEN=30) errinfo

! aslan added
character*256 OUTDIR, timeStr
INTEGER LENOUTDIR
        INTEGER :: ipart, ivec
        !integer stressfile
! aslan added
! aslan added
	REAL(KIND = 8) :: fip(NPart, 4), fip_n(NPart, 4)
		REAL(KIND = 8) :: sumofgamdot(NPart, 4), sumofgamdot_n(NPart, 4)
! aslan added
!
!---
!    write(*,*) 'ROMSolve is called!'
!----------------------------------------------------------------          
!
    if (debug==1) then
        WRITE(FILE_E,*) '====***====     Begin ROMSolve    ====***===='
        Write(FILE_E,*) 'Dstrain  and strain is:'
        Write(FILE_E,'(6e18.8)') dstrain
        Write(FILE_E,'(6e18.8)') strain
        write(FILE_E,1000) time(1),dtime,kinc   
!
        write(FILE_E, '(A16, /, 12e24.16)')  'STATEV(13:24)=', STATEV(13:24)    
    endif
! ---------------------------------------------------------------------72
!-------- fetch state variables from abaqus vector array
! 

CALL RecoverStateVars(                                             &
             statev, nstatv, gstress_n, gestran_n,                     &
             gkappa_n, ggamdot, grss, gmudot, geqvalues, gmu ,gmu_n, Plap%PhID,    &
             gfip_n, gsumofgamdot_n                  )  !!! ASLAN ADDED



   
!

!!
!!!-------- compute state
statusFlag = XTAL_CONVERGED
!
!---   Initialize the streess of each part and system Jacobian
savg_ij=pzero
cavg_ijkl=pzero           
!
!--- Initialize some local arrays
!
    stress=pzero
    PartStress=pzero
    estran=pzero
    kappa=pzero
    mu=pzero

c_ijkl=pzero
rss=pzero
mudot=pzero
gamdot=pzero

    stress_n=pzero
    kappa_n=pzero
    mu_n=pzero
! ASLAN ADDED
fip = pzero
fip_n = pzero
sumofgamdot = pzero
sumofgamdot_n = pzero
! ASLAN ADDED
!    
!---Fetch Crystal Variables At Part  
!
    call FetchCrystalVariablesAtPart(                                  &
         gstress_n, gestran_n, gkappa_n, geqvalues,                    &
        stress_n, estran_n, kappa_n, eqvalues, mu_n, gmu_n            &                                
        ,gfip_n, gsumofgamdot_n, fip_n, sumofgamdot_n) !!! aslan added

!
!-----
!
!
CALL IntegrateCrystalEqns(                                         &
                    s_ij, partstress, estran, kappa,                       &
                    mu, mu_n, eqvalues, gamdot, rss, mudot,            &
                    stress_n, estran_n, kappa_n, dtime,                &
                    iterP%maxIterstate, iterP%MaxIterNewt, ierr, iterState, &
                    iterP%tolerState, iterP%tolerNewt,                 & 
                    c_ijkl, dstrain ,kinc, kstep, argmin, argmax,      &
              SlipG%VecM0, SlipG%VecS0, SlipG%ZVEC0, Plap%PhID,    &
  fip, fip_n, sumofgamdot, sumofgamdot_n           )  ! ASLAN ADDED
 


!
    if (ierr==XTAL_CONVERGED) then
        errinfo='Converged'
    elseif (ierr==XTAL_SING_JACOBIAN) then
        errinfo= 'SING_JACOBIAN'
    elseif (ierr==XTAL_MAX_ITERS_HIT) then   
        errinfo= 'MAX_ITERS_HIT'
    else
        errinfo= 'unknown error message, impossible!'  
        call RuntimeError(FILE_O, 'unknown error message, impossible!')
    endif     
        

!
    if (iterprint==1) write(iter_O, '(3(I6,6x), 6x, A20)'), kstep, kinc, iterstate, errinfo
!
    if (debug ==1) write(FILE_E,*) 'ierr= ', ierr
!
        IF (ierr .ne. XTAL_CONVERGED) then
                CALL WriteMessage(FILE_O, 'Resetting xtal quantities')


        WRITE(FILE_O, *) ' ** Umat did not converged       **'
        WRITE(FILE_O, *) ' ** re-scaling time step by 0.75 **'
        pnewdt = 0.75

RETURN        
ENDIF
!
CALL SaveCrystalVariablesAtPart(partstress, estran, kappa,             &
            gamdot, rss, mudot, eqvalues,  gstress, gestran, gkappa,   &
            gmu,mu, ggamdot, grss, gmudot, geqvalues                   &
, gfip, fip, gsumofgamdot, sumofgamdot)  !!! ASLAN
    
   
!------Check is this right or not!
        savg_ij=s_ij
        cavg_ijkl=c_ijkl

!!
        CALL SaveStateVars(statev, nstatv, gstress, gestran, gkappa,   &
                ggamdot, grss, gmudot, gmu, geqvalues, Plap%PhID              &
                , gfip, gsumofgamdot) !!! ASLAN


    if (debug==1) write(file_e,'(A30, 6(1x, e12.5))') 'kappa(1) in SDV after save is:', statev(13:24)    

!
!!-------- stresses and algorithmic moduli in abaqus format
!!
CALL SaveStressModuli(stress, ddsdde, savg_ij,  &
                           cavg_ijkl)
   
      
1000    FORMAT(/'*---time(1)---dtime---kinc---*'/,     &
             7x, F16.8,F16.8,i8)  

write(*,*) 'time(2) is', time(2)
        

        ! ASLAN ADDED FOR PRINTING
        !do ipart=1,Npart
               ! write(stressfile,'(F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6, F15.6)') gstress(1,ipart),gstress(2,ipart),gstress(3,ipart),gstress(4,ipart),gstress(5,ipart),gstress(6,ipart), gestran(1,ipart),gestran(2,ipart),gestran(3,ipart),gestran(4,ipart),gestran(5,ipart),gestran(6,ipart), gmu(1,ipart),gmu(2,ipart),gmu(3,ipart),gmu(4,ipart),gmu(5,ipart),gmu(6,ipart), gfip(ipart,1), gfip(ipart,2), gfip(ipart,3), gfip(ipart,4)
        !end do
        ! ASLAN ADDED FOR PRINTING 
RETURN
END
!
!=====================================================================72
!
!==============================================================
SUBROUTINE RecoverStateVars(                                           &
        statev, nstatv, gstress_n, gestran_n,                          &
        gkappa_n, ggamdot, grss, gmudot,geqvalues, gmu ,gmu_n, PhID         &
        ,gfip_n, gsumofgamdot_n)

      
    USE NumType
USE FILEIO
IMPLICIT NONE

      
INTEGER nstatv,kinc, PhID(NPart)
	REAL(KIND=8) ::  statev(nstatv)

	REAL(KIND=8) ::  gstress_n (NVEC, NPart)
	REAL(KIND=8) ::  gestran_n (NVEC,  NPart)
	REAL(KIND=8) ::  gkappa_n  (maxNumSlip,  NPart)
	REAL(KIND=8) ::  gcrot_n   (NSD, NSD,  NPart)
	REAL(KIND=8) ::  grrot_n   (NSD, NSD,  NPart)
	REAL(KIND=8) ::  gcrot0    (NSD, NSD,  NPart)


	REAL(KIND=8) ::  geqvalues (NEQVA,  NPart)
	REAL(KIND=8) ::  gmu        (NVEC,NPart)  
	REAL(KIND=8) ::  gmu_n        (NVEC,NPart)  
	REAL(KIND=8) ::  ggamdot   (maxNumSlip,  NPart)
	REAL(KIND=8) ::  grss        (maxNumSlip,NPart)
	REAL(KIND=8) ::  gmudot(NVEC,NPart)

! ASLAN ADDED
	REAL(KIND = 8) :: gfip_n(NPart, 4)
	REAL(KIND = 8) :: gsumofgamdot_n(NPart, 4)  
! ASLAN ADDED      

INTEGER varsPerPart1, varsPerPart2, ip, dex, id
      
	REAL(KIND=8) :: time(2)
!
!---------------------------------------------------------------------72

!---- number of state variables per ip
!
    varsPerPart1 =  2*NVEC        & ! stress, estran                1-12
                  + PhSlip(1)      & ! kappa                         13-24 
                  + NEQVA/2      & ! eqps, eqstr, gam_star, gamtot 25-27              
                  + PhSlip(1)      & ! gamdot                        28-39              
                  + NVEC         & ! mu                            40-45
                  + PhSlip(1)      & ! rss                           46-57
          + NVEC           &! mudot                         58-63  
  + 4   &! fip 						64-67
                  + 4           ! sumgamdot                   68-71

!    
    varsPerPart2 =  2*NVEC        & ! stress, estran                1-12
                  + PhSlip(2)      & ! kappa                         13-24 
                  + NEQVA/2      & ! eqps, eqstr, gam_star, gamtot 25-27              
                  + PhSlip(2)      & ! gamdot                        28-39              
                  + NVEC         & ! mu                            40-45
                  + PhSlip(2)      & ! rss                           46-57
          + NVEC           &! mudot                         58-63
  + 4   &! fip 							64-67
                  + 4           ! sumgamdot                   68-71
      
    if (debug == 1) write(FILE_E,*) '*----Initialize StateVars----*'
    if (debug == 1) write(FILE_E,*) 'varsPerPart1 is:', varsPerPart1  

!
!!---- recover state variables from abaqus vector array
!!
    dex=0
!    
DO ip = 1, NPart
      
dex = dex + 1  
gstress_n(:,ip)=statev(dex:(dex+NVEC-1))
!
        if (debug==1) then   
            write(FILE_E,*) 'gstress_n:'
            write(FILE_E,'(6E18.6)')  gstress_n(:,ip)
        endif     
!      
  dex = dex + NVEC                      ! estrain
gestran_n(:,ip)=statev(dex:(dex+NVEC-1))
        if (debug==1) then 
            write(FILE_E,*) 'estran_n:'
            write(FILE_E,'(6E18.6)')  gestran_n(:,ip) 
        end if
!
     dex = dex + NVEC                      ! kappa
gkappa_n(1:PhSlip(PhID(ip)),ip)=statev(dex:(dex+PhSlip(PhID(ip))-1))
!
        if (debug==1) then  
            write(FILE_E,*) 'gkappa_n:'
            write(FILE_E,'(12E18.6)') gkappa_n(:,ip)
        endif
!
        dex = dex + PhSlip(PhID(ip))    ! VonMises, Shearate, gamtot
        statev(dex:(dex+NEQVA/2-1))=pzero
!        
        if (debug==1) then  
            write(FILE_E,*) '*----Initial VonMises, Shearate, gamtot in ', ip, '-th part----'
            write(FILE_E,'(6F12.5)')  statev(dex:(dex+NEQVA/2-1))
        endif         
        
        dex = dex + NEQVA/2   
ggamdot(1:PhSlip(PhID(ip)),ip)=statev(dex:(dex+PhSlip(PhID(ip))-1))
!		
        if (debug==1) then  
            write(FILE_E,*) 'ggamdot:'
            write(FILE_E,'(9E18.6)')  ggamdot(:,ip)
        endif

!
        dex = dex + PhSlip(PhID(ip))  ! mu
gmu_n(:,ip)=statev(dex:(dex+NVEC-1)) 
        if (debug==1) then   
            write(FILE_E,*) 'gmu_n:'
            write(FILE_E,'(9E18.6)')  gmu_n(:,ip)  
        endif
!        
        dex = dex + NVEC                 ! rss  
grss(1:PhSlip(PhID(ip)),ip)=statev(dex:(dex+PhSlip(PhID(ip))-1))
!
        if (debug==1) then   
            write(FILE_E,*) 'grss:'
            write(FILE_E,'(9E18.6)')  grss(:,ip)  
        endif

        dex = dex + PhSlip(PhID(ip))                  ! mudot
gmudot(:,ip)=statev(dex:(dex+NVEC-1))
!		
        if (debug==1) then   
            write(FILE_E,*) 'gmu_n:'
            write(FILE_E,'(9E18.6)')  gmudot(:,ip)  
        endif

dex = dex + NVEC                  ! fip              
    gfip_n(ip,:) = statev(dex:(dex+4-1))

        dex = dex + 4                  ! sumgamdot              
    gsumofgamdot_n(ip,:) = statev(dex:(dex+4-1))

dex = dex + 4 -1
ENDDO
!
          
    RETURN
END
!
!=====================================================================72
!

!
!=====================================================================72
!
SUBROUTINE SaveStateVars(                                              &
        statev, nstatv, gstress, gestran, gkappa,                      &
         ggamdot, grss, gmudot, gmu, geqvalues, PhID                      &
         , gfip, gsumofgamdot)
         
    USE NumType
USE FILEIO
! aslan added
USE COEFTENS
! aslan added

     
IMPLICIT NONE
      
INTEGER nstatv
	REAL(KIND=8) ::  statev(nstatv)

	REAL(KIND=8) ::  gstress   (NVEC, NPart)
	REAL(KIND=8) ::  gestran   (NVEC, NPart)
	REAL(KIND=8) ::  gkappa    (maxNumSlip, NPart)

	REAL(KIND=8) ::  gmu       (NVEC,NPart)
	REAL(KIND=8) ::  ggamdot   (maxNumSlip, NPart)
	REAL(KIND=8) ::  grss        (maxNumSlip,NPart)

    REAL(KIND=8) ::  geqvalues (NEQVA, NPart)

	REAL(KIND=8) ::  gmudot(NVEC,NPart)

INTEGER varsPerPart1,varsPerPart2, ip, dex, id, PhID(NPart)


! ASLAN ADDED
	REAL(KIND = 8) :: gfip(NPart, 4)
	REAL(KIND = 8) :: gsumofgamdot(NPart, 4)  
! ASLAN ADDED      

!
!---------------------------------------------------------------------72
!

!---- number of state variables per ip
!
    !write(FILE_E, *) '*----Save State Variables!----*'
    
!---- initialize state variable vector per ip
!
!
    varsPerPart1 =  2*NVEC        & ! stress, estran                1-12
                  + PhSlip(1)      & ! kappa                         13-24 
                  + NEQVA/2      & ! eqps, eqstr, gam_star, gamtot 25-27              
                  + PhSlip(1)      & ! gamdot                        28-39              
                  + NVEC         & ! mu                            40-45
                  + PhSlip(1)      & ! rss                           46-57
          + NVEC           &! mudot                         58-63  
  + 4&! fip 						64-67
                  + 4           ! sumgamdot                   68-71
!    
    varsPerPart2 =  2*NVEC        & ! stress, estran                1-12
                  + PhSlip(2)      & ! kappa                         13-24 
                  + NEQVA/2      & ! eqps, eqstr, gam_star, gamtot 25-27              
                  + PhSlip(2)      & ! gamdot                        28-39              
                  + NVEC         & ! mu                            40-45
                  + PhSlip(2)      & ! rss                           46-57
          + NVEC           &! mudot                         58-63	
  + 4 &! fip						64-67
                  + 4           ! sumgamdot                   68-71
!
    if (nstatv .ne. ParPerPh(1)*varsPerPart1+ParPerPh(2)*varsPerPart2)      &
        call RunTimeError(FILE_O, 'nstatv .ne. NPart*varsPerPart')
!---- save state variables in abaqus vector array
!
    dex=0
    do ip = 1, NPart
         
        dex = dex + 1           ! stress
statev(dex:(dex+NVEC-1))=gstress(:,ip)
         
dex = dex + NVEC                      ! estrain
statev(dex:(dex+NVEC-1))=gestran(:,ip)
         
dex = dex + NVEC                      ! kappa
statev(dex:(dex+PhSlip(PhID(ip))-1))=gkappa(1:PhSlip(PhID(ip)),ip)

! ASLAN FIXED 484
        dex = dex + PhSlip(PhID(ip))                 ! VonMises, Shearate, gamtot
        statev(dex:(dex+NEQVA/2-1))=geqvalues ((NEQVA/2+1):NEQVA, ip)
! ASLAN FIXED 484         
dex = dex + NEQVA/2           ! gamdot
statev(dex:(dex+PhSlip(PhID(ip))-1))=ggamdot(1:PhSlip(PhID(ip)),ip)
         
dex = dex + PhSlip(PhID(ip))                 ! mu
statev(dex:(dex+NVEC-1))=gmu(:,ip) 

dex = dex + NVEC                  ! rss
statev(dex:(dex+PhSlip(PhID(ip))-1))=grss(1:PhSlip(PhID(ip)),ip)

dex = dex + PhSlip(PhID(ip))                  ! mudot
statev(dex:(dex+NVEC-1))=gmudot(:,ip)

dex = dex + NVEC                  ! fip              
    statev(dex:(dex+4-1))=gfip(ip,:)

        dex = dex + 4                  ! sumgamdot              
    statev(dex:(dex+4-1))=gsumofgamdot(ip,1:4)
       dex = dex + 4 -1
ENDDO


RETURN
END
!
!=====================================================================72
!
!=====================================================================72
!
SUBROUTINE SaveStressModuli(                              &
         stress, ddsdde, savg_ij, cavg_ijkl               &
         )                                                

    USE NumType
    USE COEFTENS
USE FILEIO
      
IMPLICIT NONE


	REAL(KIND=8) ::  stress(NVEC), ddsdde(NVEC, NVEC)
	REAL(KIND=8) ::  savg_ij(NSD, NSD,NPart), cavg_ijkl(NVEC, NVEC,NPart)

INTEGER i, j, ip
!
!---------------------------------------------------------------------72
!
!---- stresses
!
stress=pzero
	!Assign Volume Fraction to each part. They are assumed to be identical first, now read this info from CoefTens.dat


DO ip=1,NPart     
stress(1) = stress(1)+savg_ij(1,1,ip)*VolFrac(ip)
stress(2) = stress(2)+savg_ij(2,2,ip)*VolFrac(ip)
stress(3) = stress(3)+savg_ij(3,3,ip)*VolFrac(ip)
stress(4) = stress(4)+savg_ij(1,2,ip)*VolFrac(ip)
stress(5) = stress(5)+savg_ij(1,3,ip)*VolFrac(ip)
stress(6) = stress(6)+savg_ij(2,3,ip)*VolFrac(ip)
ENDDO

!
!---- algorithmic moduli
!
ddsdde=pzero
DO ip=1,NPart
ddsdde=ddsdde+ cavg_ijkl(:,:,ip)*VolFrac(ip) 
ENDDO
RETURN
END
!=====================================================================72
!
!
!
!=====================================================================72
!
SUBROUTINE ResetCrystalQnts(                                           &
stress, estran, kappa, statev, eqvalues,                       &
            gamdot, rss, mudot, mu, mu_n, s_ij, c_ijkl,            &
                stress_n, estran_n, kappa_n,                           &
                dtime, VecM, VecS,ZVec, argmin,argmax, PhID,  &
fip, sumofgamdot, fip_n, sumofgamdot_n             )  !!!

	USE FILEIO
    USE NumType

	IMPLICIT NONE 
	INTEGER  is, ip, PhID(NPart)
	REAL(KIND=8) ::  dtime, argmin,argmax

	REAL(KIND=8) ::  stress(NVEC, NPart), estran(NVEC, NPart), kappa(maxNumSlip, NPart)
	REAL(KIND=8) ::  statev(NSTAV,NPart )
	REAL(KIND=8) ::  s_ij(NSD,NSD,NPart), c_ijkl(NVEC,NVEC,NPart)

	REAL(KIND=8) ::  stress_n(NVEC,NPart), estran_n(NVEC,NPart), kappa_n(maxNumSlip,NPart)
	REAL(KIND=8) ::  gamdot(maxNumSlip,NPart), mu(NVEC,NPart),mu_n(NVEC,NPart)
	REAL(KIND=8) ::  rss(maxNumSlip,NPart), crss(maxNumSlip,NPart) 
	REAL(KIND=8) ::  InnerProductVec, SSKineticEqn
	REAL(KIND=8)    mudot(NVEC,NPart)
	REAL(KIND=8)    eqvalues(NEQVA,NPart)



	REAL(KIND = 8) :: VecM(3,maxNumSlip,NPart)
	REAL(KIND = 8) :: VecS(3,maxNumSlip,NPart)
	REAL(KIND = 8) :: ZTen(3,3,maxNumSlip,NPart),ZVec(6,maxNumSlip,NPart)
!---------------------------------------------------------------------------------


! ASLAN ADDED ----------------
! ASLAN ADDED
	REAL(KIND = 8) :: fip(NPart, 4)
    REAL(KIND = 8) :: sumofgamdot(NPart, 4)
! ASLAN ADDED
	REAL(KIND = 8) :: fip_n(NPart, 4)
	REAL(KIND = 8) :: sumofgamdot_n(NPart, 4)  
! ASLAN ADDED  
! ASLAN ADDED ----------------

!---------------------------------------------------------------------72
!---- Reset crystal quantities
!
stress=stress_n
estran=estran_n
kappa = kappa_n

! ASLAN ADDED ----------------
! ASLAN ADDED
fip(:,:) = fip_n(:,:)
sumofgamdot(:,:) = sumofgamdot_n(:,:)  
! ASLAN ADDED  
! ASLAN ADDED ----------------

eqvalues(kMISES,:)  = eqvalues(kMISES_n,:)             
eqvalues(kSHRATE,:) = eqvalues(kSHRATE_n,:)
eqvalues(kGAMTOT,:) = eqvalues(kGAMTOT_n,:)

!------- Cauchy stress (tensor form)
!
DO ip=1,NPart
CALL Vec6x1ToMat3x3Symm(stress(1,ip), s_ij(1,1,ip), NSD)
ENDDO
!
DO ip=1,NPart
mudot(:,ip)=pzero
DO is = 1, PhSlip(PhID(ip))
!------- Resolve shear stresses
rss(is,ip) = InnerProductVec(stress(1,ip), ZVec(1,is,ip), NVEC)  
crss(is,ip)= kappa(is,ip)
write(file_e,*) 'ResetQun: ', rss(is,ip), crss(is,ip)
!------- Shear strain rate
gamdot(is,ip) = SSKineticEqn(rss(is,ip),crss(is,ip),kGAMDOT, argmin, argmax, ip, PhID)
mudot(:,ip)=mudot(:, ip)+gamdot(is,ip)*ZVec(:, is, ip) 
ENDDO
ENDDO

!------- consistent tangent
CALL PlasticModuli(                                &
                        c_ijkl, stress, kappa,                &
                        dtime, VecM, VecS,ZVec, PhID               ) 
RETURN
END
!
!=====================================================================72
!
!=====================================================================72
!
      SUBROUTINE SaveCrystalVariablesAtPart(                           &
         stress, estran, kappa,                                        &
         gamdot, rss, mudot, eqvalues, gstress, gestran, gkappa,       &
         gmu, mu,ggamdot, grss, gmudot, geqvalues   &
 , gfip, fip, gsumofgamdot, sumofgamdot)
 

    Use NumType
USE FILEIO
IMPLICIT NONE


INTEGER igrn, iqpt, ielem,ip

	REAL(KIND=8) ::   gstress  (NVEC, NPart)
	REAL(KIND=8) ::   gestran  (NVEC, NPart)
	REAL(KIND=8) ::   gkappa   (maxNumSlip, NPart)
	REAL(KIND=8) ::   eqvalues(NEQVA, NPart)
	REAL(KIND=8) ::   geqvalues(NEQVA, NPart)
	REAL(KIND=8) ::   gmu    (NVEC,NPart)
	REAL(KIND=8) ::   ggamdot  (maxNumSlip, NPart)
	REAL(KIND=8) ::   grss        (maxNumSlip,NPart)

	REAL(KIND=8) ::   stress(NVEC, NPart), estran(NVEC, NPart), kappa(maxNumSlip, NPart)
	REAL(KIND=8) ::   gamdot(maxNumSlip, NPart),mu(NVEC,NPart), rss(maxNumSlip,NPart)

	REAL(KIND=8)    gmudot(NVEC,NPart),mudot(NVEC,NPart)



! ASLAN ADDED ----------------
! ASLAN ADDED
	REAL(KIND = 8) :: gfip(NPart, 4)
    REAL(KIND = 8) :: gsumofgamdot(NPart, 4)
! ASLAN ADDED
	REAL(KIND = 8) :: fip(NPart, 4)
	REAL(KIND = 8) :: sumofgamdot(NPart, 4)  
! ASLAN ADDED  
! ASLAN ADDED ----------------
!
!---------------------------------------------------------------------72
!


! ASLAN ADDED ----------------
! ASLAN ADDED
gfip(:,:) = fip(:,:)
    gsumofgamdot(:,:) = sumofgamdot(:,:)
! ASLAN ADDED  
! ASLAN ADDED ----------------


    gstress=stress
    gestran=estran
    gkappa =kappa
    ggamdot=gamdot
    grss   =rss
    gmudot =mudot
    gmu=mu
    
    geqvalues(kMISES,:) = eqvalues(kMISES,:)
    geqvalues(kSHRATE,:) = eqvalues(kSHRATE,:)
    geqvalues(kGAMTOT,:) = eqvalues(kGAMTOT,:)    

     
RETURN 
END
!
!=====================================================================72


!=====================================================================72
!
SUBROUTINE FetchCrystalVariablesAtPart(                                &
         gstress_n, gestran_n, gkappa_n, geqvalues,                    &
        stress_n, estran_n, kappa_n, eqvalues, mu_n, gmu_n             &                                
        ,gfip_n, gsumofgamdot_n, fip_n, sumofgamdot_n)
         
         

      
    USE       NumType
USE       FILEIO
      
      
IMPLICIT NONE

INTEGER   ielem, is, ip

	REAL(KIND=8) ::   gstress_n (NVEC, NPart)
	REAL(KIND=8) ::   gestran_n (NVEC, NPart)
	REAL(KIND=8) ::   gkappa_n  (maxNumSlip, NPart)
	REAL(KIND=8) ::   geqvalues (NEQVA, NPart)
	REAL(KIND=8) ::   gmu_n    (NVEC,NPart) 

	REAL(KIND=8) ::   stress_n(NVEC,NPart), estran_n(NVEC,NPart), kappa_n(maxNumSlip,NPart),mu_n(NVEC,NPart) 
	REAL(KIND=8) ::   eqvalues(NEQVA,NPart)

! ASLAN ADDED ----------------
! ASLAN ADDED
	REAL(KIND = 8) :: fip_n(NPart, 4)
    REAL(KIND = 8) :: sumofgamdot_n(NPart, 4)
! ASLAN ADDED
	REAL(KIND = 8) :: gfip_n(NPart, 4)
	REAL(KIND = 8) :: gsumofgamdot_n(NPart, 4)  
! ASLAN ADDED  
! ASLAN ADDED ----------------

!
!---------------------------------------------------------------------72

! ASLAN ADDED ----------------
! ASLAN ADDED
fip_n(:,:) = gfip_n(:,:)
    sumofgamdot_n(:,:) = gsumofgamdot_n(:,:)
! ASLAN ADDED  
! ASLAN ADDED ----------------


stress_n= gstress_n
estran_n= gestran_n
kappa_n = gkappa_n
mu_n=gmu_n

eqvalues(kMISES_n,:)  = geqvalues(kMISES_n,:)
eqvalues(kSHRATE_n,:) = geqvalues(kSHRATE_n,:)
eqvalues(kGAMTOT_n,:) = geqvalues(kGAMTOT_n,:)
    
!
RETURN
END
!
!=====================================================================72
!

