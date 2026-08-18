INCLUDE 'mkl_pardiso.f90'
!-----------------------------------------------------------------------
!***************
 module NumType
!***************

    ! This module defines the kind of integer and real numbers.
    ! Every module, subroutine or func must use this module.
    ! To change the precision from double to single,
    ! only this module needs to be changed.
    implicit none

    integer(kind(1)),parameter :: Ikind=kind(1),Rkind=kind(0.0d0)

    real(Rkind), parameter      :: pzero=0.d0,pone=1.d0,ptwo=2.d0,pthree=3.d0
    real(Rkind), parameter      :: pfour=4.d0,pfive=5.d0,psix=6.d0,pseven=7.d0
    real(Rkind), parameter      :: peight=8.d0,pnine=9.d0,pten=10.d0
    real(Rkind), parameter      :: pthird= 0.333333333333333d0, phalf = 0.5d0, ptwothrd= 0.666666666666667d0    
    real(Rkind), parameter      :: sqr23 = 0.816496580927726d0
    real(Rkind), parameter      :: sqr32 = 1.224744871391589d0
    real(Rkind), parameter      :: sqr2  = 1.414213562373095d0
    real(Rkind), parameter      :: sqr3  = 1.732050807568877d0 
    real(Rkind), parameter      :: Ident2nd(9) = (/1.d0, 0.d0, 0.d0,                &
                                                   0.d0, 1.d0, 0.d0,                &
                                                   0.d0, 0.d0, 1.d0/)
!                                              
    real(Rkind), parameter      ::   Ident(6) = (/1.d0, 1.d0, 1.d0, 0.d0,  0.d0,  0.d0/)
!
    real(Rkind), parameter      ::   Ident4th(6,6)= (/1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0,             &
                                                      0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0,             &
                                                      0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0,             &
                                                      0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0,             &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0,             &
                                                      0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1.d0/)  
!    
    real(Rkind), parameter      ::   factoTwo(6) = (/1.d0, 1.d0, 1.d0, 2.d0,  2.d0,  2.d0/)
    real(Rkind), parameter      ::   factoHalf(6) = (/1.d0, 1.d0, 1.d0, 0.5d0, 0.5d0, 0.5d0/) 
!
    real(ikind), parameter      ::   debug=10, iterprint=1     
!
    integer(ikind), parameter   :: NSTAV=5, NEQVA=6, NSD=3   
    integer(ikind), parameter   :: NPh=2, NPart=13, PhSlip(NPH)=(/12,12/), ParPerPh(NPh)=(/0,13/), maxnumslip=12
    integer(ikind), parameter   :: nvec=6
    integer(ikind), parameter   :: kFCC=1, kBCC=2, kHCP=3 
    integer(ikind), parameter   :: XTAL_CONVERGED=0, XTAL_SING_JACOBIAN=1,  &
                                   XTAL_LS_FAILED=2, XTAL_MAX_ITERS_HIT=3 
!
    integer(ikind), parameter   :: kMISES_n=1, kSHRATE_n=2, kGAMTOT_n=3, &
                                   KMISES  =4, kSHRATE  =5, kGAMTOT  =6                                   
!
    integer(ikind), parameter   :: kGAMDOT=1, kdGAMdTAU=2, kdGAMdKAPP=3, kdGamDOTdTau=4   
    integer(ikind), parameter   :: kHARD_EXPL=1, kHARD_MIDP=2, kHARD_ANAL=3 
!
    integer(ikind) :: print_bool
!                        
!
!--- A flag to indicate it is the first increment of an analysis of of an restarted analysis
!    
    logical :: FirstIncr                       
!
end module NumType
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
module AbaData
    ! print the results at element NO. numqpt_aba integration point at NO. numel_aba
integer, parameter  :: numel_aba=-1, numqpt_aba=-1
! 
end
!-----
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
module PardisoVar
!
    USE Numtype
! variables needed for pardiso
!..     Internal solver memory pointer 
    INTEGER(ikind) pt(64)
!
!..     All other variables 
    INTEGER(Ikind) maxfct, mnum, mtype, phase, nrhs, errorpardiso, msglvl, nnz
    INTEGER(Ikind) iparm(64)
    INTEGER(Ikind), ALLOCATABLE :: ia( : ), ja(:) 
    REAL(Rkind)  dparm(64) 
    REAL(Rkind), ALLOCATABLE :: amatrix( : )
    INTEGER(Ikind) idum, solver, ddum

!.. predifined data

    DATA  nrhs /1/, maxfct /1/, mnum /1/
! 
end module PardisoVar
!-----
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
module WorkDir
!
!
!----Define the file path and root name
    character(len=255) :: jobname, outdir
    character(len=255), parameter :: filepre="test"
    character(len=255), parameter :: CoefTensFile = 'CoefTens'
    character(len=255), parameter :: textureFile='Texture'
    character(len=255), parameter :: PhasePart='ph'
! ASLAN ADDED 
    character(len=255), parameter :: GNFile='GNMap'
    character(len=255), parameter :: StressStrain='ss'
! ASLAN ADDED 
!  
end module WorkDir
!
!-----------------------------------------------------------------------
MODULE COEFTENS
!
!
    INTEGER :: MTDIM
!
    INTEGER, ALLOCATABLE :: MT(:)
!
    REAL(KIND=8), ALLOCATABLE :: VolFrac(:)    
!
    REAL(KIND=8), ALLOCATABLE :: COEFA(:,:,:)
    ! aslan commented out
    REAL(KIND=8), ALLOCATABLE :: COEFP(:,:,:)
    REAL(KIND=8), ALLOCATABLE :: COEFM(:,:,:)    
! 
END
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
Module FILEIO
    ! aslan changed
    INTEGER :: FILE_I, FILE_E, FILE_O, TXT_I, CoefTens_I, iter_O, Ph1PartFile, Ph2PartFile, NbDataFile, OutF(2), GN_I, stressfile
    ! aslan changed
    parameter (FILE_I=80, FILE_O=81,FILE_E=82,                       &
                    TXT_I=83, CoefTens_I=84, iter_O=85, Ph1PartFile=86, Ph2PartFile=87, &
! ASLAN ADDED 
                    NbDataFile=88, OutF=(/89, 90/), GN_I=92, stressfile=93)
! ASLAN ADDED 
END
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
module timing
  integer c1, c2, cm, cr
  real(kind=8) :: tss, ths, tsj, longest_step,printed_step
  integer :: most_psis,printed_psis
  integer :: most_updates,printed_updates
  logical :: printing

end module
!
!-----------------------------------------------------------------------
module DataType
    use numtype
! ASLAN ADDED 
!-----------------------------------------------------------------------
    type NBG
        integer, ALLOCATABLE :: NBGID(:)      
        integer  NNG
    end type NBG
!-----------------------------------------------------------------------
! ASLAN ADDED 
!-----------------------------------------------------------------------
    type POWERLOW_BOUNDPar
        REAL(KIND=8) :: argMin
        REAL(KIND=8) :: argMax         
    end type POWERLOW_BOUNDPar  
!-----------------------------------------------------------------------
    type XtalPar
        REAL(KIND = 8) :: h0(NPh),xm(NPh),gam0(NPh),tausi(NPh),taus0(NPh),xms(NPh), gamss0(NPh), crss0(NPh)
        !REAL(KIND = 8), ALLOCATABLE :: kappa0(:,:)
        REAL(KIND = 8) :: kappa0(maxNumSlip,NPh)    
        integer(kind=ikind)   PhID(NPart), CrystalID(NPart)    
    end type XtalPar  
!-----------------------------------------------------------------------
    type SlipSys
        !real(kind=rkind), allocatable ::  VecM0(:, :, :), VecS0(:, :, :)
        real(kind=rkind)    ::  VecM0(NSD, maxNumSlip, NPart), VecS0(NSD, maxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  ZTen0(:, :, :,:)
        real(kind=rkind)    ::  ZTen0(NSD, NSD, maxNumSlip,NPart)
        !real(kind=rkind), allocatable ::  ZVec0(:, :, :)
        real(kind=rkind)    ::  ZVec0(NVEC, maxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  ZZT0(:,:,:,:)
        real(kind=rkind)    ::  ZZT0(NVEC,NVEC,maxNumSlip,NPart)
    end type SlipSys 
!-----------------------------------------------------------------------
    type OriData
        !real (kind=rkind), allocatable :: euler(:, :)
        real (kind=rkind)   :: euler(NSD, NPart)
        !real (kind=rkind), allocatable :: gcrot0(:, :, :)
        real (kind=rkind)   :: gcrot0(NSD, NSD, NPart)          
    end type OriData 
!-----------------------------------------------------------------------
    type xtalVars
        !real(kind=rkind), allocatable ::  gstress    (:, :)
        real(kind=rkind)    ::  gstress    (NVEC, NPart)
        !real(kind=rkind), allocatable ::  gestran    (:, :)
        real(kind=rkind)    ::  gestran    (NVEC, NPart)
        !real(kind=rkind), allocatable ::  gkappa     (:, :)
        real(kind=rkind)    ::  gkappa     (maxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  geqvalues  (:, :)
        real(kind=rkind)    ::  geqvalues  (NEQVA, NPart)
        !real(kind=rkind), allocatable ::  ggamdot    (:, :)
        real(kind=rkind)    ::  ggamdot    (NVEC, NPart)
        !real(kind=rkind), allocatable ::  gmu        (:, :) 
        real(kind=rkind)    ::  gmu        (NVEC, NPart)   
        ! ASLAN ADDED
        REAL(KIND = rkind) :: gfip(NPart, 4)
        REAL(KIND = rkind) :: gsumofgamdot(NPart, 4)  
        ! ASLAN ADDED
    end type xtalVars
!-----------------------------------------------------------------------
    type xtalVars_n
        !real(kind=rkind), allocatable ::  gstress_n    (:, :)
        real(kind=rkind)    ::  gstress_n    (NVEC, NPart)
        !real(kind=rkind), allocatable ::  gestran_n    (:, :)
        real(kind=rkind)    ::  gestran_n    (NVEC, NPart)
        !real(kind=rkind), allocatable ::  gkappa_n     (:, :)
        real(kind=rkind)    ::  gkappa_n     (maxNumSlip, NPart)
        !real(kind=rkind), allocatable ::  gmu_n        (:, :)    
        real(kind=rkind)    ::  gmu_n        (NVEC, NPart)   
        ! ASLAN ADDED
        REAL(KIND = rkind) :: gfip_n(NPart, 4)
        REAL(KIND = rkind) :: gsumofgamdot_n(NPart, 4)  
        ! ASLAN ADDED               
    end  type xtalVars_n
!-----------------------------------------------------------------------    
    type  IterData
        integer maxIterstate, MaxIterNewt
        real*8  tolerState, tolerNewt 
    end type  IterData
! ASLAN ADDED 
!--- this type stores the nonzero block positions for the jacobian 
    type SparseMap
        integer(ikind),pointer :: col(:) ! Each row could have different number of columns   
    end type SparseMap  
! ASLAN ADDED 
end module DataType 
!-----------------------------------------------------------------------
!
!
module IterPar
    use datatype
    type(IterData) iterP
end module IterPar
!-----------------------------------------------------------------------
module PlaPar
    use datatype
    type(XtalPar) Plap
end module PlaPar
!-----------------------------------------------------------------------
module SlipGeo
    use datatype
    type(SlipSys) SlipG
end module SlipGeo
!-----------------------------------------------------------------------
module OriPar
    use datatype
    type(OriData) OriP
end module OriPar
!-----------------------------------------------------------------------
module CPVars
!
    use datatype
    type(xtalVars) CPV
!
end module CPVars
!-----------------------------------------------------------------------
module CPVars_n
!
    use datatype
    type(xtalVars_n) CPV_n
!
end module CPVars_n
!-----------------------------------------------------------------------
module PowerLawBoundPar
!
    use datatype
    type(POWERLOW_BOUNDPar) PLBP
!
end module PowerLawBoundPar
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! ASLAN ADDED !-----------------------------------------------------------------------
!-----------------------------------------------------------------------
module NbData

    use numtype
    use Datatype
    type (NBG) Elset_nb(NPart)
  
end module NbData
!-----
!-----------------------------------------------------------------------
!-----
module PardisoROMVar
!
    use Numtype
    use DataType    
    use mkl_pardiso
! variables needed for pardiso
!..     Internal solver memory pointer 
    INTEGER(ikind) pt(64)
    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE::pt_mkl(:)
!
!..     All other variables 
    INTEGER(Ikind) maxfct, mnum, mtype, phase, nrhs, errorpardiso, msglvl, nnz
    INTEGER(Ikind) iparm(64)
    INTEGER(Ikind), ALLOCATABLE :: ia( : ), ja(:) 
    integer(ikind), ALLOCATABLE ::  perm( : )
!    REAL(Rkind), ALLOCATABLE ::  perm( : ) 
    REAL(Rkind), allocatable ::  ddum(:)

    REAL(Rkind)  dparm(64)
    REAL(Rkind), ALLOCATABLE :: amatrix( : )
    INTEGER(Ikind) idum, solver!, ddum
    type(SparseMap), allocatable :: Krow(:)   
!.. predifined data

    DATA  nrhs /1/, maxfct /1/, mnum /1/
! 
end module PardisoROMVar
!-----
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
module PardisoDDSDDEVar
!
    use Numtype
    use DataType    
    use mkl_pardiso
! variables needed for pardiso
!..     Internal solver memory pointer 
    INTEGER(ikind) pt_ddsdde(64)
    TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE::pt_mkl_ddsdde(:)
!
!..     All other variables 
    INTEGER(Ikind) maxfct_ddsdde, mnum_ddsdde, mtype_ddsdde, phase_ddsdde, &
                nrhs_ddsdde, errorpardiso_ddsdde, msglvl_ddsdde, nnz_ddsdde
    INTEGER(Ikind) iparm_ddsdde(64)
    INTEGER(Ikind), ALLOCATABLE :: ia_ddsdde( : ), ja_ddsdde(:) 
!    REAL(Rkind), ALLOCATABLE ::  perm_ddsdde( : ) 
    REAL(Rkind), ALLOCATABLE :: amatrix_ddsdde( : )
    integer(ikind), ALLOCATABLE ::  perm_ddsdde( : )
    REAL(Rkind), allocatable ::  ddum_ddsdde(:)


    INTEGER(Ikind) idum_ddsdde, solver_ddsdde!, ddum_ddsdde
    type(SparseMap), allocatable :: Krow_DDSDDE(:)
!.. predifined data

    DATA  nrhs_ddsdde /6/, maxfct_ddsdde /1/, mnum_ddsdde /1/
! 
end module PardisoDDSDDEVar
!-----
!-----------------------------------------------------------------------
