!=====================================================================72
!
      real*8 FUNCTION SSKineticEqn(                             &
         rss1, crss1, kflag, argmin, argmax, ip, PhID           &
         )
      
      use NumType
      use PlaPar
      use FILEIO

      implicit none


      integer kflag, ip, PhID(NPart)
      real*8  crss1, rss1, pow,arg, argmin, argmax



      real*8  CheckArgPowLaw, SignOf, Power
!
!---------------------------------------------------------------------72
!


!
!------- check argument of power law
    !write(*,*) 'PhID= ', PhID
!
    arg = rss1/crss1
    if (dabs(arg) .ne. 0) arg = CheckArgPowLaw(arg, SignOf(arg), argmin, argmax)

      pow = Power(dabs(arg), 1./PlaP%xm(PhID(ip))-1.)
      ! write(        file_e,*), 'pow= ', pow
       !write(File_E,*) 'pow is:', pow
      if (kflag .eq. kGAMDOT) then 
         SSKineticEqn = PlaP%gam0(PhID(ip)) * arg * pow
      else if (kflag .eq. kdGAMdTAU) then
         SSKineticEqn = PlaP%gam0(PhID(ip)) / (PlaP%xm(PhID(ip))*crss1) * pow
      else if (kflag .eq. kdGAMdKAPP) then
         SSKineticEqn = -PlaP%gam0(PhID(ip)) / (PlaP%xm(PhID(ip))*crss1) * arg * pow
!---- This is used to caculate the dGamDOTdTau for caculating Jacobian when sloving the first set of equations
      else 
         call RunTimeError(File_o, 'SSKineticEqn: unknown kflag!')
      endif

      return
      END
!
!=====================================================================72
!
!=====================================================================72
!
      logical FUNCTION NConverged( &
         res, toler, n &
         )
      
      use NumType

      implicit none


      real*8  toler
      real*8  res(NVEC)

      integer i,n
!
!---------------------------------------------------------------------72
!     
!------- check convergence on residual
!
      NConverged = .true.
      do i = 1, n
         NConverged = ( (dabs(res(i)) .lt. toler) .and. NConverged)
      enddo

      return
      END
!     
!=====================================================================72
!
!=====================================================================72
!
      SUBROUTINE BoundForArgPowLaw( xm, argmin, argmax &
        )

      use NumType      
      implicit none

      real*8  xm, xmm, argmin, argmax

      real*8  EXPON
      data EXPON /280.d0/
!
!---------------------------------------------------------------------72
!
!---- limits (bounds) on the value of the argument for power law
!---- note: In general: xm <= 1.0 (1/xm >= 1.0)
!
      xmm = pone/xm - pone
      if (dabs(xm-pone) .lt. TINY(1.d0)) xmm = pone

      argMin = dexp(-EXPON/(xmm)*dlog(10.d0))
      argMax = dexp( EXPON/(xmm)*dlog(10.d0))

      return
      END
!
!=====================================================================72
!=====================================================================72
!

    real*8 FUNCTION CheckArgPowLaw(   &
         arg, sgn, argmin, argmax     &
         )
      
      use NumType
 
      implicit none
      
      real*8  arg, sgn, argmin, argmax


!
!---------------------------------------------------------------------72
!
!---- check range of argument for power law
!
!write(File_E,*) 'argmin and argmax are:', argmin,argmax
      if (dabs(arg) .lt. argMin) then
         CheckArgPowLaw = argMin * sgn
      else if (dabs(arg) .ge. argMax) then
         CheckArgPowLaw = argMax * sgn
      else
         CheckArgPowLaw = arg
      endif

      return
      END
!=====================================================================72

