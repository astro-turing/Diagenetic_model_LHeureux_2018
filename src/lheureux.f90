!     Last change:  IL   11 Jan 2023
! 29 nov 2023, updated code to modern standards by
!              - turning common blocks into modules
!              - made the array soide of a 1000 more flexible
!              - update outdated PAUSE statements by read statements.

!  Rythmite limestone/marl - semiimplicit diffusion on c,po and explicit advection on CA,ARA with upwind scheme
module vecsize
   implicit none
   integer, parameter :: vsize = 1000
end module

module general
   use vecsize
   use, intrinsic :: iso_fortran_env, only: real64, int64
   implicit none
   real(kind=8) pho(0:vsize), cao(0:vsize), coo(0:vsize), ARAo(0:vsize), CALo(0:vsize)
   integer outx, N, outt
   integer(int64) :: tmax
end module

module par
   real(kind=8) :: P(35)
end module

!  Main code
program marl_PDE

   use vecsize
   use general
   use par

   implicit none

   real(kind=8) ph(0:vsize), ca(0:vsize), co(0:vsize), ARA(0:vsize), CAL(0:vsize), U(0:vsize), W(0:vsize)
   real(kind=8) phalf(0:vsize), ARAhalf(0:vsize), CALhalf(0:vsize), cahalf(0:vsize), cohalf(0:vsize), whalf(0:vsize)
   real(kind=8) Rca(0:vsize), Rco(0:vsize), RAR(0:vsize), RCAL(0:vsize), RARdum(0:vsize), RCALdum(0:vsize)
   real(kind=8) OC(0:vsize), OA(0:vsize), t, dt, dx, dca(0:vsize), dco(0:vsize), sigca(0:vsize), sigco(0:vsize)
   real(kind=8) Rp(0:vsize), te, S(0:vsize), Sdum(0:vsize), Udum(0:vsize), sigpo(0:vsize)

   integer i, j

   ! Open output files. 8-11: spatial profiles at four different times; OC, OA and Rp are oversaturations
   ! and porosity reaction reaction. 12: time series taken at four positions; U and W taken at bottom of the system
   ! 13: miscellaneous variables
   open (8, file='amarlt1')
   OPEN (9, FILE='amarlt2')
   OPEN (10, FILE='amarlt3')
   OPEN (11, FILE='amarlt4')
   open (12, FILE='amarlx')
   OPEN (13, FILE='marlstuff')

! Prepare the run. INIT gets all parameters and initial conditions
   call init
   dt = P(15)
   dx = P(16)

! ph=porosity, ca = calcium ion, co=carbonate ion, ARA=aragonite, CAL=calcite
   do i = 0, N
      ph(i) = pho(i)
      ca(i) = cao(i)
      co(i) = coo(i)
      ARA(i) = ARAo(i)
      CAL(i) = CALo(i)
   end do
   t = 0.0

   WRITE(12,101)'t ','AR','AR','AR','AR','CA','CA','CA','CA','P ','P ','P ','P ','ca','ca','ca','ca','co','co','co','co'&
               , 'U ', 'W ', 'OC', 'OA', 'RA', 'RC'
   WRITE (8, 103) 'x ', 'AR', 'CA', 'Po', 'Ca', 'CO', 'TE', 'U ', 'W ', 'OC', 'OA', 'Rp'
   WRITE (9, 103) 'x ', 'AR', 'CA', 'Po', 'Ca', 'CO', 'TE', 'U ', 'W ', 'OC', 'OA', 'Rp'
   WRITE (10, 103) 'x ', 'AR', 'CA', 'Po', 'Ca', 'CO', 'TE', 'U ', 'W ', 'OC', 'OA', 'Rp'
   WRITE (11, 103) 'x ', 'AR', 'CA', 'Po', 'Ca', 'CO', 'TE', 'U ', 'W ', 'OC', 'OA', 'Rp'
   write (13, 110) 't ', 'p1', 'k1', 'u1', 'p2', 'k2', 'u2', 'p3', 'k3', 'u3'

! auxf gets all auxiliairy functions: the velocities U, W, the grid Peclet numbers, the reaction rates,
! the diffusion coefficients
   call auxf(t, n, ph, ca, co, ARA, CAL, U, S, W, OC, OA, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, RAR, RCAL)

   write (12, 100) t, P(18), P(18), P(18), P(18), P(19), P(19), P(19), P(19), P(8), P(8), P(8), P(8), &
      ca(N/4), ca(N/2), ca(3*N/4), ca(N), co(N/4), co(N/2), co(3*N/4), co(N), U(N), W(N), &
      OC(N), OA(N), RAR(N), RCAL(N)

! projectX evaluates ARA, CAL, U at t=t+dt/2
   call projectX(n, ARA, CAL, U, S, RAR, RCAL, ARAhalf, CALhalf)

! projectX evaluates ca, co,ph, W at t=t+dt/2
   call projectY(n, ph, ca, co, W, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, phalf, cahalf, cohalf)

! Start the procedure
   do j = 1, tmax
      IF (j/50000*50000 .eq. j) WRITE (6, *) 'doing t=', j*dt
      call auxf(t, n, phalf, cahalf, cohalf, ARAhalf, CALhalf, Udum, Sdum, Whalf, OC, OA, dca, dco, sigpo,&
                      sigca, sigco, Rp, Rca, Rco, RARdum, RCALdum)

      ! CRANK uses the Crank-Nicholson and Fiadeiro-Veronis algorithm to solve the diffusive variables
      call CRANK(t, n, ph, ca, co, Whalf, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, phalf)

      ! upwind uses the upwind algorithm to solve the advective equations for CAL, ARA
      call upwind(n, ARA, CAL, U, S, RAR, RCAL)
      t = t + dt
      call auxf(t, n, ph, ca, co, ARA, CAL, U, S, W, OC, OA, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, RAR, RCAL)

      !  Generates time series at every outt
      IF (j/outt*outt .eq. j) then
         write (12, 100) t, ara(N/4), ara(N/2), ara(3*N/4), ara(N), cal(N/4), cal(N/2), cal(3*N/4), cal(N), &
            ph(N/4), ph(N/2), ph(3*N/4), ph(N), ca(N/4), ca(N/2), ca(3*N/4), ca(N), &
            co(N/4), co(N/2), co(3*N/4), co(N), U(N), W(N), OC(N), OA(N), RAR(N), RCAL(N)
      end if

      IF (j .eq. outx) then
         !  Output profiles at every multiple of outx
         do i = 0, N
            ! TE=terrigeneous concentration
            TE = 1.-ara(i) - cal(i)
            write (8, 102) i*dx, ara(i), cal(i), ph(i), ca(i), co(i), TE, U(i), W(i), OC(i), OA(i), Rp(i)
         end do
      end if

      IF (j .eq. 2*outx) then
         do i = 0, N
            TE = 1 - ara(i) - cal(i)
            write (9, 102) i*dx, ara(i), cal(i), ph(i), ca(i), co(i), TE, U(i), W(i), OC(i), OA(i), Rp(i)
         end do
      end if

      IF (j .eq. 3*outx) then
         do i = 0, N
            TE = 1 - ara(i) - cal(i)
            write (10, 102) i*dx, ara(i), cal(i), ph(i), ca(i), co(i), TE, U(i), W(i), OC(i), OA(i), Rp(i)
         end do
      end if

      IF (j .eq. 4*outx) then
         do i = 0, N
            TE = 1 - ara(i) - cal(i)
            write (11, 102) i*dx, ara(i), cal(i), ph(i), ca(i), co(i), TE, U(i), W(i), OC(i), OA(i), Rp(i)
         end do
      end if

      call projectX(n, ARA, CAL, U, S, RAR, RCAL, ARAhalf, CALhalf)
      call projectY(n, ph, ca, co, W, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, phalf, cahalf, cohalf)

   end do

100 format(27(d12.5, 1x))
101 format(13x, 27(a2, 11x))
102 format(12(d12.5, 1x))
103 format(13x, 12(a2, 11x))
110 format(13x, 10(a2, 11x))

   close (8)
   close (9)
   close (10)
   close (11)
   close (12)
   close (13)
   write (6, *) 'fini'

   ! Wait for an ackowledgment from the screen
   write (*, *) "pausing execution. Press enter to continue"
   read (*, *)
   stop 'marl'
end
!
!
!-----------------------------------------------------------------------------------------
!
subroutine auxf(t, n, ph, ca, co, ARA, CAL, U, S, W, OC, OA, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, RAR, RCAL)

   use par
   use vecsize

! This routine calculates all auxiliairy functions: Peclet numbers, oversaturations, reaction rates and others
   implicit none
   real(kind=8) ph(0:vsize), ca(0:vsize), co(0:vsize), ARA(0:vsize), CAL(0:vsize), U(0:vsize), W(0:vsize)
   real(kind=8) OC(0:vsize), OA(0:vsize), k(0:vsize), dca(0:vsize), dco(0:vsize), sigca(0:vsize), sigco(0:vsize)
   real(kind=8) Rp(0:vsize), betasV, S(0:vsize), sigpo(0:vsize)
   real(kind=8) Rca(0:vsize), Rco(0:vsize), RAR(0:vsize), RCAL(0:vsize)
   real(kind=8) Pe1, Pe2, Pe3, sc, sa, t, rhosw0
   integer idis, idisf
   integer i, n
   t = t
   rhosw0 = P(20)
! Node of the begnning of the ADZ
   idis = NINT(P(22)/P(16))
! Node of the end of the ADZ
   idisf = NINT((P(22) + P(14))/P(16))
   betasV = P(32)
   u(0) = P(23)
! Sign of the u velocity
   S(0) = 1.
!   k is  permeability
   k(0) = ph(0)**3*betasV/(1 - ph(0))**2
   k(0) = k(0)*(1 - dexp(-10.d0*(1 - ph(0))/ph(0)))
   w(0) = u(0) - k(0)*(rhosw0 - 1)*(1 - ph(0))/ph(0)
! Oversaturations at the boundary
   sa = 1 - ca(0)*co(0)*P(24)
   sc = ca(0)*co(0) - 1.
! ARAGONITE:

   ! Inside the ADZ
   IF ((0 .gt. idis) .and. (0 .lt. idisf)) then
      if (sa .gt. 0.) then
         OA(0) = sa**P(1)
      else
         OA(0) = -(-sa)**P(1)
      end if
   end if

   ! Just at the ADZ boundaries. This special case allows the possibility to give
   ! a different weight to the contribution at the boundaries. Here, the weight is 1
   IF ((0 .eq. idis) .or. (0 .eq. idisf)) then
      if (sa .gt. 0.) then
         OA(0) = sa**P(1)
      else
         OA(0) = -(-sa)**P(1)
      end if
   end if

   ! Outside the ADZ
   IF ((0 .lt. idis) .or. (0 .gt. idisf)) then
      if (sa .gt. 0.) then
         OA(0) = 0.d0
      else
         OA(0) = -(-sa)**P(1)
      end if
   end if

! Calcite:
   if (sc .gt. 0.) then
      OC(0) = sc**P(2)
   else
      OC(0) = -(-sc)**P(2)
   end if

! Diffusion coefficients corrected for porosity
   dca(0) = 1./(1 - 2*dlog(ph(0)))
   dco(0) = dca(0)*P(25)

! Reaction rates
   Rca(0) = P(26)*(P(6)/P(7) - ca(0))*(OA(0)*ARA(0) - P(3)*OC(0)*CAL(0))*(1 - ph(0))/ph(0)
   Rco(0) = P(26)*(P(6)/P(7) - co(0))*(OA(0)*ARA(0) - P(3)*OC(0)*CAL(0))*(1 - ph(0))/ph(0)

! Grid Peclet numbers and sigma coefficients in the Fiaideiro/Veronis scheme
   Pe1 = w(0)*P(16)/(2*dca(0))
   Pe2 = w(0)*P(16)/(2*dco(0))
   Pe3 = w(0)*P(16)/(2*P(35))

   if (Pe1 .gt. 0.) then
      sigca(0) = (1.+dexp(-2*Pe1))/(1.-dexp(-2*Pe1)) - 1./Pe1
   else
      sigca(0) = (1.+dexp(2*Pe1))/(-1.+dexp(2*Pe1)) - 1./Pe1
   end if

   if (Pe2 .gt. 0.) then
      sigco(0) = (1.+dexp(-2*Pe2))/(1.-dexp(-2*Pe2)) - 1./Pe2
   else
      sigco(0) = (1.+dexp(2*Pe2))/(-1.+dexp(2*Pe2)) - 1./Pe2
   end if

   if (Pe3 .gt. 0.) then
      sigpo(0) = (1.+dexp(-2*Pe3))/(1.-dexp(-2*Pe3)) - 1./Pe3
   else
      sigpo(0) = (1.+dexp(2*Pe3))/(-1.+dexp(2*Pe3)) - 1./Pe3
   end if

   !  other reaction rates
   ! No terrigeneous material
   if (1 - ARA(0) - CAL(0) .lt. 1.d-70) then
      RAR(0) = -P(26)*OA(0)*(1 - CAL(0))*CAL(0) - P(26)*P(3)*OC(0)*(1 - CAL(0))*CAL(0)
      RCAL(0) = P(3)*P(26)*OC(0)*CAL(0)*(1 - CAL(0)) + P(26)*OA(0)*CAL(0)*(1 - CAL(0))
      Rp(0) = P(26)*((1 - CAL(0))*OA(0) - P(3)*CAL(0)*OC(0))*(1 - ph(0))
   else
      ! With terrigeneous material
      RAR(0) = -P(26)*OA(0)*ARA(0)*(1 - ARA(0)) - P(26)*P(3)*OC(0)*ARA(0)*CAL(0)
      RCAL(0) = P(3)*P(26)*OC(0)*CAL(0)*(1 - CAL(0)) + P(26)*OA(0)*CAL(0)*ARA(0)
      Rp(0) = P(26)*(ARA(0)*OA(0) - P(3)*CAL(0)*OC(0))*(1 - ph(0))
   end if

   !Oversaturations at other nodes
   do i = 1, N
      sa = 1 - ca(i)*co(i)*P(24)
      sc = ca(i)*co(i) - 1

! Aragonite:
      IF ((i .gt. idis) .and. (i .lt. idisf)) then
         if (sa .gt. 0.) then
            OA(i) = sa**P(1)
         else
            OA(i) = -(-sa)**P(1)
         end if
      end if

      IF ((i .eq. idis) .or. (i .eq. idisf)) then
         if (sa .gt. 0.) then
            OA(i) = sa**P(1)
         else
            OA(i) = -(-sa)**P(1)
         end if
      end if

      IF ((i .lt. idis) .or. (i .gt. idisf)) then
         if (sa .gt. 0.) then
            OA(i) = 0.d0
         else
            OA(i) = -(-sa)**P(1)
         end if
      end if

      ! Calcite:
      if (sc .gt. 0.) then
         OC(i) = sc**P(2)
      else
         OC(i) = -(-sc)**P(2)
      end if
! velocities, Peclet numbers, diffusion coefficients and reaction rates
! In case the porosity is close to one
      if ((1 - ph(i)) .le. 0.05) then
         ! write(6,*) 'p>0.95 for t,i',t,i
         ! permeability
         k(i) = betasV*10.*ph(i)**2/(1 - ph(i))
         u(i) = k(i)*(1 - ph(i))*(rhosw0 - 1) + P(23) - P(13)*k(0)*(1 - ph(0))
         w(i) = u(i) - k(i)*(rhosw0 - 1)*(1 - ph(i))/ph(i)
         dca(i) = 1./(1 - 2*dlog(ph(i)))
         dco(i) = dca(i)*P(25)
         Rca(i) = P(26)*(P(6)/P(7) - ca(i))*(OA(i)*ARA(i) - P(3)*OC(i)*CAL(i))*(1 - ph(i))/ph(i)
         Rco(i) = P(26)*(P(6)/P(7) - co(i))*(OA(i)*ARA(i) - P(3)*OC(i)*CAL(i))*(1 - ph(i))/ph(i)
         Pe1 = w(i)*P(16)/(2*dca(i))
         Pe2 = w(i)*P(16)/(2*dco(i))
         Pe3 = w(i)*P(16)/(2*P(35))

         if (Pe1 .gt. 0.) then
            sigca(i) = (1.+dexp(-2*Pe1))/(1.-dexp(-2*Pe1)) - 1./Pe1
         else
            sigca(i) = (1.+dexp(2*Pe1))/(-1.+dexp(2*Pe1)) - 1./Pe1
         end if

         if (Pe2 .gt. 0.) then
            sigco(i) = (1.+dexp(-2*Pe2))/(1.-dexp(-2*Pe2)) - 1./Pe2
         else
            sigco(i) = (1.+dexp(2*Pe2))/(-1.+dexp(2*Pe2)) - 1./Pe2
         end if

         If (Pe3 .gt. 0.) then
            sigpo(i) = (1.+dexp(-2*Pe3))/(1.-dexp(-2*Pe3)) - 1./Pe3
         else
            sigpo(i) = (1.+dexp(2*Pe3))/(-1.+dexp(2*Pe3)) - 1./Pe3
         end if

      else
         ! Special case if the porosity is nearly 0
         if (ph(i) .le. P(34)) THEN
            write (6, *) 'phi = eps at i, t=', ph(i), i, t
            k(i) = 0.
            dca(i) = 0.
            dco(i) = 0.
            Rca(i) = 0.
            Rco(i) = 0.
            u(i) = k(i)*(rhosw0 - 1)*(1 - ph(i)) + P(23) - P(13)*k(0)*(1 - ph(0))
            ! Fluid flows with the solid in this case
            w(i) = u(i)

            if (w(i) .gt. 0.) then
               sigca(i) = 1.
               sigco(i) = 1.
            else
               sigca(i) = -1.
               sigco(i) = -1.
            end if

            if (w(i) .gt. 0.) then
               sigpo(i) = 1.
            else
               sigpo(i) = -1.
            end if
         else
            ! Normal case
            ! Permeability
            k(i) = betasV*ph(i)**3/(1 - ph(i))**2
            k(i) = k(i)*(1 - dexp(-10.*(1 - ph(i))/ph(i)))
            ! diffusion coefficient
            dca(i) = 1./(1 - 2*dlog(ph(i)))
            dco(i) = dca(i)*P(25)
            ! reaction rates
            Rca(i) = P(26)*(P(6)/P(7) - ca(i))*(OA(i)*ARA(i) - P(3)*OC(i)*CAL(i))*(1 - ph(i))/ph(i)
            Rco(i) = P(26)*(P(6)/P(7) - co(i))*(OA(i)*ARA(i) - P(3)*OC(i)*CAL(i))*(1 - ph(i))/ph(i)
!  velocities, Peclet numbers and sigma factors
            u(i) = k(i)*(rhosw0 - 1)*(1 - ph(i)) + P(23) - P(13)*k(0)*(1 - ph(0))
            w(i) = u(i) - k(i)*(rhosw0 - 1)*(1 - ph(i))/ph(i)
            Pe1 = w(i)*P(16)/(2*dca(i))
            Pe2 = w(i)*P(16)/(2*dco(i))
            Pe3 = w(i)*P(16)/(2*P(35))

            if (Pe1 .gt. 0.) then
               sigca(i) = (1.+dexp(-2*Pe1))/(1.-dexp(-2*Pe1)) - 1./Pe1
            else
               sigca(i) = (1.+dexp(2*Pe1))/(-1.+dexp(2*Pe1)) - 1./Pe1
            end if

            if (Pe2 .gt. 0.) then
               sigco(i) = (1.+dexp(-2*Pe2))/(1.-dexp(-2*Pe2)) - 1./Pe2
            else
               sigco(i) = (1.+dexp(2*Pe2))/(-1.+dexp(2*Pe2)) - 1./Pe2
            end if

            if (Pe3 .gt. 0.) then
               sigpo(i) = (1.+dexp(-2*Pe3))/(1.-dexp(-2*Pe3)) - 1./Pe3
            else
               sigpo(i) = (1.+dexp(2*Pe3))/(-1.+dexp(2*Pe3)) - 1./Pe3
            end if
         end if
      end if

      !  other reaction rates
      if (1 - ARA(i) - CAL(i) .lt. 1.d-70) then
         RAR(i) = -P(26)*OA(i)*(1 - CAL(i))*CAL(i) - P(26)*P(3)*OC(i)*(1 - CAL(i))*CAL(i)
         RCAL(i) = P(3)*P(26)*OC(i)*CAL(i)*(1 - CAL(i)) + P(26)*OA(i)*CAL(i)*(1 - CAL(i))
         Rp(i) = P(26)*((1 - CAL(i))*OA(i) - P(3)*CAL(i)*OC(i))*(1 - ph(i))
      else
         RAR(i) = -P(26)*OA(i)*ARA(i)*(1 - ARA(i)) - P(26)*P(3)*OC(i)*ARA(i)*CAL(i) ! ALways negative  -> solid to liquid
         RCAL(i) = P(3)*P(26)*OC(i)*CAL(i)*(1 - CAL(i)) + P(26)*OA(i)*CAL(i)*ARA(i)
         Rp(i) = P(26)*(ARA(i)*OA(i) - P(3)*CAL(i)*OC(i))*(1 - ph(i))
      end if
! Sign of U
      S(i) = u(i)/dabs(u(i))
   end do
   return
end
!
!-----------------------------------------------------------------------------------------
!

subroutine projectX(n, ARA, CAL, U, S, RAR, RCAL, ARAhalf, CALhalf)

   use par
   use vecsize

   implicit none

! calculates the projected solid phase fields at half-time in order to treat the non-linearities
! see Rosenberg p. 58 as cited in L'Heureux (2018)
   real(kind=8) U(0:vsize), ARA(0:vsize), CAL(0:vsize), ARAhalf(0:vsize), CALhalf(0:vsize)
   real(kind=8) S(0:vsize)
   real(kind=8) dt, dx, a
   real(kind=8) RAR(0:vsize), RCAL(0:vsize)
   integer N, i

! n: grid point index, where 0 = surface, n = bottom
! U: solid phase velocity
! ARA, CAL: aragonite, calcite
! RAR: Reactino rate aragonite
! RCal: Reactino rate calcite
! S: sign of the solid phase velocity

   dx = P(16)
   dt = P(15)
   a = dt/(2*dx) ! coefficient used in upwind

! SURFACE
! boundary conditions at surface: Aragonite and Calcite value remain unchanged (and thus are fixed to the values set in init)
! Eq. 35 fro L'Heureux (2018)
   ARAhalf(0) = ARA(0)
   CALhalf(0) = CAL(0)

! INTERIOR
   do i = 1, n - 1
      ARAhalf(i) = ARA(i) - a*u(i)*(ARA(i)*S(i) - ARA(i - 1)*0.5*(S(i) + 1.) + ARA(i + 1)*0.5*(1.-S(i)))&
              + 0.5*dt*RAR(i)
      CALhalf(i) = CAL(i) - a*u(i)*(CAL(i)*S(i) - CAL(i - 1)*0.5*(S(i) + 1.) + CAL(i + 1)*0.5*(1.-S(i)))&
               + 0.5*dt*RCAL(i)
      if (1 - ARAhalf(i) - CALhalf(i) .lt. 1.d-70) ARAhalf(i) = 1 - CALhalf(i)
      if (1 - CALhalf(i) .lt. 1.d-10) CALhalf(i) = 1.
      if (1 - ARAhalf(i) .lt. 1.d-10) ARAhalf(i) = 1.
      if (ARAhalf(i) .lt. 1.d-70) ARAhalf(i) = 0.
      if (CALhalf(i) .lt. 1.d-70) CALhalf(i) = 0.
   end do

! BOTTOM
   ARAhalf(n) = ARA(n) - a*u(n)*S(n)*(ARA(n) - ARA(n - 1)) + 0.5*dt*RAR(n)
   CALhalf(n) = CAL(n) - a*u(n)*S(n)*(CAL(n) - CAL(n - 1)) + 0.5*dt*RCAL(n)
   if (1 - ARAhalf(n) - CALhalf(n) .lt. 1.d-70) ARAhalf(n) = 1 - CALhalf(n)
   if (1 - CALhalf(n) .lt. 1.d-10) CALhalf(n) = 1.
   if (1 - ARAhalf(n) .lt. 1.d-10) ARAhalf(n) = 1.
   if (ARAhalf(n) .lt. 1.d-70) ARAhalf(n) = 0.
   if (CALhalf(n) .lt. 1.d-70) CALhalf(n) = 0.
   return
end
!
!-----------------------------------------------------------------------------------------
!
SUBROUTINE projectY(n, ph, ca, co, W, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, phalf, cahalf, cohalf)

   use par
   use vecsize

   implicit none

! This routine calculates the projected aqueous phase fields at half-time in order to treat the non-linearities
! see Rosenberg p. 58 as cited in L'Heureux (2018)
   real(kind=8) ph(0:vsize), W(0:vsize), cahalf(0:vsize), cohalf(0:vsize), phalf(0:vsize), sigpo(0:vsize)
   real(kind=8) ca(0:vsize), co(0:vsize), dca(0:vsize), dco(0:vsize), sigca(0:vsize), sigco(0:vsize), Rca(0:vsize), Rco(0:vsize)
   real(kind=8) dt, dx, a, b, eps, difpor, Rp(0:vsize)
   integer N, i
! N: number of grid points
! ph: porosity
! ca: dissolved Ca
! co: dissolved Carbonate
! W: velocity of liquid phase
! dca
! dco
! sigca: peclet no of ca
! sigco: peclet no of co
! Rp: reaction rate of ph (porosity)
! Rca: reaction rate of dissolved ca
! Rco: reaction rate of dissolved carbonate
! phalf: projected porosity at half time
! cahalf: projected dissolved ca at half time
! cohalf: projected dissolved carbonate at half time

! depth and time incrementrs
   dx = P(16)
   dt = P(15)
! small value to cap porosity
   eps = P(29)
! coefficients used in the projection method
   a = dt/(4*dx)
   b = dt/(4*dx*dx)
! Constant porosity
   difpor = P(35)
! SURFACE
! projected ca, co, & porosity at surface: use constant in accordance of initial conditions (L'Heureux (2018), eq. 35)
   cahalf(0) = ca(0)
   cohalf(0) = co(0)
   phalf(0) = ph(0)

! INTERIOR
   do i = 1, n - 1
! Project porosity according to Rosenberg p. 58 as cited in L'Heureux (2018)
! uses two peclet numbers here. WHY?
      phalf(i) = ph(i) - a*((1 - sigpo(i + 1))*ph(i + 1)*w(i + 1) + 2*sigpo(i)*ph(i)*w(i) -&
             (1 + sigpo(i - 1))*ph(i - 1)*w(i - 1))&
             + 2*b*difpor*(ph(i - 1) - 2*ph(i) + ph(i + 1)) + 0.5*dt*Rp(i)
! Cap projected porosity: minimum projected value is eps
      if (phalf(i) .lt. eps) phalf(i) = eps
! Cap projected porosity: maximum projected value is 1-eps
      if (1 - phalf(i) .lt. eps) phalf(i) = 1.-eps
! Case: Porosity is small -> neglect diffusion term. Use spatial discretisation vie FV
      if (ph(i) .le. eps) then
! Project ca & co
         cahalf(i) = ca(i) - a*w(i)*((1 - sigca(i))*ca(i + 1) + 2*sigca(i)*ca(i) - (1 + sigca(i))*ca(i - 1)) &
               + 0.5*dt*Rca(i)
         cohalf(i) = co(i) - a*w(i)*((1 - sigco(i))*co(i + 1) + 2*sigco(i)*co(i) - (1 + sigco(i))*co(i - 1)) &
               + 0.5*dt*Rco(i)
      else
! Case: Porosity is not small -> include diffusion with central differences
         cahalf(i) = ca(i) - a*w(i)*((1 - sigca(i))*ca(i + 1) + 2*sigca(i)*ca(i) - (1 + sigca(i))*ca(i - 1)) &
               + b*(ph(i + 1)*dca(i + 1) + ph(i)*dca(i))*(ca(i + 1) - ca(i))/ph(i) &
               - b*(ph(i - 1)*dca(i - 1) + ph(i)*dca(i))*(ca(i) - ca(i - 1))/ph(i) + 0.5*dt*Rca(i)
         cohalf(i) = co(i) - a*w(i)*((1 - sigco(i))*co(i + 1) + 2*sigco(i)*co(i) - (1 + sigco(i))*co(i - 1)) &
               + b*(ph(i + 1)*dco(i + 1) + ph(i)*dco(i))*(co(i + 1) - co(i))/ph(i) &
               - b*(ph(i - 1)*dco(i - 1) + ph(i)*dco(i))*(co(i) - co(i - 1))/ph(i) + 0.5*dt*Rco(i)
      end if
! Case: if projected ca & co are small, set to 0
      if (cahalf(i) .lt. 1.d-15) cahalf(i) = 0.
      if (cohalf(i) .lt. 1.d-15) cohalf(i) = 0.
   end do

! BOTTOM
! What is happening here?
! Uses more than one peclet number. Why?
! Why is there only sigpo, and not (1-sigpo) as in Boudreax (1996), eq. 96
   phalf(n) = ph(n) + 2*a*(sigpo(n - 1)*ph(n - 1)*w(n - 1) - sigpo(n)*ph(n)*w(n)) &
             + 4*b*difpor*(ph(n - 1) - ph(n)) + 0.5*dt*Rp(n)
! Cap projected porosity: minimum projected value is eps
   if (phalf(n) .lt. eps) phalf(n) = eps
! Cap projected porosity: maximum projected value is 1-eps
   if (1 - phalf(n) .lt. eps) phalf(n) = 1.-eps
! Case: porosity is small: use F-V
   if (ph(n) .le. eps) then
      cahalf(n) = ca(n) - 2*a*w(n)*sigca(n)*(ca(n) - ca(n - 1)) + 0.5*dt*Rca(n)
      cohalf(n) = co(n) - 2*a*w(n)*sigco(n)*(co(n) - co(n - 1)) + 0.5*dt*Rco(n)
! Case: porosity is not small -> use FV and central diff for diffusion
   else
      cahalf(n) = ca(n) - 2*a*w(n)*sigca(n)*(ca(n) - ca(n - 1)) &
               + 2*b*(ph(n - 1)*dca(n - 1) + ph(n)*dca(n))*(ca(n - 1) - ca(n))/ph(n) + 0.5*dt*Rca(n)
      cohalf(n) = co(n) - 2*a*w(n)*sigco(n)*(co(n) - co(n - 1)) &
               + 2*b*(ph(n - 1)*dco(n - 1) + ph(n)*dco(n))*(co(n - 1) - co(n))/ph(n) + 0.5*dt*Rco(n)
   end if
! Cap ca & co if they are too small
   if (cahalf(n) .lt. 1.d-15) cahalf(n) = 0.
   if (cohalf(n) .lt. 1.d-15) cohalf(n) = 0.
   return
end
!
!-----------------------------------------------------------------------------------------

subroutine upwind(n, ARA, CAL, U, S, RAR, RCAL)

   use par
   use vecsize

   implicit none

! Implementation of the upwind scheme for the advection equations (Calcite, Aragonite)
   real(kind=8) U(0:vsize), ARA(0:vsize), CAL(0:vsize)
   real(kind=8) RAR(0:vsize), RCAL(0:vsize), ARAnew(0:vsize), CALnew(0:vsize)
   real(kind=8) S(0:vsize)
   real(kind=8) dt, dx, a
   integer N, i
! U: solid phase velocity
! Ara, Cal aragonite, calcite
! RAR: Reactino rate aragonite
! RCal: Reactino rate calcite
! S: sign of the solid phase velocity
   dx = P(16)
   dt = P(15)
   a = dt/dx
! SURFACE
! boundary conditions at surface: Aragonite and Calcite value remain unchanged (and thus are fixed to the values set in init)
! Eq. 35 fro L'Heureux (2018)
   ARAnew(0) = ARA(0)
   CALnew(0) = CAL(0)

! INTERIOR
   do i = 1, n - 1
! Analytical expression for upwind that automatically switches directions dependent on the sign of U (stored in variable S)
      ARANEW(i) = ARA(i) - a*u(i)*(ARA(i)*S(i) - ARA(i - 1)*0.5*(S(i) + 1.) + ARA(i + 1)*0.5*(1.-S(i)))&
&              + dt*RAR(i)
      CALNEW(i) = CAL(i) - a*u(i)*(CAL(i)*S(i) - CAL(i - 1)*0.5*(S(i) + 1.) + CAL(i + 1)*0.5*(1.-S(i)))&
&               + dt*RCAL(i)
! Case: No terestrial material -> calcite and aragonite add to 1
      if (1 - ARAnew(i) - CALnew(i) .lt. 1.d-70) ARAnew(i) = 1 - CALnew(i)
! Case: Calcite dominates
      if (1 - CALnew(i) .lt. 1.d-10) CALnew(i) = 1.
! Case: Aragonite dominates
      if (1 - ARAnew(i) .lt. 1.d-10) ARAnew(i) = 1.
! Case: Aragonite vanishes
      if (ARAnew(i) .lt. 1.d-70) ARAnew(i) = 0.
! Case: Calcite vanishes
      if (CALnew(i) .lt. 1.d-70) CALnew(i) = 0.
   end do

! BOTTOM
! Use upwind with velocity |u(n)| (= u(n)*S(n)) for Aragonite and Calcite
! Not sure this makes sense ??? -Nick
   ARAnew(n) = ARA(n) - a*u(n)*S(n)*(ARA(n) - ARA(n - 1)) + dt*RAR(n)
   CALnew(n) = CAL(n) - a*u(n)*S(n)*(CAL(n) - CAL(n - 1)) + dt*RCAL(n)
! Case: No terrestrial material
   if (1 - ARAnew(n) - CALnew(n) .lt. 1.d-70) ARAnew(n) = 1 - CALnew(n)
! Cases: Aragonite resp. Calcite dominate
   if (1 - CALnew(n) .lt. 1.d-10) CALnew(n) = 1.
   if (1 - ARAnew(n) .lt. 1.d-10) ARAnew(n) = 1.
! Cases: Aragonite resp. Calcite vanish
   if (ARAnew(n) .lt. 1.d-70) ARAnew(n) = 0.
   if (CALnew(n) .lt. 1.d-70) CALnew(n) = 0.
! update Values
   do i = 0, n
      ARA(i) = ARAnew(i)
      CAL(i) = CALnew(i)
   end do
   return
end
!
!-----------------------------------------------------------------------------------------
!
subroutine CRANK(t, n, ph, ca, co, Whalf, dca, dco, sigpo, sigca, sigco, Rp, Rca, Rco, phalf)

   use par
   use vecsize

   implicit none

! implementation of the Crank-Nicholson scheme for the advection-diffusion equations: X(n+1)=(A^-1).B.X(n)
   real(kind=8) ph(0:vsize), ca(0:vsize), co(0:vsize), Whalf(0:vsize), sigpo(0:vsize), phalf(0:vsize)
   real(kind=8) dt, dx, eps, sigca(0:vsize), sigco(0:vsize), dca(0:vsize), dco(0:vsize), Rca(0:vsize), Rco(0:vsize)
   real(kind=8) aA1(0:vsize), bA1(0:vsize), cA1(0:vsize), aA2(0:vsize), bA2(0:vsize), cA2(0:vsize), Byca(0:vsize), Byco(0:vsize)
   real(kind=8) aA3(0:vsize), bA3(0:vsize), cA3(0:vsize), Bypo(0:vsize), Rp(0:vsize)
   real(kind=8) a, b, t
   integer N, i

   dx = P(16)
   dt = P(15)
   eps = P(29)
   a = dt/(4*dx)
   b = dt/(4*dx*dx)
   call matA(n, a, b, phalf, whalf, sigpo, sigca, sigco, dca, dco, aA1, bA1, cA1, aA2, bA2, cA2, aA3, bA3, cA3)
   call matBy(n, dt, a, b, ph, ca, co, phalf, whalf, sigpo, sigca, sigco, dca, dco, Rp, Rca, Rco, Bypo, Byca, Byco)
!  find solution
   call tridag(t, n, aA1, bA1, cA1, Byca, ca)
   call tridag(t, n, aA2, bA2, cA2, Byco, co)
   call tridag(t, n, aA3, bA3, cA3, Bypo, ph)
   do i = 0, n
      IF (ca(i) .lt. 1.d-15) ca(i) = 0.
      IF (co(i) .lt. 1.d-15) co(i) = 0.
      if (ph(i) .lt. eps) ph(i) = eps
      if (1 - ph(i) .lt. eps) ph(i) = 1.-eps
   end do
   return
end
!
!-----------------------------------------------------------------------------------------
!
!
subroutine matA(n, a, b, phalf, whalf, sigpo, sigca, sigco, dca, dco, aA1, bA1, cA1, aA2, bA2, cA2, aA3, bA3, cA3)

   use par
   use vecsize

   implicit none

   ! Determine the matrix elements of the tridiagonal matrix A in the Crank-Nicholson scheme
   real(kind=8) Whalf(0:vsize), phalf(0:vsize)
   real(kind=8) sigpo(0:vsize), sigca(0:vsize), sigco(0:vsize), dca(0:vsize), dco(0:vsize)
   real(kind=8) aA1(0:vsize), bA1(0:vsize), cA1(0:vsize), aA2(0:vsize), bA2(0:vsize), cA2(0:vsize)
   Real(kind=8) aA3(0:vsize), bA3(0:vsize), cA3(0:vsize)
   real(kind=8) a, b, eps, difpor
   integer N, i
   eps = P(29)
   difpor = P(35)
   aA1(0) = 0.
   bA1(0) = 1.
   cA1(0) = 0.
   aA2(0) = 0.
   bA2(0) = 1.
   cA2(0) = 0.
   aA3(0) = 0.
   bA3(0) = 1.
   cA3(0) = 0.
   do i = 1, n - 1
      aA3(i) = -a*whalf(i - 1)*(1 + sigpo(i - 1)) - 2*difpor*b
      bA3(i) = 1.+2*a*whalf(i)*sigpo(i) + 4*difpor*b
      cA3(i) = a*whalf(i + 1)*(1 - sigpo(i + 1)) - 2*difpor*b
      if (phalf(i) .le. eps) then
         aA1(i) = -a*whalf(i)*(1 + sigca(i))
         bA1(i) = 1.+a*2*whalf(i)*sigca(i)
         cA1(i) = a*whalf(i)*(1 - sigca(i))
         aA2(i) = -a*whalf(i)*(1 + sigco(i))
         bA2(i) = 1.+a*2*whalf(i)*sigco(i)
         cA2(i) = a*whalf(i)*(1 - sigco(i))
      else
         aA1(i) = -a*whalf(i)*(1 + sigca(i)) - b*(dca(i - 1)*phalf(i - 1) + dca(i)*phalf(i))/phalf(i)
         bA1(i) = 1.+a*2*whalf(i)*sigca(i) + b*(2*phalf(i)*dca(i) + phalf(i + 1)*dca(i + 1) + phalf(i - 1)*dca(i - 1))/phalf(i)
         cA1(i) = a*whalf(i)*(1 - sigca(i)) - b*(dca(i + 1)*phalf(i + 1) + dca(i)*phalf(i))/phalf(i)
         aA2(i) = -a*whalf(i)*(1 + sigco(i)) - b*(dco(i - 1)*phalf(i - 1) + dco(i)*phalf(i))/phalf(i)
         bA2(i) = 1.+a*2*whalf(i)*sigco(i) + b*(2*phalf(i)*dco(i) + phalf(i + 1)*dco(i + 1) + phalf(i - 1)*dco(i - 1))/phalf(i)
         cA2(i) = a*whalf(i)*(1 - sigco(i)) - b*(dco(i + 1)*phalf(i + 1) + dco(i)*phalf(i))/phalf(i)
      end if
   end do
   aA3(n) = -2*a*whalf(n - 1)*sigpo(n - 1) - 4*difpor*b
   bA3(n) = 1.+2*a*whalf(n)*sigpo(n) + 4*difpor*b
   cA3(n) = 1.
   if (phalf(n) .le. eps) then
      aA1(n) = -2*a*whalf(n)*sigca(n)
      bA1(n) = 1.+2*a*whalf(n)*sigca(n)
      cA1(n) = 1.
      aA2(n) = -2*a*whalf(n)*sigco(n)
      bA2(n) = 1.+2*a*whalf(n)*sigco(n)
      cA2(n) = 1.
   else
      aA1(n) = -2*a*whalf(n)*sigca(n) - b*(dca(n - 1)*phalf(n - 1) + 3*dca(n)*phalf(n))/phalf(n)
      bA1(n) = 1.+2*a*whalf(n)*sigca(n) + b*(3*dca(n)*phalf(n) + phalf(n - 1)*dca(n - 1))/phalf(n)
      cA1(n) = 1.
      aA2(n) = -2*a*whalf(n)*sigco(n) - b*(dco(n - 1)*phalf(n - 1) + 3*dco(n)*phalf(n))/phalf(n)
      bA2(n) = 1.+2*a*whalf(n)*sigco(n) + b*(3*dco(n)*phalf(n) + phalf(n - 1)*dco(n - 1))/phalf(n)
      cA2(n) = 1.
   end if
   return
end
!
!-----------------------------------------------------------------------------------------
!
subroutine matby(n, dt, a, b, ph, ca, co, phalf, whalf, sigpo, sigca, sigco, dca, dco, Rp, Rca, Rco, Bypo, Byca, Byco)

   use par
   use vecsize

   implicit none

! Determine the matrix elements of the tridiagonal matrix B in the Crank-Nicholson scheme
   real(kind=8) phalf(0:vsize), whalf(0:vsize), sigpo(0:vsize), sigca(0:vsize), sigco(0:vsize), dca(0:vsize)
   real(kind=8) dco(0:vsize), Rca(0:vsize), Rco(0:vsize), Rp(0:vsize)
   real(kind=8) Byca(0:vsize), Byco(0:vsize), Bypo(0:vsize), co(0:vsize), ca(0:vsize), ph(0:vsize)
   real(kind=8) a, b, dt, eps, difpor
   integer N, i

   eps = P(29)
   difpor = P(35)
   Byca(0) = ca(0)
   Byco(0) = co(0)
   Bypo(0) = ph(0)
   do i = 1, n - 1
      Bypo(i) = (a*whalf(i - 1)*(1 + sigpo(i - 1)) + 2*b*difpor)*ph(i - 1) &
                + (1.-a*2*whalf(i)*sigpo(i) - 4*b*difpor)*ph(i) &
                + (-a*whalf(i + 1)*(1 - sigpo(i + 1)) + 2*b*difpor)*ph(i + 1) + dt*Rp(i)
      if (phalf(i) .le. eps) then
         Byca(i) = (a*whalf(i)*(1 + sigca(i)))*ca(i - 1) &
                   + (1.-a*2*whalf(i)*sigca(i))*ca(i) &
                   + (-a*whalf(i)*(1 - sigca(i)))*ca(i + 1) + dt*Rca(i)
         Byco(i) = (a*whalf(i)*(1 + sigco(i)))*co(i - 1) &
                   + (1.-a*2*whalf(i)*sigco(i))*co(i) &
                   + (-a*whalf(i)*(1 - sigco(i)))*co(i + 1) + dt*Rco(i)
      else
         Byca(i) = (a*whalf(i)*(1 + sigca(i)) + b*(dca(i - 1)*phalf(i - 1) + dca(i)*phalf(i))/phalf(i))*ca(i - 1) &
             +(1.-a*2*whalf(i)*sigca(i)-b*(2*dca(i)*phalf(i)+dca(i+1)*phalf(i+1)+dca(i-1)*phalf(i-1))/phalf(i))*ca(i) &
             + (-a*whalf(i)*(1 - sigca(i)) + b*(dca(i + 1)*phalf(i + 1) + dca(i)*phalf(i))/phalf(i))*ca(i + 1) + dt*Rca(i)
         Byco(i) = (a*whalf(i)*(1 + sigco(i)) + b*(dco(i - 1)*phalf(i - 1) + dco(i)*phalf(i))/phalf(i))*co(i - 1) &
             +(1.-a*2*whalf(i)*sigco(i)-b*(2*dco(i)*phalf(i)+dco(i+1)*phalf(i+1)+dco(i-1)*phalf(i-1))/phalf(i))*co(i) &
             + (-a*whalf(i)*(1 - sigco(i)) + b*(dco(i + 1)*phalf(i + 1) + dco(i)*phalf(i))/phalf(i))*co(i + 1) + dt*Rco(i)
      end if
   end do
   Bypo(n) = (2*a*whalf(n - 1)*sigpo(n - 1) + 4*b*difpor)*ph(n - 1) &
             + (1.-a*2*whalf(n)*sigpo(n) - 4*b*difpor)*ph(n) + dt*Rp(n)
   if (phalf(n) .le. eps) then
      Byca(n) = 2*a*whalf(n)*sigca(n)*ca(n - 1) &
             + (1.-a*2*whalf(n)*sigca(n))*ca(n) + dt*Rca(n)
      Byco(n) = 2*a*whalf(n)*sigco(n)*co(n - 1) &
             + (1.-a*2*whalf(n)*sigco(n))*co(n) + dt*Rco(n)
   else
      Byca(n) = (2*a*whalf(n)*sigca(n) + b*(dca(n - 1)*phalf(n - 1) + 3*dca(n)*phalf(n))/phalf(n))*ca(n - 1) &
             + (1.-a*2*whalf(n)*sigca(n) - b*(3*dca(n)*phalf(n) + dca(n - 1)*phalf(n - 1))/phalf(n))*ca(n) + dt*Rca(n)
      Byco(n) = (2*a*whalf(n)*sigco(n) + b*(dco(n - 1)*phalf(n - 1) + 3*dco(n)*phalf(n))/phalf(n))*co(n - 1) &
             + (1.-a*2*whalf(n)*sigco(n) - b*(3*dco(n)*phalf(n) + dco(n - 1)*phalf(n - 1))/phalf(n))*co(n) + dt*Rco(n)
   end if
   return
end
!
!-----------------------------------------------------------------------------------------
!
subroutine tridag(t, m, a, b, c, r, sol)

   use vecsize

! Find the inverse of the tridiagonal matrix A by Thomas algorithm
   implicit none
   real(kind=8) a(0:vsize), b(0:vsize), c(0:vsize), r(0:vsize), sol(0:vsize)
   real(kind=8) h(0:vsize), g(0:vsize), denom, t
   INTEGER m, i
! t: time, used for debugging only
! m: no of grid points
! a,b,c,r,sol
   g(m - 1) = -a(m)/b(m)
   h(m - 1) = r(m)/b(m)

   do i = m - 2, 0, -1
      denom = b(i + 1) + c(i + 1)*g(i + 1)
      if (dabs(denom) .le. 1.d-70) then
         write (6, *) 'tridag a echoue a i et t=', i, t
         WRITE (6, *) 'denom,b,c,g=', denom, b(i + 2), c(i + 2), g(i + 2), b(i + 1), c(i + 1), g(i + 1)
         write (*, *) "pausing execution. Press enter to continue"
         read (*, *)

      end if
      g(i) = -a(i + 1)/denom
      h(i) = (r(i + 1) - c(i + 1)*h(i + 1))/denom
   end do

   sol(0) = (r(0) - c(0)*h(0))/(b(0) + c(0)*g(0))

   do i = 1, m
      sol(i) = sol(i - 1)*g(i - 1) + h(i - 1)
   end do
   return
end

!
!-----------------------------------------------------------------------------------------
!
subroutine init
   use general
   use par
! This routine defines all parameter values and initial conditions, to be passed in vector P.
   implicit none
   real(kind=8) dt, dx, Th, mua, rhoa, rhoc, rhot, rhow, D0ca, D0co3, Ka, Kc, beta, k2, k3, length, xdis
   real(kind=8) m, nn, S, phi0, ca0, co30, ccal0, cara0, Vscale, rhos0, Xs, Ts, eps
   real(kind=8) phi00, ca00, co300, ccal00, cara00, phiinf, bb
   integer i

!   DIMENSIONLESS TIMES, tmax=max time index, outx = time index for graphic x output
!   outt= time index for graphic t outputs at four points
!  dt, dx = discrete steps (dimensionless) ; N +1 = total number of nodes;
! N must be a multiple of 4
!  Th = thickness of active dissolution zone (cm); length = system size in cm, xdis=position of the beginning of the active dissolution zone (cm)
! mua in g/mol, densities in g/cm3, diff (in pure water) in cm2/a, Ka, Kc in M2
! beta in cm/a, k2,k3 in 1/a, m, n dimensionless
! S in cm/a, phio,cal0,cara0 dimensionless, ca0, co30 in M
! phiinf=porosity at the bottom of the compressed system
!  bb=sediment compressibility in (kPa)-1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   dt = 1.d-6
   xdis = 50.
   length = 500.
   Th = 100.
   eps = 1.d-2
! dt * tmax = for how many T* the model runs
   tmax=30000000 
   outt = 1000
   outx = tmax/4
   N = 200
   mua = 100.09
   rhoa = 2.95
   rhoc = 2.71
   rhot = 2.8
   rhow = 1.023
   D0ca = 131.9
   D0co3 = 272.6
   Ka = 10.**(-6.19)
   Kc = 10.**(-6.37)
   beta = 0.1
   k2 = 1.0
   k3 = 0.1
   m = 2.48
   nn = 2.8
   S = 0.1
   bb = 5.0d0
   phiinf = eps
!  new incoming sediment
! surface porosity
   phi0 = 0.8 
   ca0 = 0.326d-3
   co30 = 0.326d-3
   ccal0 = 0.3
   cara0 = 0.6
!  old uniform sediment
! initial porosity profile
   phi00 = 0.8 
   ca00 = 0.326d-3
   co300 = 0.326d-3
   ccal00 = 0.3
   cara00 = 0.6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In other versions, beta was used to scale the velocity
! Scale velocity with beta
   Vscale = beta
! Scale velocity with S
   Vscale = S
   rhos0 = rhot/(1 - cara0*(1 - rhot/rhoa) - ccal0*(1 - rhot/rhoc))
   Xs = D0ca/Vscale
   Ts = D0ca/Vscale/Vscale
   dx = (length/Xs)/N
! Assign the array of parameters. The missing members refer to parameters that are no longer used in this version
   P(1) = m
   P(2) = nn
   P(3) = k3/k2
   P(5) = cara0/ccal0
   P(6) = rhos0*ccal0/(mua*dsqrt(Kc))
   P(7) = rhos0*ccal0/rhoc
   P(8) = phi0
   P(9) = 1.-rhot/rhoc
   P(10) = 1.-rhot/rhoa
   P(11) = rhot/rhow
   P(12) = rhoc/rhow
   P(13) = rhos0/rhow - 1.
   P(14) = Th/Xs
   P(15) = dt
   P(16) = dx
   P(17) = 1 - cara0 - ccal0
   P(18) = cara0
   P(19) = ccal0
   P(20) = rhos0/rhow
   P(22) = xdis/Xs
   P(23) = S/Vscale
   P(24) = Kc/Ka
   P(25) = D0co3/D0ca
   P(26) = k2*Ts
   P(27) = rhoc/rhoa
   P(29) = eps
   P(32) = beta/Vscale
   P(34) = phiinf
! Porosity diffusion.  9.8 is g; 100 converts 1/cm to 1/m
   P(35) = 100.*beta*phi00**3/((phi00 - phiinf)*bb*rhow*9.8*(1 - phi00))/D0ca
!  case where all densities are equal to rhos
   P(7) = ccal0
   P(9) = 0.
   P(10) = 0.
   P(11) = rhos0/rhow
   P(12) = P(11)
   P(27) = 1.

! Initial conditions (dimensionless)

   pho(0) = phi0
   cao(0) = ca0/dsqrt(Kc)
   coo(0) = co30/dsqrt(Kc)
   ARAo(0) = cara0
   CALo(0) = ccal0

   do i = 1, N
      pho(i) = phi00
      cao(i) = ca00/dsqrt(Kc)
      coo(i) = co300/dsqrt(Kc)
      ARAo(i) = cara00
      CALo(i) = ccal00
   end do

   write (6, *) 'Xs (cm), Ts (a)', Xs, Ts
   write (6, *) 'dt,dx,dtS/dx =', dt, dx, P(15)*(P(23)/P(16))
   WRITE (6, *) 'dx^2/2d=', P(16)**2/(2*P(25))

   IF (P(15)*P(23)/P(16) .gt. 1./5.) then
      WRITE (6, *) ' problem: possible instability'
      write (*, *) "pausing execution. Press enter to continue"
      read (*, *)
   end if

   write (6, *) 'scale for MA, MC =', rhos0*cara0*(1 - phi0), rhos0*ccal0*(1 - phi0)
   write (6, *) 'scale for c=', dsqrt(Kc)
   write (6, *) 'Damkohler number Da=, scaled sedimentation rate,rhos0/rhow=', P(26), P(23), rhos0/rhow
   WRITE (6, *) 'scaled length, positions of dissolution zone=', length/Xs, xdis/Xs, (xdis + Th)/Xs
   WRITE (6, *) 'rhosw-1=', P(13)
   write (6, *) 'Dpor/Dca=', P(35)

   ! export parameters needed for further analysis
   open (1000, FILE='params')
   write (1000, *) Xs, Ts, S, phi0, bb, rhow, rhos0
   close (1000)

   ! Wait for an ackowledgment from the screen before starting
   write (*, *) "If parameters are OK, press [enter] to start run."
   read (*, *)
   return
end

