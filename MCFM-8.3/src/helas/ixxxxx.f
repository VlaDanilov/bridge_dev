c
c	Subroutine returns the desired fermion or
c	anti-fermion spinor. ie., |f>
c	A replacement for the HELAS routine IXXXXX
c
c	Adam Duff,  1992 August 31
c	<duff@phenom.physics.wisc.edu>
c
      subroutine ixxxxx(
     &              p,		!in: four vector momentum
     &              fmass,	!in: fermion mass
     &              nhel,	!in: spinor helicity, -1 or 1
     &              nsf,	!in: -1=antifermion, 1=fermion
     &              fi		!out: fermion wavefunction
     &                 )
      implicit none
c
c declare input/output variables
c
      complex*16 fi(6)
      integer*4 nhel, nsf
      real*8 p(0:3), fmass
c
c declare local variables
c
      real*8 r_zero, r_one, r_two
      parameter( r_zero=0.0d0, r_one=1.0d0, r_two=2.0d0 )
      complex*16 c_zero
      parameter( c_zero=(0d0,0d0) )
c
      real*8 plat, pabs, omegap, omegam, rs2pa, spaz
c
c define kinematic parameters
c
      fi(5) = dcmplx( p(0), p(3) ) * nsf
      fi(6) = dcmplx( p(1), p(2) ) * nsf
      plat = sqrt( p(1)**2 + p(2)**2 )
      pabs = sqrt( p(1)**2 + p(2)**2 + p(3)**2 )
      omegap = sqrt( p(0) + pabs )
c
c do massive fermion case
c
      if ( fmass .ne. r_zero ) then
         omegam = fmass / omegap
         if ( nsf .eq. 1 ) then
            if ( nhel .eq. 1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = dcmplx( omegam, r_zero )
                     fi(2) = c_zero
                     fi(3) = dcmplx( omegap, r_zero )
                     fi(4) = c_zero
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = omegam * rs2pa
     &                     * dcmplx( spaz, r_zero )
                     fi(2) = omegam * rs2pa / spaz
     &                     * dcmplx( p(1), p(2) )
                     fi(3) = omegap * rs2pa
     &                     * dcmplx( spaz, r_zero )
                     fi(4) = omegap * rs2pa / spaz
     &                     * dcmplx( p(1), p(2) )
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = dcmplx( omegam, r_zero )
                     fi(3) = c_zero
                     fi(4) = dcmplx( omegap, r_zero )
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = omegam * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(2) = omegam * rs2pa * spaz / plat
     &                     * dcmplx( p(1), p(2) )
                     fi(3) = omegap * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(4) = omegap * rs2pa * spaz / plat
     &                     * dcmplx( p(1), p(2) )
                  end if
               end if
            else if ( nhel .eq. -1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = dcmplx( omegap, r_zero )
                     fi(3) = c_zero
                     fi(4) = dcmplx( omegam, r_zero )
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = omegap * rs2pa / spaz
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = omegap * rs2pa
     &                     * dcmplx( spaz, r_zero )
                     fi(3) = omegam * rs2pa / spaz
     &                     * dcmplx( -p(1), p(2) )
                     fi(4) = omegam * rs2pa
     &                     * dcmplx( spaz, r_zero )
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = dcmplx( -omegap, r_zero )
                     fi(2) = c_zero
                     fi(3) = dcmplx( -omegam, r_zero )
                     fi(4) = c_zero
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = omegap * rs2pa * spaz / plat
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = omegap * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(3) = omegam * rs2pa * spaz / plat
     &                     * dcmplx( -p(1), p(2) )
                     fi(4) = omegam * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                  end if
               end if
            else
               stop 'ixxxxx:  fermion helicity must be +1,-1'
            end if
         else if ( nsf .eq. -1 ) then
            if ( nhel .eq. 1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = dcmplx( -omegap, r_zero )
                     fi(3) = c_zero
                     fi(4) = dcmplx( omegam, r_zero )
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = -omegap * rs2pa / spaz
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = -omegap * rs2pa
     &                     * dcmplx( spaz, r_zero )
                     fi(3) = omegam * rs2pa / spaz
     &                     * dcmplx( -p(1), p(2) )
                     fi(4) = omegam * rs2pa
     &                     * dcmplx( spaz, r_zero )
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = dcmplx( omegap, r_zero )
                     fi(2) = c_zero
                     fi(3) = dcmplx( -omegam, r_zero )
                     fi(4) = c_zero
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = -omegap * rs2pa * spaz / plat
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = -omegap * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(3) = omegam * rs2pa * spaz / plat
     &                     * dcmplx( -p(1), p(2) )
                     fi(4) = omegam * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                  end if
               end if
            else if ( nhel .eq. -1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = dcmplx( omegam, r_zero )
                     fi(2) = c_zero
                     fi(3) = dcmplx( -omegap, r_zero )
                     fi(4) = c_zero
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = omegam * rs2pa
     &                     * dcmplx( spaz, r_zero )
                     fi(2) = omegam * rs2pa / spaz
     &                     * dcmplx( p(1), p(2) )
                     fi(3) = -omegap * rs2pa
     &                     * dcmplx( spaz, r_zero )
                     fi(4) = -omegap * rs2pa / spaz
     &                     * dcmplx( p(1), p(2) )
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = dcmplx( omegam, r_zero )
                     fi(3) = c_zero
                     fi(4) = dcmplx( -omegap, r_zero )
                  else
                     rs2pa = r_one / sqrt( r_two * pabs )
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = omegam * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(2) = omegam * rs2pa * spaz / plat
     &                     * dcmplx( p(1), p(2) )
                     fi(3) = -omegap * rs2pa / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(4) = -omegap * rs2pa * spaz / plat
     &                     * dcmplx( p(1), p(2) )
                  end if
               end if
            else
               stop 'ixxxxx:  fermion helicity must be +1,-1'
            end if
         else
            stop 'ixxxxx:  fermion type must be +1,-1'
         end if
c
c do massless fermion case
c
      else
         if ( nsf .eq. 1 ) then
            if ( nhel .eq. 1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = dcmplx( omegap, r_zero )
                     fi(4) = c_zero
                  else
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = dcmplx( spaz, r_zero )
                     fi(4) = r_one / spaz
     &                     * dcmplx( p(1), p(2) )
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = c_zero
                     fi(4) = dcmplx( omegap, r_zero )
                  else
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = r_one / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(4) = spaz / plat
     &                     * dcmplx( p(1), p(2) )
                  end if
               end if
            else if ( nhel .eq. -1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = dcmplx( omegap, r_zero )
                     fi(3) = c_zero
                     fi(4) = c_zero
                  else
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = r_one / spaz
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = dcmplx( spaz, r_zero )
                     fi(3) = c_zero
                     fi(4) = c_zero
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = dcmplx( -omegap, r_zero )
                     fi(2) = c_zero
                     fi(3) = c_zero
                     fi(4) = c_zero
                  else
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = spaz / plat
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = r_one / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(3) = c_zero
                     fi(4) = c_zero
                  end if
               end if
            else
               stop 'ixxxxx:  fermion helicity must be +1,-1'
            end if
         else if ( nsf .eq. -1 ) then
            if ( nhel .eq. 1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = dcmplx( -omegap, r_zero )
                     fi(3) = c_zero
                     fi(4) = c_zero
                  else
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = -r_one / spaz
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = dcmplx( -spaz, r_zero )
                     fi(3) = c_zero
                     fi(4) = c_zero
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = dcmplx( omegap, r_zero )
                     fi(2) = c_zero
                     fi(3) = c_zero
                     fi(4) = c_zero
                  else
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = -spaz / plat
     &                     * dcmplx( -p(1), p(2) )
                     fi(2) = -r_one / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(3) = c_zero
                     fi(4) = c_zero
                  end if
               end if
            else if ( nhel .eq. -1 ) then
               if ( p(3) .ge. r_zero ) then
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = dcmplx( -omegap, r_zero )
                     fi(4) = c_zero
                  else
                     spaz = sqrt( pabs + p(3) )
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = dcmplx( -spaz, r_zero )
                     fi(4) = -r_one / spaz
     &                     * dcmplx( p(1), p(2) )
                  end if
               else
                  if ( plat .eq. r_zero ) then
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = c_zero
                     fi(4) = dcmplx( -omegap, r_zero )
                  else
                     spaz = sqrt( pabs - p(3) )
                     fi(1) = c_zero
                     fi(2) = c_zero
                     fi(3) = -r_one / spaz
     &                     * dcmplx( plat, r_zero )
                     fi(4) = -spaz / plat
     &                     * dcmplx( p(1), p(2) )
                  end if
               end if
            else
               stop 'ixxxxx:  fermion helicity must be +1,-1'
            end if
         else
            stop 'ixxxxx:  fermion type must be +1,-1'
         end if
      end if
c
c done
c
      return
      end
