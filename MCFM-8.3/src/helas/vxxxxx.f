c
c
c	Subroutine returns the value of evaluated
c	helicity basis boson polarisation wavefunction.
c	Replaces the HELAS routine VXXXXX
c
c	Adam Duff,  1992 September 3
c	<duff@phenom.physics.wisc.edu>
c
      subroutine vxxxxx(
     &              p,		!in: boson four momentum
     &              vmass,	!in: boson mass
     &              nhel,	!in: boson helicity
     &              nsv,	!in: incoming (-1) or outgoing (+1)
     &              vc		!out: boson wavefunction
     &                 )
      implicit none
c
c declare input/output variables
c
      complex*16 vc(6)
      integer*4 nhel, nsv
      real*8 p(0:3), vmass
c
c declare local variables
c
      real*8 r_zero, r_one, r_two
      parameter( r_zero=0.0d0, r_one=1.0d0, r_two=2.0d0 )
      complex*16 c_zero
      parameter( c_zero=(0d0,0d0) )
c
      real*8 plat, pabs, rs2, rplat, rpabs, rden
c
c define internal/external momenta
c
      if ( nsv**2 .ne. 1 ) then
         stop 'vxxxxx:  nsv is not one of -1, +1'
      end if
c
      rs2 = sqrt( r_one / r_two )
      vc(5) = dcmplx( p(0), p(3) ) * nsv
      vc(6) = dcmplx( p(1), p(2) ) * nsv
      plat = sqrt( p(1)**2 + p(2)**2 )
      pabs = sqrt( p(1)**2 + p(2)**2 + p(3)**2 )
c
c calculate polarisation four vectors
c
      if ( nhel**2 .eq. 1 ) then
         if ( (pabs .eq. r_zero) .or. (plat .eq. r_zero) ) then
            vc(1) = c_zero
            vc(2) = dcmplx( -nhel * rs2 * dsign( r_one, p(3) ), r_zero )
            vc(3) = dcmplx( r_zero, nsv * rs2 )
            vc(4) = c_zero
         else
            rplat = r_one / plat
            rpabs = r_one / pabs
            vc(1) = c_zero
            vc(2) = dcmplx( -nhel * rs2 * rpabs * rplat * p(1) * p(3),
     &                     -nsv * rs2 * rplat * p(2) )
            vc(3) = dcmplx( -nhel * rs2 * rpabs * rplat * p(2) * p(3),
     &                     nsv * rs2 * rplat * p(1) )
            vc(4) = dcmplx( nhel * rs2 * rpabs * plat,
     &                     r_zero )
         end if
      else if ( nhel .eq. 0 ) then
         if ( vmass .gt. r_zero ) then
            if ( pabs .eq. r_zero ) then
               vc(1) = c_zero
               vc(2) = c_zero
               vc(3) = c_zero
               vc(4) = dcmplx( r_one, r_zero )
            else
               rden = p(0) / ( vmass * pabs )
               vc(1) = dcmplx( pabs / vmass, r_zero )
               vc(2) = dcmplx( rden * p(1), r_zero )
               vc(3) = dcmplx( rden * p(2), r_zero )
               vc(4) = dcmplx( rden * p(3), r_zero )
            end if
         else
            stop  'vxxxxx: nhel = 0 is only for massive bosons'
         end if
      else if ( nhel .eq. 4 ) then
         if ( vmass .gt. r_zero ) then
            rden = r_one / vmass
            vc(1) = dcmplx( rden * p(0), r_zero )
            vc(2) = dcmplx( rden * p(1), r_zero )
            vc(3) = dcmplx( rden * p(2), r_zero )
            vc(4) = dcmplx( rden * p(3), r_zero )
         elseif (vmass .eq. r_zero) then
            rden = r_one / p(0)
            vc(1) = dcmplx( rden * p(0), r_zero )
            vc(2) = dcmplx( rden * p(1), r_zero )
            vc(3) = dcmplx( rden * p(2), r_zero )
            vc(4) = dcmplx( rden * p(3), r_zero )
         else
            stop 'vxxxxx: nhel = 4 is only for m>=0'
         end if
      else
         stop 'vxxxxx:  nhel is not one of -1, 0, 1 or 4'
      end if
c
c done
c
      return
      end
