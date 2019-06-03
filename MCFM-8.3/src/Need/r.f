      function rpure(pi,pj)
          use types
          implicit none
          
          real(dp) :: rpure
          real(dp), intent(in) :: pi(4),pj(4)

          real(dp) :: pti2, ptj2, ei, ej, r1, r2, biti, bitj, delphi, dely
          real(dp), parameter :: tiny = 1.e-6_dp

          pti2=pi(1)**2+pi(2)**2
          ptj2=pj(1)**2+pj(2)**2

          ei=sqrt(pti2+pi(3)**2)
          ej=sqrt(ptj2+pj(3)**2)

          biti=pi(3)/ei
          bitj=pj(3)/ej

          if  ((abs(1d0+biti) < tiny) .or. (abs(1d0-biti) < tiny)
     &     .or.(abs(1d0+bitj) < tiny) .or. (abs(1d0-bitj) < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
              dely=100d0
          else
              r1=(1._dp+biti)*(1._dp-bitj)/((1._dp+bitj)*(1._dp-biti))
              dely=0.5_dp*log(r1)
          endif

          r2= (pi(1)*pj(1)+pi(2)*pj(2))/sqrt(pti2*ptj2)
          if (r2 > +0.9999999_dp) r2=+1._dp
          if (r2 < -0.9999999_dp) r2=-1._dp
          delphi=acos(r2)

          rpure = sqrt(dely**2+delphi**2)
      end function
      
      function r(p,i,j)
      implicit none
      include 'types.f'
      real(dp):: r
c----calculate the jets separation between p(i) and p(j)
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4),r1,r2,dely,delphi,ei,ej,pti2,ptj2,
     & biti,bitj
      integer:: i,j
      real(dp), parameter:: tiny=1.e-6

      pti2=p(i,1)**2+p(i,2)**2
      ptj2=p(j,1)**2+p(j,2)**2

      ei=sqrt(pti2+p(i,3)**2)
      ej=sqrt(ptj2+p(j,3)**2)

c      r1= (ei+p(i,3))*(ej-p(j,3))/
c     &   ((ej+p(j,3))*(ei-p(i,3)))
      biti=p(i,3)/ei
      bitj=p(j,3)/ej
      if  ((abs(1d0+biti) < tiny) .or. (abs(1d0-biti) < tiny)
     & .or.(abs(1d0+bitj) < tiny) .or. (abs(1d0-bitj) < tiny)) then
c-- set to 100 if any of these is very close to or less than zero
c-- rapidities of 100 will be rejected by any sensible cuts
        dely=100d0
      else
        r1=(one+biti)*(one-bitj)/((one+bitj)*(one-biti))
        dely=0.5_dp*dlog(r1)
      endif

      r2= (p(i,1)*p(j,1)+p(i,2)*p(j,2))/sqrt(pti2*ptj2)
      if (r2 > +0.9999999_dp) r2=+1._dp
      if (r2 < -0.9999999_dp) r2=-1._dp
      delphi=acos(r2)

      r=sqrt(dely**2+delphi**2)
      
      return
      end
      
      function delphi(pi,pj)
        include 'types.f'

        real(dp) :: delphi
        real(dp), intent(in) :: pi(4), pj(4)

        real(dp) :: pti2, ptj2

        pti2 = pi(1)**2+pi(2)**2
        ptj2 = pj(1)**2+pj(2)**2
        delphi = acos((pi(1)*pj(1)+pi(2)*pj(2))/sqrt(pti2*ptj2))

      end function

