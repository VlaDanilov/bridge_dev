      subroutine scaleset_ddis(p,mu0)
        use singletop2_scale_m
      implicit none
      include 'types.f'
c--- subroutine to calculate dynamic scale equal to
c--- "double deep-inelastic scale" defined for the single top t-channel by:
c---        Q^2=-(pt-pb)^2 on the light quark line
c---        Q^2+mt^2 on the heavy quark line
c---  (c.f. Z. Sullivan, hep-ph/0408049)
c--- AVAILABLE AT LEADING ORDER ONLY
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'kprocess.f'
      include 'masses.f'
      include 'kpart.f'
      real(dp):: p(mxpart,4),mu0
      include 'bqscale.f'
      real(dp) :: dotvec

      if((kcase==kbq_tpq) .or.
     &   (kcase==kqg_tbq) .or. kcase==ktopanom) then
        if (kcase==kbq_tpq .or. kcase==ktopanom) then

            q1scale = -dotvec(p(1,:)+p(6,:), p(1,:)+p(6,:))
            q2scale = -dotvec(p(2,:)+p(6,:), p(2,:)+p(6,:))
            b1scale = (q2scale + mt**2)
            b2scale = (q1scale + mt**2)

        else
          q1scale=-(
     &         +(p(5,4)+p(1,4))**2
     &         -(p(5,1)+p(1,1))**2
     &         -(p(5,2)+p(1,2))**2
     &         -(p(5,3)+p(1,3))**2)
          b2scale=-(
     &         +(p(5,4)+p(1,4))**2
     &         -(p(5,1)+p(1,1))**2
     &         -(p(5,2)+p(1,2))**2
     &         -(p(5,3)+p(1,3))**2)+mt**2
          q2scale=-(
     &         +(p(5,4)+p(2,4))**2
     &         -(p(5,1)+p(2,1))**2
     &         -(p(5,2)+p(2,2))**2
     &         -(p(5,3)+p(2,3))**2)
          b1scale=-(
     &         +(p(5,4)+p(2,4))**2
     &         -(p(5,1)+p(2,1))**2
     &         -(p(5,2)+p(2,2))**2
     &         -(p(5,3)+p(2,3))**2)+mt**2
      endif
      q1scale=max(sqrt(q1scale),one) ! min. of 1 GeV for safety
      b2scale=max(sqrt(b2scale),one)
      q2scale=max(sqrt(q2scale),one)
      b1scale=max(sqrt(b1scale),one)
      mu0=q1scale   ! for safety
      else
        write(6,*) 'dynamicscale DDIS not supported for this process.'
        stop
      endif
      
c     if (kpart.ne.klord) then
c       write(6,*) 'dynamicscale DDIS only available at LO.'
c       stop      
c     endif
      
      return
      end
      
