      subroutine hjetmass_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'hdecaymode.f'
      include 'first.f'
C  in is the label of the momentum contracted with n
      integer j,k,in
      double precision dot
      double precision msq(-nf:nf,-nf:nf)
      double precision n(4),p(mxpart,4),hdecay,s34,fac,
     . p1p2(-1:1,-1:1),msqgamgam
      double complex hjetmass_gg_gvec_n1
      double complex hjetmass_gg_gvec_n2
      double complex hjetmass_gg_gvec_n3
      double precision hjetmass_qqghn
      double precision qqghn,ggghn
      

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      enddo
      enddo

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

C   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
           hdecay=msqgamgam(sqrt(s34))
      else
      write(6,*) 'Unimplemented process in gg_hgg_gvec'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      fac=hdecay

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      if (first) then
          call sushi_bernini(18)
          first = .false.
      end if

c     write (*,*) "check subtractions HTL"
      if (in .eq. 1) then
      p1p2(0,-1)=-aveqg*fac*hjetmass_qqghn(2,5,1,p,n)
c     write (*,*) p1p2(0,-1)/(-aveqg*fac*qqghn(2,5,1,p,n))
      p1p2(0,+1)=-aveqg*fac*hjetmass_qqghn(2,5,1,p,n)
c     write (*,*) p1p2(0,+1)/(-aveqg*fac*qqghn(2,5,1,p,n))
      p1p2(0,0)=dreal(hjetmass_gg_gvec_n1(p,n))*fac/256d0/16d0
c     write (*,*) p1p2(0,0)/(avegg*fac*ggghn(5,2,1,p,n))
      elseif (in .eq. 2) then
      p1p2(+1,0)=-aveqg*fac*hjetmass_qqghn(1,5,2,p,n)
c     write (*,*) p1p2(+1,0)/(-aveqg*fac*qqghn(1,5,2,p,n))
      p1p2(-1,0)=-aveqg*fac*hjetmass_qqghn(5,1,2,p,n)
c     write (*,*) p1p2(-1,0)/(-aveqg*fac*qqghn(5,1,2,p,n))
      p1p2(0,0)=dreal(hjetmass_gg_gvec_n2(p,n))*fac/256d0/16d0
c     write (*,*) p1p2(0,0)/(avegg*fac*ggghn(1,5,2,p,n))
      elseif (in .eq. 5) then     
      p1p2(1,-1)=+aveqq*fac*hjetmass_qqghn(1,2,5,p,n)
c     write (*,*) p1p2(1,-1)/(aveqq*fac*qqghn(1,2,5,p,n))
      p1p2(-1,1)=+aveqq*fac*hjetmass_qqghn(2,1,5,p,n)
c     write (*,*) p1p2(-1,1)/(aveqq*fac*qqghn(2,1,5,p,n))
      p1p2(0,0)=dreal(hjetmass_gg_gvec_n3(p,n))*fac/256d0/16d0
c     write (*,*) p1p2(0,0)/(avegg*fac*ggghn(1,2,5,p,n))
      endif
c     write (*,*) "end check subtractions HTL"

      do j=-nf,nf
      do k=-nf,nf
      if     ((j .gt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j .lt. 0) .and. (k .eq. -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j .eq. 0) .and. (k .eq. 0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j .gt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(+1,0)
      elseif ((j .lt. 0) .and. (k .eq. 0)) then
          msq(j,k)=
     &    p1p2(-1,0)
      elseif ((j .eq. 0) .and. (k .gt. 0)) then
          msq(j,k)=
     &    p1p2(0,+1)
      elseif ((j .eq. 0) .and. (k .lt. 0)) then
          msq(j,k)=
     &    p1p2(0,-1)
      endif
      enddo
      enddo
 
      return
      end

      ! taken from ggHg/gg_hg_gvec.f and augmented by mass factor
      double precision function hjetmass_qqghn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j1,j2,j5
      double precision Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u
      double precision mDot, aex
      double complex asmh

      nDp1 = mDot(p(j1,:), n)
      nDp2 = mDot(p(j2,:), n)
      nDn  = mDot(n, n)

      call checkndotp(p,n,j5)

      Asq=(as/(3d0*pi))**2/vevsq

      s=2d0*Dot(p,j1,j2)
      t=2d0*Dot(p,j1,j5)
      u=2d0*Dot(p,j2,j5)

      ! this becomes 1 in HTL:
      aex = (abs(asmh(s,t+u,mt**2) * mt**2))**2

      hjetmass_qqghn=-Asq*gsq*V/2d0*(2d0*(nDp1*u-nDp2*t)**2/s**2
     .                    +0.5d0*nDn*(u+t)**2/s) * aex

      return
      end
      
      double precision function mDot(p,q)
          double precision p(4), q(4)
          mDot = p(4)*q(4) - sum(p(1:3)*q(1:3))
      end function

      double complex function hjetmass_gq_gvec_n3(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      double precision p(mxpart,4), next(4)
      double complex squaredamp, asmh, aex
      double precision sman, tman, uman
      double precision aex_prefactor
      double precision p1_next, p2_next, p3_next, next_next
      double precision mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      aex_prefactor = sqrt(gsq)**3 / 32d0 / pi**2 / sqrt(vevsq)
      aex = 1/sman*aex_prefactor*8d0/3d0
      ! this (squared) becomes 1 in HTL:
      !aex = aex * asmh(sman,tman+uman,mt**2) * mt**2

      squaredamp =
     &  + dconjg(aex)*aex * ( 4.D0*p1_next*p2_next*tman*uman - 2.D0*
     &    p1_next*p3_next*sman*uman - 2.D0*p1_next*p3_next*sman*tman - 
     &    2.D0*p1_next**2*uman**2 - 2.D0*p2_next*p3_next*sman*uman - 2.D
     &    0*p2_next*p3_next*sman*tman - 2.D0*p2_next**2*tman**2 - 1.D0/
     &    2.D0*next_next*sman*uman**2 - next_next*sman*tman*uman - 1.D0/
     &    2.D0*next_next*sman*tman**2 )

      hjetmass_gq_gvec_n3 = squaredamp * 8d0 / 2d0

      end function

      double complex function hjetmass_gg_gvec_n1(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      double precision p(mxpart,4), next(4)
      double complex squaredamp
      double precision cex_prefactor
      double complex c1ex, c2ex, c3ex, c4ex
      double complex c1smh, c2smh
      double precision sman, tman, uman
      double precision p1_next, p2_next, p3_next, next_next
      double precision mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

      c1ex = c1smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c2ex = c2smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c3ex = c2smh(tman,sman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c4ex = c2smh(uman,sman,tman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      squaredamp =
     &  + dconjg(c1ex)*c4ex * (  - 16.D0*p1_next*p2_next*sman**(-1) + 
     &    16.D0*p1_next*p3_next*tman**(-1) + 8.D0*p1_next**2*sman**(-1)
     &    *tman**(-1)*uman + 32.D0*p2_next*p3_next*uman**(-1) + 16.D0*
     &    p2_next**2*sman**(-1)*tman*uman**(-1) + 16.D0*p3_next**2*sman
     &    *tman**(-1)*uman**(-1) + 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c1ex)*c3ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c1ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c4ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c3ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman + 32.D0*p2_next*p3_next
     &    *uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*uman**(-1) + 
     &    16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0*next_next
     &     )
      squaredamp = squaredamp + dconjg(c2ex)*c2ex * ( 32.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next**2*sman**(-1)*tman**(-1)*
     &    uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c1ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c4ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c3ex * (  - 32.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p1_next**2*sman**(-1)*tman**(-1)*
     &    uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c2ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman + 32.D0*p2_next*p3_next
     &    *uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*uman**(-1) + 
     &    16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0*next_next
     &     )
      squaredamp = squaredamp + dconjg(c3ex)*c1ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c4ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c3ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c1ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman + 32.D0*p2_next*p3_next
     &    *uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*uman**(-1) + 
     &    16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0*next_next
     &     )

       hjetmass_gg_gvec_n1 = squaredamp * 24 * 16 / sman / tman / uman

      end function

      double complex function hjetmass_gg_gvec_n2(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      double precision p(mxpart,4), next(4)
      double complex squaredamp
      double precision cex_prefactor
      double complex c1ex, c2ex, c3ex, c4ex
      double complex c1smh, c2smh
      double precision sman, tman, uman
      double precision p1_next, p2_next, p3_next, next_next
      double precision mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

      c1ex = c1smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c2ex = c2smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c3ex = c2smh(tman,sman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c4ex = c2smh(uman,sman,tman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      squaredamp =
     &  + dconjg(c1ex)*c4ex * (  - 8.D0*p1_next*p2_next*sman**(-1) - 8.D
     &    0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c3ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c1ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c1ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p2_next**2*sman**(-1)*tman*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c4ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) - 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c3ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c2ex * ( 32.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c1ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c4ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c3ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p2_next**2*sman**(-1)*tman*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c1ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c4ex * (  - 32.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c3ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c2ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) - 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c1ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )

       hjetmass_gg_gvec_n2 = squaredamp * 24 * 16 / sman / tman / uman

      end function

      double complex function hjetmass_gg_gvec_n3(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      double precision p(mxpart,4), next(4)
      double complex squaredamp
      double precision cex_prefactor
      double complex c1ex, c2ex, c3ex, c4ex
      double complex c1smh, c2smh
      double precision sman, tman, uman
      double precision p1_next, p2_next, p3_next, next_next
      double precision mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

      c1ex = c1smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c2ex = c2smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c3ex = c2smh(tman,sman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c4ex = c2smh(uman,sman,tman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      squaredamp =
     &  + dconjg(c1ex)*c4ex * ( 8.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c3ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c2ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c1ex)*c1ex * (  - 16.D0*p1_next*
     &    p3_next*tman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p3_next**2*sman*tman**(-1)*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c4ex * ( 8.D0*p1_next*
     &    p3_next*tman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c3ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c2ex * (  - 16.D0*p1_next*
     &    p3_next*tman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p3_next**2*sman*tman**(-1)*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c1ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c4ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c3ex * (  - 32.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c2ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c1ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c4ex * (  - 32.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c3ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c2ex * ( 8.D0*p1_next*
     &    p3_next*tman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c1ex * ( 8.D0*p1_next*
     &    p3_next*tman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )

       hjetmass_gg_gvec_n3 = squaredamp * 24 * 16 / sman / tman / uman

      end function

