      subroutine qqb_QQb_mix_g(p,msq)
      implicit none
c---Matrix element squared averaged over initial colors and spins
c--- for mixed QCD/Z matrix element for top pair production,
c     q(-p1)+qbar(-p2) -->  Q(p3)+Q~(P4) + g(p5)

c--- all momenta are incoming
c

c--- Overall normalization agrees with the massless limit
c--- given in Eq.(III.19) of Kuhn et al., arXiv:0909.0059
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     .     gvq(nf),gaq(nf),gvt,gat,T3(nf),sw2,cw2,s12,s34,sz12,sz34,
     .     t1,t2,u1,u2,genfac,p1dp2,p3dp4,p1dp5,p2dp5,p3dp5,p4dp5,
     .     fac1,fac2,ttv1,ttv2,tta1,tta2,mz,qqb(nf),qbq(nf),cs
      real(dp):: p1Dp3,p1Dp4,p2Dp3,p2Dp4
      real(dp):: pt(4),ptt(4),p5(4)

      msq=0._dp

      sw2 = xw
!      sw2 = 1._dp - wmass**2/zmass**2
      cw2 = 1._dp - sw2

C -- check at one phase space point with Doreen
      if(1 == 2) then
         p(1,1) = 0._dp
         p(1,2) =  0._dp
         p(1,3) = - 0.25d4
         p(1,4) = - 0.25d4

         p(2,1) = 0._dp
         p(2,3) = 0._dp
         p(2,3) = + 0.25d4
         p(2,4) = - 0.25d4

         p(3,1) = + 0.127853110965d4
         p(3,2) = + 0.121447514622d4
         p(3,3) = + 0.151351146245d4
         p(3,4) = + 0.233026477382d4

         p(4,1) = - 0.254032886288d3
         p(4,2) = - 0.152708137180d3
         p(4,3) = + 0.209002671061d3
         p(4,4) = + 0.401696264592d3

         p(5,1) = - 0.102449822336d4
         p(5,2) = - 0.106176700904d4
         p(5,3) = - 0.172251413351d4
         p(5,4) = + 0.226803896159d4

         sw2 = 1._dp - wmass**2/zmass**2
         cw2 = 1._dp - sw2
      end if


      T3 = (/-0.5_dp,0.5_dp,-0.5_dp,0.5_dp,-0.5_dp/)
      gvq(1:5) = (T3(1:5)-2._dp*sw2*Q(1:5))/2._dp/sqrt(sw2*cw2)
      gaq(1:5) = T3(1:5)/2._dp/sqrt(sw2*cw2)
      gw = 0.5_dp/sqrt(2._dp*sw2)

      gvt = gvq(2)
      gat = gaq(2)

      mz = zmass

      call dotem(5,p,s)

      s12 = s(1,2)
      t1 = s(1,3)
      t2 = s(2,4)
      u1 = s(1,4)
      u2 = s(2,3)
      s34 = s(3,4) + 2._dp*mt**2
      p1dp2 = s12                 ! 2*p1.p2 = s12
      p3dp4 = s(3,4)              ! 2*p3.p4 = -(s12+t1+t2+u1+u2+2*mt**2)
      p1dp5 = s(1,5)              ! 2*p1.p5 = -(s12+t1+u1)
      p2dp5 = s(2,5)              ! 2*p2.p5 = -(s12+t2+u2)
      p3dp5 = s(3,5)              ! 2*p3.p5 = +(s12+t2+u1)
      p4dp5 = s(4,5)              ! 2*p4.p5 = +(s12+t1+u2)
      sz12 = s12 - mz**2
      sz34 = s34 - mz**2          ! sz34 = -(s12+t1+t2+u1+u2+mz**2)
! fac of 16 is due to the part (multiply,p1.p5*p2.p5*p3.p5*p4.p5*s12*s34*sz12*sz34;) 
! refer to "qqbttbg_new.frm"; there is a factor of 2 comes from each pi.pj; hence 2^4 = 16.
! and the fac 2 multiplied at the end is because we calculate 
! Amp(qqb->g->ttb+g)*conj(Amp(qqb->z->ttb+g)) while the full tree level xsec has another 
! part which is the conjugate of which we calucate

      genfac = 16._dp/p1dp5/p2dp5/p3dp5/p4dp5/s12/s34/sz12/sz34*2._dp

      fac1 = t1**2 + t2**2 + u1**2 + u2**2
      fac2 = t1**2 + t2**2 - u1**2 - u2**2

      ttv1 = 
     .     + 0.5_dp*(-mz**2*(t1 + t2 + u1 + u2) + 2._dp*s12*(s12 + t1 + t2 
     .     + u1 + u2))*(s12**2*(t1 + t2 - u1 - u2) + s12*(t1 + t2 - u1 
     .     - u2)*(t1 + t2 + u1 + u2) + (t1 + t2 + u1 + u2)*(t1*t2 
     .     - u1*u2))

      ttv2 = 
     .     + mt**2*(-mz**2*(t1 + t2 + u1 + u2) + 2._dp*s12*(s12 + t1 + t2 
     .     + u1 + u2))*(2._dp*s12**3*(t1 + t2 - u1 - u2) - (t1 - t2 + u1 
     .     - u2)**2*(t1*t2 - u1*u2) + s12**2*(t1**2 + 6._dp*t1*t2 + t2**2 
     .     - u1**2 - 6._dp*u1*u2 - u2**2) + s12*(-t1**3 - t2**3 
     .     + t2**2*(u1 - u2) + t1**2*(3._dp*t2 - u1 + u2) 
     .     + (u1 + u2)*(u1**2 - 4._dp*u1*u2 + u2**2) + t2*(-u1**2 
     .     - 2._dp*u1*u2 + u2**2) + t1*(3._dp*t2**2 + u1**2 - 2._dp*u1*u2 
     .     - u2**2 + 2._dp*t2*(u1 + u2))))


 !     tta1 = 
 !    .     - 0.5_dp*(-mz**2*(t1 + t2 + u1 + u2) + 2._dp*s12*(s12 + t1 + t2 
 !    .     + u1 + u2))*(s12**2*(t1 + t2 - u1 - u2) + s12*(t1 + t2 - u1 
 !    .     - u2)*(t1 + t2 + u1 + u2) + (t1 + t2 + u1 + u2)*(t1*t2 
 !    .     - u1*u2))

      tta1 = - ttv1

      tta2 = 
     .     -2._dp*mt**2*s12*(s12 + t1 + u1)*(s12 + t2 + u2)*(t1 + t2 + u1 
     .     + u2)*(mz**2 + s12 + t1 + t2 + u1 + u2)*(2._dp*s12 + t1 + t2 
     .     + u1 + u2)

      ttv1 = fac1*ttv1
      tta1 = fac2*tta1
      
      qqb = genfac*((ttv1+ttv2)*gvq*gvt + (tta1+tta2)*gaq*gat)
c--- for qbq, the signs of ttv1 and ttv2 flip, tta1 and tta2 remain the same
      qbq = genfac*(-(ttv1+ttv2)*gvq*gvt + (tta1+tta2)*gaq*gat)

      do j=-nf,nf
         k=-j
         if((j == 0) .and. (k == 0)) then
            msq(j,k)=0._dp
         else if((j > 0) .and. (k < 0)) then
            msq(j,k) = qqb(j)
         else if((j < 0 ) .and. (k > 0)) then
            msq(j,k) = qbq(k)
         end if
      end do

      cs  = (Nc**2 -1._dp)/4._dp

! cs/Nc**2/4 is the color and spin average
      msq = msq*esq*gsq**2*cs/Nc**2/4._dp
!     .     *fourpi**3/126.3_dp/esq/gsq**2*0.01_dp

!      write(6,*) 'msq(1,-1): ', msq(1,-1)
!      stop

!      pt = -p(1,:)-p(2,:)
!      ptt = (/0._dp,0._dp,0._dp,sqrt(s12)/)
!      call boostx(-p(5,:),pt,ptt,p5)
!      call boostx(-p(1,:),pt,ptt,p(1,:))
!      call boostx(-p(2,:),pt,ptt,p(2,:))
!      write(6,*) "pt: ", pt
!      write(6,*) "ptt: ", ptt
!      write(6,*) "p5: ", p(5,:), s(5,5)
!      write(6,*) "p5_new: ", p5, p5(4)**2-p5(3)**2-p5(2)**2-p5(1)**2
!      write(6,*) "p1: ", p(1,:)
!      write(6,*) "p2: ", p(2,:)
!      pause
!      if (dabs(p5(4)) .lt. 1._dp) msq = 0._dp
!      if (p(5,4) .lt. 1._dp) msq = 0._dp

      
      end subroutine qqb_QQb_mix_g
