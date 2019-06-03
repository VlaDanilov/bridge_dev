      subroutine qqb_QQb_g_ew(p,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      integer j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),rl(4,4,nf),
     &     gvq(nf),gaq(nf),gvt,gat,T3(nf),sw2,cw2,s12,s34,sz12,sz34,
     &     t1,t2,u1,u2,genfac,p1dp2,p1dp5,p2dp5,p3dp5,p4dp5,
     &     fac1,fac2,ttv1,ttv2,tta1,tta2,mz,qqb(nf),cs


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
      s34 = s(3,4)   ! s34 = -(s12+t1+t2+u1+u2)
      p1dp2 = s12    ! 2*p1.p2 = s12
      p1dp5 = s(1,5) ! 2*p1.p5 = -(s12+t1+u1)
      p2dp5 = s(2,5) ! 2*p2.p5 = -(s12+t2+u2)
      p3dp5 = s(3,5) ! 2*p3.p5 = +(s12+t2+u1)
      p4dp5 = s(4,5) ! 2*p4.p5 = +(s12+t1+u2)
      sz12 = s12 - mz**2
      sz34 = s34 - mz**2 ! sz34 = -(s12+t1+t2+u1+u2+mz**2)
! fac of 32 is from p1.p5,...,p4.p5, as well as s12 which is used p1.p2 in form 
! so there is a factor of 2; hence 2^5 = 32.
      genfac = 32._dp/p1dp5/p2dp5/p3dp5/p4dp5/s12/s34/sz12/sz34

      fac1 = t1**2 + t2**2 + u1**2 + u2**2

      fac2 = t1**2 + t2**2 - u1**2 - u2**2


      ttv1 = 
     .     + 0.125_dp*(-mz**2*(t1 + t2 + u1 + u2) 
     .     + 2._dp*s12*(s12 + t1 + t2 + u1 + u2))*
     .     (s12**2*(t1 + t2 - u1 - u2) 
     .     + s12*(t1 + t2 - u1 - u2)*(t1 + t2 + u1 + u2) 
     .     + (t1 + t2 + u1 + u2)*(t1*t2 - u1*u2))



      ttv2 = 
     .     + 0.125_dp*((t1**2 + t2**2 + u1**2 + u2**2)*
     .     (mz**2*(t1 + t2 + u1 + u2) 
     .     - 2._dp*s12*(s12 + t1 + t2 + u1 + u2))*
     .     (s12**2*(t1 + t2 - u1 - u2) 
     .     + s12*(t1 + t2 - u1 - u2)*(t1 + t2 + u1 + u2) 
     .     + (t1 + t2 + u1 + u2)*(t1*t2 - u1*u2)) 
     .     - (2._dp*mt**2*(mz**2 - s12) + mz**2*(t1 + t2 + u1 + u2) 
     .     - 2._dp*s12*(s12 + t1 + t2 + u1 + u2))*((t1**2 + t2**2 
     .     + u1**2 + u2**2)*(s12**2*(t1 + t2 - u1 - u2) 
     .     + s12*(t1 + t2 - u1 - u2)*(t1 + t2 + u1 + u2) 
     .     + (t1 + t2 + u1 + u2)*(t1*t2 - u1*u2)) 
     .     + 2._dp*mt**2*(2._dp*s12**3*(t1 + t2 - u1 - u2) 
     .     - (t1 - t2 + u1 - u2)**2*(t1*t2 - u1*u2) 
     .     + s12**2*(t1**2 + 6._dp*t1*t2 + t2**2 - u1**2 - 6._dp*u1*u2 
     .     - u2**2) + s12*(-t1**3 - t2**3 + t2**2*(u1 - u2) 
     .     + t1**2*(3._dp*t2 - u1 + u2) + (u1 + u2)*(u1**2 - 4._dp*u1*u2 
     .     + u2**2) + t2*(-u1**2 - 2._dp*u1*u2 + u2**2) 
     .     + t1*(3._dp*t2**2 + u1**2 - 2._dp*u1*u2 - u2**2 
     .     + 2._dp*t2*(u1 + u2))))))



      tta1 = 
     .     - 0.125_dp*(-mz**2*(t1 + t2 + u1 + u2) 
     .     + 2._dp*s12*(s12 + t1 + t2 + u1 + u2))*
     .     (s12**2*(t1 + t2 - u1 - u2) 
     .     + s12*(t1 + t2 - u1 - u2)*(t1 + t2 + u1 + u2) 
     .     + (t1 + t2 + u1 + u2)*(t1*t2 - u1*u2))



      tta2 = 
     .     - 0.25_dp*mt**2*(s12*(4._dp*s12**5 
     .     + 4._dp*s12**4*(3._dp*t1 + 4._dp*(t2 + u1) + 3._dp*u2) 
     .     + 2._dp*s12**3*(t1 + t2 + u1 + u2)*(7._dp*t1 + 11._dp*(t2 + u1) 
     .     + 7._dp*u2) + s12*(t1 + t2 + u1 + u2)*(3._dp*t1**3 + 3._dp*t2**3 
     .     + 7._dp*t2**2*u1 + 7._dp*t2*u1**2 + 3._dp*u1**3 
     .     + (t2 + 3._dp*u1)*(5._dp*t2 + 3._dp*u1)*u2 
     .     + (5._dp*t2 + 11._dp*u1)*u2**2 + 3._dp*u2**3 
     .     + t1*(3._dp*t2 + u1 + 3._dp*u2)*(3._dp*t2 + 5._dp*u1 + 3._dp*u2) 
     .     + t1**2*(11._dp*t2 + 5._dp*u1 + 9._dp*u2)) + s12**2*(9._dp*t1**3 
     .     + 13._dp*t2**3 + (u1 + u2)**2*(13._dp*u1 + 9._dp*u2) 
     .     + t1**2*(31._dp*t2 + 27._dp*u1 + 25._dp*u2) + t2**2*(37._dp*u1 
     .     + 31._dp*u2) + t1*(35._dp*t2**2 + 68._dp*t2*u1 + 31._dp*u1**2 
     .     + 60._dp*(t2 + u1)*u2 + 25._dp*u2**2) + t2*(37._dp*u1**2 
     .     + 68._dp*u1*u2 + 27._dp*u2**2)) 
     .     + (t1 + t2 + u1 + u2)*(t1**3*(3._dp*t2 + 2._dp*u2) 
     .     + u1*u2*((t2 + u1)**2 + 2._dp*(2._dp*t2 + u1)*u2 + 3._dp*u2**2) 
     .     + t1**2*(2._dp*t2*(t2 + 2._dp*u1) + 3._dp*(2._dp*t2 + u1)*u2 
     .     + 4._dp*u2**2) + t1*(t2*(t2 + u1)**2 
     .     + 2._dp*(t2**2 + 4._dp*t2*u1 + u1**2)*u2 
     .     + 3._dp*(t2 + 2._dp*u1)*u2**2 + 2._dp*u2**3))) 
     .     + mz**2*(4._dp*s12**5 + 4._dp*s12**4*(2._dp*t1 + 3._dp*(t2 + u1) 
     .     + 2._dp*u2) + 2._dp*s12**3*(t1 + t2 + u1 + u2)*(3._dp*t1 
     .     + 5._dp*(t2 + u1) + 3._dp*u2) - (t1 + t2 + u1 + u2)*(t1*t2 
     .     - u1*u2)*(t1**2 + t2**2 - u1**2 - u2**2) 
     .     - s12*(t1 + t2 + u1 + u2)*(t1**3 + t2**3 + u1**3 + u1**2*u2 
     .     - u1*u2**2 + u2**3 - t2**2*(u1 + u2) - t2*(u1 + u2)**2 
     .     - t1**2*(t2 + u1 + 3._dp*u2) + t1*(t2**2 - u1**2 - 2._dp*u1*u2 
     .     - 3._dp*u2**2 - 2._dp*t2*(u1 + u2))) + s12**2*(t1**3 + t2**3 
     .     + t2**2*(9._dp*u1 + 7._dp*u2) 
     .     + t2*(u1 + u2)*(9._dp*u1 + 7._dp*u2) + t1**2*(7._dp*(t2 + u1) 
     .     + 9._dp*u2) + (u1 + u2)*(u1**2 + 6._dp*u1*u2 + u2**2) 
     .     + t1*(7._dp*t2**2 + 16._dp*t2*(u1 + u2) 
     .     + (u1 + u2)*(7._dp*u1 + 9._dp*u2)))))


      ttv1 = fac1*ttv1
      tta1 = fac2*tta1

      cs  = (Nc**2 -1._dp)

      qqb = 0._dp

      qqb = genfac*((ttv1+ttv2)*gvq*gvt + (tta1+tta2)*gaq*gat)

      do j=-nf,nf
         k=-j
         if((j == 0) .and. (k == 0)) then
            msq(j,k)=0._dp
         else if((j > 0) .and. (k < 0)) then
            msq(j,k) = qqb(j)
         else if((j < 0 ) .and. (k > 0)) then
            msq(j,k) = qqb(k)
         end if
      end do

      msq = msq*esq*gsq**2*cs/Nc**2

      end subroutine qqb_QQb_g_ew
