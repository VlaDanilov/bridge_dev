      subroutine qqb_QQb_mix_sft(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'sprods_com.f'
      include 'qcdcouple.f'
      include 'masses.f'
      integer j,k,ii,ff
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),pp(mxpart,4),
     .     msqb(-nf:nf,-nf:nf),ss(mxpart,mxpart)

      pp(1:4,:)=p(1:4,:)
      call qqb_QQb_mix(pp,msqb)

      call dotem(5,p,s)

      ss(:,:) = s(:,:)/2._dp

      msq = 0._dp

      do j = -nf,nf
         k = -j
         if (j .ne. 0) then
            msq(j,k) = + msq(j,k)
     .           + msqb(j,k)*(
     .           + (s(1,3)/(s(1,5)+s(3,5)))/s(1,5)
     .           - (s(2,3)/(s(2,5)+s(3,5)))/s(2,5)
     .           - (s(1,4)/(s(1,5)+s(4,5)))/s(1,5)
     .           + (s(2,4)/(s(2,5)+s(4,5)))/s(2,5)
     .           + (s(3,1)/(s(3,5)+s(1,5)) - mt**2/2._dp/s(3,5))/s(3,5)
     .           - (s(4,1)/(s(4,5)+s(1,5)) - mt**2/2._dp/s(4,5))/s(4,5)
     .           - (s(3,2)/(s(3,5)+s(2,5)) - mt**2/2._dp/s(3,5))/s(3,5)
     .           + (s(4,2)/(s(4,5)+s(2,5)) - mt**2/2._dp/s(4,5))/s(4,5)
!     .           + s(1,3)/s(1,5)/s(3,5)
!     .           - s(2,3)/s(2,5)/s(3,5)
!     .           - s(1,4)/s(1,5)/s(4,5)
!     .           + s(2,4)/s(2,5)/s(4,5)
     .           )    

c--- overall factor here is canonical one for eikonal;
c--- combined with 1/16/pi^2 from phase-space it gives (as/2/pi)    
         msq(j,k)=2._dp*gsq*(
     .           + s(1,3)/s(1,5)/s(3,5)
     .           - s(2,3)/s(2,5)/s(3,5)
     .           - s(1,4)/s(1,5)/s(4,5)
     .           + s(2,4)/s(2,5)/s(4,5))
         msq(j,k)=msq(j,k)*msqb(j,k)
         end if
      end do

      msq = 4._dp*msq

c--- soft limit does not affect the born assuming sum_{i=1,4}p_i = 0

      end subroutine qqb_QQb_mix_sft
            

