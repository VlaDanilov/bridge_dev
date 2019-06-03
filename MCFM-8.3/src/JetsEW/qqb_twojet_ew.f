      subroutine qqb_twojet_ew(p,msq)
c--- LO routine for twojet production via weak or weak-strong interaction
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'sprods_com.f'
      include 'msq_mix.f'
      include 'nflav.f'
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),ss,tt,uu,
     &     qa_iijj(nf,nf),qa_ijij(nf,nf),
     &     qq_ijij(nf,nf),aq_iijj(nf,nf),
     &     aq_ijij(nf,nf),aa_ijij(nf,nf)
      integer j,k

      call dotem(4,p,s)
      ss = s(1,2)
      tt = s(1,3)
      uu = s(2,3)

C --- LO EW Born - all the crossing relations
      call twojet_ew_tree(qa_iijj,ss,tt,uu)
      call twojet_ew_tree(qa_ijij,tt,ss,uu)
      call twojet_ew_tree(qq_ijij,tt,uu,ss)

c      call twojet_ew_tree(aq_iijj,ss,tt,uu)
c      call twojet_ew_tree(aq_ijij,tt,ss,uu)
c      call twojet_ew_tree(aa_ijij,tt,uu,ss)
c--- Use symmetry relations instead of calls, for speed
      aq_iijj(:,:)=qa_iijj(:,:)
      aq_ijij(:,:)=qa_ijij(:,:)
      aa_ijij(:,:)=qq_ijij(:,:)

      msq = 0._dp

      do j = -nflav,nflav
         do k = -nflav,nflav

C --- qa
            if((j > 0) .and. (k < 0)) then
               if(j == -k) then
                  msq(j,k) = qa_iijj(j,1) + qa_iijj(j,2) 
     &                     + qa_iijj(j,3) + qa_iijj(j,4)
     &                     + qa_iijj(j,5)
               else
                  msq(j,k) = qa_ijij(j,-k)
               end if

C --- qq
            else if((j > 0) .and. (k > 0)) then
               if(j == k) then
                  msq(j,k) = qq_ijij(j,k)*half
               else
                  msq(j,k) = qq_ijij(j,k)
               end if

C --- aq
            else if((j < 0) .and. (k > 0)) then
               if(k == -j) then
                  msq(j,k) = aq_iijj(k,1) + aq_iijj(k,2)
     &                     + aq_iijj(k,3) + aq_iijj(k,4)
     &                     + aq_iijj(k,5)
               else
                  msq(j,k) = aq_ijij(-j,k)
               end if

C --- aa
            else if((j < 0) .and. (k < 0)) then
               if(j == k) then
                  msq(j,k) = aa_ijij(-j,-k)*half
               else
                  msq(j,k) = aa_ijij(-j,-k)
               end if

            end if

         enddo
      enddo

      end subroutine qqb_twojet_ew

