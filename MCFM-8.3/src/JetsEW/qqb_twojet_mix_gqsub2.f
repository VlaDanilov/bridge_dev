      subroutine qqb_twojet_mix_gqsub2(p,msq)
c--- Routine that provides LO matrix elements for some of the subtraction
c--- terms with a (gluon,quark) initial state
c--- Differs from the normal LO routine in the following respects:
c---    1.  Only computes (q,qbar) and (qbar,q) matrix elements
c---    2.  For annihilation contribution, sum over final state
c---         includes a factor of two if identical initial and final states
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'sprods_com.f'
      include 'msq_mix.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),ss,tt,uu,
     &     qaii_jj(0:3,nf,nf),qaij_ij(0:3,nf,nf),
     &     qqij_ij(0:3,nf,nf),aqii_jj(0:3,nf,nf),
     &     aqij_ij(0:3,nf,nf),aaij_ij(0:3,nf,nf)
      integer j,k
      real(dp):: mss(2),mts,mst,mtt(2),sxs

      call dotem(4,p,s)
      ss = s(1,2)
      tt = s(1,3)
      uu = s(2,3)

      call qqb_twojet_ii_jj_mix(qaii_jj,ss,tt,uu)
      call qqb_twojet_ii_jj_mix(qaij_ij,tt,ss,uu)
      
      aqii_jj(:,:,:)=qaii_jj(:,:,:)
      aqij_ij(:,:,:)=qaij_ij(:,:,:)

c-- for now, no return value; real purpose of this routine is
c-- to fill the common block msq_mix
      msq=0._dp

      msq_mix = 0._dp
      do j = fn,nf
         do k = fn,nf

C --- qa
            if((j > 0) .and. (k < 0)) then
               if(j == -k) then
                  msq_mix(:,j,k) = qaii_jj(:,j,1) + qaii_jj(:,j,2) 
     .                 + qaii_jj(:,j,3) + qaii_jj(:,j,4)
     .                 - qaii_jj(:,j,j) + qaii_jj(:,j,j)*0.5_dp
               else
                  msq_mix(:,j,k) = qaij_ij(:,j,-k)
               end if

C --- qq
c            else if((j > 0) .and. (k > 0)) then
c               if(j == k)  then
c                  msq_mix(:,j,k) = qqij_ij(:,j,k)*half
c               else
c                  msq_mix(:,j,k) = qqij_ij(:,j,k)
c               end if

C --- aq
            else if((j < 0) .and. (k > 0)) then
               if(k == -j) then
                  msq_mix(:,j,k) = aqii_jj(:,k,1) + aqii_jj(:,k,2) 
     .                 + aqii_jj(:,k,3) + aqii_jj(:,k,4)
     .                 - qaii_jj(:,k,k) + qaii_jj(:,k,k)*0.5_dp
               else
                  msq_mix(:,j,k) = aqij_ij(:,-j,k)
               end if

C --- aa
c            else if((j < 0) .and. (k < 0)) then
c               if(j == k) then
c                  msq_mix(:,j,k) = aaij_ij(:,-j,-k)*half
c               else
c                  msq_mix(:,j,k) = aaij_ij(:,-j,-k)
c               end if
            end if
         end do
      end do
      
C--- multiply by factor of N*CF=4, which is the natural factor
c--- for the sxt interference terms that have a collinear singularity
      msq_mix=msq_mix*(xn*Cf)

      return
      end
