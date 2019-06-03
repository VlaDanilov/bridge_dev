      subroutine qqb_zaj_g(p,msq)
        implicit none
        include 'types.f'
*******************************************************************
*  Original author: H. Hartanto                                            *
*                                                                 *
*  return matrix element squared for                              *
*  0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7)  *
*******************************************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        real(dp) :: msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
        real(dp) :: msq3(-nf:nf,-nf:nf)
        real(dp) :: qqb_gg(2),qbq_gg(2),gg_qqb(2)
        real(dp) :: qg_qg(2),qbg_qbg(2),gq_gq(2),gqb_gqb(2)
        real(dp) :: iib_jjb(2,2),ibi_jjb(2,2)
        real(dp) :: iib_iib(2),ibi_iib(2),ii_ii(2),ibib_ibib(2)
        real(dp) :: qiqj(2,2),qbiqbj(2,2),qiqbj(2,2),qbiqj(2,2)
        integer :: j,k
        integer, parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
        integer, parameter :: kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

        call msq_zqqbgamgg(p,qqb_gg,qbq_gg,gg_qqb,
     &         qg_qg,qbg_qbg,gq_gq,gqb_gqb)
        call msq_zqqbQQbgam(p,iib_jjb,ibi_jjb,
     &         iib_iib,ibi_iib,ii_ii,ibib_ibib,qiqj,qbiqbj,qiqbj,qbiqj)

        msq = 0._dp
        msq1 = 0._dp
        msq2 = 0._dp
        msq3 = 0._dp

        do j=-nf,nf
        do k=-nf,nf
            if ((j == 0) .and. (k == 0)) then
              msq1(j,k)=(two*gg_qqb(2)+three*gg_qqb(1))
            elseif ((j == 0) .and. (k < 0)) then
              msq1(j,k)=gqb_gqb(-kk(k))
            elseif ((j == 0) .and. (k > 0)) then
              msq1(j,k)=gq_gq(kk(k))
            elseif ((j > 0) .and. (k == -j)) then
              msq1(j,k)=qqb_gg(jj(j))
            elseif ((j < 0) .and. (k == -j)) then
              msq1(j,k)=qbq_gg(kk(k))
            elseif ((j > 0) .and. (k == 0)) then
              msq1(j,k)=qg_qg(jj(j))
            elseif ((j < 0) .and. (k == 0)) then
              msq1(j,k)=qbg_qbg(-jj(j))
            endif

            if ((j>0).and.(k>0).and.(j/=k)) then
              msq3(j,k)=qiqj(jj(j),kk(k))
            elseif ((j<0).and.(k<0).and.(j/=k)) then
              msq3(j,k)=qbiqbj(-jj(j),-kk(k))
            elseif ((j>0).and.(k<0).and.(j/=-k)) then
              msq3(j,k)=qiqbj(jj(j),-kk(k))
            elseif ((j<0).and.(k>0).and.(j/=-k)) then
              msq3(j,k)=qbiqj(-jj(j),kk(k))
            endif
        enddo
        enddo

        msq2(2,-2) = iib_iib(2)+iib_jjb(2,2)+3._dp*iib_jjb(2,1)
        msq2(4,-4) = msq2(2,-2)
        msq2(1,-1) = iib_iib(1)+2._dp*iib_jjb(1,2)+2._dp*iib_jjb(1,1)
        msq2(3,-3) = msq2(1,-1)
        msq2(5,-5) = msq2(1,-1)

        msq2(-2,2) = ibi_iib(2)+ibi_jjb(2,2)+3._dp*ibi_jjb(2,1)
        msq2(-4,4) = msq2(-2,2)
        msq2(-1,1) = ibi_iib(1)+2._dp*ibi_jjb(1,2)+2._dp*ibi_jjb(1,1)
        msq2(-3,3) = msq2(-1,1)
        msq2(-5,5) = msq2(-1,1)

        msq2(2,2) = ii_ii(2)
        msq2(4,4) = msq2(2,2)
        msq2(1,1) = ii_ii(1)
        msq2(3,3) = msq2(1,1)
        msq2(5,5) = msq2(1,1)

        msq2(-2,-2) = ibib_ibib(2)
        msq2(-4,-4) = msq2(-2,-2)
        msq2(-1,-1) = ibib_ibib(1)
        msq2(-3,-3) = msq2(-1,-1)
        msq2(-5,-5) = msq2(-1,-1)

        msq(:,:) = msq1(:,:) + msq2(:,:) + msq3(:,:)

      end

