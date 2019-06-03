      subroutine qqb_zaj_g_vdecay(p,msq)
        implicit none
        include 'types.f'
*******************************************************************
*  return matrix element squared for                              *
*  0 -> f(p1) + f(p2) + d(p3) + db(p4) + gam(p5) + f(p6) + f(p7)  *
*******************************************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        real(dp) :: msq2(-nf:nf,-nf:nf)

        real(dp) :: qqb_gg(2,2),qbq_gg(2,2)
        integer :: j,k
        integer, parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
        integer, parameter :: kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

        real(dp) :: qqb_ij(2,2,2), qbq_ij(2,2,2)
        real(dp) :: qqb_ii(2,2), qbq_ii(2,2)

        call msq_zqqbgamgg_vdecay(p,qqb_gg,qbq_gg)
        call msq_zqqbQQbgam_vdecay(p,qqb_ij,qbq_ij,qqb_ii,qbq_ii)

        msq = 0._dp

        do j=-nf,nf
        do k=-nf,nf
            if ((j > 0) .and. (k == -j)) then
              msq(j,k) = 0._dp
     &                  + n_light_down*qqb_gg(1,jj(j))
     &                  + n_light_up*qqb_gg(2,jj(j))
            elseif ((j < 0) .and. (k == -j)) then
              msq(j,k) =  0._dp
     &                  + n_light_down*qbq_gg(1,kk(k))
     &                  + n_light_up*qbq_gg(2,kk(k))
            endif
        enddo
        enddo

        ! four quark contributions.
        ! number of contributions adds up to nf*nf = (n_light_up + n_light_down)**2

        msq2 = 0._dp

        ! same flavor contributions
        msq2(-1,1) = (msq2(-1,1) + n_light_up*qbq_ii(1, 2) + n_light_down*qbq_ii(1, 1))
        msq2(1,-1) = (msq2(1,-1) + n_light_up*qqb_ii(1, 2) + n_light_down*qqb_ii(1, 1))

        msq2(-2,2) = (msq2(-2,2) + n_light_up*qbq_ii(2, 2) + n_light_down*qbq_ii(2, 1))
        msq2(2,-2) = (msq2(2,-2) + n_light_up*qqb_ii(2, 2) + n_light_down*qqb_ii(2, 1))

        ! different flavor contributions
        msq2(-1,1) = msq2(-1,1)
     &      + n_light_down*(n_light_down-1)*qbq_ij(1, 1,1)
     &      + n_light_up*(n_light_up-1)*qbq_ij(1, 2,2)
     &      + n_light_up*n_light_down*(qbq_ij(1, 1,2) + qbq_ij(1, 2,1))

        msq2(1,-1) = msq2(1,-1)
     &      + n_light_down*(n_light_down-1)*qqb_ij(1, 1,1)
     &      + n_light_up*(n_light_up-1)*qqb_ij(1, 2,2)
     &      + n_light_up*n_light_down*(qqb_ij(1, 1,2) + qqb_ij(1, 2,1))


        msq2(-2,2) = msq2(-2,2)
     &      + n_light_down*(n_light_down-1)*qbq_ij(2, 1,1)
     &      + n_light_up*(n_light_up-1)*qbq_ij(2, 2,2)
     &      + n_light_up*n_light_down*(qbq_ij(2, 1,2) + qbq_ij(2, 2,1))

        msq2(2,-2) = msq2(2,-2)
     &      + n_light_down*(n_light_down-1)*qqb_ij(2, 1,1)
     &      + n_light_up*(n_light_up-1)*qqb_ij(2, 2,2)
     &      + n_light_up*n_light_down*(qqb_ij(2, 1,2) + qqb_ij(2, 2,1))


        msq2(-3,3) = msq2(-1,1)
        msq2(-5,5) = msq2(-1,1)

        msq2(-4,4) = msq2(-2,2)

        msq2(3,-3) = msq2(1,-1)
        msq2(5,-5) = msq2(1,-1)

        msq2(4,-4) = msq2(2,-2)

        msq = msq + msq2

      end

