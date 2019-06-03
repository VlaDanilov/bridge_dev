      subroutine msq_zqqbgamgg_vdecay(p, qqb_gg,qbq_gg)
        use VVconfig_m
        implicit none
*****************************************************************
* return averaged matrix element squared for                    *
* 0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7) *
* coming from (0->q+qb+gam+g+g+lb+l) amplitude                  *
*****************************************************************
      
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'nf.f'
        include 'zprods_decl.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: qqb_gg(2,2),qbq_gg(2,2)

        integer (kind(decayElAntiEl)) :: vDecay

        real(dp) :: pnt(mxpart,4)
        integer :: j, flavor
        complex(dp) :: a70h1(2,2,2,2,2,2),a70h3(2,2,2,2,2,2)

        if (decayChannel() /= decayQuarks) call abort

c-----convert to Nagy-Trocsanyi momentum convention
c--   and switch mcfm momenta 2 with 3 and 1 with 4 for crossing!
        pnt(1,:) = p(3,:)
        pnt(2,:) = p(5,:)
        pnt(3,:) = p(6,:)
        pnt(4,:) = p(7,:)
        pnt(5,:) = p(4,:)
        pnt(6,:) = p(1,:)
        pnt(7,:) = p(2,:)

        qqb_gg = 0._dp
        qbq_gg = 0._dp

        call spinoru(7,pnt,za,zb)

        a70h1 = 0._dp
        a70h3 = 0._dp

c----- qqb_gg
        call xzqqagg_qq(1,2,3,4,5,6,7,za,zb,a70h1)
        call xzqqagg_ql(1,2,3,4,5,6,7,za,zb,a70h3)
        do j=1,2; do flavor=1,2
           call zajj_qqbgg_sq_vdecay(flavor,j,za,zb,a70h1,a70h3,qqb_gg(flavor,j))
        enddo; enddo

c----- qbq_gg
        call xzqqagg_qq(1,2,3,4,5,7,6,za,zb,a70h1)
        call xzqqagg_ql(1,2,3,4,5,7,6,za,zb,a70h3)
        do j=1,2; do flavor=1,2
           call zajj_qqbgg_sq_vdecay(flavor,j,za,zb,a70h1,a70h3,qbq_gg(flavor,j))
        enddo; enddo

        qqb_gg(:,:) = qqb_gg(:,:)*half*aveqq
        qbq_gg(:,:) = qbq_gg(:,:)*half*aveqq
      end
      
      subroutine zajj_qqbgg_sq_vdecay(flavor,qi,za,zb,a70h1,a70h3,msq)
        use VVconfig_m
        implicit none
        include 'types.f'
********************************************************************
* 0 -> q(-p1) + qb(-p5) + a(p2) + g(p3) + g(p4) + lb(p6) + l(p7)
* return matrix element squared, for initial flavor qi
* given the helicity amplitudes from each channel
********************************************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'ewcouple.f'
        include 'ewcharge.f'
        include 'qcdcouple.f'
        include 'zcouple.f'
        include 'masses.f'

        integer, intent(in) :: flavor,qi
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: a70h1(2,2,2,2,2,2), a70h3(2,2,2,2,2,2)
        real(dp), intent(out) :: msq

        integer:: ic,h1,h2,h3,h4,h5
        complex(dp) :: amp(2,2,2,2,2,2)
        complex(dp):: propQQ,propQL
        complex(dp) :: propZ_s67, propZ_s267, propA_s67, propA_s267

        real(dp) :: vq, vl, qq, ql

        integer :: i1,i2,i3
        real(dp) :: s,t
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        propQQ = s(6,7)/(s(6,7)-zmass**2 + im*zwidth*zmass)
        propQL = t(2,6,7)/(t(2,6,7)-zmass**2 + im*zwidth*zmass)

        propZ_s67 = 1/(s(6,7) - zmass**2 + im*zwidth*zmass)
        propZ_s267 = 1/(t(2,6,7) - zmass**2 + im*zwidth*zmass)

        propA_s67 = 1/s(6,7)
        propA_s267 = 1/t(2,6,7)

        do ic=1,2
          do h1=1,2
            if (h1==1) then
              vq = l(qi)
            else
              vq = r(qi)
            endif

            do h2=1,2
              do h3=1,2
                do h4=1,2
                  do h5=1,2
                    if (h5==1) then
                      vl = l(flavor)
                    else
                      vl = r(flavor)
                    endif

                    if (h5 == h1) then
                      ql = Q(flavor)
                    else
                      ql = -Q(flavor)
                    endif

                    qq = Q(qi)

                    amp(ic,h1,h2,h3,h4,h5) = 4*(sqrt(esq))**3*gsq* (
c XXX comment next two lines
     &                ql*(q(flavor)*Q(qi) + vl*vq*propQQ) * a70h1(ic,h1,h2,h3,h4,h5)
     &              + qq*(q(flavor)*Q(qi) + vl*vq*propQL) * a70h3(ic,h1,h2,h3,h4,h5)
     &              )
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        msq =     sum((abs(amp(1,:,:,:,:,:)))**2)
     &         +  sum((abs(amp(2,:,:,:,:,:)))**2)
     &         -  sum((abs(amp(1,:,:,:,:,:) + amp(2,:,:,:,:,:)))**2)/Nc**2
        msq = (Nc**2-1)*Nc*msq/2*cA

      end subroutine
