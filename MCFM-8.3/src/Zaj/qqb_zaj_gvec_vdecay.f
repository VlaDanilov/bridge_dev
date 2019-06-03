      subroutine qqb_zaj_gvec_vdecay(p_mcfm,n,in,msq)
        implicit none
        include 'types.f'
C*********************************************************************** 
c     Matrix element for Z+gamma+jet production                        *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p6)              *
c     0 -> q(p1)+qbar(p2)+l(p3)+lb(p4)+gam(p5)+glu(p6)                 *
C*********************************************************************** 
      
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'masses.f'

        real(dp), intent(in) :: p_mcfm(mxpart,4), n(4)
        integer, intent(in) :: in
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zab(mxpart,mxpart),zba(mxpart,mxpart)

        complex(dp) :: aqqbn_qq(2,2,2), aqqbn_ql(2,2,2)
        complex(dp) :: aqbqn_qq(2,2,2), aqbqn_ql(2,2,2)

        real(dp) :: p_crossed(mxpart,4)

        real(dp) :: qqbn(2,2),qbqn(2,2)
        integer :: j,k
        integer, parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
        integer, parameter :: kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

        complex(dp) :: htZa,htZb,htGa,htGb

        integer :: j1,j2,j3
        integer :: flavor
        real(dp) :: s,t
        s(j1,j2) = real(za(j1,j2)*zb(j2,j1))
        t(j1,j2,j3) = s(j1,j2) + s(j2,j3) + s(j3,j1)

        qqbn = 0._dp
        qbqn = 0._dp

        !!! crossing
        if (in == 6) then
          p_crossed(1,:) = p_mcfm(4,:)
          p_crossed(2,:) = p_mcfm(3,:)
          p_crossed(3,:) = p_mcfm(2,:)
          p_crossed(4,:) = p_mcfm(1,:)
          p_crossed(5,:) = p_mcfm(5,:)
          p_crossed(6,:) = p_mcfm(6,:)
        elseif (in==4) then
          ! 6 and 4 exchanged here
          p_crossed(1,:) = p_mcfm(6,:)
          p_crossed(2,:) = p_mcfm(3,:)
          p_crossed(3,:) = p_mcfm(2,:)
          p_crossed(4,:) = p_mcfm(1,:)
          p_crossed(5,:) = p_mcfm(5,:)
          p_crossed(6,:) = p_mcfm(4,:)
        endif

        call spinoru(6,p_crossed,za,zb)
        call spinork(6,p_crossed,zab,zba,n)

        if (in==6) then
           call zajn_a60h(1,2,3,4,5,6,p_crossed,n,za,zb,zab,zba,aqbqn_qq,aqbqn_ql) !qbq
           call zajn_a60h(2,1,3,4,5,6,p_crossed,n,za,zb,zab,zba,aqqbn_qq,aqqbn_ql) !qqb
           do j=1,2; do flavor=1,2
              call zajn_m60sq_vdecay(flavor,j,1,2,3,4,5,6,za,zb,aqbqn_qq,aqbqn_ql,qbqn(flavor,j))
              call zajn_m60sq_vdecay(flavor,j,2,1,3,4,5,6,za,zb,aqqbn_qq,aqqbn_ql,qqbn(flavor,j))
           enddo; enddo
        elseif (in==4) then
           call zajn_a60h(1,2,3,4,5,6,p_crossed,n,za,zb,zab,zba,aqbqn_qq,aqbqn_ql) !qbq
           call zajn_a60h(2,1,3,4,5,6,p_crossed,n,za,zb,zab,zba,aqqbn_qq,aqqbn_ql) !qqb
           do j=1,2; do flavor=1,2
              call zajn_m60sq_vdecay(flavor,j,1,2,3,4,5,6,za,zb,aqbqn_qq,aqbqn_ql,qbqn(flavor,j))
              call zajn_m60sq_vdecay(flavor,j,2,1,3,4,5,6,za,zb,aqqbn_qq,aqqbn_ql,qqbn(flavor,j))
           enddo; enddo
        else
            write (*,*) "in = ", in
            call abort
        endif

        msq = 0._dp

        do j=-nf,nf
        do k=-nf,nf
            if ((j > 0) .and. (k == -j)) then
              msq(j,k) = 
     &            + n_light_down*aveqq*qqbn(1,jj(j))
     &            + n_light_up*aveqq*qqbn(2,jj(j))
            elseif ((j < 0) .and. (k == -j)) then
              msq(j,k) = 
     &            + n_light_down*aveqq*qbqn(1,kk(k))
     &            + n_light_up*aveqq*qbqn(2,kk(k))
            endif
        enddo
        enddo
      end

      subroutine zajn_m60sq_vdecay(flavor,qi,i1,i2,i3,i4,i5,i6,za,zb,amp_qq,amp_ql,msqn)
          use VVconfig_m
        implicit none
        include 'types.f'
***********************************************
* squared matrix element
* qqb_zag, for initial state flavor qi
* gluon line is contracted with a vector n(mu)
***********************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'ewcouple.f'
        include 'ewcharge.f'
        include 'qcdcouple.f'
        include 'zcouple.f'
        include 'masses.f'

        integer, intent(in) :: flavor,qi
        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: amp_qq(2,2,2), amp_ql(2,2,2)
        real(dp), intent(out) :: msqn

        integer :: h1,h2,h3
        complex(dp) :: ampsDress(2,2,2)
        real(dp) :: vq,vl,Ql,qq

        complex(dp):: propzQ,propzL
        complex(dp) :: propZ_s34, propZ_s345, propA_s34, propA_s345

        real(dp) :: s,t
        integer :: j1,j2,j3
        s(j1,j2) = real(za(j1,j2)*zb(j2,j1))
        t(j1,j2,j3) = s(j1,j2) + s(j2,j3) + s(j3,j1)

        propzQ=s(i3,i4)/(s(i3,i4)-zmass**2 + im*zwidth*zmass)
        propzL=t(i3,i4,i5)/(t(i3,i4,i5)-zmass**2 + im*zwidth*zmass)

        propZ_s34 = 1/(s(i3,i4) - zmass**2 + im*zwidth*zmass)
        propZ_s345 = 1/(t(i3,i4,i5) - zmass**2 + im*zwidth*zmass)
        propA_s34 = 1/s(i3,i4)
        propA_s345 = 1/t(i3,i4,i5)

        do h1=1,2
          do h2=1,2
            if (h2==1) then
              vq = l(qi)
            else
              vq = r(qi)
            endif

            do h3=1,2
              if (h3==1) then
                vl = l(flavor)
              else
                vl = r(flavor)
              endif

              if (h2==h3) then
                qq = +Q(flavor)
              else
                qq = -Q(flavor)
              endif

              Ql = Q(qi)


              ampsDress(h1,h2,h3) = twort2*sqrt(esq)**3*sqrt(gsq)*(
     &           amp_qq(h1,h2,h3) * qq*(q(flavor)*Q(qi) + vl*vq*propzQ)
     &         + amp_ql(h1,h2,h3) * Ql*(q(flavor)*Q(qi)*propA_s345 + vl*vq*propZ_s345)
     &        )
            enddo
          enddo
        enddo

        msqn = sum(abs(ampsDress)**2)*cA

      end
