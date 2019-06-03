      subroutine qqb_zaj_gvec(p,n,in,msq)
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
        include 'anomcoup.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4), n(4)
        integer, intent(in) :: in
        real(dp), intent(out) :: msq(-nf:nf,-nf:nf)

        complex(dp) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp) :: zab(mxpart,mxpart),zba(mxpart,mxpart)

        complex(dp) :: aqqbn_qq(2,2,2), aqqbn_ql(2,2,2)
        complex(dp) :: aqbqn_qq(2,2,2), aqbqn_ql(2,2,2)
        complex(dp) :: aqgn_qq(2,2,2), aqgn_ql(2,2,2)
        complex(dp) :: aqbgn_qq(2,2,2), aqbgn_ql(2,2,2)
        complex(dp) :: agqbn_qq(2,2,2), agqbn_ql(2,2,2)
        complex(dp) :: agqn_qq(2,2,2), agqn_ql(2,2,2)

        real(dp) :: qqbn(2),qbqn(2),qgn(2),qbgn(2),gqn(2),gqbn(2)
        integer :: j,k
        integer, parameter :: jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
        integer, parameter :: kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

        complex(dp) :: aqbqn_anomZZ(2,2,2), aqbqn_anomZa(2,2,2), aqbqn_anomaZ(2,2,2)
        complex(dp) :: aqqbn_anomZZ(2,2,2), aqqbn_anomZa(2,2,2), aqqbn_anomaZ(2,2,2)
        complex(dp) :: aqbgn_anomZZ(2,2,2), aqbgn_anomZa(2,2,2), aqbgn_anomaZ(2,2,2)
        complex(dp) :: aqgn_anomZZ(2,2,2), aqgn_anomZa(2,2,2), aqgn_anomaZ(2,2,2)
        complex(dp) :: agqbn_anomZZ(2,2,2), agqbn_anomZa(2,2,2), agqbn_anomaZ(2,2,2)
        complex(dp) :: agqn_anomZZ(2,2,2), agqn_anomZa(2,2,2), agqn_anomaZ(2,2,2)

        complex(dp) :: htZa,htZb,htGa,htGb

        integer :: j1,j2,j3
        real(dp) :: s,t
        s(j1,j2) = real(za(j1,j2)*zb(j2,j1))
        t(j1,j2,j3) = s(j1,j2) + s(j2,j3) + s(j3,j1)

        aqbqn_anomZZ = 0._dp
        aqqbn_anomZZ = 0._dp
        aqbgn_anomZZ = 0._dp
        aqgn_anomZZ = 0._dp
        agqbn_anomZZ = 0._dp
        agqn_anomZZ = 0._dp

        aqbqn_anomZa = 0._dp
        aqqbn_anomZa = 0._dp
        aqbgn_anomZa = 0._dp
        aqgn_anomZa = 0._dp
        agqbn_anomZa = 0._dp
        agqn_anomZa = 0._dp

        aqbqn_anomaZ = 0._dp
        aqqbn_anomaZ = 0._dp
        aqbgn_anomaZ = 0._dp
        aqgn_anomaZ = 0._dp
        agqbn_anomaZ = 0._dp
        agqn_anomaZ = 0._dp

        qqbn = 0._dp
        qbqn = 0._dp
        qgn = 0._dp
        qbgn = 0._dp
        gqn = 0._dp
        gqbn = 0._dp

        call spinoru(6,p,za,zb)
        call spinork(6,p,zab,zba,n)

        if (anomtgc) then
          htZa = im*hitZ(1) + hitZ(3)
          htZb = im*hitZ(2) + hitZ(4)
          htGa = im*hitgam(1) + hitgam(3)
          htGb = im*hitgam(2) + hitgam(4)
        endif


        if (in==6) then
           call zajn_a60h(1,2,3,4,5,6,p,n,za,zb,zab,zba,aqbqn_qq,aqbqn_ql) !qbq
           call zajn_a60h(2,1,3,4,5,6,p,n,za,zb,zab,zba,aqqbn_qq,aqqbn_ql) !qqb
           if (anomtgc) then
             call zajn_anom(1,2,3,4,5,6,p,n,za,zb,zab,zba,htZa,htZb,htGa,htGb,
     &                        aqbqn_anomZZ,aqbqn_anomZa,aqbqn_anomaZ)
             call zajn_anom(2,1,3,4,5,6,p,n,za,zb,zab,zba,htZa,htZb,htGa,htGb,
     &                        aqqbn_anomZZ,aqqbn_anomZa,aqqbn_anomaZ)
           endif
           do j=1,2
              call zajn_m60sq(j,za,zb,aqqbn_qq,aqqbn_ql,aqqbn_anomZZ,aqqbn_anomZa,aqqbn_anomaZ,qqbn(j))
              call zajn_m60sq(j,za,zb,aqbqn_qq,aqbqn_ql,aqbqn_anomZZ,aqbqn_anomZa,aqbqn_anomaZ,qbqn(j))
           enddo
        elseif (in==2) then
           call zajn_a60h(1,6,3,4,5,2,p,n,za,zb,zab,zba,aqbgn_qq,aqbgn_ql)  !qbg
           call zajn_a60h(6,1,3,4,5,2,p,n,za,zb,zab,zba,aqgn_qq,aqgn_ql) !qg
           if (anomtgc) then
             call zajn_anom(1,6,3,4,5,2,p,n,za,zb,zab,zba,htZa,htZb,htGa,htGb,
     &                        aqbgn_anomZZ,aqbgn_anomZa,aqbgn_anomaZ)
             call zajn_anom(6,1,3,4,5,2,p,n,za,zb,zab,zba,htZa,htZb,htGa,htGb,
     &                        aqgn_anomZZ,aqgn_anomZa,aqgn_anomaZ)
           endif
           do j=1,2
              call zajn_m60sq(j,za,zb,aqgn_qq,aqgn_ql,aqgn_anomZZ,aqgn_anomZa,aqgn_anomaZ,qgn(j))
              call zajn_m60sq(j,za,zb,aqbgn_qq,aqbgn_ql,aqbgn_anomZZ,aqbgn_anomZa,aqbgn_anomaZ,qbgn(j))
           enddo
        elseif (in==1) then
           call zajn_a60h(2,6,3,4,5,1,p,n,za,zb,zab,zba,agqbn_qq,agqbn_ql) !gqb
           call zajn_a60h(6,2,3,4,5,1,p,n,za,zb,zab,zba,agqn_qq,agqn_ql)  !gq
           if (anomtgc) then
             call zajn_anom(2,6,3,4,5,1,p,n,za,zb,zab,zba,htZa,htZb,htGa,htGb,
     &                        agqbn_anomZZ,agqbn_anomZa,agqbn_anomaZ)
             call zajn_anom(6,2,3,4,5,1,p,n,za,zb,zab,zba,htZa,htZb,htGa,htGb,
     &                        agqn_anomZZ,agqn_anomZa,agqn_anomaZ)
           endif
           do j=1,2
              call zajn_m60sq(j,za,zb,agqn_qq,agqn_ql,agqn_anomZZ,agqn_anomZa,agqn_anomaZ,gqn(j))
              call zajn_m60sq(j,za,zb,agqbn_qq,agqbn_ql,agqbn_anomZZ,agqbn_anomZa,agqbn_anomaZ,gqbn(j))
           enddo
        endif

        msq = 0._dp

        do j=-nf,nf
        do k=-nf,nf
            if ((j == 0) .and. (k == 0)) then
              msq(j,k)=zip
            elseif ((j == 0) .and. (k < 0)) then
              msq(j,k)=aveqg*gqbn(-kk(k))
            elseif ((j == 0) .and. (k > 0)) then
              msq(j,k)=aveqg*gqn(kk(k))
            elseif ((j > 0) .and. (k == -j)) then
              msq(j,k)=aveqq*qqbn(jj(j))
            elseif ((j < 0) .and. (k == -j)) then
              msq(j,k)=aveqq*qbqn(kk(k))
            elseif ((j > 0) .and. (k == 0)) then
              msq(j,k)=aveqg*qgn(jj(j))
            elseif ((j < 0) .and. (k == 0)) then
              msq(j,k)=aveqg*qbgn(-jj(j))
            else
              msq(j,k)=zip
            endif
        enddo
        enddo
      end

      subroutine zajn_m60sq(qi,za,zb,amp_qq,amp_ql,anomZZ,anomZa,anomaZ,msqn)
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

        integer, intent(in) :: qi
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: amp_qq(2,2,2), amp_ql(2,2,2)
        complex(dp), intent(in) :: anomZZ(2,2,2), anomZa(2,2,2), anomaZ(2,2,2)
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

        propzQ=s(3,4)/(s(3,4)-zmass**2 + im*zwidth*zmass)
        propzL=t(3,4,5)/(t(3,4,5)-zmass**2 + im*zwidth*zmass)

        propZ_s34 = 1/(s(3,4) - zmass**2 + im*zwidth*zmass)
        propZ_s345 = 1/(t(3,4,5) - zmass**2 + im*zwidth*zmass)
        propA_s34 = 1/s(3,4)
        propA_s345 = 1/t(3,4,5)

        do h1=1,2
          do h2=1,2
            if (h2==1) then
              vq = l(qi)
            else
              vq = r(qi)
            endif

            do h3=1,2
              if (h3==1) then
                vl = l1
              else
                vl = r1
              endif

              if (decayChannel() == decayQuarks) then
                if (h2==h3) then
                  qq = -abs(q1)
                else
                  qq = abs(q1)
                endif
                Ql = Q(qi)
              else
                if (h2==h3) then
                  Ql = -abs(q1)
                else
                  Ql = abs(q1)
                endif
                qq = Q(qi)
              endif


              
              ampsDress(h1,h2,h3) = twort2*sqrt(esq)**3*sqrt(gsq)*(
c XXX comment next two lines
     &           amp_qq(h1,h2,h3) * qq*(q1*Q(qi) + vl*vq*propzQ)
     &         + amp_ql(h1,h2,h3) * Ql*(q1*Q(qi)*propA_s345 + vl*vq*propZ_s345)
     &         + anomZZ(h1,h2,h3) * vl*vq * propZ_s345 * propZ_s34
     &         + anomZa(h1,h2,h3) * q1*vq * propZ_s345 * propA_s34
     &         + anomaZ(h1,h2,h3) * vl*Q(qi) * propA_s345 * propZ_s34
     &        )
            enddo
          enddo
        enddo

        msqn = sum(abs(ampsDress)**2)

      end

      subroutine zajn_anom(i1,i2,i3,i4,i5,i6,p,n,za,zb,zab,zba,
     &       htZa,htZb,htGa,htGb,anomZZ,anomZa,anomaZ)
        implicit none
        include 'types.f'
        include 'mxpart.f'

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        real(dp), intent(in) :: p(mxpart,4), n(4)
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: zab(mxpart,mxpart), zba(mxpart,mxpart)
        complex(dp), intent(out) :: anomZZ(2,2,2), anomZa(2,2,2), anomaZ(2,2,2)
        complex(dp), intent(in) :: htZa,htZb,htGa,htGb

        complex(dp) :: zaj_gvec_anomZZ_photMinus
        complex(dp) :: zaj_gvec_anomZa_photMinus
        complex(dp) :: zaj_gvec_anomaZ_photMinus
        real(dp) :: isgn

        call checkndotp(p,n,i6)

        !helicity label ordering: hgamma, hquark, hlepton

        anomZZ(1,1,1) = zaj_gvec_anomZZ_photMinus(i1,i2,i3,i4,i5,i6,
     &                    za,zb,zab,zba,htZa,htZb)
        anomZZ(1,1,2) = zaj_gvec_anomZZ_photMinus(i1,i2,i4,i3,i5,i6,
     &                    za,zb,zab,zba,htZa,htZb) * isgn(i6)
        anomZZ(1,2,1) = zaj_gvec_anomZZ_photMinus(i2,i1,i3,i4,i5,i6,
     &                    za,zb,zab,zba,htZa,htZb)
        anomZZ(1,2,2) = zaj_gvec_anomZZ_photMinus(i2,i1,i4,i3,i5,i6,
     &                    za,zb,zab,zba,htZa,htZb) * isgn(i6)

        anomZZ(2,1:2,1:2) = conjg(anomZZ(1,2:1:-1,2:1:-1))

c       anomZZ(2,2,2) = zaj_gvec_anomZZ_photMinus(i1,i2,i3,i4,i5,i6,
c    &                    zb,za,zba,zab,conjg(htZa),conjg(htZb)) * (-isgn(i6))
c       anomZZ(2,2,1) = zaj_gvec_anomZZ_photMinus(i1,i2,i4,i3,i5,i6,
c    &                    zb,za,zba,zab,conjg(htZa),conjg(htZb)) * (-1)
c       anomZZ(2,1,2) = zaj_gvec_anomZZ_photMinus(i2,i1,i3,i4,i5,i6,
c    &                    zb,za,zba,zab,conjg(htZa),conjg(htZb)) * (-isgn(i6))
c       anomZZ(2,1,1) = zaj_gvec_anomZZ_photMinus(i2,i1,i4,i3,i5,i6,
c    &                    zb,za,zba,zab,conjg(htZa),conjg(htZb)) * (-1)

        ! anomZa

        anomZa(1,1,1) = zaj_gvec_anomZa_photMinus(i1,i2,i3,i4,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb)
        anomZa(1,1,2) = zaj_gvec_anomZa_photMinus(i1,i2,i4,i3,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb) * isgn(i6)
        anomZa(1,2,1) = zaj_gvec_anomZa_photMinus(i2,i1,i3,i4,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb)
        anomZa(1,2,2) = zaj_gvec_anomZa_photMinus(i2,i1,i4,i3,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb) * isgn(i6)

        anomZa(2,1:2,1:2) = conjg(anomZa(1,2:1:-1,2:1:-1))

c       anomZa(2,2,2) = zaj_gvec_anomZa_photMinus(i1,i2,i3,i4,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-isgn(i6))
c       anomZa(2,2,1) = zaj_gvec_anomZa_photMinus(i1,i2,i4,i3,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-1)
c       anomZa(2,1,2) = zaj_gvec_anomZa_photMinus(i2,i1,i3,i4,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-isgn(i6))
c       anomZa(2,1,1) = zaj_gvec_anomZa_photMinus(i2,i1,i4,i3,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-1)

        ! anomaZ

        anomaZ(1,1,1) = zaj_gvec_anomaZ_photMinus(i1,i2,i3,i4,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb)
        anomaZ(1,1,2) = zaj_gvec_anomaZ_photMinus(i1,i2,i4,i3,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb) * isgn(i6)
        anomaZ(1,2,1) = zaj_gvec_anomaZ_photMinus(i2,i1,i3,i4,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb)
        anomaZ(1,2,2) = zaj_gvec_anomaZ_photMinus(i2,i1,i4,i3,i5,i6,
     &                    za,zb,zab,zba,htGa,htGb) * isgn(i6)

        anomaZ(2,1:2,1:2) = conjg(anomaZ(1,2:1:-1,2:1:-1))

c       anomaZ(2,2,2) = zaj_gvec_anomaZ_photMinus(i1,i2,i3,i4,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-isgn(i6))
c       anomaZ(2,2,1) = zaj_gvec_anomaZ_photMinus(i1,i2,i4,i3,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-1)
c       anomaZ(2,1,2) = zaj_gvec_anomaZ_photMinus(i2,i1,i3,i4,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-isgn(i6))
c       anomaZ(2,1,1) = zaj_gvec_anomaZ_photMinus(i2,i1,i4,i3,i5,i6,
c    &                    zb,za,zba,zab,conjg(htGa),conjg(htGb)) * (-1)

      end subroutine
      
      subroutine zajn_a60h(j1,j2,j3,j4,j5,j6,p,n,za,zb,zab,zba,a6nh_qq,a6nh_ql)
        implicit none
        include 'types.f'
****************************************************************
*  Color ordered amplitudes for:
*  0 -> q(p1) + qb(p2) + l(p3) + lb(p4) + gam(p5) + glu(p6)
*  where line 6 is contracted with the vector n(mu)
****************************************************************
c-----the order of momentum in the argument
c-----(q,qb,l,lb,ph,g,p,n,za,zb,zab,zba,a6nh)
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'cplx.h'

        integer, intent(in) :: j1,j2,j3,j4,j5,j6
        real(dp), intent(in) :: p(mxpart,4), n(4)
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: zab(mxpart,mxpart),zba(mxpart,mxpart)
        complex(dp), intent(out) :: a6nh_qq(2,2,2), a6nh_ql(2,2,2)

        real(dp):: nDp5,isgn
        complex(dp):: qcdabn(2,2,2),qcdban(2,2,2)
        !complex(dp) :: a6nql(2,2,2)

        complex(dp) :: zaj_gvec_photMinus

        real(dp) :: s,t
        s(j1,j2) = real(za(j1,j2)*zb(j2,j1))
        t(j1,j2,j3) = s(j1,j2) + s(j2,j3) + s(j3,j1)

        nDp5=n(4)*p(j5,4)-n(3)*p(j5,3)-n(2)*p(j5,2)-n(1)*p(j5,1)
        call checkndotp(p,n,j6)
        call subqcdn(j2,j1,j3,j4,j5,j6,nDp5,za,zb,zab,zba,qcdabn,qcdban)

        ! all signs adjusted such that amplitudes with opposite photon
        ! helicities are obtained by plain complex conjugation

        a6nh_qq(2,2,2) = qcdabn(2,2,2) + qcdban(2,2,2)
        a6nh_qq(1,2,2) = qcdabn(1,2,2) + qcdban(1,2,2)
        a6nh_qq(2,2,1) = qcdabn(2,2,1) + qcdban(2,2,1)
        a6nh_qq(1,2,1) = qcdabn(1,2,1) + qcdban(1,2,1)
        a6nh_qq(2,1,1) = qcdabn(2,1,1) + qcdban(2,1,1)
        a6nh_qq(1,1,1) = qcdabn(1,1,1) + qcdban(1,1,1)
        a6nh_qq(2,1,2) = qcdabn(2,1,2) + qcdban(2,1,2)
        a6nh_qq(1,1,2) = qcdabn(1,1,2) + qcdban(1,1,2)

        a6nh_ql(1,1,1) = zaj_gvec_photMinus(j1,j2,j3,j4,j5,j6,za,zb,zab,zba)
        a6nh_ql(1,1,2) = isgn(j6) * zaj_gvec_photMinus(j1,j2,j4,j3,j5,j6,za,zb,zab,zba)
        a6nh_ql(1,2,1) = zaj_gvec_photMinus(j2,j1,j3,j4,j5,j6,za,zb,zab,zba)
        a6nh_ql(1,2,2) = isgn(j6) * zaj_gvec_photMinus(j2,j1,j4,j3,j5,j6,za,zb,zab,zba)

        a6nh_ql(2,1:2,1:2) = conjg(a6nh_ql(1,2:1:-1,2:1:-1))

c       a6nh_ql(2,2,2) = -isgn(j6) * zaj_gvec_photMinus(j1,j2,j3,j4,j5,j6,zb,za,zba,zab)
c       a6nh_ql(2,2,1) = -zaj_gvec_photMinus(j1,j2,j4,j3,j5,j6,zb,za,zba,zab)
c       a6nh_ql(2,1,2) = -isgn(j6) * zaj_gvec_photMinus(j2,j1,j3,j4,j5,j6,zb,za,zba,zab)
c       a6nh_ql(2,1,1) = -zaj_gvec_photMinus(j2,j1,j4,j3,j5,j6,zb,za,zba,zab)

      end

      function zaj_gvec_photPlus(i1,i2,i3,i4,i5,i6,za,zb,zab,zba)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex(dp) :: zaj_gvec_photPlus

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: zab(mxpart,mxpart), zba(mxpart,mxpart)

        complex(dp) :: s
        s(i1,i2) = za(i1,i2)*zb(i2,i1)

        zaj_gvec_photPlus = 
     -  2*(s(i2,i6)*(za(i1,i3)*zab(i1,i1) - za(i3,i6)*zab(i1,i6))*
     -     (s(i4,i5)*za(i3,i4)*zb(i4,i2)*zb(i5,i3) + 
     -       s(i3,i5)*za(i4,i5)*zb(i5,i2)*zb(i5,i4)) + 
     -    s(i1,i6)*za(i1,i3)*
     -     (s(i4,i5)*za(i3,i4)*zb(i5,i3)*
     -        (-(zab(i2,i2)*zb(i4,i2)) + zab(i6,i2)*zb(i6,i4)) + 
     -       s(i3,i5)*za(i4,i5)*zb(i5,i4)*
     -        (-(zab(i2,i2)*zb(i5,i2)) + zab(i6,i2)*zb(i6,i5))))/
     -  (s(i1,i6)*s(i2,i6)*s(i3,i5)*s(i4,i5)*za(i4,i5))

      end function

      function zaj_gvec_photMinus(i1,i2,i3,i4,i5,i6,za,zb,zab,zba)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex(dp) :: zaj_gvec_photMinus

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: zab(mxpart,mxpart), zba(mxpart,mxpart)

        complex(dp) :: s
        s(i1,i2) = za(i1,i2)*zb(i2,i1)

        zaj_gvec_photMinus = 
     -  (-2*za(i3,i5)*(s(i2,i6)*zb(i4,i2)*
     -       (-(za(i1,i3)*zab(i1,i1)*zb(i4,i3)) + 
     -         za(i3,i6)*zab(i1,i6)*zb(i4,i3) + 
     -         (za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*zb(i5,i4))
     -       + s(i1,i6)*(za(i1,i3)*zb(i4,i3) - za(i1,i5)*zb(i5,i4))*
     -       (zab(i2,i2)*zb(i4,i2) - zab(i6,i2)*zb(i6,i4))))/
     -  (s(i1,i6)*s(i2,i6)*s(i3,i5)*zb(i5,i4))
      end function

      function zaj_gvec_anomZZ_photMinus(i1,i2,i3,i4,i5,i6,za,zb,zab,zba,ha,hb)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex(dp) :: zaj_gvec_anomZZ_photMinus

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: zab(mxpart,mxpart), zba(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s
        s(i1,i2) = za(i1,i2)*zb(i2,i1)

        zaj_gvec_anomZZ_photMinus =
     -  (ha*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (-2*s(i2,i6)*(za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*
     -        zb(i4,i2) + 2*s(i1,i6)*za(i1,i5)*zab(i2,i2)*zb(i4,i2) - 
     -       2*s(i1,i6)*za(i1,i5)*zab(i6,i2)*zb(i6,i4)))/
     -   (2._dp*s(i1,i6)*s(i2,i6)) + 
     -  (hb*(s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i2,i6)*s(i3,i5)*
     -        (za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*zb(i4,i2) + 
     -       s(i2,i6)*s(i4,i5)*
     -        (za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*zb(i4,i2) - 
     -       s(i1,i6)*s(i3,i5)*za(i1,i5)*zab(i2,i2)*zb(i4,i2) - 
     -       s(i1,i6)*s(i4,i5)*za(i1,i5)*zab(i2,i2)*zb(i4,i2) - 
     -       s(i2,i6)*za(i3,i5)*
     -        (za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*zb(i4,i3)*
     -        zb(i5,i2) + s(i1,i6)*za(i1,i5)*za(i3,i5)*zab(i2,i2)*
     -        zb(i4,i3)*zb(i5,i2) + 
     -       s(i1,i6)*s(i3,i5)*za(i1,i5)*zab(i6,i2)*zb(i6,i4) + 
     -       s(i1,i6)*s(i4,i5)*za(i1,i5)*zab(i6,i2)*zb(i6,i4) - 
     -       s(i1,i6)*za(i1,i5)*za(i3,i5)*zab(i6,i2)*zb(i4,i3)*zb(i6,i5)
     -       ))/(2._dp*s(i1,i6)*s(i2,i6))

      end function

      function zaj_gvec_anomZa_photMinus(i1,i2,i3,i4,i5,i6,za,zb,zab,zba,ha,hb)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex(dp) :: zaj_gvec_anomZa_photMinus

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: zab(mxpart,mxpart), zba(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s
        s(i1,i2) = za(i1,i2)*zb(i2,i1)

        zaj_gvec_anomZa_photMinus =
     -  (ha*s(i3,i4)*za(i3,i5)*
     -     (2*s(i2,i6)*za(i1,i5)*zab(i1,i1)*zb(i4,i2) - 
     -       2*s(i2,i6)*za(i5,i6)*zab(i1,i6)*zb(i4,i2) - 
     -       2*s(i1,i6)*za(i1,i5)*zab(i2,i2)*zb(i4,i2) + 
     -       2*s(i1,i6)*za(i1,i5)*zab(i6,i2)*zb(i6,i4)))/
     -   (2._dp*s(i1,i6)*s(i2,i6)) + 
     -  (hb*s(i3,i4)*za(i3,i5)*
     -     (s(i2,i6)*za(i3,i5)*
     -        (-(za(i1,i3)*zab(i1,i1)*zb(i3,i2)) + 
     -          za(i3,i6)*zab(i1,i6)*zb(i3,i2) + 
     -          (-(za(i1,i4)*zab(i1,i1)) + za(i4,i6)*zab(i1,i6))*
     -           zb(i4,i2))*zb(i4,i3) + 
     -       s(i1,i6)*za(i3,i5)*zb(i4,i3)*
     -        (za(i1,i3)*(zab(i2,i2)*zb(i3,i2) - 
     -             zab(i6,i2)*zb(i6,i3)) + 
     -          za(i1,i4)*(zab(i2,i2)*zb(i4,i2) - zab(i6,i2)*zb(i6,i4)))
     -       ))/(2._dp*s(i1,i6)*s(i2,i6))

      end function

      function zaj_gvec_anomaZ_photMinus(i1,i2,i3,i4,i5,i6,za,zb,zab,zba,ha,hb)
        implicit none
        include 'types.f'
        include 'mxpart.f'
        complex(dp) :: zaj_gvec_anomaZ_photMinus

        integer, intent(in) :: i1,i2,i3,i4,i5,i6
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: zab(mxpart,mxpart), zba(mxpart,mxpart)
        complex(dp), intent(in) :: ha,hb

        complex(dp) :: s
        s(i1,i2) = za(i1,i2)*zb(i2,i1)

        zaj_gvec_anomaZ_photMinus =
     -  (ha*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (-2*s(i2,i6)*(za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*
     -        zb(i4,i2) + 2*s(i1,i6)*za(i1,i5)*zab(i2,i2)*zb(i4,i2) - 
     -       2*s(i1,i6)*za(i1,i5)*zab(i6,i2)*zb(i6,i4)))/
     -   (2._dp*s(i1,i6)*s(i2,i6)) + 
     -  (hb*(s(i3,i4) + s(i3,i5) + s(i4,i5))*za(i3,i5)*
     -     (s(i2,i6)*s(i3,i5)*
     -        (za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*zb(i4,i2) + 
     -       s(i2,i6)*s(i4,i5)*
     -        (za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*zb(i4,i2) - 
     -       s(i1,i6)*s(i3,i5)*za(i1,i5)*zab(i2,i2)*zb(i4,i2) - 
     -       s(i1,i6)*s(i4,i5)*za(i1,i5)*zab(i2,i2)*zb(i4,i2) - 
     -       s(i2,i6)*za(i3,i5)*
     -        (za(i1,i5)*zab(i1,i1) - za(i5,i6)*zab(i1,i6))*zb(i4,i3)*
     -        zb(i5,i2) + s(i1,i6)*za(i1,i5)*za(i3,i5)*zab(i2,i2)*
     -        zb(i4,i3)*zb(i5,i2) + 
     -       s(i1,i6)*s(i3,i5)*za(i1,i5)*zab(i6,i2)*zb(i6,i4) + 
     -       s(i1,i6)*s(i4,i5)*za(i1,i5)*zab(i6,i2)*zb(i6,i4) - 
     -       s(i1,i6)*za(i1,i5)*za(i3,i5)*zab(i6,i2)*zb(i4,i3)*zb(i6,i5)
     -       ))/(2._dp*s(i1,i6)*s(i2,i6))

      end function

      function isgn(j1)
        implicit none
        include 'types.f'

        real(dp) :: isgn
        integer :: j1

        if (j1<=2) then
           isgn=-1._dp
        else
           isgn=1._dp
        endif
      end

