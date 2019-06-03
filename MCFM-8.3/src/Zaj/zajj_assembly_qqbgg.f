      subroutine msq_zqqbgamgg(p,qqb_gg,qbq_gg,gg_qqb,
     & qg_qg,qbg_qbg,gq_gq,gqb_gqb)
        use zajj_anomcoup_m
        use VVconfig_m
        implicit none
*****************************************************************
* return averaged matrix element squared for                    *
* 0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7) *
* coming from (0->q+qb+gam+g+g+lb+l) amplitude                  *
*****************************************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'zprods_decl.f'
        include 'anomcoup.f'
        include 'masses.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: qqb_gg(2),qbq_gg(2),gg_qqb(2)
        real(dp), intent(out) :: qg_qg(2),qbg_qbg(2),gq_gq(2),gqb_gqb(2)

        integer (kind(decayElAntiEl)) :: vDecay

        real(dp) :: pnt(mxpart,4)
        integer :: j
        complex(dp) :: a70h1(2,2,2,2,2,2),a70h3(2,2,2,2,2,2)
        complex(dp) :: anomZZ(2,2,2,2,2,2), anomZa(2,2,2,2,2,2)
        complex(dp) :: anomaZ(2,2,2,2,2,2)
        complex(dp) :: htZa, htZb, htGa, htGb

        if (anomtgc) then
          htZa = im*hitZ(1) + hitZ(3)
          htZb = im*hitZ(2) + hitZ(4)
          htGa = im*hitgam(1) + hitgam(3)
          htGb = im*hitgam(2) + hitgam(4)
        endif

        vDecay = decayChannel()

c-----convert to Nagy-Trocsanyi momentum convention
        pnt(1,:) = p(2,:)
        pnt(2,:) = p(5,:)
        pnt(3,:) = p(6,:)
        pnt(4,:) = p(7,:)
        pnt(5,:) = p(1,:)
        pnt(6,:) = p(4,:)
        pnt(7,:) = p(3,:)

        qqb_gg = 0._dp
        qbq_gg = 0._dp
        gg_qqb = 0._dp
        qg_qg = 0._dp
        qbg_qbg = 0._dp
        gq_gq = 0._dp
        gqb_gqb = 0._dp

        call spinoru(7,pnt,za,zb)


c----- In determining assignments, recall translation from NT to MCFM:
c-----  1_NT -> 2_MCFM
c-----  5_NT -> 1_MCFM
c-----  3_NT -> 6_MCFM
c-----  4_NT -> 7_MCFM

        a70h1 = 0._dp
        a70h3 = 0._dp
        anomZZ = 0._dp
        anomZa = 0._dp
        anomaZ = 0._dp

c----- qqb_gg
        call xzqqagg_qq(1,2,3,4,5,6,7,za,zb,a70h1)
        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call xzqqagg_ql(1,2,3,4,5,6,7,za,zb,a70h3)
        endif
        if (anomtgc) then
          call xzqqagg_anom(1,2,3,4,5,6,7,za,zb,amp_qqagg_anomZZ,htZa,htZb,anomZZ)
          call xzqqagg_anom(1,2,3,4,5,6,7,za,zb,amp_qqagg_anomZa,htGa,htGb,anomZa)
          call xzqqagg_anom(1,2,3,4,5,6,7,za,zb,amp_qqagg_anomaZ,htGa,htGb,anomaZ)
        endif
        do j=1,2
           call zajj_qqbgg_sq(j,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,qqb_gg(j))
        enddo

c----- qbq_gg
        call xzqqagg_qq(5,2,3,4,1,6,7,za,zb,a70h1)
        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call xzqqagg_ql(5,2,3,4,1,6,7,za,zb,a70h3)
        endif
        if (anomtgc) then
          call xzqqagg_anom(5,2,3,4,1,6,7,za,zb,amp_qqagg_anomZZ,htZa,htZb,anomZZ)
          call xzqqagg_anom(5,2,3,4,1,6,7,za,zb,amp_qqagg_anomZa,htGa,htGb,anomZa)
          call xzqqagg_anom(5,2,3,4,1,6,7,za,zb,amp_qqagg_anomaZ,htGa,htGb,anomaZ)
        endif
        do j=1,2
           call zajj_qqbgg_sq(j,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,qbq_gg(j))
        enddo

c----- gg_qqb
        call xzqqagg_qq(3,2,1,5,4,6,7,za,zb,a70h1)
        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call xzqqagg_ql(3,2,1,5,4,6,7,za,zb,a70h3)
        endif
        if (anomtgc) then
          call xzqqagg_anom(3,2,1,5,4,6,7,za,zb,amp_qqagg_anomZZ,htZa,htZb,anomZZ)
          call xzqqagg_anom(3,2,1,5,4,6,7,za,zb,amp_qqagg_anomZa,htGa,htGb,anomZa)
          call xzqqagg_anom(3,2,1,5,4,6,7,za,zb,amp_qqagg_anomaZ,htGa,htGb,anomaZ)
        endif
        do j=1,2
           call zajj_qqbgg_sq(j,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,gg_qqb(j))
        enddo

c----- qg_qg
        call xzqqagg_qq(3,2,1,4,5,6,7,za,zb,a70h1)
        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call xzqqagg_ql(3,2,1,4,5,6,7,za,zb,a70h3)
        endif
        if (anomtgc) then
          call xzqqagg_anom(3,2,1,4,5,6,7,za,zb,amp_qqagg_anomZZ,htZa,htZb,anomZZ)
          call xzqqagg_anom(3,2,1,4,5,6,7,za,zb,amp_qqagg_anomZa,htGa,htGb,anomZa)
          call xzqqagg_anom(3,2,1,4,5,6,7,za,zb,amp_qqagg_anomaZ,htGa,htGb,anomaZ)
        endif
        do j=1,2
           call zajj_qqbgg_sq(j,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,qg_qg(j))
        enddo

c----- qbg_qbg
        call xzqqagg_qq(5,2,1,4,3,6,7,za,zb,a70h1)
        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call xzqqagg_ql(5,2,1,4,3,6,7,za,zb,a70h3)
        endif
        if (anomtgc) then
          call xzqqagg_anom(5,2,1,4,3,6,7,za,zb,amp_qqagg_anomZZ,htZa,htZb,anomZZ)
          call xzqqagg_anom(5,2,1,4,3,6,7,za,zb,amp_qqagg_anomZa,htGa,htGb,anomZa)
          call xzqqagg_anom(5,2,1,4,3,6,7,za,zb,amp_qqagg_anomaZ,htGa,htGb,anomaZ)
        endif
        do j=1,2
           call zajj_qqbgg_sq(j,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,qbg_qbg(j))
        enddo

c----- gq_gq
        call xzqqagg_qq(3,2,5,4,1,6,7,za,zb,a70h1)
        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call xzqqagg_ql(3,2,5,4,1,6,7,za,zb,a70h3)
        endif
        if (anomtgc) then
          call xzqqagg_anom(3,2,5,4,1,6,7,za,zb,amp_qqagg_anomZZ,htZa,htZb,anomZZ)
          call xzqqagg_anom(3,2,5,4,1,6,7,za,zb,amp_qqagg_anomZa,htGa,htGb,anomZa)
          call xzqqagg_anom(3,2,5,4,1,6,7,za,zb,amp_qqagg_anomaZ,htGa,htGb,anomaZ)
        endif
        do j=1,2
           call zajj_qqbgg_sq(j,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,gq_gq(j))
        enddo

c----- gqb_gqb
        call xzqqagg_qq(1,2,5,4,3,6,7,za,zb,a70h1)
        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call xzqqagg_ql(1,2,5,4,3,6,7,za,zb,a70h3)
        endif
        if (anomtgc) then
          call xzqqagg_anom(1,2,5,4,3,6,7,za,zb,amp_qqagg_anomZZ,htZa,htZb,anomZZ)
          call xzqqagg_anom(1,2,5,4,3,6,7,za,zb,amp_qqagg_anomZa,htGa,htGb,anomZa)
          call xzqqagg_anom(1,2,5,4,3,6,7,za,zb,amp_qqagg_anomaZ,htGa,htGb,anomaZ)
        endif
        do j=1,2
           call zajj_qqbgg_sq(j,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,gqb_gqb(j))
        enddo

        qqb_gg(:) = qqb_gg(:)*half*aveqq
        qbq_gg(:) = qbq_gg(:)*half*aveqq
        gg_qqb(:) = gg_qqb(:)*avegg
        qg_qg(:) = qg_qg(:)*aveqg
        qbg_qbg(:) = qbg_qbg(:)*aveqg
        gq_gq(:) = gq_gq(:)*aveqg
        gqb_gqb(:) = gqb_gqb(:)*aveqg

      end

      subroutine zajj_qqbgg_sq(qi,za,zb,a70h1,a70h3,anomZZ,anomZa,anomaZ,msq)
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

        integer, intent(in) :: qi
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(in) :: a70h1(2,2,2,2,2,2), a70h3(2,2,2,2,2,2)
        complex(dp), intent(in) :: anomZZ(2,2,2,2,2,2), anomZa(2,2,2,2,2,2),
     &                             anomaZ(2,2,2,2,2,2)
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
                      vl = l1
                    else
                      vl = r1
                    endif


                    if (decayChannel() == decayQuarks) then
                      if (h5 == h1) then
                        ql = -abs(q1)
                      else
                        ql = abs(q1)
                      endif

                      qq = Q(qi)

                    else
                      if (h5 == h1) then
                        qq = -abs(q1)
                      else
                        qq = abs(q1)
                      endif

                      ql = Q(qi)
                    endif

                    amp(ic,h1,h2,h3,h4,h5) = 4*(sqrt(esq))**3*gsq* (
c XXX comment next two lines
     &                ql*(q1*Q(qi) + vl*vq*propQQ) * a70h1(ic,h1,h2,h3,h4,h5)
     &              + qq*(q1*Q(qi) + vl*vq*propQL) * a70h3(ic,h1,h2,h3,h4,h5)
     &              + vl*vq * propZ_s267 * propZ_s67 * anomZZ(ic,h1,h2,h3,h4,h5)
     &              + q1*vq * propZ_s267 * propA_s67 * anomZa(ic,h1,h2,h3,h4,h5)
     &              + Q(qi)*vl * propA_s267 * propZ_s67 * anomaZ(ic,h1,h2,h3,h4,h5)
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
        msq = (Nc**2-1)*Nc*msq/2

      end subroutine

      subroutine xzqqagg_qq(i1,i2,i3,i4,i5,i6,i7,za,zb,a70h)
        use zajj_treeamps_m, only: zajj_tree_qqgg_ppp,
     &          zajj_tree_qqgg_ppm, zajj_tree_qqgg_pmp,
     &          zajj_tree_qqgg_mpp
        implicit none
        include 'types.f'
        include 'mxpart.f'
      
************************************************************************
*     return the helicity amplitudes for                               *
*     0 ---> q(p1)+ph(p2)+g(p3)+g(p4)+qbar(p5)+lb(p6)+l(p7)            *
*     photon is coming from quark line                                 *
************************************************************************

      complex(dp):: amp_qqgga
      integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
      complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
      complex(dp), intent(out) :: a70h(2,2,2,2,2,2)

      a70h = 0._dp

      a70h(1, 2, 2,2,2, 1) = zajj_tree_qqgg_ppp(i1,i5,i4,i3,i2,i7,i6,za,zb)
      a70h(1, 2, 1,2,2, 1) = zajj_tree_qqgg_ppm(i1,i5,i4,i3,i2,i7,i6,za,zb)
      a70h(1, 2, 2,2,1, 1) = zajj_tree_qqgg_mpp(i1,i5,i4,i3,i2,i7,i6,za,zb)
      a70h(1, 2, 2,1,2, 1) = zajj_tree_qqgg_pmp(i1,i5,i4,i3,i2,i7,i6,za,zb)

      a70h(1, 2, 2,2,2, 2) = zajj_tree_qqgg_ppp(i1,i5,i4,i3,i2,i6,i7,za,zb)
      a70h(1, 2, 1,2,2, 2) = zajj_tree_qqgg_ppm(i1,i5,i4,i3,i2,i6,i7,za,zb)
      a70h(1, 2, 2,2,1, 2) = zajj_tree_qqgg_mpp(i1,i5,i4,i3,i2,i6,i7,za,zb)
      a70h(1, 2, 2,1,2, 2) = zajj_tree_qqgg_pmp(i1,i5,i4,i3,i2,i6,i7,za,zb)

      a70h(1, 1, 2,2,2, 1) = -zajj_tree_qqgg_ppp(i5,i1,i3,i4,i2,i7,i6,za,zb)
      a70h(1, 1, 1,2,2, 1) = -zajj_tree_qqgg_ppm(i5,i1,i3,i4,i2,i7,i6,za,zb)
      a70h(1, 1, 2,2,1, 1) = -zajj_tree_qqgg_pmp(i5,i1,i3,i4,i2,i7,i6,za,zb)
      a70h(1, 1, 2,1,2, 1) = -zajj_tree_qqgg_mpp(i5,i1,i3,i4,i2,i7,i6,za,zb)

      a70h(1, 1, 2,2,2, 2) = -zajj_tree_qqgg_ppp(i5,i1,i3,i4,i2,i6,i7,za,zb)
      a70h(1, 1, 1,2,2, 2) = -zajj_tree_qqgg_ppm(i5,i1,i3,i4,i2,i6,i7,za,zb)
      a70h(1, 1, 2,2,1, 2) = -zajj_tree_qqgg_pmp(i5,i1,i3,i4,i2,i6,i7,za,zb)
      a70h(1, 1, 2,1,2, 2) = -zajj_tree_qqgg_mpp(i5,i1,i3,i4,i2,i6,i7,za,zb)

      ! CC
      a70h(1, 1, 1,1,1, 2) = -zajj_tree_qqgg_ppp(i1,i5,i4,i3,i2,i7,i6,zb,za)
      a70h(1, 1, 2,1,1, 2) = -zajj_tree_qqgg_ppm(i1,i5,i4,i3,i2,i7,i6,zb,za)
      a70h(1, 1, 1,1,2, 2) = -zajj_tree_qqgg_mpp(i1,i5,i4,i3,i2,i7,i6,zb,za)
      a70h(1, 1, 1,2,1, 2) = -zajj_tree_qqgg_pmp(i1,i5,i4,i3,i2,i7,i6,zb,za)

      a70h(1, 1, 1,1,1, 1) = -zajj_tree_qqgg_ppp(i1,i5,i4,i3,i2,i6,i7,zb,za)
      a70h(1, 1, 2,1,1, 1) = -zajj_tree_qqgg_ppm(i1,i5,i4,i3,i2,i6,i7,zb,za)
      a70h(1, 1, 1,1,2, 1) = -zajj_tree_qqgg_mpp(i1,i5,i4,i3,i2,i6,i7,zb,za)
      a70h(1, 1, 1,2,1, 1) = -zajj_tree_qqgg_pmp(i1,i5,i4,i3,i2,i6,i7,zb,za)

      a70h(1, 2, 1,1,1, 2) = zajj_tree_qqgg_ppp(i5,i1,i3,i4,i2,i7,i6,zb,za)
      a70h(1, 2, 2,1,1, 2) = zajj_tree_qqgg_ppm(i5,i1,i3,i4,i2,i7,i6,zb,za)
      a70h(1, 2, 1,1,2, 2) = zajj_tree_qqgg_pmp(i5,i1,i3,i4,i2,i7,i6,zb,za)
      a70h(1, 2, 1,2,1, 2) = zajj_tree_qqgg_mpp(i5,i1,i3,i4,i2,i7,i6,zb,za)

      a70h(1, 2, 1,1,1, 1) = zajj_tree_qqgg_ppp(i5,i1,i3,i4,i2,i6,i7,zb,za)
      a70h(1, 2, 2,1,1, 1) = zajj_tree_qqgg_ppm(i5,i1,i3,i4,i2,i6,i7,zb,za)
      a70h(1, 2, 1,1,2, 1) = zajj_tree_qqgg_pmp(i5,i1,i3,i4,i2,i6,i7,zb,za)
      a70h(1, 2, 1,2,1, 1) = zajj_tree_qqgg_mpp(i5,i1,i3,i4,i2,i6,i7,zb,za)


c     ! second color ordering i3 <> i4
      a70h(2, 2, 2,2,2, 1) = zajj_tree_qqgg_ppp(i1,i5,i3,i4,i2,i7,i6,za,zb)
      a70h(2, 2, 1,2,2, 1) = zajj_tree_qqgg_ppm(i1,i5,i3,i4,i2,i7,i6,za,zb)
      a70h(2, 2, 2,2,1, 1) = zajj_tree_qqgg_pmp(i1,i5,i3,i4,i2,i7,i6,za,zb)
      a70h(2, 2, 2,1,2, 1) = zajj_tree_qqgg_mpp(i1,i5,i3,i4,i2,i7,i6,za,zb)

      a70h(2, 2, 2,2,2, 2) = zajj_tree_qqgg_ppp(i1,i5,i3,i4,i2,i6,i7,za,zb)
      a70h(2, 2, 1,2,2, 2) = zajj_tree_qqgg_ppm(i1,i5,i3,i4,i2,i6,i7,za,zb)
      a70h(2, 2, 2,2,1, 2) = zajj_tree_qqgg_pmp(i1,i5,i3,i4,i2,i6,i7,za,zb)
      a70h(2, 2, 2,1,2, 2) = zajj_tree_qqgg_mpp(i1,i5,i3,i4,i2,i6,i7,za,zb)

      a70h(2, 1, 2,2,2, 1) = -zajj_tree_qqgg_ppp(i5,i1,i4,i3,i2,i7,i6,za,zb)
      a70h(2, 1, 1,2,2, 1) = -zajj_tree_qqgg_ppm(i5,i1,i4,i3,i2,i7,i6,za,zb)
      a70h(2, 1, 2,2,1, 1) = -zajj_tree_qqgg_mpp(i5,i1,i4,i3,i2,i7,i6,za,zb)
      a70h(2, 1, 2,1,2, 1) = -zajj_tree_qqgg_pmp(i5,i1,i4,i3,i2,i7,i6,za,zb)

      a70h(2, 1, 2,2,2, 2) = -zajj_tree_qqgg_ppp(i5,i1,i4,i3,i2,i6,i7,za,zb)
      a70h(2, 1, 1,2,2, 2) = -zajj_tree_qqgg_ppm(i5,i1,i4,i3,i2,i6,i7,za,zb)
      a70h(2, 1, 2,2,1, 2) = -zajj_tree_qqgg_mpp(i5,i1,i4,i3,i2,i6,i7,za,zb)
      a70h(2, 1, 2,1,2, 2) = -zajj_tree_qqgg_pmp(i5,i1,i4,i3,i2,i6,i7,za,zb)

      ! CC
      a70h(2, 1, 1,1,1, 2) = -zajj_tree_qqgg_ppp(i1,i5,i3,i4,i2,i7,i6,zb,za)
      a70h(2, 1, 2,1,1, 2) = -zajj_tree_qqgg_ppm(i1,i5,i3,i4,i2,i7,i6,zb,za)
      a70h(2, 1, 1,1,2, 2) = -zajj_tree_qqgg_pmp(i1,i5,i3,i4,i2,i7,i6,zb,za)
      a70h(2, 1, 1,2,1, 2) = -zajj_tree_qqgg_mpp(i1,i5,i3,i4,i2,i7,i6,zb,za)
 
      a70h(2, 1, 1,1,1, 1) = -zajj_tree_qqgg_ppp(i1,i5,i3,i4,i2,i6,i7,zb,za)
      a70h(2, 1, 2,1,1, 1) = -zajj_tree_qqgg_ppm(i1,i5,i3,i4,i2,i6,i7,zb,za)
      a70h(2, 1, 1,1,2, 1) = -zajj_tree_qqgg_pmp(i1,i5,i3,i4,i2,i6,i7,zb,za)
      a70h(2, 1, 1,2,1, 1) = -zajj_tree_qqgg_mpp(i1,i5,i3,i4,i2,i6,i7,zb,za)

      a70h(2, 2, 1,1,1, 2) = zajj_tree_qqgg_ppp(i5,i1,i4,i3,i2,i7,i6,zb,za)
      a70h(2, 2, 2,1,1, 2) = zajj_tree_qqgg_ppm(i5,i1,i4,i3,i2,i7,i6,zb,za)
      a70h(2, 2, 1,1,2, 2) = zajj_tree_qqgg_mpp(i5,i1,i4,i3,i2,i7,i6,zb,za)
      a70h(2, 2, 1,2,1, 2) = zajj_tree_qqgg_pmp(i5,i1,i4,i3,i2,i7,i6,zb,za)
 
      a70h(2, 2, 1,1,1, 1) = zajj_tree_qqgg_ppp(i5,i1,i4,i3,i2,i6,i7,zb,za)
      a70h(2, 2, 2,1,1, 1) = zajj_tree_qqgg_ppm(i5,i1,i4,i3,i2,i6,i7,zb,za)
      a70h(2, 2, 1,1,2, 1) = zajj_tree_qqgg_mpp(i5,i1,i4,i3,i2,i6,i7,zb,za)
      a70h(2, 2, 1,2,1, 1) = zajj_tree_qqgg_pmp(i5,i1,i4,i3,i2,i6,i7,zb,za)
      end


************************************************************************
*     return the helicity amplitudes for                               *
*     0 ---> q(p1)+ph(p2)+g(p3)+g(p4)+qbar(p5)+lb(p6)+l(p7)            *
*     photon is coming from lepton line                                *
************************************************************************
      subroutine xzqqagg_ql(j1,j2,j3,j4,j5,j6,j7,za,zb,a70h3)
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'

        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(out) :: a70h3(2,2,2,2,2,2)

        integer :: lh,h2,h3,h4,hq
        integer :: hqc,lhc,h2c,h3c,h4c
        complex(dp) :: amp_qqagg_ql

        a70h3 = 0._dp

        hq = 2
        hqc = 1

        do lh=1,2; do h2=1,2; do h3=1,2; do h4=1,2
          a70h3(1,hq,h2,h3,h4,lh) =
     &          amp_qqagg_ql(j1,hq,j2,h2,j3,h3,j4,h4,j5,lh,j6,j7,za,zb)
          a70h3(2,hq,h2,h3,h4,lh) =
     &          amp_qqagg_ql(j1,hq,j2,h2,j4,h4,j3,h3,j5,lh,j6,j7,za,zb)

          lhc=mod(lh,2)+1
          h2c=mod(h2,2)+1
          h3c=mod(h3,2)+1
          h4c=mod(h4,2)+1

          a70h3(1,hqc,h2c,h3c,h4c,lhc) = conjg(a70h3(1,hq,h2,h3,h4,lh))
          a70h3(2,hqc,h2c,h3c,h4c,lhc) = conjg(a70h3(2,hq,h2,h3,h4,lh))
! manuallyfix signs due to crossing that are ruined by conjg relation above
          if ((j4 == 4) .and. (j3 /= 3)) then
            a70h3(1,hqc,h2c,h3c,h4c,lhc)=-a70h3(1,hqc,h2c,h3c,h4c,lhc)
            a70h3(2,hqc,h2c,h3c,h4c,lhc)=-a70h3(2,hqc,h2c,h3c,h4c,lhc)
          endif
        enddo; enddo; enddo; enddo

      end
 
