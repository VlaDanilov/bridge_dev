      subroutine msq_zqqbQQbgam(p,iib_jjb,ibi_jjb,
     & iib_iib,ibi_iib,ii_ii,ibib_ibib,qiqj,qbiqbj,qiqbj,qbiqj)
        implicit none
        include 'types.f'
*****************************************************************
* return matrix element squared for                             *
* 0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7) *
* coming from (0->q+qb+Q+Qb+gam+lb+l) amplitude                 *
*****************************************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'mxpart.f'
        include 'zprods_decl.f'
        include 'sprods_decl.f'

        real(dp), intent(in) :: p(mxpart,4)
        real(dp), intent(out) :: iib_jjb(2,2),ibi_jjb(2,2)
        real(dp), intent(out) :: iib_iib(2),ibi_iib(2),ii_ii(2),ibib_ibib(2)
        real(dp), intent(out) :: qiqj(2,2),qbiqbj(2,2), qiqbj(2,2), qbiqj(2,2)

        real(dp) :: pnt(mxpart,4)

        ! convert to Nagy-Trocsanyi momentum convention
        pnt(1,:) = p(2,:)
        pnt(2,:) = p(1,:)
        pnt(3,:) = p(6,:)
        pnt(4,:) = p(7,:)
        pnt(5,:) = p(5,:)
        pnt(6,:) = p(4,:)
        pnt(7,:) = p(3,:)

        iib_jjb = 0._dp
        ibi_jjb = 0._dp
        iib_iib = 0._dp
        ibi_iib = 0._dp
        ii_ii = 0._dp
        ibib_ibib = 0._dp

        qiqj = 0._dp
        qbiqbj = 0._dp
        qiqbj = 0._dp
        qbiqj = 0._dp

        call spinoru(7,pnt,za,zb)

        call zqqbQQba_msqij(1,2,3,4,5,6,7,za,zb,iib_jjb) 
        call zqqbQQba_msqij(2,1,3,4,5,6,7,za,zb,ibi_jjb) 
        call zqqbQQba_msqii(1,2,3,4,5,6,7,za,zb,iib_iib) 
        call zqqbQQba_msqii(2,1,3,4,5,6,7,za,zb,ibi_iib) 
        call zqqbQQba_msqii(3,1,4,2,5,6,7,za,zb,ii_ii) 
        call zqqbQQba_msqii(1,3,2,4,5,6,7,za,zb,ibib_ibib)

        ! In determining assignments, recall translation from NT to MCFM:
        ! 1_NT -> 2_MCFM
        ! 2_NT -> 1_MCFM
        ! 3_NT -> 6_MCFM
        ! 4_NT -> 7_MCFM

        call zqqbQQba_msqij(3,2,4,1,5,6,7,za,zb,qiqj)        ! q + Q -> q + Q
        call zqqbQQba_msqij(2,3,1,4,5,6,7,za,zb,qbiqbj)      ! qb + Qb -> qb + Qb
        call zqqbQQba_msqij(3,2,1,4,5,6,7,za,zb,qiqbj)       ! q + Qb -> q + Qb
        call zqqbQQba_msqij(2,4,3,1,5,6,7,za,zb,qbiqj)       ! qb + Q -> Q + qb

        iib_jjb(:,:) = iib_jjb(:,:)*aveqq
        ibi_jjb(:,:) = ibi_jjb(:,:)*aveqq

        iib_iib(:) = iib_iib(:)*aveqq
        ibi_iib(:) = ibi_iib(:)*aveqq
        ii_ii(:) = ii_ii(:)*aveqq*half
        ibib_ibib(:) = ibib_ibib(:)*aveqq*half

        qiqj(:,:) = qiqj(:,:)*aveqq
        qbiqbj(:,:) = qbiqbj(:,:)*aveqq
        qiqbj(:,:) = qiqbj(:,:)*aveqq
        qbiqj(:,:) = qbiqj(:,:)*aveqq
      end

      subroutine zqqbQQba_msqij(j1,j2,j3,j4,j5,j6,j7,za,zb,msq)
        use zajj_anomcoup_m
        use VVconfig_m
        implicit none
********************************************************************
* return squared matrix element for                                *
* 0 -> q(p1)+qb(p2)+Q(p3)+Qb(p4)+gamma(p5)+lb(p6)+l(p7)            *
* q and Q have different flavour                                   *
* above is Nagy-Trocsanyi momentum assignment                      *
* NO average, NO identical particle factor included                *                                        
********************************************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'cplx.h'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'zcouple.f'
        include 'masses.f'
        include 'ewcharge.f'
        include 'anomcoup.f'

        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        real(dp), intent(out) :: msq(2,2)

        integer :: qi,qj
        real(dp) :: ix,Ql
        integer :: hq,Qh,lh,ha
        real(dp) :: vq, Qv, vl
        complex(dp) :: propQQ,propQL
        complex(dp) :: m70hA(2,2, 2,2,2,2),m70hB(2,2, 2,2,2,2)
        complex(dp) :: m70Anom(2,2, 2,2,2,2)
        complex(dp) :: Ai(2,2,2,2,2,2),Bi(2,2,2,2,2,2)
        complex(dp) :: Aii(2,2,2,2),Bii(2,2,2,2)

        integer (kind(decayElAntiEl)) :: vDecay
        
        complex(dp) :: htZa, htZb, htGa, htGb

        complex(dp) :: anomZZ_p1(2,2,2,2), anomZZ_p2(2,2,2,2)
        complex(dp) :: anomZa_p1(2,2,2,2), anomZa_p2(2,2,2,2)
        complex(dp) :: anomaZ_p1(2,2,2,2), anomaZ_p2(2,2,2,2)
        complex(dp) :: propZ_s67, propZ_s567
        complex(dp) :: propA_s67, propA_s567

        real(dp) :: s,t
        s(j1,j2) = real(za(j1,j2)*zb(j2,j1))
        t(j1,j2,j3) = s(j1,j2) + s(j2,j3) + s(j3,j1)

        Ai = 0._dp
        Bi = 0._dp
        Aii = 0._dp
        Bii = 0._dp

        vDecay = decayChannel()

        call amp_zqqQQa_qq_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Ai,Bi)

        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call amp_zqqQQa_ql_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Aii)
          call amp_zqqQQa_ql_new(j3,j4,j1,j2,j5,j6,j7,za,zb,Bii)
        endif

        if (anomtgc) then
          htZa = im*hitZ(1) + hitZ(3)
          htZb = im*hitZ(2) + hitZ(4)
          htGa = im*hitgam(1) + hitgam(3)
          htGb = im*hitgam(2) + hitgam(4)

          call amp_zqqQQa_anomZZ(j1,j2,j3,j4,j5,j6,j7,za,zb,htZa,htZb,anomZZ_p1)
          call amp_zqqQQa_anomZZ(j3,j4,j1,j2,j5,j6,j7,za,zb,htZa,htZb,anomZZ_p2)

          call amp_zqqQQa_anomZa(j1,j2,j3,j4,j5,j6,j7,za,zb,htGa,htGb,anomZa_p1)
          call amp_zqqQQa_anomZa(j3,j4,j1,j2,j5,j6,j7,za,zb,htGa,htGb,anomZa_p2)

          call amp_zqqQQa_anomaZ(j1,j2,j3,j4,j5,j6,j7,za,zb,htGa,htGb,anomaZ_p1)
          call amp_zqqQQa_anomaZ(j3,j4,j1,j2,j5,j6,j7,za,zb,htGa,htGb,anomaZ_p2)

        endif

        propQQ = s(6,7)/(s(6,7)-zmass**2 + im*zwidth*zmass)
        propQL = t(5,6,7)/(t(5,6,7)-zmass**2 + im*zwidth*zmass)

        propZ_s67 = 1/(s(6,7) - zmass**2 + im*zmass*zwidth)
        propZ_s567 = 1/(t(5,6,7) - zmass**2 + im*zmass*zwidth)
        propA_s67 = 1/s(6,7)
        propA_s567 = 1/t(5,6,7)

        msq = 0._dp

        do qi=1,2; do qj=1,2
          m70hA = 0._dp
          m70hB = 0._dp
          m70Anom = 0._dp

          do hq=1,2; 
            if (hq == 1) then
              vq = l(qi)
            else
              vq = r(qi)
            endif
          
            do Qh=1,2; 
              if (Qh == 1) then
                Qv = l(qj)
              else
                Qv = r(qj)
              endif

              if (hq /= Qh) then
                 ix=-1._dp
              else
                 ix=1._dp
              endif
        
              do ha=1,2; do lh=1,2
                if (lh == 1) then
                  vl = l1
                else
                  vl = r1
                endif

                if (hq == lh) then
                   Ql = -abs(q1)
                else
                   Ql = abs(q1)
                endif

                m70hA(qi,qj, hq,Qh,ha,lh) = 
     &               (q1*Q(qi) + vq*vl*propQQ)*Ai(qi,qj,hq,Qh,ha,lh)
     &             + (q1*Q(qj) + Qv*vl*propQQ)*Bi(qi,qj,Qh,hq,ha,lh)
                m70hB(qi,qj, hq,Qh,ha,lh) = 
     &             + Ql*( (q1*Q(qi)*propA_s567 + vq*vl*propZ_s567)*Aii(hq,Qh,ha,lh)
     &               + ix*(q1*Q(qj)*propA_s567 + Qv*vl*propZ_s567)*Bii(Qh,hq,ha,lh))

                if (anomtgc) then
                  m70Anom(qi,qj, hq,Qh,ha,lh) = m70Anom(qi,qj, hq,Qh,ha,lh)
     &               + ( vq*vl*propZ_s567*propZ_s67 * anomZZ_p1(hq,Qh,ha,lh)
     &              + Qv*vl*propZ_s567*propZ_s67 * anomZZ_p2(Qh,hq,ha,lh))

                  m70Anom(qi,qj, hq,Qh,ha,lh) = m70Anom(qi,qj, hq,Qh,ha,lh)
     &               + ( q1*vq*propZ_s567*propA_s67 * anomZa_p1(hq,Qh,ha,lh)
     &              + q1*Qv*propZ_s567*propA_s67 * anomZa_p2(Qh,hq,ha,lh))
                  m70Anom(qi,qj, hq,Qh,ha,lh) = m70Anom(qi,qj, hq,Qh,ha,lh)
     &               + ( Q(qi)*vl*propA_s567*propZ_s67 * anomaZ_p1(hq,Qh,ha,lh)
     &              + Q(qj)*vl*propA_s567*propZ_s67 * anomaZ_p2(Qh,hq,ha,lh))
                endif

              enddo
            enddo
          enddo; enddo

          ! XXX uncomment next two lines
          !m70hA = 0._dp
          !m70hB = 0._dp

          msq(qi,qj) = sum(abs(m70hA(qi,qj,:,:,:,:) +
     &            m70hB(qi,qj,:,:,:,:) + m70Anom(qi,qj,:,:,:,:))**2)
          ! include color factor, gauge couplings, but no averaging
          msq(qi,qj)=8._dp*8._dp*esq**3*gsq**2*msq(qi,qj)

        enddo; enddo

      end


      subroutine zqqbQQba_msqii(j1,j2,j3,j4,j5,j6,j7,za,zb,msq)
        use zajj_anomcoup_m
        use VVconfig_m
        implicit none
********************************************************************
* return squared matrix element for                                *
* 0 -> q(p1)+qb(p2)+q(p3)+qb(p4)+gamma(p5)+lb(p6)+l(p7)            *
* q and Q have the same flavour                                    *           
* above is Nagy-Trocsanyi momentum assignment                      *
* NO average, NO identical particle factor included                *                                        
********************************************************************
      
        include 'constants.f'
        include 'nf.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'zcouple.f'
        include 'masses.f'
        include 'ewcharge.f'
        include 'anomcoup.f'

        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        real(dp), intent(out) :: msq(2)

        integer :: qi,qj
        real(dp) :: ix,Ql
        real(dp) :: ampsq(2,2,2,2)
        integer :: hq,Qh,lh,ha
        complex(dp) :: m70hA(2,2,2,2,2),m70hB(2,2,2,2,2),m70Anom(2,2,2,2,2)
        complex(dp) :: Ai(2,2,2,2,2,2),Bi(2,2,2,2,2,2)
        complex(dp) :: Ci(2,2,2,2,2,2),Di(2,2,2,2,2,2)
        complex(dp) :: Aii(2,2,2,2),Bii(2,2,2,2)
        complex(dp) :: Cii(2,2,2,2),Dii(2,2,2,2)

        integer (kind(decayElAntiEl)) :: vDecay

        complex(dp) :: htZa, htZb, htGa, htGb
        real(dp) :: vq, Qv, vl

        complex(dp) :: anomZZ_c1_p1(2,2,2,2), anomZZ_c1_p2(2,2,2,2)
        complex(dp) :: anomZa_c1_p1(2,2,2,2), anomZa_c1_p2(2,2,2,2)
        complex(dp) :: anomaZ_c1_p1(2,2,2,2), anomaZ_c1_p2(2,2,2,2)

        complex(dp) :: anomZZ_c2_p1(2,2,2,2), anomZZ_c2_p2(2,2,2,2)
        complex(dp) :: anomZa_c2_p1(2,2,2,2), anomZa_c2_p2(2,2,2,2)
        complex(dp) :: anomaZ_c2_p1(2,2,2,2), anomaZ_c2_p2(2,2,2,2)

        complex(dp) :: propQQ,propQL
        complex(dp) :: propZ_s567, propZ_s67
        complex(dp) :: propA_s567, propA_s67

        real(dp) :: s,t
        s(j1,j2) = real(za(j1,j2)*zb(j2,j1))
        t(j1,j2,j3) = s(j1,j2) + s(j2,j3) + s(j3,j1)

        Ai = 0._dp
        Bi = 0._dp
        Ci = 0._dp
        Di = 0._dp
        Aii = 0._dp
        Bii = 0._dp
        Cii = 0._dp
        Dii = 0._dp

        vDecay = decayChannel()

        call amp_zqqQQa_qq_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Ai,Bi)
        call amp_zqqQQa_qq_new(j1,j4,j3,j2,j5,j6,j7,za,zb,Ci,Di)

        if (vDecay == decayElAntiEl .or. vDecay == decayQuarks) then
          call amp_zqqQQa_ql_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Aii)
          call amp_zqqQQa_ql_new(j3,j4,j1,j2,j5,j6,j7,za,zb,Bii)

          call amp_zqqQQa_ql_new(j1,j4,j3,j2,j5,j6,j7,za,zb,Cii)
          call amp_zqqQQa_ql_new(j3,j2,j1,j4,j5,j6,j7,za,zb,Dii)
        endif


        if (anomtgc) then
          htZa = im*hitZ(1) + hitZ(3)
          htZb = im*hitZ(2) + hitZ(4)
          htGa = im*hitgam(1) + hitgam(3)
          htGb = im*hitgam(2) + hitgam(4)

          call amp_zqqQQa_anomZZ(j1,j2,j3,j4,j5,j6,j7,za,zb,htZa,htZb,anomZZ_c1_p1)
          call amp_zqqQQa_anomZZ(j3,j4,j1,j2,j5,j6,j7,za,zb,htZa,htZb,anomZZ_c1_p2)

          call amp_zqqQQa_anomZa(j1,j2,j3,j4,j5,j6,j7,za,zb,htGa,htGb,anomZa_c1_p1)
          call amp_zqqQQa_anomZa(j3,j4,j1,j2,j5,j6,j7,za,zb,htGa,htGb,anomZa_c1_p2)

          call amp_zqqQQa_anomaZ(j1,j2,j3,j4,j5,j6,j7,za,zb,htGa,htGb,anomaZ_c1_p1)
          call amp_zqqQQa_anomaZ(j3,j4,j1,j2,j5,j6,j7,za,zb,htGa,htGb,anomaZ_c1_p2)


          call amp_zqqQQa_anomZZ(j1,j4,j3,j2,j5,j6,j7,za,zb,htZa,htZb,anomZZ_c2_p1)
          call amp_zqqQQa_anomZZ(j3,j2,j1,j4,j5,j6,j7,za,zb,htZa,htZb,anomZZ_c2_p2)

          call amp_zqqQQa_anomZa(j1,j4,j3,j2,j5,j6,j7,za,zb,htGa,htGb,anomZa_c2_p1)
          call amp_zqqQQa_anomZa(j3,j2,j1,j4,j5,j6,j7,za,zb,htGa,htGb,anomZa_c2_p2)

          call amp_zqqQQa_anomaZ(j1,j4,j3,j2,j5,j6,j7,za,zb,htGa,htGb,anomaZ_c2_p1)
          call amp_zqqQQa_anomaZ(j3,j2,j1,j4,j5,j6,j7,za,zb,htGa,htGb,anomaZ_c2_p2)

        endif

        propQQ = s(6,7)/(s(6,7)-zmass**2 + im*zwidth*zmass)
        propQL = t(5,6,7)/(t(5,6,7)-zmass**2 + im*zwidth*zmass)

        propZ_s67 = 1/(s(6,7) - zmass**2 + im*zmass*zwidth)
        propZ_s567 = 1/(t(5,6,7) - zmass**2 + im*zmass*zwidth)
        propA_s67 = 1/s(6,7)
        propA_s567 = 1/t(5,6,7)

        msq = 0._dp
      
        do qi=1,2
          qj=qi

          m70hA = 0._dp
          m70hB = 0._dp
          m70Anom = 0._dp

          do hq=1,2; 
            if (hq == 1) then
              vq = l(qi)
            else
              vq = r(qi)
            endif

              do Qh=1,2
              if (Qh == 1) then
                Qv = l(qj)
              else
                Qv = r(qj)
              endif

              if (hq /= Qh) then
                 ix = -1._dp
              else
                 ix = 1._dp
              endif

              do ha=1,2; do lh=1,2
                if (lh == 1) then
                  vl = l1
                else
                  vl = r1
                endif

                if (hq == lh) then
                   Ql = -abs(q1)
                else
                   Ql = abs(q1)
                endif

                m70hA(1,hq,Qh,ha,lh) =
     &            (q1*Q(qi) + vq*vl*propQQ)*Ai(qi,qj,hq,Qh,ha,lh)
     &           +(q1*Q(qj) + Qv*vl*propQQ)*Bi(qi,qj,Qh,hq,ha,lh)
                m70hA(2,hq,Qh,ha,lh) =
     &            ( (q1*Q(qi) + vq*vl*propQQ)*Ci(qi,qj,hq,Qh,ha,lh)
     &             +(q1*Q(qj) + Qv*vl*propQQ)*Di(qi,qj,Qh,hq,ha,lh))

                m70hB(1,hq,Qh,ha,lh) =
     &            Ql*( (q1*Q(qi)*propA_s567 + vq*vl*propZ_s567)*Aii(hq,Qh,ha,lh)
     &             +ix*(q1*Q(qj)*propA_s567 + Qv*vl*propZ_s567)*Bii(Qh,hq,ha,lh))
                m70hB(2,hq,Qh,ha,lh) =
     &            Ql*( (q1*Q(qi)*propA_s567 + vq*vl*propZ_s567)*Cii(hq,Qh,ha,lh)
     &             +ix*(q1*Q(qj)*propA_s567 + Qv*vl*propZ_s567)*Dii(Qh,hq,ha,lh))

                if (anomtgc) then
                  m70Anom(1,hq,Qh,ha,lh) = m70Anom(1,hq,Qh,ha,lh) +
     &              ( vq*vl*propZ_s567*propZ_s67 * anomZZ_c1_p1(hq,Qh,ha,lh)
     &               +Qv*vl*propZ_s567*propZ_s67 * anomZZ_c1_p2(Qh,hq,ha,lh))
                  m70Anom(2,hq,Qh,ha,lh) = m70Anom(2,hq,Qh,ha,lh) +
     &              ( vq*vl*propZ_s567*propZ_s67 * anomZZ_c2_p1(hq,Qh,ha,lh)
     &               +Qv*vl*propZ_s567*propZ_s67 * anomZZ_c2_p2(Qh,hq,ha,lh))

                  m70Anom(1,hq,Qh,ha,lh) = m70Anom(1,hq,Qh,ha,lh) +
     &              ( q1*vq*propZ_s567*propA_s67 * anomZa_c1_p1(hq,Qh,ha,lh)
     &               +q1*Qv*propZ_s567*propA_s67 * anomZa_c1_p2(Qh,hq,ha,lh))
                  m70Anom(2,hq,Qh,ha,lh) = m70Anom(2,hq,Qh,ha,lh) +
     &              ( q1*vq*propZ_s567*propA_s67 * anomZa_c2_p1(hq,Qh,ha,lh)
     &               +q1*Qv*propZ_s567*propA_s67 * anomZa_c2_p2(Qh,hq,ha,lh))

                  m70Anom(1,hq,Qh,ha,lh) = m70Anom(1,hq,Qh,ha,lh) +
     &              ( Q(qi)*vl*propA_s567*propZ_s67 * anomaZ_c1_p1(hq,Qh,ha,lh)
     &               +Q(qi)*vl*propA_s567*propZ_s67 * anomaZ_c1_p2(Qh,hq,ha,lh))
                  m70Anom(2,hq,Qh,ha,lh) = m70Anom(2,hq,Qh,ha,lh) +
     &              ( Q(qi)*vl*propA_s567*propZ_s67 * anomaZ_c2_p1(hq,Qh,ha,lh)
     &               +Q(qi)*vl*propA_s567*propZ_s67 * anomaZ_c2_p2(Qh,hq,ha,lh))
                endif

              enddo;
            enddo
          enddo; enddo
      
          ampsq = 0._dp

          ! XXX uncomment next two lines
          !m70hA = 0._dp
          !m70hB = 0._dp

          do hq=1,2; do Qh=1,2; do ha=1,2; do lh=1,2;
            ampsq(hq,Qh,ha,lh) = 8._dp * (
     &           abs(m70hA(1,hq,Qh,ha,lh) +
     &                m70hB(1,hq,Qh,ha,lh) + m70Anom(1,hq,Qh,ha,lh))**2
     &         + abs(m70hA(2,hq,Qh,ha,lh) +
     &                m70hB(2,hq,Qh,ha,lh) + m70Anom(2,hq,Qh,ha,lh))**2 )

            if (hq == Qh) then
              ampsq(hq,Qh,ha,lh) = ampsq(hq,Qh,ha,lh) + 
     &          + 8._dp/3._dp*2._dp*real(
     &                (m70hA(1,hq,Qh,ha,lh) +
     &                  m70hB(1,hq,Qh,ha,lh) + m70Anom(1,hq,Qh,ha,lh))*
     &               conjg(m70hA(2,hq,Qh,ha,lh) +
     &                      m70hB(2,hq,Qh,ha,lh) + m70Anom(2,hq,Qh,ha,lh)))
            endif
          enddo; enddo; enddo; enddo

          msq(qi) = sum(ampsq(:,:,:,:))
        enddo

        !include color factor, gauge couplings, but no averaging
        msq = msq * 8*esq**3*gsq**2

      end

      subroutine amp_zqqQQa_ql_new(j1,j2,j3,j4,j5,j6,j7,za,zb,ai)
        use zajj_treeamps_m
        implicit none
        
        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(out) :: ai(2,2,2,2)

        ! hq,hQ,ha,hl

        ai(2,2,2,2) = zajj_tree_qqQQ_fsr_plus(j2,j1,j6,j7,j5,j3,j4,za,zb)
        ai(2,2,1,2) = zajj_tree_qqQQ_fsr_minus(j2,j1,j6,j7,j5,j3,j4,za,zb)

        ! lepton exchange
        ai(2,2,2,1) = zajj_tree_qqQQ_fsr_plus(j2,j1,j7,j6,j5,j3,j4,za,zb)
        ai(2,2,1,1) = zajj_tree_qqQQ_fsr_minus(j2,j1,j7,j6,j5,j3,j4,za,zb)

        ! Q exchange
        ai(2,1,2,2) = zajj_tree_qqQQ_fsr_plus(j2,j1,j6,j7,j5,j4,j3,za,zb)
        ai(2,1,1,2) = zajj_tree_qqQQ_fsr_minus(j2,j1,j6,j7,j5,j4,j3,za,zb)

        ! lepton +  Q exchange
        ai(2,1,2,1) = zajj_tree_qqQQ_fsr_plus(j2,j1,j7,j6,j5,j4,j3,za,zb)
        ai(2,1,1,1) = zajj_tree_qqQQ_fsr_minus(j2,j1,j7,j6,j5,j4,j3,za,zb)


        ai(1,1,1,1) = conjg(ai(2,2,2,2))
        ai(1,1,2,1) = conjg(ai(2,2,1,2))
        
        ai(1,1,1,2) = conjg(ai(2,2,2,1))
        ai(1,1,2,2) = conjg(ai(2,2,1,1))

        ai(1,2,1,1) = conjg(ai(2,1,2,2))
        ai(1,2,2,1) = conjg(ai(2,1,1,2))

        ai(1,2,1,2) = conjg(ai(2,1,2,1))
        ai(1,2,2,2) = conjg(ai(2,1,1,1))

c       ai(1,1,1,1) = conjg(zajj_tree_qqQQ_fsr_plus(j2,j1,j6,j7,j5,j3,j4,za,zb))
c       ai(1,1,2,1) = conjg(zajj_tree_qqQQ_fsr_minus(j2,j1,j6,j7,j5,j3,j4,za,zb))

c       ai(1,1,1,2) = conjg(zajj_tree_qqQQ_fsr_plus(j2,j1,j7,j6,j5,j3,j4,za,zb))
c       ai(1,1,2,2) = conjg(zajj_tree_qqQQ_fsr_minus(j2,j1,j7,j6,j5,j3,j4,za,zb))

c       ai(1,2,1,1) = conjg(zajj_tree_qqQQ_fsr_plus(j2,j1,j6,j7,j5,j4,j3,za,zb))
c       ai(1,2,2,1) = conjg(zajj_tree_qqQQ_fsr_minus(j2,j1,j6,j7,j5,j4,j3,za,zb))

c       ai(1,2,1,2) = conjg(zajj_tree_qqQQ_fsr_plus(j2,j1,j7,j6,j5,j4,j3,za,zb))
c       ai(1,2,2,2) = conjg(zajj_tree_qqQQ_fsr_minus(j2,j1,j7,j6,j5,j4,j3,za,zb))

      end
      
************************************************************************
*     0 ---> q(p1)+qb(p2)+Q(p3)+Qb(p4)+ph(p5)+lb(p6)+l(p7)              *
************************************************************************
      
      subroutine amp_zqqQQa_qq_new(i1,i2,i3,i4,i5,i6,i7,za,zb,ai,bi)
        implicit none
        include 'types.f'
        include 'constants.f'
        include 'mxpart.f'
        include 'nf.f'
        include 'ewcharge.f'

        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        complex(dp), intent(out) :: ai(2,2,2,2,2,2),bi(2,2,2,2,2,2)

        complex(dp)::a1(2,2,2,2),a2(2,2,2,2),a3(2,2,2,2),a4(2,2,2,2),
     &               b1(2,2,2,2),b2(2,2,2,2),b3(2,2,2,2),b4(2,2,2,2)
        integer:: q12,q34

        call nagyqqQQg_clean(i1,i2,i3,i4,i5,i6,i7,za,zb,a1,a2,a3,a4)
        call nagyqqQQg_clean(i3,i4,i1,i2,i5,i6,i7,za,zb,b1,b2,b3,b4)

        do q12=1,2; do q34=1,2
          ai(q12,q34,:,:,:,:)= Q(q12)*(a1+a2) + Q(q34)*(a3+a4)
          bi(q12,q34,:,:,:,:)= Q(q34)*(b1+b2) + Q(q12)*(b3+b4)
        enddo; enddo
      end

      ! hep-ph/9806317
      pure function nagyqqQQg_a1_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A50
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a1_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a1_ppp = 
     & -zb(i1,i5)*t2a(i4,i1,i5,i3)*t2a(i4,i2,i6,i7)*za(i6,i2)
     &  /(za(i4,i5)*s(i1,i5)*s(i3,i4)*t(i2,i6,i7))  
     & -zb(i1,i7)*t2a(i6,i1,i7,i5)*za(i4,i2)**2*zb(i2,i3)
     &  /(za(i4,i5)*s(i3,i4)*t(i2,i3,i4)*t(i1,i6,i7))  
     & -t2a(i4,i1,i5,i7)*t2a(i6,i2,i4,i3)*za(i4,i2)
     &  /(za(i1,i5)*za(i5,i4)*s(i3,i4)*t(i2,i3,i4))  
     & +zb(i5,i3)*t2a(i4,i3,i5,i1)*t2a(i4,i2,i6,i7)*za(i6,i2)
     &  /(za(i4,i5)*s(i3,i4)*t(i3,i4,i5)*t(i2,i6,i7))  
     & +zb(i1,i7)
     & *(t2a(i6,i1,i7,i3)*za(i3,i4)+t2a(i6,i1,i7,i5)*za(i5,i4))
     & *zb(i3,i5)*za(i4,i2)/(za(i4,i5)*s(i3,i4)*t(i3,i4,i5)*t(i1,i6,i7))  

      end function

      pure function nagyqqQQg_a1_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A51
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a1_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a1_ppm = 
     & +zb(i1,i3)**2*za(i5,i1)*t2a(i4,i2,i6,i7)*za(i6,i2)
     &  /(zb(i3,i5)*s(i1,i5)*s(i3,i4)*t(i2,i6,i7))  
     & +zb(i1,i7)*t2a(i6,i1,i7,i3)*t2a(i5,i2,i4,i3)*za(i4,i2)
     &  /(zb(i3,i5)*s(i3,i4)*t(i2,i3,i4)*t(i1,i6,i7))  
     & -zb(i1,i3)*zb(i1,i7)*t2a(i6,i2,i4,i3)*za(i4,i2)
     &  /(zb(i1,i5)*zb(i5,i3)*s(i3,i4)*t(i2,i3,i4))  
     & +zb(i1,i3)*za(i5,i4)
     &  *(+zb(i3,i4)*za(i4,i2)*zb(i2,i7)+zb(i3,i5)*za(i5,i2)*zb(i2,i7)
     &    +zb(i3,i4)*za(i4,i6)*zb(i6,i7)+zb(i3,i5)*za(i5,i6)*zb(i6,i7))
     & *za(i6,i2)/(zb(i3,i5)*s(i3,i4)*t(i3,i4,i5)*t(i2,i6,i7))  
     & -zb(i1,i7)*t2a(i6,i1,i7,i3)*za(i5,i4)*t2a(i2,i4,i5,i3)
     &  /(zb(i3,i5)*s(i3,i4)*t(i3,i4,i5)*t(i1,i6,i7))  

      end function

      pure function nagyqqQQg_a2_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A53
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a2_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a2_ppp = 
     & -zb(i5,i3)*t2a(i4,i3,i5,i1)*t2a(i4,i2,i6,i7)*za(i6,i2)
     &  /(za(i4,i5)*s(i3,i4)*t(i3,i4,i5)*t(i2,i6,i7))  
     & -zb(i1,i7)*t2a(i6,i1,i7,i3)*zb(i2,i5)*za(i4,i2)**2
     &  /(za(i4,i5)*s(i2,i5)*s(i3,i4)*t(i1,i6,i7))  
     & -zb(i1,i3)*t2a(i4,i1,i3,i5)*t2a(i4,i2,i6,i7)*za(i6,i2)
     &  /(za(i4,i5)*s(i3,i4)*t(i1,i3,i4)*t(i2,i6,i7))  
     & -zb(i1,i3)*t2a(i4,i1,i3,i7)*za(i6,i2)*za(i4,i2)
     &  /(za(i4,i5)*za(i5,i2)*s(i3,i4)*t(i1,i3,i4))  
     & -zb(i1,i7)
     &  *(t2a(i6,i1,i7,i3)*za(i3,i4)+t2a(i6,i1,i7,i5)*za(i5,i4))
     & *zb(i3,i5)*za(i4,i2)/(za(i4,i5)*s(i3,i4)*t(i3,i4,i5)*t(i1,i6,i7))  

      end function

      pure function nagyqqQQg_a2_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A54
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a2_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a2_ppm = 
     & -zb(i1,i3)*za(i5,i4)*
     & (+zb(i3,i4)*za(i4,i2)*zb(i2,i7)+zb(i3,i5)*za(i5,i2)*zb(i2,i7)
     &  +zb(i3,i4)*za(i4,i6)*zb(i6,i7)+zb(i3,i5)*za(i5,i6)*zb(i6,i7))
     & *za(i6,i2)/(zb(i3,i5)*s(i3,i4)*t(i3,i4,i5)*t(i2,i6,i7))  
     & +zb(i1,i3)**2*za(i4,i1)*t2a(i5,i2,i6,i7)*za(i6,i2)
     &  /(zb(i3,i5)*s(i3,i4)*t(i1,i3,i4)*t(i2,i6,i7))  
     & +zb(i1,i7)*t2a(i6,i1,i7,i3)*za(i5,i4)*t2a(i2,i4,i5,i3)
     &  /(zb(i3,i5)*s(i3,i4)*t(i3,i4,i5)*t(i1,i6,i7))  
     & -zb(i1,i3)*t2a(i4,i1,i3,i7)*t2a(i6,i2,i5,i3)
     &  /(zb(i3,i5)*zb(i5,i2)*s(i3,i4)*t(i1,i3,i4))  
     & +zb(i1,i7)*t2a(i6,i1,i7,i3)*za(i5,i2)*t2a(i4,i2,i5,i3)
     &  /(zb(i3,i5)*s(i2,i5)*s(i3,i4)*t(i1,i6,i7))  

      end function

      pure function nagyqqQQg_a3_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A56
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a3_ppp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a3_ppp = 
     & -zb(i5,i3)*t2a(i4,i3,i5,i1)*t2a(i4,i2,i6,i7)*za(i6,i2)
     & /(za(i4,i5)*s(i3,i5)*t(i3,i4,i5)*t(i2,i6,i7))
     & -zb(i1,i7)
     & *(t2a(i6,i1,i7,i3)*za(i3,i4)+t2a(i6,i1,i7,i5)*za(i5,i4))
     & *zb(i3,i5)*za(i4,i2)/(za(i4,i5)*s(i3,i5)*t(i3,i4,i5)*t(i1,i6,i7))

      end function

      pure function nagyqqQQg_a3_pmm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A57
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a3_pmm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a3_pmm = 
     &  -zb(i1,i4)*za(i5,i3)*
     & (+zb(i4,i3)*za(i3,i2)*zb(i2,i7)+zb(i4,i5)*za(i5,i2)*zb(i2,i7)
     &  +zb(i4,i3)*za(i3,i6)*zb(i6,i7)+zb(i4,i5)*za(i5,i6)*zb(i6,i7))
     & *za(i6,i2)/(zb(i4,i5)*s(i3,i5)*t(i4,i3,i5)*t(i2,i6,i7))
     & +zb(i1,i7)*t2a(i6,i1,i7,i4)*za(i5,i3)*t2a(i2,i3,i5,i4)
     & /(zb(i4,i5)*s(i3,i5)*t(i4,i3,i5)*t(i1,i6,i7))

      end function

      pure function nagyqqQQg_a4_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A58
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a4_ppm
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a4_ppm = 
     & +zb(i1,i3)*za(i5,i4)*
     & (+zb(i3,i4)*za(i4,i2)*zb(i2,i7)+zb(i3,i4)*za(i4,i6)*zb(i6,i7)
     &  +zb(i3,i5)*za(i5,i2)*zb(i2,i7)+zb(i3,i5)*za(i5,i6)*zb(i6,i7))
     & *za(i6,i2)/(zb(i3,i5)*s(i4,i5)*t(i3,i4,i5)*t(i2,i6,i7))
     & -zb(i1,i7)*t2a(i6,i1,i7,i3)*za(i5,i4)*t2a(i2,i4,i5,i3)
     &  /(zb(i3,i5)*s(i4,i5)*t(i3,i4,i5)*t(i1,i6,i7))

      end function

      pure function nagyqqQQg_a4_pmp(i1,i2,i3,i4,i5,i6,i7,za,zb) ! A59
        implicit none
        include 'types.f'
        include 'mxpart.f'

        complex(dp) :: nagyqqQQg_a4_pmp
        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)

        real(dp) :: s,t
        complex(dp) :: t2a

        t2a(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))
        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

         nagyqqQQg_a4_pmp = 
     &  +zb(i5,i4)*t2a(i3,i4,i5,i1)*t2a(i3,i2,i6,i7)
     & *za(i6,i2)/(za(i3,i5)*s(i4,i5)*t(i4,i3,i5)*t(i2,i6,i7))
     &  +zb(i1,i7)
     & *(t2a(i6,i1,i7,i4)*za(i4,i3)+t2a(i6,i1,i7,i5)*za(i5,i3))
     & *zb(i4,i5)*za(i3,i2)/(za(i3,i5)*s(i4,i5)*t(i4,i3,i5)*t(i1,i6,i7))

      end function

      subroutine nagyqqQQg_clean(i1,i2,i3,i4,i5,i6,i7,za,zb,a1,a2,a3,a4)
        implicit none
        include 'types.f'
        include 'mxpart.f'

        integer, intent(in) :: i1,i2,i3,i4,i5,i6,i7
        complex(dp), intent(in) :: za(mxpart,mxpart),zb(mxpart,mxpart)
        complex(dp), intent(out) :: a1(2,2,2,2),a2(2,2,2,2),
     &                              a3(2,2,2,2),a4(2,2,2,2)
        real(dp) :: s
        s(i1,i2) = real(za(i1,i2)*zb(i2,i1))

        complex(dp) :: nagyqqQQg_a1_ppp, nagyqqQQg_a1_ppm
        complex(dp) :: nagyqqQQg_a2_ppp, nagyqqQQg_a2_ppm
        complex(dp) :: nagyqqQQg_a3_ppp, nagyqqQQg_a3_pmm
        complex(dp) :: nagyqqQQg_a4_ppm, nagyqqQQg_a4_pmp

        a1 = 0._dp
        a2 = 0._dp
        a3 = 0._dp
        a4 = 0._dp

        ! hq,Qh,hg,lh

c a1
        a1(2,2,2,2) = nagyqqQQg_a1_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)
        a1(2,2,1,2) = nagyqqQQg_a1_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)
        ! via 3 <-> 4
        a1(2,1,2,2) = nagyqqQQg_a1_ppp(i1,i2,i4,i3,i5,i6,i7,za,zb)/s(i6,i7)
        a1(2,1,1,2) = nagyqqQQg_a1_ppm(i1,i2,i4,i3,i5,i6,i7,za,zb)/s(i6,i7)

        ! via 6 <-> 7
        a1(2,2,2,1) = nagyqqQQg_a1_ppp(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)
        a1(2,2,1,1) = nagyqqQQg_a1_ppm(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)
        a1(2,1,2,1) = nagyqqQQg_a1_ppp(i1,i2,i4,i3,i5,i7,i6,za,zb)/s(i6,i7)
        a1(2,1,1,1) = nagyqqQQg_a1_ppm(i1,i2,i4,i3,i5,i7,i6,za,zb)/s(i6,i7)


        a1(1,1,1,1) = conjg(a1(2,2,2,2))
        a1(1,1,2,1) = conjg(a1(2,2,1,2))
        a1(1,2,1,1) = conjg(a1(2,1,2,2))
        a1(1,2,2,1) = conjg(a1(2,1,1,2))
        a1(1,1,1,2) = conjg(a1(2,2,2,1))
        a1(1,1,2,2) = conjg(a1(2,2,1,1))
        a1(1,2,1,2) = conjg(a1(2,1,2,1))
        a1(1,2,2,2) = conjg(a1(2,1,1,1))

c a2
        a2(2,2,2,2) = nagyqqQQg_a2_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)
        a2(2,2,1,2) = nagyqqQQg_a2_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)
        ! via 3 <-> 4
        a2(2,1,2,2) = nagyqqQQg_a2_ppp(i1,i2,i4,i3,i5,i6,i7,za,zb)/s(i6,i7)
        a2(2,1,1,2) = nagyqqQQg_a2_ppm(i1,i2,i4,i3,i5,i6,i7,za,zb)/s(i6,i7)

        ! via 6 <-> 7
        a2(2,2,2,1) = nagyqqQQg_a2_ppp(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)
        a2(2,2,1,1) = nagyqqQQg_a2_ppm(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)
        a2(2,1,2,1) = nagyqqQQg_a2_ppp(i1,i2,i4,i3,i5,i7,i6,za,zb)/s(i6,i7)
        a2(2,1,1,1) = nagyqqQQg_a2_ppm(i1,i2,i4,i3,i5,i7,i6,za,zb)/s(i6,i7)

        a2(1,1,1,1) = conjg(a2(2,2,2,2))
        a2(1,1,2,1) = conjg(a2(2,2,1,2))
        a2(1,2,1,1) = conjg(a2(2,1,2,2))
        a2(1,2,2,1) = conjg(a2(2,1,1,2))
        a2(1,1,1,2) = conjg(a2(2,2,2,1))
        a2(1,1,2,2) = conjg(a2(2,2,1,1))
        a2(1,2,1,2) = conjg(a2(2,1,2,1))
        a2(1,2,2,2) = conjg(a2(2,1,1,1))

c a3
        a3(2,2,2,2) = nagyqqQQg_a3_ppp(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)
        a3(2,1,1,2) = nagyqqQQg_a3_pmm(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)

        ! via 6 <> 7
        a3(2,2,2,1) = nagyqqQQg_a3_ppp(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)
        a3(2,1,1,1) = nagyqqQQg_a3_pmm(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)

        a3(1,1,1,1) = conjg(a3(2,2,2,2))
        a3(1,2,2,1) = conjg(a3(2,1,1,2))
        a3(1,1,1,2) = conjg(a3(2,2,2,1))
        a3(1,2,2,2) = conjg(a3(2,1,1,1))

c a4
        a4(2,2,1,2) = nagyqqQQg_a4_ppm(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)
        a4(2,1,2,2) = nagyqqQQg_a4_pmp(i1,i2,i3,i4,i5,i6,i7,za,zb)/s(i6,i7)

        ! via 6 <-> 7
        a4(2,2,1,1) = nagyqqQQg_a4_ppm(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)
        a4(2,1,2,1) = nagyqqQQg_a4_pmp(i1,i2,i3,i4,i5,i7,i6,za,zb)/s(i6,i7)

        a4(1,1,2,1) = conjg(a4(2,2,1,2))
        a4(1,2,1,1) = conjg(a4(2,1,2,2))
        a4(1,1,2,2) = conjg(a4(2,2,1,1))
        a4(1,2,1,2) = conjg(a4(2,1,2,1))

      end

