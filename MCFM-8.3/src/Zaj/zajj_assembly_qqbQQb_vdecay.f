      subroutine msq_zqqbQQbgam_vdecay(p,qqb_ij,qbq_ij,qqb_ii,qbq_ii)
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
        real(dp), intent(out) :: qqb_ij(2, 2,2), qbq_ij(2, 2,2)
        real(dp), intent(out) :: qqb_ii(2, 2), qbq_ii(2, 2)

        integer :: flavor
        real(dp) :: pnt(mxpart,4)

        ! convert to Nagy-Trocsanyi momentum convention
        ! and do crossing manually in calling amplitudes
        pnt(1,:) = p(2,:)
        pnt(2,:) = p(1,:)
        pnt(3,:) = p(6,:)
        pnt(4,:) = p(7,:)
        pnt(5,:) = p(5,:)
        pnt(6,:) = p(4,:)
        pnt(7,:) = p(3,:)

        ! NT convention, but crossed for vdecay:
c       pnt(1,:) = p(4,:)
c       pnt(2,:) = p(3,:)
c       pnt(3,:) = p(6,:)
c       pnt(4,:) = p(7,:)
c       pnt(5,:) = p(5,:)
c       pnt(6,:) = p(2,:)
c       pnt(7,:) = p(1,:)

        qqb_ij = 0._dp
        qbq_ij = 0._dp
        qqb_ii = 0._dp
        qbq_ii = 0._dp

        call spinoru(7,pnt,za,zb)

c       call zqqbQQba_msqij(1,2,3,4,5,6,7,za,zb,qqb_ij)
c       call zqqbQQba_msqij(2,1,3,4,5,6,7,za,zb,qbq_ij)

c       call zqqbQQba_msqii(1,2,3,4,5,6,7,za,zb,qqb_ii)
c       call zqqbQQba_msqii(2,1,3,4,5,6,7,za,zb,qbq_ii)

        do flavor=1,2
            call zqqbQQba_msqij_vdecay(flavor, 6,7,3,4,5,1,2,za,zb,qqb_ij(flavor,:,:)) 
            call zqqbQQba_msqij_vdecay(flavor, 6,7,3,4,5,2,1,za,zb,qbq_ij(flavor,:,:)) 

            call zqqbQQba_msqii_vdecay(flavor, 6,7,3,4,5,1,2,za,zb,qqb_ii(flavor,:)) 
            call zqqbQQba_msqii_vdecay(flavor, 6,7,3,4,5,2,1,za,zb,qbq_ii(flavor,:)) 
        enddo

        qqb_ij = qqb_ij * aveqq
        qbq_ij = qbq_ij * aveqq

        qqb_ii = qqb_ii * aveqq*half
        qbq_ii = qbq_ii * aveqq*half

      end

      subroutine zqqbQQba_msqij_vdecay(flavor,j1,j2,j3,j4,j5,j6,j7,za,zb,msq)
        implicit none
********************************************************************
* return squared matrix element for                                *
* 0 -> q(p1)+qb(p2)+Q(p3)+Qb(p4)+gamma(p5)+lb(p6)+l(p7)            *
* q and Q have different flavour                                   *
* above is Nagy-Trocsanyi momentum assignment                      *
* NO average, NO identical particle factor included                *                                        
********************************************************************
      
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'nf.f'
        include 'cplx.h'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'zcouple.f'
        include 'masses.f'
        include 'ewcharge.f'

        integer, intent(in) :: flavor
        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        real(dp), intent(out) :: msq(2,2)

        integer :: qi,qj
        real(dp) :: ix,Ql
        integer :: hq,Qh,lh,ha
        real(dp) :: vq, Qv, vl
        complex(dp) :: propQQ,propQL
        complex(dp) :: m70hA(2,2, 2,2,2,2),m70hB(2,2, 2,2,2,2)
        complex(dp) :: Ai(2,2,2,2,2,2),Bi(2,2,2,2,2,2)
        complex(dp) :: Aii(2,2,2,2),Bii(2,2,2,2)

        complex(dp) :: propZ_s67, propZ_s567
        complex(dp) :: propA_s67, propA_s567

        real(dp) :: s,t
        s(j1,j2) = real(za(j1,j2)*zb(j2,j1))
        t(j1,j2,j3) = s(j1,j2) + s(j2,j3) + s(j3,j1)

        Ai = 0._dp
        Bi = 0._dp
        Aii = 0._dp
        Bii = 0._dp

        call amp_zqqQQa_qq_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Ai,Bi)
        call amp_zqqQQa_ql_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Aii)
        call amp_zqqQQa_ql_new(j3,j4,j1,j2,j5,j6,j7,za,zb,Bii)

        propQQ = s(j6,j7)/(s(j6,j7)-zmass**2 + im*zwidth*zmass)
        propQL = t(j5,j6,j7)/(t(j5,j6,j7)-zmass**2 + im*zwidth*zmass)

        propZ_s67 = 1/(s(j6,j7) - zmass**2 + im*zmass*zwidth)
        propZ_s567 = 1/(t(j5,j6,j7) - zmass**2 + im*zmass*zwidth)
        propA_s67 = 1/s(j6,j7)
        propA_s567 = 1/t(j5,j6,j7)

        msq = 0._dp

        do qi=1,2; do qj=1,2
          m70hA = 0._dp
          m70hB = 0._dp

          do hq=1,2; 
            do Qh=1,2; 
              do ha=1,2; do lh=1,2
                if (hq == 1) then
                  vq = l(qi)
                else
                  vq = r(qi)
                endif

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

                if (lh == 1) then
                  vl = l(flavor)
                else
                  vl = r(flavor)
                endif

                if (hq == lh) then
                   Ql = +Q(flavor)
                else
                   Ql = -Q(flavor)
                endif

                m70hA(qi,qj, hq,Qh,ha,lh) = 
     &               (Q(flavor)*Q(qi) + vq*vl*propQQ)*Ai(qi,qj,hq,Qh,ha,lh)
     &             + (Q(flavor)*Q(qj) + Qv*vl*propQQ)*Bi(qi,qj,Qh,hq,ha,lh)
                m70hB(qi,qj, hq,Qh,ha,lh) = 
     &             + Ql*( (Q(flavor)*Q(qi)*propA_s567 + vq*vl*propZ_s567)*Aii(hq,Qh,ha,lh)
     &               + ix*(Q(flavor)*Q(qj)*propA_s567 + Qv*vl*propZ_s567)*Bii(Qh,hq,ha,lh))

              enddo
            enddo
          enddo; enddo

          ! XXX uncomment next two lines
          !m70hA = 0._dp
          !m70hB = 0._dp

          msq(qi,qj) = sum(abs(m70hA(qi,qj,:,:,:,:) +
     &            m70hB(qi,qj,:,:,:,:))**2)
          ! include color factor, gauge couplings, but no averaging
          msq(qi,qj)=8._dp*8._dp*esq**3*gsq**2*msq(qi,qj)*cA

        enddo; enddo

      end


      subroutine zqqbQQba_msqii_vdecay(flavor,j1,j2,j3,j4,j5,j6,j7,za,zb,msq)
        implicit none
********************************************************************
* return squared matrix element for                                *
* 0 -> q(p1)+qb(p2)+q(p3)+qb(p4)+gamma(p5)+lb(p6)+l(p7)            *
* q and Q have the same flavour                                    *           
* above is Nagy-Trocsanyi momentum assignment                      *
* NO average, NO identical particle factor included                *                                        
********************************************************************
        include 'types.f'
        include 'mxpart.f'
        include 'constants.f'
        include 'nf.f'
        include 'qcdcouple.f'
        include 'ewcouple.f'
        include 'zcouple.f'
        include 'masses.f'
        include 'ewcharge.f'

        integer, intent(in) :: flavor
        integer, intent(in) :: j1,j2,j3,j4,j5,j6,j7
        complex(dp), intent(in) :: za(mxpart,mxpart), zb(mxpart,mxpart)
        real(dp), intent(out) :: msq(2)

        integer :: qi,qj
        real(dp) :: ix,Ql
        real(dp) :: ampsq(2,2,2,2)
        integer :: hq,Qh,lh,ha
        complex(dp) :: m70hA(2,2,2,2,2),m70hB(2,2,2,2,2)
        complex(dp) :: Ai(2,2,2,2,2,2),Bi(2,2,2,2,2,2)
        complex(dp) :: Ci(2,2,2,2,2,2),Di(2,2,2,2,2,2)
        complex(dp) :: Aii(2,2,2,2),Bii(2,2,2,2)
        complex(dp) :: Cii(2,2,2,2),Dii(2,2,2,2)

        real(dp) :: vq, Qv, vl

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

        call amp_zqqQQa_qq_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Ai,Bi)
        call amp_zqqQQa_qq_new(j1,j4,j3,j2,j5,j6,j7,za,zb,Ci,Di)

        call amp_zqqQQa_ql_new(j1,j2,j3,j4,j5,j6,j7,za,zb,Aii)
        call amp_zqqQQa_ql_new(j3,j4,j1,j2,j5,j6,j7,za,zb,Bii)

        call amp_zqqQQa_ql_new(j1,j4,j3,j2,j5,j6,j7,za,zb,Cii)
        call amp_zqqQQa_ql_new(j3,j2,j1,j4,j5,j6,j7,za,zb,Dii)

        propQQ = s(j6,j7)/(s(j6,j7)-zmass**2 + im*zwidth*zmass)
        propQL = t(j5,j6,j7)/(t(j5,j6,j7)-zmass**2 + im*zwidth*zmass)

        propZ_s67 = 1/(s(j6,j7) - zmass**2 + im*zmass*zwidth)
        propZ_s567 = 1/(t(j5,j6,j7) - zmass**2 + im*zmass*zwidth)
        propA_s67 = 1/s(j6,j7)
        propA_s567 = 1/t(j5,j6,j7)

        msq = 0._dp
      
        do qi=1,2
          qj=qi

          m70hA = 0._dp
          m70hB = 0._dp

          do hq=1,2; 
            do Qh=1,2
              do ha=1,2; do lh=1,2
                if (hq == 1) then
                  vq = l(qi)
                else
                  vq = r(qi)
                endif

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

                if (lh == 1) then
                  vl = l(flavor)
                else
                  vl = r(flavor)
                endif

                if (hq == lh) then
                   Ql = +Q(flavor)
                else
                   Ql = -Q(flavor)
                endif

                ! A: original q_q type
                ! m70hA contributions 1,2: crossing qbar and Qbar
                m70hA(1,hq,Qh,ha,lh) =
                ! Z/gamma coupling to different quark line (q,Q)
     &            (Q(flavor)*Q(qi) + vq*vl*propQQ)*Ai(qi,qj,hq,Qh,ha,lh)
     &           +(Q(flavor)*Q(qj) + Qv*vl*propQQ)*Bi(qi,qj,Qh,hq,ha,lh)
                m70hA(2,hq,Qh,ha,lh) =
     &            ( (Q(flavor)*Q(qi) + vq*vl*propQQ)*Ci(qi,qj,hq,Qh,ha,lh)
     &             +(Q(flavor)*Q(qj) + Qv*vl*propQQ)*Di(qi,qj,Qh,hq,ha,lh))

                ! B: original q_l type
                m70hB(1,hq,Qh,ha,lh) =
     &            Ql*( (Q(flavor)*Q(qi)*propA_s567 + vq*vl*propZ_s567)*Aii(hq,Qh,ha,lh)
     &             +ix*(Q(flavor)*Q(qj)*propA_s567 + Qv*vl*propZ_s567)*Bii(Qh,hq,ha,lh))
                m70hB(2,hq,Qh,ha,lh) =
     &            Ql*( (Q(flavor)*Q(qi)*propA_s567 + vq*vl*propZ_s567)*Cii(hq,Qh,ha,lh)
     &             +ix*(Q(flavor)*Q(qj)*propA_s567 + Qv*vl*propZ_s567)*Dii(Qh,hq,ha,lh))

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
     &                m70hB(1,hq,Qh,ha,lh))**2
     &         + abs(m70hA(2,hq,Qh,ha,lh) +
     &                m70hB(2,hq,Qh,ha,lh))**2 )

            if (hq == Qh) then
              ampsq(hq,Qh,ha,lh) = ampsq(hq,Qh,ha,lh) + 
     &          + 8._dp/3._dp*2._dp*real(
     &                (m70hA(1,hq,Qh,ha,lh) +
     &                  m70hB(1,hq,Qh,ha,lh))*
     &               conjg(m70hA(2,hq,Qh,ha,lh) +
     &                      m70hB(2,hq,Qh,ha,lh)))
            endif
          enddo; enddo; enddo; enddo

          msq(qi) = sum(ampsq(:,:,:,:))
        enddo

        !include color factor, gauge couplings, but no averaging
        msq = msq * 8*esq**3*gsq**2*cA

      end

