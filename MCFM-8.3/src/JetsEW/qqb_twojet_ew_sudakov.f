      subroutine qqb_twojet_ew_sudakov(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'sprods_com.f'
      include 'scheme.f'
      real(dp):: msq(fn:nf,fn:nf),p(mxpart,4),ss,tt,uu,
     &     msq_qa_iijj(0:nf,1:nf),msq_qa_ijij(0:nf,1:nf),
     &     msq_qq_ijij(0:nf,1:nf),msq_aq_iijj(0:nf,1:nf),
     &     msq_aq_ijij(0:nf,1:nf),msq_aa_ijij(0:nf,1:nf),
     &     msq_gg_qa(nf),msq_gq_gq(nf),msq_qg_qg(nf),
     &     msq_qa_gg(nf),msq_aq_gg(nf),msq_ga_ga(nf),
     &     msq_ag_ag(nf)
      integer j,k

!      scheme='dred' ! dummy value of scheme, in case routine called from virtint
 
      call dotem(4,p,s)
      ss = s(1,2)
      tt = s(1,3)
      uu = s(2,3)

      call dijet_qqb_ew_sud(msq_qa_iijj,ss,tt,uu)
      call dijet_qqb_ew_sud(msq_qa_ijij,tt,ss,uu)
      call dijet_qqb_ew_sud(msq_qq_ijij,tt,uu,ss)
      call dijet_qqb_ew_sud(msq_aq_iijj,ss,tt,uu)
      call dijet_qqb_ew_sud(msq_aq_ijij,tt,ss,uu)
      call dijet_qqb_ew_sud(msq_aa_ijij,tt,uu,ss)

      msq_gg_qa(1:nf)=+avegg*msq_qa_iijj(0,1:nf)
      msq_gq_gq(1:nf)=-aveqg*msq_qq_ijij(0,1:nf)
      msq_qg_qg(1:nf)=-aveqg*msq_qa_ijij(0,1:nf)
      msq_qa_gg(1:nf)=+aveqq*msq_qa_iijj(0,1:nf)
      msq_aq_gg(1:nf)=+aveqq*msq_qa_iijj(0,1:nf)
      msq_ga_ga(1:nf)=-aveqg*msq_qa_ijij(0,1:nf)
      msq_ag_ag(1:nf)=-aveqg*msq_qq_ijij(0,1:nf)

      msq=0._dp
      do j=fn+1,nf-1
         do k=fn+1,nf-1
C --- qa
            if((j > 0) .and. (k < 0)) then
               if(j == -k) then
                  msq(j,k) = msq_qa_iijj(j,1) + msq_qa_iijj(j,2) 
     &                 + msq_qa_iijj(j,3) + msq_qa_iijj(j,4)
     &                 + msq_qa_gg(j)*half
                else
                  msq(j,k) = msq_qa_ijij(j,-k)
               end if

C --- qq
            else if((j > 0) .and. (k > 0)) then
               if(j == k) then
                  msq(j,k) = msq_qq_ijij(j,k)*half
               else
                  msq(j,k) = msq_qq_ijij(j,k)
               end if

C --- aq
            else if((j < 0) .and. (k > 0)) then
               if(k == -j) then
                  msq(j,k) = msq_aq_iijj(k,1) + msq_aq_iijj(k,2) 
     &                 + msq_aq_iijj(k,3) + msq_aq_iijj(k,4)
     &                 + msq_aq_gg(k)*half
               else
                  msq(j,k) = msq_aq_ijij(-j,k)
               end if

C --- aa
            else if((j < 0) .and. (k < 0)) then
               if(j == k) then
                  msq(j,k) = msq_aa_ijij(-j,-k)*half
               else
                  msq(j,k) = msq_aa_ijij(-j,-k)
               end if

C --- gg
            else if((j == 0) .and. (k == 0)) then
               msq(j,k) = msq_gg_qa(1) + msq_gg_qa(2) 
     &              + msq_gg_qa(3) + msq_gg_qa(4)

C --- gq
            else if((j == 0) .and. (k > 0)) then
               msq(j,k) = msq_gq_gq(k)

C --- qg
            else if((j > 0) .and. (k == 0)) then
               msq(j,k) = msq_qg_qg(j)

C --- ga
            else if((j == 0) .and. (k < 0)) then
               msq(j,k) = msq_ga_ga(-k)

C --- ag
            else if((j < 0) .and. (k == 0)) then
               msq(j,k) = msq_ag_ag(-j)
            end if

         end do
      end do
      
      return
      end
      
      
      subroutine dijet_qqb_ew_sud(msq,ss,tt,uu) 
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'breit.f'      
      include 'cplx.h'      
      integer j,k,cs
      real(dp):: msq(0:nf,1:nf),ss,tt,uu
      real(dp):: wtgg_up,wtgg_dn,t1,t2,ro,cdiff,csame,
     & wta_up,wta_dn,wtb_up,wtb_dn,wti_up,wti_dn,
     & DLs,SLs,SLt,SLt1,SLt2,Lm,sw2,cw2,alpha_em,
     & Cew_uL,Cew_ubL,Cew_uR,Cew_ubR,
     & Cew_dL,Cew_dbL,Cew_dR,Cew_dbR,
     & Iz2_uL,Iz2_ubL,Iz2_uR,Iz2_ubR,
     & Iz2_dL,Iz2_dbL,Iz2_dR,Iz2_dbR,
     & Iz_uL,Iz_ubL,Iz_uR,Iz_ubR,
     & Iz_dL,Iz_dbL,Iz_dR,Iz_dbR,
     & csame_up,csame_dn,cdiff_up,cdiff_dn,
     & csame_up_up_ssc,csame_up_dn_ssc,csame_dn_up_ssc,csame_dn_dn_ssc,
     & cdiff_up_up_ssc,cdiff_up_dn_ssc,cdiff_dn_up_ssc,cdiff_dn_dn_ssc,
     & wtqqb_up_up,wtqqb_up_dn,wtqqb_dn_up,wtqqb_dn_dn,
     & wtqqb_id_up_up,wtqqb_id_dn_dn
      real(dp):: vsq(nf,nf),sscw(0:nf,1:nf)
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      vsq = 0._dp

      vsq(1,2) = Vud**2
      vsq(2,1) = vsq(1,2)
      vsq(1,4) = Vcd**2
      vsq(4,1) = vsq(1,4)
      vsq(2,3) = Vus**2
      vsq(3,2) = vsq(2,3)
      vsq(3,4) = Vcs**2
      vsq(4,3) = vsq(3,4)

C----set all elements to zero and compute invariants
      msq(:,:)=0._dp

      t1=-tt/ss
      t2=-uu/ss

c--- Fine-structure constant
c      alpha_em=1._dp/128._dp
c--- Fine-structure constant to reproduce Kuhn result
c     alpha_em=1._dp/126.3_dp
      alpha_em=esq/fourpi

c--- Logarithms from Eq. (2.8)
      DLs=alpha_em/(4._dp*pi)*real(log(cplx1(ss/wmass**2))**2,dp)
      SLs=alpha_em/(4._dp*pi)*real(log(cplx1(ss/wmass**2)),dp)
c--- Logarithm that enters LSC to account for Z/W mass difference
      Lm=log(zmass**2/wmass**2)
c--- Angular logarithm that enter SSC 
      SLt1=real(log(cplx1(t1)),dp)
      SLt2=real(log(cplx1(t2)),dp)
      
c--- change xw to sw2 for ease of translating formulae to code
      sw2=xw
      cw2=1._dp-sw2

c--- Definition of Casimir factor, without photon contribution
c--- i.e. Cew = (Iz)^2 + (Iw)^2
      Cew_uL=(3._dp*cw2-sw2)**2/(36._dp*sw2*cw2)+1._dp/(2._dp*sw2)
      Cew_ubL=Cew_uL
      Cew_uR=4._dp*sw2/(9._dp*cw2)
      Cew_ubR=Cew_uR
      Cew_dL=(3._dp*cw2+sw2)**2/(36._dp*sw2*cw2)+1._dp/(2._dp*sw2)
      Cew_dbL=Cew_dL
      Cew_dR=sw2/(9._dp*cw2)
      Cew_dbR=Cew_dR
c--- Definition of (Iz) alone
      Iz_uL=(3._dp*cw2-sw2)/sqrt(36._dp*sw2*cw2)
      Iz_ubL=Iz_uL
      Iz_uR=-sqrt(4._dp*sw2/(9._dp*cw2))
      Iz_ubR=Iz_uR
      Iz_dL=-(3._dp*cw2+sw2)/sqrt(36._dp*sw2*cw2)
      Iz_dbL=Iz_dL
      Iz_dR=sqrt(sw2/(9._dp*cw2))
      Iz_dbR=Iz_dR
c--- Calculation of (Iz)^2 alone
      Iz2_uL=Iz_uL**2
      Iz2_ubL=Iz_ubL**2
      Iz2_uR=Iz_uR**2
      Iz2_ubR=Iz_ubR**2
      Iz2_dL=Iz_dL**2
      Iz2_dbL=Iz_dbL**2
      Iz2_dR=Iz_dR**2
      Iz2_dbR=Iz_dbR**2

c--- GLUON-GLUON CHANNEL
c--- csame corresponds to the amplitude where t and t~ share the same chirality
c--- cdiff corresponds to the amplitude where t and t~ chiralities are different

c--- Leading soft-collinear logarithms, Eq. (3.6) and (3.7)
      csame_up=-0.5_dp*(Cew_uL+Cew_ubL+Cew_uR+Cew_ubR)*DLs
     &      +(Iz2_uL+Iz2_ubL+Iz2_uR+Iz2_ubR)*Lm*SLs
      cdiff_up=-0.5_dp*(Cew_uL+Cew_ubR+Cew_uR+Cew_ubL)*DLs
     &      +(Iz2_uL+Iz2_ubR+Iz2_uR+Iz2_ubL)*Lm*SLs
      csame_dn=-0.5_dp*(+Cew_dL+Cew_dbL+Cew_dR+Cew_dbR)*DLs
     &      +(+Iz2_dL+Iz2_dbL+Iz2_dR+Iz2_dbR)*Lm*SLs
      cdiff_dn=-0.5_dp*(+Cew_dL+Cew_dbR+Cew_dR+Cew_dbL)*DLs
     &      +(+Iz2_dL+Iz2_dbR+Iz2_dR+Iz2_dbL)*Lm*SLs

c--- Collinear and soft single logs, Eq. (4.2) and (4.6)
      csame_up=csame_up
     & +3._dp/2._dp*(Cew_uL+Cew_ubL)*SLs
     & +3._dp/2._dp*(Cew_uR+Cew_ubR)*SLs
      cdiff_up=cdiff_up
     & +3._dp/2._dp*(Cew_uL+Cew_ubR)*SLs
     & +3._dp/2._dp*(Cew_uR+Cew_ubL)*SLs
      csame_dn=csame_dn
     & +3._dp/2._dp*(Cew_dL+Cew_dbL)*SLs
     & +3._dp/2._dp*(Cew_dR+Cew_dbR)*SLs
     
      cdiff_dn=cdiff_dn
     & +3._dp/2._dp*(Cew_dL+Cew_dbR)*SLs
     & +3._dp/2._dp*(Cew_dR+Cew_dbL)*SLs
      
c--- @ 5 TEV: DLs --> -0.132 ,  SLs --> -0.032

c--- note that factor of two for (LO+EW)^2 is already included, so to reproduce
c--- LO result need to set csame = cdiff = 1
c      csame=1._dp
c      cdiff=1._dp
      
c--- color structure TaTb*TaTb
      wta_up=
     &+csame_up*(4._dp*t2/t1-8._dp*t2**2)
     
      wta_dn=
     &+csame_dn*(4._dp*t2/t1-8._dp*t2**2)
      
c--- color structure TbTa*TbTa
      wtb_up=
     &+csame_up*(4._dp*t1/t2-8._dp*t1**2)

      wtb_dn=
     &+csame_dn*(4._dp*t1/t2-8._dp*t1**2)
     
c--- color structure (TaTb+TbTa)*(TaTb+TbTa)
      wti_up=
     &+csame_up*(4._dp*t1/t2+4._dp*t2/t1)
     
      wti_dn=
     &+csame_dn*(4._dp*t1/t2+4._dp*t2/t1)      

      wtgg_up=2._dp*avegg*V/4._dp*gsq**2*(xn*(wta_up+wtb_up)-1._dp/xn*wti_up)
      wtgg_dn=2._dp*avegg*V/4._dp*gsq**2*(xn*(wta_dn+wtb_dn)-1._dp/xn*wti_dn)
      

c--- ANTIQUARK-QUARK CHANNEL(aqaq_iijj <=> qa_qa_iijj <=> a_iq_iq_ja_j -> 0)
c--- csame corresponds to the amplitude where q and t share the same chirality
c--- cdiff corresponds to the amplitude where q and t chiralities are different
c      DLs=0._dp
c      SLs=0._dp
c      SLt1=0._dp
c     SLt2=0._dp
      
c--- Subleading soft-collinear logarithms, Eq. (3.9) and (3.10)
c--- modified by log(s/Mw^2) -> log(s/Mz^2);
c--- to recover canonical result set Lm=0 here
      csame_up_up_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_uL*Iz_uL+Iz_ubL*Iz_ubL+Iz_uR*Iz_uR+Iz_ubR*Iz_ubR)
     &  -SLt2*(Iz_uL*Iz_ubL+Iz_ubL*Iz_uL+Iz_uR*Iz_ubR+Iz_ubR*Iz_uR))
      cdiff_up_up_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_uL*Iz_uR+Iz_ubL*Iz_ubR+Iz_uR*Iz_uL+Iz_ubR*Iz_ubL)
     &  -SLt2*(Iz_uL*Iz_ubR+Iz_ubL*Iz_uR+Iz_uR*Iz_ubL+Iz_ubR*Iz_uL))
      csame_up_dn_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_uL*Iz_dL+Iz_ubL*Iz_dbL+Iz_uR*Iz_dR+Iz_ubR*Iz_dbR)
     &  -SLt2*(Iz_uL*Iz_dbL+Iz_ubL*Iz_dL+Iz_uR*Iz_dbR+Iz_ubR*Iz_dR))
      cdiff_up_dn_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_uL*Iz_dR+Iz_ubL*Iz_dbR+Iz_uR*Iz_dL+Iz_ubR*Iz_dbL)
     &  -SLt2*(Iz_uL*Iz_dbR+Iz_ubL*Iz_dR+Iz_uR*Iz_dbL+Iz_ubR*Iz_dL))
      csame_dn_up_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_dL*Iz_uL+Iz_dbL*Iz_ubL+Iz_dR*Iz_uR+Iz_dbR*Iz_ubR)
     &  -SLt2*(Iz_dL*Iz_ubL+Iz_dbL*Iz_uL+Iz_dR*Iz_ubR+Iz_dbR*Iz_uR))
      cdiff_dn_up_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_dL*Iz_uR+Iz_dbL*Iz_ubR+Iz_dR*Iz_uL+Iz_dbR*Iz_ubL)
     &  -SLt2*(Iz_dL*Iz_ubR+Iz_dbL*Iz_uR+Iz_dR*Iz_ubL+Iz_dbR*Iz_uL))
      csame_dn_dn_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_dL*Iz_dL+Iz_dbL*Iz_dbL+Iz_dR*Iz_dR+Iz_dbR*Iz_dbR)
     &  -SLt2*(Iz_dL*Iz_dbL+Iz_dbL*Iz_dL+Iz_dR*Iz_dbR+Iz_dbR*Iz_dR))
      cdiff_dn_dn_ssc=-2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_dL*Iz_dR+Iz_dbL*Iz_dbR+Iz_dR*Iz_dL+Iz_dbR*Iz_dbL)
     &  -SLt2*(Iz_dL*Iz_dbR+Iz_dbL*Iz_dR+Iz_dR*Iz_dbL+Iz_dbR*Iz_dL))
      
      wtqqb_up_up=gsq**2*4._dp/9._dp
     & *((2._dp*csame_up+csame_up_up_ssc)*t2**2
     &  +(2._dp*cdiff_up+cdiff_up_up_ssc)*t1**2)

      wtqqb_up_dn=gsq**2*4._dp/9._dp
     & *((csame_up+csame_dn+csame_up_dn_ssc)*t2**2
     &  +(cdiff_up+cdiff_dn+cdiff_up_dn_ssc)*t1**2)


      wtqqb_dn_up=gsq**2*4._dp/9._dp
     & *((csame_dn+csame_up+csame_dn_up_ssc)*t2**2 
     &  +(cdiff_dn+cdiff_up+cdiff_dn_up_ssc)*t1**2)

      wtqqb_dn_dn=gsq**2*4._dp/9._dp
     & *((2._dp*csame_dn+csame_dn_dn_ssc)*t2**2 
     &  +(2._dp*cdiff_dn+cdiff_dn_dn_ssc)*t1**2)

      wtqqb_id_up_up=gsq**2*4._dp/9._dp
     & *((2._dp*csame_up+csame_up_up_ssc)*(t2/t1)**2
     &  +(2._dp*cdiff_up+cdiff_up_up_ssc)/t1**2
     &  -2._dp/xn*(2._dp*csame_up+csame_up_up_ssc)*t2**2/t1)

      wtqqb_id_dn_dn=gsq**2*4._dp/9._dp
     & *((2._dp*csame_dn+csame_dn_dn_ssc)*(t2/t1)**2
     &  +(2._dp*cdiff_dn+cdiff_dn_dn_ssc)/t1**2
     &  -2._dp/xn*(2._dp*csame_dn+csame_dn_dn_ssc)*t2**2/t1)

C --- ssc in W-exchange; it occurs when QCD(or EW) born has t-channel exchanged boson    
      sscw(:,:)=0._dp
      do j=1,nf-1
         do k=1,nf-1
            if(j == k) then
               sscw(j,k)=-SLs*SLt1/sw2
               sscw(j,k)=gsq**2*4._dp/9._dp
     &              *(2._dp*sscw(j,k)*((t2/t1)**2-2._dp/xn*t2**2/t1))
            else
               sscw(j,k)=+SLs*SLt2/sw2*vsq(j,k)
               sscw(j,k)=gsq**2*4._dp/9._dp
     &              *(2._dp*sscw(j,k)*(-2._dp/xn*t2**2/t1))
            end if
         end do
      end do

C---fill qb-q, gg and q-qb elements
      do j=0,nf
         do k=1,nf
            if ((j .eq. 0)) then
               if(mod(k,2) .eq. 0) then
                  msq(j,k)=wtgg_up/avegg
               else
                  msq(j,k)=wtgg_dn/avegg
               endif
            elseif ((mod(j,2) .eq. 0) .and. (mod(k,2) .eq. 0)) then
               if (j .eq. k) then
                  msq(j,k)=wtqqb_up_up+wtqqb_id_up_up
               else
                  msq(j,k)=wtqqb_up_up
               endif
            elseif ((mod(j,2) .eq. 0) .and. (mod(k,2) .ne. 0)) then
               msq(j,k)=wtqqb_up_dn
            elseif ((mod(j,2) .ne. 0) .and. (mod(k,2) .eq. 0)) then
               msq(j,k)=wtqqb_dn_up
            else
               if (j .eq. k) then
                  msq(j,k)=wtqqb_dn_dn+wtqqb_id_dn_dn
               else
                  msq(j,k)=wtqqb_dn_dn
               endif
            endif
            msq(j,k)=msq(j,k)+sscw(j,k)
         enddo
      enddo

      return
      end
