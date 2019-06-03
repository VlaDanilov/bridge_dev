      subroutine qqb_QQb_ew_sudakov(p,msq) 
      implicit none

************************************************************************
*     Author: J. M. Campbell                                           *
*     March, 2014.                                                     *
*                                                                      *
*     Calculate the leading EW Sudakov logs                            *
*     for the processes                                                *
*                                                                      *
*      g(-P1) + g(-P2)    --> Q(P3) + Qbar(P4)                         *
*      q(-P1) + qbar(-P2) --> Q(P3) + Qbar(P4)                         *
*                                                                      *
*     Results taken from                                               *
*   "One-loop leading logarithms in electroweak radiative corrections" *
*     A. Denner, S. Pozzorini                                          *
*     Eur. Phys. J. C 18, 461-480 (2001)                               *
*                                                                      *
************************************************************************
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcouple.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'breit.f'      
      integer j,k,cs
      logical first
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: wtgg,t1,t2,ro,cdiff,csame,wta,wtb,wti,
     & DLs,SLs,SLt,SLt1,SLt2,Lm,sw2,cw2,alpha_em,
     & Cew_uL,Cew_ubL,Cew_uR,Cew_ubR,
     & Cew_dL,Cew_dbL,Cew_dR,Cew_dbR,
     & Iz2_uL,Iz2_ubL,Iz2_uR,Iz2_ubR,
     & Iz2_dL,Iz2_dbL,Iz2_dR,Iz2_dbR,
     & Iz_uL,Iz_ubL,Iz_uR,Iz_ubR,
     & Iz_dL,Iz_dbL,Iz_dR,Iz_dbR,
     & csame_up,csame_dn,cdiff_up,cdiff_dn,
     & csame_up_ssc,csame_dn_ssc,cdiff_up_ssc,cdiff_dn_ssc,
     & wtqqb_up,wtqqb_dn,wtqbq_up,wtqbq_dn,
     & csame_bextra,cdiff_bextra,wt_bbbar,wt_bbarb,aemmz
      common/em/aemmz
      data first/.true./
      save first

      if (first) then
        first=.false.
        write(6,*) 'Heavy Quark mass:',mass2
      endif 

C----set all elements to zero and compute invariants
      msq(:,:)=0._dp
      call dotem(4,p,s)

      t1=-s(1,3)/s(1,2)
      t2=-s(2,3)/s(1,2)
      ro=4._dp*mass2**2/s(1,2)

c--- Fine-structure constant
c      alpha_em=1._dp/128._dp
c--- Fine-structure constant to reproduce Kuhn result
c      alpha_em=1._dp/126.3_dp
c--- alpha -> standard MCFM value
      alpha_em=aemmz

c--- Logarithms from Eq. (2.8)
      DLs=alpha_em/(4._dp*pi)*log(s(1,2)/wmass**2)**2
      SLs=alpha_em/(4._dp*pi)*log(s(1,2)/wmass**2)
c--- Logarithm that enters LSC to account for Z/W mass difference
      Lm=log(zmass**2/wmass**2)
c--- Angular logarithm that enter SSC 
      SLt1=log(t1)
      SLt2=log(t2)
      
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
      csame=-0.5_dp*(Cew_uL+Cew_ubL+Cew_uR+Cew_ubR)*DLs
     &      +(Iz2_uL+Iz2_ubL+Iz2_uR+Iz2_ubR)*Lm*SLs
      cdiff=-0.5_dp*(Cew_uL+Cew_ubR+Cew_uR+Cew_ubL)*DLs
     &      +(Iz2_uL+Iz2_ubR+Iz2_uR+Iz2_ubL)*Lm*SLs
c      csame=0._dp
c      cdiff=0._dp

c--- Collinear and soft single logs, Eq. (4.2) and (4.6)
c      csame=csame
c     & +(3._dp/2._dp*(Cew_uL+Cew_ubL)-2._dp*(mt/wmass)**2/8._dp/sw2)*SLs
c     & +(3._dp/2._dp*(Cew_uR+Cew_ubR)-4._dp*(mt/wmass)**2/8._dp/sw2)*SLs
c      cdiff=cdiff
c     & +(3._dp/2._dp*(Cew_uL+Cew_ubR)-3._dp*(mt/wmass)**2/8._dp/sw2)*SLs
c     & +(3._dp/2._dp*(Cew_uR+Cew_ubL)-3._dp*(mt/wmass)**2/8._dp/sw2)*SLs

      SLt=alpha_em/(4._dp*pi)*log(s(1,2)/mt**2) ! c.f. arXiv:0809.0800
c--- to recover "canonical" expression of the original reference, set SLt=SLs
c      SLt=SLs
c--- Collinear and soft single logs, Eq. (4.2) and (4.6)
      csame=csame
     & +(3._dp/2._dp*(Cew_uL+Cew_ubL)*SLs-2._dp*(mt/wmass)**2/8._dp/sw2*SLt)
     & +(3._dp/2._dp*(Cew_uR+Cew_ubR)*SLs-4._dp*(mt/wmass)**2/8._dp/sw2*SLt)
      cdiff=cdiff
     & +(3._dp/2._dp*(Cew_uL+Cew_ubR)*SLs-3._dp*(mt/wmass)**2/8._dp/sw2*SLt)
     & +(3._dp/2._dp*(Cew_uR+Cew_ubL)*SLs-3._dp*(mt/wmass)**2/8._dp/sw2*SLt)

c--- @ 5 TEV: DLs --> -0.132 ,  SLs --> -0.032

c--- note that factor of two for (LO+EW)^2 is already included, so to reproduce
c--- LO result need to set csame = cdiff = 1
c      csame=1._dp
c      cdiff=1._dp
      
c--- color structure TaTb*TaTb
      wta=
     &+csame*(4._dp*t2/t1-ro*(-6._dp*t2/t1+1._dp/t1**2)
     &       -8._dp*t2**2-0.5_dp*(ro/t1)**2)
     &+cdiff*(-ro*(2._dp*t2/t1-1._dp/t1**2)-0.5_dp*(ro/t1)**2)
     
c--- color structure TbTa*TbTa
      wtb=
     &+csame*(4._dp*t1/t2-ro*(-6._dp*t1/t2+1._dp/t2**2)
     &       -8._dp*t1**2-0.5_dp*(ro/t2)**2)
     &+cdiff*(-ro*(2._dp*t1/t2-1._dp/t2**2)-0.5_dp*(ro/t2)**2)
     
c--- color structure (TaTb+TbTa)*(TaTb+TbTa)
      wti=
     &+csame*(4._dp*t1/t2+4._dp*t2/t1+ro*(4._dp-t2/t1-t1/t2)/t1/t2
     &       -0.5_dp*(ro/t1/t2)**2)
     &+cdiff*(ro*(t2/t1+t1/t2)/t1/t2-0.5_dp*(ro/t1/t2)**2)

      wtgg=2._dp*avegg*V/4._dp*gsq**2*(xn*(wta+wtb)-1._dp/xn*wti)

c--- ANTIQUARK-QUARK CHANNEL
c--- csame corresponds to the amplitude where q and t share the same chirality
c--- cdiff corresponds to the amplitude where q and t chiralities are different
c      DLs=0._dp
c      SLs=0._dp
c      SLt1=0._dp
c      SLt2=0._dp
c--- Leading soft-collinear logarithms, Eq. (3.6) and (3.7)
      csame_up=-0.5_dp*(Cew_uL+Cew_ubL+Cew_uR+Cew_ubR)*2._dp*DLs
     &      +(Iz2_uL+Iz2_ubL+Iz2_uR+Iz2_ubR)*2._dp*Lm*SLs
      cdiff_up=-0.5_dp*(Cew_uL+Cew_ubR+Cew_uR+Cew_ubL)*2._dp*DLs
     &      +(Iz2_uL+Iz2_ubR+Iz2_uR+Iz2_ubL)*2._dp*Lm*SLs
      csame_dn=-0.5_dp*(Cew_uL+Cew_ubL+Cew_uR+Cew_ubR
     &                +Cew_dL+Cew_dbL+Cew_dR+Cew_dbR)*DLs
     &      +(Iz2_uL+Iz2_ubL+Iz2_uR+Iz2_ubR
     &       +Iz2_dL+Iz2_dbL+Iz2_dR+Iz2_dbR)*Lm*SLs
      cdiff_dn=-0.5_dp*(Cew_uL+Cew_ubR+Cew_uR+Cew_ubL
     &                +Cew_dL+Cew_dbR+Cew_dR+Cew_dbL)*DLs
     &      +(Iz2_uL+Iz2_ubR+Iz2_uR+Iz2_ubL
     &       +Iz2_dL+Iz2_dbR+Iz2_dR+Iz2_dbL)*Lm*SLs

      SLt=alpha_em/(4._dp*pi)*log(s(1,2)/mt**2) ! c.f. arXiv:0809.0800
c--- to recover "canonical" expression of the original reference, set SLt=SLs
c      SLt=SLs
c--- Collinear and soft single logs, Eq. (4.2) and (4.6)
      csame_up=csame_up
     & +(3._dp/2._dp*(Cew_uL+Cew_ubL)*2._dp*SLs-2._dp*(mt/wmass)**2/8._dp/sw2*SLt)
     & +(3._dp/2._dp*(Cew_uR+Cew_ubR)*2._dp*SLs-4._dp*(mt/wmass)**2/8._dp/sw2*SLt)
      cdiff_up=cdiff_up
     & +(3._dp/2._dp*(Cew_uL+Cew_ubR)*2._dp*SLs-3._dp*(mt/wmass)**2/8._dp/sw2*SLt)
     & +(3._dp/2._dp*(Cew_uR+Cew_ubL)*2._dp*SLs-3._dp*(mt/wmass)**2/8._dp/sw2*SLt)
      csame_dn=csame_dn
     & +(3._dp/2._dp*(Cew_uL+Cew_ubL+Cew_dL+Cew_dbL)*SLs
     &  -2._dp*(mt/wmass)**2/8._dp/sw2*SLt)
     & +(3._dp/2._dp*(Cew_uR+Cew_ubR+Cew_dR+Cew_dbR)*SLs
     &  -4._dp*(mt/wmass)**2/8._dp/sw2*SLt)
      cdiff_dn=cdiff_dn
     & +(3._dp/2._dp*(Cew_uL+Cew_ubR+Cew_dL+Cew_dbR)*SLs
     &  -3._dp*(mt/wmass)**2/8._dp/sw2*SLt)
     & +(3._dp/2._dp*(Cew_uR+Cew_ubL+Cew_dR+Cew_dbL)*SLs
     &  -3._dp*(mt/wmass)**2/8._dp/sw2*SLt)

c--- extra contribution (to _dn) if b-quark in initial state
      csame_bextra=-2._dp*(mt/wmass)**2/8._dp/sw2*SLt
      cdiff_bextra=-2._dp*(mt/wmass)**2/8._dp/sw2*SLt

c--- Subleading soft-collinear logarithms, Eq. (3.9) and (3.10)
c--- modified by log(s/Mw^2) -> log(s/Mz^2)
c--- to recover canonical result set Lm=0 here
      csame_up_ssc=2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_uL*Iz_uL+Iz_ubL*Iz_ubL+Iz_uR*Iz_uR+Iz_ubR*Iz_ubR)
     &  -Slt2*(Iz_uL*Iz_ubL+Iz_ubL*Iz_uL+Iz_uR*Iz_ubR+Iz_ubR*Iz_uR))
      cdiff_up_ssc=2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_uL*Iz_uR+Iz_ubL*Iz_ubR+Iz_uR*Iz_uL+Iz_ubR*Iz_ubL)
     &  -Slt2*(Iz_uL*Iz_ubR+Iz_ubL*Iz_uR+Iz_uR*Iz_ubL+Iz_ubR*Iz_uL))
      csame_dn_ssc=2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_dL*Iz_uL+Iz_dbL*Iz_ubL+Iz_dR*Iz_uR+Iz_dbR*Iz_ubR)
     &  -Slt2*(Iz_dL*Iz_ubL+Iz_dbL*Iz_uL+Iz_dR*Iz_ubR+Iz_dbR*Iz_uR))
      cdiff_dn_ssc=2._dp*(SLs-alpha_em/(4._dp*pi)*Lm)
     & *(SLt1*(Iz_dL*Iz_uR+Iz_dbL*Iz_ubR+Iz_dR*Iz_uL+Iz_dbR*Iz_ubL)
     &  -Slt2*(Iz_dL*Iz_ubR+Iz_dbL*Iz_uR+Iz_dR*Iz_ubL+Iz_dbR*Iz_uL))

      wtqbq_up=gsq**2*4._dp/9._dp
     & *((csame_up+csame_up_ssc)*(t1**2+ro/4._dp)
     &  +(cdiff_up+cdiff_up_ssc)*(t2**2+ro/4._dp))
      wtqbq_dn=gsq**2*4._dp/9._dp
     & *((csame_dn+csame_dn_ssc)*(t1**2+ro/4._dp)
     &  +(cdiff_dn+cdiff_dn_ssc)*(t2**2+ro/4._dp))

c--- to obtain quark-antiquark, interchange t1 and t2, which flips sign of _ssc terms
      wtqqb_up=gsq**2*4._dp/9._dp
     & *((csame_up-csame_up_ssc)*(t2**2+ro/4._dp)
     &  +(cdiff_up-cdiff_up_ssc)*(t1**2+ro/4._dp))
      wtqqb_dn=gsq**2*4._dp/9._dp
     & *((csame_dn-csame_dn_ssc)*(t2**2+ro/4._dp)
     &  +(cdiff_dn-cdiff_dn_ssc)*(t1**2+ro/4._dp))

c--- add extra contribution to b-quark initial states
      wt_bbbar=wtqqb_dn+gsq**2*4._dp/9._dp
     & *((csame_bextra)*(t2**2+ro/4._dp)
     &  +(cdiff_bextra)*(t1**2+ro/4._dp))
      wt_bbarb=wtqbq_dn+gsq**2*4._dp/9._dp
     & *((csame_bextra)*(t2**2+ro/4._dp)
     &  +(cdiff_bextra)*(t1**2+ro/4._dp))

C---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      k=-j
      if ((j .eq. 0) .and. (k.eq.0)) then
          msq(j,k)=wtgg
      elseif ((j .gt. 0) .and. (k.lt.0)) then
          if (mod(j,2) .eq. 0) then
            msq(j,k)=wtqqb_up
          else
            if (j .eq. 5) then
              msq(j,k)=wt_bbbar
            else
              msq(j,k)=wtqqb_dn
            endif
          endif
      elseif ((j .lt. 0) .and. (k.gt.0)) then
          if (mod(k,2) .eq. 0) then
            msq(j,k)=wtqbq_up
          else
            if (k .eq. 5) then
              msq(j,k)=wt_bbarb
            else
            msq(j,k)=wtqbq_dn
            endif
          endif
      endif
      enddo

      return
      end
