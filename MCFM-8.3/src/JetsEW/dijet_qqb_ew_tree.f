      subroutine dijet_qqb_ew_tree(msq,ss,tt,uu)
c---- Jia Zhou, Feb. 2015
c---- computes matrix element for LO dijet production via
c---- weak or weak-strong interaction:
c----  q(p1) + q~(p2) -> q(p3) + q~(p4)
c----  ss=(p1+p2)^2, tt=(p1+p3)^2, uu=(p2+p3)^2
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      real(dp):: msq(nf,nf),ss,tt,uu
      real(dp):: brn_mix(0:3,nf,nf),brn_wk(nf,nf),ai(nf),af(nf),
     &     vi(nf),vf(nf),vsq(nf,nf),cwsq,mz,mw
      real(dp):: prop_zxz_sxs,prop_zxz_sxt,prop_zxz_txt,
     &     prop_zxr_sxs,prop_zxr_sxt,prop_zxr_txs,prop_zxr_txt,
     &     prop_rxr_sxs,prop_rxr_sxt,prop_rxr_txt,
     &     prop_zxw_sxt,prop_wxw_txt,prop_rxw_sxt
      integer j,k
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb
      real(dp):: sxtfac
      parameter(sxtfac=-1._dp/3._dp)

      vi(1:nf) = (l(1:nf) + r(1:nf))/2._dp
      vf(1:nf) = vi(1:nf)
      ai(1:nf) = (l(1:nf) - r(1:nf))/2._dp
      af(1:nf) = ai(1:nf)

      cwsq = 1._dp/(2._dp*xw)

      vsq = 0._dp
      vsq(1,2) = Vud**2
      vsq(2,1) = vsq(1,2)
      vsq(1,4) = Vcd**2
      vsq(4,1) = vsq(1,4)
      vsq(2,3) = Vus**2
      vsq(3,2) = vsq(2,3)
      vsq(3,4) = Vcs**2
      vsq(4,3) = vsq(3,4)

      mz = zmass
      mw = wmass

      prop_zxz_sxs = 1._dp/((ss - mz**2)**2 + mz**2*zwidth**2)
      prop_zxz_sxt = (ss - mz**2)*(tt - mz**2) + mz**2*zwidth**2
      prop_zxz_sxt = prop_zxz_sxt/(prop_zxz_sxt**2 
     &     + mz**2*zwidth**2*(ss - tt)**2)
      prop_zxz_txt = 1._dp/((tt - mz**2)**2 + mz**2*zwidth**2)

      prop_zxr_sxs = (ss - mz**2)/((ss - mz**2)**2 + mz**2*zwidth**2)/ss
      prop_zxr_sxt = (ss - mz**2)/((ss - mz**2)**2 + mz**2*zwidth**2)/tt
      prop_zxr_txs = (tt - mz**2)/((tt - mz**2)**2 + mz**2*zwidth**2)/ss
      prop_zxr_txt = (tt - mz**2)/((tt - mz**2)**2 + mz**2*zwidth**2)/tt

      prop_rxr_sxs = 1._dp/ss**2
      prop_rxr_sxt = 1._dp/ss/tt
      prop_rxr_txt = 1._dp/tt**2

      prop_zxw_sxt = (ss - mz**2)*(tt - mw**2) + mz*mw*zwidth*wwidth
      prop_zxw_sxt = prop_zxw_sxt/(prop_zxw_sxt**2 
     &     + (mz*zwidth*(tt - mw**2) - mw*wwidth*(ss - mz**2))**2)
      prop_wxw_txt = 1._dp/((tt - mw**2)**2 + mw**2*wwidth**2)
      prop_rxw_sxt = (tt - mw**2)/((tt - mw**2)**2 + mw**2*wwidth**2)/ss

c--- Added by JC: accounts for 1 fermion loop (-1) and
c--- less factors of Nc (one vs. two in the square) ---> -1/Nc
      prop_zxz_sxt=prop_zxz_sxt*sxtfac
      prop_zxr_sxt=prop_zxr_sxt*sxtfac
      prop_zxr_txs=prop_zxr_txs*sxtfac
      prop_rxr_sxt=prop_rxr_sxt*sxtfac
      prop_zxw_sxt=prop_zxw_sxt*sxtfac
      prop_rxw_sxt=prop_rxw_sxt*sxtfac
c--- Added by JC

C --- weak-QCD mixed born; only interference between s- and t-exchange amp
!do NOT use this subroutine when flg = .false. in 'dijet_qqb_mix.f'
c      call dijet_qqb_ii_jj_mix(brn_mix,ss,tt,uu)       
C --- '4' is the color structure
c      brn_mix(1:2,:,:) = 4._dp*brn_mix(1:2,:,:)

      brn_mix = 0._dp
C --- pure weak mediated born
      brn_wk = 0._dp

      do j = 1,nf
         do k = 1,nf

C --- sxs: (1) ZxZ (2) Zxr(rxZ) (3) rxr, where 'r' denotes photon exchange

C --- in (1) ZxZ
      brn_wk(j,k) = prop_zxz_sxs*(
     &           -8._dp*(
     &           - (tt**2 + uu**2)*(ai(j)**2 + vi(j)**2)*(af(k)**2 
     &           + vf(k)**2) 
     &           + 4._dp*(tt**2 - uu**2)*(ai(j)*vi(j)*af(k)*vf(k))
     &           )
     &           )

C --- in (2) Zxr .and. rxZ
      brn_wk(j,k) = brn_wk(j,k)
     &     + 2._dp*prop_zxr_sxs*(
     &     - 8._dp*Q(j)*Q(k)*(
     &     - (tt**2 + uu**2)*vi(j)*vf(k)
     &     + (tt**2 - uu**2)*ai(j)*af(k)
     &     )
     &     )

C --- in (3) rxr
      brn_wk(j,k) = brn_wk(j,k)
     &     + prop_rxr_sxs*Q(j)**2*Q(k)**2*8._dp*(tt**2 + uu**2)


C --- sxt .and. txs: (1) ZxZ (2) Zxr(rxZ) (3) rxr (4) ZxW(WxZ) (5) rxW(Wxr)
C --- txt: (1) ZxZ (2) Zxr(rxZ) (3) rxr (4) WxW
c      if(.false.) then
      if(j == k) then
         brn_wk(j,k) = brn_wk(j,k)
C --- ZxZ (sxt) + (txs) + (txt)
     &        + 2._dp*prop_zxz_sxt*(
     &        - 8._dp*uu**2*(
     &        + (ai(j)**2 + vi(j)**2)*(af(k)**2 + vf(k)**2)
     &        + 4._dp*ai(j)*af(k)*vi(j)*vf(k)
     &        )
     &        )
     &        + prop_zxz_txt*(
     &        - 8._dp*(
     &        - (ss**2 + uu**2)*(ai(j)**2 + vi(j)**2)*(af(k)**2 
     &        + vf(k)**2)
     &        + 4._dp*(ss**2 - uu**2)*(ai(j)*vi(j)*af(k)*vf(k))
     &        )
     &        )
C --- Zxr(sxt) + rxZ(txs)
     &        + 2._dp*prop_zxr_sxt*(
     &        - 8._dp*uu**2*Q(j)*Q(k)*(ai(j)*af(k) + vi(j)*vf(k))
     &        )
C --- rxZ(sxt) + Zxr(txs)
     &        + 2._dp*prop_zxr_txs*(
     &        - 8._dp*uu**2*Q(j)*Q(k)*(ai(j)*af(k) + vi(j)*vf(k))
     &        )
C --- Zxr(txt) + rxZ(txt)
     &        + 2._dp*prop_zxr_txt*(
     &        - 8._dp*Q(j)*Q(k)*(
     &        - (ss**2 + uu**2)*vi(j)*vf(k)
     &        + (ss**2 - uu**2)*ai(j)*af(k)
     &        )
     &        )
C --- rxr (sxt) + (txs) + (txt)
     &        + 2._dp*prop_rxr_sxt*(-8._dp*uu**2)*Q(j)**2*Q(k)**2
     &        + prop_rxr_txt*8._dp*(ss**2 + uu**2)*Q(j)**2*Q(k)**2


C --- add in mixed born as the above subroutine 'dijet_qqb_ii_jj_mix' uncalled
         brn_mix(1,j,k) = - 8._dp*uu**2*(
     &        + (ai(j)*af(k) + vi(j)*vf(k))*prop_zxr_sxt
     &        + Q(j)*Q(k)*prop_rxr_sxt
     &        )

         brn_mix(2,j,k) = - 8._dp*uu**2*(
     &        +  (ai(j)*af(k) + vi(j)*vf(k))*prop_zxr_txs
     &        + Q(j)*Q(k)*prop_rxr_sxt
     &        )

      else
         brn_wk(j,k) = brn_wk(j,k)
C --- ZxW + WxZ
     &        + 2._dp*prop_zxw_sxt*(
     &        - 4._dp*cwsq*vsq(j,k)*(ai(j) + vi(j))*(af(k) + vf(k))*uu**2
     &        )
C --- rxW + Wxr
     &        + 2._dp*prop_rxw_sxt*(
     &        - 4._dp*cwsq*vsq(j,k)*Q(j)*Q(k)*uu**2
     &        )
C --- WxW
     &        + prop_wxw_txt*(
     &        + 4._dp*cwsq**2*vsq(j,k)**2*uu**2
     &        )

C --- add in mixed born
         brn_mix(2,j,k) = - 4._dp*uu**2*cwsq*vsq(j,k)*prop_rxw_sxt
      end if
c      end if
c         msq(j,k) = brn_mix(1,j,k) + brn_mix(2,j,k) 
c     &        + 2._dp*brn_wk(j,k)*esq**2*aveqq
c--- John:  I have modified the factors in the lines below;
c---   I believe the second factor (xn**2/2._dp) is correct
c---   but am not sure about the factor on the first line
          msq(j,k) = aveqq*(
     &     + 4._dp*gsq*esq*(brn_mix(1,j,k) + brn_mix(2,j,k))*2._dp*xn
     &     + esq**2*brn_wk(j,k)*(xn**2)
     &     )
      end do
      end do

      end subroutine dijet_qqb_ew_tree



C******************************************************************C
C** use product of spinor to calculate the amp as a double-check **C
C******************************************************************C
      subroutine dijet_qqb_ew_tree2(msq,p)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'cplx.h'
      real(dp):: msq(nf,nf),p(mxpart,4),ss,tt,uu,mz,mw,vsq(nf,nf),
     &     msq_mix(nf,nf),msq_wk(nf,nf),cwsq
      complex(dp):: prop_z_s,prop_z_t,prop_w,
     &     LL_s,RR_s,LR_s,RL_s,
     &     LL_t,RR_t,LR_t,RL_t,
     &     LL_s_g(nf,nf),RR_s_g(nf,nf),LR_s_g(nf,nf),RL_s_g(nf,nf),
     &     LL_t_g(nf,nf),RR_t_g(nf,nf),LR_t_g(nf,nf),RL_t_g(nf,nf),
     &     LL_s_z(nf,nf),RR_s_z(nf,nf),LR_s_z(nf,nf),RL_s_z(nf,nf),
     &     LL_t_z(nf,nf),RR_t_z(nf,nf),LR_t_z(nf,nf),RL_t_z(nf,nf),
     &     LL_s_r(nf,nf),RR_s_r(nf,nf),LR_s_r(nf,nf),RL_s_r(nf,nf),
     &     LL_t_r(nf,nf),RR_t_r(nf,nf),LR_t_r(nf,nf),RL_t_r(nf,nf),
     &     LL_t_w(nf,nf),RR_t_w(nf,nf),LR_t_w(nf,nf),RL_t_w(nf,nf)
      integer j,k
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

      cwsq = 1._dp/(2._dp*xw)

      call dotem(4,p,s)

      ss = s(1,2)
      tt = s(1,3)
      uu = s(1,4)

      mz = zmass
      mw = wmass

      prop_z_s = 1._dp/cplx2((ss - mz**2),mz*zwidth)
      prop_z_t = 1._dp/cplx2((tt - mz**2),mz*zwidth)
      prop_w = 1._dp/cplx2((tt - mw**2),mw*wwidth) 

      call spinoru(4,p,za,zb)

C -- helicity correspondence: 
C --- s-channel: 1, 3; i.e., LL_s denotes L(1) & L(3)
      LL_s = 2._dp*za(3,2)*zb(1,4)
      RR_s = 2._dp*zb(3,2)*za(1,4)
      LR_s = 2._dp*zb(3,1)*za(2,4)
      RL_s = 2._dp*za(3,1)*zb(2,4)
C --- t-channel: 1, 2
C --- crossing relation w/ s-channel: 2 <-> 3
      LL_t = 2._dp*za(2,3)*zb(1,4)
      RR_t = 2._dp*zb(2,3)*za(1,4)
      LR_t = 2._dp*zb(2,1)*za(3,4)
      RL_t = 2._dp*za(2,1)*zb(3,4)

C --- separate mixed and pure weak processes
      LL_s_g = 0._dp
      RR_s_g = 0._dp
      LR_s_g = 0._dp
      RL_s_g = 0._dp
      LL_s_z = 0._dp
      RR_s_z = 0._dp
      LR_s_z = 0._dp
      RL_s_z = 0._dp
      LL_s_r = 0._dp
      RR_s_r = 0._dp
      LR_s_r = 0._dp
      RL_s_r = 0._dp
      LL_t_g = 0._dp
      RR_t_g = 0._dp
      LR_t_g = 0._dp
      RL_t_g = 0._dp
      LL_t_z = 0._dp
      RR_t_z = 0._dp
      LR_t_z = 0._dp
      RL_t_z = 0._dp
      LL_t_r = 0._dp
      RR_t_r = 0._dp
      LR_t_r = 0._dp
      RL_t_r = 0._dp
      LL_t_w = 0._dp
      RR_t_w = 0._dp
      LR_t_w = 0._dp
      RL_t_w = 0._dp
      do j = 1,nf
         do k = 1,nf

            LL_s_g(j,k) = gsq*LL_s/ss
            RR_s_g(j,k) = gsq*RR_s/ss
            LR_s_g(j,k) = gsq*LR_s/ss
            RL_s_g(j,k) = gsq*RL_s/ss

            LL_s_z(j,k) = esq*prop_z_s*LL_s*L(j)*L(k)
            RR_s_z(j,k) = esq*prop_z_s*RR_s*R(j)*R(k)
            LR_s_z(j,k) = esq*prop_z_s*LR_s*L(j)*R(k)
            RL_s_z(j,k) = esq*prop_z_s*RL_s*R(j)*L(k)

            LL_s_r(j,k) = esq/ss*Q(j)*Q(k)*LL_s
            RR_s_r(j,k) = esq/ss*Q(j)*Q(k)*RR_s
            LR_s_r(j,k) = esq/ss*Q(j)*Q(k)*LR_s
            RL_s_r(j,k) = esq/ss*Q(j)*Q(k)*RL_s

C --- identical initial and final states
            if(j == k) then
               LL_t_g(j,k) = gsq*LL_t/tt
               RR_t_g(j,k) = gsq*RR_t/tt
               LR_t_g(j,k) = gsq*LR_t/tt
               RL_t_g(j,k) = gsq*RL_t/tt

               LL_t_z(j,k) = esq*prop_z_t*LL_t*L(j)*L(k)
               RR_t_z(j,k) = esq*prop_z_t*RR_t*R(j)*R(k)
               LR_t_z(j,k) = esq*prop_z_t*LR_t*L(j)*R(k)
               RL_t_z(j,k) = esq*prop_z_t*RL_t*R(j)*L(k)

               LL_t_r(j,k) = esq/tt*Q(j)*Q(k)*LL_t
               RR_t_r(j,k) = esq/tt*Q(j)*Q(k)*RR_t
               LR_t_r(j,k) = esq/tt*Q(j)*Q(k)*LR_t
               RL_t_r(j,k) = esq/tt*Q(j)*Q(k)*RL_t

            else
C --- no right-handed portion
               LL_t_w(j,k) = esq*prop_w*cwsq*vsq(j,k)*LL_t
c               RR_t_w(j,k) = esq*prop_w*cwsq*vsq(j,k)*RR_t
c               LR_t_w(j,k) = esq*prop_w*cwsq*vsq(j,k)*LR_t
c               RL_t_w(j,k) = esq*prop_w*cwsq*vsq(j,k)*RL_t
            end if

            msq_mix(j,k) = 2._dp*4._dp*(
     &           + real((LL_t_z(j,k) + LL_t_r(j,k)
     &           + LL_t_w(j,k))*conjg(LL_s_g(j,k)),dp)
     &           + real((RR_t_z(j,k) + RR_t_r(j,k)
     &           )*conjg(RR_s_g(j,k)),dp)
     &           + real((LL_s_z(j,k) + LL_s_r(j,k)
     &           )*conjg(LL_t_g(j,k)),dp)
     &           + real((RR_s_z(j,k) + RR_s_r(j,k)
     &           )*conjg(RR_t_g(j,k)),dp)
     &           )

            msq_wk(j,k) = 2._dp*(
C --- sxs
     &           + real((LL_s_z(j,k) + LL_s_r(j,k))*conjg(LL_s_z(j,k) 
     &           + LL_s_r(j,k)),dp)
     &           + real((RR_s_z(j,k) + RR_s_r(j,k))*conjg(RR_s_z(j,k) 
     &           + RR_s_r(j,k)),dp)
     &           + real((LR_s_z(j,k) + LR_s_r(j,k))*conjg(LR_s_z(j,k) 
     &           + LR_s_r(j,k)),dp)
     &           + real((RL_s_z(j,k) + RL_s_r(j,k))*conjg(RL_s_z(j,k) 
     &           + RL_s_r(j,k)),dp)
C --- sxt + txs
     &           + 2._dp*(
     &           + real((LL_s_z(j,k) + LL_s_r(j,k))*conjg(LL_t_z(j,k)
     &           + LL_t_r(j,k) + LL_t_w(j,k)),dp)
     &           + real((RR_s_z(j,k) + RR_s_r(j,k))*conjg(RR_t_z(j,k)
     &           + RR_t_r(j,k)),dp)
     &           )
C --- txt
     &           + real((LL_t_z(j,k) + LL_t_r(j,k) 
     &           + LL_t_w(j,k))*conjg(LL_t_z(j,k) + LL_t_r(j,k) 
     &           + LL_t_w(j,k)),dp)
     &           + real((RR_t_z(j,k) + RR_t_r(j,k))*conjg(RR_t_z(j,k) 
     &           + RR_t_r(j,k)),dp)
     &           + real((LR_t_z(j,k) + LR_t_r(j,k))*conjg(LR_t_z(j,k) 
     &           + LR_t_r(j,k)),dp)
     &           + real((RL_t_z(j,k) + RL_t_r(j,k))*conjg(RL_t_z(j,k) 
     &           + RL_t_r(j,k)),dp)
     &           )

            msq(j,k) = msq_mix(j,k) + msq_wk(j,k)

         end do
      end do

      msq = msq*aveqq

      end subroutine dijet_qqb_ew_tree2
      
