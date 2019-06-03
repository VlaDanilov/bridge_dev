      subroutine qqb_twojet_mix(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'sprods_com.f'
      include 'msq_mix.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'masses.f'
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),ss,tt,uu,
     .     qaii_jj(0:3,nf,nf),qaij_ij(0:3,nf,nf),
     .     qqij_ij(0:3,nf,nf),aqii_jj(0:3,nf,nf),
     .     aqij_ij(0:3,nf,nf),aaij_ij(0:3,nf,nf)
      integer j,k
      real(dp):: mss(2),mts,mst,mtt(2),sxs

C --- cross check with Doreen at one phase space point
      if(.false.) then
      p(1,:) = (/0._dp, 0._dp, -50.123897978046799_dp, 
     &     -50.123897978046799_dp/)
      p(2,:) = (/0._dp, 0._dp, 6880.9140123772686_dp, -6880.9140123772686_dp/)
      p(3,:) = (/-18.352460460194017_dp, -9.1444724432301498_dp, 
     &     48.011003809174781_dp, 52.206232079832716_dp/)
      p(4,:) = (/18.352460460194035_dp, 9.1444724432301534_dp,  
     &     -6878.8011182083965_dp, 6878.8316782754819_dp/)
      end if

      call dotem(4,p,s)
      ss = s(1,2)
      tt = s(1,3)
      uu = s(2,3)

c      print*, "ss,tt,uu: ", ss,tt,uu

      call qqb_twojet_ii_jj_mix(qaii_jj,ss,tt,uu)
      call qqb_twojet_ii_jj_mix(qaij_ij,tt,ss,uu)
      call qqb_twojet_ii_jj_mix(qqij_ij,tt,uu,ss)

c      call qqb_twojet_ii_jj_mix(aqii_jj,ss,tt,uu)
c      call qqb_twojet_ii_jj_mix(aqij_ij,tt,ss,uu)
c      call qqb_twojet_ii_jj_mix(aaij_ij,tt,uu,ss)
c--- Use symmetry relations instead of calls, for speed
      aqii_jj(:,:,:)=qaii_jj(:,:,:)
      aqij_ij(:,:,:)=qaij_ij(:,:,:)
      aaij_ij(:,:,:)=qqij_ij(:,:,:)

c-- for now, no return value; real purpose of this routine is
c-- to fill the common block msq_mix
      msq=0._dp

      msq_mix = 0._dp

      do j = fn,nf
         do k = fn,nf

C --- qa
            if((j > 0) .and. (k < 0)) then
               if(j == -k) then
                  msq_mix(:,j,k) = qaii_jj(:,j,1) + qaii_jj(:,j,2) 
     &                           + qaii_jj(:,j,3) + qaii_jj(:,j,4)
               else
                  msq_mix(:,j,k) = qaij_ij(:,j,-k)
               end if

C --- qq
            else if((j > 0) .and. (k > 0)) then
               if(j == k)  then
                  msq_mix(:,j,k) = qqij_ij(:,j,k)*half
               else
                  msq_mix(:,j,k) = qqij_ij(:,j,k)
               end if

C --- aq
            else if((j < 0) .and. (k > 0)) then
               if(k == -j) then
                  msq_mix(:,j,k) = aqii_jj(:,k,1) + aqii_jj(:,k,2) 
     &                           + aqii_jj(:,k,3) + aqii_jj(:,k,4)
               else
                  msq_mix(:,j,k) = aqij_ij(:,-j,k)
               end if

C --- aa
            else if((j < 0) .and. (k < 0)) then
               if(j == k) then
                  msq_mix(:,j,k) = aaij_ij(:,-j,-k)*half
               else
                  msq_mix(:,j,k) = aaij_ij(:,-j,-k)
               end if
            end if
         end do
      end do

C--- multiply by factor of N*CF=4, which is the natural factor
c--- for the sxt interference terms that have a collinear singularity
      msq_mix=msq_mix*(xn*Cf)


c--- debugging code
      if(.false.) then
        print*, "phase point: "
        print*, "p1: ", p(1,:)
        print*, "p2: ", p(2,:)
        print*, "p3: ", p(3,:)
        print*, "p4: ", p(4,:)
        print*, "mixed: ", msq_mix(0,2,1)/gsq*(4._dp*pi*0.118_dp),
     &       msq_mix(2,2,1)/gsq*(4._dp*pi*0.118_dp), esq
        print*, "strong: ", sxs(tt,uu,ss)*(4._dp*pi*0.118_dp)**2*aveqq*2._dp
        pause
      end if

      end subroutine qqb_twojet_mix




C --- summing the following four subroutines up
      subroutine qqb_twojet_ii_jj_mix(msq_cs,ss,tt,uu)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      real(dp):: msq_cs(0:3,nf,nf),ss,tt,uu,
     .     msq_cs1(0:3,1:nf,fn:-1),msq_cs2(0:3,1:nf,fn:-1),
     .     msq_cs3(0:3,1:nf,fn:-1),msq_cs4(0:3,1:nf,fn:-1)
      integer j,k
      logical flg/.false./

      call qqb_twojet_uub_mix(msq_cs2,ss,tt,uu,flg)
      call qqb_twojet_ddb_mix(msq_cs1,ss,tt,uu,flg)
      call qqb_twojet_ccb_mix(msq_cs4,ss,tt,uu,flg)
      call qqb_twojet_ssb_mix(msq_cs3,ss,tt,uu,flg)

      msq_cs = 0._dp
      do j = 1,nf-1
         k = -j
         msq_cs(:,j,1) = msq_cs1(:,j,k)
         msq_cs(:,j,2) = msq_cs2(:,j,k)
         msq_cs(:,j,3) = msq_cs3(:,j,k)
         msq_cs(:,j,4) = msq_cs4(:,j,k)
      end do

      msq_cs = msq_cs*aveqq

c--- DEBUG: Jia check
c      msq_cs(1:2,:,:)=0._dp

c      msq_cs(0,:,:)=0._dp ! DEBUG
c      msq_cs(1,:,:)=0._dp ! DEBUG
c      msq_cs(2,:,:)=0._dp ! DEBUG
c      msq_cs(3,:,:)=0._dp ! DEBUG

      
      end subroutine qqb_twojet_ii_jj_mix





      subroutine qqb_twojet_uub_mix(msq_cs,ss,tt,uu,flg)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      real(dp):: msq_cs(0:3,1:nf,fn:-1),ss,tt,uu,
     .     vi(nf),vf(nf),ai(nf),af(nf),Qf(nf),cwsq,mz,mw,
     .     propz_s,propz_t,propw,mss(2),mst,mts,mtt(2)
      integer j,k
      logical flg
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      mz = zmass
      mw = wmass

C --- cwsq is W coupling to fermion
      cwsq = 1._dp/4._dp/xw


      propz_s = ss*(ss-mz**2)/((ss-mz**2)**2 + zwidth**2*mz**2)
      propz_t = tt*(tt-mz**2)/((tt-mz**2)**2 + zwidth**2*mz**2)
      propw = tt*(tt-mw**2)/((tt-mw**2)**2 + wwidth**2*mw**2)

      vi(1:nf) = (l(1:nf) + r(1:nf))/2._dp
      vf(1:nf) = vi(1:nf)
      ai(1:nf) = (l(1:nf) - r(1:nf))/2._dp
      af(1:nf) = ai(1:nf)

C --- flag determines whether photon mediated amplitude is included

      Qf(1:nf) = 0._dp

      if(flg) Qf(1:nf) = Q(1:nf)

C --- color structure cs1, cs2, cs3
C --- cs1 = (Nc^2 - 1)/4; cs2 = -(Nc^2 -1)/(4*Nc); cs3 = (Nc^2 - 1)^2/(4*Nc)

      call qqb_twojet_born_mix(mss,mst,mts,mtt,ss,tt,uu)

C --- msq(0:3,:,:)
C --- first component corresponds to 0 = sxs, 1 = sxt, 2 = txs, 3 = txt
      msq_cs = 0._dp
C --- non-zero components
c      msq_cs(1,1,-1) = gwsq*mst*propw
c      msq_cs(1,2,-2) = (vi(2)**2 + ai(2)**2)*(mst*propz_t + mts*propz_s)
c     .     + Q(2)**2*(mst + mts)
c      msq_cs(1,3,-3) = msq_cs(1,1,-1)
c      msq_cs(2,2,-2) = (vi(2)**2*mtt(1) + ai(2)**2*mtt(2))*propz_t
c     .     + Q(2)**2*mss(1)

C --- add one more column 
      msq_cs(1,2,-2) = (vi(2)**2 + ai(2)**2)*mst*propz_s + Qf(2)**2*mst
      msq_cs(2,1,-1) = cwsq*mts*propw*Vud**2
      msq_cs(2,3,-3) = cwsq*mts*propw*Vus**2
      msq_cs(2,2,-2) = (vi(2)**2 + ai(2)**2)*mts*propz_t + Qf(2)**2*mts
      msq_cs(3,2,-2) = (vi(2)**2*mtt(1) + ai(2)**2*mtt(2))*propz_t
     .     + Qf(2)**2*mtt(1)

      do j = 1,nf-1
         k = -j
         msq_cs(0,j,k) = 
     .        + (vi(j)*vf(2)*mss(1) + ai(j)*af(2)*mss(2))*propz_s
     .        + Qf(j)*Qf(2)*mss(1)
      end do

C --- msq_cs has order of alpha_s*alpha
      msq_cs = 2._dp*msq_cs*gsq*esq

      end subroutine qqb_twojet_uub_mix





      subroutine qqb_twojet_ddb_mix(msq_cs,ss,tt,uu,flg)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      real(dp):: msq_cs(0:3,1:nf,fn:-1),ss,tt,uu,
     .     vi(nf),vf(nf),ai(nf),af(nf),Qf(nf),cwsq,mz,mw,
     .     propz_s,propz_t,propw,mss(2),mst,mts,mtt(2)
      integer j,k
      logical flg
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      mz = zmass
      mw = wmass

C --- cwsq is W coupling to fermion
      cwsq = 1._dp/4._dp/xw


      propz_s = ss*(ss-mz**2)/((ss-mz**2)**2 + zwidth**2*mz**2)
      propz_t = tt*(tt-mz**2)/((tt-mz**2)**2 + zwidth**2*mz**2)
      propw = tt*(tt-mw**2)/((tt-mw**2)**2 + wwidth**2*mw**2)
! DEBUG
c      propw = propw + 1._dp

      vi(1:nf) = (l(1:nf) + r(1:nf))/2._dp
      vf(1:nf) = vi(1:nf)
      ai(1:nf) = (l(1:nf) - r(1:nf))/2._dp
      af(1:nf) = ai(1:nf)

C --- flag determines whether photon mediated amplitude is included

      Qf(1:nf) = 0._dp

      if(flg) Qf(1:nf) = Q(1:nf)


C --- color structure cs1, cs2, cs3
C --- cs1 = (Nc^2 - 1)/4; cs2 = -(Nc^2 -1)/(4*Nc); cs3 = (Nc^2 - 1)^2/(4*Nc)

      call qqb_twojet_born_mix(mss,mst,mts,mtt,ss,tt,uu)

C --- msq_cs(0:2,:,:)
C --- first component corresponds to 0 = sxs, 1 = sxt, 2 = txt
      msq_cs = 0._dp
C --- non-zero components
c      msq_cs(1,2,-2) = gwsq*mst*propw
c      msq_cs(1,2,-2) = mst*propw*0.222_dp**4 ! removed gwsq, added Vsq(c,d)
c      msq_cs(1,1,-1) = (vi(1)**2 + ai(1)**2)*(mst*propz_t + mts*propz_s)
c     .     + Q(1)**2*(mst + mts)
c      msq_cs(1,4,-4) = msq_cs(1,2,-2)
c      msq_cs(2,1,-1) = (vi(1)**2*mtt(1) + ai(1)**2*mtt(2))*propz_t
c     .     + Q(1)**2*mss(1)


      msq_cs(1,1,-1) = (vi(1)**2 + ai(1)**2)*mst*propz_s + Qf(1)**2*mst
      msq_cs(2,2,-2) = cwsq*mts*propw*Vud**2
      msq_cs(2,4,-4) = cwsq*mts*propw*Vcd**2
      msq_cs(2,1,-1) = (vi(1)**2 + ai(1)**2)*mts*propz_t + Qf(1)**2*mts
      msq_cs(3,1,-1) = (vi(1)**2*mtt(1) + ai(1)**2*mtt(2))*propz_t
     .     + Qf(1)**2*mtt(1)

      do j = 1,nf-1
         k = -j
         msq_cs(0,j,k) = 
     .        + (vi(j)*vf(1)*mss(1) + ai(j)*af(1)*mss(2))*propz_s
     .        + Qf(j)*Qf(1)*mss(1)
      end do

c      print*, 'v,a', vi(1:2), ai(1:2)

C--- msq_cs has order of alpha_s*alpha
      msq_cs = 2._dp*msq_cs*gsq*esq

      end subroutine qqb_twojet_ddb_mix





      subroutine qqb_twojet_ccb_mix(msq_cs,ss,tt,uu,flg)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      real(dp):: msq_cs(0:3,1:nf,fn:-1),ss,tt,uu,
     .     vi(nf),vf(nf),ai(nf),af(nf),Qf(nf),cwsq,mz,mw,
     .     propz_s,propz_t,propw,mss(2),mst,mts,mtt(2)
      integer j,k
      logical flg
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      mz = zmass
      mw = wmass

C --- cwsq is W coupling to fermion
      cwsq = 1._dp/4._dp/xw


      propz_s = ss*(ss-mz**2)/((ss-mz**2)**2 + zwidth**2*mz**2)
      propz_t = tt*(tt-mz**2)/((tt-mz**2)**2 + zwidth**2*mz**2)
      propw = tt*(tt-mw**2)/((tt-mw**2)**2 + wwidth**2*mw**2)

      vi(1:nf) = (l(1:nf) + r(1:nf))/2._dp
      vf(1:nf) = vi(1:nf)
      ai(1:nf) = (l(1:nf) - r(1:nf))/2._dp
      af(1:nf) = ai(1:nf)


C --- flag determines whether photon mediated amplitude is included

      Qf(1:nf) = 0._dp

      if(flg) Qf(1:nf) = Q(1:nf)

C --- color structure cs1, cs2, cs3
C --- cs1 = (Nc^2 - 1)/4; cs2 = -(Nc^2 -1)/(4*Nc); cs3 = (Nc^2 - 1)^2/(4*Nc)

      call qqb_twojet_born_mix(mss,mst,mts,mtt,ss,tt,uu)

C --- msq_cs(0:2,:,:)
C --- first component corresponds to 0 = sxs, 1 = sxt, 2 = txt
      msq_cs = 0._dp
C --- non-zero components
c      msq_cs(1,1,-1) = gwsq*mst*propw
c      msq_cs(1,4,-4) = (vi(4)**2 + ai(4)**2)*(mst*propz_t + mts*propz_s)
c     .     + Q(4)**2*(mst + mts)
c      msq_cs(1,3,-3) = msq_cs(1,1,-1)
c      msq_cs(2,4,-4) = (vi(4)**2*mtt(1) + ai(4)**2*mtt(2))*propz_t
c     .     + Q(4)**2*mss(1)


      msq_cs(1,4,-4) = (vi(4)**2 + ai(4)**2)*mst*propz_s + Qf(4)**2*mst
      msq_cs(2,1,-1) = cwsq*mts*propw*Vcd**2
      msq_cs(2,3,-3) = cwsq*mts*propw*Vcs**2
!      print*, msq_cs(2,3,-3),mts,propw,cwsq,vcs,mss,mst
!      stop
      msq_cs(2,4,-4) = (vi(4)**2 + ai(4)**2)*mts*propz_t + Qf(4)**2*mts
      msq_cs(3,4,-4) = (vi(4)**2*mtt(1) + ai(4)**2*mtt(2))*propz_t
     .     + Qf(4)**2*mtt(1)

      do j = 1,nf-1
         k = -j
         msq_cs(0,j,k) = 
     .        + (vi(j)*vf(4)*mss(1) + ai(j)*af(4)*mss(2))*propz_s
     .        + Qf(j)*Qf(4)*mss(1)
      end do

C --- msq_cs has order of alpha_s*alpha
      msq_cs = 2._dp*msq_cs*gsq*esq

      end subroutine qqb_twojet_ccb_mix





      subroutine qqb_twojet_ssb_mix(msq_cs,ss,tt,uu,flg)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      real(dp):: msq_cs(0:3,1:nf,fn:-1),ss,tt,uu,
     .     vi(nf),vf(nf),ai(nf),af(nf),Qf(nf),cwsq,mz,mw,
     .     propz_s,propz_t,propw,mss(2),mst,mts,mtt(2)
      integer j,k
      logical flg
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb

      mz = zmass
      mw = wmass

C --- cwsq is W coupling to fermion
      cwsq = 1._dp/4._dp/xw


      propz_s = ss*(ss-mz**2)/((ss-mz**2)**2 + zwidth**2*mz**2)
      propz_t = tt*(tt-mz**2)/((tt-mz**2)**2 + zwidth**2*mz**2)
      propw = tt*(tt-mw**2)/((tt-mw**2)**2 + wwidth**2*mw**2)

      vi(1:nf) = (l(1:nf) + r(1:nf))/2._dp
      vf(1:nf) = vi(1:nf)
      ai(1:nf) = (l(1:nf) - r(1:nf))/2._dp
      af(1:nf) = ai(1:nf)

C --- flag determines whether photon mediated amplitude is included

      Qf(1:nf) = 0._dp

      if(flg) Qf(1:nf) = Q(1:nf)


C --- color structure cs1, cs2, cs3
C --- cs1 = (Nc^2 - 1)/4; cs2 = -(Nc^2 -1)/(4*Nc); cs3 = (Nc^2 - 1)^2/(4*Nc)

      call qqb_twojet_born_mix(mss,mst,mts,mtt,ss,tt,uu)

C --- msq_cs(0:2,:,:)
C --- first component corresponds to 0 = sxs, 1 = sxt, 2 = txt
      msq_cs = 0._dp
C --- non-zero components
c      msq_cs(1,2,-2) = gwsq*mst*propw
c      msq_cs(1,3,-3) = (vi(3)**2 + ai(3)**2)*(mst*propz_t + mts*propz_s)
c     .     + Q(3)**2*(mst + mts)
c      msq_cs(1,4,-4) = msq_cs(1,2,-2)
c      msq_cs(2,3,-3) = (vi(3)**2*mtt(1) + ai(3)**2*mtt(2))*propz_t
c     .     + Q(3)**2*mss(1)


      msq_cs(1,3,-3) = (vi(3)**2 + ai(3)**2)*mst*propz_s + Qf(3)**2*mst
      msq_cs(2,2,-2) = cwsq*mts*propw*Vus**2
      msq_cs(2,4,-4) = cwsq*mts*propw*Vcs**2
      msq_cs(2,3,-3) = (vi(3)**2 + ai(3)**2)*mts*propz_t + Qf(3)**2*mts
      msq_cs(3,3,-3) = (vi(3)**2*mtt(1) + ai(3)**2*mtt(2))*propz_t
     .     + Qf(3)**2*mtt(1)


      do j = 1,nf-1
         k = -j
         msq_cs(0,j,k) = 
     .        + (vi(j)*vf(3)*mss(1) + ai(j)*af(3)*mss(2))*propz_s
     .        + Qf(j)*Qf(3)*mss(1)
      end do


C --- msq_cs has order of alpha_s*alpha
      msq_cs = 2._dp*msq_cs*gsq*esq


      end subroutine qqb_twojet_ssb_mix





      subroutine qqb_twojet_born_mix(mss,mst,mts,mtt,ss,tt,uu)
      implicit none
      include 'types.f'
      real(dp):: mss(2),mst,mts,mtt(2),ss,tt,uu

C --- the vector and axial vector part in sxs are different, 
C --- '1' denotes vector portion and '2' axial
      mss(1) = + 8._dp*(tt**2 + uu**2)
      mss(2) = - 8._dp*(tt**2 - uu**2)
      mss = mss/ss**2

C --- mtt = mss w/ ss <-> tt
      mtt(1) = + 8._dp*(ss**2 + uu**2)
      mtt(2) = - 8._dp*(ss**2 - uu**2)
      mtt = mtt/tt**2

C --- the vector and axial portion in sxt are the same
c      mst = - 8._dp*uu**2
c      mts = - 8._dp*uu**2
C --- add '-' due to single fermion-loop in s- and t-channel interference
      mst = + 8._dp*uu**2
      mts = + 8._dp*uu**2
      mst = mst/ss/tt
      mts = mts/tt/ss

      end subroutine qqb_twojet_born_mix

      


      

