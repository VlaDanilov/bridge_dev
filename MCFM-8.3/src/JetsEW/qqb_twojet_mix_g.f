      subroutine qqb_twojet_mix_g(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),
     .     s12,s34,t1,t2,u1,u2,s15,s25,s35,s45,
     .     qaii_jj(nf,nf),qqij_ij(nf,nf),aqii_jj(nf,nf),aaij_ij(nf,nf),
     .     aqij_ij(nf,nf),qaij_ij(nf,nf),qg(nf,nf),gq(nf,nf),ga(nf,nf),
     .     ag(nf,nf),old,new
      integer j,k
      real(dp):: pt,aetarap,ayrap
      logical flg/.false./

      call dotem(5,p,s)

      s12 = s(1,2)
      t1 = s(1,3)
      t2 = s(2,4)
      u1 = s(1,4)
      u2 = s(2,3)
      s34 = s(3,4)

      s15 = s(1,5)
      s25 = s(2,5)
      s35 = s(3,5)
      s45 = s(4,5)

c      s34 = -(s12 + t1 + t2 + u1 + u2)

      msq = 0._dp
C --- quark induced processes
      call qaii_jj_mix(qaii_jj,s12,t1,t2,u1,u2,flg)
      call qaii_jj_mix(qaij_ij,t1,s12,s34,u1,u2,flg)
      call qaii_jj_mix(qqij_ij,t1,u1,u2,s12,s34,flg)

c      call qaii_jj_mix(aqii_jj,s12,t2,t1,u2,u1,flg)
c      call qaii_jj_mix(aqij_ij,t1,s34,s12,u2,u1,flg)
c      call qaii_jj_mix(aaij_ij,t1,u2,u1,s34,s12,flg)
c--- Use symmetry relations instead of calls, for speed
      aqii_jj(:,:)=qaii_jj(:,:)
      aqij_ij(:,:)=qaij_ij(:,:)
      aaij_ij(:,:)=qqij_ij(:,:)     

C --- gluon induced processes
c      call qaii_jj_mix(ga,s35,u2,s45,s34,s25,flg) ! gq~ -> r~q~ r
      call qaii_jj_mix(gq,s35,s45,u2,s25,s34,flg) ! gq  -> rq r~

      call qaii_jj_mix(qg,t1,u1,s35,s15,s34,flg) ! qg  -> qr r~
c      call qaii_jj_mix(ag,t1,s35,u1,s34,s15,flg) ! q~g  -> q~r~ r
c--- flip sign to account for crossing
c      ga(:,:)=-ga(:,:)
      gq(:,:)=-gq(:,:)
      qg(:,:)=-qg(:,:)
c      ag(:,:)=-ag(:,:)

c--- Use symmetry relations instead of calls, for speed
      ga(:,:)=gq(:,:)
      ag(:,:)=qg(:,:)
      
      do j = -4,4
         do k = -4,4

C --- initial states qa

            if((j > 0) .and. (k < 0)) then
               if(j == -k) then
                  msq(j,k) =
     &             aveqq*(qaii_jj(j,1) + qaii_jj(j,2) 
     &                  + qaii_jj(j,3) + qaii_jj(j,4))
               else
                  msq(j,k) = aveqq*qaij_ij(j,-k)
               end if

C --- qq
            else if((j > 0) .and. (k > 0)) then
               if(j == k) then
                  msq(j,k) = aveqq*qqij_ij(j,k)*half
               else
                  msq(j,k) = aveqq*qqij_ij(j,k)
               end if


C --- aa
            else if((j < 0) .and. (k < 0)) then
               if(j == k) then
                  msq(j,k) = aveqq*aaij_ij(-j,-k)*half
               else
                  msq(j,k) = aveqq*aaij_ij(-j,-k)
               end if

C --- aq
            else if((j < 0) .and. (k > 0)) then
               if(k == -j) then
                  msq(j,k) = 
     &             aveqq*(aqii_jj(k,1) + aqii_jj(k,2)
     &                  + aqii_jj(k,3) + aqii_jj(k,4))
               else
                  msq(j,k) = aveqq*aqij_ij(-j,k)
               end if

C --- qg
            else if((j > 0) .and. (k == 0)) then
               msq(j,k) = aveqg*(qg(j,1)+qg(j,2)+qg(j,3)+qg(j,4)
     &                          -qg(j,j)+qg(j,j)*half) ! account for identical q

C --- ag
            else if((j < 0) .and. (k == 0)) then
               msq(j,k) = aveqg*(ag(-j,1)+ag(-j,2)+ag(-j,3)+ag(-j,4)
     &                          -ag(-j,-j)+ag(-j,-j)*half) ! account for identical q

C --- gq
            else if((j == 0) .and. (k > 0)) then
               msq(j,k) = aveqg*(gq(k,1)+gq(k,2)+gq(k,3)+gq(k,4)
     &                          -gq(k,k)+gq(k,k)*half) ! account for identical q

C --- ga
            else if((j == 0) .and. (k < 0)) then
               msq(j,k) = aveqg*(ga(-k,1)+ga(-k,2)+ga(-k,3)+ga(-k,4)
     &                          -ga(-k,-k)+ga(-k,-k)*half) ! account for identical q

            end if

         end do
      end do

c--- apply overall factor
      msq = msq*gsq**2*esq

c--- Cross-check with JC calculation of real matrix elements
c      old=qaii_jj(2,4)*gsq**2*esq
c      write(6,*) 'Old uu~ -> cc~',old
c      call realmix_uubccb(1,2,3,4,5,new)
c      write(6,*) 'New uu~ -> cc~',new,new/old
c      write(6,*)
c      old=qaii_jj(2,2)*gsq**2*esq
c      write(6,*) 'Old uu~ -> uu~',old
c      call realmix_uubuub(1,2,3,4,5,new)
c      write(6,*) 'New uu~ -> uu~',new,new/old
c      write(6,*)      
c      old=qaii_jj(2,1)*gsq**2*esq
c      write(6,*) 'Old uu~ -> dd~',old
c      call realmix_uubddb(1,2,3,4,5,new)
c      write(6,*) 'New uu~ -> dd~',new,new/old
c      write(6,*)      
c      pause


! DEBUGGING CODE
c      msq(:,:) = 0._dp
c      msq(2,2) = qqij_ij(2,2)*half
c      msq(2,1) = qqij_ij(2,1)
c      msq(1,2) = qqij_ij(1,2)

      if(.false.) then
        print*, "phase point: "
        print*, "p1: ", p(1,:)
        print*, "p2: ", p(2,:)
        print*, "p3: ", p(3,:)
        print*, "p4: ", p(4,:)
        print*, "p5: ", p(5,:)
        write(6,*) "s12,s34,t1,t2,u1,u2: ", s12, s34, t1, t2, u1, u2 
        write(6,*) "ud -> ud: ", qqij_ij(2,1)*(4._dp*pi*0.118_dp)**2*esq
        pause
      end if

C --- for cuts on emitted gluon
c      if ((p(5,4) < 2._dp) .or. (sqrt(dabs(s(3,5))) < 5._dp) .or. 
c     .     (sqrt(dabs(s(4,5))) < 5._dp) .or. (sqrt(dabs(s(1,5))) < 5._dp) 
c     .     .or. (sqrt(dabs(s(2,5))) < 5._dp)) then
c         msq(2,1) = 0._dp
c      else
c         msq(2,1) = qqij_ij(2,1)
c      end if

c      if((pt(3,p) > 1d2) .and. (pt(4,p) > 1d2) .and. 
c     .     (aetarap(3,p) < 0.5_dp) .and. (aetarap(4,p) < 0.5_dp)) then
c     .     (ayrap(3,p) < 2.5_dp) .and. (ayrap(4,p) < 2.5_dp)) then
c         msq(2,1) = qqij_ij(2,1)
c      endif
      
c      msq = 0._dp
c      msq = msq*esq*(4._dp*pi*0.118_dp)**2

c      print*, "p(5,4): ",p(5,4),s(1,5),s(2,5),s(3,5),s(4,5)
c      print*, msq(2,1)
c      stop


      end subroutine qqb_twojet_mix_g



      subroutine qaii_jj_mix(msq,s12,t1,t2,u1,u2,flg)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'ewcharge.f'
      include 'zcouple.f'
      real(dp):: msq(nf,nf),s12,t1,t2,u1,u2,cc(nf,2),
     .     vsq(nf,nf),Qf(nf)
      integer j,k
      logical flg
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

      cc(1:nf,1) = (l(1:nf) + r(1:nf))/2._dp 
      cc(1:nf,2) = (l(1:nf) - r(1:nf))/2._dp 

C --- flag determines whether photon mediated amplitude is included
c      flg = .false.

      Qf(1:nf) = 0._dp

      if(flg) Qf(1:nf) = Q(1:nf)

      msq = 0._dp
      do j = 1,nf-1
         do k = j,nf-1
            if(j == k) then
               call qaii_iig(msq(j,k),s12,t1,t2,u1,u2,cc(j,:),Qf(j))
            else if(abs(j-k) == 2) then
               call qaii_jjg1(msq(j,k),s12,t1,t2,u1,u2,cc(j,:),cc(k,:),
     .              Qf(j),Qf(k))
            else
               call qaii_jjg2(msq(j,k),s12,t1,t2,u1,u2,cc(j,:),cc(k,:),
     .              Qf(j),Qf(k),vsq(j,k))
            endif
            msq(k,j)=msq(j,k)
         end do
      end do
         
      end subroutine qaii_jj_mix




C --- for u_i ub_i -> u_j + ub_j + g or d_i db_i -> d_j db_j + g
      subroutine qaii_jjg1(msq,s12,t1,t2,u1,u2,cci,ccf,ccqi,ccqf)
C --- ccz is Z coupling to fermion
C --- ccq is electric charge
      implicit none
      include 'types.f'
      include 'masses.f'
      real(dp):: msq,mss(2),s12,s34,t1,t2,u1,u2,cci(2),ccf(2),
     .     ccqi,ccqf,cqi(2),cqf(2),mz(2),mph(2)

      mz(1) = zmass
      mz(2) = mz(1)*zwidth
      mph = 0._dp

C --- if we include decay width in Z/W and neglect the part supressed by width
c      propz = (s12 - mz**2)**2/((s12 - mz**2)**2 + zwidth**2*mz**2)

      cqi = (/ccqi, 0._dp/)
      cqf = (/ccqf, 0._dp/)

      s34 = -(s12 + t1 + t2 + u1 + u2)
      
      call qaii_jjg_sxs(mss(1),s12,t1,t2,u1,u2,mz,cci,ccf)
      call qaii_jjg_sxs(mss(2),s12,t1,t2,u1,u2,mph,cqi,cqf)

      msq = mss(1) + mss(2)

      end subroutine qaii_jjg1



C --- u_i ub_i -> d_j db_j or d_i db_i -> u_j ub_j
      subroutine qaii_jjg2(msq,s12,t1,t2,u1,u2,cci,ccf,ccqi,ccqf,vsqij)
C --- vsqij is CKM matrix element square
      implicit none
      include 'types.f'
      include 'masses.f'
      include 'ewcouple.f'
      real(dp):: msq,s12,t1,t2,u1,u2,s34,cci(2),ccf(2),ccqi,ccqf,
     .     vsqij,cwi(2),cwf(2),mss,mts,mw(2),ccw

      mw(1) = wmass
      mw(2) = mw(1)*wwidth

      ccw = 1._dp/2._dp/sqrt(2._dp*xw)
      cwi = (/ccw,ccw/)
      cwf = cwi


      s34 = -(s12 + t1 + t2 + u1 + u2)

      call qaii_jjg1(mss,s12,t1,t2,u1,u2,cci,ccf,ccqi,ccqf)
c      mss=0._dp ! DEBUG
      call qaii_jjg_sxt(mts,t1,s12,s34,u1,u2,mw,cwi,cwf)
c      mts=0._dp ! DEBUG
c--- DEBUG: Jia check
c      mts=0._dp

      msq = mss + mts*vsqij

      end subroutine qaii_jjg2



C --- u_i ub_i -> u_i ub_i or d_i db_i -> d_i db_i
      subroutine qaii_iig(msq,s12,t1,t2,u1,u2,ccz,ccq)
      implicit none
      include 'types.f' 
      include 'masses.f'
      real(dp):: msq,s12,s34,t1,t2,u1,u2,ccz(2),cci(2),ccf(2),ccq,
     .     cqi(2),cqf(2),mz(2),mph(2),mss(2),mtt(2),mst(2),mts(2)

      mz(1) = zmass
      mz(2) = mz(1)*zwidth
C --- photon mass and width
      mph = 0._dp

      cci = ccz
      ccf = cci
      cqi = (/ccq, 0._dp/)
      cqf = cqi

      s34 = -(s12 + t1 + t2 + u1 + u2)

      call qaii_jjg_sxs(mss(1),s12,t1,t2,u1,u2,mz,cci,ccf)
      call qaii_jjg_sxs(mss(2),s12,t1,t2,u1,u2,mph,cqi,cqf)
      call qaii_jjg_sxs(mtt(1),t1,s12,s34,u1,u2,mz,cci,ccf)
      call qaii_jjg_sxs(mtt(2),t1,s12,s34,u1,u2,mph,cqi,cqf)
      call qaii_jjg_sxt(mst(1),s12,t1,t2,u1,u2,mz,cci,ccf)
      call qaii_jjg_sxt(mst(2),s12,t1,t2,u1,u2,mph,cqi,cqf)
      call qaii_jjg_sxt(mts(1),t1,s12,s34,u1,u2,mz,cci,ccf)
      call qaii_jjg_sxt(mts(2),t1,s12,s34,u1,u2,mph,cqi,cqf)

      msq = 
     .     + mss(1) + mss(2) 
     .     + mtt(1) + mtt(2) 
     .     + mst(1) + mst(2) 
     .     + mts(1) + mts(2)

C --- DEBUG
c      msq = mtt(1) + mtt(2) 
c      msq = mst(1) + mst(2) 
c      msq = mss(1) + mss(2)
c      msq = mts(1) + mts(2)

      end subroutine qaii_iig
      


      subroutine qaii_jjg_sxs(msq,s12,t1,t2,u1,u2,mv,cci,ccf)
      implicit none
      include 'types.f'
      real(dp):: msq,s12,t1,t2,u1,u2,mv(2),cci(2),ccf(2),
     .     vi,ai,vf,af,s34,sv12,sv34,genfac,p1dp5,p2dp5,p3dp5,p4dp5,
     .     prop12,prop34,mss12,mss34

C --- complex mass mv = (mv(1), mv(2)) w/ mv(2) = Gamma_v*mv(1)
      s34 = -(s12 + t1 + t2 + u1 + u2)
      sv12 = s12 - mv(1)**2
      sv34 = s34 - mv(1)**2
C --- ratio of propagator w/o width to that w/ width (real)
      prop12 = sv12**2/(sv12**2 + mv(2)**2)
      prop34 = sv34**2/(sv34**2 + mv(2)**2)
C --- removed; since there is no massive t-channel propagator needed
c      tv13 = t1 - mv**2
c      tv24 = t2 - mv**2
      p1dp5 = -(s12 + t1 + u1)/2._dp
      p2dp5 = -(s12 + t2 + u2)/2._dp
      p3dp5 = +(s12 + t2 + u1)/2._dp
      p4dp5 = +(s12 + t1 + u2)/2._dp

c      genfac = 2._dp/p1dp5/p2dp5/p3dp5/p4dp5/s12/s34/sv12/sv34
c     .     /t1/t2/tv13/tv24

      genfac = 2._dp/p1dp5/p2dp5/p3dp5/p4dp5/s12/s34/sv12/sv34

      vi = cci(1)
      vf = ccf(1)
      ai = cci(2)
      af = ccf(2)

      mss12 = 
     .     + s12*sv34*(s12**2*(t1 + t2 - u1 - u2) 
     .     + s12*(t1 + t2 - u1 - u2)*(t1 + t2 + u1 + u2) + (t1 + t2 + u1 
     .     + u2)*(t1*t2 - u1*u2))*(af*ai*(t1**2 + t2**2 - u1**2 - u2**2) 
     .     - (t1**2 + t2**2 + u1**2 + u2**2)*vf*vi)
      mss12 = mss12*prop12

      mss34 = 
     .     + sv12*s34*(s12**2*(t1 + t2 - u1 - u2)
     .     + s12*(t1 + t2 - u1 - u2)*(t1 + t2 + u1 + u2) + (t1 + t2 + u1 
     .     + u2)*(t1*t2 - u1*u2))*(af*ai*(t1**2 + t2**2 - u1**2 - u2**2)
     .     - (t1**2 + t2**2 + u1**2 + u2**2)*vf*vi)
      mss34 = mss34*prop34

      msq = genfac*(mss12 + mss34)


      end subroutine qaii_jjg_sxs




      subroutine qaii_jjg_sxt(msq,s12,t1,t2,u1,u2,mv,cci,ccf)
      implicit none
      include 'types.f'
      real(dp):: msq,s12,t1,t2,u1,u2,mv(2),cci(2),ccf(2),
     .     vi,ai,vf,af,s34,sv12,sv34,genfac,p1dp5,p2dp5,p3dp5,p4dp5,
     .     prop12,prop34,mst12,mst34

      s34 = -(s12 + t1 + t2 + u1 + u2)
      sv12 = s12 - mv(1)**2
      sv34 = s34 - mv(1)**2
C --- removed
c      tv13 = t1 - mv**2
c      tv24 = t2 - mv**2

      prop12 = sv12**2/(sv12**2 + mv(2)**2)
      prop34 = sv34**2/(sv34**2 + mv(2)**2)


      p1dp5 = -(s12 + t1 + u1)/2._dp
      p2dp5 = -(s12 + t2 + u2)/2._dp
      p3dp5 = +(s12 + t2 + u1)/2._dp
      p4dp5 = +(s12 + t1 + u2)/2._dp

c      genfac = 2._dp/p1dp5/p2dp5/p3dp5/p4dp5/s12/s34/sv12/sv34
c     .     /t1/t2/tv13/tv24

      genfac = 2._dp/p1dp5/p2dp5/p3dp5/p4dp5/sv12/sv34/t1/t2

      vi = cci(1)
      vf = ccf(1)
      ai = cci(2)
      af = ccf(2)



      mst12 = 
     .     - 1._dp/3._dp*sv34*(u1**2 + u2**2)*(9._dp*s12**4 - (t1*t2 
     .     - u1*u2)*(7._dp*t1*t2 + 8._dp*t2*u1 + 8._dp*t1*u2 
     .     + 9._dp*u1*u2) + s12*(t2*(t1**2 + 8._dp*u1*(t2 + u1) 
     .     + t1*(t2 + 8._dp*u1)) + (8._dp*t1*(t1 + t2) 
     .     + 25._dp*(t1 + t2)*u1 + 18._dp*u1**2)*u2 + 2._dp*(4._dp*t1 
     .     + 9._dp*u1)*u2**2) + s12**3*(17._dp*t1 + 17._dp*t2 + 18._dp*(u1 
     .     + u2)) + s12**2*(8._dp*t1**2 + 8._dp*t2**2 + t2*(25._dp*u1 
     .     + 17._dp*u2) + t1*(18._dp*t2 + 17._dp*u1 + 25._dp*u2) 
     .     + 9._dp*(u1**2 + 4._dp*u1*u2 + u2**2)))*(af*ai + vf*vi)

      mst12 = mst12*prop12


      mst34 = 
     .     - 1._dp/3._dp*sv12*(u1**2 + u2**2)*(9._dp*s12**4 + t1*t2*(t1**2 
     .     - 5._dp*t1*(t2 + u1) + (t2 + u1)*(t2 + 2._dp*u1)) 
     .     + (-t1*(5._dp*t2 - 9._dp*u1)*(t2 + u1) + t2*u1*(t2 + u1) 
     .     + t1**2*(3._dp*t2 + u1))*u2 + (9._dp*u1*(t2 + u1) + t1*(2._dp*t2 
     .     + u1))*u2**2 + s12**3*(19._dp*t1 + 19._dp*t2 + 18._dp*(u1 + u2)) 
     .     + s12**2*(11._dp*t1**2 + 11._dp*t2**2 + 4._dp*t1*(6._dp*t2 
     .     + 7._dp*u1 + 5._dp*u2) + 4._dp*t2*(5._dp*u1 + 7._dp*u2) 
     .     + 9._dp*(u1**2 + 4._dp*u1*u2 + u2**2)) + s12*(t1**3 
     .     + t2*(t2 + u1)**2 + (t2 + 2._dp*u1)*(10._dp*t2 + 9._dp*u1)*u2 
     .     + 9._dp*(t2 + 2._dp*u1)*u2**2 + 2._dp*t1**2*(3._dp*t2 + 5._dp*u1 
     .     + u2) + t1*(6._dp*t2**2 + 9._dp*u1**2 + 29._dp*u1*u2 + u2**2 
     .     + 16._dp*t2*(u1 + u2))))*(af*ai + vf*vi)

      mst34 = mst34*prop34


C --- for the contribution missing trial-gluon coupling
      if(.false.) then
      mst12 = 
     .     + 1._dp/3._dp*sv34*(9._dp*s12**4*(t1 - t2)*(u1 - u2) 
     .     + s12*(t2*u1*(18._dp*t1*(t1 - t2)*(t1 + t2) + (8._dp*t1**2 
     .     - 19._dp*t1*t2 - 9._dp*t2**2)*u1 - 8._dp*(t1 + t2)*u1**2 
     .     + u1**3) - (18._dp*t1*(t1 - t2)*t2*(t1 + t2) - 9._dp*(t1 
     .     - t2)**2*(t1 + t2)*u1 + 8._dp*t1*(t1 + t2)*u1**2 
     .     + 2._dp*(8._dp*t1 - t2)*u1**3)*u2 + (t1*(-9._dp*t1**2 
     .     - 19._dp*t1*t2 + 8._dp*t2**2) - 8._dp*t2*(t1 + t2)*u1 + (t1 
     .     + t2)*u1**2)*u2**2 - 2._dp*(4._dp*t1*(t1 + t2) - (t1 
     .     - 8._dp*t2)*u1)*u2**3 + t1*u2**4) + s12**3*(18._dp*t1**2*(u1 
     .     - u2) + t1*(u1**2 - 17._dp*u2**2) + t2*(-17._dp*u1**2 + u2**2 
     .     + 18._dp*t2*(-u1 + u2))) + s12**2*(9._dp*t1**3*(u1 - u2) 
     .     + t1**2*(u1**2 + 27._dp*t2*(u1 - u2) + 9._dp*u1*u2 
     .     - 26._dp*u2**2) - t1*(27._dp*t2**2*(u1 - u2) + 9._dp*t2*(u1 
     .     + u2)**2 + (8._dp*u1 + 7._dp*u2)*(u1**2 + u2**2)) 
     .     + t2*(9._dp*t2**2*(-u1 + u2) - (7._dp*u1 + 8._dp*u2)*(u1**2 
     .     + u2**2) + t2*(-26._dp*u1**2 + 9._dp*u1*u2 + u2**2))) + (t1*t2 
     .     - u1*u2)*(9._dp*t1**2*(u1 - u2)*(t2 + u2) - t2*u1*(u1**2 
     .     + 9._dp*t2*(u1 - u2) + 9._dp*u1*u2 - 8._dp*u2**2) 
     .     - t1*(9._dp*t2**2*(u1 - u2) + t2*(-7._dp*u1**2 + 18._dp*u1*u2 
     .     - 7._dp*u2**2) + u2*(-8._dp*u1**2 + 9._dp*u1*u2 
     .     + u2**2))))*(af*ai + vf*vi)


      mst12 = mst12*prop12


      mst34 = 
     .     - 1._dp/3._dp*sv12*(9._dp*s12**4*(t1 - t2)*(u1 - u2) - t2*u1*(t2 
     .     + u1)*u2*(18._dp*t2*u1 + 17._dp*u1**2 + 9._dp*u1*u2 
     .     + 8._dp*u2**2) + t1**3*(t2*u1*(9._dp*t2 + 10._dp*u1) 
     .     - 9._dp*t2*(t2 + u1)*u2 - (17._dp*t2 + 18._dp*u1)*u2**2) 
     .     + s12**2*(u1*(9._dp*t1**3 - t1*(9._dp*t2 - 2._dp*u1)*(3._dp*t2 
     .     + 5._dp*u1) + t1**2*(27._dp*t2 + 29._dp*u1) - t2*(t2 
     .     + u1)*(9._dp*t2 + 43._dp*u1)) + (-9._dp*(t1 - t2)*(t1**2 
     .     + 4._dp*t1*t2 + t2**2) - 9._dp*(t1 - t2)**2*u1 + 11._dp*(t1 
     .     - 4._dp*t2)*u1**2)*u2 + (-52._dp*t1**2 - 39._dp*t1*t2 
     .     + 29._dp*t2**2 - 44._dp*t1*u1 + 11._dp*t2*u1)*u2**2 + (-43._dp*t1 
     .     + 10._dp*t2)*u2**3) - t1**2*(t2*u1*(t2 + u1)*(9._dp*t2 
     .     + 14._dp*u1) - (9._dp*t2**3 + 18._dp*t2**2*u1 + 3._dp*t2*u1**2 
     .     - 8._dp*u1**3)*u2 + (t2 + u1)*(23._dp*t2 + 9._dp*u1)*u2**2 
     .     + (33._dp*t2 + 35._dp*u1)*u2**3) + s12*(u1*(-17._dp*t2*u1*(t2 
     .     + u1)**2 + 2._dp*t1**2*u1*(3._dp*t2 + 5._dp*u1) 
     .     + 2._dp*t1**3*(9._dp*t2 + 5._dp*u1) - t1*t2*(18._dp*t2**2 
     .     + 75._dp*t2*u1 + 56._dp*u1**2)) + (18._dp*t1*t2*(-t1**2 + t2**2) 
     .     - 9._dp*(t1 - t2)**2*(t1 + t2)*u1 + 2._dp*(t1**2 - 10._dp*t1*t2 
     .     - 31._dp*t2**2)*u1**2 + 2._dp*(t1 - 26._dp*t2)*u1**3)*u2 
     .     - (17._dp*t1**3 + 75._dp*t1**2*t2 - 6._dp*t1*t2**2 - 10._dp*t2**3 
     .     + 2._dp*(31._dp*t1**2 + 10._dp*t1*t2 - t2**2)*u1 + 17._dp*(t1 
     .     + t2)*u1**2)*u2**2 + 2._dp*(-17._dp*t1**2 + t2*(5._dp*t2 + u1) 
     .     - 2._dp*t1*(14._dp*t2 + 13._dp*u1))*u2**3 - 17._dp*t1*u2**4) 
     .     + t1*(-t2*u1**2*(t2 + u1)*(17._dp*t2 + 16._dp*u1) - t2*u1*(t2 
     .     + u1)*(9._dp*t2 + 23._dp*u1)*u2 + (5._dp*t2 
     .     + 4._dp*u1)*(2._dp*t2**2 - t2*u1 - 2._dp*u1**2)*u2**2 - (t2 
     .     + u1)*(14._dp*t2 + 9._dp*u1)*u2**3 - (16._dp*t2 
     .     + 17._dp*u1)*u2**4) + s12**3*(18._dp*t1**2*(u1 - u2) 
     .     + t1*(19._dp*u1**2 - 35._dp*u2**2) + t2*(-35._dp*u1**2 
     .     + 19._dp*u2**2 + 18._dp*t2*(-u1 + u2))))*(af*ai + vf*vi)


      mst34 = mst34*prop34
      end if

C --- add '-' due to single fermion-loop in s- and t-channel interference
      msq = - genfac*(mst12 + mst34)

C --- Debug

c      mst12 = 0._dp
c      mst34 = 
c     .   + 1._dp/3._dp*sv12*t1*(s12 + t1 + u2)*(s12 + t2 + u2)*(s12**2*(u1 
c     .   + u2) + u1*u2*(u1 + u2) + s12*(u1 + u2)*(t1 + t2 + u1 + u2) 
c     .   + t1*(2._dp*u1*u2 + t2*(u1 + u2)))*(af*ai + vf*vi)

c      msq = 0._dp

      end subroutine qaii_jjg_sxt
