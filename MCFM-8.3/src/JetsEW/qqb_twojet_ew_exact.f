      subroutine qqb_twojet_ew_exact(p,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'nflav.f'
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf),ss,tt,uu,
     &     msq_qa_iijj(nf,nf),msq_qa_ijij(nf,nf),
     &     msq_qq_ijij(nf,nf),msq_aq_iijj(nf,nf),
     &     msq_aq_ijij(nf,nf),msq_aa_ijij(nf,nf),
     &     msq_gg_qa(nf),msq_gq_gq(nf),msq_qg_qg(nf),
     &     msq_qa_gg(nf),msq_aq_gg(nf),msq_ga_ga(nf),
     &     msq_ag_ag(nf),
C --- variables for the LO mixed and pure weak Born
     &     qa_iijj(nf,nf),qa_ijij(nf,nf),
     &     qq_ijij(nf,nf),aq_iijj(nf,nf),
     &     aq_ijij(nf,nf),aa_ijij(nf,nf),
     &     qa_iijj2(nf,nf),qa_ijij2(nf,nf),
     &     qq_ijij2(nf,nf),aq_iijj2(nf,nf),
     &     aq_ijij2(nf,nf),aa_ijij2(nf,nf),
     &     p_qa_iijj(mxpart,4),p_qa_ijij(mxpart,4),
     &     p_qq_ijij(mxpart,4),p_aq_iijj(mxpart,4),
     &     p_aq_ijij(mxpart,4),p_aa_ijij(mxpart,4)
      integer j,k
      logical, parameter:: ew_tree=.false.

      scheme = 'dred'

      call dotem(4,p,s)
      ss = s(1,2)
      tt = s(1,3)
      uu = s(2,3)

      call dijet_qqb_ii_jj_v(msq_qa_iijj,ss,tt,uu)
      call dijet_qqb_ii_jj_v(msq_qa_ijij,tt,ss,uu)
      call dijet_qqb_ii_jj_v(msq_qq_ijij,tt,uu,ss)
c      call dijet_qqb_ii_jj_v(msq_aq_iijj,ss,tt,uu)
c      call dijet_qqb_ii_jj_v(msq_aq_ijij,tt,ss,uu)
c      call dijet_qqb_ii_jj_v(msq_aa_ijij,tt,uu,ss)
c--- Use symmetry relations instead of calls, for speed
      msq_aq_iijj(:,:)=msq_qa_iijj(:,:)
      msq_aq_ijij(:,:)=msq_qa_ijij(:,:)
      msq_aa_ijij(:,:)=msq_qq_ijij(:,:)

C --- call for the routine returns g-fusion contribution
      call dijet_gg_ew(msq_gg_qa,ss,tt,uu)
      call dijet_gg_ew(msq_qg_qg,tt,ss,uu)
c      call dijet_gg_ew(msq_gq_gq,tt,uu,ss)
c      call dijet_gg_ew(msq_qa_gg,ss,uu,tt)
c      call dijet_gg_ew(msq_aq_gg,ss,tt,uu)
c      call dijet_gg_ew(msq_ga_ga,tt,ss,uu)
c      call dijet_gg_ew(msq_ag_ag,tt,uu,ss)
            
c--- Use symmetry relations instead of calls, for speed
c--- Note that there is a symmetry in the last two arguments of
c--- the function call that gives a high degree of symmetry
      msq_qa_gg(:)=msq_gg_qa(:)
      msq_aq_gg(:)=msq_gg_qa(:)
      msq_gq_gq(:)=msq_qg_qg(:)
      msq_ga_ga(:)=msq_qg_qg(:)
      msq_ag_ag(:)=msq_qg_qg(:)
      
C --- add initial state average here for the gluon-induced processes
      msq_gg_qa = + avegg*msq_gg_qa

      msq_qa_gg = + aveqq*msq_qa_gg
      msq_aq_gg = + aveqq*msq_aq_gg

      msq_gq_gq = - aveqg*msq_gq_gq
      msq_qg_qg = - aveqg*msq_qg_qg
      msq_ga_ga = - aveqg*msq_ga_ga
      msq_ag_ag = - aveqg*msq_ag_ag
      

      qa_iijj = 0._dp
      qa_ijij = 0._dp
      qq_ijij = 0._dp
      aq_iijj = 0._dp
      aq_ijij = 0._dp
      aa_ijij = 0._dp
c      qa_iijj2 = 0._dp
c      qa_ijij2 = 0._dp
c      qq_ijij2 = 0._dp
c      aq_iijj2 = 0._dp
c      aq_ijij2 = 0._dp
c      aa_ijij2 = 0._dp
      if(ew_tree) then
C --- LO EW Born - all the crossing relations
      call dijet_qqb_ew_tree(qa_iijj,ss,tt,uu)
      call dijet_qqb_ew_tree(qa_ijij,tt,ss,uu)
      call dijet_qqb_ew_tree(qq_ijij,tt,uu,ss)
      call dijet_qqb_ew_tree(aq_iijj,ss,tt,uu)
      call dijet_qqb_ew_tree(aq_ijij,tt,ss,uu)
      call dijet_qqb_ew_tree(aa_ijij,tt,uu,ss)

      p_qa_iijj = p
C --- qa_ij -> qa_ij: 2 <-> 3
      p_qa_ijij(1,:) = p(1,:)
      p_qa_ijij(2,:) = p(3,:)
      p_qa_ijij(3,:) = p(2,:)
      p_qa_ijij(4,:) = p(4,:)
      p_qa_ijij(5:mxpart,:) = p(5:mxpart,:)
C --- qq_ij -> qq_ij: 2 -> 3, 3 -> 4, 4 -> 2
      p_qq_ijij(1,:) = p(1,:)
      p_qq_ijij(2,:) = p(3,:)
      p_qq_ijij(3,:) = p(4,:)
      p_qq_ijij(4,:) = p(2,:)
      p_qq_ijij(5:mxpart,:) = p(5:mxpart,:)
C --- aq_ii -> aq_jj: 1 < -> 2, 3 <-> 4
      p_aq_iijj(1,:) = p(2,:)
      p_aq_iijj(2,:) = p(1,:)
      p_aq_iijj(3,:) = p(4,:)
      p_aq_iijj(4,:) = p(3,:)
      p_aq_iijj(5:mxpart,:) = p(5:mxpart,:)
C --- aq_ij -> aq_ij: 1 -> 3, 3 -> 4, 4 -> 2, 2 -> 1
      p_aq_ijij(1,:) = p(3,:)
      p_aq_ijij(2,:) = p(1,:)
      p_aq_ijij(3,:) = p(4,:)
      p_aq_ijij(4,:) = p(2,:)
      p_aq_ijij(5:mxpart,:) = p(5:mxpart,:)
C --- aa_ij -> aa_ij: 1 -> 3, 3 -> 2, 2 -> 1
      p_aa_ijij(1,:) = p(3,:)
      p_aa_ijij(2,:) = p(1,:)
      p_aa_ijij(3,:) = p(2,:)
      p_aa_ijij(4,:) = p(4,:)
      p_aa_ijij(5:mxpart,:) = p(5:mxpart,:)

      call dijet_qqb_ew_tree2(qa_iijj2,p_qa_iijj)
      call dijet_qqb_ew_tree2(qa_ijij2,p_qa_ijij)
      call dijet_qqb_ew_tree2(qq_ijij2,p_qq_ijij)
      call dijet_qqb_ew_tree2(aq_iijj2,p_aq_iijj)
      call dijet_qqb_ew_tree2(aq_ijij2,p_aq_ijij)
      call dijet_qqb_ew_tree2(aa_ijij2,p_aa_ijij)
      

C --- testing the ew_tree
c      call dijet_qqb_ew_tree(qa_iijj,ss,tt,uu)
c      call dijet_qqb_ew_tree2(qa_iijj2,p)

      write(6,*), 'testing: '
      write(6,*),qa_iijj(1,1),qa_iijj2(1,1),qa_iijj(1,1)/qa_iijj2(1,1)
      write(6,*),qq_ijij(2,3),qq_ijij2(2,3),qq_ijij(2,3)/qq_ijij2(2,3)
      write(6,*),aa_ijij(4,4),aa_ijij2(4,4),aa_ijij(4,4)/aa_ijij2(4,4)
      pause

      end if


      msq = 0._dp
      do j = -nflav,nflav
         do k = -nflav,nflav

C --- qa
            if((j > 0) .and. (k < 0)) then
               if(j == -k) then
                  msq(j,k) = msq_qa_iijj(j,1) + msq_qa_iijj(j,2) 
     &                 + msq_qa_iijj(j,3) + msq_qa_iijj(j,4)
     &                 + msq_qa_gg(j)*half
C --- LO EW
     &                 + qa_iijj(j,1) + qa_iijj(j,2) 
     &                 + qa_iijj(j,3) + qa_iijj(j,4)
                  
               else
                  msq(j,k) = msq_qa_ijij(j,-k)
C --- LO EW
     &                 + qa_ijij(j,-k)
               end if

C --- qq
            else if((j > 0) .and. (k > 0)) then
               if(j == k) then
                  msq(j,k) = msq_qq_ijij(j,k)*half
C --- LO EW
     &                 + qq_ijij(j,k)*half
               else
                  msq(j,k) = msq_qq_ijij(j,k)
C --- LO EW
     &                 + qq_ijij(j,k)
               end if
      
C --- aq
            else if((j < 0) .and. (k > 0)) then
               if(k == -j) then
                  msq(j,k) = msq_aq_iijj(k,1) + msq_aq_iijj(k,2) 
     &                 + msq_aq_iijj(k,3) + msq_aq_iijj(k,4)
     &                 + msq_aq_gg(k)*half
C --- LO EW
     &                 + aq_iijj(k,1) + aq_iijj(k,2)
     &                 + aq_iijj(k,3) + aq_iijj(k,4)
               else
                  msq(j,k) = msq_aq_ijij(-j,k)
C --- LO EW
     &                 + aq_ijij(-j,k)
               end if

C --- aa
            else if((j < 0) .and. (k < 0)) then
               if(j == k) then
                  msq(j,k) = msq_aa_ijij(-j,-k)*half
C --- LO EW
     &                 + aa_ijij(-j,-k)*half
               else
                  msq(j,k) = msq_aa_ijij(-j,-k)
C --- LO EW
     &                 + aa_ijij(-j,-k)
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


! DEBUG
c      msq(:,:) = 0._dp
c      msq(2,4) = msq_qq_ijij(2,4)
c      msq(2,1) = msq_qq_ijij(2,1)
c      msq(1,2) = msq_qq_ijij(1,2)
c      msq(2,-2) = msq_qa_iijj(2,1)
c      msq(2,2) = msq_qq_ijij(2,2)*half
c      msq(2,-3) = msq_qa_ijij(2,3)
c      msq(2,3) = msq_qq_ijij(2,3)
c      msq(3,2) = msq_qq_ijij(3,2)
c      msq(1,-1) = msq_qa_iijj(1,2)
c      msq(1,1) = msq_qq_ijij(1,1)*half
c      msq(3,3) = msq_qq_ijij(3,3)*half
c      msq(-4,-1) = msq_aa_ijij(4,1)
c      msq(-1,-4) = msq_aa_ijij(1,4)
c      msq(1,4) = msq_qq_ijij(1,4)
c-      msq(4,1) = msq_qq_ijij(4,1)
c-      msq(4,-4) = msq_qa_iijj(4,1)
c-      msq(-4,4) = msq_aq_iijj(4,1)
c      msq(1,-1) = msq_qa_iijj(1,4)
c      msq(-1,1) = msq_aq_iijj(1,4)
c      msq(2,-4) = msq_qa_ijij(2,4)
c      msq(4,-4) = msq_qa_iijj(1,4) + msq_qa_iijj(2,4) 
c     &     + msq_qa_iijj(3,4) + msq_qa_iijj(4,4)



      end subroutine qqb_twojet_ew_exact

      




      subroutine dijet_qqb_ii_jj_v(msq,ss,tt,uu)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'qcdcouple.f'
      include 'scale.f'
      real(dp):: msq(nf,nf),ss,tt,uu,mz,mw,vrt1(nf,nf),
     &     vrt2(nf,nf),vrt3(nf,nf),vrt(nf,nf),box(nf,nf),
     &     czuu(2,0:1),czdd(2,0:1),czud(2,0:1),cruu(2,0:1),crdd(2,0:1),
     &     crud(2,0:1),ccw(2,0:1),Qf(nf),
     &     propz_s,propz_t,propw,cs0,cs1,cs2,msq_cs(0:3,nf,nf)!,sxt
      integer j,k,ep,epp
      logical, parameter:: flg=.false., chk=.false.

      real(dp):: 
     &     bx1z_uu(-2:0),bx1z_uu_x(-2:0),bx1z_dd(-2:0),bx1z_dd_x(-2:0),
     &     bx1z_ud(-2:0),
     &     bx2z_uu(-2:0),bx2z_uu_x(-2:0),bx2z_dd(-2:0),bx2z_dd_x(-2:0),
     &     bx2z_ud(-2:0),
     &     bx3z_uu(-2:0),bx3z_uu_x(-2:0),bx3z_dd(-2:0),bx3z_dd_x(-2:0),
     &     bx3w_ud_x(-2:0),bx3g_ud(-2:0),
     &     bx4z_uu(-2:0),bx4z_uu_x(-2:0),bx4z_dd(-2:0),bx4z_dd_x(-2:0),
     &     bx4w_ud_x(-2:0),bx4g_ud(-2:0)

      real(dp)::
     &     bx1g_uu(2,-2:0),bx1g_uu_x(2,-2:0),bx1g_dd(2,-2:0),
     &     bx1g_dd_x(2,-2:0),bx1g_ud(2,-2:0),
     &     bx2g_uu(2,-2:0),bx2g_uu_x(2,-2:0),bx2g_dd(2,-2:0),
     &     bx2g_dd_x(2,-2:0),bx2g_ud(2,-2:0),
     &     bx3g_uu(2,-2:0),bx3g_uu_x(2,-2:0),bx3g_dd(2,-2:0),
     &     bx3g_dd_x(2,-2:0),
     &     bx4g_uu(2,-2:0),bx4g_uu_x(2,-2:0),bx4g_dd(2,-2:0),
     &     bx4g_dd_x(2,-2:0)

C --- variables for new routines return box contribution
      real(dp):: 
     &     bx1z_uu_new,bx1z_uu_x_new,bx1z_dd_new,bx1z_dd_x_new,
     &     bx1z_ud_new,
     &     bx2z_uu_new,bx2z_uu_x_new,bx2z_dd_new,bx2z_dd_x_new,
     &     bx2z_ud_new,
     &     bx3z_uu_new,bx3z_uu_x_new,bx3z_dd_new,bx3z_dd_x_new,
     &     bx3w_ud_x_new,bx3g_ud_new,
     &     bx4z_uu_new,bx4z_uu_x_new,bx4z_dd_new,bx4z_dd_x_new,
     &     bx4w_ud_x_new,bx4g_ud_new

      real(dp)::
     &     bx1g_uu_new(2),bx1g_uu_x_new(2),bx1g_dd_new(2),
     &     bx1g_dd_x_new(2),bx1g_ud_new(2),
     &     bx2g_uu_new(2),bx2g_uu_x_new(2),bx2g_dd_new(2),
     &     bx2g_dd_x_new(2),bx2g_ud_new(2),
     &     bx3g_uu_new(2),bx3g_uu_x_new(2),bx3g_dd_new(2),
     &     bx3g_dd_x_new(2),
     &     bx4g_uu_new(2),bx4g_uu_x_new(2),bx4g_dd_new(2),
     &     bx4g_dd_x_new(2)

      real(dp)::
     &     bx_type_I_uu,bx_type_I_dd,bx_type_II_ud,bx_type_III_uu,
     &     bx_type_III_dd

      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,Vcd,Vcs,Vcb


      mz = zmass
      mw = wmass

C --- initial and final Z couplings used in box
C --- initial and final are both up-type quarks
      czuu(1:2,0) = (l(2) - r(2))/2._dp
      czuu(1:2,1) = (l(2) + r(2))/2._dp
C --- I & F down-type quarks
      czdd(1:2,0) = (l(1) - r(1))/2._dp
      czdd(1:2,1) = (l(1) + r(1))/2._dp
C --- I & F up- and down-type quarks
      czud(1,0) = (l(2) - r(2))/2._dp
      czud(2,0) = (l(1) - r(1))/2._dp
      czud(1,1) = (l(2) + r(2))/2._dp
      czud(2,1) = (l(1) + r(1))/2._dp


C --- flag determines whether photon mediated amplitude is included

      Qf(1:nf) = 0._dp

      if(flg) Qf(1:nf) = Q(1:nf)

C --- photon coupling
      cruu(1:2,0) = 0._dp
      cruu(1:2,1) = (/Qf(2),Qf(2)/)

      crdd(1:2,0) = 0._dp
      crdd(1:2,1) = (/Qf(1),Qf(1)/)

      crud(1:2,0) = 0._dp
      crud(1:2,1) = (/Qf(2),Qf(1)/)

C --- W coupling
      ccw = 1._dp/2._dp/sqrt(2._dp*xw)

      msq = 0._dp


C --- vertex corrections
      call vertex1(vrt1,ss,tt,uu)
      call vertex2(vrt2,ss,tt,uu)
      call vertex3(vrt3,ss,tt,uu)

      vrt = 0._dp
      do j = 1,nf-1
         do k = 1,nf-1
            vrt(j,k) = vrt1(j,k) + vrt2(j,k) + vrt3(j,k)
         end do
      end do

C --- box corrections initialized
      box = 0._dp

C --- finite and single pole; double pole cancelled when add together
C --- for Weak-QCD boxes, the double pole cancellation is complete;
C --- however, pure QCD boxes have double pole remainings due to diff cs.
      if(1 == 2) then
      epp = 0
      if(chk) epp = -2
      do ep = epp,0
         call dijet_bx1(bx1z_uu(ep),ss,tt,uu,czuu,mz,ep)
         call dijet_bx1(bx1g_uu(1,ep),ss,tt,uu,czuu,0._dp,ep)
         call dijet_bx1(bx1g_uu(2,ep),ss,tt,uu,cruu,0._dp,ep)
         call dijet_bx1(bx1z_dd(ep),ss,tt,uu,czdd,mz,ep)
         call dijet_bx1(bx1g_dd(1,ep),ss,tt,uu,czdd,0._dp,ep)
         call dijet_bx1(bx1g_dd(2,ep),ss,tt,uu,crdd,0._dp,ep)
         call dijet_bx2(bx2z_uu(ep),ss,tt,uu,czuu,mz,ep)
         call dijet_bx2(bx2g_uu(1,ep),ss,tt,uu,czuu,0._dp,ep)
         call dijet_bx2(bx2g_uu(2,ep),ss,tt,uu,cruu,0._dp,ep)
         call dijet_bx2(bx2z_dd(ep),ss,tt,uu,czdd,mz,ep)
         call dijet_bx2(bx2g_dd(1,ep),ss,tt,uu,czdd,0._dp,ep)
         call dijet_bx2(bx2g_dd(2,ep),ss,tt,uu,crdd,0._dp,ep)

         call dijet_bx1(bx1z_ud(ep),ss,tt,uu,czud,mz,ep)
         call dijet_bx1(bx1g_ud(1,ep),ss,tt,uu,czud,0._dp,ep)
         call dijet_bx1(bx1g_ud(2,ep),ss,tt,uu,crud,0._dp,ep)
         call dijet_bx2(bx2z_ud(ep),ss,tt,uu,czud,mz,ep)
         call dijet_bx2(bx2g_ud(1,ep),ss,tt,uu,czud,0._dp,ep)
         call dijet_bx2(bx2g_ud(2,ep),ss,tt,uu,crud,0._dp,ep)
         call dijet_bx3(bx3w_ud_x(ep),tt,ss,uu,ccw,mw,ep)
         call dijet_bx3(bx3g_ud(ep),ss,tt,uu,ccw,0._dp,ep)
         call dijet_bx4(bx4w_ud_x(ep),tt,ss,uu,ccw,mw,ep)
         call dijet_bx4(bx4g_ud(ep),ss,tt,uu,ccw,0._dp,ep)

         call dijet_bx1(bx1z_uu_x(ep),tt,ss,uu,czuu,mz,ep)
         call dijet_bx1(bx1g_uu_x(1,ep),tt,ss,uu,czuu,0._dp,ep)
         call dijet_bx1(bx1g_uu_x(2,ep),tt,ss,uu,cruu,0._dp,ep)
         call dijet_bx2(bx2z_uu_x(ep),tt,ss,uu,czuu,mz,ep)
         call dijet_bx2(bx2g_uu_x(1,ep),tt,ss,uu,czuu,0._dp,ep)
         call dijet_bx2(bx2g_uu_x(2,ep),tt,ss,uu,cruu,0._dp,ep)
         call dijet_bx3(bx3z_uu(ep),ss,tt,uu,czuu,mz,ep)
         call dijet_bx3(bx3z_uu_x(ep),tt,ss,uu,czuu,mz,ep)
         call dijet_bx3(bx3g_uu(1,ep),ss,tt,uu,czuu,0._dp,ep)
         call dijet_bx3(bx3g_uu(2,ep),ss,tt,uu,cruu,0._dp,ep)
         call dijet_bx3(bx3g_uu_x(1,ep),tt,ss,uu,czuu,0._dp,ep)
         call dijet_bx3(bx3g_uu_x(2,ep),tt,ss,uu,cruu,0._dp,ep)
         call dijet_bx4(bx4z_uu(ep),ss,tt,uu,czuu,mz,ep)
         call dijet_bx4(bx4z_uu_x(ep),tt,ss,uu,czuu,mz,ep)
         call dijet_bx4(bx4g_uu(1,ep),ss,tt,uu,czuu,0._dp,ep)
         call dijet_bx4(bx4g_uu(2,ep),ss,tt,uu,cruu,0._dp,ep)
         call dijet_bx4(bx4g_uu_x(1,ep),tt,ss,uu,czuu,0._dp,ep)
         call dijet_bx4(bx4g_uu_x(2,ep),tt,ss,uu,cruu,0._dp,ep)

         call dijet_bx1(bx1z_dd_x(ep),tt,ss,uu,czdd,mz,ep)
         call dijet_bx1(bx1g_dd_x(1,ep),tt,ss,uu,czdd,0._dp,ep)
         call dijet_bx1(bx1g_dd_x(2,ep),tt,ss,uu,crdd,0._dp,ep)
         call dijet_bx2(bx2z_dd_x(ep),tt,ss,uu,czdd,mz,ep)
         call dijet_bx2(bx2g_dd_x(1,ep),tt,ss,uu,czdd,0._dp,ep)
         call dijet_bx2(bx2g_dd_x(2,ep),tt,ss,uu,crdd,0._dp,ep)
         call dijet_bx3(bx3z_dd(ep),ss,tt,uu,czdd,mz,ep)
         call dijet_bx3(bx3z_dd_x(ep),tt,ss,uu,czdd,mz,ep)
         call dijet_bx3(bx3g_dd(1,ep),ss,tt,uu,czdd,0._dp,ep)
         call dijet_bx3(bx3g_dd(2,ep),ss,tt,uu,crdd,0._dp,ep)
         call dijet_bx3(bx3g_dd_x(1,ep),tt,ss,uu,czdd,0._dp,ep)
         call dijet_bx3(bx3g_dd_x(2,ep),tt,ss,uu,crdd,0._dp,ep)
         call dijet_bx4(bx4z_dd(ep),ss,tt,uu,czdd,mz,ep)
         call dijet_bx4(bx4z_dd_x(ep),tt,ss,uu,czdd,mz,ep)
         call dijet_bx4(bx4g_dd(1,ep),ss,tt,uu,czdd,0._dp,ep)
         call dijet_bx4(bx4g_dd(2,ep),ss,tt,uu,crdd,0._dp,ep)
         call dijet_bx4(bx4g_dd_x(1,ep),tt,ss,uu,czdd,0._dp,ep)
         call dijet_bx4(bx4g_dd_x(2,ep),tt,ss,uu,crdd,0._dp,ep)
      end do

C --- the routine dijet_bxn(n=1,...,4) includes color structure
C --- when use dijet_bx3 to get pure QCD box contributions, 
C --- the color structure changes(dijet_bx4 happen to have same cs). 
C --- The following is a rescaling of this cs 
C --- (-8._dp due to conversion of cs1 = -2/3 to cs2 = 16/3)
      bx3g_uu = -8._dp*bx3g_uu
      bx3g_uu_x = -8._dp*bx3g_uu_x
      bx3g_dd = -8._dp*bx3g_dd
      bx3g_dd_x = -8._dp*bx3g_dd_x
      bx3g_ud = -8._dp*bx3g_ud

C --- the mixed QCD-Eweak box contribution should be doubled for the reason 
C --- that the two gauge boson propagators in such box do not have a symmetry
C --- factor '1/2'. In other words, these contributions with the massive gauge 
C --- boson up(left) and gluon down(right), and vice versa should be differentiated.
      bx1z_uu = 2._dp*bx1z_uu
      bx2z_uu = 2._dp*bx2z_uu
      bx1z_dd = 2._dp*bx1z_dd
      bx2z_dd = 2._dp*bx2z_dd
      bx1z_ud = 2._dp*bx1z_ud
      bx2z_ud = 2._dp*bx2z_ud
      bx1z_uu_x = 2._dp*bx1z_uu_x
      bx2z_uu_x = 2._dp*bx2z_uu_x
      bx1z_dd_x = 2._dp*bx1z_dd_x
      bx2z_dd_x = 2._dp*bx2z_dd_x
      bx3w_ud_x = 2._dp*bx3w_ud_x
      bx4w_ud_x = 2._dp*bx4w_ud_x
      bx3z_uu = 2._dp*bx3z_uu
      bx3z_uu_x = 2._dp*bx3z_uu_x
      bx3z_dd = 2._dp*bx3z_dd
      bx3z_dd_x = 2._dp*bx3z_dd_x
      bx4z_uu = 2._dp*bx4z_uu
      bx4z_uu_x = 2._dp*bx4z_uu_x
      bx4z_dd = 2._dp*bx4z_dd
      bx4z_dd_x = 2._dp*bx4z_dd_x

      propz_s = ss*(ss - mz**2)/((ss - mz**2)**2 + zwidth**2*mz**2)
      propz_t = tt*(tt - mz**2)/((tt - mz**2)**2 + zwidth**2*mz**2)
      propw = tt*(tt - mw**2)/((tt - mw**2)**2 + wwidth**2*mw**2)


      call qqb_twojet_ii_jj_mix(msq_cs,ss,tt,uu)
      cs0 = + 2._dp
      cs1 = -2._dp/3._dp
      cs2 = + 16._dp/3._dp

c      epinv = 0._dp !DEBUG
c      epinv2 = 0._dp

C --- reaction type: u_i u_i~ -> u_j u_j~ || d_i d_i~ -> d_j d_j~
      bx_type_I_uu = 
     &     + 2._dp*(
     &     + bx1z_uu(0) + bx2z_uu(0) 
     &     + (bx1g_uu(1,0) + bx2g_uu(1,0))*propz_s
     &     + bx1g_uu(2,0) + bx2g_uu(2,0)
c     &     + epinv*(
c     &     + bx1z_uu(-1) + bx2z_uu(-1) 
c     &     + (bx1g_uu(1,-1) + bx2g_uu(1,-1))*propz_s
c     &     + bx1g_uu(2,-1) + bx2g_uu(2,-1)
c     &     )
     &     )

      bx_type_I_uu = 
     &     + bx_type_I_uu 
     &     + ason4pi*8._dp*msq_cs(0,2,4)*(
     &     + epinv*( 
     &     + cs0*log(abs(tt/uu))
     &     )
     &     )


c      write(6,*) "bx_I --- 1: ", bx_type_I_uu
c      write(6,*) "bx_I --- 2: ", bx_type_I_uu
c      pause

      bx_type_I_dd = 
     &     + 2._dp*(
     &     + bx1z_dd(0) + bx2z_dd(0)
     &     + (bx1g_dd(1,0) + bx2g_dd(1,0))*propz_s
     &     + bx1g_dd(2,0) + bx2g_dd(2,0)
c     &     + epinv*(
c     &     + bx1z_dd(-1) + bx2z_dd(-1)
c     &     + (bx1g_dd(1,-1) + bx2g_dd(1,-1))*propz_s
c     &     + bx1g_dd(2,-1) +bx2g_dd(2,-1)
c     &     )
     &     )

      bx_type_I_dd = 
     &     + bx_type_I_dd
     &     + ason4pi*8._dp*msq_cs(0,1,3)*(
     &     + epinv*(
     &     + cs0*log(abs(tt/uu))
     &     )
     &     )

      box(1,3) = bx_type_I_dd
      box(2,4) = bx_type_I_uu
      box(3,1) = bx_type_I_dd
      box(4,2) = bx_type_I_uu


      if(chk) then
      write(6,*) "test the poles === type I"
      write(6,*) "uu: "
      write(6,*) 2._dp*(+ bx1z_uu(-1) + bx2z_uu(-1) 
     &     + (bx1g_uu(1,-1) + bx2g_uu(1,-1))*propz_s
     &     + bx1g_uu(2,-1) + bx2g_uu(2,-1))/
     &     (ason4pi*8._dp*msq_cs(0,2,4)*cs0*log(abs(tt/uu)))
      write(6,*) "dd: "
      write(6,*) 2._dp*(+ bx1z_dd(-1) + bx2z_dd(-1)
     &     + (bx1g_dd(1,-1) + bx2g_dd(1,-1))*propz_s
     &     + bx1g_dd(2,-1) +bx2g_dd(2,-1))/
     &     (ason4pi*8._dp*msq_cs(0,1,3)*cs0*log(abs(tt/uu)))
      pause
      end if

C --- u_i u_i~ -> d_j d_j~ || d_i d_i~ -> u_j u_j~
      bx_type_II_ud = 
     &     + 2._dp*(
     &     + bx3w_ud_x(0) + bx4w_ud_x(0)
     &     + (bx3g_ud(0) + bx4g_ud(0))*propw
c     &     + epinv*(
c     &     + bx3w_ud_x(-1) + bx4w_ud_x(-1)
c     &     + (bx3g_ud(-1) + bx4g_ud(-1))*propw
c     &     )
c     &     + epinv2**2*(
c     &     + (bx3g_ud(-2) + bx4g_ud(-2))*propw
c     &     )
     &     )

      bx_type_II_ud = 
     &     + bx_type_II_ud
     &     + ason4pi*4._dp*msq_cs(2,2,1)/Vud**2*(
     &     + epinv*(
     &     + cs1*log(abs(ss/uu)) 
     &     + cs2*log(abs(tt/musq)) 
     &     - cs1*log(abs(uu/musq))
     &     )
     &     + epinv2**2*(
     &     - cs2 + cs1
     &     )
     &     )

      box(1,2) = bx_type_II_ud*Vud**2
      box(1,4) = bx_type_II_ud*Vcd**2
      box(2,1) = bx_type_II_ud*Vud**2
      box(2,3) = bx_type_II_ud*Vus**2
      box(3,2) = bx_type_II_ud*Vus**2
      box(3,4) = bx_type_II_ud*Vcs**2
      box(4,3) = bx_type_II_ud*Vcs**2
      box(4,1) = bx_type_II_ud*Vcd**2

      bx_type_II_ud = 
     &     + 2._dp*(
     &     + bx1z_ud(0) + bx2z_ud(0)
     &     + (bx1g_ud(1,0) + bx2g_ud(1,0))*propz_s
     &     + bx1g_ud(2,0) + bx2g_ud(2,0)
c     &     + epinv*(
c     &     + bx1z_ud(-1) + bx2z_ud(-1)
c     &     + (bx1g_ud(1,-1) + bx2g_ud(1,-1))*propz_s
c     &      + bx1g_ud(2,-1) + bx2g_ud(2,-1)
c     &     )
     &     )
c      bx_type_II_ud = 0._dp

      bx_type_II_ud = 
     &     + bx_type_II_ud
     &     + ason4pi*8._dp*msq_cs(0,2,1)*(
     &     + epinv*(
     &     + cs0*log(abs(tt/uu))
     &     )
     &     )

      box(1,2) = box(1,2) + bx_type_II_ud
      box(1,4) = box(1,4) + bx_type_II_ud
      box(2,1) = box(2,1) + bx_type_II_ud
      box(2,3) = box(2,3) + bx_type_II_ud
      box(3,2) = box(3,2) + bx_type_II_ud
      box(3,4) = box(3,4) + bx_type_II_ud
      box(4,3) = box(4,3) + bx_type_II_ud
      box(4,1) = box(4,1) + bx_type_II_ud

      if(chk) then
      write(6,*) "test the poles === type II"
      write(6,*) "single pole: "
      write(6,*) "W: "
      write(6,*) 2._dp*(+ bx3w_ud_x(-1) + bx4w_ud_x(-1)
     &     + (bx3g_ud(-1) + bx4g_ud(-1))*propw)*Vud**2/
     &     (ason4pi*4._dp*msq_cs(2,2,1)*
     &     (cs1*log(abs(ss/uu)) + cs2*log(abs(tt/musq)) 
     &     - cs1*log(abs(uu/musq))))
      write(6,*) "Z: "
      write(6,*) 2._dp*(+ bx1z_ud(-1) + bx2z_ud(-1)
     &     + (bx1g_ud(1,-1) + bx2g_ud(1,-1))*propz_s
     &      + bx1g_ud(2,-1) + bx2g_ud(2,-1))/
     &     (ason4pi*8._dp*msq_cs(0,2,1)*cs0*log(abs(tt/uu)))
      write(6,*) "double pole: "
      write(6,*) 2._dp*(bx3g_ud(-2) + bx4g_ud(-2))*propw*Vud**2/
     &     (ason4pi*4._dp*msq_cs(2,2,1)*(-cs2 + cs1))
      pause
      end if

C --- u_i u_i~ -> u_i u_i~ || d_i d_i~ -> d_i d_i~
      bx_type_III_uu = 
     &     + 2._dp*(
     &     + bx1z_uu(0) + bx2z_uu(0) 
     &     + bx3z_uu(0) + bx4z_uu(0)
     &     + bx1z_uu_x(0) + bx2z_uu_x(0) 
     &     + bx3z_uu_x(0) + bx4z_uu_x(0)
     &     + (bx1g_uu(1,0) + bx2g_uu(1,0))*propz_s 
     &     + bx1g_uu(2,0) + bx2g_uu(2,0)
     &     + (bx3g_uu(1,0) + bx4g_uu(1,0))*propz_t
     &     + bx3g_uu(2,0) + bx4g_uu(2,0)
     &     + (bx1g_uu_x(1,0) + bx2g_uu_x(1,0))*propz_t
     &     + bx1g_uu_x(2,0) + bx2g_uu_x(2,0)
     &     + (bx3g_uu_x(1,0) + bx4g_uu_x(1,0))*propz_s
     &     + bx3g_uu_x(2,0) + bx4g_uu_x(2,0)
c     &     + epinv*(
c     &     + bx1z_uu(-1) + bx2z_uu(-1) 
c     &     + bx3z_uu(-1) + bx4z_uu(-1)
c     &     + bx1z_uu_x(-1) + bx2z_uu_x(-1) 
c     &     + bx3z_uu_x(-1) + bx4z_uu_x(-1)
c     &     + (bx1g_uu(1,-1) + bx2g_uu(1,-1))*propz_s
c     &     + bx1g_uu(2,-1) + bx2g_uu(2,-1)
c     &     + (bx3g_uu(1,-1) + bx4g_uu(1,-1))*propz_t
c     &     + bx3g_uu(2,-1) + bx4g_uu(2,-1)
c     &     + (bx1g_uu_x(1,-1) + bx2g_uu_x(1,-1))*propz_t
c     &     + bx1g_uu_x(2,-1) + bx2g_uu_x(2,-1)
c     &     + (bx3g_uu_x(1,-1) + bx4g_uu_x(1,-1))*propz_s
c     &     + bx3g_uu_x(2,-1) + bx4g_uu_x(2,-1)
c     &     )
c     &     + epinv2**2*(
c     &     + (bx3g_uu(1,-2) + bx4g_uu(1,-2))*propz_t
c     &     + bx3g_uu(2,-2) + bx4g_uu(2,-2)
c     &     + (bx3g_uu_x(1,-2) + bx4g_uu_x(1,-2))*propz_s
c     &     + bx3g_uu_x(2,-2) + bx4g_uu_x(2,-2)
c     &     )
     &     )

      bx_type_III_uu = 
     &     + bx_type_III_uu
     &     + epinv*(
     &     + ason4pi*8._dp*msq_cs(0,2,2)*cs0*log(abs(tt/uu))
     &     + ason4pi*8._dp*msq_cs(3,2,2)*cs0*log(abs(ss/uu))
     &     + ason4pi*4._dp*msq_cs(1,2,2)*(
     &     + cs1*log(abs(tt/uu)) 
     &     + cs2*log(abs(ss/musq))
     &     - cs1*log(abs(uu/musq))
     &     )
     &     + ason4pi*4._dp*msq_cs(2,2,2)*(
     &     + cs1*log(abs(ss/uu))
     &     + cs2*log(abs(tt/musq))
     &     - cs1*log(abs(uu/musq))
     &     )
     &     )
     &     + epinv2**2*(
     &     + ason4pi*4._dp*msq_cs(1,2,2)*(- cs2 + cs1)
     &     + ason4pi*4._dp*msq_cs(2,2,2)*(- cs2 + cs1)
     &     )

      bx_type_III_dd = 
     &     + 2._dp*(
     &     + bx1z_dd(0) + bx2z_dd(0) 
     &     + bx3z_dd(0) + bx4z_dd(0)
     &     + bx1z_dd_x(0) + bx2z_dd_x(0) 
     &     + bx3z_dd_x(0) + bx4z_dd_x(0)
     &     + (bx1g_dd(1,0) + bx2g_dd(1,0))*propz_s 
     &     + bx1g_dd(2,0) + bx2g_dd(2,0)
     &     + (bx3g_dd(1,0) + bx4g_dd(1,0))*propz_t
     &     + bx3g_dd(2,0) + bx4g_dd(2,0)
     &     + (bx1g_dd_x(1,0) + bx2g_dd_x(1,0))*propz_t
     &     + bx1g_dd_x(2,0) + bx2g_dd_x(2,0)
     &     + (bx3g_dd_x(1,0) + bx4g_dd_x(1,0))*propz_s
     &     + bx3g_dd_x(2,0) + bx4g_dd_x(2,0)
c     &     + epinv*(
c     &     + bx1z_dd(-1) + bx2z_dd(-1) 
c     &     + bx3z_dd(-1) + bx4z_dd(-1)
c     &     + bx1z_dd_x(-1) + bx2z_dd_x(-1) 
c     &     + bx3z_dd_x(-1) + bx4z_dd_x(-1)
c     &     + (bx1g_dd(1,-1) + bx2g_dd(1,-1))*propz_s
c     &     + bx1g_dd(2,-1) + bx2g_dd(2,-1)
c     &     + (bx3g_dd(1,-1) + bx4g_dd(1,-1))*propz_t
c     &     + bx3g_dd(2,-1) + bx4g_dd(2,-1)
c     &     + (bx1g_dd_x(1,-1) + bx2g_dd_x(1,-1))*propz_t
c     &     + bx1g_dd_x(2,-1) + bx2g_dd_x(2,-1)
c     &     + (bx3g_dd_x(1,-1) + bx4g_dd_x(1,-1))*propz_s
c     &     + bx3g_dd_x(2,-1) + bx4g_dd_x(2,-1)
c     &     )
c     &     + epinv2**2*(
c     &     + (bx3g_dd(1,-2) + bx4g_dd(1,-2))*propz_t
c     &     + bx3g_dd(2,-2) + bx4g_dd(2,-2)
c     &     + (bx3g_dd_x(1,-2) + bx4g_dd_x(1,-2))*propz_s
c     &     + bx3g_dd_x(2,-2) + bx4g_dd_x(2,-2)      
c     &     )
     &     )

      bx_type_III_dd = 
     &     + bx_type_III_dd
     &     + epinv*(
     &     + ason4pi*8._dp*msq_cs(0,1,1)*cs0*log(abs(tt/uu))
     &     + ason4pi*8._dp*msq_cs(3,1,1)*cs0*log(abs(ss/uu))
     &     + ason4pi*4._dp*msq_cs(1,1,1)*(
     &     + cs1*log(abs(tt/uu)) 
     &     + cs2*log(abs(ss/musq))
     &     - cs1*log(abs(uu/musq))
     &     )
     &     + ason4pi*4._dp*msq_cs(2,1,1)*(
     &     + cs1*log(abs(ss/uu))
     &     + cs2*log(abs(tt/musq))
     &     - cs1*log(abs(uu/musq))
     &     )
     &     )
     &     + epinv2**2*(
     &     + ason4pi*4._dp*msq_cs(1,1,1)*(- cs2 + cs1)
     &     + ason4pi*4._dp*msq_cs(2,1,1)*(- cs2 + cs1)
     &     )

      box(1,1) = bx_type_III_dd
      box(2,2) = bx_type_III_uu
      box(3,3) = bx_type_III_dd
      box(4,4) = bx_type_III_uu

      if(chk) then
      write(6,*) "test the poles === type III"
      write(6,*) "==========================="
      write(6,*) "uu: "
      write(6,*) "==========================="
      write(6,*) "single pole: "
      write(6,*) "sxs: "
      write(6,*) 2._dp*( + bx1z_uu(-1) + bx2z_uu(-1)
     &     + (bx1g_uu(1,-1) + bx2g_uu(1,-1))*propz_s
     &     + bx1g_uu(2,-1) + bx2g_uu(2,-1))/
     &     (ason4pi*8._dp*msq_cs(0,2,2)*cs0*log(abs(tt/uu)))
      write(6,*) "txt: "
      write(6,*) 2._dp*( + bx1z_uu_x(-1) + bx2z_uu_x(-1) 
     &     + (bx1g_uu_x(1,-1) + bx2g_uu_x(1,-1))*propz_t
     &     + bx1g_uu_x(2,-1) + bx2g_uu_x(2,-1))/
     &     (ason4pi*8._dp*msq_cs(3,2,2)*cs0*log(abs(ss/uu)))
      write(6,*) "sxt: "
      write(6,*) 2._dp*( + bx3z_uu(-1) + bx4z_uu(-1)
     &     + (bx3g_uu_x(1,-1) + bx4g_uu_x(1,-1))*propz_s
     &     + bx3g_uu_x(2,-1) + bx4g_uu_x(2,-1))/
     &     ( + ason4pi*4._dp*msq_cs(1,2,2)*(
     &     + cs1*log(abs(tt/uu)) 
     &     + cs2*log(abs(ss/musq))
     &     - cs1*log(abs(uu/musq))
     &     ))
      write(6,*) "txs: "
      write(6,*) 2._dp*( + bx3z_uu_x(-1) + bx4z_uu_x(-1)
     &     + (bx3g_uu(1,-1) + bx4g_uu(1,-1))*propz_t
     &     + bx3g_uu(2,-1) + bx4g_uu(2,-1))/
     &     (+ ason4pi*4._dp*msq_cs(2,2,2)*(
     &     + cs1*log(abs(ss/uu))
     &     + cs2*log(abs(tt/musq))
     &     - cs1*log(abs(uu/musq))
     &     ))
      write(6,*) "double pole: "
      write(6,*) "sxt: "
      write(6,*) 2._dp*( + (bx3g_uu_x(1,-2) + bx4g_uu_x(1,-2))*propz_s
     &     + bx3g_uu_x(2,-2) + bx4g_uu_x(2,-2))/
     &     (+ ason4pi*4._dp*msq_cs(1,2,2)*(- cs2 + cs1))
      write(6,*) "txs: "
      write(6,*) 2._dp*( + (bx3g_uu(1,-2) + bx4g_uu(1,-2))*propz_t
     &     + bx3g_uu(2,-2) + bx4g_uu(2,-2))/
     &     (+ ason4pi*4._dp*msq_cs(2,2,2)*(- cs2 + cs1))
      write(6,*) "==========================="
      write(6,*) "dd: "
      write(6,*) "==========================="
      write(6,*) "single pole: "
      write(6,*) "sxs: "
      write(6,*) 2._dp*( + bx1z_dd(-1) + bx2z_dd(-1)
     &     + (bx1g_dd(1,-1) + bx2g_dd(1,-1))*propz_s
     &     + bx1g_dd(2,-1) + bx2g_dd(2,-1))/
     &     (ason4pi*8._dp*msq_cs(0,1,1)*cs0*log(abs(tt/uu)))
      write(6,*) "txt: "
      write(6,*) 2._dp*( + bx1z_dd_x(-1) + bx2z_dd_x(-1) 
     &     + (bx1g_dd_x(1,-1) + bx2g_dd_x(1,-1))*propz_t
     &     + bx1g_dd_x(2,-1) + bx2g_dd_x(2,-1))/
     &     (ason4pi*8._dp*msq_cs(3,1,1)*cs0*log(abs(ss/uu)))
      write(6,*) "sxt: "
      write(6,*) 2._dp*( + bx3z_dd(-1) + bx4z_dd(-1)
     &     + (bx3g_dd_x(1,-1) + bx4g_dd_x(1,-1))*propz_s
     &     + bx3g_dd_x(2,-1) + bx4g_dd_x(2,-1))/
     &     ( + ason4pi*4._dp*msq_cs(1,1,1)*(
     &     + cs1*log(abs(tt/uu)) 
     &     + cs2*log(abs(ss/musq))
     &     - cs1*log(abs(uu/musq))
     &     ))
      write(6,*) "txs: "
      write(6,*) 2._dp*( + bx3z_dd_x(-1) + bx4z_dd_x(-1)
     &     + (bx3g_dd(1,-1) + bx4g_dd(1,-1))*propz_t
     &     + bx3g_dd(2,-1) + bx4g_dd(2,-1))/
     &     (+ ason4pi*4._dp*msq_cs(2,1,1)*(
     &     + cs1*log(abs(ss/uu))
     &     + cs2*log(abs(tt/musq))
     &     - cs1*log(abs(uu/musq))
     &     ))
      write(6,*) "double pole: "
      write(6,*) "sxt: "
      write(6,*) 2._dp*( + (bx3g_dd_x(1,-2) + bx4g_dd_x(1,-2))*propz_s
     &     + bx3g_dd_x(2,-2) + bx4g_dd_x(2,-2))/
     &     (+ ason4pi*4._dp*msq_cs(1,1,1)*(- cs2 + cs1))
      write(6,*) "txs: "
      write(6,*) 2._dp*( + (bx3g_dd(1,-2) + bx4g_dd(1,-2))*propz_t
     &     + bx3g_dd(2,-2) + bx4g_dd(2,-2))/
     &     (+ ason4pi*4._dp*msq_cs(2,1,1)*(- cs2 + cs1))
      pause
      end if

C --- begin w/ if(1 == 2) then
      end if




c      if(.false.) then

C****************************************************************C
C********************* add new routines *************************C
C****************************************************************C
      call dijet_bx1_new(bx1z_uu_new,ss,tt,uu,czuu,mz)
      call dijet_bx1_new(bx1g_uu_new(1),ss,tt,uu,czuu,0._dp)
c      call dijet_bx1_new(bx1g_uu_new(2),ss,tt,uu,cruu,0._dp)
      call dijet_bx1_new(bx1z_dd_new,ss,tt,uu,czdd,mz)
      call dijet_bx1_new(bx1g_dd_new(1),ss,tt,uu,czdd,0._dp)
c      call dijet_bx1_new(bx1g_dd_new(2),ss,tt,uu,crdd,0._dp)
      call dijet_bx2_new(bx2z_uu_new,ss,tt,uu,czuu,mz)
      call dijet_bx2_new(bx2g_uu_new(1),ss,tt,uu,czuu,0._dp)
c      call dijet_bx2_new(bx2g_uu_new(2),ss,tt,uu,cruu,0._dp)
      call dijet_bx2_new(bx2z_dd_new,ss,tt,uu,czdd,mz)
      call dijet_bx2_new(bx2g_dd_new(1),ss,tt,uu,czdd,0._dp)
c      call dijet_bx2_new(bx2g_dd_new(2),ss,tt,uu,crdd,0._dp)


      call dijet_bx1_new(bx1z_ud_new,ss,tt,uu,czud,mz)
      call dijet_bx1_new(bx1g_ud_new(1),ss,tt,uu,czud,0._dp)
c      call dijet_bx1_new(bx1g_ud_new(2),ss,tt,uu,crud,0._dp)
      call dijet_bx2_new(bx2z_ud_new,ss,tt,uu,czud,mz)
      call dijet_bx2_new(bx2g_ud_new(1),ss,tt,uu,czud,0._dp)
c      call dijet_bx2_new(bx2g_ud_new(2),ss,tt,uu,crud,0._dp)
      call dijet_bx3_new(bx3w_ud_x_new,tt,ss,uu,ccw,mw)
      call dijet_bx3_new(bx3g_ud_new,ss,tt,uu,ccw,0._dp)
      call dijet_bx4_new(bx4w_ud_x_new,tt,ss,uu,ccw,mw)
      call dijet_bx4_new(bx4g_ud_new,ss,tt,uu,ccw,0._dp)


      call dijet_bx1_new(bx1z_uu_x_new,tt,ss,uu,czuu,mz)
      call dijet_bx1_new(bx1g_uu_x_new(1),tt,ss,uu,czuu,0._dp)
c      call dijet_bx1_new(bx1g_uu_x_new(2),tt,ss,uu,cruu,0._dp)
      call dijet_bx2_new(bx2z_uu_x_new,tt,ss,uu,czuu,mz)
      call dijet_bx2_new(bx2g_uu_x_new(1),tt,ss,uu,czuu,0._dp)
c      call dijet_bx2_new(bx2g_uu_x_new(2),tt,ss,uu,cruu,0._dp)
      call dijet_bx3_new(bx3z_uu_new,ss,tt,uu,czuu,mz)
      call dijet_bx3_new(bx3z_uu_x_new,tt,ss,uu,czuu,mz)
      call dijet_bx3_new(bx3g_uu_new(1),ss,tt,uu,czuu,0._dp)
c      call dijet_bx3_new(bx3g_uu_new(2),ss,tt,uu,cruu,0._dp)
      call dijet_bx3_new(bx3g_uu_x_new(1),tt,ss,uu,czuu,0._dp)
c      call dijet_bx3_new(bx3g_uu_x_new(2),tt,ss,uu,cruu,0._dp)
      call dijet_bx4_new(bx4z_uu_new,ss,tt,uu,czuu,mz)
      call dijet_bx4_new(bx4z_uu_x_new,tt,ss,uu,czuu,mz)
      call dijet_bx4_new(bx4g_uu_new(1),ss,tt,uu,czuu,0._dp)
c      call dijet_bx4_new(bx4g_uu_new(2),ss,tt,uu,cruu,0._dp)
      call dijet_bx4_new(bx4g_uu_x_new(1),tt,ss,uu,czuu,0._dp)
c      call dijet_bx4_new(bx4g_uu_x_new(2),tt,ss,uu,cruu,0._dp)

      call dijet_bx1_new(bx1z_dd_x_new,tt,ss,uu,czdd,mz)
      call dijet_bx1_new(bx1g_dd_x_new(1),tt,ss,uu,czdd,0._dp)
c      call dijet_bx1_new(bx1g_dd_x_new(2),tt,ss,uu,crdd,0._dp)
      call dijet_bx2_new(bx2z_dd_x_new,tt,ss,uu,czdd,mz)
      call dijet_bx2_new(bx2g_dd_x_new(1),tt,ss,uu,czdd,0._dp)
c      call dijet_bx2_new(bx2g_dd_x_new(2),tt,ss,uu,crdd,0._dp)
      call dijet_bx3_new(bx3z_dd_new,ss,tt,uu,czdd,mz)
      call dijet_bx3_new(bx3z_dd_x_new,tt,ss,uu,czdd,mz)
      call dijet_bx3_new(bx3g_dd_new(1),ss,tt,uu,czdd,0._dp)
c      call dijet_bx3_new(bx3g_dd_new(2),ss,tt,uu,crdd,0._dp)
      call dijet_bx3_new(bx3g_dd_x_new(1),tt,ss,uu,czdd,0._dp)
c      call dijet_bx3_new(bx3g_dd_x_new(2),tt,ss,uu,crdd,0._dp)
      call dijet_bx4_new(bx4z_dd_new,ss,tt,uu,czdd,mz)
      call dijet_bx4_new(bx4z_dd_x_new,tt,ss,uu,czdd,mz)
      call dijet_bx4_new(bx4g_dd_new(1),ss,tt,uu,czdd,0._dp)
c      call dijet_bx4_new(bx4g_dd_new(2),ss,tt,uu,crdd,0._dp)
      call dijet_bx4_new(bx4g_dd_x_new(1),tt,ss,uu,czdd,0._dp)
c      call dijet_bx4_new(bx4g_dd_x_new(2),tt,ss,uu,crdd,0._dp)

      bx3g_uu_new = -8._dp*bx3g_uu_new
      bx3g_uu_x_new = -8._dp*bx3g_uu_x_new
      bx3g_dd_new = -8._dp*bx3g_dd_new
      bx3g_dd_x_new = -8._dp*bx3g_dd_x_new
      bx3g_ud_new = -8._dp*bx3g_ud_new

      bx1z_uu_new = 2._dp*bx1z_uu_new
      bx2z_uu_new = 2._dp*bx2z_uu_new
      bx1z_dd_new = 2._dp*bx1z_dd_new
      bx2z_dd_new = 2._dp*bx2z_dd_new
      bx1z_ud_new = 2._dp*bx1z_ud_new
      bx2z_ud_new = 2._dp*bx2z_ud_new
      bx1z_uu_x_new = 2._dp*bx1z_uu_x_new
      bx2z_uu_x_new = 2._dp*bx2z_uu_x_new
      bx1z_dd_x_new = 2._dp*bx1z_dd_x_new
      bx2z_dd_x_new = 2._dp*bx2z_dd_x_new
      bx3w_ud_x_new = 2._dp*bx3w_ud_x_new
      bx4w_ud_x_new = 2._dp*bx4w_ud_x_new
      bx3z_uu_new = 2._dp*bx3z_uu_new
      bx3z_uu_x_new = 2._dp*bx3z_uu_x_new
      bx3z_dd_new = 2._dp*bx3z_dd_new
      bx3z_dd_x_new = 2._dp*bx3z_dd_x_new
      bx4z_uu_new = 2._dp*bx4z_uu_new
      bx4z_uu_x_new = 2._dp*bx4z_uu_x_new
      bx4z_dd_new = 2._dp*bx4z_dd_new
      bx4z_dd_x_new = 2._dp*bx4z_dd_x_new

      propz_s = ss*(ss - mz**2)/((ss - mz**2)**2 + zwidth**2*mz**2)
      propz_t = tt*(tt - mz**2)/((tt - mz**2)**2 + zwidth**2*mz**2)
      propw = tt*(tt - mw**2)/((tt - mw**2)**2 + wwidth**2*mw**2)

      bx_type_I_uu = 
     &     + 2._dp*(
     &     + bx1z_uu_new + bx2z_uu_new
     &     + (bx1g_uu_new(1) + bx2g_uu_new(1))*propz_s
c     &     + bx1g_uu_new(2) + bx2g_uu_new(2)
     &     )


      bx_type_I_dd = 
     &     + 2._dp*(
     &     + bx1z_dd_new + bx2z_dd_new
     &     + (bx1g_dd_new(1) + bx2g_dd_new(1))*propz_s
c     &     + bx1g_dd_new(2) + bx2g_dd_new(2)
     &     )

      box(1,3) = bx_type_I_dd
      box(2,4) = bx_type_I_uu
      box(3,1) = bx_type_I_dd
      box(4,2) = bx_type_I_uu

      bx_type_II_ud = 
     &     + 2._dp*(
     &     + bx3w_ud_x_new + bx4w_ud_x_new
     &     + (bx3g_ud_new + bx4g_ud_new)*propw
     &     ) 

      box(1,2) = bx_type_II_ud*Vud**2
      box(1,4) = bx_type_II_ud*Vcd**2
      box(2,1) = bx_type_II_ud*Vud**2
      box(2,3) = bx_type_II_ud*Vus**2
      box(3,2) = bx_type_II_ud*Vus**2
      box(3,4) = bx_type_II_ud*Vcs**2
      box(4,3) = bx_type_II_ud*Vcs**2
      box(4,1) = bx_type_II_ud*Vcd**2

      bx_type_II_ud = 
     &     + 2._dp*(
     &     + bx1z_ud_new + bx2z_ud_new
     &     + (bx1g_ud_new(1) + bx2g_ud_new(1))*propz_s
c     &     + bx1g_ud_new(2) + bx2g_ud_new(2)
     &     )

      box(1,2) = box(1,2) + bx_type_II_ud
      box(1,4) = box(1,4) + bx_type_II_ud
      box(2,1) = box(2,1) + bx_type_II_ud
      box(2,3) = box(2,3) + bx_type_II_ud
      box(3,2) = box(3,2) + bx_type_II_ud
      box(3,4) = box(3,4) + bx_type_II_ud
      box(4,3) = box(4,3) + bx_type_II_ud
      box(4,1) = box(4,1) + bx_type_II_ud

      bx_type_III_uu = 
     &     + 2._dp*(
     &     + bx1z_uu_new + bx2z_uu_new 
     &     + bx3z_uu_new + bx4z_uu_new
     &     + bx1z_uu_x_new + bx2z_uu_x_new 
     &     + bx3z_uu_x_new + bx4z_uu_x_new
     &     + (bx1g_uu_new(1) + bx2g_uu_new(1))*propz_s 
c     &     + bx1g_uu_new(2) + bx2g_uu_new(2)
     &     + (bx3g_uu_new(1) + bx4g_uu_new(1))*propz_t
c     &     + bx3g_uu_new(2) + bx4g_uu_new(2)
     &     + (bx1g_uu_x_new(1) + bx2g_uu_x_new(1))*propz_t
c     &     + bx1g_uu_x_new(2) + bx2g_uu_x_new(2)
     &     + (bx3g_uu_x_new(1) + bx4g_uu_x_new(1))*propz_s
c     &     + bx3g_uu_x_new(2) + bx4g_uu_x_new(2)
     &     )

      bx_type_III_dd = 
     &     + 2._dp*(
     &     + bx1z_dd_new + bx2z_dd_new 
     &     + bx3z_dd_new + bx4z_dd_new
     &     + bx1z_dd_x_new + bx2z_dd_x_new 
     &     + bx3z_dd_x_new + bx4z_dd_x_new
     &     + (bx1g_dd_new(1) + bx2g_dd_new(1))*propz_s 
c     &     + bx1g_dd_new(2) + bx2g_dd_new(2)
     &     + (bx3g_dd_new(1) + bx4g_dd_new(1))*propz_t
c     &     + bx3g_dd_new(2) + bx4g_dd_new(2)
     &     + (bx1g_dd_x_new(1) + bx2g_dd_x_new(1))*propz_t
c     &     + bx1g_dd_x_new(2) + bx2g_dd_x_new(2)
     &     + (bx3g_dd_x_new(1) + bx4g_dd_x_new(1))*propz_s
c     &     + bx3g_dd_x_new(2) + bx4g_dd_x_new(2)
     &     )

      box(1,1) = bx_type_III_dd
      box(2,2) = bx_type_III_uu
      box(3,3) = bx_type_III_dd
      box(4,4) = bx_type_III_uu

c      end if
c      msq = vrt

      msq = vrt + box

C --- DEBUG
c      msq = 0._dp
c      box(2,2) = epinv*(
c     &     + (bx1g_uu_x(1,-1) + bx2g_uu_x(1,-1))*propz_t
c     &     + bx1g_uu_x(2,-1) + bx2g_uu_x(2,-1)
c     &     + (bx3g_uu_x(1,-1) + bx4g_uu_x(1,-1))*propz_s
c     &     + bx3g_uu_x(2,-1) + bx4g_uu_x(2,-1))
c     &     epinv2**2*(
c     &     + (bx3g_uu_x(1,-2) + bx4g_uu_x(1,-2))*propz_s
c     &     + bx3g_uu_x(2,-2) + bx4g_uu_x(2,-2))
c      msq(2,2) = (vrt(2,2) + box(2,2))
c      msq(2,2) =  vrt3(2,2) + box(2,2)
c      box(1,1) = 
c     &     + epinv*(
c     &     + bx1z_dd(-1) + bx2z_dd(-1) 
c     &     + bx1z_dd_x(-1) + bx2z_dd_x(-1) 
c     &     + (bx1g_dd(1,-1) + bx2g_dd(1,-1))*propz_s
c     &     + bx1g_dd(2,-1) + bx2g_dd(2,-1)
c     &     + (bx1g_dd_x(1,-1) + bx2g_dd_x(1,-1))*propz_t
c     &     + bx1g_dd_x(2,-1) + bx2g_dd_x(2,-1)
c     &     )
c      msq(1,1) = box(1,2)

      end subroutine dijet_qqb_ii_jj_v
