      subroutine qqb_waj_g(p,msq)
      implicit none
      include 'types.f'

c------------------------------------------------------------------------------
c---- Author: R. Mondini
c---- January 2017
c----
c---- return matrix element squared for Wgamjj production
c---- averaged over initial colors and spins
c----
c---- for nwz=+1
c----    u(-p1)+dbar(-p2)-->W^+(n(p3)+e^+(p4)) + gamma(p5) + j(p6) + j(p7)
c---- for nwz=-1
c----    ubar(-p1)+d(-p2)-->W^-(e^-(p3)+nbar(p4)) + gamma(p5) + j(p6) + j(p7)
c---- where jj=gg,QQB
c---- and crossed channels
c------------------------------------------------------------------------------

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ckm.f'
      include 'nwz.f'

      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      real(dp):: msq3(-nf:nf,-nf:nf)
      real(dp):: qqb_gg,qbq_gg,gg_qqb
      real(dp):: qg_qg,qbg_qbg,gq_gq,gqb_gqb

c      real(dp):: iib_jjb(4),ibi_jjb(4)
c      real(dp):: iib_iib(2),ibi_iib(2),ii_ii(2),ibib_ibib(2)
c      real(dp):: qiqj(2,2),qjqi(2,2),qbiqbj(2,2),qbjqbi(2,2)
c      real(dp):: qiqbj(2,2),qjqbi(2,2),qbiqj(2,2),qbjqi(2,2)

      integer:: i,j,k

c      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
c      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c-----call squared matrix elements from each subprocesses
      call msq_wqqbgamgg(p,qqb_gg,qbq_gg,gg_qqb,
     .qg_qg,qbg_qbg,gq_gq,gqb_gqb)

c      call msq_zqqbQQbgam(p,iib_jjb,ibi_jjb,
c     .iib_iib,ibi_iib,ii_ii,ibib_ibib,
c     .qiqj,qbiqbj,qiqbj,qbiqj,qjqi,qbjqbi,qjqbi,qbjqi)

c-----initialize msq
      do j=-nf,nf
      do i=-nf,nf
         msq(i,j)=zip
         msq1(i,j)=zip
         msq2(i,j)=zip
         msq3(i,j)=zip
      enddo
      enddo

c-----fill msq1 and msq3
      do k=-nf,nf
      do j=-nf,nf

c---------fill msq1
          if ((j == 0) .and. (k == 0) .and. (nwz == +1)) then
            msq1(j,k)=(Vsum(2)+Vsum(4))*gg_qqb
          elseif ((j == 0) .and. (k == 0) .and. (nwz == -1)) then
            msq1(j,k)=(Vsum(-2)+Vsum(-4))*gg_qqb
          elseif ((j == 0) .and. (k < 0)) then
            msq1(j,k)=Vsum(k)*gqb_gqb
          elseif ((j == 0) .and. (k > 0)) then
            msq1(j,k)=Vsum(k)*gq_gq
          elseif ((j > 0) .and. (k < 0)) then
            msq1(j,k)=Vsq(j,k)*qqb_gg
          elseif ((j < 0) .and. (k > 0)) then
            msq1(j,k)=Vsq(j,k)*qbq_gg
          elseif ((j > 0) .and. (k == 0)) then
            msq1(j,k)=Vsum(j)*qg_qg
          elseif ((j < 0) .and. (k == 0)) then
            msq1(j,k)=Vsum(j)*qbg_qbg
          else
            msq1(j,k)=0._dp
          endif

c---------fill msq3
c          if ((j>0).and.(k>0)) then
c            if (j<k) then
c               msq3(j,k)=qiqj(jj(j),kk(k))
c            elseif (k<j) then
c               msq3(j,k)=qjqi(kk(k),jj(j))
c            endif
c          elseif ((j<0).and.(k<0)) then
c            if (j<k) then
c               msq3(j,k)=qbjqbi(-kk(k),-jj(j))
c            elseif (k<j) then
c               msq3(j,k)=qbiqbj(-jj(j),-kk(k))
c            endif
c          elseif ((j>0).and.(k<0)) then
c            if (j<abs(k)) then
c               msq3(j,k)=qiqbj(jj(j),-kk(k))
c            elseif (abs(k)<j) then
c               msq3(j,k)=qjqbi(-kk(k),jj(j))
c            endif
c          elseif ((j<0).and.(k>0)) then
c            if (abs(j)<k) then
c               msq3(j,k)=qbiqj(-jj(j),kk(k))
c            elseif (k<abs(j)) then
c               msq3(j,k)=qbjqi(kk(k),-jj(j))
c            endif
c          endif

      enddo
      enddo

c-----fill msq2
c      msq2(2,-2)=iib_iib(2)+iib_jjb(3)+3._dp*iib_jjb(1)
c      msq2(4,-4)=msq2(2,-2)
c      msq2(1,-1)=iib_iib(1)+2._dp*iib_jjb(2)+2._dp*iib_jjb(4)
c      msq2(3,-3)=msq2(1,-1)
c      msq2(5,-5)=msq2(1,-1)
c-----
c      msq2(-2,2)=ibi_iib(2)+ibi_jjb(3)+3._dp*ibi_jjb(1)
c      msq2(-4,4)=msq2(-2,2)
c      msq2(-1,1)=ibi_iib(1)+2._dp*ibi_jjb(2)+2._dp*ibi_jjb(4)
c      msq2(-3,3)=msq2(-1,1)
c      msq2(-5,5)=msq2(-1,1)
c-----
c      msq2(2,2)=ii_ii(2)
c      msq2(4,4)=msq2(2,2)
c      msq2(1,1)=ii_ii(1)
c      msq2(3,3)=msq2(1,1)
c      msq2(5,5)=msq2(1,1)
c-----
c      msq2(-2,-2)=ibib_ibib(2)
c      msq2(-4,-4)=msq2(-2,-2)
c      msq2(-1,-1)=ibib_ibib(1)
c      msq2(-3,-3)=msq2(-1,-1)
c      msq2(-5,-5)=msq2(-1,-1)

c-----sum msq1,2,3
      do j=-nf,nf
      do i=-nf,nf
         msq(i,j)=msq1(i,j)+msq2(i,j)+msq3(i,j)
      enddo
      enddo

      return
      end


      subroutine msq_wqqbgamgg(p,qqb_gg,qbq_gg,gg_qqb,
     & qg_qg,qbg_qbg,gq_gq,gqb_gqb)
      implicit none
      include 'types.f'

*****************************************************************
* return averaged matrix element squared for                    *
* 0 -> f(p1) + f(p2) + l(p3) + vb(p4) + gam(p5) + f(p6) + f(p7) *
* coming from (0->q+qb+l+vb+gam+g+g) amplitude                  *
*****************************************************************
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nwz.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      real(dp):: p(mxpart,4)
      integer:: i,j,k
      complex(dp):: a70qq(2,2,2,2),a70ql(2,2,2,2)
      real(dp):: qqb_gg,qbq_gg,gg_qqb
      real(dp):: qg_qg,qbg_qbg,gq_gq,gqb_gqb

c-----initialize matelem
      qqb_gg=zip
      qbq_gg=zip
      gg_qqb=zip
      qg_qg=zip
      qbg_qbg=zip
      gq_gq=zip
      gqb_gqb=zip

c-----calculate spinor products, invariants
      call spinoru(7,p,za,zb)

      if (nwz == -1) then
c-----call color ordered amplitudes and compute msq for qbq_gg
      call wagg_a70h(1,2,3,4,5,6,7,za,zb,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qbq_gg)
c-----call color ordered amplitudes and compute msq for qqb_gg
      call wagg_a70h(2,1,3,4,5,6,7,za,zb,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qqb_gg)
c-----call color ordered amplitudes and compute msq for gq_gq
      call wagg_a70h(6,2,3,4,5,1,7,za,zb,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,gq_gq)
c-----call color ordered amplitudes and compute msq for qg_qg
      call wagg_a70h(6,1,3,4,5,2,7,za,zb,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qg_qg)
c-----call color ordered amplitudes and compute msq for gqb_gqb
      call wagg_a70h(2,6,3,4,5,1,7,za,zb,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,gqb_gqb)
c-----call color ordered amplitudes and compute msq for qbg_qbg
      call wagg_a70h(1,6,3,4,5,2,7,za,zb,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qbg_qbg)
c-----call color ordered amplitudes and compute msq for gg_qqb
      call wagg_a70h(6,7,3,4,5,1,2,za,zb,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,gg_qqb)

      elseif (nwz == +1) then 
c-----call color ordered amplitudes and compute msq for qbq_gg
      call wagg_a70h(2,1,4,3,5,6,7,zb,za,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qbq_gg)
c-----call color ordered amplitudes and compute msq for qqb_gg
      call wagg_a70h(1,2,4,3,5,6,7,zb,za,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qqb_gg)
c-----call color ordered amplitudes and compute msq for gq_gq
      call wagg_a70h(2,6,4,3,5,1,7,zb,za,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,gq_gq)
c-----call color ordered amplitudes and compute msq for qg_qg
      call wagg_a70h(1,6,4,3,5,2,7,zb,za,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qg_qg)
c-----call color ordered amplitudes and compute msq for gqb_gqb
      call wagg_a70h(6,2,4,3,5,1,7,zb,za,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,gqb_gqb)
c-----call color ordered amplitudes and compute msq for qbg_qbg
      call wagg_a70h(6,1,4,3,5,2,7,zb,za,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,qbg_qbg)
c-----call color ordered amplitudes and compute msq for gg_qqb
      call wagg_a70h(7,6,4,3,5,1,2,zb,za,a70qq,a70ql)
      call wagg_m70sq(a70qq,a70ql,gg_qqb)

      endif

c-----averaging and put identical particle factor
      qqb_gg=qqb_gg*half*aveqq
      qbq_gg=qbq_gg*half*aveqq
      gg_qqb=gg_qqb*avegg
      qg_qg=qg_qg*aveqg
      qbg_qbg=qbg_qbg*aveqg
      gq_gq=gq_gq*aveqg
      gqb_gqb=gqb_gqb*aveqg

      return
      end


c      subroutine msq_zqqbQQbgam(p,iib_jjb,ibi_jjb,
c     & iib_iib,ibi_iib,ii_ii,ibib_ibib,
c     & qiqj,qbiqbj,qiqbj,qbiqj,qjqi,qbjqbi,qjqbi,qbjqi)
c      implicit none
c      include 'types.f'
c*****************************************************************
c* return matrix element squared for                             *
c* 0 -> f(p1) + f(p2) + l(p3) + lb(p4) + gam(p5) + f(p6) + f(p7) *
c* coming from (0->q+qb+Q+Qb+gam+lb+l) amplitude                 *
c*****************************************************************
c      
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'qcdcouple.f'
c      include 'ewcouple.f'
c      include 'zcouple.f'
c      include 'zprods_com.f'
c      include 'sprods_com.f'
c      include 'masses.f'
c      include 'ewcharge.f'
c      real(dp):: p(mxpart,4),pnt(mxpart,4)
c      real(dp):: iib_jjb(4),ibi_jjb(4)
c      real(dp):: iib_iib(2),ibi_iib(2),ii_ii(2),ibib_ibib(2)
c      real(dp):: qiqj(2,2),qjqi(2,2),qbiqbj(2,2),qbjqbi(2,2)
c      real(dp):: qiqbj(2,2),qjqbi(2,2),qbiqj(2,2),qbjqi(2,2)
c      real(dp):: xqiqj(4),xqjqi(4),xqbiqbj(4),xqbjqbi(4)
c      real(dp):: xqiqbj(4),xqjqbi(4),xqbiqj(4),xqbjqi(4)
c      integer:: i,j,k
cc-----convert to Nagy-Trocsanyi momentum convention
c      do i=1,4
c         pnt(1,i)=p(2,i)
c         pnt(2,i)=p(1,i)
c         pnt(3,i)=p(6,i)
c         pnt(4,i)=p(7,i)
c         pnt(5,i)=p(5,i)
c         pnt(6,i)=p(4,i)
c         pnt(7,i)=p(3,i)
c      enddo
cc-----initialize msq
c      do i=1,4
c         iib_jjb(i)=zip
c         ibi_jjb(i)=zip
c      enddo
c      do i=1,2
c         iib_iib(i)  =zip
c         ibi_iib(i)  =zip
c         ii_ii(i)    =zip
c         ibib_ibib(i)=zip
c      enddo
c      do i=1,2
c      do j=1,2
c         qiqj(i,j)=zip
c         qbiqbj(i,j)=zip
c         qiqbj(i,j)=zip
c         qbiqj(i,j)=zip
c         qjqi(i,j)=zip
c         qbjqbi(i,j)=zip
c         qjqbi(i,j)=zip
c         qbjqi(i,j)=zip
c      enddo
c      enddo
cc-----evaluate spinor products, invariants
c      call spinoru(7,pnt,za,zb)
cc-----
c      call zqqbQQba_msqij(1,2,3,4,5,6,7,iib_jjb) 
cc-----
c      call zqqbQQba_msqij(2,1,3,4,5,6,7,ibi_jjb) 
cc-----
c      call zqqbQQba_msqii(1,2,3,4,5,6,7,iib_iib) 
cc-----
c      call zqqbQQba_msqii(2,1,3,4,5,6,7,ibi_iib) 
cc-----
c      call zqqbQQba_msqii(3,1,4,2,5,6,7,ii_ii) 
cc-----
c      call zqqbQQba_msqii(1,3,2,4,5,6,7,ibib_ibib)
cc-----
c      call zqqbQQba_msqij(3,2,4,1,5,6,7,xqiqj)
c      call zqqbQQba_msqij(4,1,3,2,5,6,7,xqjqi)
c      qiqj(1,1)=xqiqj(4)
c      qiqj(1,2)=xqjqi(1)
c      qiqj(2,1)=xqiqj(1)
c      qiqj(2,2)=xqiqj(3)
c      qjqi(1,1)=xqjqi(4)
c      qjqi(1,2)=xqiqj(1)
c      qjqi(2,1)=xqjqi(1)
c      qjqi(2,2)=xqjqi(3)
cc-----
c      call zqqbQQba_msqij(2,3,1,4,5,6,7,xqbiqbj)
c      call zqqbQQba_msqij(1,4,2,3,5,6,7,xqbjqbi)
c      qbiqbj(1,1)=xqbiqbj(4)
c      qbiqbj(1,2)=xqbjqbi(1)
c      qbiqbj(2,1)=xqbiqbj(1)
c      qbiqbj(2,2)=xqbiqbj(3)
c      qbjqbi(1,1)=xqbjqbi(4)
c      qbjqbi(1,2)=xqbiqbj(1)
c      qbjqbi(2,1)=xqbjqbi(1)
c      qbjqbi(2,2)=xqbjqbi(3)
cc-----
c      call zqqbQQba_msqij(3,2,1,4,5,6,7,xqiqbj) 
c      call zqqbQQba_msqij(1,4,3,2,5,6,7,xqjqbi) 
c      qiqbj(1,1)=xqiqbj(4)
c      qiqbj(1,2)=xqjqbi(1)
c      qiqbj(2,1)=xqiqbj(1)
c      qiqbj(2,2)=xqiqbj(3)
c      qjqbi(1,1)=xqjqbi(4)
c      qjqbi(1,2)=xqiqbj(1)
c      qjqbi(2,1)=xqjqbi(1)
c      qjqbi(2,2)=xqjqbi(3)
cc-----
c      call zqqbQQba_msqij(2,4,3,1,5,6,7,xqbiqj) 
c      call zqqbQQba_msqij(3,1,2,4,5,6,7,xqbjqi)
c      qbiqj(1,1)=xqbiqj(4)
c      qbiqj(1,2)=xqbjqi(1)
c      qbiqj(2,1)=xqbiqj(1)
c      qbiqj(2,2)=xqbiqj(3)
c      qbjqi(1,1)=xqbjqi(4)
c      qbjqi(1,2)=xqbiqj(1)
c      qbjqi(2,1)=xqbjqi(1)
c      qbjqi(2,2)=xqbjqi(3)
cc-----averaging and put identical particle factors
c      do i=1,4
c         iib_jjb(i)=iib_jjb(i)*aveqq
c         ibi_jjb(i)=ibi_jjb(i)*aveqq
c      enddo
c      do i=1,2
c         iib_iib(i)  =iib_iib(i)*aveqq
c         ibi_iib(i)  =ibi_iib(i)*aveqq
c         ii_ii(i)    =ii_ii(i)*aveqq*half
c         ibib_ibib(i)=ibib_ibib(i)*aveqq*half
c      enddo
c      do i=1,2
c      do j=1,2
c         qiqj(i,j)=qiqj(i,j)*aveqq
c         qjqi(i,j)=qjqi(i,j)*aveqq
c         qbiqbj(i,j)=qbiqbj(i,j)*aveqq
c         qbjqbi(i,j)=qbjqbi(i,j)*aveqq
c         qiqbj(i,j)=qiqbj(i,j)*aveqq
c         qbiqj(i,j)=qbiqj(i,j)*aveqq
c         qjqbi(i,j)=qjqbi(i,j)*aveqq
c         qbjqi(i,j)=qbjqi(i,j)*aveqq
c      enddo
c      enddo
cc-----done
c      return
c      end

