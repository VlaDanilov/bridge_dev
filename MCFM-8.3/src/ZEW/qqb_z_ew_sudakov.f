      subroutine qqb_z_ew_sudakov(p,msq)
      implicit none
C---- Author Jia Zhou
C---- December 2013
c----Matrix element for Z production
C----averaged over initial colours and spins
c     q(-p1)+qbar(-p2)-->(e^-(p3)+e^+(p4))
c---
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'cplx.h'
      integer j,k,i,m,legs(4),chs(4,4)
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),s,fac,s34,
     & ss,tt,uu,ddl,ssl,sudakov(4,2),DL,SL,dd(7),Rff2,lrr(4,4,2),
     & aem0,alpha,aemmz
      complex(dp):: prop,qqb,qbq
      common/em/aemmz
      
c---statement function
      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))

c--set msq=0 to initialize
      msq(:,:)=0._dp
      s34=s(3,4)
      ss=s34
      tt=-s(1,3)
      uu=-s(2,3)

      fac=4._dp*esq**2*xn

      call mandelstam_var(lrr(:,:,1),ss,tt,uu,s34)
      call mandelstam_var(lrr(:,:,2),ss,uu,tt,s34)

c--- PDG code for legs 3 and 4 (muons)
      legs(3:4)=(/13,-13/)
c--- Chirality for legs (+1 left, -1 right)
      chs(:,1)=(/ 1, 1, 1, 1/)
      chs(:,2)=(/-1,-1,-1,-1/)
      chs(:,3)=(/ 1, 1,-1,-1/)
      chs(:,4)=(/-1,-1, 1, 1/)

c--- alpha_em(0)
c      aem0=1._dp/127.9_dp
c--- alpha -> alpha_em(0)
c      alpha=aem0
c--- alpha -> standard MCFM value
      alpha=aemmz
    
      DL=alpha/4._dp/pi*log(ss/wmass**2)**2
      SL=alpha/4._dp/pi*log(ss/wmass**2)
c      SL=0._dp ! DEBUG: remove single logs
      
c--   calculate propagators
      fac=aveqq*fac/s34**2
      prop=s34/cplx2((s34-zmass**2),zmass*zwidth)
      prop=cplx2(one,zip)

      call spinoru(4,p,za,zb)

c---case qbar-q or q-qbar
c      qqb=fac*s(1,4)**2
c      qbq=fac*s(2,4)**2
      qqb=za(2,3)*zb(4,1)
      qbq=za(1,3)*zb(4,2)

c--- correction to amplitude squared is  2 * delta * |M|^2
c--- where delta is given by sudakov(i,j)
      fac=fac*2._dp
      
      do j=-nf,nf
      k=-j
      legs(1:2)=(/-abs(j),abs(j)/)
      do i=1,4
      do m=1,2
        call sudakov_DY(dd,Rff2,legs,chs(:,i),lrr(:,:,m))
        ddl=(dd(1)-dd(2))*DL
        ssl=(dd(3)-dd(4)+dd(6)-dd(7))*SL
        sudakov(i,m)=ddl+ssl
      end do
      end do
         if ((j .eq. 0) .and. (k .eq. 0)) then
           msq(j,k)=0._dp
         elseif ((j .gt. 0) .and. (k .lt. 0)) then
           msq(j,k)=+sudakov(1,1)*abs((Q(j)*q1+L(j)*l1*prop)*qqb)**2
     .              +sudakov(2,1)*abs((Q(j)*q1+R(j)*r1*prop)*qqb)**2
     .              +sudakov(3,1)*abs((Q(j)*q1+L(j)*r1*prop)*qbq)**2
     .              +sudakov(4,1)*abs((Q(j)*q1+R(j)*l1*prop)*qbq)**2
         elseif ((j .lt. 0) .and. (k .gt. 0)) then
           msq(j,k)=+sudakov(1,2)*abs((Q(k)*q1+L(k)*l1*prop)*qbq)**2
     .              +sudakov(2,2)*abs((Q(k)*q1+R(k)*r1*prop)*qbq)**2
     .              +sudakov(3,2)*abs((Q(k)*q1+L(k)*r1*prop)*qqb)**2
     .              +sudakov(4,2)*abs((Q(k)*q1+R(k)*l1*prop)*qqb)**2
         endif
         msq(j,k)=msq(j,k)*fac
c      write(6,*) 'j,k,msq(j,k)',j,k,msq(j,k)
      enddo
c      pause

      return
      end
