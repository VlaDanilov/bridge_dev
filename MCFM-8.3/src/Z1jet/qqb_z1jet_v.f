      subroutine qqb_z1jet_v(p,msq)
      implicit none
************************************************************************
*     Authors: R.K. Ellis and John Campbell                            *
*     May, 2001.                                                       *
*     Matrix element for Z + jet production                            *
*     in order alpha_s^2                                               *
*     averaged over initial colours and spins                          *
*     q(-p1)+qbar(-p2)-->Z^+(l(p3)+a(p4))+g(p5)                        *
************************************************************************
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'epinv.f'
      include 'scheme.f'
      include 'nflav.f'
      include 'cplx.h'
      integer:: j,k,nu,om
      real(dp):: msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf),
     & p(mxpart,4),fac,sz,virt5,subuv,sin2winv
      complex(dp):: qqbZgLL(2),qqbZgRR(2),qqbZgLR(2),qqbZgRL(2)
      complex(dp):: gqZqLL(2),gqZqRR(2),gqZqLR(2),gqZqRL(2)
      complex(dp):: qgZqLL(2),qgZqRR(2),qgZqLR(2),qgZqRL(2)
      complex(dp):: qbqZgLL(2),qbqZgRR(2),qbqZgLR(2),qbqZgRL(2)
      complex(dp):: gqbZqbLL(2),gqbZqbRR(2),gqbZqbLR(2),gqbZqbRL(2)
      complex(dp):: qbgZqbLL(2),qbgZqbRR(2),qbgZqbLR(2),qbgZqbRL(2)
      complex(dp):: prop,virt5ax,cLL(nf),cLR(nf),cRL(nf),cRR(nf)
      real(dp):: phi,muk,rho,ssig,csig,theta,mtsq,musq,
     & p1true(4),p2true(4),p3true(4),p4true(4),p5true(4)
      logical, parameter:: includeanom=.true.
      integer,parameter::
     & iqqbgLL(5)=(/1,2,3,4,5/),iqqbgRR(5)=(/2,1,4,3,5/),
     & iqqbgRL(5)=(/2,1,3,4,5/),iqqbgLR(5)=(/1,2,4,3,5/),
     & iqgqLL(5)=(/1,5,3,4,2/),iqgqRR(5)=(/5,1,4,3,2/),
     & iqgqRL(5)=(/5,1,3,4,2/),iqgqLR(5)=(/1,5,4,3,2/),
     & igqqLL(5)=(/2,5,3,4,1/),igqqRR(5)=(/5,2,4,3,1/),
     & igqqRL(5)=(/5,2,3,4,1/),igqqLR(5)=(/2,5,4,3,1/)
      common/virt5ax/virt5ax
!$omp threadprivate(/virt5ax/)

      scheme='dred'

c--set msq=0 to initialize
      msq(:,:)=zero
      
c--calculate spinor and dot-products (using BDK type notation)
      call spinoru(5,p,za,zb)

! BEGIN checking code
!        include 'kinpoint5.f'
!        mu=1._dp
!        musq=mu**2
!        mt=0.4255266775_dp
!!        mtsq=mt**2
!        write(6,*) 'mt=',mt
!        do nu=1,4
!        om=nu-1
!        if (nu==1) om=4
!        p(1,om)=p1true(nu)
!        p(2,om)=p2true(nu)
!        p(3,om)=p4true(nu)
!        p(4,om)=p5true(nu)
!        p(5,om)=p3true(nu)
!        enddo
!        do nu=1,5
!        write(6,'(a4,i2,4f18.12)') 'p_',nu,
!     &   p(nu,4),p(nu,1),p(nu,2),p(nu,3)
!        enddo
!        call spinorz(5,p,za,zb)
! END checking code

c--- calculate lowest order
      call qqb_z1jet(p,msq0)
      
c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      subuv=ason2pi*xn
     & *(epinv*(11._dp-2._dp*real(nflav,dp)/xn)-1._dp)/6._dp   ! * zip !!!! DEBUG: compare with Madgraph

c--   calculate propagator
      sz=s(3,4)
      prop=sz/cplx2((sz-zmass**2),zmass*zwidth)

      fac=8._dp*cf*xnsq*esq**2*gsq

!     first letter L/R is the index of the fermion line
!     second letter L/R is the index of the Z/gamma
!     already summed over gluons
      qqbZgLL(1)=aveqq*fac*cplx1(virt5(iqqbgLL,za,zb))
      qqbZgLL(2)=aveqq*fac*virt5ax
      qqbZgLR(1)=aveqq*fac*cplx1(virt5(iqqbgLR,za,zb))
      qqbZgLR(2)=aveqq*fac*virt5ax
      qqbZgRL(1)=aveqq*fac*cplx1(virt5(iqqbgRL,za,zb))
      qqbZgRL(2)=aveqq*fac*virt5ax
      qqbZgRR(1)=aveqq*fac*cplx1(virt5(iqqbgRR,za,zb))
      qqbZgRR(2)=aveqq*fac*virt5ax


      qbqZgLL(:)=qqbZgRL(:)
      qbqZgLR(:)=qqbZgRR(:)
      qbqZgRL(:)=qqbZgLL(:)
      qbqZgRR(:)=qqbZgLR(:)

      gqZqLL(1)=aveqg*fac*cplx1(virt5(igqqLL,za,zb))
      gqZqLL(2)=aveqg*fac*virt5ax
      gqZqLR(1)=aveqg*fac*cplx1(virt5(igqqLR,za,zb))
      gqZqLR(2)=aveqg*fac*virt5ax
      gqZqRL(1)=aveqg*fac*cplx1(virt5(igqqRL,za,zb))
      gqZqRL(2)=aveqg*fac*virt5ax
      gqZqRR(1)=aveqg*fac*cplx1(virt5(igqqRR,za,zb))
      gqZqRR(2)=aveqg*fac*virt5ax

      gqbZqbRL(:)=gqZqLL(:)
      gqbZqbRR(:)=gqZqLR(:)
      gqbZqbLL(:)=gqZqRL(:)
      gqbZqbLR(:)=gqZqRR(:)


      qgZqLL(1)=aveqg*fac*cplx1(virt5(iqgqLL,za,zb))
      qgZqLL(2)=aveqg*fac*virt5ax
      qgZqLR(1)=aveqg*fac*cplx1(virt5(iqgqLR,za,zb))
      qgZqLR(2)=aveqg*fac*virt5ax
      qgZqRL(1)=aveqg*fac*cplx1(virt5(iqgqRL,za,zb))
      qgZqRL(2)=aveqg*fac*virt5ax
      qgZqRR(1)=aveqg*fac*cplx1(virt5(iqgqRR,za,zb))
      qgZqRR(2)=aveqg*fac*virt5ax

      qbgZqbRL(:)=qgZqLL(:)
      qbgZqbRR(:)=qgZqLR(:)
      qbgZqbLL(:)=qgZqRL(:)
      qbgZqbLR(:)=qgZqRR(:)
      
      sin2winv=L(2)-R(2)
      do j=1,nf
      cLL(j)=Q(j)*q1+L(j)*l1*prop
      cRR(j)=Q(j)*q1+R(j)*r1*prop
      cLR(j)=Q(j)*q1+L(j)*r1*prop
      cRL(j)=Q(j)*q1+R(j)*l1*prop
      enddo
      
      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j .ne. 0 .and. k .ne. 0 .and. j .ne. -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
         msq(j,k)=0._dp
      elseif ((j > 0) .and. (k < 0)) then
         msq(j,k)=+real(cLL(j)*conjg(cLL(j))*qqbZgLL(1),dp)
     &            +real(cRR(j)*conjg(cRR(j))*qqbZgRR(1),dp)
     &            +real(cLR(j)*conjg(cLR(j))*qqbZgLR(1),dp)
     &            +real(cRL(j)*conjg(cRL(j))*qqbZgRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(conjg(sin2winv*l1*prop)*(qqbZgLL(2)*cLL(j)+qqbZgRL(2)*cRL(j)),dp)
     &            +real(conjg(sin2winv*r1*prop)*(qqbZgRR(2)*cRR(j)+qqbZgLR(2)*cLR(j)),dp)
         endif
      elseif ((j < 0) .and. (k > 0)) then
         msq(j,k)=+real(cLL(k)*conjg(cLL(k))*qbqZgLL(1),dp)
     &            +real(cRR(k)*conjg(cRR(k))*qbqZgRR(1),dp)
     &            +real(cLR(k)*conjg(cLR(k))*qbqZgLR(1),dp)
     &            +real(cRL(k)*conjg(cRL(k))*qbqZgRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(conjg(sin2winv*l1*prop)*(qbqZgLL(2)*cLL(k)+qbqZgRL(2)*cRL(k)),dp)
     &            +real(conjg(sin2winv*r1*prop)*(qbqZgRR(2)*cRR(k)+qbqZgLR(2)*cLR(k)),dp)
         endif
      elseif ((j > 0) .and. (k == 0)) then
         msq(j,k)=+real(cLL(j)*conjg(cLL(j))*qgZqLL(1),dp)
     &            +real(cRR(j)*conjg(cRR(j))*qgZqRR(1),dp)
     &            +real(cLR(j)*conjg(cLR(j))*qgZqLR(1),dp)
     &            +real(cRL(j)*conjg(cRL(j))*qgZqRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(conjg(sin2winv*l1*prop)*(qgZqLL(2)*cLL(j)+qgZqRL(2)*cRL(j)),dp)
     &            +real(conjg(sin2winv*r1*prop)*(qgZqRR(2)*cRR(j)+qgZqLR(2)*cLR(j)),dp)
         endif
      elseif ((j < 0) .and. (k == 0)) then
         msq(j,k)=+real(cLL(-j)*conjg(cLL(-j))*qbgZqbLL(1),dp)
     &            +real(cRR(-j)*conjg(cRR(-j))*qbgZqbRR(1),dp)
     &            +real(cLR(-j)*conjg(cLR(-j))*qbgZqbLR(1),dp)
     &            +real(cRL(-j)*conjg(cRL(-j))*qbgZqbRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(conjg(sin2winv*l1*prop)*(qbgZqbLL(2)*cLL(-j)+qbgZqbRL(2)*cRL(-j)),dp)
     &            +real(conjg(sin2winv*r1*prop)*(qbgZqbRR(2)*cRR(-j)+qbgZqbLR(2)*cLR(-j)),dp)
         endif
      elseif ((j == 0) .and. (k > 0)) then
         msq(j,k)=+real(cLL(k)*conjg(cLL(k))*gqZqLL(1),dp)
     &            +real(cRR(k)*conjg(cRR(k))*gqZqRR(1),dp)
     &            +real(cLR(k)*conjg(cLR(k))*gqZqLR(1),dp)
     &            +real(cRL(k)*conjg(cRL(k))*gqZqRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(conjg(sin2winv*l1*prop)*(gqZqLL(2)*cLL(k)+gqZqRL(2)*cRL(k)),dp)
     &            +real(conjg(sin2winv*r1*prop)*(gqZqRR(2)*cRR(k)+gqZqLR(2)*cLR(k)),dp)
         endif
      elseif ((j == 0) .and. (k < 0)) then
         msq(j,k)=+real(cLL(-k)*conjg(cLL(-k))*gqbZqbLL(1),dp)
     &            +real(cRR(-k)*conjg(cRR(-k))*gqbZqbRR(1),dp)
     &            +real(cLR(-k)*conjg(cLR(-k))*gqbZqbLR(1),dp)
     &            +real(cRL(-k)*conjg(cRL(-k))*gqbZqbRL(1),dp)
     &            -subuv*msq0(j,k)
         if (includeanom) then
         msq(j,k)=msq(j,k)
     &            +real(conjg(sin2winv*l1*prop)*(gqbZqbLL(2)*cLL(-k)+gqbZqbRL(2)*cRL(-k)),dp)
     &            +real(conjg(sin2winv*r1*prop)*(gqbZqbRR(2)*cRR(-k)+gqbZqbLR(2)*cLR(-k)),dp)
         endif
      endif
      
   19 continue
      enddo
      enddo

      return
      end
