      subroutine qqbZgg_floop(p,mtin,ampTree,ampVec,ampAx)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'qqbggintnames.f'
c--- Returns a series of arrays representing the amp[itudes
c--- for vector and axial couplings of internal particles

c--- for the process 0 -> qqb Z gg
c--- The overall factor on the amplitude is:
c---
c---      4._dp*esq*gsq/(16._dp*pisq)*esq * delta(a,b)
c---
c--- apportioned as follows:
c--- 8*e^2*gsq from gg->ZZ, extra esq from Z decays, 1/(16*pisq) from
c--- QCDLoop normalization and 1/2 from color factor
c---
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
!      include 'masses.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scale.f'
      include 'docheck.f'
      complex(dp)::
     & ampTree(2,2,2,2,2),ampVec(2,2,2,2,2),ampAx(0:2,2,2,2,2),amp0(2,2),
     & resKC(0:2,2,2,2,2),resKCvec(2,2,2,2),resKClo(2,2,2,2,2),
     & zavec(mxpart,mxpart),zbvec(mxpart,mxpart),
     & Coeffs0(15,2,2),Coeffs2(15,2,2),Coeffsl(15,2,2),xInt(15),ABDK(2,2)
      real(dp):: s12,s34,s56,s123,s124
      logical:: ggZZuse6d
      integer:: h1,h2,h3,h4,h5,h6,h12,h34,h56,j1,j2,j3,j4,j5,j6,j34,up,dn,om,nu,
     & k34,k3,k4,k
      real(dp):: p(mxpart,4),pvec(mxpart,4)
      real(dp):: phi,muk,rho,ssig,csig,theta,mtin,mt,mtsq,mu,
     & p1true(4),p2true(4),p3true(4),p4true(4),p5true(4),p6true(4),sign
      complex(dp):: AmtLL(2,2,2,2),AmtLR(2,2,2,2),Amp
c     & ,Avec(2,2,2,2),a64v,,Fax,Faxsl,Amp,sum(2,2),App,Apm,Appnew,Apmnew
c     & AmtLL_new(2,2,2,2),AmtLR_new(2,2,2,2)
      character*7 junk
      parameter(up=2,dn=1)
      logical, parameter:: writecoeffs=.false.
      logical, parameter:: useBDK=.true.     ! use A6axBDK or not
      
      docheck=.false.   ! compare with KC at phase space point

!!! WARNING: do not set ggZZuse6d=.true. without first upgrading ZZmassiveboxtri
      ggZZuse6d=.false.
      if (ggZZuse6d) then
        write(6,*) 'Must check gq channels before proceeding'
        stop
      endif
!!! WARNING

! DEBUG
!      mt=0._dp
      mt=mtin
      mtsq=mt**2

c--- if this variable (passed via common block) is true then
c--- set the kinematics to a special point for numerical check
      if (docheck) then
        include 'kinpoint.f'
        mu=1._dp
        musq=mu**2
        mt=0.4255266775_dp
!        mt=1.e6_dp
!        mt=0._dp
        mtsq=mt**2
        write(6,*) 'mt=',mt
        do nu=1,4
        om=nu-1
        if (nu==1) om=4
        p(1,om)=p1true(nu)
        p(2,om)=p2true(nu)
        p(3,om)=p3true(nu)
        p(4,om)=p4true(nu)
        p(5,om)=p5true(nu)
        p(6,om)=p6true(nu)
! for qqbggVec comparison
!        p(1,om)=p3true(nu)
!        p(2,om)=p4true(nu)
!        p(3,om)=p1true(nu)
!        p(4,om)=p2true(nu)
!        p(5,om)=p6true(nu)
!        p(6,om)=p5true(nu)
! end for qqbggVec comparison
        enddo
        do nu=1,6
        write(6,'(a4,i2,4f18.12)') 'p_',nu,
     &   p(nu,4),p(nu,1),p(nu,2),p(nu,3)
        enddo
!        call spinorz(6,p,za,zb)
      endif
c--- fill amplitudes
c--- labels are: helicity of gluons, lepton 3 and lepton 5

c--- compute amplitudes with massless quarks
!      if (docheck) then
!        call spinoru(6,p,za,zb)
!        AmpVec(1,2,2,1)=a64v('q+qb-g-g-',1,2,3,4,5,6,zb,za)
!        AmpVec(1,2,1,1)=a64v('q+qb-g-g+',1,2,3,4,5,6,zb,za)
!        AmpVec(1,1,2,1)=a64v('q+qb-g+g-',1,2,3,4,5,6,zb,za)
!        AmpVec(1,1,1,1)=a64v('q+qb-g+g+',1,2,3,4,5,6,zb,za)
!
!        AmpVec(2,2,2,1)=a64v('q+qb-g+g+',1,2,3,4,6,5,za,zb)
!        AmpVec(2,2,1,1)=a64v('q+qb-g+g-',1,2,3,4,6,5,za,zb)
!        AmpVec(2,1,2,1)=a64v('q+qb-g-g+',1,2,3,4,6,5,za,zb)
!        AmpVec(2,1,1,1)=a64v('q+qb-g-g-',1,2,3,4,6,5,za,zb)
!
!        AmpVec(1,2,2,2)=a64v('q+qb-g-g-',1,2,3,4,6,5,zb,za)
!        AmpVec(1,2,1,2)=a64v('q+qb-g-g+',1,2,3,4,6,5,zb,za)
!        AmpVec(1,1,2,2)=a64v('q+qb-g+g-',1,2,3,4,6,5,zb,za)
!        AmpVec(1,1,1,2)=a64v('q+qb-g+g+',1,2,3,4,6,5,zb,za)
!
!        AmpVec(2,2,2,2)=a64v('q+qb-g+g+',1,2,3,4,5,6,za,zb)
!        AmpVec(2,2,1,2)=a64v('q+qb-g+g-',1,2,3,4,5,6,za,zb)
!        AmpVec(2,1,2,2)=a64v('q+qb-g-g+',1,2,3,4,5,6,za,zb)
!        AmpVec(2,1,1,2)=a64v('q+qb-g-g-',1,2,3,4,5,6,za,zb)
!        write(6,*) 'massless loops: Avec/(-im)'
!        do h1=1,2
!        do h2=1,2
!        do h34=1,2
!        do h56=1,2
!        write(6,*) h1,h2,h34,h56,Avec(h1,h2,h34,h56)/(-im),abs(Avec(h1,h2,h34,h56))
!        enddo
!        enddo
!        enddo
!        enddo
!        write(6,*)
c--- a64v amplitudes have overall color factor delta(A,B);
c--- put in factor of two to extract delta(A,B)/2 as in massive case
!        Avec(:,:,:,:)=two*Avec(:,:,:,:)
!      endif

c--- compute amplitudes with a massive top quark internal loop
      pvec(1,:)=p(3,:)
      pvec(2,:)=p(4,:)
      pvec(3,:)=p(1,:)
      pvec(4,:)=p(2,:)
      pvec(5,:)=p(6,:)
      pvec(6,:)=p(5,:)
c--- set up spinor products
      call spinoru(6,pvec,zavec,zbvec)
      if (docheck) call spinorz(6,pvec,zavec,zbvec)
      
!      if (docheck) write(6,*) '>>>>>> TOP AMPLITUDE <<<<<'
      if (ggZZuse6d) then ! scalar integrals including 6-d box
        call ggZZmassamp_new(pvec,zavec,zbvec,mt,AmtLL,AmtLR)
      else                ! all scalar integrals in 4-d
        call ggZZmassamp(pvec,zavec,zbvec,mt,AmtLL,AmtLR)
      endif
      do h12=1,2
      do h3=1,2
      do h4=1,2
      do h56=1,2
      ampVec(1,h12,h3,h4,h56)=-(AmtLL(h3,h4,h12,h56)+AmtLR(h3,h4,h12,h56))
      enddo
      enddo
      enddo
      enddo
      
      ampVec(2,:,:,:,:)=ampVec(1,:,:,:,:)
      
      if (docheck) then
c--- this block of code is to check A6s primitive
         pvec(1,:)=p(1,:)
         pvec(2,:)=p(3,:)
         pvec(3,:)=p(4,:)
         pvec(4,:)=p(2,:)
         pvec(5,:)=p(5,:)
         pvec(6,:)=p(6,:)
         call spinorz(6,pvec,zavec,zbvec)
         j1=1
         j2=2
         j3=3
         j4=4
         j5=5
         j6=6
         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(one/three*zip+mtsq*(
     &            two*qlI3(s(j2,j3),zip,zip,mtsq,mtsq,mtsq,musq,-1)
     &           +four/s(j2,j3)*(qlI2(s(j2,j3),mtsq,mtsq,musq,-1)
     &                          -qlI2(zip,mtsq,mtsq,musq,-1))))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic (-1)',amp,abs(amp)
         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(one/three+mtsq*(
     &            two*qlI3(s(j2,j3),zip,zip,mtsq,mtsq,mtsq,musq,0)
     &           +four/s(j2,j3)*(qlI2(s(j2,j3),mtsq,mtsq,musq,0)
     &                          -qlI2(zip,mtsq,mtsq,musq,0))))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic ( 0)',amp,abs(amp)
         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(+mtsq*(
     &            two
     &                          ))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic C0(s34) coeff',amp,abs(amp)
         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(mtsq*(
     &           +four/s(j2,j3)
     &                          ))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic B0(s34) coeff',amp,abs(amp)

c--- bubble diagrams only
         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(-two/nine*zip+two/three*(qlI2(zip,mtsq,mtsq,musq,-1))
     &           +(two/three+four/three/s(j2,j3)*mtsq)*(qlI2(s(j2,j3),mtsq,mtsq,musq,-1)
     &                                                 -qlI2(zip,mtsq,mtsq,musq,-1)))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic bubbles (-1)',amp,abs(amp)
         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(-two/nine+two/three*(qlI2(zip,mtsq,mtsq,musq,0))
     &           +(two/three+four/three/s(j2,j3)*mtsq)*(qlI2(s(j2,j3),mtsq,mtsq,musq,0)
     &                                                 -qlI2(zip,mtsq,mtsq,musq,0)))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic bubbles ( 0)',amp,abs(amp)

         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(five/nine*zip
     &           +two*mtsq*qlI3(s(j2,j3),zip,zip,mtsq,mtsq,mtsq,musq,-1)
     &           -two/three*qlI2(zip,mtsq,mtsq,musq,-1)
     &           +(eight/three/s(j2,j3)*mtsq-two/three)
     &               *(qlI2(s(j2,j3),mtsq,mtsq,musq,-1)-qlI2(zip,mtsq,mtsq,musq,-1)))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic triangles (-1)',amp,abs(amp)
         Amp= im/zavec(j2,j3)**2/s(j5,j6)*(
     &     -zavec(j4,j5)*(zbvec(j6,j1)*zavec(j1,j3)+zbvec(j6,j2)*zavec(j2,j3))*zbvec(j3,j1)
     &        /(s(j1,j2)+s(j1,j3)+s(j2,j3))
     &     +zbvec(j1,j6)*(zavec(j5,j4)*zbvec(j4,j3)+zavec(j5,j2)*zbvec(j2,j3))*zavec(j3,j4)
     &        /(s(j2,j3)+s(j2,j4)+s(j3,j4)))
     &        *(five/nine
     &           +two*mtsq*qlI3(s(j2,j3),zip,zip,mtsq,mtsq,mtsq,musq,0)
     &           -two/three*qlI2(zip,mtsq,mtsq,musq,0)
     &           +(eight/three/s(j2,j3)*mtsq-two/three)
     &               *(qlI2(s(j2,j3),mtsq,mtsq,musq,0)-qlI2(zip,mtsq,mtsq,musq,0)))
        Amp=Amp*im ! why this factor?
         write(6,*) 'analytic triangles ( 0)',amp,abs(amp)
!         stop
! end this block of code

!        call ggZZwritetable
!        write(6,*) 'mt=',mt
!        write(6,*)
!        write(6,'(a15,SP,2e20.10)') 'AmtLR(2,2,1,1)',AmtLR(2,2,1,1)
!        write(6,'(a15,SP,2e20.10)') 'AmtLR(1,2,1,1)',AmtLR(1,2,1,1)
!        write(6,*)
!        write(6,'(a15,SP,2e20.10)') 'AmtLL(2,2,1,1)',AmtLL(2,2,1,1)
!        write(6,'(a15,SP,2e20.10)') 'AmtLL(1,2,1,1)',AmtLL(1,2,1,1)
!        write(6,*)
          open(unit=11,file='resvec',status='old')
          do nu=1,16
          read(11,*) junk,h1,h2,h3,h4,h5,h6,Amp
          resKCvec((3+h1)/2,(3+h3)/2,(3+h4)/2,(3+h6)/2)=Amp
          enddo
          close(11)
          do h12=2,1,-1
          do h3=1,2
          do h4=1,2
          do h56=1,2
          write(6,*) 'v :',2*h12-3,3-2*h12,h3,h4,3-2*h56,2*h56-3,
     &     ampVec(1,h12,h3,h4,h56)/resKCvec(h12,h3,h4,h56)
          enddo
          enddo
          enddo
          enddo
!          stop
!          stop
!        write(6,*) 'massive loops: -(AmtLL+AmtLR)'
!        do h1=1,2
!        do h2=1,2
!        do h34=1,2
!        do h56=1,2
!        write(6,*) h1,h2,h34,h56,
!     &   -(AmtLL(h1,h2,h34,h56)+AmtLR(h1,h2,h34,h56)),abs(AmtLL(h1,h2,h34,h56)+AmtLR(h1,h2,h34,h56))
!        enddo
!        enddo
!        enddo
!        enddo
!        stop
      endif

!      call triangleassemble(1,2,3,4,5,6,za,zb,mtsq,Ccoeff)
!      pause
!      call qqbggAxtri12x3x456(1,2,3,4,5,6,1,1,za,zb,mtsq,tri12_3)
!      write(6,*) 'tri12_3(1,1)',tri12_3(1,1)
!      write(6,*) 'tri12_3(1,2)',tri12_3(1,2)
!      write(6,*) 'tri12_3(2,1)',tri12_3(2,1)
!      write(6,*) 'tri12_3(2,2)',tri12_3(2,2)
!      call qqbggAxtri12x3x456(1,2,3,4,5,6,2,2,za,zb,mtsq,tri12_3)
!      write(6,*) 'tri12_3(1,1)',tri12_3(1,1)
!      write(6,*) 'tri12_3(1,2)',tri12_3(1,2)
!      write(6,*) 'tri12_3(2,1)',tri12_3(2,1)
!      write(6,*) 'tri12_3(2,2)',tri12_3(2,2)
!      pause

c--- set up spinor products
      call spinoru(6,p,za,zb)
      if (docheck) call spinorz(6,p,za,zb)

! debug      
!         App=Fax('q+qb-g+g+',1,2,3,4,5,6,za,zb)
!         Apm=Fax('q+qb-g+g-',1,2,3,4,5,6,za,zb)
!         Amp=Fax('q+qb-g-g+',1,2,3,4,5,6,za,zb)
         
!         call BDKqqbggAx(1,2,3,4,5,6,za,zb,sum)
         
!         write(6,*) 'Fax pp',sum(2,2),App,sum(2,2)/App
!         write(6,*) 'Fax pm',sum(2,1),Apm,sum(2,1)/Apm
!         write(6,*) 'Fax mp',sum(1,2),Amp,sum(1,2)/Amp
!         write(6,*)

!         App=Faxsl('q+qb-g+g+',1,2,3,4,5,6,za,zb)
!         Apm=Faxsl('q+qb-g+g-',1,2,3,4,5,6,za,zb) 
!         call BDKqqbggAxsla(1,2,3,4,5,6,za,zb,Appnew,Apmnew)
!         write(6,*) 'Faxsl pp',Appnew,App,Appnew/App
!         write(6,*) 'Faxsl pm',Apmnew,Apm,Apmnew/Apm
!         pause
! debug
         
! loop over color orderings (j34=1 T34mdel, j34=2 T43mdel)
      do j34=1,2
      if (j34 == 1) then
        j3=3
        j4=4
      else
        j3=4
        j4=3
      endif

c--- compute basis integrals
      s12=s(1,2)
      s34=s(j3,j4)
      s56=s(5,6)
      s123=s(1,2)+s(1,j3)+s(2,j3)
      s124=s(1,2)+s(1,j4)+s(2,j4)
      xInt(d3_12_4)=qlI4(zip,s56,zip,s12,s124,s123,mtsq,mtsq,mtsq,mtsq,musq,0)
      xInt(d4_3_12)=qlI4(s56,zip,zip,s12,s123,s34,mtsq,mtsq,mtsq,mtsq,musq,0)
      xInt(d3_4_12)=qlI4(s56,zip,zip,s12,s124,s34,mtsq,mtsq,mtsq,mtsq,musq,0)
      xInt(c12_34)=qlI3(s56,s34,s12,mtsq,mtsq,mtsq,musq,0)
      xInt(c12_3)=qlI3(s123,s12,zip,mtsq,mtsq,mtsq,musq,0)
      xInt(c12_4)=qlI3(s124,s12,zip,mtsq,mtsq,mtsq,musq,0)
      xInt(c3_124)=qlI3(s56,s124,zip,mtsq,mtsq,mtsq,musq,0)
      xInt(c4_123)=qlI3(s56,s123,zip,mtsq,mtsq,mtsq,musq,0)
      xInt(c3_4)=qlI3(s34,zip,zip,mtsq,mtsq,mtsq,musq,0)
      xInt(b12)=qlI2(s12,mtsq,mtsq,musq,0)
      xInt(b34)=qlI2(s34,mtsq,mtsq,musq,0)
      xInt(b56)=qlI2(s56,mtsq,mtsq,musq,0)
      xInt(b123)=qlI2(s123,mtsq,mtsq,musq,0)
      xInt(b124)=qlI2(s124,mtsq,mtsq,musq,0)
      xInt(rat)=cone

c--- default ordering corresponds to (outgoing) q(+) qb(-) e+(-) e+(-) e-(+)
c--- h12 , h56 thus refer to helicities of quark(1), electron(6) respectively
      do h12=1,2
      do h56=1,2
      if (h12 == 2) then
        j1=1
        j2=2
      else
        j1=2
        j2=1
      endif
      if (h56 == 2) then
        j5=5
        j6=6
      else
        j5=6
        j6=5
      endif

      if (useBDK) then
        call BDKqqbggAxAmp(j1,j2,j3,j4,j5,j6,za,zb,mtsq,xInt,ABDK)
      else
        call BDKqqbggAxCoeffs(j1,j2,j3,j4,j5,j6,za,zb,mtsq,Coeffs0,Coeffs2)
      endif
      
      if (j34 == 1) then
        call qqbggAxslCoeffs(j1,j2,j3,j4,j5,j6,za,zb,mtsq,Coeffsl)
      endif
     
      do h3=2,1,-1
      do h4=2,1,-1
!      if (h12 == 3) then
!        ampAx(j34,h12,h3,h4,3-h56)=
!     &   +conjg(ebox0(h3,h4)+mtsq*ebox2(h3,h4))*D0(d3_12_4)
!     &   +conjg(h3box0(h3,h4)+mtsq*h3box2(h3,h4))*D0(d4_3_12)
!     &   +conjg(h4box0(h3,h4)+mtsq*h4box2(h3,h4))*D0(d3_4_12)
!     &   +conjg(tri12_34_0(h3,h4)+mtsq*tri12_34_2(h3,h4))*C0(c12_34)
!     &   +conjg(tri12_4(h3,h4))*C0(c12_4)
!     &   +conjg(tri4_123_0(h3,h4)+mtsq*tri4_123_2(h3,h4))*C0(c4_123)
!     &   +conjg(tri3_124_0(h3,h4)+mtsq*tri3_124_2(h3,h4))*C0(c3_124)
!     &   +conjg(tri12_3(h3,h4))*C0(c12_3)
!     &   +conjg(tri3_4(h3,h4))*C0(c3_4)
!     &   +conjg(bub12(h3,h4))*B0(b12)
!     &   +conjg(bub34(h3,h4))*B0(b34)
!     &   +conjg(bub56(h3,h4))*B0(b56)
!     &   +conjg(bub123(h3,h4))*B0(b123)
!     &   +conjg(bub124(h3,h4))*B0(b124)
!     &   +conjg(rat(h3,h4))
!      else

! interchange labels appropriately: j34 -> k34, h3 -> k3, h4 -> k4
        if (h12 == 1) then
           k34=3-j34
           sign=-one
        else
          k34=j34
          k3=h3
          k4=h4
          sign=one
        endif
        if (j34 == 2) then
           k3=h4
           k4=h3
        else
           k3=h3
           k4=h4
        endif
        
        if (useBDK) then 
          ampAx(k34,h12,k3,k4,h56)=ABDK(h3,h4)
        else
          ampAx(k34,h12,k3,k4,h56)=czip
          do k=1,15
          ampAx(k34,h12,k3,k4,h56)=ampAx(k34,h12,k3,k4,h56)
     &     +xInt(k)*(Coeffs0(k,h3,h4)+mtsq*Coeffs2(k,h3,h4))
          enddo
        endif
        
        ampAx(k34,h12,k3,k4,h56)=ampAx(k34,h12,k3,k4,h56)*sign

!        ampAx(k34,h12,k3,k4,h56)=sign*(
!     &   +(ebox0(h3,h4)+mtsq*ebox2(h3,h4))*D0(d3_12_4)
!     &   +(h3box0(h3,h4)+mtsq*h3box2(h3,h4))*D0(d4_3_12)
!     &   +(h4box0(h3,h4)+mtsq*h4box2(h3,h4))*D0(d3_4_12)
!     &   +(tri12_34_0(h3,h4)+mtsq*tri12_34_2(h3,h4))*C0(c12_34)
!     &   +(tri12_4(h3,h4))*C0(c12_4)
!     &   +(tri4_123_0(h3,h4)+mtsq*tri4_123_2(h3,h4))*C0(c4_123)
!     &   +(tri3_124_0(h3,h4)+mtsq*tri3_124_2(h3,h4))*C0(c3_124)
!     &   +(tri12_3(h3,h4))*C0(c12_3)
!     &   +(tri3_4(h3,h4))*C0(c3_4)
!     &   +bub12(h3,h4)*B0(b12)
!     &   +bub34(h3,h4)*B0(b34)
!     &   +bub56(h3,h4)*B0(b56)
!     &   +bub123(h3,h4)*B0(b123)
!     &   +bub124(h3,h4)*B0(b124)
!     &   +rat(h3,h4))
     
        if (j34 == 1) then
        ampAx(0,h12,k3,k4,h56)=sign*(
     &   +coeffsl(c3_124,h3,h4)*xInt(c3_124)
     &   +coeffsl(c4_123,h3,h4)*xInt(c4_123)
     &   +coeffsl(b56,h3,h4)*xInt(b56)
     &   +coeffsl(b123,h3,h4)*xInt(b123)
     &   +coeffsl(b124,h3,h4)*xInt(b124)
     &   +coeffsl(rat,h3,h4))
        endif

        call qqbZggtree(j1,j3,j4,j2,j5,j6,za,zb,amp0)
        ampTree(k34,h12,k3,k4,h56)=amp0(h3,h4)
        
!        write(6,*) 'Ccoeffsl(c3_124,2,2)',Ccoeffsl(c3_124,2,2)
!        write(6,*) 'Ccoeffsl(c4_123,2,2)',Ccoeffsl(c4_123,2,2)
!        write(6,*) 'Bcoeffsl(b56,2,2)',Bcoeffsl(b56,2,2)
!        write(6,*) 'Bcoeffsl(b123,2,2)',Bcoeffsl(b123,2,2)
!        write(6,*) 'Bcoeffsl(b124,2,2)',Bcoeffsl(b124,2,2)
!        write(6,*) 'Ratsl(2,2)',Ratsl(2,2)

        
!      endif
        if ((writecoeffs) .and. (h12 == 2) .and. (h56 == 2)) then
!          write(6,*) 'easy ',h3,h4,Coeffs0(d3_12_4,h3,h4)+mtsq*Coeffs2(d3_12_4,h3,h4)
!          write(6,*) 'h123 ',h3,h4,Coeffs0(d4_3_12,h3,h4)+mtsq*Coeffs2(d4_3_12,h3,h4)
!          write(6,*) 'h124 ',h3,h4,Coeffs0(d3_4_12,h3,h4)+mtsq*Coeffs2(d3_4_12,h3,h4)
!          write(6,*) '12_34',h3,h4,Coeffs0(c12_34,h3,h4)+mtsq*Coeffs2(c12_34,h3,h4)
!          write(6,*) '12_4 ',h3,h4,Coeffs0(c12_4,h3,h4)+mtsq*Coeffs2(c12_4,h3,h4)
!          write(6,*) '123_4',h3,h4,Coeffs0(c4_123,h3,h4)+mtsq*Coeffs2(c4_123,h3,h4)
!          write(6,*) '124_3',h3,h4,Coeffs0(c3_124,h3,h4)+mtsq*Coeffs2(c3_124,h3,h4)
!          write(6,*) '12_3 ',h3,h4,Coeffs0(c12_3,h3,h4)+mtsq*Coeffs2(c12_3,h3,h4)
!          write(6,*) '3_4  ',h3,h4,Coeffs0(c3_4,h3,h4)+mtsq*Coeffs2(c3_4,h3,h4)
!          write(6,*) 'b34  ',h3,h4,Coeffs0(b34,h3,h4)+mtsq*Coeffs2(b34,h3,h4)
!          write(6,*) 'b56  ',h3,h4,Coeffs0(b56,h3,h4)+mtsq*Coeffs2(b56,h3,h4)
!          write(6,*) 'b12  ',h3,h4,Coeffs0(b12,h3,h4)+mtsq*Coeffs2(b12,h3,h4)
!          write(6,*) 'b124 ',h3,h4,Coeffs0(b124,h3,h4)+mtsq*Coeffs2(b124,h3,h4)
!          write(6,*) 'b123 ',h3,h4,Coeffs0(b123,h3,h4)+mtsq*Coeffs2(b123,h3,h4)
!          write(6,*) 'rat  ',h3,h4,Coeffs0(rat,h3,h4)
!          write(6,*)
          write(6,51) 'easy ',h3,h4,Coeffs0(d3_12_4,h3,h4),Coeffs2(d3_12_4,h3,h4)
          write(6,51) 'h123 ',h3,h4,Coeffs0(d4_3_12,h3,h4),Coeffs2(d4_3_12,h3,h4)
          write(6,51) 'h124 ',h3,h4,Coeffs0(d3_4_12,h3,h4),Coeffs2(d3_4_12,h3,h4)
          write(6,*)
          write(6,51) '12_34',h3,h4,Coeffs0(c12_34,h3,h4),Coeffs2(c12_34,h3,h4)
          write(6,51) '12_3 ',h3,h4,Coeffs0(c12_3,h3,h4),Coeffs2(c12_3,h3,h4)
          write(6,51) '12_4 ',h3,h4,Coeffs0(c12_4,h3,h4),Coeffs2(c12_4,h3,h4)
          write(6,51) '124_3',h3,h4,Coeffs0(c3_124,h3,h4),Coeffs2(c3_124,h3,h4)
          write(6,51) '123_4',h3,h4,Coeffs0(c4_123,h3,h4),Coeffs2(c4_123,h3,h4)
          write(6,51) '3_4  ',h3,h4,Coeffs0(c3_4,h3,h4),Coeffs2(c3_4,h3,h4)
          write(6,*)
          write(6,51) 'b12  ',h3,h4,Coeffs0(b12,h3,h4),Coeffs2(b12,h3,h4)
          write(6,51) 'b34  ',h3,h4,Coeffs0(b34,h3,h4),Coeffs2(b34,h3,h4)
          write(6,51) 'b56  ',h3,h4,Coeffs0(b56,h3,h4),Coeffs2(b56,h3,h4)
          write(6,51) 'b123 ',h3,h4,Coeffs0(b123,h3,h4),Coeffs2(b123,h3,h4)
          write(6,51) 'b124 ',h3,h4,Coeffs0(b124,h3,h4),Coeffs2(b124,h3,h4)
          write(6,51) 'rat  ',h3,h4,Coeffs0(rat,h3,h4)
          write(6,51)
        endif
      enddo
      enddo

      enddo
      enddo
      
      if (writecoeffs) stop
   51 format(a8,2i5,4(f16.10,'   & '))
      enddo

c--- write-out comparison with KC results if required
      if (docheck) then
        open(unit=11,file='res34',status='old')
        open(unit=12,file='res43',status='old')
        open(unit=13,file='resdel',status='old')
        do nu=1,16
        read(11,*) junk,h1,h2,h3,h4,h5,h6,Amp
        resKC(1,(3+h1)/2,(3+h3)/2,(3+h4)/2,(3+h6)/2)=Amp
        read(12,*) junk,h1,h2,h3,h4,h5,h6,Amp
        resKC(2,(3+h1)/2,(3+h3)/2,(3+h4)/2,(3+h6)/2)=Amp
        read(13,*) junk,h1,h2,h3,h4,h5,h6,Amp
        resKC(0,(3+h1)/2,(3+h3)/2,(3+h4)/2,(3+h6)/2)=Amp
        enddo
        close(11)
        close(12)
        close(13)
        do j34=2,0,-1
        write(6,*)
        do h12=2,1,-1
        do h3=1,2
        do h4=1,2
        do h56=1,2
        write(6,*) j34,':',2*h12-3,3-2*h12,h3,h4,3-2*h56,2*h56-3,
     &   ampAx(j34,h12,h3,h4,h56)/resKC(j34,h12,h3,h4,h56)
        enddo
        enddo
        enddo
        enddo
        enddo
        open(unit=11,file='reslo34',status='old')
        open(unit=12,file='reslo43',status='old')
        do nu=1,16
        read(11,*) junk,h1,h2,h3,h4,h5,h6,Amp
        resKClo(1,(3+h1)/2,(3+h3)/2,(3+h4)/2,(3+h6)/2)=Amp
        read(12,*) junk,h1,h2,h3,h4,h5,h6,Amp
        resKClo(2,(3+h1)/2,(3+h3)/2,(3+h4)/2,(3+h6)/2)=Amp
        enddo
        close(11)
        close(12)
        do j34=1,2
        write(6,*)
        do h12=2,1,-1
        do h3=1,2
        do h4=1,2
        do h56=1,2
c--- recall that ampTree = -i*A6tree
!        write(6,*) j34,' lo :',2*h12-3,3-2*h12,h3,h4,3-2*h56,2*h56-3,
!     &   ampTree(j34,h12,h3,h4,h56)/(-im*resKClo(j34,h12,h3,h4,h56))
        enddo
        enddo
        enddo
        enddo
        enddo
        stop
      endif
      
      return
      end
