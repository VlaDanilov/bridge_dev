      subroutine qqbggHamp(j1,j2,j3,j4,flip,p,mtsq,Amp)
      implicit none
!--   This routine calculates the amplitude
!--   for the process q^+(1) qb^-(2) g(3) g(4) H(p1234)
!--   with all momenta outgoing.
!--   The amplitude returned is the T^A T^B colour ordering
!--   with the overall factor i g^4 m^2/(4 pi^2 v) removed
!--   Normalization such that Tr {T^A T^B}=1/2*delta^{AB}
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'first.f'
      integer,parameter::a=1,b=2,c=3
      integer::i1,i2,i3,i4,j1,j2,j3,j4,jAB,j12
      real(dp):: p(mxpart,4),p12(4),p34(4),p123(4),p124(4),p3(4),p4(4),
     & s12,s13,s14,s23,s24,s34,s123,s124,mtsq
      complex(dp)::Amp(2,2,2,2),zab2,A21,A12,BA(2,2),tamp(2,2,2,2),
     & B3x4x12(3),B4x3x12(3),
     & B12x3x4(3),B12x4x3(3),
     & B3x12x4(3),B4x12x3(3),
     & FT3x124,FT4x123,FT12x34
      logical:: flip
!--statement function
      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
!--end statement function


      if (first) then
        call pvsetmudim(scale)
        call TRsetmaxindex(3,4,0)
        call pvArraysetup
        first=.false.
      endif

!-- this can be removed when call is through qqbggHampsq
!      include 'specialmom.f'
      call spinoru(4,p,za,zb)
      do j12=1,2
      if (j12 == 2) then
      i1=j1
      i2=j2
      else
      i1=j2
      i2=j1
      endif

      do jAB=1,2
      if (jAB == 1) then
      i3=j3
      i4=j4
      else
      i3=j4
      i4=j3
      endif
      p3(:)=p(i3,:)
      p4(:)=p(i4,:)
      p12(:)=p(i1,:)+p(i2,:)
      p123(:)=p(i1,:)+p(i2,:)+p(i3,:)
      p124(:)=p(i1,:)+p(i2,:)+p(i4,:)
      p34(:)=p(i3,:)+p(i4,:)

      s12=s(i1,i2)
      s13=s(i1,i3)
      s14=s(i1,i4)
      s23=s(i2,i3)
      s24=s(i2,i4)
      s34=s(i3,i4)
      s123=s12+s13+s23
      s124=s12+s14+s24

      call DelDBfunctions(p4,p3,p12,mtsq,B4x3x12)
      call DelDBfunctions(p3,p4,p12,mtsq,B3x4x12)
      call DelDBfunctions(p12,p3,p4,mtsq,B12x3x4)
      call DelDBfunctions(p12,p4,p3,mtsq,B12x4x3)
      call DelDBfunctions(p3,p12,p4,mtsq,B3x12x4)
      call DelDBfunctions(p4,p12,p3,mtsq,B4x12x3)
      call DelDTfunctions(p3,p124,mtsq,FT3x124)
      call DelDTfunctions(p4,p123,mtsq,FT4x123)
      call DelDTfunctions(p12,p34,mtsq,FT12x34)

      if ((jab == 1) .and. (j12 == 2)) then
! Explicit check of other color order
      BA(2,2)=
     &  -zab2(i2,i3,i4,i1)**2/za(i3,i4)
     & /(zb(i1,i2)*za(i2,i3)*za(i4,i2))*FT12x34
     &  +zab2(i2,i1,i4,i3)**2*zb(i1,i4)
     & /(za(i2,i4)*s124)*(one/s12+one/s14)*FT3x124
     &  -zab2(i2,i1,i3,i4)**2*zb(i1,i3)/(za(i2,i3)*s123*s12)*FT4x123
     & -zb(i1,i3)**2/(two*zb(i1,i2)*za(i2,i4))
     & *(B3x4x12(b)*za(i1,i2)*zb(i1,i4)-B12x4x3(b)*za(i2,i3)*zb(i4,i3))
     & +zb(i1,i4)**2/(two*zb(i1,i2)*za(i2,i3))
     & *(B4x3x12(b)*za(i1,i2)*zb(i1,i3)-B12x3x4(b)*za(i2,i4)*zb(i3,i4))
     & -zb(i3,i4)**2/(two*s12)
     & *(B4x12x3(b)*za(i2,i3)*zb(i1,i3)-B3x12x4(b)*za(i2,i4)*zb(i1,i4))
     & -zb(i1,i3)*zb(i1,i4)*zb(i3,i4)
     & /(two*zb(i1,i2))*(B3x12x4(c)+B12x3x4(c))

      BA(1,1)=
     & -zab2(i2,i3,i4,i1)**2
     & /(za(i1,i2)*zb(i1,i4)*zb(i3,i4)*zb(i1,i3))*FT12x34
     & -zab2(i4,i2,i3,i1)**2*za(i2,i3)
     & /(zb(i1,i3)*s123)*(one/s23+one/s12)*FT4x123
     & +zab2(i3,i2,i4,i1)**2*za(i2,i4)/(zb(i1,i4)*s124*s12)*FT3x124
     & -za(i2,i3)**2/(two*za(i1,i2)*zb(i1,i4))
     & *(B3x4x12(b)*za(i2,i4)*zb(i1,i2)-B12x4x3(b)*za(i3,i4)*zb(i1,i3))
     & +za(i2,i4)**2/(two*za(i1,i2)*zb(i1,i3))
     & *(B4x3x12(b)*za(i2,i3)*zb(i1,i2)+B12x3x4(b)*za(i3,i4)*zb(i1,i4))
     & -za(i3,i4)**2/(two*s12)
     & *(B4x12x3(b)*za(i2,i3)*zb(i1,i3)-B3x12x4(b)*za(i2,i4)*zb(i1,i4))
     & +za(i2,i3)*za(i2,i4)*za(i3,i4)
     & /(two*za(i1,i2))*(B3x12x4(c)+B12x3x4(c))

      BA(1,2)=
     & -zb(i2,i4)*zb(i1,i4)**2/(zb(i1,i2)*zb(i2,i3)*zb(i3,i4))*FT4x123
     & -za(i3,i1)*za(i2,i3)**2/(za(i3,i4)*za(i1,i4)*za(i1,i2))*FT3x124
     & -zab2(i3,i1,i2,i4)**2/(two*s34*s12)
     & *(za(i2,i3)*zb(i1,i3)*B3x4x12(b)-za(i2,i4)*zb(i1,i4)*B4x3x12(b))
     & -za(i2,i3)*zb(i1,i4)*zab2(i3,i1,i2,i4)
     & /(s12*s34)*(B12x3x4(a)-B12x4x3(a))

      BA(2,1)=
     & -za(i4,i2)*za(i2,i4)**2/(za(i3,i4)*za(i1,i2)*za(i2,i3))*FT4x123
     & -zb(i1,i3)*zb(i1,i3)**2/(zb(i3,i4)*zb(i1,i2)*zb(i1,i4))*FT3x124
     & -zab2(i4,i1,i2,i3)**2/(two*s34*s12)
     & *(za(i2,i3)*zb(i1,i3)*B3x4x12(b)-za(i2,i4)*zb(i1,i4)*B4x3x12(b))
     & -za(i2,i4)*zb(i1,i3)*zab2(i4,i1,i2,i3)
     & /s34/s12*(B12x3x4(a)-B12x4x3(a))
      endif

!      Bmp(2,2)=
!     & -FT3x124*zab2(i2,i1,i4,i3)**2*zb(i1,i4)
!     & /(za(i2,i4)*s124*s12)
!     & +FT4x123*zab2(i2,i1,i3,i4)**2*zb(i1,i3)
!     * /(za(i2,i3)*s123)*(one/s12+one/s13)
!     & -FT12x34*zab2(i2,i3,i4,i1)**2
!     & /(za(i3,i4)*zb(i1,i2)*za(i2,i3)*za(i2,i4))
!     & +half/(zb(i1,i2)*za(i2,i4)*za(i2,i3))
!     & *(za(i2,i3)*zb(i1,i3)**2
!     & *(B3x4x12(b)*za(i1,i2)*zb(i1,i4)-B12x4x3(b)*za(i2,i3)*zb(i4,i3))
!     & -za(i2,i4)*zb(i1,i4)**2
!     & *(B4x3x12(b)*za(i1,i2)*zb(i1,i3)-B12x3x4(b)*za(i2,i4)*zb(i3,i4)))
!     & + half*zb(i3,i4)**2/s12
!     & *(B4x12x3(b)*za(i2,i3)*zb(i1,i3)-B3x12x4(b)*za(i2,i4)*zb(i1,i4))
!     & +half*(B3x12x4(c)+B12x3x4(c))*zb(i1,i3)*zb(i1,i4)*zb(i3,i4)
!     & /zb(i1,i2)

!      Bmp(1,1)=
!     & +FT4x123*zab2(i4,i2,i3,i1)**2*za(i2,i3)
!     & /(zb(i1,i3)*s12*s123)
!     & -FT3x124*zab2(i3,i2,i4,i1)**2*za(i2,i4)/(zb(i1,i4)*s124)
!     & *(one/s24+one/s12)
!     & +FT12x34*zab2(i2,i3,i4,i1)**2
!     & /(za(i1,i2)*zb(i1,i4)*zb(i3,i4)*zb(i1,i3))
!     & +half*za(i2,i3)**2/(za(i1,i2)*zb(i1,i4))
!     & *(B3x4x12(b)*za(i2,i4)*zb(i1,i2)-B12x4x3(b)*za(i3,i4)*zb(i1,i3))
!     & -half*za(i2,i4)**2/(za(i1,i2)*zb(i1,i3))
!     & *(B4x3x12(b)*za(i2,i3)*zb(i1,i2)+B12x3x4(b)*za(i3,i4)*zb(i1,i4))
!     & +half*za(i3,i4)**2/s12
!     & *(B4x12x3(b)*za(i2,i3)*zb(i1,i3)-B3x12x4(b)*za(i2,i4)*zb(i1,i4))
!     & -half*(B3x12x4(c)+B12x3x4(c))*za(i2,i3)*za(i2,i4)*za(i3,i4)
!     & /za(i1,i2)

!      Bmp(2,1)=
!     & +FT3x124*zb(i2,i3)*zb(i1,i3)**2/(zb(i1,i2)*zb(i2,i4)*zb(i3,i4))
!     & +FT4x123*za(i4,i1)*za(i2,i4)**2/(za(i3,i4)*za(i1,i3)*za(i1,i2))
!     & +half*B3x4x12(b)/(s12*s34)
!     & *za(i2,i3)*zb(i1,i3)*zab2(i4,i1,i2,i3)**2
!     & -half*B4x3x12(b)/(s12*s34)
!     & *za(i2,i4)*zb(i1,i4)*zab2(i4,i1,i2,i3)**2
!     & +(B12x3x4(a)-B12x4x3(a))
!     & /(s12*s34)*za(i2,i4)*zb(i1,i3)*zab2(i4,i1,i2,i3)

!      Bmp(1,2)=
!     & +FT3x124*za(i3,i2)*za(i2,i3)**2/(za(i3,i4)*za(i1,i2)*za(i2,i4))
!     & +FT4x123*zb(i1,i4)*zb(i1,i4)**2/(zb(i3,i4)*zb(i1,i2)*zb(i1,i3))
!     & +half*B3x4x12(b)/(s12*s34)
!     & *za(i2,i3)*zb(i1,i3)*zab2(i3,i1,i2,i4)**2
!     & -half*B4x3x12(b)/(s12*s34)
!     & *za(i2,i4)*zb(i1,i4)*zab2(i3,i1,i2,i4)**2
!     & +(B12x3x4(a)-B12x4x3(a))/(s12*s34)
!     & *za(i2,i3)*zb(i1,i4)*zab2(i3,i1,i2,i4)
!       
!      write(6,*) 'old:AB(2,2)',AB(2,2)
!      write(6,*) 'old:AB(1,1)',AB(1,1)
!      write(6,*) 'old:AB(1,2)',AB(1,2)
!      write(6,*) 'old:AB(2,1)',AB(2,1)

       Amp(jAB,j12,2,2)=
     &  zab2(i2,i3,i4,i1)**2
     & /(za(i3,i4)*zb(i1,i2)*za(i2,i3)*za(i4,i2))*FT12x34
     & -zab2(i2,i1,i4,i3)**2*zb(i1,i4)
     & /(za(i2,i4)*s124*s12)*FT3x124
     & +zab2(i2,i1,i3,i4)**2*zb(i1,i3)
     & /(za(i2,i3)*s123)*(one/s12+one/s13)*FT4x123
     & +zb(i1,i3)**2/(two*zb(i1,i2)*za(i2,i4))
     & *(B3x4x12(b)*za(i1,i2)*zb(i1,i4)-B12x4x3(b)*za(i2,i3)*zb(i4,i3))
     & -zb(i1,i4)**2/(two*zb(i1,i2)*za(i2,i3))
     & *(B4x3x12(b)*za(i1,i2)*zb(i1,i3)-B12x3x4(b)*za(i2,i4)*zb(i3,i4))
     & +zb(i3,i4)**2/(two*s12)
     & *(B4x12x3(b)*za(i2,i3)*zb(i1,i3)-B3x12x4(b)*za(i2,i4)*zb(i1,i4))
     & +zb(i1,i3)*zb(i1,i4)*zb(i3,i4)
     & /(two*zb(i1,i2))*(B3x12x4(c)+B12x3x4(c))

      Amp(jAB,j12,1,1)=
     &  zab2(i2,i3,i4,i1)**2
     & /(za(i1,i2)*zb(i1,i4)*zb(i3,i4)*zb(i1,i3))*FT12x34
     & +zab2(i4,i2,i3,i1)**2*za(i2,i3)/zb(i1,i3)/s12/s123*FT4x123
     & -zab2(i3,i2,i4,i1)**2*za(i2,i4)
     & /zb(i1,i4)/s124*(one/s24+one/s12)*FT3x124
     & +za(i2,i3)**2/(two*za(i1,i2)*zb(i1,i4))
     & *(B3x4x12(b)*za(i2,i4)*zb(i1,i2)-B12x4x3(b)*za(i3,i4)*zb(i1,i3))
     & -za(i2,i4)**2/(two*za(i1,i2)*zb(i1,i3))
     & *(B4x3x12(b)*za(i2,i3)*zb(i1,i2)+B12x3x4(b)*za(i3,i4)*zb(i1,i4))
     & +za(i3,i4)**2/(two*s12)
     & *(B4x12x3(b)*za(i2,i3)*zb(i1,i3)-B3x12x4(b)*za(i2,i4)*zb(i1,i4))
     & -za(i2,i3)*za(i2,i4)*za(i3,i4)
     & /(two*za(i1,i2))*(B3x12x4(c)+B12x3x4(c))

!      Amp(jAB,j12,2,1)=
      A21=
     &  zb(i2,i3)*zb(i1,i3)**2/(zb(i1,i2)*zb(i2,i4)*zb(i3,i4))*FT3x124
     & +za(i4,i1)*za(i2,i4)**2/(za(i3,i4)*za(i1,i3)*za(i1,i2))*FT4x123
     & +zab2(i4,i1,i2,i3)**2/(two*s34*s12)
     & *(za(i2,i3)*zb(i1,i3)*B3x4x12(b)-za(i2,i4)*zb(i1,i4)*B4x3x12(b))
     & +za(i2,i4)*zb(i1,i3)*zab2(i4,i1,i2,i3)
     & /(s12*s34)*(B12x3x4(a)-B12x4x3(a))

!      Amp(jAB,j12,1,2)=
      A12=
     &  za(i3,i2)*za(i2,i3)**2/(za(i3,i4)*za(i1,i2)*za(i2,i4))*FT3x124
     & +zb(i1,i4)*zb(i1,i4)**2/(zb(i3,i4)*zb(i1,i2)*zb(i1,i3))*FT4x123
     & +zab2(i3,i1,i2,i4)**2/(two*s34*s12)
     & *(za(i2,i3)*zb(i1,i3)*B3x4x12(b)-za(i2,i4)*zb(i1,i4)*B4x3x12(b))
     & +za(i2,i3)*zb(i1,i4)*zab2(i3,i1,i2,i4)
     & /s34/s12*(B12x3x4(a)-B12x4x3(a))
      if (jab == 1) then
       Amp(jAB,j12,2,1)=A21
       Amp(jAB,j12,1,2)=A12
      else
       Amp(jAB,j12,1,2)=A21
       Amp(jAB,j12,2,1)=A12
      endif
! The counter-intuitive labelling below is what the old code uses
      if (jab == 1) then
       Amp(jAB,j12,2,1)=A12
       Amp(jAB,j12,1,2)=A21
      else
       Amp(jAB,j12,1,2)=A12
       Amp(jAB,j12,2,1)=A21
      endif
      enddo
      enddo
! Other color ordering

! Fix flippery involved with computing j12=1
      tAmp(:,:,:,:)=Amp(:,:,:,:)
      do i1=1,2
      Amp(i1,1,:,:)=tAmp(3-i1,1,:,:)
      enddo
      
! Normalization for comparison with old code
! eventually to be removed
! debug
      Amp(:,:,:,:)=-im/two*Amp(:,:,:,:)
      BA(:,:)=-im/two*BA(:,:)
! debug
      
!      write(6,*) 'AB(1,1)',Amp(1,2,1,1)
!      write(6,*) 'AB(1,2)',Amp(1,2,1,2)
!      write(6,*) 'AB(2,1)',Amp(1,2,2,1)
!      write(6,*) 'AB(2,2)',Amp(1,2,2,2)
!      write(6,*)
!      write(6,*) 'BA(1,1)',Amp(2,2,1,1)
!      write(6,*) 'BA(1,2)',Amp(2,2,1,2)
!      write(6,*) 'BA(2,1)',Amp(2,2,2,1)
!      write(6,*) 'BA(2,2)',Amp(2,2,2,2)
c      write(6,*)
c      write(6,*) 'BA(1,1)',BA(1,1)
c      write(6,*) 'BA(2,1)',BA(2,1)
c      write(6,*) 'BA(1,2)',BA(1,2)
c      write(6,*) 'BA(2,2)',BA(2,2)

!      if (flip) then
!      BA(1,1)=AB(1,1)
!      BA(2,2)=AB(2,2)
!      BA(1,2)=AB(2,1)
!      BA(2,1)=AB(1,2)
!      AB(:,:)=-BA(:,:)
!      endif

      return
      end
