      subroutine BDKqqbggAxAmp(p1,p2,p3,p4,p5,p6,za,zb,mtsq,xInt,ABDK)
      implicit none
C                               |6
C     3                         ^
C      ccccccc---->-----cccccccc|
C            |          |       |5
C            |          |
C            |          |
C     1|     |          |
C      |     |          |
C      |cccccc-----<-----cccccccc 4
C      |
C     2|


      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'color34.f'
      include 'qqbggintnames.f'
      integer p1,p2,p3,p4,p5,p6,h3,h4,k
      complex(dp)::zab2,izab2,zba2,izba2,iza,izb,
     & coeffs0(15,2,2),coeffs2(15,2,2),
     & box0(2,2),box2(2,2),tri0(2,2),tri2(2,2),xInt(15),xIntflip2(15),ABDK(2,2),
     & A6axBDKpp,A6axBDKpm,A6axBDKmp
      real(dp)::s3,s12,s34,s56,s123,s124,mtsq,s123s,s124s,
     & is,is3,del3,idel3,d12,d34,d56
      logical, parameter:: writecoeffs=.false.
      
C--- begin statement function
      s3(p1,p2,p3)=s(p1,p2)+s(p2,p3)+s(p3,p1)
      is3(p1,p2,p3)=1d0/s3(p1,p2,p3)
      is(p1,p2)=1d0/s(p1,p2)
      zab2(p1,p2,p3,p4)=za(p1,p2)*zb(p2,p4)+za(p1,p3)*zb(p3,p4)
      zba2(p1,p2,p3,p4)=zb(p1,p2)*za(p2,p4)+zb(p1,p3)*za(p3,p4)
      izab2(p1,p2,p3,p4)=cone/zab2(p1,p2,p3,p4)
      izba2(p1,p2,p3,p4)=cone/zba2(p1,p2,p3,p4)
      izb(p1,p2)=cone/zb(p1,p2)
      iza(p1,p2)=cone/za(p1,p2)
C--- end statement functions
      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s123=s3(p1,p2,p3)
      s124=s3(p1,p2,p4)
      s124s=s124-s12
      s123s=s123-s12
      d12=s12-s34-s56
      d34=s34-s56-s12
      d56=s56-s12-s34
      del3=s12**2+s34**2+s56**2
     & -2*s12*s34-2*s34*s56-2*s56*s12
      idel3=1d0/del3

      coeffs0(:,:,:)=czip
      coeffs2(:,:,:)=czip

c--- natural ordering of hardbox gives (s12,s56,s124,s34) -> h4box
      call qqbggAxbox3x12x4(p1,p2,p3,p4,p5,p6,za,zb,box0,box2)
      coeffs0(d3_12_4,:,:)=box0(:,:)
      coeffs2(d3_12_4,:,:)=box2(:,:)
      call qqbggAxbox3x4x12(p1,p2,p3,p4,p5,p6,za,zb,box0,box2)
      coeffs0(d3_4_12,:,:)=box0(:,:)
      coeffs2(d3_4_12,:,:)=box2(:,:)
      call qqbggAxbox3x4x12(p1,p2,p4,p3,p5,p6,za,zb,box0,box2)
      do h3=1,2
      do h4=1,2
      coeffs0(d4_3_12,h3,h4)=-box0(h4,h3)
      coeffs2(d4_3_12,h3,h4)=-box2(h4,h3)
      enddo
      enddo

      call qqbggAxtri12x3x456(p1,p2,p3,p4,p5,p6,p4,p3,za,zb,mtsq,tri0)
      coeffs0(c12_3,:,:)=tri0(:,:)
      call qqbggAxtri12x3x456(p1,p2,p4,p3,p5,p6,p3,p4,za,zb,mtsq,tri0)
      do h3=1,2
      do h4=1,2
      coeffs0(c12_4,h3,h4)=-tri0(h4,h3)
      enddo
      enddo

      call qqbggAxtri123x4x56(p1,p2,p3,p4,p5,p6,za,zb,tri0,tri2)
      coeffs0(c4_123,:,:)=tri0(:,:)
      coeffs2(c4_123,:,:)=tri2(:,:)
! permutations related to color ordering being wrong under naive 3<->4
      call qqbggAxtri123x4x56(p2,p1,p4,p3,p6,p5,zb,za,tri0,tri2)
      do h3=1,2
      do h4=1,2
        coeffs0(c3_124,h3,h4)=tri0(3-h4,3-h3)
        coeffs2(c3_124,h3,h4)=tri2(3-h4,3-h3)
      enddo
      enddo

c--- need rational parts separately to determine c2_12|34
      coeffs0(rat,2,2)=  + za(p2,p3)*za(p2,p5)*zb(p3,p1)*zb(p3,p6)*iza(
     &    p2,p4)*iza(p3,p4)*is(p5,p6)*izab2(p3,p5,p6,p3) + za(p2,p3)*
     &    za(p2,p5)*zb(p3,p6)*iza(p1,p2)*iza(p4,p3)**2*s(p3,p4)*is(p5,
     &    p6)*izab2(p3,p5,p6,p3) - za(p2,p4)*za(p2,p5)*zb(p4,p6)*iza(p1
     &    ,p2)*iza(p3,p4)**2*s(p3,p4)*is(p5,p6)*izab2(p4,p5,p6,p4) + 
     &    za(p2,p5)*zb(p4,p6)*iza(p1,p3)*iza(p3,p4)*s(p1,p4)*is(p5,p6)*
     &    izab2(p4,p5,p6,p4) + za(p2,p5)*zb(p4,p6)*iza(p1,p3)*iza(p3,p4
     &    )*s(p3,p4)*is(p5,p6)*izab2(p4,p5,p6,p4)

      coeffs0(rat,2,1)=  + za(p1,p4)*za(p2,p4)*zb(p4,p6)*iza(p1,p2)*
     &    iza(p1,p3)*izb(p5,p6)*zab2(p2,p1,p3,p6)*izab2(p3,p1,p2,p4)*
     &    izab2(p4,p5,p6,p4) - za(p2,p4)**2*za(p3,p5)*zb(p3,p6)*iza(p1,
     &    p2)*iza(p3,p4)*is(p5,p6)*izab2(p3,p1,p2,p4) + za(p2,p4)**2*
     &    za(p3,p5)*zb(p3,p6)*iza(p1,p2)*iza(p3,p4)*izab2(p3,p1,p2,p4)*
     &    idel3*d56 - za(p2,p4)**2*zb(p3,p6)*zb(p4,p6)*iza(p1,p2)*izb(
     &    p5,p6)*izab2(p3,p1,p2,p4)*idel3*d34 + za(p2,p4)*za(p3,p5)*za(
     &    p4,p5)*zb(p1,p3)*iza(p4,p3)*iza(p6,p5)*izab2(p3,p2,p1,p4)*
     &    idel3*d12 - za(p2,p4)*za(p3,p5)*iza(p1,p3)*iza(p3,p4)*zab2(p4
     &    ,p1,p3,p6)*is(p5,p6)*izab2(p3,p1,p2,p4) + 2.D0*za(p2,p4)*za(
     &    p5,p3)*zb(p3,p1)*zb(p3,p6)*izab2(p3,p1,p2,p4)*idel3 + za(p2,
     &    p4)*zb(p1,p3)*zb(p3,p6)*zb(p4,p6)*izb(p3,p4)*izb(p5,p6)*
     &    izab2(p3,p1,p2,p4)*idel3*d12 - za(p3,p5)*za(p4,p5)*zb(p1,p3)
     &    **2*iza(p6,p5)*izb(p2,p1)*izab2(p3,p2,p1,p4)*idel3*d34 + za(
     &    p3,p5)*zb(p1,p3)*zb(p2,p3)*iza(p6,p5)*izb(p2,p1)*izb(p2,p4)*
     &    zab2(p5,p2,p4,p1)*izab2(p3,p2,p1,p4)*izab2(p3,p5,p6,p3)
      coeffs0(rat,2,1) = coeffs0(rat,2,1) + 2.D0*za(p4,p2)*za(p4,p5)*
     & zb(p1,p3)*zb(p6,p4)*izab2(p3,p2,p1,p4)*idel3 - za(p4,p5)*zb(p1,
     &    p3)**2*zb(p4,p6)*izb(p2,p1)*izb(p4,p3)*is(p5,p6)*izab2(p3,p2,
     &    p1,p4) + za(p4,p5)*zb(p1,p3)**2*zb(p4,p6)*izb(p2,p1)*izb(p4,
     &    p3)*izab2(p3,p2,p1,p4)*idel3*d56 - zb(p1,p3)*zb(p4,p6)*izb(p2
     &    ,p4)*izb(p4,p3)*zab2(p5,p2,p4,p3)*is(p5,p6)*izab2(p3,p2,p1,p4
     &    )

      coeffs0(rat,1,2)=  - za(p2,p3)**2*za(p4,p5)*zb(p1,p6)*iza(p2,p4)*
     &    iza(p4,p3)*is(p5,p6)*izab2(p4,p2,p1,p3) + za(p2,p3)**2*za(p4,
     &    p5)*zb(p4,p6)*iza(p1,p2)*iza(p4,p3)*is(p5,p6)*izab2(p4,p1,p2,
     &    p3) - za(p2,p3)**2*za(p4,p5)*zb(p4,p6)*iza(p1,p2)*iza(p4,p3)*
     &    izab2(p4,p1,p2,p3)*idel3*d56 + za(p2,p3)**2*zb(p3,p6)*zb(p4,
     &    p6)*iza(p1,p2)*izb(p5,p6)*izab2(p4,p1,p2,p3)*idel3*d34 - za(
     &    p2,p3)**2*zb(p3,p6)*iza(p2,p1)*iza(p2,p4)*izb(p6,p5)*zab2(p2,
     &    p1,p4,p6)*izab2(p3,p5,p6,p3)*izab2(p4,p2,p1,p3) - za(p2,p3)*
     &    za(p3,p5)*za(p4,p5)*zb(p1,p4)*iza(p3,p4)*iza(p6,p5)*izab2(p4,
     &    p2,p1,p3)*idel3*d12 - 2.D0*za(p2,p3)*za(p5,p4)*zb(p4,p1)*zb(
     &    p4,p6)*izab2(p4,p1,p2,p3)*idel3 - za(p2,p3)*zb(p1,p4)*zb(p3,
     &    p6)*zb(p4,p6)*izb(p4,p3)*izb(p5,p6)*izab2(p4,p1,p2,p3)*idel3*
     &    d12 - za(p2,p5)*zb(p1,p4)**2*zb(p3,p6)*izb(p1,p3)*izb(p3,p4)*
     &    is(p5,p6)*izab2(p4,p1,p2,p3) - 2.D0*za(p3,p2)*za(p3,p5)*zb(p1
     &    ,p4)*zb(p6,p3)*izab2(p4,p2,p1,p3)*idel3 + za(p3,p5)*za(p4,p5)
     &    *zb(p1,p4)**2*iza(p6,p5)*izb(p2,p1)*izab2(p4,p2,p1,p3)*idel3*
     &    d34
      coeffs0(rat,1,2) = coeffs0(rat,1,2) + za(p3,p5)*zb(p1,p4)**2*zb(
     & p3,p6)*izb(p2,p1)*izb(p3,p4)*is(p5,p6)*izab2(p4,p2,p1,p3) - za(
     &    p3,p5)*zb(p1,p4)**2*zb(p3,p6)*izb(p2,p1)*izb(p3,p4)*izab2(p4,
     &    p2,p1,p3)*idel3*d56 - za(p4,p5)*zb(p1,p4)**2*iza(p5,p6)*izb(
     &    p1,p2)*izb(p1,p3)*zab2(p5,p2,p3,p1)*izab2(p4,p1,p2,p3)*izab2(
     &    p4,p5,p6,p4)

      coeffs0(rat,1,1)=  - za(p3,p5)*zb(p1,p3)*zb(p1,p6)*izb(p2,p1)*
     &    izb(p4,p3)**2*s(p3,p4)*is(p5,p6)*izab2(p3,p5,p6,p3) + za(p3,
     &    p5)*zb(p1,p6)*izb(p2,p4)*izb(p4,p3)*s(p2,p3)*is(p5,p6)*izab2(
     &    p3,p5,p6,p3) + za(p3,p5)*zb(p1,p6)*izb(p2,p4)*izb(p4,p3)*s(p3
     &    ,p4)*is(p5,p6)*izab2(p3,p5,p6,p3) + za(p4,p2)*za(p4,p5)*zb(p1
     &    ,p4)*zb(p1,p6)*izb(p1,p3)*izb(p4,p3)*is(p5,p6)*izab2(p4,p5,p6
     &    ,p4) + za(p4,p5)*zb(p1,p4)*zb(p1,p6)*izb(p2,p1)*izb(p3,p4)**2
     &    *s(p3,p4)*is(p5,p6)*izab2(p4,p5,p6,p4)

c--- integrals that use call with p3,p4 flipped
c--- must fix helicities and insert minus sign for color

c--- triangle integral from infrared relation
      coeffs0(c4_123,:,:)=(s56-s123)*(
     &  coeffs0(d4_3_12,:,:)/s34/s123-coeffs0(d3_4_12,:,:)/s34/s124
     & -coeffs0(c12_4,:,:)/(s124-s12))

      coeffs0(c3_124,:,:)=(s56-s124)*(
     &  coeffs0(d3_4_12,:,:)/s34/s124-coeffs0(d4_3_12,:,:)/s34/s123
     & -coeffs0(c12_3,:,:)/(s123-s12))

c--- triangle integral from infrared relation
      coeffs0(c3_4,:,:)=
     & -coeffs0(d4_3_12,:,:)/s123-coeffs0(d3_4_12,:,:)/s124

c--- 3-mass triangle massive coefficient from rational relation
      coeffs2(c12_34,:,:)=two*coeffs0(rat,:,:)
     &  -coeffs2(c4_123,:,:)-coeffs2(c3_124,:,:)

c--- BDK pieces
      ABDK(:,:)=czip
      
c--- fill bubble integral arrays for flip2 operation (effectively just 3 and 4 exchanged)
      xIntflip2(:)=xInt(:)
      xIntflip2(b123)=xInt(b124)
      xIntflip2(b124)=xInt(b123)
      
      ABDK(2,2)=A6axBDKpp(p1,p2,p3,p4,p5,p6,za,zb,xInt)
      ABDK(1,1)=A6axBDKpp(p2,p1,p4,p3,p6,p5,zb,za,xIntflip2) ! this is flip2
      ABDK(2,1)=A6axBDKpm(p1,p2,p3,p4,p5,p6,za,zb,xInt)
     &         +A6axBDKpm(p2,p1,p4,p3,p6,p5,zb,za,xIntflip2) ! this is flip2
      ABDK(1,2)=A6axBDKmp(p1,p2,p3,p4,p5,p6,za,zb,xInt)
     &         +A6axBDKmp(p2,p1,p4,p3,p6,p5,zb,za,xIntflip2) ! this is flip2

c--- rational terms may also be obtained this way
!      write(6,*) 'old 2,2',coeffs0(rat,2,2)
!      write(6,*) 'old 1,1',coeffs0(rat,1,1)
!      write(6,*) 'old 2,1',coeffs0(rat,2,1)
!      write(6,*) 'old 1,2',coeffs0(rat,1,2)
!      xIntzip(:)=czip            
!      coeffs0(rat,2,2)=A6axBDKpp(p1,p2,p3,p4,p5,p6,za,zb,xIntzip)
!      coeffs0(rat,1,1)=A6axBDKpp(p2,p1,p4,p3,p6,p5,zb,za,xIntzip) ! this is flip2
!      coeffs0(rat,2,1)=A6axBDKpm(p1,p2,p3,p4,p5,p6,za,zb,xIntzip)
!     &                +A6axBDKpm(p2,p1,p4,p3,p6,p5,zb,za,xIntzip) ! this is flip2
!      coeffs0(rat,1,2)=A6axBDKmp(p1,p2,p3,p4,p5,p6,za,zb,xIntzip)
!     &                +A6axBDKmp(p2,p1,p4,p3,p6,p5,zb,za,xIntzip) ! this is flip2
!      write(6,*) 'new 2,2',coeffs0(rat,2,2)
!      write(6,*) 'new 1,1',coeffs0(rat,1,1)
!      write(6,*) 'new 2,1',coeffs0(rat,2,1)
!      write(6,*) 'new 1,2',coeffs0(rat,1,2)
!      pause

      if (writecoeffs) then
        if ((p1 == 1) .and. (p2 == 2) .and. (p3 == 3) .and. (p4 == 4) .and. (p5 == 5) .and. (p6 == 6)) then
        write(6,*)
        write(6,*) 'BDK coefficients'
        write(6,*)
        do h4=2,1,-1
        do h3=2,1,-1
        write(6,51) 'easy ',h3,h4,Coeffs0(d3_12_4,h3,h4),Coeffs2(d3_12_4,h3,h4)
        write(6,51) 'h123 ',h3,h4,Coeffs0(d4_3_12,h3,h4),Coeffs2(d4_3_12,h3,h4)
        write(6,51) 'h124 ',h3,h4,Coeffs0(d3_4_12,h3,h4),Coeffs2(d3_4_12,h3,h4)
        write(6,*)
        write(6,51) '12_34',h3,h4,1.e30,Coeffs2(c12_34,h3,h4)
        write(6,51) '12_3 ',h3,h4,Coeffs0(c12_3,h3,h4),Coeffs2(c12_3,h3,h4)
        write(6,51) '12_4 ',h3,h4,Coeffs0(c12_4,h3,h4),Coeffs2(c12_4,h3,h4)
        write(6,51) '124_3',h3,h4,Coeffs0(c3_124,h3,h4),Coeffs2(c3_124,h3,h4)
        write(6,51) '123_4',h3,h4,Coeffs0(c4_123,h3,h4),Coeffs2(c4_123,h3,h4)
        write(6,51) '3_4  ',h3,h4,Coeffs0(c3_4,h3,h4),Coeffs2(c3_4,h3,h4)
        write(6,*)
        write(6,51) 'ABDK ',h3,h4,ABDK(h3,h4)
        write(6,51)
        enddo
        enddo
        stop
        endif
   51   format(a8,2i5,4(f16.10,'   & '))
      endif

c--- Now add contributions of boxes and triangles
      do k=1,c3_4
        ABDK(:,:)=ABDK(:,:)+xInt(k)*(Coeffs0(k,:,:)+mtsq*Coeffs2(k,:,:))
      enddo

      return
      end
