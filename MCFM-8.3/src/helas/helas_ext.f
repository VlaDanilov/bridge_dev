!     ------------------------------
!     Routines to compute amplitudes
!     ------------------------------


      SUBROUTINE IOVCXX(fi,fo,vc,g,vertex)
      IMPLICIT NONE
      COMPLEX*16 fi(6),fo(6),vc(6),vertex,g(2)
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))

      vertex=(0d0,0d0)
      IF ( g(1) .NE. 0d0 ) THEN
      vertex =  g(1)*( (fo(3)*fi(1)+fo(4)*fi(2))*vc(1)
     &                +(fo(3)*fi(2)+fo(4)*fi(1))*vc(2)
     &                -(fo(3)*fi(2)-fo(4)*fi(1))*vc(3)*c_imag
     &                +(fo(3)*fi(1)-fo(4)*fi(2))*vc(4))
      ENDIF

      IF ( g(2) .NE. 0d0 ) THEN
      vertex = vertex+ g(2)*(
     &                   (fo(1)*fi(3)+fo(2)*fi(4))*vc(1)
     &                  -(fo(1)*fi(4)+fo(2)*fi(3))*vc(2)
     &                  +(fo(1)*fi(4)-fo(2)*fi(3))*vc(3)*c_imag
     &                  -(fo(1)*fi(3)-fo(2)*fi(4))*vc(4))
      ENDIF

      RETURN
      END

!     =====================================================

      SUBROUTINE IOVSMX(fi,fo,vc,g,ISIGN,vertex)
      IMPLICIT NONE
      INTEGER ISIGN
      COMPLEX*16 fi(6),fo(6),vc(6),vertex,g(2)
      REAL*8 Q(4)
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))

!      Q(1)=DREAL(FI(5)-FO(5))*FLOAT(ISIGN) ! Instable numerically!!!
!      Q(2)=DREAL(FI(6)-FO(6))*FLOAT(ISIGN) !
!      Q(3)=DIMAG(FI(6)-FO(6))*FLOAT(ISIGN) !
!      Q(4)=DIMAG(FI(5)-FO(5))*FLOAT(ISIGN) !
      Q(1)=DREAL(VC(5))*FLOAT(ISIGN)
      Q(2)=DREAL(VC(6))*FLOAT(ISIGN)
      Q(3)=DIMAG(VC(6))*FLOAT(ISIGN)
      Q(4)=DIMAG(VC(5))*FLOAT(ISIGN)

      vertex=(0d0,0d0)
      IF ( g(1) .NE. 0d0 ) THEN
      vertex = g(1)*(
     .         -(fo(1)*fi(2)+fo(2)*fi(1))*(vc(1)*q(2)-vc(2)*q(1))
     .  +c_imag*(fo(1)*fi(2)-fo(2)*fi(1))*(vc(1)*q(3)-vc(3)*q(1))
     .         -(fo(1)*fi(1)-fo(2)*fi(2))*(vc(1)*q(4)-vc(4)*q(1))
     .  +c_imag*(fo(1)*fi(1)-fo(2)*fi(2))*(vc(2)*q(3)-vc(3)*q(2))
     .         -(fo(1)*fi(2)-fo(2)*fi(1))*(vc(2)*q(4)-vc(4)*q(2))
     .  +c_imag*(fo(1)*fi(2)+fo(2)*fi(1))*(vc(3)*q(4)-vc(4)*q(3)))
      ENDIF

      IF ( g(2) .NE. 0d0 ) THEN
      vertex = vertex + g(2)*(
     .         +(fo(3)*fi(4)+fo(4)*fi(3))*(vc(1)*q(2)-vc(2)*q(1))
     .  -c_imag*(fo(3)*fi(4)-fo(4)*fi(3))*(vc(1)*q(3)-vc(3)*q(1))
     .         +(fo(3)*fi(3)-fo(4)*fi(4))*(vc(1)*q(4)-vc(4)*q(1))
     .  +c_imag*(fo(3)*fi(3)-fo(4)*fi(4))*(vc(2)*q(3)-vc(3)*q(2))
     .         -(fo(3)*fi(4)-fo(4)*fi(3))*(vc(2)*q(4)-vc(4)*q(2))
     .  +c_imag*(fo(3)*fi(4)+fo(4)*fi(3))*(vc(3)*q(4)-vc(4)*q(3)))
      ENDIF

      RETURN
      END

!     =====================================================

      SUBROUTINE IOIOVX(fi1,fo1,fi2,fo2,g1,g2,vertex)
      IMPLICIT NONE

      COMPLEX*16 fi1(6),fo1(6),fi2(6),fo2(6),vertex
      REAL*8 g1(2),g2(2)
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))
      COMPLEX*16 J1(4),J2(4)

      J1(1) = g1(1) * (fo1(3)*fi1(1)+fo1(4)*fi1(2))
     &         + g1(2) * (fo1(1)*fi1(3)+fo1(2)*fi1(4))
      J1(2) = g1(1) * (-fo1(3)*fi1(2)-fo1(4)*fi1(1))
     &         + g1(2) * (fo1(1)*fi1(4)+fo1(2)*fi1(3))
      J1(3) = g1(1) * (fo1(3)*fi1(2)-fo1(4)*fi1(1))*c_imag
     &         + g1(2) * (-fo1(1)*fi1(4)+fo1(2)*fi1(3))*c_imag
      J1(4) = g1(1) * (-fo1(3)*fi1(1)+fo1(4)*fi1(2))
     &         + g1(2) * (fo1(1)*fi1(3)-fo1(2)*fi1(4))

      J2(1) = g2(1) * (fo2(3)*fi2(1)+fo2(4)*fi2(2))
     &         + g2(2) * (fo2(1)*fi2(3)+fo2(2)*fi2(4))
      J2(2) = g2(1) * (-fo2(3)*fi2(2)-fo2(4)*fi2(1))
     &         + g2(2) * (fo2(1)*fi2(4)+fo2(2)*fi2(3))
      J2(3) = g2(1) * (fo2(3)*fi2(2)-fo2(4)*fi2(1))*c_imag
     &         + g2(2) * (-fo2(1)*fi2(4)+fo2(2)*fi2(3))*c_imag
      J2(4) = g2(1) * (-fo2(3)*fi2(1)+fo2(4)*fi2(2))
     &         + g2(2) * (fo2(1)*fi2(3)-fo2(2)*fi2(4))

      vertex = J1(1)*J2(1)-J1(2)*J2(2)-J1(3)*J2(3)-J1(4)*J2(4)
      RETURN
      END

!     =====================================================

      SUBROUTINE IOIOSX(fi1,fo1,fi2,fo2,g1,g2,vertex)
      IMPLICIT NONE

      COMPLEX*16 fi1(6),fo1(6),fi2(6),fo2(6),vertex
      REAL*8 g1(2),g2(2)
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))
      COMPLEX*16 S1,S2

      S1 = g1(1)*(fi1(1)*fo1(1)+fi1(2)*fo1(2))
     &   + g1(2)*(fi1(3)*fo1(3)+fi1(4)*fo1(4))
      S2 = g2(1)*(fi2(1)*fo2(1)+fi2(2)*fo2(2))
     &   + g2(2)*(fi2(3)*fo2(3)+fi2(4)*fo2(4))
      
      vertex = S1*S2
      RETURN
      END


!     -------------------------------------------
!     Routines to compute in lines from Out lines
!     -------------------------------------------


      SUBROUTINE FVOCXX(fo,vc,g,fmass,fwidth,fvo)
      IMPLICIT NONE
      COMPLEX*16 fo(6),vc(6),fvo(6),g(2)

      REAL*8 pf(0:3),fmass,fwidth,pf2
      COMPLEX*16 FVO1(4),FVO2(4),sl1,sl2,sr1,sr2,d
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))

      fvo(5) = fo(5)+vc(5)
      fvo(6) = fo(6)+vc(6)

      pf(0)=DREAL(fvo(5))
      pf(1)=DREAL(fvo(6))
      pf(2)=DIMAG(fvo(6))
      pf(3)=DIMAG(fvo(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d=-1d0/dcmplx( pf2-fmass**2,max(sign(fmass*fwidth,pf2),0d0))

      IF (g(1) .NE. 0d0) THEN
        sl1= (vc(1)+vc(4))*fo(3)
     &    +(vc(2)+c_imag*vc(3))*fo(4)
        sl2= (vc(2)-c_imag*vc(3))*fo(3)
     &    +(vc(1)-vc(4))*fo(4)
        fvo1(1) = g(1)*fmass*sl1*d
        fvo1(2) = g(1)*fmass*sl2*d
        fvo1(3) = g(1)*( (pf(0)-pf(3))*sl1-fvo(6)*sl2)*d
        fvo1(4) = g(1)*(-DCONJG(fvo(6))*sl1+(pf(0)+pf(3))*sl2)*d
      ELSE
        fvo1(1)=(0d0,0d0)
        fvo1(2)=(0d0,0d0)
        fvo1(3)=(0d0,0d0)
        fvo1(4)=(0d0,0d0)
      ENDIF

      IF (g(2) .NE. 0d0) THEN
        sr1= (vc(1)-vc(4))*fo(1)
     &    -(vc(2)+c_imag*vc(3))*fo(2)
        sr2=-(vc(2)-c_imag*vc(3))*fo(1)
     &    +(vc(1)+vc(4))*fo(2)
        fvo2(1) = g(2)*( (pf(0)+pf(3))*sr1+fvo(6)*sr2)*d
        fvo2(2) = g(2)*( dconjg(fvo(6))*sr1 +(pf(0)-pf(3))*sr2)*d
        fvo2(3) = g(2)*fmass*sr1*d
        fvo2(4) = g(2)*fmass*sr2*d
      ELSE
        fvo2(1)=(0d0,0d0)
        fvo2(2)=(0d0,0d0)
        fvo2(3)=(0d0,0d0)
        fvo2(4)=(0d0,0d0)
      ENDIF
      FVO(1)=FVO1(1)+FVO2(1)
      FVO(2)=FVO1(2)+FVO2(2)
      FVO(3)=FVO1(3)+FVO2(3)
      FVO(4)=FVO1(4)+FVO2(4)
      RETURN
      END

!     =====================================================

      SUBROUTINE FVOSmX(fo,vc,g,fmass,fwidth,ISIGN,fvo)
      IMPLICIT NONE
      COMPLEX*16 fo(6),vc(6),fvo(6),g(2)
      REAL*8 pf(0:3),fmass,fwidth,pf2
      COMPLEX*16 sl1,sl2,sr1,sr2,d
      INTEGER ISIGN
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))
      REAL*8 Q(4)
      COMPLEX*16 VQ01,VQ02,VQ03,VQ12,VQ13,VQ23,FVO1(4),FVO2(4)

      Q(1)=DREAL(VC(5))
      Q(2)=DREAL(VC(6))
      Q(3)=DIMAG(VC(6))
      Q(4)=DIMAG(VC(5))
      VQ01=(VC(1)*Q(2)-VC(2)*Q(1))*DFLOAT(ISIGN)
      VQ02=(VC(1)*Q(3)-VC(3)*Q(1))*DFLOAT(ISIGN)
      VQ03=(VC(1)*Q(4)-VC(4)*Q(1))*DFLOAT(ISIGN)
      VQ12=(VC(2)*Q(3)-VC(3)*Q(2))*DFLOAT(ISIGN)
      VQ13=(VC(2)*Q(4)-VC(4)*Q(2))*DFLOAT(ISIGN)
      VQ23=(VC(3)*Q(4)-VC(4)*Q(3))*DFLOAT(ISIGN)

      fvo(5) = fo(5)+vc(5)
      fvo(6) = fo(6)+vc(6)

      pf(0)=DBLE(fvo(5))
      pf(1)=DBLE(fvo(6))
      pf(2)=DIMAG(fvo(6))
      pf(3)=DIMAG(fvo(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d=-1d0/dcmplx( pf2-fmass**2,max(sign(fmass*fwidth,pf2),0d0))

      IF (g(1) .NE. 0d0) THEN
        sl1 = (-VQ03+c_imag*VQ12)*FO(1) +
     .        (-VQ01-c_imag*VQ02+VQ13+c_imag*VQ23)*FO(2)
        sl2 = (-VQ01+c_imag*VQ02-VQ13+c_imag*VQ23)*FO(1) +
     .        (VQ03-c_imag*VQ12)*FO(2)
        fvo1(1) = g(1)*fmass*sl1*d
	fvo1(2) = g(1)*fmass*sl2*d
	fvo1(3) = g(1)*d*
     .          ( (pf(0)-pf(3)) * sl1 - fvo(6) * sl2 )
        fvo1(4) = g(1)*d*
     .          ( -dconjg(fvo(6)) * sl1 + (pf(0)+pf(3)) * sl2 )
      ELSE
        fvo1(1)=(0d0,0d0)
        fvo1(2)=(0d0,0d0)
        fvo1(3)=(0d0,0d0)
        fvo1(4)=(0d0,0d0)
      ENDIF
      IF (g(2) .NE. 0d0) THEN
        sr1 = (VQ03+c_imag*VQ12)*FO(3) +
     .        (VQ01+c_imag*VQ02+VQ13+c_imag*VQ23)*FO(4)
        sr2 = (VQ01-c_imag*VQ02-VQ13+c_imag*VQ23)*FO(3) +
     .        (-VQ03-c_imag*VQ12)*FO(4)
        fvo2(1) = g(2)*d*
     .          ( (pf(0)+pf(3)) * sr1 + fvo(6) * sr2 )
        fvo2(2) = g(2)*d*
     .          ( dconjg(fvo(6)) * sr1 + (pf(0)-pf(3)) * sr2 )
        fvo2(3) = g(2)*fmass*sr1*d
        fvo2(4) = g(2)*fmass*sr2*d
      ELSE
        fvo2(1)=(0d0,0d0)
        fvo2(2)=(0d0,0d0)
        fvo2(3)=(0d0,0d0)
        fvo2(4)=(0d0,0d0)
      ENDIF
      FVO(1)=FVO1(1)+FVO2(1)
      FVO(2)=FVO1(2)+FVO2(2)
      FVO(3)=FVO1(3)+FVO2(3)
      FVO(4)=FVO1(4)+FVO2(4)
      RETURN
      END

!     =====================================================

      SUBROUTINE FVOIOX(fo1,fi2,fo2,g1,g2,fmass,fwidth,fi1)
      IMPLICIT NONE
      COMPLEX*16 fo1(6),fi2(6),fo2(6),fi1(6),sl1,sl2,sr1,sr2,d
      REAL*8 g1(2),g2(2),pf(0:3),fmass,fwidth,pf2
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))
      COMPLEX*16 VC(4),FL(4),FR(4)

      VC(1) = g2(1) * (fo2(3)*fi2(1)+fo2(4)*fi2(2))
     &         + g2(2) * (fo2(1)*fi2(3)+fo2(2)*fi2(4))
      VC(2) = g2(1) * (-fo2(3)*fi2(2)-fo2(4)*fi2(1))
     &         + g2(2) * (fo2(1)*fi2(4)+fo2(2)*fi2(3))
      VC(3) = g2(1) * (fo2(3)*fi2(2)-fo2(4)*fi2(1))*c_imag
     &         + g2(2) * (-fo2(1)*fi2(4)+fo2(2)*fi2(3))*c_imag
      VC(4) = g2(1) * (-fo2(3)*fi2(1)+fo2(4)*fi2(2))
     &         + g2(2) * (fo2(1)*fi2(3)-fo2(2)*fi2(4))

      fi1(5) = fo1(5)+fo2(5)-fi2(5)
      fi1(6) = fo1(6)+fo2(6)-fi2(6)

      pf(0)=DREAL(fi1(5))
      pf(1)=DREAL(fi1(6))
      pf(2)=DIMAG(fi1(6))
      pf(3)=DIMAG(fi1(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d=-1d0/DCMPLX( pf2-fmass**2,MAX(SIGN(fmass*fwidth,pf2),0d0))

      IF (g1(1) .NE. 0d0) THEN
        sl1 = (VC(1)+VC(4))*fo1(3)
     &      + (VC(2)+c_imag*VC(3))*fo1(4)
        sl2 = (VC(2)-c_imag*VC(3))*fo1(3)
     &      + (VC(1)-VC(4))*fo1(4)
        FL(1) = g1(1)*fmass*sl1*d
        FL(2) = g1(1)*fmass*sl2*d
        FL(3) = g1(1)*( (pf(0)-pf(3))*sl1-fi1(6)*sl2)*d
        FL(4) = g1(1)*(-dconjg(fi1(6))*sl1 +(pf(0)+pf(3))*sl2)*d
      ELSE
        FL(1) = (0d0,0d0)
        FL(2) = (0d0,0d0)
        FL(3) = (0d0,0d0)
        FL(4) = (0d0,0d0)
      ENDIF

      IF (g1(2) .NE. 0d0) THEN
        sr1 = (VC(1)-VC(4))*fo1(1)
     &      - (VC(2)+c_imag*VC(3))*fo1(2)
        sr2 = -(VC(2)-c_imag*VC(3))*fo1(1)
     &      + (VC(1)+       VC(4))*fo1(2)
        FR(1) = g1(2)*( (pf(0)+pf(3))*sr1 + fi1(6)*sr2)*d
        FR(2) = g1(2)*( dconjg(fi1(6))*sr1 + (pf(0)-pf(3))*sr2)*d
        FR(3) = g1(2)*fmass*sr1*d
        FR(4) = g1(2)*fmass*sr2*d
      ELSE
        FR(1) = (0d0,0d0)
        FR(2) = (0d0,0d0)
        FR(3) = (0d0,0d0)
        FR(4) = (0d0,0d0)
      ENDIF
      fi1(1)=FL(1)+FR(1)
      fi1(2)=FL(2)+FR(2)
      fi1(3)=FL(3)+FR(3)
      fi1(4)=FL(4)+FR(4)
      RETURN
      END

!     =====================================================

      SUBROUTINE FSOIOX(fo1,fi2,fo2,g1,g2,fmass,fwidth,fi1)
      IMPLICIT NONE
      COMPLEX*16 fo1(6),fi2(6),fo2(6),fi1(6),sl1,sl2,sr1,sr2,ds
      REAL*8 g1(2),g2(2),pf(0:3),fmass,fwidth,pf2,p0p3,p0m3
      COMPLEX*16 S2

      S2 = g2(1)*(fi2(1)*fo2(1)+fi2(2)*fo2(2))
     &   + g2(2)*(fi2(3)*fo2(3)+fi2(4)*fo2(4))

      fi1(5) = fo1(5)+fo2(5)-fi2(5)
      fi1(6) = fo1(6)+fo2(6)-fi2(6)

      pf(0)=DREAL(fi1(5))
      pf(1)=DREAL(fi1(6))
      pf(2)=DIMAG(fi1(6))
      pf(3)=DIMAG(fi1(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      ds=-S2/DCMPLX(pf2-fmass**2,MAX(DSIGN(fmass*fwidth,pf2),0d0))
      p0p3=pf(0)+pf(3)
      p0m3=pf(0)-pf(3)
      sl1=g1(2)*(p0p3*fo1(3)+fi1(6) *fo1(4))
      sl2=g1(2)*(p0m3*fo1(4)+DCONJG(fi1(6))*fo1(3))
      sr1=g1(1)*(p0m3*fo1(1)-fi1(6) *fo1(2))
      sr2=g1(1)*(p0p3*fo1(2)-DCONJG(fi1(6))*fo1(1))

      fi1(1) = ( g1(1)*fmass*fo1(1) + sl1 )*ds
      fi1(2) = ( g1(1)*fmass*fo1(2) + sl2 )*ds
      fi1(3) = ( g1(2)*fmass*fo1(3) + sr1 )*ds
      fi1(4) = ( g1(2)*fmass*fo1(4) + sr2 )*ds
      RETURN
      END


!     -------------------------------------------
!     Routines to compute out lines from In lines
!     -------------------------------------------


      SUBROUTINE FVICXX(fi,vc,g,fmass,fwidth,fvi)
      IMPLICIT NONE
      COMPLEX*16 fi(6),vc(6),fvi(6),g(2)
      REAL*8 pf(0:3),fmass,fwidth,pf2

      COMPLEX*16 FVI1(4),FVI2(4),sl1,sl2,sr1,sr2,d
      COMPLEX*16 c_imag
      parameter (c_imag=(0d0,1d0))

      fvi(5) = fi(5)-vc(5)
      fvi(6) = fi(6)-vc(6)

      pf(0)=DREAL(fvi(5))
      pf(1)=DREAL(fvi(6))
      pf(2)=DIMAG(fvi(6))
      pf(3)=DIMAG(fvi(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d=-1d0/dcmplx( pf2-fmass**2,max(sign(fmass*fwidth,pf2),0d0))

      IF (g(1) .NE. 0d0) THEN
        sl1= (vc(1)+vc(4))*fi(1)
     &    +(vc(2)-c_imag*vc(3))*fi(2)
        sl2= (vc(2)+c_imag*vc(3))*fi(1)
     &    +(vc(1)-vc(4))*fi(2)
        fvi1(1) = g(1)*((pf(0)-pf(3))*sl1 -dconjg(fvi(6))*sl2)*d
        fvi1(2) = g(1)*(-fvi(6)*sl1 +(pf(0)+pf(3))*sl2)*d
        fvi1(3) = g(1)*fmass*sl1*d
        fvi1(4) = g(1)*fmass*sl2*d
      ELSE
        fvi1(1)=(0d0,0d0)
        fvi1(2)=(0d0,0d0)
        fvi1(3)=(0d0,0d0)
        fvi1(4)=(0d0,0d0)
      ENDIF

      IF ( g(2) .NE. 0d0 ) THEN
        sr1= (vc(1)-vc(4))*fi(3)
     &    -(vc(2)-c_imag*vc(3))*fi(4)
        sr2=-(vc(2)+c_imag*vc(3))*fi(3)
     &    +(vc(1)+vc(4))*fi(4)
        fvi2(1) = g(2)*fmass*sr1*d
        fvi2(2) = g(2)*fmass*sr2*d
        fvi2(3) = g(2)*((pf(0)+pf(3))*sr1 +dconjg(fvi(6))*sr2)*d
        fvi2(4) = g(2)*(fvi(6)*sr1 +(pf(0)-pf(3))*sr2)*d
      ELSE
        fvi2(1)=(0d0,0d0)
        fvi2(2)=(0d0,0d0)
        fvi2(3)=(0d0,0d0)
        fvi2(4)=(0d0,0d0)
      ENDIF
      FVI(1)=FVI1(1)+FVI2(1)
      FVI(2)=FVI1(2)+FVI2(2)
      FVI(3)=FVI1(3)+FVI2(3)
      FVI(4)=FVI1(4)+FVI2(4)
      RETURN
      END

!     =====================================================

      SUBROUTINE FVISMX(fi,vc,g,fmass,fwidth,ISIGN,fvi)
      IMPLICIT NONE
      COMPLEX*16 fi(6),vc(6),fvi(6),g(2)
      INTEGER ISIGN
      REAL*8 pf(0:3),fmass,fwidth,pf2

      COMPLEX*16 c_imag
      PARAMETER(c_imag=(0d0,1d0))
      REAL*8 Q(4)
      COMPLEX*16 VQ01,VQ02,VQ03,VQ12,VQ13,VQ23,FVI1(4),FVI2(4)
      COMPLEX*16 sl1,sl2,sr1,sr2,d

      Q(1)=DREAL(VC(5))
      Q(2)=DREAL(VC(6))
      Q(3)=DIMAG(VC(6))
      Q(4)=DIMAG(VC(5))
      VQ01=(VC(1)*Q(2)-VC(2)*Q(1))*DFLOAT(ISIGN)
      VQ02=(VC(1)*Q(3)-VC(3)*Q(1))*DFLOAT(ISIGN)
      VQ03=(VC(1)*Q(4)-VC(4)*Q(1))*DFLOAT(ISIGN)
      VQ12=(VC(2)*Q(3)-VC(3)*Q(2))*DFLOAT(ISIGN)
      VQ13=(VC(2)*Q(4)-VC(4)*Q(2))*DFLOAT(ISIGN)
      VQ23=(VC(3)*Q(4)-VC(4)*Q(3))*DFLOAT(ISIGN)

      fvi(5) = fi(5)-vc(5)
      fvi(6) = fi(6)-vc(6)

      pf(0)=DBLE(fvi(5))
      pf(1)=DBLE(fvi(6))
      pf(2)=DIMAG(fvi(6))
      pf(3)=DIMAG(fvi(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d=-1d0/dcmplx( pf2-fmass**2,max(sign(fmass*fwidth,pf2),0d0))

      IF (g(1) .NE. 0d0) THEN
        sl1 = (-VQ03+c_imag*VQ12)*FI(1) +
     .        (-VQ01+c_imag*VQ02-VQ13+c_imag*VQ23)*FI(2)
        sl2 = (-VQ01-c_imag*VQ02+VQ13+c_imag*VQ23)*FI(1) +
     .        (VQ03-c_imag*VQ12)*FI(2)
        fvi1(1) = g(1)*fmass*sl1*d
        fvi1(2) = g(1)*fmass*sl2*d
        fvi1(3) = g(1)*d*
     .          ( (pf(0)+pf(3)) * sl1 + dconjg(fvi(6)) * sl2 )
        fvi1(4) = g(1)*d*
     .          ( fvi(6) * sl1 + (pf(0)-pf(3)) * sl2 )
      ELSE
        fvi1(1)=(0d0,0d0)
        fvi1(2)=(0d0,0d0)
        fvi1(3)=(0d0,0d0)
        fvi1(4)=(0d0,0d0)
      ENDIF
      IF (g(2) .NE. 0d0) THEN
        sr1 = (VQ03+c_imag*VQ12)*FI(3) +
     .        (VQ01-c_imag*VQ02-VQ13+c_imag*VQ23)*FI(4)
        sr2 = (VQ01+c_imag*VQ02+VQ13+c_imag*VQ23)*FI(3) +
     .        (-VQ03-c_imag*VQ12)*FI(4)
        fvi2(1) = g(2)*d*
     .         ( (pf(0)-pf(3)) * sr1 - dconjg(fvi(6))*sr2 )
        fvi2(2) = g(2)*d*
     .         ( -fvi(6) * sr1 + (pf(0)+pf(3)) * sr2)
        fvi2(3) = g(2)*fmass*sr1*d
        fvi2(4) = g(2)*fmass*sr2*d
      ELSE
        fvi2(1)=(0d0,0d0)
        fvi2(2)=(0d0,0d0)
        fvi2(3)=(0d0,0d0)
        fvi2(4)=(0d0,0d0)
      ENDIF
      FVI(1)=FVI1(1)+FVI2(1)
      FVI(2)=FVI1(2)+FVI2(2)
      FVI(3)=FVI1(3)+FVI2(3)
      FVI(4)=FVI1(4)+FVI2(4)
      RETURN
      END

!     =====================================================

      SUBROUTINE FVIIOX(fi1,fi2,fo2,g1,g2,fmass,fwidth,fo1)
      IMPLICIT NONE
      COMPLEX*16 fi1(6),fi2(6),fo2(6),fo1(6),sl1,sl2,sr1,sr2,d
      REAL*8 g1(2),g2(2),pf(0:3),fmass,fwidth,pf2
      COMPLEX*16 c_imag
      PARAMETER (c_imag=(0d0,1d0))
      COMPLEX*16 VC(4),FL(4),FR(4)

      VC(1) = g2(1) * (fo2(3)*fi2(1)+fo2(4)*fi2(2))
     &         + g2(2) * (fo2(1)*fi2(3)+fo2(2)*fi2(4))
      VC(2) = g2(1) * (-fo2(3)*fi2(2)-fo2(4)*fi2(1))
     &         + g2(2) * (fo2(1)*fi2(4)+fo2(2)*fi2(3))
      VC(3) = g2(1) * (fo2(3)*fi2(2)-fo2(4)*fi2(1))*c_imag
     &         + g2(2) * (-fo2(1)*fi2(4)+fo2(2)*fi2(3))*c_imag
      VC(4) = g2(1) * (-fo2(3)*fi2(1)+fo2(4)*fi2(2))
     &         + g2(2) * (fo2(1)*fi2(3)-fo2(2)*fi2(4))

      fo1(5) = fi1(5)-fo2(5)+fi2(5)
      fo1(6) = fi1(6)-fo2(6)+fi2(6)

      pf(0)=DREAL(fo1(5))
      pf(1)=DREAL(fo1(6))
      pf(2)=DIMAG(fo1(6))
      pf(3)=DIMAG(fo1(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      d=-1d0/dcmplx( pf2-fmass**2,max(sign(fmass*fwidth,pf2),0d0))

      IF (g1(1) .NE. 0d0) THEN
        sl1 = (VC(1)+VC(4))*fi1(1)
     &      + (VC(2)-c_imag*VC(3))*fi1(2)
        sl2 = (VC(2)+c_imag*VC(3))*fi1(1)
     &      + (VC(1)-VC(4))*fi1(2)
        FL(1) = g1(1)*((pf(0)-pf(3))*sl1 - DCONJG(fo1(6))*sl2)*d
        FL(2) = g1(1)*(-fo1(6)*sl1 +(pf(0)+pf(3))*sl2)*d
        FL(3) = g1(1)*fmass*sl1*d
        FL(4) = g1(1)*fmass*sl2*d
      ELSE
        FL(1) = (0d0,0d0)
        FL(2) = (0d0,0d0)
        FL(3) = (0d0,0d0)
        FL(4) = (0d0,0d0)
      ENDIF

      IF (g1(2) .NE. 0d0) THEN
        sr1 = (VC(1)-VC(4))*fi1(3)
     &      - (VC(2)-c_imag*VC(3))*fi1(4)
        sr2 = -(VC(2)+c_imag*VC(3))*fi1(3)
     &       + (VC(1)+VC(4))*fi1(4)

        FR(1) = g1(2)*fmass*sr1*d
        FR(2) = g1(2)*fmass*sr2*d
        FR(3) = g1(2)*((pf(0)+pf(3))*sr1 + DCONJG(fo1(6))*sr2)*d
        FR(4) = g1(2)*(fo1(6)*sr1 +(pf(0)-pf(3))*sr2)*d
      ELSE
        FR(1) = (0d0,0d0)
        FR(2) = (0d0,0d0)
        FR(3) = (0d0,0d0)
        FR(4) = (0d0,0d0)
      ENDIF
      fo1(1)=FL(1)+FR(1)
      fo1(2)=FL(2)+FR(2)
      fo1(3)=FL(3)+FR(3)
      fo1(4)=FL(4)+FR(4)
      RETURN
      END

!     =====================================================

      SUBROUTINE FSIIOX(fi1,fi2,fo2,g1,g2,fmass,fwidth,fo1)
      IMPLICIT NONE
      COMPLEX*16 fi1(6),fi2(6),fo2(6),fo1(6),sl1,sl2,sr1,sr2,ds
      REAL*8 g1(2),g2(2),pf(0:3),fmass,fwidth,pf2,p0p3,p0m3
      COMPLEX*16 S2

      S2 = g2(1)*(fi2(1)*fo2(1)+fi2(2)*fo2(2))
     &   + g2(2)*(fi2(3)*fo2(3)+fi2(4)*fo2(4))

      fo1(5) = fi1(5)-fo2(5)+fi2(5)
      fo1(6) = fi1(6)-fo2(6)+fi2(6)

      pf(0)=DREAL(fo1(5))
      pf(1)=DREAL(fo1(6))
      pf(2)=DIMAG(fo1(6))
      pf(3)=DIMAG(fo1(5))
      pf2=pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

      ds=-S2/DCMPLX(pf2-fmass**2,MAX(DSIGN(fmass*fwidth,pf2),0d0))
      p0p3=pf(0)+pf(3)
      p0m3=pf(0)-pf(3)
      sl1=g1(1)*(p0p3*fi1(1)+DCONJG(fo1(6))*fi1(2))
      sl2=g1(1)*(p0m3*fi1(2)+fo1(6)*fi1(1))
      sr1=g1(2)*(p0m3*fi1(3)-DCONJG(fo1(6))*fi1(4))
      sr2=g1(2)*(p0p3*fi1(4)-fo1(6)*fi1(3))

      fo1(1) = ( g1(1)*fmass*fi1(1) + sr1 )*ds
      fo1(2) = ( g1(1)*fmass*fi1(2) + sr2 )*ds
      fo1(3) = ( g1(2)*fmass*fi1(3) + sl1 )*ds
      fo1(4) = ( g1(2)*fmass*fi1(4) + sl2 )*ds
      RETURN
      END
