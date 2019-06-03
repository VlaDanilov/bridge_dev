*********************************************************
*AUTHOR: FABIO MALTONI                                  *
*DATE  : 2/15/2003                                     *
*NOTES : PROGRAM GENERATED BY ToFortran.m               *
*        AMPLITUDES CALCULATED BY ALBERTO FRIZZO        *
*********************************************************

      function  Appppp(I1,I2,I3,I4,I5)
* ---------------------------------------------------------------------
*                            1+ 2+ 3+ 4+ 5+                            
* ---------------------------------------------------------------------


      IMPLICIT NONE
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      complex(dp) Appppp
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer::I1,I2,I3,I4,I5
      real(dp)::MSQ      

      MSQ=
     & +s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i1,i5)
     &          +s(i2,i3)+s(i2,i4)+s(i2,i5)
     &                   +s(i3,i4)+s(i3,i5)
     &                            +s(i4,i5)

      Appppp=MSQ**2/(za(i1,i2)*za(i2,i3)*za(i3,i4)*za(i4,i5)*za(i5,i1))

      END

