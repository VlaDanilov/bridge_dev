      subroutine hjetmass_qqbgg_bubints_init_dp
     &   (i1,i2,i3,i4,bubints,ord,s)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_decl.f'
      include 'scale.f'
      include 'masses.f'
      integer i1,i2,i3,i4
      include 'qqbggnames.f'
      integer ord
      complex(dp)::bubints(Nbubints)
      real*8 t
      real*8 mt2, mhsq,s1234

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      mt2 = mt**2
      s1234=s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
      bubints(:) = (0d0,0d0)

c   "(mhsq;mtsq,mtsq)",
        bubints(b1234) = qlI2(s1234, mt2, mt2, musq, ord)

c   "(s123;mtsq,mtsq)",
        bubints(b123) = qlI2(t(i1,i2,i3), mt2, mt2, musq, ord)

c   "(s124;mtsq,mtsq)",
        bubints(b124) = qlI2(t(i1,i2,i4), mt2, mt2, musq, ord)

c   "( s12;mtsq,mtsq)",
        bubints(b12) = qlI2(s(i1,i2), mt2, mt2, musq, ord)

c   "( s34;mtsq,mtsq)"]
        bubints(b34) = qlI2(s(i3,i4), mt2, mt2, musq, ord)

      end subroutine
