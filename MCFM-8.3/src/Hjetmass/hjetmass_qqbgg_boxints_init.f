      subroutine hjetmass_qqbgg_boxints_init(i1,i2,i3,i4,boxints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'

          include 'constants.f'
          real*16 s(mxpart,mxpart)
          include 'scale.f'
          include 'masses.f'

          integer i1,i2,i3,i4
          integer ord
          complex*32 boxints(3)
          real*16 t
          real*16 mt2, mhsq, qmusq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)
        qmusq = musq
      

        boxints(:) = (0q0,0q0)

c    "(   0,   0, s12,mhsq, s34,s123;mtsq,mtsq,mtsq,mtsq)",
        boxints(1) = qlI4q(0q0, 0q0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i3),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

c    "(   0,   0, s12,mhsq, s34,s124;mtsq,mtsq,mtsq,mtsq)",
        boxints(2) = qlI4q(0q0, 0q0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

c    "(   0, s12,   0,mhsq,s123,s124;mtsq,mtsq,mtsq,mtsq)"
        boxints(3) = qlI4q(0q0, s(i1,i2), 0q0, mhsq, t(i1,i2,i3), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, qmusq, ord)

      end subroutine
