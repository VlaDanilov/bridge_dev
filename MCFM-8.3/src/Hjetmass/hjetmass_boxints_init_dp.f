      subroutine hjetmass_boxints_init_dp(i1,i2,i3,i4,boxints,ord,s)
        use mod_qcdloop_c
          implicit none
          include 'types.f'
          include 'mxpart.f'

          include 'constants.f'
          include 'scale.f'
          include 'masses.f'

          integer i1,i2,i3,i4
          integer ord
          complex*16 boxints(16)
          real*8 t, s(mxpart,mxpart)
          real*8 mt2, mhsq

        t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

        mt2 = mt**2
        mhsq = s(i1,i2)+s(i1,i3)+s(i1,i4)+s(i2,i3)+s(i2,i4)+s(i3,i4)

        boxints(:) = (0d0,0d0)

!      (   0,   0,   0,s123, s12, s23;mts,mts,mts,mts)
        boxints(1) =
     &  qlI4(0d0, 0d0, 0d0, t(i1,i2,i3), s(i1,i2), s(i2,i3),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0,   0,s124, s14, s12;mts,mts,mts,mts)
        boxints(2) =
     &  qlI4(0d0, 0d0, 0d0, t(i1,i2,i4), s(i1,i4), s(i1,i2),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0,   0,s134, s34, s14;mts,mts,mts,mts)
        boxints(3) =
     &  qlI4(0d0, 0d0, 0d0, t(i1,i3,i4), s(i3,i4), s(i1,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0,   0,s234, s23, s34;mts,mts,mts,mts)
        boxints(4) =
     &  qlI4(0d0, 0d0, 0d0, t(i2,i3,i4), s(i2,i3), s(i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s12,mhsq, s34,s123;mts,mts,mts,mts)
        boxints(5) =
     &  qlI4(0d0, 0d0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i3),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s12,mhsq, s34,s124;mts,mts,mts,mts)
        boxints(6) =
     &  qlI4(0d0, 0d0, s(i1,i2), mhsq, s(i3,i4), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s14,mhsq, s23,s124;mts,mts,mts,mts)
        boxints(7) =
     &  qlI4(0d0, 0d0, s(i1,i4), mhsq, s(i2,i3), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s14,mhsq, s23,s134;mts,mts,mts,mts)
        boxints(8) =
     &  qlI4(0d0, 0d0, s(i1,i4), mhsq, s(i2,i3), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s23,mhsq, s14,s123;mts,mts,mts,mts)
        boxints(9) =
     &  qlI4(0d0, 0d0, s(i2,i3), mhsq, s(i1,i4), t(i1,i2,i3),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s23,mhsq, s14,s234;mts,mts,mts,mts)
        boxints(10) =
     &  qlI4(0d0, 0d0, s(i2,i3), mhsq, s(i1,i4), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s34,mhsq, s12,s134;mts,mts,mts,mts)
        boxints(11) =
     &  qlI4(0d0, 0d0, s(i3,i4), mhsq, s(i1,i2), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0,   0, s34,mhsq, s12,s234;mts,mts,mts,mts)
        boxints(12) =
     &  qlI4(0d0, 0d0, s(i3,i4), mhsq, s(i1,i2), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0, s12,   0,mhsq,s123,s124;mts,mts,mts,mts)
        boxints(13) =
     &  qlI4(0d0, s(i1,i2), 0d0, mhsq, t(i1,i2,i3), t(i1,i2,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0, s14,   0,mhsq,s124,s134;mts,mts,mts,mts)
        boxints(14) =
     &  qlI4(0d0, s(i1,i4), 0d0, mhsq, t(i1,i2,i4), t(i1,i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0, s23,   0,mhsq,s123,s234;mts,mts,mts,mts)
        boxints(15) =
     &  qlI4(0d0, s(i2,i3), 0d0, mhsq, t(i1,i2,i3), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

!      (   0, s34,   0,mhsq,s134,s234;mts,mts,mts,mts)
        boxints(16) =
     &  qlI4(0d0, s(i3,i4), 0d0, mhsq, t(i1,i3,i4), t(i2,i3,i4),
     &          mt2, mt2, mt2, mt2, musq, ord)

      end subroutine
