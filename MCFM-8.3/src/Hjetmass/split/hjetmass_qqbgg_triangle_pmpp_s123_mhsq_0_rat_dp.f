
      double complex function hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_rat_dp
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i2, i4)
      t2 = zb(i4, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i3)
      t5 = zb(i4, i3)
      t6 = t2 * t3 - t4 * t5
      t7 = za(i1, i4)
      t8 = zb(i4, i2)
      t9 = za(i3, i4)
      t10 = t7 * t2
      t11 = t1 * t8
      t12 = t9 * t5
      t13 = t12 + t10 + t11
      if ( dreal(t13) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t13 ** 2) + t10 + t11 + t12
      t12 = zb(i2, i1)
      t13 = zb(i3, i1)
      t14 = za(i1, i3)
      t15 = zb(i3, i2)
      t16 = t3 * t12
      t17 = t14 * t13
      t18 = t15 * t4 + t16 + t17
      t11 = 0.1D1 / t11
      t19 = -0.2D1 * t1 * t18 * t2 * t11 + t13 * t4
      t20 = -0.2D1 * t9 * t18 * t2 * t11 - t12 * t4
      t10 = -0.2D1 * t10 * t18 * t11 + t16 + t17
      t16 = t7 * (t10 * t2 + t19 * t8)
      t17 = -0.2D1 * t7 * t18
      t21 = t2 * (t1 * (t17 * t8 * t11 + t14 * t15) + t10 * t7)
      t15 = t17 * t5 * t11 - t15 * t3
      t17 = 0.1D1 / t21
      t9 = 0.1D1 / t9
      t12 = 0.1D1 / t12
      t20 = 0.1D1 / t20
      t16 = 0.1D1 / t16
      t21 = 0.1D1 / t4
      t22 = 0.1D1 / t15
      t23 = 0.1D1 / t3
      t24 = t2 * t17
      t18 = 0.1D1 / t18
      ret = -0.32D2 * t23 * t11 * t21 * t12 * (t6 * (t2 * (t24 * t15 - t
     #9) + (-t22 * t9 + t24) * t10 * t5) * t1 + t2 * t6 * (t5 * t7 * t16
     # - t20) * t19 + t10 * (t5 * (t14 * t2 + t4 * t8) * t9 ** 2 * t22 +
     # t2 ** 2 * t15 * t17 ** 2 * (t2 * (t14 * t5 + t3 * t8) - t6 * t8))
     # * t1 ** 2) - 0.128D3 * t13 * t6 * t11 * t12 * t18 * (t2 * t21 - t
     #23 * t5)

      hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0_rat_dp = ret/32d0*(0,1d0)
      return

      end function
