
      double complex function hjetmass_triangle_pmpm_s14_s134_0_rat_dp 
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
      t2 = zb(i3, i2)
      t3 = za(i2, i3)
      t4 = zb(i3, i1)
      t5 = zb(i4, i1)
      t6 = t1 * t5 + t3 * t4
      t7 = za(i1, i2)
      t8 = zb(i4, i2)
      t9 = t3 * t2
      t10 = t1 * t8
      t11 = t9 + t10
      if ( dreal(t11) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t11 ** 2) + t10 + t9
      t11 = za(i3, i4)
      t12 = zb(i2, i1)
      t13 = zb(i4, i3)
      t14 = t4 * za(i1, i3) + t5 * za(i1, i4)
      t15 = t7 * t11 * t12 * t13
      t16 = -t9 * t14 / 2 + t15
      t17 = -t9 / 2
      t18 = t3 * t12 * t13
      t19 = t11 * (t17 * t5 - t18)
      t20 = t1 * t12
      t21 = t20 * t13
      t17 = t11 * (t17 * t4 + t21)
      t9 = 0.1D1 / t9
      t21 = t11 * (2 * t21 * t9 - t4)
      t14 = -2 * t15 * t9 + t14
      t5 = t11 * (2 * t18 * t9 + t5)
      t8 = (t12 * t14 - t2 * t5 - t21 * t8) * t7
      t8 = 0.1D1 / t8
      t5 = 0.1D1 / t5
      t15 = t6 ** 2
      t18 = 0.1D1 / t19
      t22 = 0.1D1 / t7
      t14 = 0.1D1 / t14
      t23 = 0.1D1 / t12
      t16 = 0.1D1 / t16
      t24 = 0.1D1 / t11
      t25 = 0.1D1 / t13
      t19 = -0.1D1 / t19
      ret = -64 * t20 * t5 * t9 ** 2 * t8 * t15 * t2 * (t3 * t21 * (t7 *
     # t2 * t8 - t5) + t1) - 16 * t6 * (t1 * t4 * t21 * t24 * t25 * t5 *
     # (t22 * t23 * t14 - t8) + (-t17 ** 2 * t23 * t24 * t25 * t19 * t16
     # ** 2 + t22 * t1 * (t18 * (-t3 * t17 * t18 - t1) + t17 * t23 * t24
     # * t25 * t19) * t16) * t6 * t2) - 32 * t5 * t9 * t25 * t24 * t8 * 
     #t21 * t15 * t2 * (t7 ** 2 * t12 * t21 * t8 + t1 - t20 * t7 * t11 *
     # t13 * (-2 * t10 * t9 + 1) * t8)

      hjetmass_triangle_pmpm_s14_s134_0_rat_dp = ret/32d0/(0,1d0)
      return

      end function
