
      double complex function hjetmass_triangle_pppm_s123_0_mhsq_rat_dp 
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
      t2 = zb(i2, i1)
      t3 = za(i3, i4)
      t4 = zb(i3, i1)
      t5 = t1 * t2 + t3 * t4
      t6 = zb(i4, i2)
      t7 = za(i1, i4)
      t8 = zb(i4, i1)
      t9 = zb(i4, i3)
      t10 = t7 * t8
      t11 = t1 * t6
      t12 = t3 * t9
      t13 = t11 + t12 + t10
      if ( dreal(t13) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i2, i3)
      t13 = za(i1, i2)
      t14 = zb(i3, i2)
      t15 = t13 * t2
      t16 = za(i1, i3) * t4
      t17 = t12 * t14 + t15 + t16
      t11 = 0.1D1 / t11
      t10 = -2 * t10 * t17 * t11 + t15 + t16
      t15 = t1 * t14 + t4 * t7
      t16 = -2 * t3 * t17 * t8 * t11 - t12 * t2
      t4 = (-2 * t1 * t17 * t8 * t11 + t12 * t4) * t6
      t17 = t7 * (t10 * t8 + t4)
      t2 = -t14 * t3 + t2 * t7
      t4 = t7 * (t16 * t9 + t4)
      t14 = 0.1D1 / t1
      t13 = 0.1D1 / t13
      t17 = 0.1D1 / t17
      t18 = 0.1D1 / t8
      t10 = 0.1D1 / t10
      t4 = 0.1D1 / t4
      t19 = 0.1D1 / t9
      t12 = 0.1D1 / t12
      t20 = 0.1D1 / t3
      t16 = 0.1D1 / t16
      t21 = 0.1D1 / t7
      t22 = t5 ** 2
      t23 = t17 ** 2
      t24 = t10 ** 2
      t25 = t7 ** 2
      t26 = t15 * t14
      t27 = t2 * t20
      t28 = t6 * t20
      t29 = t6 * t5
      t30 = t9 * t14
      ret = -32 * t11 * t22 * (t7 * (t28 * t17 * t9 * (t22 * t14 * t18 *
     # t24 + (t1 * t13 - t10 * t5) * t17 * t2) + (t12 * (t6 * (-t27 * t1
     # - t15) + t9 * (-t26 * t3 - t2)) - t6 * (t27 + t26) * t16 * t5) * 
     #t4 ** 2 * t8) + t16 ** 2 * t15 * (-t6 * t13 * t21 + (-t29 * t4 - t
     #13) * t14 * t8) * t19 - t30 * t7 * t25 * t5 * t8 * t13 * t20 * t23
     # - t5 * t12 * t18 * t24 * (t28 + t30) + t20 * t23 * t9 * (t22 * t6
     # * t14 * t10 + t8 * t2 * t13 - t29 * t13) * t25)

      hjetmass_triangle_pppm_s123_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
