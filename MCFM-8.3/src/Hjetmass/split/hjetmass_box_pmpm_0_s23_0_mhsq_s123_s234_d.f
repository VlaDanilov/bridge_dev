
      double complex function hjetmass_box_pmpm_0_s23_0_mhsq_s123_s234_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i4)
      t2 = za(i2, i3)
      t3 = zb(i3, i2)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = za(i2, i4)
      t7 = zb(i2, i1)
      t8 = zb(i4, i2)
      t9 = za(i1, i3)
      t10 = zb(i3, i1)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t6 * t8
      t14 = t11 * t12
      t15 = t1 * t2 * t3 * t4
      t16 = t9 * t10
      t17 = t1 * t4
      t18 = t17 * ((-t14 - t13) * t7 * t5 + t15 - t16 * (t14 + t13))
      t19 = t10 * t11 + t6 * t7
      t20 = mt ** 2
      t21 = t5 * t8
      t9 = t9 * t12
      t22 = t9 + t21
      t23 = t1 ** 2 * t4 ** 2
      t24 = t23 * t22
      t25 = 16 * t17 * t20 * t19 * t24 + 4 * t18 ** 2
      t25 = cdsqrt(t25)
      t18 = 2 * t18
      t26 = -t18 - t25
      t27 = t14 + t13
      t24 = 0.1D1 / t24
      t28 = 0.1D1 / t4
      t28 = t27 * t28
      t29 = t28 * t10
      t30 = t26 * t24 / 4
      t31 = -t30 * t1 * t12 + t29
      t32 = t9 + t21
      t32 = t5 * t32 * t7 + t16 * t32
      t16 = -t27 * t7 * t5 - t16 * t27 + t15
      t9 = -t23 * t24 * (t13 * t32 + t15 * (-t9 - t21) + t14 * t32) / 2
      t13 = t17 * t25 * t24 * t22 / 4
      t14 = t17 * t19
      t15 = 0.1D1 / t1
      t14 = 0.1D1 / t14
      t15 = t28 * t15 * t5
      t17 = (-t16 + t9 + t13) * t14 * t6
      t19 = t17 - t15
      t21 = t30 * t5
      t22 = -t10 * t19 - t21 * t12
      t19 = -t7 * t19 - t21 * t8
      t18 = -t18 + t25
      t21 = t18 * t24 / 4
      t23 = -t21 * t1 * t12 + t29
      t9 = (-t16 + t9 - t13) * t14 * t6
      t13 = t9 - t15
      t14 = t21 * t5
      t15 = -t10 * t13 - t14 * t12
      t8 = -t7 * t13 - t14 * t8
      t13 = 0.1D1 / t19
      t14 = 0.1D1 / t5
      t12 = 0.1D1 / t12
      t3 = 0.1D1 / t3
      t8 = 0.1D1 / t8
      t2 = 0.1D1 / t2
      t16 = t22 ** 2
      t19 = t31 * t3
      t21 = t10 * t12
      t25 = t22 * t13
      t27 = t16 * t13
      t28 = t1 * t15
      t29 = t15 * t8
      t4 = t2 * t24 * t4
      t24 = t15 ** 2
      t13 = t13 * t22 * t16
      t8 = t24 * t8
      t22 = t8 * t18
      t30 = t27 * t26
      t32 = t29 + t25
      t33 = t10 ** 2
      t34 = t10 * t14
      ret = -3 * t4 * t6 * t12 * (t22 + t30) - t4 * (t26 * (t27 * t1 * t
     #3 * (-t14 * t6 + t21) + t25 * (t6 * (t19 - t6) + t21 * (-t19 + t6)
     # * t5)) + t29 * t18 * (t6 * (t3 * (-t28 * t14 + t23) - t6) + t21 *
     # (t3 * (-t23 * t5 + t28) + t5 * t6))) + 2 * t4 * t3 * t12 * (t30 *
     # t31 + (t26 * (t27 * t17 * t10 + t13) + t22 * (t9 * t10 + t15)) * 
     #t14 * t1 + t23 * t22) + 4 * t2 * t3 * (-t34 * t20 * t32 * t6 ** 2 
     #+ t12 * (t7 * (-t13 * t14 + t8 * (-t14 * t15 + t10) + t27 * t10) +
     # t34 * (t24 + t16) + t33 * t32 * t20) * t6 + t33 * t12 * (t8 + t27
     #) * t11)

      hjetmass_box_pmpm_0_s23_0_mhsq_s123_s234_dp = ret/32d0/(0,1d0)
      return

      end function
