
      double complex function hjetmass_box_pmpm_0_0_s12_mhsq_s34_s123_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = t1 * t3
      t10 = t5 * t6
      t11 = t7 * t8
      t12 = t2 * t4
      t13 = t12 * (t11 + t9 + t10)
      t14 = zb(i4, i1)
      t15 = zb(i4, i2)
      t5 = t14 * t5 + t15 * t7
      t16 = mt ** 2
      t17 = za(i2, i4)
      t18 = za(i1, i4) * t6
      t19 = t17 * t8
      t20 = t19 + t18
      t21 = t12 * t20
      t22 = 16 * t16 * t5 * t21 + 4 * t13 ** 2
      t22 = cdsqrt(t22)
      t13 = 2 * t13
      t23 = -t13 - t22
      t21 = 0.1D1 / t21
      t18 = t19 + t18
      t19 = t12 * t21
      t18 = t19 * (t10 * t18 + t11 * t18 + t9 * t18)
      t24 = (0.1D1 / 0.2D1)
      t20 = t24 * t22 * t21 * t20
      t9 = 2 * t11 + 2 * t9 + 2 * t10
      t10 = -t18 - t20 + t9
      t5 = 0.1D1 / t5
      t11 = t6 * t17
      t25 = t11 / 4
      t26 = t24 * t10 * t5 * t7
      t27 = -t25 * t23 * t21 + t26 * t14
      t26 = t8 * (t23 * t17 * t21 / 4 + t7) - t26 * t15
      t13 = -t13 + t22
      t9 = -t18 + t20 + t9
      t18 = t9 * t5
      t20 = t18 * t24 * t7
      t22 = -t25 * t13 * t21 + t20 * t14
      t8 = t8 * (t17 * t13 * t21 / 4 + t7) - t20 * t15
      t3 = 0.1D1 / t3
      t1 = 0.1D1 / t1
      t15 = 0.1D1 / t7
      t8 = 0.1D1 / t8
      t20 = 0.1D1 / t26
      t21 = t17 ** 2
      t25 = t10 ** 2
      t26 = t27 ** 2
      t28 = t17 * t15
      t29 = t27 * t20
      t30 = t9 * t8
      t31 = t30 * t22
      t32 = t25 * t20
      t33 = t10 * t20
      t34 = t33 * t27
      t35 = t6 * t7
      t36 = t6 ** 2 * t7
      t37 = t3 * t6
      t9 = t9 ** 2 * t8
      t38 = t9 + t32
      t39 = -t35 + t27
      t40 = -t35 + t22
      t41 = t32 * t27
      t42 = t9 * t22
      t43 = t5 ** 2 * t3
      t44 = 0.1D1 / t4
      t45 = t22 * t8
      ret = -t24 * t43 * t11 * t19 * t1 * (t9 * t13 * t40 + t32 * t39 * 
     #t23) - 6 * t43 * t35 * t12 * t1 * (t41 + t42) - 2 * t5 * (t37 * t2
     # * (t30 * (t15 * t22 ** 2 + t36) + t33 * (t26 * t15 + t36)) + (t21
     # * (t35 * (t33 + t30) - t34 - t31) + (t10 * (-t28 * t20 * t26 + t2
     #9 * t11) + t31 * (t22 * (-t28 - t18) + t11) - t32 * t26 * t5) * t3
     # * t2) * t1 * t4) - t43 * t14 * t2 * (t6 * (t35 * t38 - t41 - t42)
     # + (t7 * t5 * (t35 - t27) * t20 * t10 * t25 + t9 * (t17 * t40 + t1
     #8 * t7 * (t35 - t22)) + t32 * t17 * t39) * t1 * t4) + 4 * t37 * (t
     #16 * (t11 * t44 * (t15 * (-t45 - t29) + t6 * (t8 + t20)) + t1 * (t
     #15 * (t45 + t29) + t6 * (-t8 - t20)) * t21) + t5 * t6 * t2 * (t38 
     #* t5 * t1 * t7 ** 2 * t4 + t31 + t34))

      hjetmass_box_pmpm_0_0_s12_mhsq_s34_s123_dp = ret/32d0/(0,1d0)
      return

      end function
