
      double complex function hjetmass_box_pmpm_0_0_0_s123_s12_s23_dp 
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
      t2 = za(i2, i3)
      t3 = zb(i2, i1)
      t4 = zb(i3, i2)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = mt ** 2
      t8 = t2 * t4
      t9 = t8 * t1 * t3
      t10 = t9 * (4 * t7 * t5 * t6 + t9)
      t10 = cdsqrt(t10)
      t11 = -t9 + t10
      t12 = za(i1, i4)
      t13 = zb(i4, i1)
      t14 = 0.1D1 / t1
      t15 = 0.1D1 / t3
      t16 = (0.1D1 / 0.2D1)
      t1 = t16 *cdsqrt(-t1 * t3 * t8 * (-8 * t7 * t5 * t6 - 2 * t9)) * d
     #sqrt(0.2D1) * t14 * t15
      t3 = t1 + t8
      t17 = za(i2, i4)
      t18 = zb(i4, i2)
      t19 = za(i3, i4)
      t2 = 0.1D1 / t2
      t5 = 0.1D1 / t5
      t20 = 0.1D1 / t6
      t4 = 0.1D1 / t4
      t21 = t3 * t2
      t22 = t11 * t14 * t5 * t15
      t23 = t17 * t18 + t19 * zb(i4, i3)
      t24 = t16 * (-t22 * t4 * t12 * t18 + t21 * t20 * t17 * t13)
      t25 = t24 - t23
      t21 = t21 * t17
      t22 = -t22 * t12 + t21
      t9 = t9 + t10
      t1 = -t1 + t8
      t8 = t1 * t2 * t17
      t10 = t9 * t14
      t26 = t10 * t5
      t27 = t26 * t15 * t12
      t28 = t8 + t27
      t8 = t16 * (t8 * t20 * t13 + t27 * t4 * t18)
      t18 = t8 - t23
      t20 = t12 * t13
      t8 = t8 + t20
      t20 = t24 + t20
      t23 = 0.1D1 / t12
      t18 = 0.1D1 / t18
      t24 = 0.1D1 / t25
      t13 = 0.1D1 / t13
      t25 = t22 * t4
      t27 = -t25 + t17
      t29 = t17 ** 2
      t30 = t5 ** 2
      t31 = t11 ** 2
      t32 = t11 * t31
      t33 = t4 ** 2
      t34 = t28 * t4
      t35 = t6 * t17
      t36 = t17 * t23
      t12 = t6 ** 2 * t12
      t11 = t11 * t24
      t37 = t9 * t18
      t38 = t4 * t13 * t5
      t20 = 0.1D1 / t20
      t8 = 0.1D1 / t8
      t39 = t9 ** 2
      t40 = t9 * t39
      t18 = (-t18 + t8) * t39
      t41 = t31 * t20
      t42 = t2 * t5
      t43 = t15 * t4 * t14
      t44 = t6 * t19
      t45 = t44 * t15
      t26 = t35 * t38 * t2 * t7 * (-t43 * t31 * t5 * t24 - t11 * t23 * (
     #t45 + t17) + t37 * (t15 * (t44 * t23 - t26 * t4) + t36))
      t44 = -t37 + t11
      t1 = -t16 * t14 * t5 * t30 * t4 * t33 * t2 ** 2 * t17 * t15 * t13 
     #* (t40 * t1 * t8 - t32 * t3 * t20) - t43 * t30 * (t13 * (t42 * (t3
     #2 * t22 * t20 - t40 * t28 * t8) * t33 + t12 * (t41 + t18) * t4) + 
     #t29 * t2 * (t31 * (t20 - t24) + t18)) - t38 * (t2 * t30 * t14 * t4
     # * t15 * t27 * t24 * t32 + t11 * t36 * t6 * (-t21 + t22) + (-t12 *
     # t5 * t14 * t4 * t15 - t25 * t36 * t2 * t5 + (t25 * t5 * t14 + (t2
     #7 * t5 * t19 - t14 * t29) * t23 * t2) * t15 * t6) * t24 * t31 + t3
     #7 * (t36 * (t2 * (-t34 * t9 * t5 + t35 * t1) - t28 * t6) + t10 * t
     #15 * (t30 * (-t17 * t9 * t2 * t4 + t9 * t28 * t2 * t33) - t6 * t29
     # * t2 * t23 + t34 * t6 * t5) + t6 * t9 * t2 * t5 * t23 * t15 * (-t
     #34 + t17) * t19)) + 4 * t26 + 2 * t38 * t6 * (t6 * (t17 * t44 + t4
     #5 * t44) + t43 * t42 * t17 * (t39 * t8 + t41) * t7)
      ret = t1

      hjetmass_box_pmpm_0_0_0_s123_s12_s23_dp = ret/32d0/(0,1d0)
      return

      end function
