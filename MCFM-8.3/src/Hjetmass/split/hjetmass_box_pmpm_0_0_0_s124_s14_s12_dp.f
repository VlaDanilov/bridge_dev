
      double complex function hjetmass_box_pmpm_0_0_0_s124_s14_s12_dp 
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
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = mt ** 2
      t8 = t1 * t2
      t9 = t3 * t4
      t10 = 0.1D1 / t4
      t11 = 0.1D1 / t3
      t12 = (0.1D1 / 0.2D1)
      t13 = t12 *cdsqrt(-t9 * t8 * (-2 * t8 * t3 * t4 - 8 * t7 * t5 * t6
     #)) * dsqrt(0.2D1) * t11 * t10
      t14 = t13 + t8
      t15 = zb(i3, i1)
      t16 = zb(i4, i3)
      t17 = za(i2, i3)
      t18 = zb(i3, i2)
      t19 = t9 * t8
      t20 = t19 * (4 * t7 * t5 * t6 + t19)
      t20 = cdsqrt(t20)
      t21 = t20 - t19
      t22 = za(i3, i4)
      t23 = za(i1, i3)
      t1 = 0.1D1 / t1
      t2 = 0.1D1 / t2
      t6 = 0.1D1 / t6
      t24 = 0.1D1 / t5
      t25 = t21 * t11
      t26 = t22 * t16
      t27 = t12 * (t25 * t24 * t2 * t10 * t22 * t15 - t14 * t1 * t6 * t2
     #3 * t16)
      t28 = -t27 + t26
      t29 = t14 * t6
      t25 = -t25 * t2 * t10 * t15 + t29 * t16
      t8 = -t13 + t8
      t13 = t20 + t19
      t11 = t13 * t11
      t19 = t12 * (t11 * t24 * t2 * t10 * t22 * t15 + t8 * t1 * t6 * t23
     # * t16)
      t20 = t19 + t26
      t22 = t8 * t6
      t10 = t11 * t2 * t10 * t15 + t22 * t16
      t11 = t15 * t23 + t17 * t18
      t19 = t19 - t11
      t11 = -t11 - t27
      t17 = 0.1D1 / t17
      t19 = 0.1D1 / t19
      t18 = 0.1D1 / t18
      t11 = 0.1D1 / t11
      t23 = t6 ** 2
      t24 = t14 ** 2
      t26 = t14 * t24
      t27 = t8 ** 2
      t30 = t8 * t27
      t31 = t27 * t19
      t32 = t31 * t13
      t33 = t11 * t24
      t34 = t33 * t21
      t20 = 0.1D1 / t20
      t28 = 0.1D1 / t28
      t35 = -t20 + t19
      t36 = (-t28 + t11) * t25
      t37 = t8 * t10
      t38 = t36 * t14 + t37 * t35
      t39 = t20 * t30 + t26 * t28
      t40 = t5 * t16 * t1
      t41 = -t40 + t15
      t42 = t15 ** 2
      t43 = t5 ** 2
      t27 = t27 * t20
      t24 = t24 * t28
      t44 = t42 * t17
      t45 = t14 * (t28 - t11) + t8 * (t20 - t19)
      t46 = t14 * t25
      t47 = t14 * t11
      t14 = t14 * t28
      t48 = t8 * t20
      t49 = t18 * t15
      t9 = t49 * t6 * t17 * (t7 * (t5 * (t1 * (-t37 * t20 - t46 * t28) +
     # t15 * t45) - t16 * t1 * (t48 + t14) * t43) - t9 * t6 * t15 * (t33
     # + t31))
      ret = t12 * t1 ** 2 * t23 * t15 * t2 * t17 * t18 * (-t32 * t10 + t
     #34 * t25) - t6 * t1 * (t17 * t43 * t4 * t38 + t18 * (t38 * t42 + (
     #t6 * (t27 * t10 * t41 + t24 * t41 * t25) + t16 * (t1 * (t30 * t10 
     #* t35 + t36 * t26) + t15 * t39 - t40 * t39) * t23) * t17 * t4) * t
     #3 + t44 * t6 * t18 * (t34 - t32) * t2) - 2 * t6 * (t17 * (t43 * (-
     #t15 * (t19 * t8 + t47) + t40 * (-t48 - t14)) + t1 * t15 * (t6 * (t
     #5 * (t24 + t27) + t6 * (-t11 * t26 - t19 * t30)) * t16 + t37 * (-t
     #22 * t19 + t20 * t5) + t46 * (-t29 * t11 + t28 * t5)) * t18 * t3) 
     #* t4 + t44 * t5 * t1 * t18 * (-t48 * t13 + t14 * t21) * t2 + t49 *
     # (t45 * t42 * t3 - t7 * t1 * t5 * t17 * (t37 * t19 + t47 * t25))) 
     #+ 4 * t9

      hjetmass_box_pmpm_0_0_0_s124_s14_s12_dp = ret/32d0/(0,1d0)
      return

      end function
