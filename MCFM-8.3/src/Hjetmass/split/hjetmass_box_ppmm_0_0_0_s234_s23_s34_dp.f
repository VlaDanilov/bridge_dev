
      double complex function hjetmass_box_ppmm_0_0_0_s234_s23_s34_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i2, i3)
      t2 = za(i3, i4)
      t3 = zb(i3, i2)
      t4 = zb(i4, i3)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = mt ** 2
      t8 = t2 * t4
      t9 = t8 * t1 * t3
      t10 = t9 * (4 * t7 * t5 * t6 + t9)
      t10 = cdsqrt(t10)
      t11 = t10 + t9
      t12 = zb(i2, i1)
      t13 = 0.1D1 / t1
      t14 = 0.1D1 / t3
      t1 = cdsqrt(-t1 * t3* t8 * (-8 * t7 * t5 * t6 - 2 * t9)) * dsqrt(0
     #.2D1) * t13 * t14 / 2
      t3 = t8 - t1
      t15 = zb(i3, i1)
      t6 = 0.1D1 / t6
      t4 = 0.1D1 / t4
      t16 = t3 * t6 * t12
      t17 = t11 * t13 * t14 * t4 * t15
      t18 = t16 + t17
      t19 = za(i1, i4)
      t20 = za(i1, i2)
      t21 = za(i1, i3)
      t5 = 0.1D1 / t5
      t22 = 0.1D1 / t2
      t16 = t16 * t22 * t21 / 2 + t17 * t5 * t20 / 2
      t17 = t20 * t12
      t23 = t16 + t17
      t24 = zb(i4, i1)
      t9 = t10 - t9
      t1 = t8 + t1
      t8 = t9 * t13
      t10 = -t8 * t5 * t14 * t4 * t20 * t15 / 2 + t1 * t22 * t6 * t21 * 
     #t12 / 2
      t17 = t17 + t10
      t8 = t8 * t14 * t4 * t15 - t1 * t6 * t12
      t20 = t15 * t21 + t19 * t24
      t16 = -t20 + t16
      t10 = -t20 + t10
      t10 = 0.1D1 / t10
      t20 = 0.1D1 / t23
      t21 = 0.1D1 / t24
      t17 = 0.1D1 / t17
      t19 = 0.1D1 / t19
      t16 = 0.1D1 / t16
      t22 = (-t20 + t16) * t11
      t23 = (t10 - t17) * t9
      t24 = t22 * t18 + t23 * t8
      t25 = -t10 + t17
      t26 = t20 - t16
      t27 = t9 ** 2
      t28 = t9 * t27
      t29 = t13 ** 2
      t30 = t12 ** 2
      t31 = t4 ** 2
      t32 = t2 ** 2
      t33 = t11 ** 2
      t34 = t11 * t33
      t35 = t14 ** 2
      t36 = t27 * t10
      t37 = t36 * t1
      t38 = t33 * t16
      t39 = t38 * t3
      t33 = t33 * t20
      t40 = t33 * t18
      t27 = t27 * t17
      t41 = t27 * t8
      t42 = t34 * t20
      t43 = t28 * t17
      t44 = t21 * t14
      t45 = t15 ** 2
      t46 = t9 * t10
      t47 = t11 * t16
      t20 = t11 * t20
      t48 = t20 * t18
      t17 = t9 * t17
      t49 = t17 * t8
      t50 = t12 * t19
      t51 = t50 * t13
      t52 = t2 * t4
      t53 = t5 ** 2
      t9 = t52 * t44 * t19 * t13 * (t7 * (t30 * (-t17 + t20) + t5 * ((t2
     #5 * t9 + t22) * t15 * t2 + t48 + t49) * t12) + t52 * t53 * t45 * (
     #t38 + t36))
      ret = t4 * t5 * (t44 * t30 * t24 + t21 * t35 * t31 * t15 * (t12 * 
     #(-t43 + t42) + t5 * (t26 * t18 * t34 + t28 * t8 * t25) + t5 * (t43
     # - t42) * t15 * t2) * t19 * t29 + (t24 * t32 + (t12 * (-t41 + t40)
     # + ((t39 + t37) * t6 * t12 - t40 + t41) * t5 * t15 * t2) * t21 * t
     #14 * t4) * t19 * t13) + 2 * t52 * (t51 * (t17 - t20) * t2 + t5 * (
     #t14 * (((t12 * (-t33 - t27) + t5 * (t38 * t18 - t36 * t8)) * t4 * 
     #t15 + t7 * t12 * (t47 * t18 + t46 * t8)) * t19 * t13 + t15 * t30 *
     # (t11 * t26 + t23)) + t45 * t29 * t5 * t31 * t19 * (-t10 * t28 + t
     #16 * t34) * t35 - t50 * (t49 + t48) + t30 * t19 * (-t17 * t1 + t20
     # * t3) * t6) * t21 + t13 * t15 * t5 * t19 * (-t47 + t46) * t32) - 
     #4 * t9 - t51 * t44 * t53 * t31 * t6 * (t39 * t18 - t37 * t8) / 2

      hjetmass_box_ppmm_0_0_0_s234_s23_s34_dp = ret/32d0/(0,1d0)
      return

      end function
