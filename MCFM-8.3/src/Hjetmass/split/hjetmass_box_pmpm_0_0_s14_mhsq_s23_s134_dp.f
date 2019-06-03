
      double complex function hjetmass_box_pmpm_0_0_s14_mhsq_s23_s134_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i3, i4)
      t6 = zb(i4, i3)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = za(i2, i4)
      t10 = za(i1, i2) * t2
      t11 = t9 * t6
      t12 = t11 + t10
      t13 = t7 * t8
      t14 = t13 * t12
      t15 = zb(i2, i1)
      t16 = zb(i4, i2)
      t17 = t1 * t15 + t16 * t5
      t18 = mt ** 2
      t1 = t1 * t2
      t19 = t3 * t4
      t20 = t5 * t6
      t21 = t13 * (t20 + t1 + t19)
      t22 = 16 * t18 * t17 * t14 + 4 * t21 ** 2
      t22 = cdsqrt(t22)
      t14 = 0.1D1 / t14
      t23 = t20 + t1
      t10 = t13 * t14 * (t10 * t23 + t11 * t23 + t19 * (t11 + t10))
      t11 = 2
      t23 = (0.1D1 / 0.2D1)
      t1 = t11 * (t20 + t1 + t19)
      t12 = t23 * t22 * t14 * t12
      t19 = t1 + t12 - t10
      t20 = t11 * t21
      t21 = -t20 + t22
      t17 = 0.1D1 / t17
      t24 = t21 * t14 * t9 / 4
      t25 = t23 * t19 * t17 * t5
      t26 = -t25 * t15 + t24 * t2
      t1 = t1 - t12 - t10
      t10 = -t20 - t22
      t12 = t10 * t14 * t9 / 4
      t20 = t23 * t1 * t17 * t5
      t22 = t12 * t2 - t20 * t15
      t12 = t6 * (t5 + t12) - t20 * t16
      t16 = t6 * (t5 + t24) - t25 * t16
      t20 = 0.1D1 / t4
      t6 = 0.1D1 / t6
      t24 = 0.1D1 / t5
      t16 = 0.1D1 / t16
      t3 = 0.1D1 / t3
      t12 = 0.1D1 / t12
      t25 = 0.1D1 / t15
      t27 = t22 ** 2
      t28 = t26 ** 2
      t29 = t1 ** 2
      t30 = t9 ** 2
      t31 = t19 ** 2
      t32 = t8 ** 2
      t33 = t2 ** 2
      t34 = t2 * t33
      t35 = t26 * t19
      t36 = t1 * t22
      t37 = t12 * t29
      t38 = t31 * t16
      t39 = t7 * t17
      t40 = t7 * t34 * t20
      t41 = t17 * t6
      t4 = t4 * t5
      t42 = t19 * t31 * t16
      t43 = t41 * t8
      t6 = t3 * t6
      t44 = t6 * t20
      t45 = t2 * t5
      t45 = t43 * t25 * t3 * t2 * (t7 * (t24 * t20 * (t1 * t27 + t19 * t
     #28) + (t37 * t22 + t38 * t26) * t17 * t8 * t5) + (-t1 * t12 * (t45
     # + t22) + (-t45 - t26) * t16 * t19) * t9 * t18)
      t6 = t6 * t17 ** 2 * t9 * t32 * t7 * t5 * (t17 * t15 * (t12 * t1 *
     # t29 + t42) + (t37 * t10 + t38 * t21) * t25 * t14 * t33)
      ret = -t11 * t41 * (t40 * (t19 + t1) + t39 * t25 * (t37 * t27 + t3
     #8 * t28) * t3 * t32 + t2 * (t30 * (t1 * (t12 * t18 + 1) + t19 * (t
     #16 * t18 + 1)) + (-t17 * (t22 * t29 + t26 * t31) + (-t36 - t35) * 
     #t24 * t9) * t20 * t7) * t3 * t8) + t23 * t6 + 8 * t44 * t33 * (t18
     # * (-t24 * t30 + t25 * (t26 + t22) * t24 * t9) - t13 * t17 * t25 *
     # (t36 + t35)) - t43 * (t17 * (t7 * t33 * t5 * (t38 + t37) + (-t42 
     #* t39 * t5 * t26 + t37 * (t4 * t30 + (-t1 * t5 * t17 - t9) * t22 *
     # t7) + t38 * t9 * (-t26 * t7 + t4 * t9)) * t3 * t8 - t2 * t7 * t9 
     #* t3 * (t31 + t29) * t20 * t15) - t40 * t9 * t3 * (t1 * t10 + t19 
     #* t21) * t25 * t14) - 4 * t45 + 16 * t44 * t18 * t9 * t34 * t25

      hjetmass_box_pmpm_0_0_s14_mhsq_s23_s134_dp = ret/32d0/(0,1d0)
      return

      end function
