
      double complex function hjetmass_box_pppm_0_0_s14_mhsq_s23_s134_dp 
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
      t2 = za(i2, i3)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = zb(i2, i1)
      t10 = t7 * zb(i4, i2)
      t11 = t1 * t9
      t12 = t11 + t10
      t13 = mt ** 2
      t14 = za(i1, i2)
      t15 = za(i2, i4)
      t16 = t14 * t3
      t17 = t15 * t8
      t18 = t17 + t16
      t19 = t2 * t4
      t20 = t19 * t18
      t21 = t1 * t3
      t22 = t5 * t6
      t23 = t7 * t8
      t24 = t19 * (t23 + t21 + t22)
      t25 = 16 * t13 * t12 * t20 + 4 * t24 ** 2
      t25 = cdsqrt(t25)
      t24 = 2 * t24
      t26 = -t24 + t25
      t20 = 0.1D1 / t20
      t17 = t17 + t16
      t17 = t19 * t20 * (t21 * t17 + t22 * t17 + t23 * t17)
      t19 = (0.1D1 / 0.2D1)
      t18 = t19 * t25 * t20 * t18
      t22 = 2 * t23 + 2 * t21 + 2 * t22
      t23 = t18 + t22 - t17
      t12 = 0.1D1 / t12
      t27 = (0.1D1 / 0.4D1)
      t28 = t15 * t26 * t20 * t27
      t10 = t10 * t19
      t29 = t8 * (t28 + t7) - t10 * t23 * t12
      t30 = t16 * t27
      t11 = t11 * t19
      t31 = t11 * t23 * t12 - t30 * t26 * t20
      t28 = -t19 * t23 * t12 * t7 * t9 + t28 * t3
      t24 = -t24 - t25
      t17 = -t18 + t22 - t17
      t18 = t15 * t24
      t22 = t18 * t20 * t27
      t8 = t8 * (t7 + t22) - t10 * t17 * t12
      t10 = t11 * t17 * t12 - t30 * t24 * t20
      t11 = -t19 * t17 * t12 * t7 * t9 + t22 * t3
      t22 = t1 * t15 - t14 * t7
      t25 = t26 * t20 * t23 * t12 * t4 * t22
      t8 = 0.1D1 / t8
      t29 = 0.1D1 / t29
      t30 = 0.1D1 / t5
      t32 = 0.1D1 / t14
      t33 = 0.1D1 / t1
      t6 = 0.1D1 / t6
      t34 = t3 * t7
      t35 = t34 + t28
      t36 = t34 + t11
      t37 = t35 * t29
      t38 = t6 * t20
      t39 = t10 * t33
      t40 = -t39 + t3
      t41 = t39 - t3
      t42 = t31 * t33
      t43 = -t42 + t3
      t44 = t15 * t31
      t45 = t34 * t33
      t46 = t10 * t32
      t47 = t33 * t36
      t48 = t29 * t26
      t49 = t37 * t26
      t50 = t24 * t11 * t8
      t51 = t9 * t14
      t52 = t33 * t3
      t53 = t38 * t4
      t54 = t3 ** 2
      t55 = t54 * t7
      t56 = t48 * (t13 - t31)
      t1 = t53 * ((t18 * t8 * (t11 * t41 + t3 * (t1 * t11 + t7 * (t21 - 
     #t10)) * t12 * t32 * t17) + t48 * (t3 * t15 * (t1 * t28 + t7 * (t21
     # - t31)) * t12 * t32 * t23 + t15 * (t42 - t3) * t28 - t55 * t15 + 
     #t54 * t7 ** 2 * t14 * t33)) * t30 * t4 + t52 * (t50 * t13 + t34 * 
     #(t56 + (-t51 - t10 + t13) * t8 * t24) + t56 * t28)) + t53 * (t52 *
     # (-t50 * t10 + t51 * (-t50 - t49)) + (t48 * (t45 * (t14 * t28 + t4
     #4) + (t7 * (t28 * t43 + t34 * t43) - t44 * t28 * t32) * t12 * t23)
     # + ((t7 * (t11 * t40 + t34 * t40) - t46 * t15 * t11) * t12 * t17 +
     # t34 * (t47 * t14 + t15 * t41)) * t8 * t24) * t30 * t4)
      t14 = t7 * t33 * t35
      t18 = t35 * t32
      t40 = t47 * t7
      t41 = t36 * t32
      t22 = t8 * t24 * t20 * t17 * t12 * t4 * t22
      t43 = t32 * t6
      t44 = t21 - t31 - t13
      t21 = t21 - t10 - t13
      t47 = t10 * t8
      t50 = t31 * t29
      t54 = t50 + t47
      t56 = t8 + t29
      t57 = t31 ** 2
      t58 = t3 * t9
      t59 = t29 * t28
      t5 = t5 * t32
      t60 = t9 * t6
      t61 = (-t60 + t5) * t13
      t62 = t7 * t4 * t30
      t5 = t4 * (t59 * (t33 * (-t60 * t31 - t43 * t57 + t61) + t6 * (t31
     # * t32 - t62) * t3) + (t33 * (-t43 * t10 ** 2 - t60 * t10 + t61) +
     # t6 * (t46 - t62) * t3) * t8 * t11) + t4 * (t6 * (t34 * (t3 * (t32
     # * t54 + t56 * t9) + t33 * (t32 * (-t29 * t57 - t47 * (t13 + t10) 
     #- t50 * t13) + t9 * (t29 * (-t13 - t31) + t8 * (-t13 - t10)))) + t
     #8 * (-t46 * t13 * t33 + t58) * t11 + t59 * (-t42 * t13 * t32 + t58
     #)) + t5 * t45 * t13 * t56 + t6 * t7 * (t33 * (t47 * t11 + t59 * t3
     #1) + t55 * (-t8 - t29) + t45 * t54) * t30 * t4)
      ret = t19 * t4 * ((t48 * t23 * (t14 * t51 * t6 * t30 + t18 * t3 * 
     #t6 * t2 + t18 * t15 - t14) + (t41 * t3 * t6 * t2 + t40 * t51 * t6 
     #* t30 + t15 * t41 - t40) * t8 * t17 * t24) * t12 * t20 * t4 + t43 
     #* t30 * (t22 * (t3 * t36 + t39 * (-t34 - t11)) + (t3 * t35 + t42 *
     # (-t34 - t28)) * t29 * t25)) - t27 * t38 * t4 ** 2 * (t16 * (t24 *
     #* 2 * t8 * t36 + t37 * t26 ** 2) * t33 * t30 * t20 * t15 + t12 ** 
     #2 * t32 * t9 * t2 * (t24 * t17 ** 2 * t8 * t36 + t37 * t23 ** 2 * 
     #t26)) - t53 * t52 * t30 * (t22 * t24 * t36 + t49 * t25) / 8 - t1 +
     # 2 * t43 * t12 * t9 * t4 * (t8 * t17 * (t11 * t21 + t34 * t21) + (
     #t44 * t28 + t34 * t44) * t29 * t23) - 4 * t5

      hjetmass_box_pppm_0_0_s14_mhsq_s23_s134_dp = ret/32d0/(0,1d0)
      return

      end function
