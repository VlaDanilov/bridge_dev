
      double complex function hjetmass_box_pppm_0_0_s14_mhsq_s23_s124_dp 
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
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t1 * t3
      t10 = t5 * t6
      t11 = t7 * t8
      t12 = t2 * t4
      t13 = t12 * (t11 + t9 + t10)
      t14 = za(i1, i3)
      t15 = za(i3, i4)
      t16 = t14 * t3
      t8 = t15 * t8
      t17 = t8 + t16
      t18 = t12 * t17
      t19 = zb(i3, i1)
      t20 = t1 * t19
      t21 = t7 * zb(i4, i3) + t20
      t22 = mt ** 2
      t23 = 16 * t18 * t22 * t21 + 4 * t13 ** 2
      t23 = cdsqrt(t23)
      t13 = 2 * t13
      t24 = -t13 + t23
      t18 = 0.1D1 / t18
      t8 = t8 + t16
      t8 = t12 * t18 * (t10 * t8 + t11 * t8 + t9 * t8)
      t12 = (0.1D1 / 0.2D1)
      t17 = t12 * t23 * t18 * t17
      t10 = 2 * t11 + 2 * t9 + 2 * t10
      t11 = t17 + t10 - t8
      t21 = 0.1D1 / t21
      t25 = (0.1D1 / 0.4D1)
      t20 = t20 * t12
      t26 = t20 * t11 * t21
      t27 = t16 * t25
      t28 = -t27 * t24 * t18 + t26
      t26 = t3 * (t24 * t14 * t18 * t25 + t1) - t26
      t13 = -t13 - t23
      t8 = -t17 + t10 - t8
      t10 = t25 * t13 * t18
      t17 = -t12 * t8 * t21 * t7 * t19 + t10 * t15 * t3
      t20 = t20 * t8 * t21
      t10 = t3 * (t10 * t14 + t1) - t20
      t23 = -t12 * t11 * t21 * t7 * t19 + t25 * t24 * t18 * t15 * t3
      t20 = -t27 * t13 * t18 + t20
      t27 = -t1 * t15 + t14 * t7
      t29 = t13 * t18 * t8 * t21 * t4 * t27
      t27 = t24 * t18 * t11 * t21 * t4 * t27
      t30 = 0.1D1 / t14
      t31 = 0.1D1 / t5
      t32 = 0.1D1 / t1
      t10 = 0.1D1 / t10
      t6 = 0.1D1 / t6
      t26 = 0.1D1 / t26
      t33 = t7 ** 2
      t34 = t14 * t32
      t35 = t20 * t32
      t36 = -t35 + t3
      t37 = t3 * t7
      t38 = t37 + t23
      t39 = t37 + t17
      t40 = t3 ** 2
      t41 = t28 * t32
      t42 = t15 * t30
      t43 = (-t15 * t7 + t34 * t33) * t40
      t44 = t15 * t28
      t45 = t15 * t17
      t46 = t10 * t13
      t47 = t46 * t39
      t48 = t38 * t26 * t24
      t16 = t16 * t32
      t49 = t6 * t18
      t50 = t49 * t4
      t51 = t7 * t32
      t52 = -t51 - t42
      t53 = -t22 + t20
      t54 = -t22 + t28
      t55 = t32 * t3
      t34 = t34 * t7 * t6
      t56 = t3 * t39
      t57 = t3 * t38
      t58 = t26 * t11
      t59 = t6 * t31
      t34 = t4 * ((t46 * (t34 * t39 * t31 * t19 - t56 * t30 * t6 * t2 + 
     #t42 * t39 - t51 * t39) * t8 + t58 * t24 * (t34 * t38 * t31 * t19 -
     # t57 * t30 * t6 * t2 + t42 * t38 - t51 * t38)) * t21 * t18 * t4 + 
     #t59 * t30 * ((t41 * t38 - t57) * t26 * t27 + t10 * t29 * (t35 * t3
     #9 - t56)))
      t38 = t10 + t26
      t56 = t17 * t10
      t57 = t23 * t26
      t60 = t37 * t38 + t56 + t57
      t61 = t22 * t32
      t62 = -t61 + t3
      t63 = t28 ** 2
      t64 = t20 * t10
      t65 = (-t22 - t20) * t10
      t66 = (-t37 - t23) * t26
      t5 = t6 * t32 * t4 * (t7 * (t66 * t28 - t64 * t39) * t31 * t4 - t2
     #2 * t30 * (t64 * t17 + t57 * t28) + t37 * (t19 * (t26 * (-t22 - t2
     #8) + t65) + t30 * (-t10 * t20 ** 2 - t26 * t63))) + t4 * (t6 * (t1
     #9 * (t3 * (t57 + t56) + t32 * (t65 * t17 - t57 * (t22 + t28)) + t3
     #8 * t7 * t40) + t30 * (-t57 * t63 * t32 + t64 * (t17 * t36 + t37 *
     # t62) + t3 * (t62 * t7 + t23) * t26 * t28)) + t59 * t37 * t60 * t4
     # + t61 * t30 * t60 * t5)
      t38 = -t9 + t28 + t22
      t22 = -t9 + t20 + t22
      t27 = t55 * t50 * t31 * (t48 * t27 + t47 * t29)
      t2 = t49 * t4 ** 2 * (t21 ** 2 * t30 * t19 * t2 * (t48 * t11 ** 2 
     #+ t47 * t8 ** 2) + t16 * (-t13 ** 2 * t10 * t39 + t66 * t24 ** 2) 
     #* t31 * t18 * t15)
      ret = t12 * t34 + t25 * t2 + t27 / 8 - t50 * ((t46 * (t3 * (t51 * 
     #(t14 * t17 + t15 * t20) - t45) + t17 * (t20 * t52 + t37) * t21 * t
     #8) + (t3 * (-t15 * t23 + t51 * (t14 * t23 + t44)) + t23 * (t28 * t
     #52 + t37) * t21 * t11) * t26 * t24) * t31 * t4 + t55 * (t46 * (t17
     # * t53 + t37 * t53) + (t23 * t54 + t37 * t54) * t26 * t24)) - t50 
     #* (t16 * t19 * (t48 + t47) + ((t44 * t23 * t32 + t43 + t3 * (t33 *
     # (-t41 + t3) + t42 * (t1 * t23 + t7 * (t9 - t28))) * t21 * t11) * 
     #t26 * t24 + t46 * (t35 * t45 + t43 + t3 * (t33 * t36 + t42 * (t1 *
     # t17 + t7 * (t9 - t20))) * t21 * t8)) * t31 * t4) + 4 * t5 + 2 * t
     #6 * t21 * t30 * t19 * t4 * (t10 * t8 * (t17 * t22 + t37 * t22) + t
     #58 * (t23 * t38 + t37 * t38))

      hjetmass_box_pppm_0_0_s14_mhsq_s23_s124_dp = ret/32d0/(0,1d0)
      return

      end function
