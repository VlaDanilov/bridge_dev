
      double complex function hjetmass_box_pmpm_0_0_s14_mhsq_s23_s124_dp 
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
      t17 = t15 * t8
      t18 = t17 + t16
      t19 = t12 * t18
      t20 = zb(i3, i1)
      t21 = zb(i4, i3)
      t22 = t1 * t20
      t23 = t21 * t7 + t22
      t24 = mt ** 2
      t25 = 16 * t19 * t24 * t23 + 4 * t13 ** 2
      t25 = cdsqrt(t25)
      t13 = 2 * t13
      t26 = -t25 - t13
      t19 = 0.1D1 / t19
      t27 = t17 + t16
      t28 = t12 * t19
      t27 = t28 * (t10 * t27 + t11 * t27 + t9 * t27)
      t18 = t25 * t19 * t18 / 2
      t9 = 2 * t11 + 2 * t9 + 2 * t10
      t10 = -t27 + t9 - t18
      t23 = 0.1D1 / t23
      t22 = t22 / 2
      t29 = t22 * t10 * t23
      t30 = t16 / 4
      t31 = -t30 * t26 * t19 + t29
      t13 = t25 - t13
      t9 = -t27 + t9 + t18
      t18 = t22 * t9 * t23
      t22 = -t30 * t13 * t19 + t18
      t25 = t3 * (t1 + t26 * t14 * t19 / 4) - t29
      t18 = t3 * (t1 + t14 * t13 * t19 / 4) - t18
      t27 = t26 * t19 * t15 * t3 / 4 - t10 * t23 * t7 * t20 / 2
      t29 = -t9 * t23 * t7 * t20 / 2 + t13 * t19 * t15 * t3 / 4
      t6 = 0.1D1 / t6
      t25 = 0.1D1 / t25
      t30 = 0.1D1 / t5
      t1 = 0.1D1 / t1
      t21 = 0.1D1 / t21
      t32 = 0.1D1 / t3
      t18 = 0.1D1 / t18
      t33 = t26 ** 2
      t34 = t26 * t33
      t35 = t13 ** 2
      t36 = t13 * t35
      t37 = t35 + t33
      t38 = -t35 - t33
      t39 = t19 ** 2
      t40 = t4 ** 2
      t9 = t35 * t9
      t10 = t33 * t10
      t41 = t1 * t37
      t42 = t30 * t37
      t43 = t42 * t15
      t44 = t1 * t2
      t45 = t20 * t35
      t46 = t36 * t18
      t47 = t34 * t25
      t48 = t2 * t19
      t49 = t20 * t21
      t50 = t3 * t7
      t51 = t7 ** 2
      t52 = t7 * t1
      t53 = t35 * t29
      t54 = t33 * t27
      t55 = t54 + t53
      t56 = t47 + t46
      t57 = t33 * t25
      t58 = t35 * t18
      t59 = t58 + t57
      t60 = t36 + t34
      t61 = t2 ** 2
      t62 = t20 ** 2
      t63 = t7 * t30
      t53 = t39 * t4 * (t50 * t61 * t19 * t1 * t21 * t56 * t4 + t49 * t2
     # * (-t53 * t18 - t54 * t25 + (t48 * t60 * t3 + t63 * (t22 * t35 + 
     #t31 * t33)) * t1 * t6) - t15 * t2 * t21 * t59 * t62 + t52 * (t2 * 
     #(-t58 * (t50 + t29) + t57 * (-t50 - t27)) - t43 * t7) + t49 * t44 
     #* t11 * t6 * t30 * t55 * t32)
      t54 = t26 + t13
      t58 = -t26 - t13
      t31 = t26 * t31
      t22 = t13 * t22
      t64 = t11 * t54
      t65 = t26 * t27
      t66 = t13 * t29
      t67 = t66 + t65
      t68 = t20 * t6
      t69 = t49 * t19
      t70 = t69 * t1
      t11 = t28 * t7 * (-t45 * t19 * t3 * t18 * t21 + t70 * t59 * t24 + 
     #t1 * (t32 * (-t68 * t67 + (t58 * t20 * t7 - t19 * t55) * t21 * t4)
     # + t69 * t11 * t6 * t37) * t30) + t19 * t20 * t7 * (t28 * t21 * (t
     #38 * t30 * t6 - t57) * t3 + t1 * t6 * (-t20 * t2 * t54 + t63 * (t1
     #2 * t58 - t24 * t54)) + (t4 * (t30 * (t52 * (t22 + t31) + t65 + t6
     #6) + t44 * ((t64 - t31 - t22) * t30 * t6 - t26 - t13) * t20) + t64
     # * t20 * t6 * t30 * t1 * t24 + t44 * t8 * t62 * t6 * t54) * t32 * 
     #t21)
      t22 = t26 * t25
      t28 = t13 * t18
      t31 = t28 + t22
      t55 = t20 * t24
      t22 = t19 * ((t20 * t62 * (t5 * (-t28 - t22) + t58 * t6) - t62 * t
     #31 * t7 * t4 - t42 * t19 * t1 * t51 * t40) * t21 * t2 + t51 * t1 *
     # (-t55 * t31 + (t32 * t67 + t54 * t7) * t30 * t4))
      t28 = t9 + t10
      t31 = t15 * t3
      t5 = t12 * t1 * t39 * (t7 * (t19 * (t31 * t56 + t49 * (t16 * t56 +
     # (-t16 * t60 + t17 * t60) * t30 * t6)) + t63 * t20 * (-t28 * t6 + 
     #(t28 * t4 + t68 * (-t9 - t10) * t8) * t32 * t21) * t23) + t48 * t4
     #9 * t3 * t56 * t5)
      t12 = 0.1D1 / t2
      t16 = (0.1D1 / 0.8D1)
      ret = t16 * (t49 * t39 * (t4 * (-t42 * t14 * t51 * t1 + (-t52 * t4
     # * t33 * t30 * t32 + (-t47 - t46) * t19 * t3) * t15 * t2) + t7 * t
     #2 * (t9 * t4 * t30 * t23 + t41 * t8) * t6 * t20) + t39 * (-t50 * t
     #44 * t20 * t6 * t37 + t44 * (t48 * (t47 * t27 + t46 * t29) + (t19 
     #* (-t36 - t34) - t45 * t32) * t30 * t15 * t7) * t21 * t40 + t49 * 
     #t7 * (t2 * (t41 + t20 * (t38 * t1 * t14 + t10 * t23) * t30 * t6 + 
     #(t10 * t25 + t9 * t18) * t23 * t20) + t43 + t44 * t43 * t20 * t6 *
     # t32 * t8) * t4)) + t53 / 4 + t5 / 16 - (0.3D1 / 0.2D1) * t70 * t2
     #4 * t51 * t4 * t30 * t32 * (t26 + t13) + t31 * t39 ** 2 * t61 * t4
     #0 * t1 * t21 * (t18 * t35 ** 2 + t25 * t33 ** 2) / 32 + t11 / 2 + 
     #t22 - 4 * t63 * t24 * t62 * t32 * t6 * (-t49 + t52) + 8 * t55 * t7
     # * t51 * t1 * t30 * t12 * t32

      hjetmass_box_pmpm_0_0_s14_mhsq_s23_s124_dp = ret/32d0/(0,1d0)
      return

      end function
