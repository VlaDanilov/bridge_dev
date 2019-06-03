
      double complex function hjetmass_triangle_pppm_s34_mhsq_s12_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret

          parameter (cg = 1d0)

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t1 * t2
      t10 = t3 * t4
      t11 = t5 * t6
      t12 = t7 * t8
      t13 = t12 + t9 + t10 + t11
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t13 = -4 * t14 * t15 * t16 * t17 + t13 ** 2
      t10 = cg * cdsqrt(t13) + t10 + t11 + t12 + t9
      t13 = t2 * t5 + t4 * t7
      t9 = t9 + t11
      t18 = t14 * t16
      t19 = (0.1D1 / 0.2D1)
      t20 = t18 * t19 * t10 - t18 * t9
      t21 = t19 * t10
      t22 = -t1 * t16 * t17 + t21 * t8
      t23 = t14 * t22
      t24 = t10 ** 2
      t25 = t18 * t15 * t17
      t26 = -t24 / 4 + t25
      t27 = t1 * t4 + t5 * t8
      t28 = t14 * t15
      t1 = -t21 * t1 + t28 * t8
      t8 = t17 * t1
      t9 = t19 * t10 * t9 - t25
      t29 = t5 * t16 * t17 + t21 * t4
      t30 = -t15 * t29
      t31 = t7 * t16 * t17
      t32 = -t21 * t2 + t31
      t33 = t15 * t32
      t34 = t2 * t3 + t6 * t7
      t35 = t15 * (t3 * t16 * t17 + t21 * t6)
      t36 = t28 * t2
      t37 = -t21 * t7 + t36
      t38 = t16 * t37
      t39 = t15 * t17
      t40 = t39 * t13
      t3 = -t16 * (t21 * t3 + t28 * t6)
      t6 = t12 + t11
      t11 = t19 * t10 * t6 - t25
      t4 = t21 * t5 + t28 * t4
      t12 = t17 * t4
      t25 = t40 * t10
      t28 = t16 * t27
      t29 = -t14 * t29
      t6 = t39 * t21 - t39 * t6
      t18 = t18 * t27
      t21 = -t10 / 4
      t36 = t10 * t16 * (t36 * t19 + t21 * t7)
      t2 = t10 * t15 * (t31 * t19 + t21 * t2)
      t19 = 0.1D1 / t35
      t5 = 0.1D1 / t5
      t21 = -0.1D1 / t26
      t31 = 0.1D1 / t34
      t34 = 0.1D1 / t14
      t39 = -0.1D1 / t3
      t41 = 0.1D1 / t27
      t42 = 0.1D1 / t20
      t43 = 0.1D1 / t9
      t44 = -0.1D1 / t23
      t45 = 0.1D1 / t15
      t23 = 0.1D1 / t23
      t7 = 0.1D1 / t7
      t46 = 0.1D1 / t16
      t47 = 0.1D1 / t3
      t8 = 0.1D1 / t8
      t48 = 0.1D1 / t17
      t49 = t33 ** 2
      t21 = t21 ** 2
      t50 = t13 ** 2
      t51 = t19 ** 2
      t52 = t38 ** 2
      t53 = t38 * t52
      t26 = 0.1D1 / t26 ** 2
      t54 = t8 ** 2
      t55 = t38 * t7
      t4 = t16 * t4 * t5
      t56 = t5 * t7 * t21
      t57 = t56 * t52
      t58 = t26 * t34
      t59 = t58 * t7
      t60 = t20 * t30 * t34 * t41 * t46
      t61 = t30 * t5
      t62 = t33 * t7
      t63 = t55 * t21 * t42 + t58
      t64 = t39 ** 2
      t65 = t56 * t53
      t6 = t24 * (t10 * t50 * t27 * t5 * t8 * t63 + (t13 * (t29 * (t58 *
     # (t23 * (t62 + t61) + t20 * t6 * t7 * t23 ** 2) + t44 * t21 * t5 *
     # (-t28 * t6 * t44 + t62 * (-t6 * t18 * t44 - t30) * t43)) - t58 * 
     #t52 * t15 * t22 * t5 * t64 + t61 * t58 * t33 * t19 - t65 * t35 * t
     #34 * t46 * t31 * t47 ** 2) + t27 * (t12 * (-t40 * t11 * t5 * t63 *
     # t54 + t57 * t40 * t42 * t43 * t8) + t65 * t42 ** 2 * t43 * t8 * t
     #12 ** 2)) * t48 * t45)
      t15 = 0.1D1 / t10
      t22 = t36 ** 2
      t2 = t2 ** 2
      t2 = t48 * t45 * ((t4 * t53 * t17 * t37 * t21 * t31 ** 2 * t34 ** 
     #2 * t7 * t47 * t46 ** 2 + t4 * t52 * t40 * t21 * t31 * t34 * t7 * 
     #t47 * t46) * t41 * t9 + t56 * (-t13 * t18 * t29 * t2 * t42 * t44 *
     # t43 + t3 * t33 * t49 * t40 * t31 * t51) + t58 * (t13 * t7 * (t35 
     #* t22 * t64 - t19 * t2) + t19 * (t4 + t55) * t25 * t33))
      t17 = -8
      ret = t17 * (t50 * t24 * (t59 * (t20 * t23 - t8 * t9) + (t62 * t18
     # * t43 + t28) * t44 * t21 * t5) + (t12 * (t59 * t9 * t11 * t25 * t
     #54 + t25 * (-t4 * t55 * t21 * t42 + t58 * (-t55 - t4)) * t8) + t56
     # * t49 * (t33 * t29 ** 2 * t18 * t42 * t43 ** 2 * t44 + (-t60 + t3
     #3) * t31 * t19 * t13) * t10 + t34 * (t38 * (t57 * t31 * t46 * t47 
     #- t4 * t26 * t39 - t55 * t26 * t39) - t26 * (-t1 * t16 * t5 + t3 *
     # t7) * t51 * t49) * t25 + t58 * t13 * t38 * t39 * (t62 + t61) * t2
     #4) * t48 * t45) - 32 * t56 * t48 * t45 * t31 * t33 * (-t60 * t14 *
     # t32 * t49 * t15 * t31 * t19 - t22 * t13 * t34 * t46 * t47 - t38 *
     # t33 * t40 * t19) - 16 * t2 + 4 * t6

      hjetmass_triangle_pppm_s34_mhsq_s12_rat_dp = ret/32d0/(0,1d0)
      return

      end function
