
      double complex function hjetmass_qqbgg_bubble_pmmp_s34_dp 
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
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = zb(i4, i1)
      t5 = zb(i4, i2)
      t6 = t3 * t4
      t7 = t1 * t5
      t8 = t7 + t6
      t9 = za(i1, i2)
      t10 = zb(i2, i1)
      t11 = zb(i3, i2)
      t12 = za(i1, i4)
      t13 = za(i2, i4)
      t14 = t12 * t2
      t15 = t13 * t11
      t16 = t14 + t15
      t17 = t3 * t2
      t18 = t1 * t11
      t19 = t12 * t4
      t5 = t13 * t5
      t20 = -t19 - t5 + t17 + t18
      t20 = t20 ** 2
      t21 = 4 * t16 * t8 + t20
      t21 = cdsqrt(t21)
      t22 = t21 - t19 - t5 + t17 + t18
      t21 = -t21 - t19 - t5 + t17 + t18
      t23 = t9 * t10
      t24 = (-t23 + t19 + t5) * t13
      t25 = t13 * (t18 + t17)
      t26 = t24 - 2 * t25
      t27 = 0.1D1 / t8
      t28 = t27 ** 2
      t29 = t27 * t28
      t30 = 2 * t4
      t31 = t13 * t16
      t32 = 0.1D1 / t22
      t33 = 0.1D1 / t21
      t34 = t32 ** 2
      t35 = t32 * t34
      t36 = t33 ** 2
      t37 = t33 * t36
      t38 = t8 ** 2
      t39 = t38 ** 2
      t40 = t38 * t39
      t41 = t31 * t4
      t42 = t8 * t39 * t34 * t36 * (t32 + t33)
      t31 = t2 * t31 * (384 * t40 * t34 * t36 * (t36 + t34) + 512 * t40 
     #* t35 * t37) * t28
      t39 = 16 * t39
      t40 = t23 + t19
      t43 = t1 * t40
      t44 = 2 * t13 * (t19 + t5) + (t23 - t17 - t18) * t13
      t45 = (t18 + t17) * t27
      t46 = 0.1D1 / t8
      t47 = (-t23 + t19 + t5) * t46
      t48 = (t7 + t6) * t28
      t49 = t21 * t13
      t50 = t23 + t17 + t18
      t51 = t21 ** 2
      t52 = t1 * t51 * t28
      t53 = (0.1D1 / 0.2D1)
      t54 = (0.1D1 / 0.4D1)
      t55 = -t1 * t21 * t51 * t29 * t8 / 8
      t56 = -t53 * t49 * (t48 * t21 + t47) + t54 * t52 * t50 + t13 * (t4
     #5 * t21 + t14 + t15) + t55
      t20 = 4 * t16 * t8 + t20
      t20 = cdsqrt(t20)
      t57 = -t20 - t19 - t5 + t17 + t18
      t57 = 0.1D1 / t57
      t58 = 2 * t16
      t59 = t53 * t21 * t46
      t60 = -t58 * t57 - t59
      t61 = t46 * (-t22 + t21)
      t62 = -t59 * t4 + t2
      t63 = -t53 * t4 * t22 * t46 + t2
      t64 = (-t5 - t19) * t27
      t65 = t22 ** 2
      t66 = t22 * t13
      t29 = -t1 * t22 * t65 * t29 * t8 / 8
      t67 = (0.3D1 / 0.4D1) * t7 * t13
      t68 = -t53 * t66 * (t46 * (t23 - t17 - t18) + t6 * t22 * t28) - t5
     #4 * t43 * t65 * t28 + t13 * (t64 * t22 + t14 + t15) + t29 - t67 * 
     #t65 * t28
      t20 = t20 - t19 - t5 + t17 + t18
      t20 = 0.1D1 / t20
      t69 = t53 * t22 * t46
      t58 = -t58 * t20 - t69
      t65 = t1 * t65 * t28
      t29 = -t53 * t66 * (t48 * t22 + t47) + t54 * t65 * t50 + t29 + t13
     # * (t45 * t22 + t14 + t15)
      t45 = t1 * t22 * t27
      t47 = (0.3D1 / 0.4D1) * t65 * t8
      t48 = -t23 - t17
      t11 = t1 ** 2 * t11
      t50 = (-t7 - t6) * t27
      t65 = t1 * t21 * t27
      t70 = (0.3D1 / 0.4D1) * t52 * t8
      t14 = t53 * t49 * (t46 * (-t23 + t17 + t18) - t6 * t21 * t28) - t5
     #4 * t52 * t40 + t55 + t13 * (t64 * t21 + t14 + t15) - t67 * t51 * 
     #t28
      t9 = 0.1D1 / t9
      t10 = 0.1D1 / t10
      t15 = t32 + t33
      t40 = t13 * t4
      t46 = t1 * t2
      t51 = 0.1D1 / t61
      t52 = 0.1D1 / t60
      t53 = 0.1D1 / t58
      t54 = t46 + t40
      t55 = t57 * t36 * t52
      t58 = t9 * t10
      t60 = t58 * t51 * t16 * t8
      t61 = -t56 + t14
      t64 = t51 ** 2
      t47 = t63 * (-t68 + t29) * t20
      t14 = (t56 - t14) * t52
      t37 = t37 * t52 * t57
      t52 = t58 * t64 * t16
      t3 = t52 * t8 * (t35 * (t47 * t53 ** 2 + (t4 * (t68 - t29) + t63 *
     # (2 * t13 * (t4 * (t3 * t22 * t27 + t12) + t5) + t19 * t45 - t25 +
     # t23 * (t45 + t13) + 3 * t66 * t7 * t27 + t11 * t22 * t27 - t45 * 
     #t48 - t24 + 2 * t13 * (t50 * t22 + t17 + t18))) * t20 * t53) + t37
     # * (t4 * t61 + t62 * (t14 + 2 * t13 * (t4 * (t3 * t21 * t27 + t12)
     # + t5) + t65 * t19 - t25 + t23 * (t65 + t13) + 3 * t49 * t7 * t27 
     #+ t11 * t21 * t27 - t65 * t48 - t24 + 2 * t13 * (t50 * t21 + t17 +
     # t18))))
      ret = -64 * t60 * (t54 * t20 * t53 * t34 - t55 * t54) + 192 * t60 
     #* (t55 * (t1 * t62 - t4 * (-t59 * t1 - t13)) + t20 * t34 * t53 * (
     #-t1 * t63 + t4 * (-t69 * t1 - t13))) + 1024 * t3 + 4096 * t58 * t5
     #1 * t64 * t16 * t8 * (t47 * t35 * t53 + t37 * t61 * t62) - 6144 * 
     #t52 * t38 * (t14 * t36 ** 2 * t57 * t62 + t47 * t34 ** 2 * t53) - 
     #8 * t58 * (-128 * t28 * (t2 * t26 + t41) * t42 + t39 * t28 * (t2 *
     # (2 * (t23 + t17 + t18) * t1 - 4 * t13 * t8) + t30 * t26) * t34 * 
     #t36 + 128 * t28 * (t2 * t44 + t41) * t42 - t39 * t28 * (t2 * (t13 
     #* (-6 * t7 - 4 * t6) - 2 * t43) + t30 * t44) * t34 * t36) - 128 * 
     #t33 * t10 * t9 * t32 * t38 * (t40 * t15 + t46 * t15)

      hjetmass_qqbgg_bubble_pmmp_s34_dp = -ret/16d0*(0,1d0)
      return

      end function
