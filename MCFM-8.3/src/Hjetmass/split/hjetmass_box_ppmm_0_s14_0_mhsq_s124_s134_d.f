
      double complex function hjetmass_box_ppmm_0_s14_0_mhsq_s124_s134_dp 
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
      t2 = za(i3, i4)
      t3 = zb(i2, i1)
      t4 = za(i1, i4)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = za(i1, i2)
      t10 = zb(i3, i1)
      t11 = za(i2, i4)
      t12 = zb(i4, i3)
      t13 = t1 * t3 + t2 * t8
      t14 = mt ** 2
      t15 = t9 * t10
      t16 = t11 * t12
      t17 = t15 + t16
      t18 = t6 ** 2
      t19 = t5 ** 2 * t18
      t20 = t19 * t17
      t21 = t3 * t9
      t22 = t8 * t11
      t23 = t22 + t21
      t24 = t4 * t5
      t25 = t24 * t6 * t7
      t26 = t2 * t12
      t27 = t1 * t10
      t28 = t5 * t6
      t23 = t28 * (t26 * t23 + t27 * t23 - t25)
      t29 = 16 * t14 * t5 * t6 * t13 * t20 + 4 * t23 ** 2
      t29 = cdsqrt(t29)
      t23 = 2 * t23
      t30 = t23 + t29
      t31 = 0.1D1 / t6
      t20 = 0.1D1 / t20
      t32 = (0.1D1 / 0.4D1)
      t33 = (t26 + t27) * t31
      t34 = t33 * t3
      t35 = t32 * t30 * t20
      t36 = t35 * t5
      t37 = t36 * t10 - t34
      t38 = t26 + t27
      t38 = t15 * t38 + t16 * t38
      t39 = t22 + t21
      t26 = t26 * t39 + t27 * t39 - t25
      t16 = t19 * t20 * (t21 * t38 + t22 * t38 + t25 * (-t16 - t15)) / 2
      t17 = t32 * t28 * t29 * t20 * t17
      t13 = t28 * t13
      t13 = 0.1D1 / t13
      t19 = 0.1D1 / t5
      t22 = t33 * t19
      t25 = t22 * t9
      t27 = (-t26 + t16 + t17) * t13
      t38 = t3 * (t27 * t1 + t25) - t15 * t35
      t23 = t23 - t29
      t22 = t22 * t11
      t13 = (-t26 + t16 - t17) * t13
      t16 = t32 * t23 * t20
      t17 = -t3 * (t13 * t2 + t22) + t16 * t11 * t10
      t13 = t3 * (t13 * t1 + t25) - t16 * t15
      t15 = -t3 * (t27 * t2 + t22) + t35 * t11 * t10
      t22 = t33 * t8
      t16 = t16 * t5
      t25 = -t16 * t12 + t22
      t12 = -t36 * t12 + t22
      t10 = t16 * t10 - t34
      t16 = 0.1D1 / t1
      t22 = 0.1D1 / t8
      t26 = 0.1D1 / t4
      t13 = 0.1D1 / t13
      t27 = 0.1D1 / t38
      t7 = 0.1D1 / t7
      t29 = t23 * t10
      t33 = t30 * t37
      t34 = t33 + t29
      t35 = t30 * t27
      t36 = t23 * t13
      t38 = t36 + t35
      t39 = t2 ** 2
      t40 = t2 * t39
      t41 = t3 ** 2
      t42 = t14 * t2
      t43 = t16 * t7
      t44 = t43 * t2
      t45 = (t12 * t30 + t23 * t25) * t11
      t46 = t26 * t7
      t47 = t46 * t2
      t48 = t36 * t10
      t49 = t36 * t17
      t50 = t15 * t35
      t51 = t39 * t8 * t16
      t52 = t39 * t16
      t53 = t16 * t2
      t54 = t35 * t37 + t48
      t55 = t30 + t23
      t56 = t30 ** 2
      t57 = t30 * t56
      t58 = t23 ** 2
      t59 = t23 * t58
      t60 = t11 * t19
      t61 = t3 * t25 * t22
      t62 = t58 * t10
      t63 = t56 * t37
      t64 = t7 * t22 * t3
      t65 = t16 * t6
      t24 = t20 * (-t65 * t26 * t55 * t40 + t64 * (t41 * t55 + (-t63 - t
     #62) * t20 * t16 * t6) * t5 + t52 * (t11 * (-t6 * t58 * t20 * t26 +
     # t54 * t19 * t8) + t21 * t19 * (t33 * (t46 + t27) + t29 * (t46 + t
     #13))) + t3 * (t30 * (t27 * (-t11 * t3 + t37 * (-t16 * t4 - t60)) +
     # t3 * (-t60 * t27 + t43) * t22 * t12) + t36 * t11 * (t19 * (-t61 -
     # t10) - t3)) * t2) + t20 * (t53 * (t45 * t47 * t19 - t48 * t4) * t
     #3 + t22 * (-t24 * t38 + (-t36 - t35) * t11 * t1) * t3 * t41 - t51 
     #* (t50 + t49) + t52 * (t30 * (-t20 * t11 * t26 * t30 + t27 * t37) 
     #+ t48) * t6 + t22 * (t23 * (t13 * (-t1 * t17 + t42) + t44 * t25) +
     # t35 * (-t1 * t15 + t42) + t47 * (t34 * t9 + t45) * t19) * t41)
      t29 = t20 * t3
      t33 = t11 * t55
      t35 = t13 + t27
      t36 = t22 * t41
      t43 = t3 * t22
      t45 = t56 * t27
      t48 = t58 * t13
      t60 = t48 + t45
      t66 = t63 * t27
      t67 = t62 * t13
      t68 = t26 * t19
      t69 = t10 * t11 + t17 * t5
      t70 = t58 + t56
      t71 = t20 ** 2
      t72 = t5 * t15
      t73 = t11 * t37
      t6 = t71 * t6
      t15 = t6 * (-t53 * t5 * (t45 * t15 + t48 * t17) + t43 * (t11 * t70
     # * t26 * t2 + t45 * (-t73 - t72) - t48 * t69))
      t6 = t6 * (t52 * t21 * t22 * t26 * t70 + t53 * t11 * (t28 * (t59 +
     # t57) * t26 * t20 * t22 + t66 + t67) + t36 * t5 * t11 * t60)
      t17 = t3 * t2 * (t20 * (t43 * t46 * t55 * t14 - t49 - t50) + t44 *
     # (-t14 * t20 * t26 * t55 + (t3 * (-t25 - t12) + t8 * (-t37 - t10))
     # * t19 * t31))
      t18 = t20 * t71 * t22 * t16 * t18 * t5 * (t59 * t13 * t69 + (t73 +
     # t72) * t27 * t57)
      ret = t32 * t6 + 3 * t29 * (-t4 * t3 * t22 * t54 + t52 * t14 * t38
     # + t44 * t34) + (0.3D1 / 0.2D1) * t29 * t26 * t19 * t39 * (t33 + t
     #53 * (-t30 - t23) * t9) - 6 * t7 * t19 * t31 * t41 * t2 * (t43 * (
     #t25 + t12) + t37 + t10) - (0.3D1 / 0.4D1) * t15 - t18 / 8 - 8 * t7
     # * t31 * t3 * (t41 * (t43 * t1 + t2) + t42 * t68 * (t36 * t1 - t51
     #)) + t24 + 4 * t3 * (t31 * (t3 * (t4 * (t2 * (t13 * (t19 * (t61 + 
     #t10) + t3) + t27 * (t19 * (t3 * t12 * t22 + t37) + t3)) + t36 * t3
     #5 * t1) - t19 * t39 * t35 * t14) - t14 * t40 * t19 * t35 * t16 * t
     #8) + t64 * t34 * t20) - t29 * t22 * (t65 * (t5 * (-t4 * (t67 + t66
     #) + t42 * t60) + t47 * (t11 * (t12 * t56 + t25 * t58) + t9 * (t63 
     #+ t62))) * t20 + t68 * t3 * t2 * (t55 * t9 * t2 - t33 * t1)) / 2 +
     # 2 * t17

      hjetmass_box_ppmm_0_s14_0_mhsq_s124_s134_dp = ret/32d0/(0,1d0)
      return

      end function
