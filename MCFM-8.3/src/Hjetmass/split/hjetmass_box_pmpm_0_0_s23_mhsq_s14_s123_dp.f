
      double complex function hjetmass_box_pmpm_0_0_s23_mhsq_s14_s123_dp 
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
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i1, i4)
      t8 = za(i2, i4)
      t9 = zb(i4, i1)
      t10 = za(i3, i4)
      t11 = t8 * t2
      t12 = t10 * t4
      t13 = t12 + t11
      t14 = t7 * t9
      t15 = t14 * t13
      t16 = t1 * t2
      t17 = t3 * t4
      t18 = t5 * t6
      t19 = t14 * (t18 + t16 + t17)
      t20 = zb(i4, i2)
      t21 = zb(i4, i3)
      t22 = t1 * t20
      t23 = t3 * t21
      t24 = t23 + t22
      t25 = mt ** 2
      t26 = 16 * t15 * t25 * t24 + 4 * t19 ** 2
      t26 = cdsqrt(t26)
      t15 = 0.1D1 / t15
      t27 = t17 + t16
      t28 = t14 * t15
      t11 = t28 * (t11 * t27 + t18 * (t12 + t11) + t12 * t27)
      t27 = 2
      t29 = (0.1D1 / 0.2D1)
      t13 = t29 * t26 * t15 * t13
      t16 = t27 * (t18 + t16 + t17)
      t18 = t16 - t11 - t13
      t11 = t16 - t11 + t13
      t13 = t27 * t19
      t16 = -t13 - t26
      t19 = 0.1D1 / t24
      t24 = t8 * t16 * t15 / 4
      t30 = t22 * t29
      t31 = t2 * (t24 + t1) - t30 * t18 * t19
      t13 = -t13 + t26
      t26 = t8 * t13 * t15 / 4
      t30 = t2 * (t26 + t1) - t30 * t11 * t19
      t32 = t29 * t18 * t19
      t24 = t32 * t1 * t21 - t24 * t4
      t33 = t29 * t11 * t19
      t26 = t33 * t1 * t21 - t26 * t4
      t12 = t12 / 4
      t32 = -t12 * t16 * t15 + t23 * t32
      t12 = -t12 * t13 * t15 + t23 * t33
      t33 = t4 * t20
      t34 = -t21 * t2 + t33
      t35 = t11 * t13
      t36 = t18 * t16
      t2 = 0.1D1 / t2
      t30 = 0.1D1 / t30
      t10 = 0.1D1 / t10
      t37 = 0.1D1 / t6
      t38 = 0.1D1 / t21
      t5 = 0.1D1 / t5
      t31 = 0.1D1 / t31
      t39 = t11 * t26
      t40 = t18 * t24
      t41 = t40 + t39
      t42 = t30 * t11 + t18 * t31
      t43 = t11 ** 2
      t44 = t11 * t43
      t45 = t18 ** 2
      t46 = t18 * t45
      t47 = t43 * t30
      t48 = t45 * t31
      t49 = t48 + t47
      t50 = -t18 - t11
      t51 = t18 + t11
      t52 = t4 ** 2
      t53 = t9 ** 2
      t54 = t8 ** 2
      t55 = t8 * t54
      t56 = t25 * t4
      t57 = t42 * t6
      t58 = t5 * t41
      t59 = t11 * t12
      t60 = t37 * t5
      t61 = t60 * t52
      t62 = t25 * t8
      t63 = t45 * t24
      t64 = t43 * t26
      t65 = t52 * t2
      t66 = t7 * t8
      t67 = t10 * t8
      t68 = t46 * t31
      t69 = t44 * t30
      t70 = t69 + t68
      t44 = t46 + t44
      t46 = t7 ** 2
      t71 = t19 ** 2
      t72 = t1 * t4
      t73 = t45 + t43
      t74 = t1 * t2
      t75 = t1 * t73
      t76 = t67 * t19
      t77 = t60 * t4
      t78 = t2 * t4
      t79 = t5 * t38
      t80 = t79 * t10
      t81 = t19 * t8
      t82 = t81 * (t1 * (-t14 * t21 * t10 * t70 * t71 + t76 * t6 * t9 * 
     #t49) + t4 * (t62 * t60 * t10 * t51 + (t50 * t10 * t7 * t4 + t58) *
     # t38 * t2 * t9))
      t12 = t82 + t19 * (t14 * (t67 * t2 * (t38 * (t48 * t22 * t52 + (t4
     #8 * (-t72 + t24) + t47 * (-t72 + t26)) * t9 * t6) + t77 * (t12 * t
     #43 + t32 * t45 + t17 * (-t45 - t43))) + t78 * (-t64 * t30 - t63 * 
     #t31 + t72 * t49) + t21 * t49 * t10 * t54) * t19 + t80 * t2 * t54 *
     # t4 * (t56 * t51 * t37 * t20 - t59 * t9) - t2 * t53 * t10 * (t68 *
     # t24 + t69 * t26) * t71 * t46) + t76 * (t7 * (t4 * (t1 * (t60 * t7
     #3 + t47 + t48) + (t43 * (t30 * (t72 - t26) - t60 * t26) + t63 * (-
     #t60 - t31)) * t38 * t2 * t20) * t19 * t9 + t52 * t37 * (t33 * t51 
     #* t38 * t2 - t11 - t18) + t74 * t5 * t44 * t71 * t53) + t5 * t9 * 
     #(t75 * t19 + (-t11 * t25 + t18 * (-t25 - t32)) * t38 * t2 * t4) * 
     #t8)
      t14 = t36 + t35
      t17 = t7 * t21 * t10
      t13 = t43 * t13
      t16 = t45 * t16
      t24 = t19 * t7
      t1 = t24 * (t2 * (t72 * t9 * (-t68 * t21 + t69 * (-t67 * t20 - t21
     #)) * t71 + t17 * t1 * t53 * (t30 * t43 ** 2 + t31 * t45 ** 2) * t1
     #9 * t71 + t76 * t4 * (t43 * (-t33 * t37 - t9) - t45 * t9)) - t54 *
     # t52 * t9 * t15 * t10 * t38 * (t35 * t30 + t36 * t31)) + t81 * (t1
     #9 * (t4 * (t17 * t37 * t73 + (t5 * (t67 * t73 * t3 - t75) + t10 * 
     #(-t68 * t22 + t60 * (t22 * t44 - t23 * t44)) * t19 * t7) * t2 * t9
     #) - t65 * t45 * t7 * t20 * t10 * t37 - t74 * t7 * t53 * t10 * t19 
     #* t70 * t6) + t79 * t28 * t4 * (t65 * t37 * t14 + t67 * (-t4 * t37
     # * t14 + (-t16 - t13) * t19 * t2 * t9)))
      t3 = t2 * t38
      t14 = 0.1D1 / t7
      t17 = t10 * t4 * (t25 * (t54 * (t19 * t42 + (-t30 - t31) * t38 * t
     #4) - t6 * t14 * t38 * (t30 + t31) * t55) + t78 * t71 * t46 * t9 * 
     #t37 * t73)
      ret = t27 * (t19 * (t62 * t61 * t2 * t51 + (t4 * (t21 * t2 * t71 *
     # t44 * t37 * t46 - t66 * t60 * t38 * t41) + t6 * t54 * t38 * (t39 
     #* t30 + t40 * t31)) * t10 * t9 - t72 * t46 * t2 * t71 * t70 * t10 
     #* t53) + t19 * (t67 * (t54 * (t5 * t51 + t57) + t8 * (-t57 * t56 *
     # t38 * t2 + t42 * t7 * t4 + t58 * t38) + t61 * t7 * t2 * t38 * (t1
     #8 * t32 + t59) + t56 * t7 * t2 * t19 * t49) * t9 + t65 * (t50 * t3
     #7 * t7 * t4 + t62 * t42) + t65 * t25 * t54 * t10 * t38 * t42 * t20
     # - t66 * t2 * t10 * t38 * t5 * t19 * (t64 + t63) * t53)) + t29 * t
     #1 + t77 * t19 * t10 * t9 * t54 * (t24 * (t3 * t33 * (t16 + t13) - 
     #t13 - t16) * t15 + t3 * (t11 * t35 * t15 * t19 * t7 * t34 + t18 * 
     #t36 * t15 * t19 * t7 * t34)) / 8 - t12 - 4 * t17 + 8 * t80 * t56 *
     # t55 * t14

      hjetmass_box_pmpm_0_0_s23_mhsq_s14_s123_dp = ret/32d0/(0,1d0)
      return

      end function
