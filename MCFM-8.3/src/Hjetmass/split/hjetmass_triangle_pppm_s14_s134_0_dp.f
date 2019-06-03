
      double complex function hjetmass_triangle_pppm_s14_s134_0_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt
          double precision cg

      t1 = za(i1, i4)
      t2 = zb(i4, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i3, i4)
      t6 = zb(i4, i3)
      t7 = t3 * t4
      t8 = t5 * t6
      t9 = t8 + t7
      if ( dreal(t9) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t7 = cg * cdsqrt(t9 ** 2) + t7 + t8
      t9 = zb(i3, i2)
      t10 = zb(i4, i2)
      t11 = 0.1D1 / t7
      t12 = 2 * t3
      t13 = t1 * (-t12 * t9 * t2 * t11 + t10)
      t12 = t12 * t4 * t11
      t14 = t1 * t2
      t15 = t14 * (1 - t12)
      t16 = za(i2, i3)
      t12 = t8 * t12 * t14
      t14 = -t12 + t4 * (t13 * t16 + t15 * t3)
      t17 = zb(i2, i1)
      t18 = za(i2, i4)
      t19 = za(i1, i2)
      t20 = t1 * t16
      t21 = t2 * (-2 * t20 * t4 * t11 + t18)
      t12 = -t12 + t3 * (t15 * t4 + t21 * t9)
      t22 = 0.1D1 / t16
      t23 = 0.1D1 / t3
      t24 = 0.1D1 / t7
      t2 = 0.1D1 / t2
      t25 = 0.1D1 / t13
      t26 = 0.1D1 / t6
      t27 = 0.1D1 / t1
      t14 = 0.1D1 / t14
      t28 = 0.1D1 / t18
      t29 = t5 * t28
      t30 = t13 * t26
      t31 = t30 * t27
      t32 = t31 + t29
      t33 = t4 * t28
      t34 = t33 * t5
      t35 = t34 + t17
      t36 = t14 ** 2
      t37 = t14 * t36
      t38 = t17 ** 2
      t39 = t4 ** 2
      t40 = t4 * t39
      t41 = t13 ** 2
      t42 = t10 * t2
      t43 = t42 * t5
      t44 = t1 * t28
      t45 = t4 * t13
      t46 = t9 * t15
      t47 = t31 * t4
      t48 = t5 * t4
      t49 = t27 * t2
      t50 = t15 * t25
      t51 = t15 * t27
      t52 = t26 * t27 ** 2
      t53 = t29 * t27
      t54 = t14 * t4
      t55 = t2 * t23
      t56 = t55 * t5
      t57 = t55 * t15
      t58 = t10 * t18
      t59 = t17 * t19
      t60 = t59 + t58
      t61 = t59 * t28 + t10
      t62 = t5 ** 2
      t63 = t15 ** 2
      t64 = t25 ** 2
      t65 = t2 ** 2
      t66 = t9 ** 2
      t67 = t17 * t26
      t68 = t62 * t7
      t69 = t26 * t27
      t70 = t22 * t23
      t71 = t17 * t14
      t72 = t69 * t22
      t73 = t65 * t23
      t31 = t73 * t15 * t7 * (t66 * (t53 * t22 * t64 * t26 * t63 + t50 *
     # t52 * t5 * t22) + t9 * (t15 * (t22 * (t29 * t17 * t64 - t69 * t35
     # * t25) + t54 * t27 * (t68 * t23 * t24 * t28 + t67 + t54 * t5 * (t
     #31 * t60 + t5 * t61) + t31 * t5 * t7 * t23 * t24)) - t70 * t53 * t
     #7 * t24 * t25 * t26 * t63) + t45 * t27 * (t5 * (t72 + t71) * t24 *
     # t23 * t7 + t54 * (-t59 * t33 * t62 * t14 - t67 + (t17 * t60 + t47
     # * (-t59 - t58)) * t14 * t5)))
      t58 = t59 * t2
      t60 = t28 * (t58 + t1) + t42
      t67 = t9 * t16
      t74 = t67 * t23 + t4
      t75 = t45 * t14
      t76 = -t75 + t22
      t77 = t23 ** 2
      t78 = mt ** 2
      t79 = t7 ** 2
      t80 = t24 ** 2
      t81 = t46 * t23
      t82 = t15 * t60 * t36 * t13
      t83 = t1 * t17
      t84 = t9 * t63
      t85 = t45 * t23
      t58 = t48 * (t7 * (t54 * t15 * (t4 * (t82 * t23 * t6 * t62 + (t81 
     #* t16 * t60 * t36 + t23 * (t28 * (-t58 - t1) - t42) * t14) * t13 *
     # t5) + t82 * t5 * t39 - t83 * t55 * t28) * t24 + t34 * t1 * (t14 *
     # (t14 * (t84 * t74 - t75 * t63 * (t66 * t16 ** 2 * t23 + t3 * t39)
     #) - t85) + t70 - t85 * t62 * t6 ** 2 * t63 * t37 + t84 * t8 * t23 
     #* t36) * t80) + t29 * t15 * t77 * t2 * t76 * t80 * t79 + t73 * t51
     # * t33 * t78 * t26 * t76)
      t76 = t13 * t23
      t82 = t76 * t5
      t86 = t9 * t28
      t87 = t46 * t14
      t88 = t14 * t23
      t89 = t2 * t24
      t90 = t83 * t6
      t16 = t56 * t24 * t7 * t4 * (-t17 * t22 + t45 * (-t69 * t66 * t63 
     #* t16 * t36 - t46 * t16 * t17 * t36 + t72) + t69 * t46 * t16 * t39
     # * t41 * t36 - t33 * t6 * t15 * t36 * (t90 + t46) * t62 + t34 * t1
     #4 * (t87 * t16 * (t45 - t46) + t90)) + t89 * t7 * t5 * (t15 * t13 
     #* t36 * t32 * t39 ** 2 + t14 * (-t84 * t32 * t14 + (t13 * (t82 * t
     #27 - t17) + t8 * (t82 - t83) * t28) * t14 * t15 - t76 * t32) * t40
     # + t84 * t70 * t28 * (t46 * t26 + t83) * t64 + t88 * (t13 * t17 + 
     #t5 * (t15 * (t86 + t71 * (-t86 * t20 - t13) * t6) - t84 * t13 * t2
     #7 * t14)) * t39 + t81 * (t69 * (t87 - t22) + t71) * t4)
      t62 = t13 * t39
      t64 = t22 * t25
      t20 = t63 * t5 * (t54 * t7 * (t39 * (-t89 * t28 * (t82 * t7 * t24 
     #+ t83) * t14 + t56 * t13 * (t59 * (t42 * t27 + t28) + t10) * t36) 
     #+ t29 * t2 * t80 * t7 * t9 * t77 + t54 * t24 * t23 * (t28 * (t1 * 
     #t5 * t9 + (t9 * (t19 * t5 - t20) - t8 * t1) * t2 * t17) - t55 * t2
     #9 * (t67 + t8) * t24 * t13 * t7 + t43 * t9)) + t86 * t73 * t69 * (
     #t4 * (-t64 + t54) + (-t25 * t22 ** 2 + t62 * t36) * t27 * (t18 * t
     #6 + t19 * t4) * t5) * t78)
      t21 = 0.1D1 / t21
      t12 = 0.1D1 / t12
      t3 = -t3 * t9 * t12 + t21
      t12 = t75 - t22
      t21 = t2 * t17
      t10 = t5 * (t7 * (t63 * (t52 * t33 * t78 * t2 * t65 * t12 * t77 + 
     #t36 * t39 * (t21 * t60 + t75 * ((t38 * t19 ** 2 * t28 + t18 * t10 
     #** 2) * t27 * t65 + t44) * t5) * t23) + t57 * t33 * t22 * t24 * (-
     #t83 * t25 + t26 * t4)) + t89 * t54 * (t21 * t28 + t75 * (t49 * t61
     # + t28) * t5) * t77 * t63 * t79 + t26 * t11 * t40 * (t33 * t3 * t1
     # + t3 * t9) + t53 * t80 * t65 * t12 * t23 * t77 * t63 * t7 * t79)
      ret = -32 * t58 + 4 * t16 + 16 * t20 + 8 * t10 - 64 * t29 * t39 * 
     #t1 * (t78 * t39 * t11 ** 2 * t26 * t3 + t68 * (-t45 * t6 * t74 * t
     #37 * t80 * t63 + t85 * t6 * t80 * t36 * t15) + (-t62 * t67 * t80 *
     # t37 * t63 + (t45 * t74 * t36 - t88 * t9) * t80 * t15) * t7 * t5) 
     #+ 2 * t31 + 2 * t57 * t7 * (t4 * (t49 * t5 * (-t30 * t39 * t28 + t
     #17 * t9 - t47 * t9) * t14 + t48 * (t17 * (t8 * (t42 + t44) + t13) 
     #+ t45 * (t27 * (-t43 - t30) - t29) + t46 * t32 + t8 * t38 * t19 * 
     #t2 * t28) * t36 + t2 * t26 * t27 * t35 * t22) + t56 * (t22 * (t17 
     #* (-t50 * t28 - t27) + t51 * (-t27 * t9 + t33) * t26) + t54 * (t4 
     #* (-t53 * t13 - t52 * t41) + t8 * t17 * t28)) * t24 * t7) - 12 * t
     #64 * t84 * t55 * t34 * t24 * t7 * t26

      hjetmass_triangle_pppm_s14_s134_0_dp = ret/32d0/(0,1d0)
      return

      end function
