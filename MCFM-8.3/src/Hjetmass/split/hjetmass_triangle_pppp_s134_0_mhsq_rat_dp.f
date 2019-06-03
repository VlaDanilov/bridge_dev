
      double complex function hjetmass_triangle_pppp_s134_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = zb(i2, i1)
      t2 = zb(i4, i1)
      t3 = za(i2, i3)
      t4 = za(i3, i4)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = za(i1, i4)
      t8 = zb(i4, i3)
      t9 = t5 * t6
      t10 = t7 * t2
      t11 = t4 * t8
      t12 = t11 + t9 + t10
      t13 = zb(i3, i2)
      t14 = za(i1, i2)
      t15 = zb(i4, i2)
      t16 = t13 * t5 + t15 * t7
      t17 = za(i2, i4)
      t18 = t14 * t1
      t19 = t3 * t13
      t20 = t17 * t15
      t21 = t20 + t18 + t19
      if ( dreal(t21) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t21 = cg * cdsqrt(t21 ** 2) + t18 + t19 + t20
      t21 = 0.1D1 / t21
      t22 = -2 * t18 * t12 * t21 + t10 + t9
      t23 = 2 * t3
      t24 = t23 * t1 * t12 * t21 + t2 * t4
      t25 = t1 * t5 + t15 * t4
      t26 = 2 * t14
      t27 = t26 * t13 * t12 * t21 + t7 * t8
      t9 = -t23 * t13 * t12 * t21 + t11 + t9
      t26 = t26 * t15 * t12 * t21 - t5 * t8
      t23 = -t23 * t15 * t12 * t21 + t2 * t5
      t28 = -2 * t20 * t12 * t21 + t10 + t11
      t29 = t1 * t7 - t13 * t4
      t30 = 2 * t17
      t31 = t15 * (t30 * t1 * t12 * t21 - t4 * t6)
      t32 = t1 * t22
      t33 = t14 * (t32 - t31)
      t31 = t14 * (t13 * t24 + t31)
      t34 = t17 * t26
      t35 = t14 * t22
      t36 = t1 * (t35 - t34)
      t37 = t3 * t27
      t38 = t1 * (t37 + t34)
      t30 = -t30 * t13 * t12 * t21 + t6 * t7
      t38 = 0.1D1 / t38
      t39 = 0.1D1 / t5
      t12 = 0.1D1 / t12
      t40 = 0.1D1 / t22
      t41 = 0.1D1 / t24
      t42 = 0.1D1 / t4
      t33 = 0.1D1 / t33
      t43 = 0.1D1 / t14
      t36 = 0.1D1 / t36
      t44 = 0.1D1 / t27
      t45 = 0.1D1 / t3
      t31 = 0.1D1 / t31
      t7 = 0.1D1 / t7
      t46 = t17 * t25
      t47 = t46 * t45 + t29
      t48 = t42 ** 2
      t49 = t22 ** 2
      t50 = t24 ** 2
      t51 = t7 ** 2
      t52 = t1 ** 2
      t53 = t8 * t39
      t54 = t2 * t12
      t55 = t15 * t39
      t56 = t55 * t40 * t33
      t57 = t1 * t45
      t58 = t52 * t26
      t59 = t15 * t24
      t60 = t3 * t28
      t61 = t13 * t7
      t62 = t16 * t17
      t63 = t62 * t43
      t64 = t45 * t39
      t65 = t55 * t21
      t66 = t42 * t43
      t67 = t66 * t3 + t7
      t68 = t45 * t7
      t69 = t8 * t16
      t70 = t61 * t22
      t71 = t6 * t42
      t72 = t38 + t31
      t73 = t64 * t20 * t14 * t25
      t74 = t43 * t39
      t75 = t14 * t29
      t76 = t12 * t8
      t77 = t6 * t12
      t78 = t14 ** 2
      t79 = t16 * t45
      t80 = (t17 * t9 + t3 * t30) * t36
      t81 = t25 * t48 * t43
      t79 = t39 * (t51 * (t15 * (-t78 * t29 * t38 * t45 * t22 + t41 * t3
     #1 * (t75 - t62) * t49) + t79 * t14) + t42 * t15 * (t16 * (t72 * t1
     #7 * t22 + t80) - t35 * t72 * t29) * t7 + t62 * t3 * t26 * t40 * t4
     #8 * t43 ** 2 + t3 ** 2 * t15 * t16 * t30 * t48 * t43 * t36) * t1 +
     # t39 * (-t62 * t14 * t23 * t44 * t45 ** 2 * t51 + t81 * t3) * t13 
     #- t20 * t39 * (t79 * t51 + t81) + t80 * t16 * t48 * t7 * t52
      t80 = t71 * t7 * t12
      t5 = t21 * t79 + t12 * (t66 * t13 * t2 + t57 * t8 * t7) + t21 * (t
     #42 * (t13 * (t73 * t24 * t33 * t7 + t32 * (t62 * t72 - t75 * t72) 
     #* t51) - t54 * t64 * t25 * (t20 + t18)) + t48 * (t13 * (t24 * (t18
     # * t47 * t7 * t33 - t25 * t39 * t40) - t73 * t40 * t33 * t50) + t7
     #4 * t16 * t3 * t1 * (t40 * (-t37 * t43 - t9) + t20 * t9 * t36)) - 
     #t35 * t64 * t13 * t51 * (t14 * t25 * t45 + t16) * t44 + t7 * t39 *
     # t16 * (-t76 * t43 * (t20 + t19) + t32 * (t20 * t14 * t38 * t45 - 
     #t41) * t7)) + t77 * t15 * (t68 + t66) + t80 * (t13 * t43 + t57) * 
     #t5 + t55 * (t10 * t66 + t11 * t68) * t12 - t75 * t56 * t13 * t50 *
     # t21 * t48 + t75 * t61 * t65 * t24 * t42 * t33
      t11 = t23 * t13
      t32 = t15 * t28 + t11
      t19 = t20 + t19
      t62 = t20 + t18
      t72 = t38 ** 2
      t73 = t36 ** 2
      t75 = t34 * t42
      t79 = t42 * t1
      t4 = t16 * (t51 * ((t38 * (t64 * t14 * t26 + t27 * t42) - t35 * t1
     #3 * t9 * t42 * t72) * t17 * t52 + t13 * t12 * (t53 * t4 + t6)) + t
     #7 * ((-t37 * t17 * t22 * t48 * t73 + t75 * t37 * t39 * t73) * t1 *
     # t52 + (-t11 * t17 * t39 * t36 - t77) * t42 * t1 + t54 * t13 * t39
     #) + t74 * t1 * t3 * t17 * t48 * t36 * (t37 * t58 * t36 - t11) - t7
     #9 * t54 * t39) + t12 * (t79 * (t42 * (t10 * t39 + t6) + t53) + t61
     # * (-t71 - t53)) * t25 + t20 * t52 * t16 * t22 * t42 * t7 * (-t7 *
     # t14 * t72 + t3 * t42 * t73) * t30
      t4 = t21 * t4 + t21 * (t39 * (t16 * (t52 * (t72 * (-t78 * t17 * t2
     #2 * t45 * t32 * t51 - t35 * t17 * t42 * t32 * t7) + t75 * t7 * t38
     #) - t18 * t13 * t17 * t22 * t44 * t45 * (t14 * t23 * t45 + t26) * 
     #t51 * t38 + t66 * t54 * t19) + t27 * (-t60 * t20 * t52 * t16 * t42
     # * t67 * t73 + t63 * t52 * t3 * t40 * t48 * (t3 * t26 * t43 + t23)
     # * t36) + t76 * t68 * t62 * t25) + t70 * t46 * t1 * t48 * t36 + t8
     #0 * (t16 * t43 * t19 + t62 * t45 * t25))
      t6 = t39 * t21 * t12 * (t13 * t2 * t25 * t42 + t69 * t1 * t7)
      t9 = t39 * (t15 * (t43 * (-t41 * t22 * t7 + t42) + t45 * (-t24 * t
     #40 * t42 + t7)) + t16 * (t66 * t2 * t40 + t68 * t8 * t44))
      t10 = 8
      ret = t10 * (t1 * (t39 * (-t36 * t67 * t26 * t13 + t16 * (t68 * t1
     #4 + t42) * t38 * t2 - t15 * t22 * t42 * t31 + t69 * t36 * t7) + t7
     #0 * t42 * (t36 - t31) + t71 * t16 * t7 * (t38 + t36)) + t42 * (t1 
     #* (t54 * t45 + t53 * t3 * t16 * t43 * t36 + (t26 * t39 + t27 * t7)
     # * t38 * t1) + (-t57 * t33 * t7 * t24 + t56 * t45 * t50) * t14 * t
     #13) + t13 * t16 * t39 * t21 * t45 * t51 * (t14 * t28 - t34) * t44 
     #+ t64 * ((-t59 * t13 * t33 + t58 * t38) * t7 * t14 + t2 * t15 * t1
     #2) + t12 * (t61 + t55) * t43 * t8 + t65 * t51 * (t63 - t29) * t41 
     #* t22 + t21 * t40 * t39 * (t59 * t47 + (t17 * t23 + t60) * t43 * t
     #16 * t1) * t48 + t55 * t1 * t49 * t7 * t41 * t31) + 16 * t5 - 32 *
     # t4 + 48 * t6 + 4 * t9

      hjetmass_triangle_pppp_s134_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
