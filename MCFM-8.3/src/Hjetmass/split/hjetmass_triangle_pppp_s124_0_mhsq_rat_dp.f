
      double complex function hjetmass_triangle_pppp_s124_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
      implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg


      t1 = zb(i3, i1)
      t2 = za(i1, i4)
      t3 = za(i2, i3)
      t4 = za(i2, i4)
      t5 = za(i1, i2)
      t6 = za(i3, i4)
      t7 = zb(i3, i2)
      t8 = t2 * t1
      t9 = t4 * t7 + t8
      t10 = zb(i4, i1)
      t11 = za(i1, i3)
      t12 = zb(i4, i3)
      t13 = t11 * t1
      t14 = t3 * t7
      t15 = t6 * t12
      t16 = t13 + t14 + t15
      if ( dreal(t16) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t16 = cg * cdsqrt(t16 ** 2) + t13 + t14 + t15
      t17 = zb(i2, i1)
      t18 = zb(i4, i2)
      t19 = t5 * t17
      t20 = t2 * t10
      t21 = t4 * t18
      t22 = t21 + t19 + t20
      t16 = 0.1D1 / t16
      t23 = -2 * t13 * t22 * t16 + t19 + t20
      t24 = 2 * t3
      t25 = -t24 * t1 * t22 * t16 + t10 * t4
      t26 = 2 * t6
      t27 = t26 * t1 * t22 * t16 - t17 * t4
      t28 = t5 * t1
      t29 = t12 * t4 + t28
      t30 = t1 * t23
      t31 = t27 * t12
      t32 = t11 * (-t31 + t30)
      t33 = 2 * t11
      t34 = -t33 * t7 * t22 * t16 + t18 * t2
      t35 = -t12 * t2 + t5 * t7
      t33 = t33 * t22 * t12 * t16 - t18 * t5
      t36 = t11 * t23
      t37 = t6 * t33
      t38 = t1 * (t36 - t37)
      t21 = -t26 * t22 * t12 * t16 + t20 + t21
      t26 = t11 * (t25 * t7 + t30)
      t30 = t3 * t34
      t36 = t1 * (t30 + t36)
      t24 = t24 * t22 * t12 * t16 + t10 * t5
      t39 = 0.1D1 / t34
      t40 = 0.1D1 / t3
      t41 = 0.1D1 / t33
      t42 = 0.1D1 / t27
      t43 = 0.1D1 / t6
      t38 = 0.1D1 / t38
      t44 = 0.1D1 / t25
      t22 = 0.1D1 / t22
      t45 = 0.1D1 / t4
      t26 = 0.1D1 / t26
      t36 = 0.1D1 / t36
      t2 = 0.1D1 / t2
      t5 = 0.1D1 / t5
      t46 = t12 * t23
      t47 = t10 * t35 - t46
      t48 = t7 * t23
      t49 = t2 ** 2
      t50 = t5 ** 2
      t51 = (t17 * t35 + t48) * t38
      t52 = t10 * t40
      t53 = t17 * t43
      t54 = t11 * t25 * t40
      t55 = t9 * t23
      t56 = t16 * t50
      t57 = t3 * t35
      t58 = t11 * t29
      t59 = t9 * t11
      t60 = t6 * t35
      t61 = t48 * t39
      t62 = t61 * t40
      t63 = t46 * t16
      t64 = t12 * t5
      t65 = t7 * t43
      t66 = t12 * t34
      t67 = t2 * t5
      t32 = 0.1D1 / t32
      t68 = t25 ** 2
      t69 = t9 * t16 * t44 * t49
      t70 = t40 * t2
      t71 = t7 * t11
      t72 = t27 * t43
      t73 = t72 * t26
      t74 = t25 * t40 * t32
      t75 = t1 * t45
      t76 = t71 * t64
      t77 = t3 ** 2
      t78 = t11 ** 2
      t79 = t43 ** 2
      t25 = t25 * t42
      t80 = t7 * t29
      t81 = t3 * t43
      t82 = t66 * t77
      t83 = t82 * t9
      t84 = t9 * t12
      t85 = t84 * t73
      t86 = t43 * t5
      t87 = t10 * t43
      t88 = t70 * t17
      t89 = t75 + t64
      t90 = t40 ** 2
      t91 = t6 ** 2
      t92 = t1 ** 2
      t93 = t1 * t92
      t94 = t92 * t11
      t95 = t7 * t39
      t96 = t60 * t21
      t97 = (t11 * t21 - t37) * t36 * t35
      t98 = t23 * t45
      t14 = t49 * (-t95 * t91 * t35 * t24 * t90 * t45 + t5 * t7 * (t97 *
     # t1 - t85 * t78) + (-t92 * t91 * t35 * t33 * t36 + t96 * (t94 * t3
     #6 - t95) + t13 * t9) * t45 * t40) + t67 * (t1 * (-t57 * t89 + t58 
     #* t89) * t38 * t34 + t76 * t29 * t32 * (t54 + t23) + t97 * t45 * t
     #92 + t98 * t13 * (t84 * t26 + t80 * t32)) + t53 * t9 * t45 * t2 * 
     #(t14 + t13) * t22 + t81 * t12 * t45 * t50 * (-t55 * t11 * t43 + t2
     #9 * t34) * t41 + t81 * t92 * t34 * t50 * t45 * t38 * (t58 - t57)
      t4 = t16 * t14 + t45 * (t5 * (((t13 + t15) * t16 * t40 * t29 - t8 
     #* t43) * t22 * t10 + t1 * t78 * t16 * t2 * (t80 * t74 - t85)) + t7
     #0 * (t6 * t16 * t2 * (t58 * t62 - t84) - t28 * t17 * t22) + t56 * 
     #(t1 * (-t80 * t68 * t40 * t42 * t32 * t78 + t58 * (-t25 * t48 * t3
     #2 + t43)) + t80 * (t25 - t81) + t83 * t41 * t79) + t31 * t69 * (1 
     #+ t13 * t26 * (t72 * t11 - t23))) + (t59 * t63 * t26 * t49 - t87 *
     # t22) * t5 * t7 + t22 * (t1 * (-t86 - t70) - t67 * (t12 * t40 + t6
     #5) * t4) * t18 - t88 * t12 * t22
      t8 = t17 * t45 + t18 * t5
      t14 = t38 ** 2
      t15 = t36 ** 2
      t28 = t60 * t7 * t24 * t15 * t2
      t31 = t9 * t10
      t72 = t17 * t29
      t78 = t70 * t6
      t85 = t10 * t45 + t18 * t2
      t89 = t35 * t24
      t92 = t98 * t29
      t19 = t22 * (t12 * t29 * t50 * (t20 * t45 + t18) + (t49 * (t19 * t
     #45 + t18) + t81 * t85 * t5) * t9 * t7) + t13 * (t5 * (t2 * (t84 * 
     #(t2 * t34 + t98) * t36 + t7 * t38 * (t89 * t5 + t92)) + t9 * t85 *
     # t22 * t43 - t82 * t55 * t41 * t79 * t38 * t45 * t5 - t81 * t92 * 
     #t66 * t41 * t38 * t5) + t89 * t61 * t91 * t90 * t45 * t36 * t49 + 
     #t96 * t62 * t45 * t36 * t49) + t98 * t94 * t5 * (t86 * t83 * t14 +
     # t28) + t67 * t34 * t35 * t33 * (-t3 * t14 * t5 + t6 * t2 * t15) *
     # t11 * t93
      t8 = t16 * t19 + t16 * (t22 * (t45 * ((-t72 + t31) * t2 * t7 + t64
     # * (t72 - t31)) - t67 * (t80 + t84) * t18 + t78 * t12 * t29 * t8) 
     #+ t13 * (t40 * (t46 * t6 * t9 * t45 * t36 * t49 + t29 * t22 * t8 *
     # t2) + t81 * t48 * t29 * t38 * t45 * t50) + t98 * t35 * (t77 * t34
     # * t43 * t14 * t50 + t30 * t67 * t14 + t37 * t15 * t2 * (t78 + t5)
     #) * t11 * t93 + t94 * t2 * (t34 * (t5 * (t46 * t3 * t9 * t14 * t45
     # + t28) + t57 * t12 * t21 * t14 * t50) + t70 * t48 * t91 * t35 * t
     #24 * t45 * t15))
      t9 = -8
      ret = t9 * (t22 * (-t65 * t17 * t2 - t64 * t52) + t45 * (t1 * (t2 
     #* (-t11 * t12 * t27 ** 2 * t43 * t44 * t26 - t51) + t22 * (-t52 - 
     #t53) + t47 * t5 * t36 + t55 * t16 * t44 * t49 - t56 * (t54 + t23) 
     #* t42 * t29) + t62 * t16 * t49 * (t59 + t60) + t63 * t43 * t50 * (
     #-t58 + t57) * t41) + t67 * (t33 * t7 * t38 - t66 * t36) * t1 + t67
     # * t1 * t35 * (-t38 + t36) * t18 + t75 * (t43 * (-t51 * t5 * t3 + 
     #t11 * t27 * (t64 * t26 - t69)) - t54 * t7 * t32 * t2 + t70 * t47 *
     # t36 * t6 + t71 * t68 * t40 * t42 * t32 * t5) + t76 * t2 * (t73 - 
     #t74)) - 16 * t4 - 32 * t8 + 48 * t16 * t45 * t22 * (t80 * t10 * t5
     # + t84 * t17 * t2) - 4 * t45 * (t1 * (t40 * (t25 * t5 - t2) + t43 
     #* (t27 * t44 * t2 - t5)) + t35 * (t87 * t41 * t5 + t88 * t39))
          



      hjetmass_triangle_pppp_s124_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
