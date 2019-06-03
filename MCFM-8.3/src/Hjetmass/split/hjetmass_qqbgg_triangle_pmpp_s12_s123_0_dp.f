
      double complex function hjetmass_qqbgg_triangle_pmpp_s12_s123_0_dp
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

      t1 = za(i1, i2)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = t2 * t3
      t7 = t4 * t5
      t8 = t7 + t6
      if ( dreal(t8) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t8 = cg * cdsqrt(t8 ** 2) + t6 + t7
      t9 = za(i2, i4)
      t10 = za(i3, i4)
      t11 = zb(i2, i1)
      t12 = 0.1D1 / t8
      t13 = 2 * t1
      t14 = t13 * t10
      t15 = t11 * (t14 * t3 * t12 - t9)
      t16 = t1 * t11
      t17 = t16 * (-2 * t6 * t12 + 1)
      t18 = zb(i4, i1)
      t19 = zb(i4, i3)
      t20 = zb(i4, i2)
      t21 = t1 * (2 * t2 * t11 * t19 * t12 - t20)
      t22 = t3 * t17
      t13 = t7 * t6 * t13 * t11 * t12
      t23 = -t13 + t2 * (-t15 * t19 + t22)
      t24 = t10 * t21
      t13 = -t13 - t3 * (-t17 * t2 + t24)
      t25 = za(i1, i4)
      t26 = t25 * t18
      t27 = t9 * t20
      t14 = -t14 * t11 * t19 * t12 + t26 + t27
      t28 = 0.1D1 / t11
      t29 = 0.1D1 / t21
      t30 = 0.1D1 / t10
      t13 = 0.1D1 / t13
      t31 = 0.1D1 / t1
      t32 = 0.1D1 / t9
      t33 = t13 ** 2
      t34 = t13 * t33
      t35 = t29 ** 2
      t36 = t3 ** 2
      t37 = t3 * t36
      t38 = t5 ** 2
      t39 = t5 * t38
      t40 = t4 ** 2
      t41 = t18 ** 2
      t42 = t1 ** 2
      t43 = t7 * t13
      t44 = t42 * t40
      t45 = t17 * t13
      t46 = t28 * t31
      t47 = t28 * t32
      t48 = t4 * t19
      t49 = t13 * t21
      t50 = t28 * t18
      t51 = t50 * t25
      t52 = t20 * t28
      t53 = t19 * t13
      t54 = t17 * t19
      t55 = t1 * t18
      t56 = t55 * t5
      t57 = t47 * t54
      t58 = t49 * t4
      t59 = t27 + t26
      t60 = t21 ** 2
      t61 = t7 * t1
      t62 = t4 * t3
      t15 = 0.1D1 / t15
      t8 = 0.1D1 / t8
      t23 = 0.1D1 / t23
      t63 = t10 * t19
      t64 = -t63 - t6 - t7
      t65 = t63 + t6
      t66 = -t62 * t32 + t18
      t67 = t19 ** 2
      t68 = t54 * t35
      t69 = t3 * t29
      t70 = t3 * t21
      t71 = t1 * t32
      t72 = t53 * t3
      t73 = t8 * t17
      t74 = t1 * t8
      t24 = t4 * (t72 * t8 * (t22 * t58 + t55) * t5 + t37 * t19 * (t71 *
     # t12 * t23 + t73 * t21 * t33) * t2 + t73 * t67 * (t71 * t17 * t30 
     #* t35 + t24 * t36 * t33) + t74 * t36 * t4 * t21 * t33 * t66 * t38 
     #- t44 * t8 * t36 * t39 * t18 * t32 * t33) + t4 * (t8 * (t72 * (t70
     # + t54) + t71 * (t56 * (t69 + t68) + t36) * t30 + t36 * ((t63 + t6
     # + t7) * t3 * t60 + t5 * t1 * (t6 * t66 + t63 * t66) * t21 + t61 *
     # t54 * t32 * t64 - t42 * t4 * t38 * t18 * t32 * t65) * t33) + t36 
     #* (t15 * (t71 * t3 + t19) + t2 * t67 * t23) * t12)
      t55 = t32 * (t16 + t26) + t20
      t66 = t8 ** 2
      t72 = mt ** 2
      t73 = t22 * t16 * t21 * t32
      t75 = t16 * t32
      t22 = t22 * t21
      t2 = t36 * t4 * t1 * (-t73 * t66 * t4 * t40 * t39 * t34 + (t22 * t
     #55 * t34 * t8 + t75 * (t70 - t54) * t33 * t66) * t38 * t40 + t7 * 
     #(t66 * (t75 * (-t54 * t65 + t70 * t65) * t33 - t73 * (t67 * t10 **
     # 2 + t36 * t2 ** 2) * t34 + t75 * t53) + t22 * (t6 * t55 + t63 * t
     #55) * t34 * t8) - t72 * t3 * t12 ** 2 * t32 * (t2 * t19 * t23 + t1
     #5))
      t10 = -t51 - t1
      t15 = t13 * t4
      t22 = t5 * t1 * (2 * t4 * t11 * t19 * t12 + t18)
      t12 = t4 * (t47 * (t69 * (t22 + t54) * t30 + t13 * t36 * (t43 * t1
     #7 * (t14 * t19 - t22 - t70) + t54 + t22) - t68 * t7 * t14 * t30 **
     # 2) * t12 * t72 + t13 * t5 * t36 * (t71 * (t18 * (t74 * (t45 * t64
     # + 1) + t15 * (t54 * t8 + t70 * (-t45 - t8)) * t25) + t15 * t16 * 
     #t8 * (-t70 + t54)) + t15 * (t74 * t54 + t70 * (t45 * t10 - t74)) *
     # t20))
      t8 = t62 * t17 * (-t71 * t8 * t19 * t30 * t29 + t7 * t21 * (t32 * 
     #(t41 * t25 ** 2 * t28 + t11 * t42) + t9 * t20 ** 2 * t28) * t34 * 
     #t36 + t56 * (t10 * t32 - t52) * t33 * t3)
      t9 = -4
      ret = t9 * (t62 * (-t46 * t60 * t59 * t33 * t36 + t57 * t30 * t29 
     #+ (t54 * (-t46 * t21 * t59 + t7 * (t32 * (t51 + t1) + t52)) + t50 
     #* t5 * (t61 * t20 - t26 * t21)) * t33 * t3) + t36 * (t47 * t1 * t4
     #0 * t38 * t41 * t25 * t33 + t47 * t4 * t30 + t13 * (t44 * t38 * t1
     #3 * t32 + t21 * t28 - t43 * t21 * (t27 * t28 + t1)) * t18 + t49 * 
     #t48 * (t46 - t45)) + t50 * (t53 * (t7 + t17) + t30) * t3 + t57 * t
     #4 * t30 * (t56 + t54) * t35 + t54 * t28 * t30 * (-t48 * t31 + t18)
     # * t29 + t58 * (-t49 + t47 + t43 * (t32 * (t51 + t1) + t52)) * t37
     #) - 8 * t24 - 64 * t2 + 128 * t44 * t32 * t34 * t66 * t21 * t17 * 
     #t11 * t5 * t37 * (t63 * t7 + t6 * (t63 + t7)) - 32 * t12 + 16 * t8

      hjetmass_qqbgg_triangle_pmpp_s12_s123_0_dp = ret/32d0/(0,1d0)
      return

      end function
