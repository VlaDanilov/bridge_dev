
      double complex function hjetmass_triangle_pppm_s14_s124_0_dp 
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

      t1 = za(i2, i4)
      t2 = zb(i2, i1)
      t3 = zb(i3, i2)
      t4 = za(i1, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = zb(i4, i2)
      t8 = t4 * t2
      t9 = t1 * t7
      t10 = t8 + t9
      if ( dreal(t10) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t8 = cg * cdsqrt(t10 ** 2) + t8 + t9
      t10 = zb(i4, i3)
      t11 = 0.1D1 / t8
      t12 = 2 * t4
      t13 = t5 * (t12 * t3 * t6 * t11 + t10)
      t12 = t12 * t2 * t11
      t14 = t5 * t6
      t15 = t14 * (1 - t12)
      t16 = za(i2, i3)
      t17 = t4 * t15
      t12 = t9 * t12 * t14
      t18 = -t12 + t2 * (-t13 * t16 + t17)
      t19 = za(i3, i4)
      t20 = zb(i3, i1)
      t21 = t6 * (2 * t5 * t16 * t2 * t11 + t19)
      t12 = -t4 * (-t15 * t2 + t21 * t3) - t12
      t22 = za(i1, i3)
      t23 = t5 * (-2 * t1 * t3 * t6 * t11 + t20)
      t24 = t22 * t20
      t25 = t19 * t10
      t26 = t3 * t16
      t27 = 0.1D1 / t4
      t28 = 0.1D1 / t13
      t18 = 0.1D1 / t18
      t16 = 0.1D1 / t16
      t29 = 0.1D1 / t6
      t30 = 0.1D1 / t19
      t31 = 0.1D1 / t5
      t32 = t24 * t29
      t33 = t32 + t5
      t32 = -t32 - t5
      t34 = t18 ** 2
      t35 = t18 * t34
      t36 = t7 ** 2
      t37 = t28 ** 2
      t38 = t2 ** 2
      t39 = t2 * t38
      t40 = t1 ** 2
      t41 = t3 * t15
      t42 = t41 * t1
      t43 = t7 * t20
      t44 = t15 * t18
      t45 = t3 * t18
      t46 = t5 * t30
      t47 = t13 * t34
      t48 = t10 * t29
      t49 = t16 * t29 * t30
      t50 = t43 * t5
      t51 = t1 * t3
      t52 = t20 * t29
      t53 = t28 * t16
      t54 = t27 * t1
      t55 = t1 * t13
      t8 = 0.1D1 / t8
      t56 = t26 * t27
      t57 = t56 + t2
      t58 = t13 * t27
      t59 = t46 * t40 * t36 * t27
      t60 = t13 * t57
      t57 = t9 * (t46 * t57 + t58) + t59 + t60
      t61 = t2 * t13
      t62 = t38 * t34
      t63 = t27 * t16
      t64 = t8 * t40
      t28 = t64 * (t45 * t27 * t2 * (t41 + t61 - t50) + t62 * (t41 * t57
     # + t61 * t57 + t50 * (t9 * (t46 * (-t56 - t2) - t58) - t59 - t60))
     # + t63 * t46 * (-t3 ** 2 * t15 ** 2 * t37 - t38 + t50 * t28 * (t41
     # * t28 + t2)))
      t12 = 0.1D1 / t12
      t21 = 0.1D1 / t21
      t56 = 0.1D1 / t15
      t4 = t4 * t3 * t12 + t21
      t12 = t11 ** 2
      t6 = t5 ** 2 * t6
      t21 = t2 * t40
      t3 = t21 * (t14 * t56 * t12 * t38 * (t46 * t4 * t2 - t4 * t3) + (-
     #t53 * t46 * t3 * t8 + t9 * t13 * (t30 * (t20 ** 2 * t22 ** 2 * t29
     # + t6) + t10 ** 2 * t19 * t29) * t35 * t38 + t50 * (t30 * t33 + t4
     #8) * t34 * t2) * t27 * t15)
      t57 = t26 + t9
      t58 = mt ** 2
      t59 = t8 ** 2
      t65 = t57 * t27
      t66 = t6 * t30
      t67 = t66 * t59
      t68 = t67 * t1
      t69 = t5 * t8
      t70 = t58 * t23 * t11 * t29
      t71 = t27 * t30
      t6 = t21 * t7 * (t2 * (-t67 * t41 * t54 * t57 * t34 + t71 * (t6 * 
     #t51 * t59 - t70) * t18) + t38 * (t69 * t54 * t15 * t13 * (t24 * t5
     #7 * t30 + t10 * t57 + t14 * t30 * (-t64 * t36 + t9 + t26 * (-t26 *
     # t8 + 1))) * t35 + t68 * (t65 * t13 - t41) * t34) + t39 * (t68 * t
     #47 + t69 * t55 * t15 * (t30 * (t14 + t24) + t10) * t35) - t17 * t6
     #7 * t55 * t38 ** 2 * t35 - t71 * t70 * t53)
      t17 = t54 * t18
      t8 = t40 * (t71 * t29 * (t1 * (t7 * (t2 * (-t47 * t38 * t15 - t61 
     #* t27 * t18 - t63) + t41 * (-t37 * t16 ** 2 + t62) * (-2 * t26 * t
     #14 * t11 + t24 + t25)) + t62 * t15 * t23 * t36) + t41 * t2 * (t18 
     #* t2 + t53)) * t11 * t58 + t18 * t7 * t38 * (t17 * (t69 * t41 + t6
     #1 * (t44 * t32 - t69)) * t10 + t46 * (t20 * (t17 * (t41 * t8 + t61
     # * (-t44 - t8)) * t22 + t69 * (t44 * (t65 + t2) - t27)) + t54 * t1
     #4 * t8 * t18 * (t41 - t61))))
      t4 = t66 * t39 * t40 * (-t58 * t2 * t11 * t12 * t56 * t4 + t61 * t
     #26 * t9 * t59 * t15 * t35 + t60 * t35 * t59 * t15 * t36 * t40)
      t11 = -4
      ret = t11 * (t54 * t2 * (t2 * (-t52 * t13 * t18 + t9 * (-t42 * (t4
     #6 + t48) + t52 * (t9 * t5 * t10 + t24 * t13)) * t34) + t38 * (-t55
     # * t29 * t30 * t18 + t13 ** 2 * t1 * (-t29 * t31 * (t25 + t24) - 1
     #) * t34) - t53 * t42 * t30 * t29) + t54 * (t49 * t42 * (-t41 + t50
     #) * t37 + t1 * (t46 * t1 * t20 * t34 * t33 * t36 - t41 * t13 * t34
     # + t43 * t34 * (-t42 * t22 * t29 * t30 + t13 * t5) + t45 * t13 * t
     #29 * (-t44 * t24 + 1) * t31 + t48 * t47 * (-t41 * t31 + t43) * t19
     # - t49) * t38 + t52 * (-t45 * (t9 + t15) - t16) * t2 - t53 * t41 *
     # t29 * (t51 * t31 + t20) + t47 * t7 * t40 * (t30 * t32 - t48) * t3
     #9)) - 8 * t28 - 16 * t3 + 64 * t6 + 32 * t8 - 128 * t4

      hjetmass_triangle_pppm_s14_s124_0_dp = ret/32d0/(0,1d0)
      return

      end function
