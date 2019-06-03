
      complex*32 function hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = zb(i4, i1)
      t5 = zb(i4, i2)
      t6 = t1 * t5
      t7 = t3 * t4 + t6
      t8 = za(i1, i2)
      t9 = zb(i2, i1)
      t10 = za(i1, i4)
      t11 = za(i2, i4)
      t12 = t8 * t9
      t13 = t10 * t4
      t14 = t11 * t5 + t12 + t13
      t15 = zb(i3, i2)
      t16 = za(i3, i4)
      t17 = zb(i4, i3)
      t18 = t3 * t2
      t19 = t1 * t15
      t20 = t16 * t17
      t21 = t18 + t19 + t20
      if ( real(t21) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t20 = cg * sqrt(t21 ** 2) + t18 + t19 + t20
      t21 = 0.1q1 / t20
      t12 = -2 * t14 * t18 * t21 + t12 + t13
      t13 = 2 * t3
      t18 = t2 * (t1 * (-t13 * t14 * t15 * t21 + t10 * t5) + t12 * t3)
      t22 = t1 * t9 - t16 * t4
      t13 = t13 * t14 * t17 * t21 - t5 * t8
      t23 = -2 * t1 * t14 * t2 * t21 + t11 * t4
      t24 = t15 * t23
      t25 = t3 * (t12 * t2 + t24)
      t26 = 2 * t14 * t16 * t2 * t21 - t11 * t9
      t27 = 0.1q1 / t8
      t28 = 0.1q1 / t9
      t29 = 0.1q1 / t16
      t30 = 0.1q1 / t17
      t11 = 0.1q1 / t11
      t18 = 0.1q1 / t18
      t25 = 0.1q1 / t25
      t26 = 0.1q1 / t26
      t31 = t3 * t25
      t32 = t30 * t26
      t33 = t32 + t31
      t34 = t25 ** 2
      t35 = t25 * t34
      t36 = t3 ** 2
      t37 = t36 ** 2
      t38 = t3 * t36
      t39 = t23 ** 2
      t40 = t1 ** 2
      t41 = t4 ** 2
      t42 = t1 * t12
      t43 = t28 * t27
      t44 = t1 * t11
      t45 = t44 * t23
      t46 = t39 * t11
      t47 = t16 * t11
      t48 = t47 * t23
      t49 = t43 * t1
      t50 = t42 * t18
      t51 = t31 * t23
      t52 = t23 * t26
      t53 = t22 * t28
      t54 = t40 * t27 * t11
      t55 = 0.1q1 / t2
      t20 = 0.1q1 / t20
      t10 = 0.1q1 / t10
      t13 = 0.1q1 / t13
      t56 = 0.1q1 / t14
      t57 = mt ** 2
      t13 = t29 * t55 * t13
      t58 = t22 * t1
      t59 = t40 * t10 * t56
      t60 = t8 * t11
      t61 = t20 ** 2
      t62 = t2 ** 2
      t63 = t26 ** 2
      t64 = t14 ** 2
      t65 = t46 * t43
      t66 = t17 * t55
      t67 = t16 * t63
      t68 = t43 * t11
      t66 = t58 * (t68 * (-t19 * t36 * t17 * t34 + t67) * t61 * t23 * t6
     #4 * t2 + t56 * t21 * t57 * (t10 * (t66 + t60) + (t66 * t27 + t11) 
     #* t28 * t4) + (-t43 * t38 * t17 * t11 * t61 * t34 * t23 + t65 * t3
     #7 * t15 * t17 * t61 * t35) * t64 * t62)
      t69 = t31 * t40
      t70 = t39 * t34
      t63 = t16 ** 2 * t39 * t63
      t71 = t1 * t57
      t19 = t68 * t4 * ((t2 * (t36 * (t24 * t40 * t12 * t34 + t1 * t23 *
     # t25) - t69 * t12 - t70 * t19 * t38) + t62 * (t42 * t23 * t34 * t3
     #8 - t70 * t37) - t40 - t63) * t20 * t14 + t71 * (t23 * (t32 + t31)
     # + t42 * (t13 + t18)))
      t24 = t53 * t20 * t14 * t4 * (t2 * (t45 * t15 * t34 * t36 - t44 * 
     #t31) + t27 * t26 * (t52 * t16 + t1) + t38 * t62 * t23 * t11 * t34)
      t32 = t53 * t27
      t15 = t2 * (t45 * t36 * t17 * t34 * (-t51 * t15 + 1) * t20 * t22 *
     # t14 + t32 * (t26 * (t63 + t40) + t40 * t38 * t15 ** 2 * t17 * t39
     # * t35 + t69 * t17) * t61 * t11 * t64) + t46 * t53 * t64 * t61 * t
     #3 * t37 * t2 * t62 * t17 * t27 * t35 - t46 * t14 * t20 * t37 * t62
     # * t22 * t17 * t35 + t71 * t4 * t55 * t56 * (t43 * t4 + t10)
      t16 = t65 * t21 * t2 * t57 * (t67 * t7 * t30 + t36 * t17 * t34 * (
     #t1 * (t16 * t5 + t3 * t9) - t22 * t3))
      t30 = -4
      t37 = 16
      ret = (t4 * (t2 * (t44 * t28 * (t18 * t27 * t42 + t22 * t25) * t3 
     #- t46 * t34 * t38 + t45 * t25 * (t12 * t25 + t43) * t36) + t49 * (
     #t26 * (-t48 - t22) + t44)) - t53 * t33 * t41 + t54 * (-t1 * t29 + 
     #t52) + t54 * t2 * t28 * (t51 + t50) * t5) * t30 + (t1 * t28 * t56 
     #* (t27 * t3 - t47) * t41 + t59 * t5 + t49 * t11 * (t58 * ((t13 + t
     #18) * t12 * t17 - t29) + (t23 * t33 + t50) * t7 * t2) * t21 * t57 
     #+ t1 * (t56 * (t10 * (-t47 * t8 + t3) + t44) + t43 * (-t48 * t14 *
     # t20 * t26 + t6 * t56)) * t4 + t60 * (-t2 * t38 * t17 * t22 * t39 
     #* t35 + t59) * t9) * t37 - 128 * t66 + 8 * t19 + 24 * t24 - 64 * t
     #15 + 32 * t16 + 48 * t32 * t45 * t57 * t21 * (t31 * t17 + t26) - 1
     #2 * t60 * t36 * t2 * t4 * t22 * t23 * t34

      hjetmass_qqbgg_triangle_pmmp_s124_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function
