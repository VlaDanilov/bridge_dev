
      complex*32 function hjetmass_triangle_pmpm_s12_s124_0
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

      t1 = za(i2, i4)
      t2 = zb(i2, i1)
      t3 = zb(i4, i1)
      t4 = za(i1, i4)
      t5 = zb(i4, i2)
      t6 = t4 * t3
      t7 = t1 * t5
      t8 = t6 + t7
      if ( real(t8) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t8 = cg * sqrt(t8 ** 2) + t6 + t7
      t9 = za(i3, i4)
      t10 = za(i1, i2)
      t11 = 0.1q1 / t8
      t12 = -2 * t6 * t11
      t13 = t10 * t2
      t14 = t13 * (1 + t12)
      t15 = za(i2, i3)
      t16 = t2 * (-2 * t10 * t9 * t3 * t11 - t15)
      t17 = zb(i3, i1)
      t18 = zb(i4, i3)
      t19 = t7 * t13
      t12 = t19 * t12 + t4 * (t14 * t3 + t16 * t18)
      t20 = za(i1, i3)
      t21 = zb(i3, i2)
      t22 = t1 * t21
      t23 = t4 * t17
      t12 = 0.1q1 / t12
      t8 = 0.1q1 / t8
      t14 = 0.1q1 / t14
      t16 = 0.1q1 / t16
      t24 = 0.1q1 / t9
      t25 = 0.1q1 / t18
      t26 = t2 * t18
      t27 = t3 * t21
      t28 = t27 - t26
      t29 = t11 ** 2
      t30 = t17 ** 2
      t31 = t1 ** 2
      t32 = t31 ** 2
      t33 = t1 * t31
      t34 = mt ** 2
      t35 = t4 ** 2
      t36 = t35 ** 2
      t37 = t4 * t35
      t38 = t12 ** 2
      t39 = t12 * t38
      t40 = t10 ** 2
      t41 = t16 ** 2
      t42 = t8 ** 2
      t43 = t3 ** 2
      t44 = t3 * t43
      t45 = t2 ** 2
      t46 = t45 ** 2
      t47 = t2 * t45
      t5 = t5 * t17
      t48 = t17 * t18
      t49 = t40 * t45
      t50 = t10 * t18
      t51 = t50 * t24
      t52 = t49 * t8
      t53 = t52 * t17
      t52 = t52 * t29
      t54 = t13 * t24 * t11
      t55 = t54 * t8
      t56 = t55 * t3
      t57 = t34 * (t23 + t22) * t3 * t29
      t58 = t9 * t41
      t59 = t34 * t29
      t60 = t59 * t24
      t61 = t4 * t12
      t62 = t14 * t3
      t63 = t62 * t2 * t33
      t64 = t10 * t14
      t65 = t10 * t30
      t66 = t14 * t16
      t67 = t11 * t3
      t31 = t67 * t2 * t31 * (t61 * (t26 * t31 * t14 * t24 + t65 * t14 +
     # t17 * (t64 * t25 * (t5 - t27) + 1) * t24 * t1) - t66 * (t31 * t2 
     #* t24 + t65 * t25))
      t65 = t15 * t21
      t68 = t17 * t20
      t69 = t68 + t65
      t70 = t12 * t8
      t71 = t8 * t41
      t72 = t37 * t39
      t73 = t12 * t18
      t74 = t67 * t69
      t75 = t35 * t18
      t76 = t75 * t38
      t59 = t64 * t11 * t43 * t45 * t33 * (t1 * (t2 * (-t48 * t72 * t27 
     #* t15 * t20 * t11 * t24 + t66 * t60 * t4 + t71 + t73 * (-t59 * t14
     # + t70 * t69) * t24 * t35) + t62 * t37 * t40 * t18 * t38 * t8 * t1
     #1 * t24 * t47 + t10 * t4 * (t4 * t18 * t38 * t8 * t24 + t71 * t62 
     #* t11 + t75 * (t62 * t8 * t11 * t69 * t24 * t38 - t74 * t24 * t39)
     #) * t45) + t57 * (-t25 * t41 + t76))
      t62 = t13 + t65 + t68
      t77 = t14 ** 2
      t78 = t18 ** 2
      t79 = t7 * t8
      t80 = t14 * t11 * t8 * t40 * t43 * t46 * t32
      t4 = t80 * (t67 * (t9 * t16 * t41 + t13 * (t78 * (-t37 * t8 * t14 
     #* t38 + t72) + t16 * t8 * (t35 * t77 * t24 + t66 * t4 + t58))) + t
     #8 * t4 * (-t4 * t78 * t38 + t66 * t24)) + t80 * (t37 * (t18 * (t67
     # * t7 * t24 * (t13 * (-t79 + 1) + t65 + t68) * t39 - t8 * t3 * t24
     # * (t19 * t11 * t14 + 1) * t38 - t56 * t77 * t12) + t74 * t39 * t7
     #8 - t67 * t13 * t9 * t39 * t8 * t18 * t78) + t71 + t75 * (-t70 * t
     #14 * t24 - t79 * t24 * t38) + (t24 * t11 * t43 * t62 * t39 - t55 *
     # t43 * t14 * t38) * t18 * t36 - t55 * t4 * t36 * t44 * t18 * t39)
      t9 = t63 * t11 * (t45 * (t3 * (t51 * t1 * t35 * t38 * t11 * t69 + 
     #t70 * t17 * t35 * (t24 * (t7 * t12 + t14) + t73) * t40) + t24 * t3
     #8 * t35 * (-t22 * t11 + t23 * (t8 - t11)) * t40 * t43) + t76 * t67
     # * t40 * t1 * t47 * t24 + t34 * t17 * t24 * (t16 * t25 - t61) + t5
     #4 * t38 * t35 * (-t22 * t69 - t23 * t69) * t43)
      t19 = t6 * t24
      t22 = t19 * t64 * t12 * t17 * t45 * t33 * (t6 * t11 * (t12 * t62 -
     # t11) - t8)
      ret = -64 * t63 * (t16 * (t1 * (t60 * t2 + t10 * (t58 * t43 * t29 
     #+ t24 * t42) * t45) - t57 * t24 * t25) + t35 * (t52 * t3 * t1 * (t
     #7 * t28 * t24 + t18 * t28) * t38 + t56 * (t13 * t1 * t11 * t14 * t
     #28 + t17) * t12) + t37 * (t52 * t43 * (t48 + (t5 - t26 + t27) * t2
     #4 * t1) * t38 + t51 * t29 * t43 * t45 * t1 * (-t15 ** 2 * t21 ** 2
     # - t30 * t20 ** 2 - t49) * t39 + t53 * t43 * t29 * t24 * t14 * t12
     #) + t53 * t36 * t44 * t29 * t24 * t38 + t61 * t24 * (t1 * (t11 * (
     #t27 * t10 * t8 - t34 * t18 * t11) * t2 - t50 * t8 * (t8 + t11) * t
     #45) + t57)) - 8 * t31 - 128 * t59 - 256 * t4 - 32 * t9 + 512 * t72
     # * t14 * t29 * t42 * t18 * t10 * t40 * t44 * t2 * t46 * t32 * (t6 
     #* t18 + t7 * (t19 + t18)) + 16 * t22

      hjetmass_triangle_pmpm_s12_s124_0 = ret/32q0/(0,1q0)
      return

      end function
