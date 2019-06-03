
      complex*32 function hjetmass_triangle_pmpm_s23_s123_0
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

      t1 = za(i1, i2)
      t2 = za(i1, i3)
      t3 = za(i2, i3)
      t4 = za(i2, i4)
      t5 = zb(i2, i1)
      t6 = zb(i3, i1)
      t7 = zb(i3, i2)
      t8 = t1 * t5
      t9 = t2 * t6
      t10 = t8 + t9
      if ( real(t10) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t10 = cg * sqrt(t10 ** 2) + t8 + t9
      t11 = za(i1, i4)
      t12 = 0.1q1 / t10
      t13 = -2 * t8 * t12
      t14 = t3 * t7
      t15 = t14 * (1 + t13)
      t16 = zb(i4, i1)
      t17 = zb(i4, i3)
      t18 = t3 * (-2 * t1 * t7 * t16 * t12 - t17)
      t19 = t9 * t14
      t13 = t5 * (t1 * t15 + t11 * t18) + t19 * t13
      t20 = zb(i4, i2)
      t21 = za(i3, i4)
      t22 = t3 * t11
      t23 = t7 * (-2 * t22 * t6 * t12 + t4)
      t18 = 0.1q1 / t18
      t13 = 0.1q1 / t13
      t15 = 0.1q1 / t15
      t10 = 0.1q1 / t10
      t24 = 0.1q1 / t16
      t25 = t4 * t20
      t26 = t17 * t21
      t27 = t26 + t25
      t28 = t5 ** 2
      t29 = t5 * t28
      t30 = t18 ** 2
      t31 = t6 ** 2
      t32 = t31 ** 2
      t33 = t6 * t31
      t34 = t1 ** 2
      t35 = t1 * t34
      t36 = t3 ** 2
      t37 = t36 ** 2
      t38 = t3 * t36
      t39 = t13 ** 2
      t40 = t13 * t39
      t41 = t11 ** 2
      t42 = t7 ** 2
      t43 = t7 * t42
      t44 = t9 * t10
      t45 = t14 * t34
      t46 = t45 * t28
      t47 = t10 * t15
      t48 = t24 * t12
      t49 = t48 * t29 * t1
      t50 = t24 * t10
      t51 = t16 * t18
      t52 = t5 * t15
      t53 = t14 * t10
      t54 = t12 * t27
      t55 = t12 * t30 * t1
      t56 = t47 * t12
      t57 = t56 * t42
      t58 = t57 * t32 * t37 * t34
      t59 = t9 * t24
      t60 = t41 * t28 * t39
      t61 = 0.1q1 / t11
      t62 = t1 * t21
      t63 = -t62 + t22
      t64 = t59 + t11
      t65 = t4 ** 2
      t66 = mt ** 2
      t67 = t2 * t4
      t68 = t4 * t11
      t69 = t66 * t11
      t70 = t24 * t6
      t71 = t15 * t12
      t72 = t71 * t33
      t73 = t12 ** 2
      t74 = t10 ** 2
      t75 = t66 * t73
      t33 = t15 * t33 * t36
      t6 = t33 * t1 * (t3 * (-t50 * t42 * t35 * t4 * t28 ** 2 * t73 * t3
     #9 + t6 * (-t34 * t16 * t18 * t30 * t73 + (t5 * t11 * t13 - t18) * 
     #t24 * t74) * t7) - t75 * (t1 * t11 * t28 * t23 * t39 + t70 * t18) 
     #+ t60 * t36 * t42 * t10 * t1 * t6 * t73 + t70 * t34 * t11 * t38 * 
     #t29 * t43 * t73 * t40) + t72 * t36 * t1 * (t28 * (t12 * t6 * t1 * 
     #(t44 * t36 * t42 * t11 * t24 - t69 * t14 * (-2 * t9 * t12 + 1) * t
     #24 - t62 * t64 * t10 * t42 * t3) * t39 + t50 * t7 * t1 * (t14 * t6
     # * t12 * t15 * t63 - t4) * t13) + t29 * (-t57 * t34 * t3 * t4 * t2
     #4 * t13 + t10 * t12 * t42 * t3 * t34 * (-t68 + (-t67 + t22 - t62) 
     #* t24 * t6) * t39 + t48 * t45 * t6 * t11 * (t17 ** 2 * t21 ** 2 + 
     #t65 * t20 ** 2) * t40) + t55 * t66 * t23 * t61 + t70 * (t63 * t10 
     #* t7 + t69 * t12) * t13 * t5)
      t17 = t71 * t7 * t32 * t38 * t34 * (t11 * (t75 * t28 * t15 * t24 *
     # t13 + t28 * t24 * (-t10 * t27 + t75 * t8 + t53 * (-t54 * t8 * t15
     # - 1) - t57 * t8 * t36) * t39 + t49 * (t14 * t27 + t26 * t25) * t4
     #0) + t18 * (-t10 * t18 + t5 * (-t56 * t14 * t1 * t18 - t75 * t15 *
     # t24)))
      t21 = t65 * t7
      t30 = t3 * t31
      t21 = t12 * t30 * t1 * ((t21 * t15 + t70 * t4 * (t7 * t15 * t61 * 
     #(t67 - t62) + 1) + t22 * t31 * t15 * t24) * t13 * t5 - t15 * t18 *
     # (t21 * t61 + t30 * t24))
      t22 = 0.1q1 / t7
      t22 = t23 * t22
      t2 = t72 * t3 * t1 * (t28 * (t47 * t1 * t36 * t4 * t42 * t24 * t13
     # + t14 * t1 * (t68 * t53 + t70 * (t4 * (t12 * t63 * t20 + t53 * t2
     #) + t14 * t12 * t63 + t26 * t12 * t63)) * t39) + t45 * t24 * t4 * 
     #(t14 * (-t12 + t10) - t54) * t39 * t29 + t66 * t61 * t24 * t18 * (
     #-t22 + t4) + t66 * t24 * (t22 - t4) * t13 * t5)
      t20 = 256
      ret = t20 * (t58 * (t10 * (t18 * (-t52 * t24 - t18) + t60) + t28 *
     # (t5 * (t41 * (t47 * t14 * t39 - t14 * t40) - t59 * t27 * t40 * t1
     #1) - t53 * t18 * t15 ** 2 * t24) * t12 * t1) + t58 * (t11 * (t50 *
     # t28 * (t46 * t12 * t15 + t9 + t8 * (t19 * t12 * t15 + 1)) * t39 +
     # t47 * t28 * t24 * (t14 * t8 * t12 * t15 + 1) * t13 + t49 * (t19 *
     # (t44 - 1) + t8 * (-t26 - t14 - t25) + t46 * t10) * t40) - t54 * t
     #1 * t29 * t40 * t41 + t53 * t1 * t11 * t41 * t29 * t16 * t12 * t40
     # + t55 * (t53 * (-t52 - t51) - t51))) + 64 * t6 + 128 * t17 - 8 * 
     #t21 + 512 * t74 * t15 * t40 * t73 * t11 * t43 * t32 * t29 * t3 * t
     #37 * t35 * (t9 * t11 + t8 * t64) - 32 * t2 + 16 * t33 * t8 * t24 *
     # t13 * t7 * t4 * (t8 * t12 * (t13 * (t26 + t14 + t25) - t12) - t10
     #)

      hjetmass_triangle_pmpm_s23_s123_0 = ret/32q0/(0,1q0)
      return

      end function
