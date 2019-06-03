
      complex*32 function hjetmass_triangle_ppmm_s14_s134_0
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
      t3 = za(i2, i4)
      t4 = za(i3, i4)
      t5 = zb(i2, i1)
      t6 = zb(i3, i1)
      t7 = zb(i3, i2)
      t8 = zb(i4, i1)
      t9 = zb(i4, i2)
      t10 = zb(i4, i3)
      t11 = t2 * t6
      t12 = t4 * t10
      t13 = t11 + t12
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t13 = cg * sqrt(t13 ** 2) + t11 + t12
      t14 = za(i2, i3)
      t15 = za(i1, i4)
      t16 = 0.1q1 / t13
      t17 = 2 * t11 * t16
      t18 = t15 * t8
      t19 = t18 * (1 - t17)
      t20 = t8 * (-2 * t15 * t14 * t6 * t16 + t3)
      t17 = -t12 * t17 * t18 + t2 * (t19 * t6 + t20 * t7)
      t21 = t4 * t7
      t22 = t15 * (2 * t21 * t8 * t16 + t5)
      t23 = t4 * t9
      t24 = t2 * t5 + t23
      t25 = 0.1q1 / t14
      t17 = 0.1q1 / t17
      t13 = 0.1q1 / t13
      t20 = 0.1q1 / t20
      t10 = 0.1q1 / t10
      t26 = t13 ** 2
      t27 = t7 ** 2
      t28 = t20 ** 2
      t29 = t4 ** 2
      t30 = t4 * t29
      t31 = mt ** 2
      t32 = t17 ** 2
      t33 = t17 * t32
      t34 = t8 ** 2
      t35 = t2 ** 2
      t36 = t35 ** 2
      t37 = t2 * t35
      t38 = t6 ** 2
      t39 = t6 * t38
      t40 = t3 * t9
      t41 = t1 * t5
      t42 = t18 * t26
      t43 = t10 * t16
      t44 = t43 * t25
      t45 = t44 * t6
      t46 = t25 * t6
      t47 = t46 * t10
      t48 = t42 * t32
      t49 = t18 + t41 + t40
      t50 = t18 * t13
      t51 = t21 * t25
      t52 = t7 * t10
      t53 = t52 * t50
      t54 = t16 * t13 * t15
      t55 = 0.1q1 / t7
      t56 = 0.1q1 / t15
      t57 = t15 ** 2
      t58 = t5 ** 2
      t59 = t4 * t25
      t60 = t7 * t8
      t61 = t55 * t20
      t62 = t38 * t29
      t63 = t40 + t41
      t64 = t6 * t9
      t65 = t2 * t17
      t66 = t8 * t10
      t1 = t66 * t6 * t29 * (t35 * (-t45 * t31 * t7 * t56 * t17 + t7 * t
     #6 * (-t31 * t22 * t56 * t16 + (-t64 * t15 * t16 + t25 * t63) * t13
     # * t8 * t4) * t32) - t23 * t50 * t37 * t39 * t16 * t25 * t32 + t59
     # * t31 * t16 * t56 * t20 - t65 * t51 * t31 * t16 * t56) + t62 * (t
     #34 * (t37 * (t45 * t21 * (-t3 ** 2 * t9 ** 2 - t58 * t1 ** 2 - t34
     # * t57) * t33 + t54 * t6 * (-t52 * t5 + t59 * (-t52 * t8 - t5)) * 
     #t32) + t4 * t28 * t10 * (t14 * t6 * t16 * t20 + t13) - t44 * t15 *
     # t13 * t36 * t5 * t38 * t32 + t13 * t15 * t4 * (-t23 * t46 * t16 +
     # t60 * t25 * (-t16 * t4 + t10) - t43 * t27 * t8) * t32 * t35) + t4
     #3 * (t51 * t35 * t8 * t18 * (-2 * t12 * t16 + 1) * t32 + t25 * (t8
     # * t20 * t10 - t17 * t24) * t2 + t61 * (t8 * t22 * t20 + t24 * t25
     #)) * t56 * t31)
      t3 = t64 + t60
      t9 = t31 * t5 * t56
      t15 = t66 * t25
      t3 = t6 * t29 * (t35 * (-t47 * t5 * t8 * t13 * t17 + t8 * t6 * (-t
     #53 * t5 + t59 * (-t50 * t5 + t43 * (t18 * t3 + t40 * t3 + t41 * t3
     #))) * t32) + t15 * t38 * t5 * (t16 * t63 + t18 * (-t13 + t16)) * t
     #32 * t37 + t65 * t10 * t25 * (-t8 * t4 * t13 * t3 + t9) - t9 * t61
     # * t25 * t10)
      t9 = t15 * t62 * t17 * t5 * t35 * (t17 * t49 - t16)
      ret = -64 * t34 * t38 * t30 * (t35 * (t48 * t10 * t27 + t48 * t21 
     #* t25) - t42 * t28 * t10 + (t47 * (t16 ** 2 * t31 + t42) * t32 + t
     #45 * (t40 * t18 + t41 * (t18 + t40)) * t33) * t7 * t37) + 128 * t5
     #4 * t8 * t34 * t39 * t30 * (t37 * (t49 * t10 * t33 * t27 + t51 * (
     #t18 * (-t12 * t13 + 1) + t41 + t40) * t33 - t50 * t14 * t33 * t10 
     #* t7 * t27) + t14 * t20 * t28 * t10 * (t50 + 1) - t53 * t2 * t36 *
     # t38 * t25 * t33 + t47 * t49 * t33 * t7 * t36) + 32 * t1 - 4 * t6 
     #* t4 * (t66 * t56 * (-t65 * t7 + t20) * t25 * t29 + t58 * t10 * (-
     #t65 + t61) + t65 * t59 * t5 * (t10 * (t19 * t56 + t64 * t55) - t5 
     #* t55)) + 16 * t3 - 256 * t16 * t26 * t33 * t57 * t34 ** 2 * t7 * 
     #t39 * t30 * t37 * (t11 * (t59 + t52) + t21) + 8 * t9

      hjetmass_triangle_ppmm_s14_s134_0 = ret/32q0/(0,1q0)
      return

      end function
