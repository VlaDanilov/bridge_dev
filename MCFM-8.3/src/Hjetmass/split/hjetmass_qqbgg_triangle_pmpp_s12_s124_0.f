
      complex*32 function hjetmass_qqbgg_triangle_pmpp_s12_s124_0
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

      t1 = za(i1, i3)
      t2 = za(i2, i4)
      t3 = zb(i3, i1)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = zb(i2, i1)
      t7 = za(i1, i4)
      t8 = zb(i4, i2)
      t9 = t7 * t4
      t10 = t2 * t8
      t11 = t9 + t10
      if ( real(t11) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t11 ** 2) + t10 + t9
      t12 = 0.1q1 / t11
      t13 = -2 * t9 * t12
      t14 = t5 * t6
      t15 = t14 * (1 + t13)
      t16 = zb(i3, i2)
      t17 = zb(i4, i3)
      t18 = t5 * (-2 * t7 * t6 * t17 * t12 - t16)
      t19 = za(i2, i3)
      t20 = za(i3, i4)
      t21 = t7 * t15
      t22 = t20 * t18
      t13 = t10 * t14 * t13
      t23 = t4 * (t21 + t22) + t13
      t24 = t6 * (-2 * t5 * t20 * t4 * t12 - t19)
      t13 = t7 * (t15 * t4 + t17 * t24) + t13
      t25 = 0.1q1 / t19
      t26 = 0.1q1 / t20
      t27 = 0.1q1 / t6
      t28 = 0.1q1 / t5
      t29 = 0.1q1 / t18
      t23 = 0.1q1 / t23
      t30 = t4 * t2
      t31 = t30 * t25
      t32 = t4 ** 2
      t33 = t18 ** 2
      t34 = t29 ** 2
      t35 = t23 ** 2
      t36 = t23 * t35
      t37 = t8 ** 2
      t38 = t17 ** 2
      t39 = t16 * t19
      t40 = t1 * t3
      t41 = t40 * t25
      t42 = t2 * t3 * t37
      t43 = t3 * t18
      t44 = t23 * t4
      t45 = t27 * (t3 + t31)
      t46 = t5 * t25
      t47 = t46 * t10
      t48 = t26 * t29
      t49 = t5 ** 2
      t50 = t15 * t23
      t51 = t8 * t3
      t52 = t15 * t17
      t53 = t4 * t18
      t24 = 0.1q1 / t24
      t13 = 0.1q1 / t13
      t11 = 0.1q1 / t11
      t54 = t2 ** 2
      t55 = t18 * t2
      t56 = t5 * t2
      t57 = t10 * t18
      t58 = t32 * t7
      t59 = t51 * t5
      t21 = t17 * t2 * (t4 * (-t58 * t46 * t12 * t13 + t4 * (-t42 * t49 
     #* t20 * t25 + t4 * t33 * t20 + (t43 * t20 + t31 * (-t22 + t21)) * 
     #t8 * t5) * t11 * t35 - t59 * t11 * t23) + t46 * t15 * t11 * (t59 -
     # t52) * t34 * t26) + t30 * (t11 * (t17 * (-t53 + t52) * t23 - t46 
     #* (t59 * t29 + t4) * t26 + t4 * (t15 * (t17 * (t46 * t54 * t37 - t
     #9 * t18 - t57) + t20 * (t47 - t18) * t38) + t8 * (t8 * (-t51 * t54
     # * t49 * t25 + t56 * ((-t5 * t3 * t7 - t55) * t25 * t4 + t43)) + t
     #53 * ((-t31 + t3) * t7 * t5 + t55)) + t58 * t33) * t35) + t4 * (t2
     #4 * (t46 * t4 - t17) + t7 * t38 * t13) * t12)
      t6 = t49 * t6 * t25
      t22 = t30 * t15 * (-t46 * t48 * t11 * t17 + t59 * (t27 * (-t41 - t
     #16) - t46) * t35 * t4 + t57 * (t27 * (t1 ** 2 * t3 ** 2 * t25 + t1
     #9 * t16 ** 2) + t6) * t36 * t32)
      t55 = t53 + t52
      t57 = mt ** 2
      t58 = t17 * t20
      t59 = t5 * t11
      t60 = t23 * t2
      t19 = t2 * (t27 * t25 * (t32 * (t44 * t18 - t26) + t52 * (t4 * (t4
     #4 - t48) + (-t32 * t18 * t35 + t29 * t26 ** 2) * t28 * (t1 * t4 + 
     #t19 * t8) * t2)) * t12 * t57 + t23 * t8 * t32 * (t46 * (t3 * (t59 
     #* (t50 * (t10 + t58 + t9) - 1) + t60 * (t52 * t11 + t53 * (t50 + t
     #11)) * t1) + t14 * t2 * t23 * t11 * t55) + t60 * (t59 * t52 + t53 
     #* (t59 + t50 * (t40 * t27 + t5))) * t16))
      t59 = t25 * (t14 + t40) + t16
      t60 = t11 ** 2
      t14 = t14 * t25
      t61 = t14 * t60
      t62 = t53 * t15
      t63 = t62 * t61
      t7 = t56 * t32 * (-t63 * t2 * t54 * t8 * t37 * t36 + t57 * t4 * t1
     #2 ** 2 * t25 * (t7 * t17 * t13 - t24) + (t62 * t59 * t11 * t36 + t
     #61 * t55 * t35) * t37 * t54 + t10 * (t36 * (t62 * (t58 * t59 + t9 
     #* t59) * t11 - t63 * (t38 * t20 ** 2 + t32 * t7 ** 2)) + t14 * (t5
     #8 * t55 + t9 * t55) * t60 * t35 - t61 * t17 * t23))
      t11 = -4
      ret = t11 * (t30 * (t27 * (-t51 * t17 * t23 + t52 * (-t40 * t10 * 
     #t4 * t35 + t48) * t25 + t44 * (t44 * t10 * t16 + (t50 * (t39 + t40
     #) - 1) * t28 * t17) * t18) + t53 * t8 * t35 * (t31 - t3) * t5 + t3
     #1 * t49 * t3 * t37 * t35) + t15 * (t17 * (-t44 * t3 * t27 + t2 * (
     #t18 - t10 * (t16 * t27 + t46)) * t35 * t32 + t48 * t3 * t27 * (t47
     # * t29 + 1)) + t48 * t2 * t27 * t28 * t38) + t4 * (t45 * t44 * t18
     # - t45 * t26 + t30 * (t43 * t27 * t8 * (t1 * (t31 - t3) - t39) + t
     #4 * (-1 - t27 * t28 * (t39 + t40)) * t33 + t42 * t27 * (t16 + t41)
     # * t5) * t35) - t2 * t15 ** 2 * t38 * t25 * t26 * t27 * t34) - 8 *
     # t21 + 16 * t22 + 32 * t19 - 64 * t7 + 128 * t6 * t60 * t36 * t18 
     #* t15 * t8 * t4 * t32 * t54 * (t58 * t9 + t10 * (t9 + t58))

      hjetmass_qqbgg_triangle_pmpp_s12_s124_0 = -ret/32q0/(0,1q0)
      return

      end function
