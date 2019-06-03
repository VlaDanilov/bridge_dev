
      complex*32 function hjetmass_box_pmpm_0_0_s23_mhsq_s14_s234
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i4)
      t2 = za(i2, i3)
      t3 = zb(i3, i2)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = zb(i2, i1)
      t10 = zb(i3, i1)
      t11 = t7 * t10
      t12 = t5 * t9 + t11
      t13 = mt ** 2
      t14 = za(i1, i2)
      t15 = za(i1, i3)
      t16 = t14 * t6
      t17 = t15 * t8
      t18 = t17 + t16
      t19 = t1 * t4
      t20 = t19 * t18
      t21 = t2 * t3
      t6 = t5 * t6
      t22 = t7 * t8
      t23 = t19 * (t22 + t21 + t6)
      t24 = 16 * t13 * t12 * t20 + 4 * t23 ** 2
      t24 = sqrt(t24)
      t23 = 2 * t23
      t25 = t24 - t23
      t20 = 0.1q1 / t20
      t26 = t22 + t6
      t27 = t19 * t20
      t16 = t27 * (t16 * t26 + t21 * (t17 + t16) + t17 * t26)
      t17 = (0.1q1 / 0.2q1)
      t6 = 2 * t22 + 2 * t21 + 2 * t6
      t18 = t17 * t24 * t20 * t18
      t21 = t18 - t16 + t6
      t12 = 0.1q1 / t12
      t22 = (0.1q1 / 0.4q1)
      t26 = t17 * t21 * t12
      t28 = t22 * t25 * t20
      t29 = t26 * t5 * t10 - t28 * t14 * t8
      t26 = t8 * (t28 * t15 + t7) - t11 * t26
      t23 = -t24 - t23
      t6 = -t18 - t16 + t6
      t7 = t8 * (t15 * t23 * t20 * t22 + t7) - t11 * t17 * t6 * t12
      t11 = t17 * t6 * t12 * t5 * t10 - t23 * t22 * t20 * t14 * t8
      t16 = 0.1q1 / t2
      t18 = 0.1q1 / t14
      t7 = 0.1q1 / t7
      t24 = 0.1q1 / t26
      t26 = 0.1q1 / t8
      t9 = 0.1q1 / t9
      t3 = 0.1q1 / t3
      t28 = t29 * t24
      t30 = t11 * t7 + t28
      t31 = t16 * t3
      t32 = t31 + t7
      t33 = t7 + t24
      t34 = t10 ** 2
      t35 = t25 * t24
      t36 = t23 * t7
      t37 = t34 * t26
      t38 = t25 * t29
      t39 = t23 * t11
      t40 = t27 * t18
      t41 = t9 * t10
      t42 = t41 * t5
      t43 = t11 * t26
      t44 = t43 + t5
      t45 = t29 * t26
      t46 = t45 + t5
      t47 = t20 ** 2
      t48 = t25 ** 2
      t49 = t25 * t48
      t50 = t23 ** 2
      t51 = t50 * t7
      t52 = t48 * t24
      t53 = t52 + t51
      t54 = -t31 - t24
      t55 = t11 ** 2
      t56 = t5 ** 2
      t57 = t5 * t56
      t58 = t29 ** 2
      t59 = t4 * t16
      t60 = t23 * t6
      t61 = t21 * t24
      t62 = t10 * t56
      t63 = t4 ** 2
      t64 = t47 * t1
      t65 = t36 + t35
      t66 = t39 * t7
      t67 = t34 * t3
      t68 = t67 * t56 * t18
      t69 = t62 * t18
      t70 = t5 * t8
      t71 = t25 + t23
      t72 = t38 * t24
      t2 = t20 * (t4 * (-t65 * t8 * t16 * t57 + (t41 * t71 + t66 + t72) 
     #* t16 * t56) + t9 * t10 * t34 * t1 * (t2 * t65 + t3 * t71)) + t20 
     #* t9 * t4 * (t1 * (t16 * (t26 * ((t25 * (t18 * t24 * t58 + t28 * t
     #10) + t66 * (t11 * t18 + t10)) * t5 * t4 + t67 * (t39 + t38)) + t6
     #8 * (t21 * t25 + t60) * t12) + t34 * (t35 * t46 + t36 * t44) + t69
     # * t65 * t8) - t37 * t5 * t16 * t65 * t14 * t13 + t69 * t16 * (t36
     # * (-t70 + t11) + t35 * (-t70 + t29)) * t15)
      t15 = t23 * t50 * t7
      t21 = t64 * t9 * t16 * t63 * (t62 * (-t52 * t21 - t51 * t6) * t12 
     #+ t20 * t14 * (t15 * t11 + t28 * t49))
      t28 = t42 * t13 * (-t59 * t56 * t18 * t33 + t18 * (t26 * t30 + t5 
     #* (-t7 - t24) + t31 * t26 * (t29 + t11)) * t10 + t31 * t37)
      t34 = t66 + t72
      t6 = t9 * t4 * (t40 * t56 * t16 * t34 + t37 * t13 * t57 * t12 * t1
     #6 * t18 * (t6 * t7 + t61) + t20 * (t5 * t13 * t16 * t34 + (-t23 * 
     #t55 * t32 + t54 * t58 * t25) * t18 * t1) * t26 * t10)
      ret = -t17 * t27 * t9 * (t20 * (t10 * (t53 * t14 * t10 - t50 * t11
     # * t32 + t54 * t29 * t48) + t16 * (t51 * t55 + t52 * t58) * t26 * 
     #t4 + t59 * t53 * t8 * t56) + t62 * (t61 * (t45 * t4 * t16 + t10) *
     # t25 + t60 * t7 * (t43 * t4 * t16 + t10)) * t18 * t12) - t22 * t41
     # * t19 * t47 * t16 * t14 * (t4 * (t51 * t44 + t52 * t46) + (t50 + 
     #t48) * t3 * t10) + (0.7q1 / 0.4q1) * t64 * t5 * t63 * t16 * t9 * (
     #t51 * t11 + t52 * t29) - t41 * t20 * t47 * t14 ** 2 * t63 * t1 * t
     #16 * (t24 * t49 + t15) / 16 + t21 / 8 + 8 * t28 + 4 * t42 * (t13 *
     # (t37 * t33 + (t20 * (-t36 - t35) + t26 * (t10 * t33 + t18 * t30))
     # * t16 * t5 * t4) + t40 * (t39 * t32 + t38 * (t31 + t24))) - t2 + 
     #2 * t6 - 16 * t68 * t13 * t16 * t9

      hjetmass_box_pmpm_0_0_s23_mhsq_s14_s234 = ret/32q0/(0,1q0)
      return

      end function
