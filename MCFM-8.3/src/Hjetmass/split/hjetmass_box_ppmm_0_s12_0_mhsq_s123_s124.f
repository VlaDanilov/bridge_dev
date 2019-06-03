
      complex*32 function hjetmass_box_ppmm_0_s12_0_mhsq_s123_s124
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i3)
      t2 = za(i1, i4)
      t3 = zb(i3, i1)
      t4 = zb(i4, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = za(i1, i2)
      t10 = za(i3, i4)
      t11 = zb(i2, i1)
      t12 = zb(i4, i3)
      t13 = t1 * t4
      t14 = t5 * t8
      t15 = t14 + t13
      t16 = t10 ** 2
      t17 = t16 * t12 ** 2
      t18 = t17 * t15
      t19 = t2 * t4
      t20 = t7 * t8
      t21 = t10 * t12
      t22 = t21 * t9 * t11
      t23 = t5 * t6
      t24 = t21 * (-t23 * (t20 + t19) + t22 + (-t20 - t19) * t3 * t1)
      t25 = t7 * t6
      t26 = t2 * t3
      t27 = t25 + t26
      t28 = mt ** 2
      t29 = 16 * t21 * t28 * t27 * t18 + 4 * t24 ** 2
      t29 = sqrt(t29)
      t18 = 0.1q1 / t18
      t30 = t14 + t13
      t30 = t1 * t30 * t3 + t23 * t30
      t31 = t20 + t19
      t23 = t31 * t3 * t1 + t23 * t31 - t22
      t15 = t21 * t29 * t18 * t15 / 4
      t13 = t17 * t18 * (t19 * t30 + t20 * t30 + t22 * (-t14 - t13)) / 2
      t17 = -t23 - t15 + t13
      t19 = t21 * t27
      t13 = -t23 + t15 + t13
      t15 = -2 * t24
      t20 = t15 - t29
      t21 = 0.1q1 / t12
      t21 = t31 * t21
      t22 = t21 * t3
      t23 = t20 * t18 / 4
      t24 = t23 * t10
      t27 = t24 * t4 - t22
      t15 = t15 + t29
      t29 = t21 * t6
      t30 = t15 * t18 / 4
      t31 = t30 * t10
      t32 = t31 * t8 - t29
      t33 = 0.1q1 / t10
      t19 = 0.1q1 / t19
      t21 = t21 * t33
      t33 = t21 * t5
      t34 = t13 * t19
      t35 = -t14 * t30 + t6 * (t34 * t7 + t33)
      t24 = t24 * t8 - t29
      t29 = t17 * t19
      t14 = t6 * (t29 * t7 + t33) - t14 * t23
      t4 = t31 * t4 - t22
      t21 = t21 * t1
      t22 = t6 * (t29 * t2 + t21) - t23 * t1 * t8
      t21 = t6 * (t34 * t2 + t21) - t1 * t30 * t8
      t23 = 0.1q1 / t35
      t5 = 0.1q1 / t5
      t8 = 0.1q1 / t8
      t9 = 0.1q1 / t9
      t14 = 0.1q1 / t14
      t30 = 0.1q1 / t1
      t31 = 0.1q1 / t6
      t33 = t19 ** 2
      t34 = t13 ** 2
      t35 = t17 ** 2
      t36 = t11 ** 2
      t37 = t11 * t14
      t38 = t11 * t23
      t39 = t34 * t32
      t40 = t35 * t24
      t41 = t35 * t14
      t42 = t34 * t23
      t43 = t6 * t11
      t44 = t28 * t5
      t45 = t17 * t24
      t46 = t45 * t14
      t47 = t13 * t32
      t48 = t7 * t11
      t24 = t3 * t24
      t49 = t3 * t32
      t50 = t25 * t11
      t51 = t13 * t34 * t23
      t52 = t17 * t35 * t14
      t53 = t30 * t5
      t54 = t53 * t16
      t55 = t44 * t10
      t56 = t17 + t13
      t57 = -t17 - t13
      t58 = t2 * t11
      t59 = t38 * t34
      t60 = t9 * t11
      t61 = t19 * t8
      t62 = t61 * t10
      t63 = t13 * t23
      t64 = t13 * t4
      t65 = t17 * t27
      t66 = t17 * t22
      t67 = t13 * t21
      t68 = t9 * t5
      t69 = t35 + t34
      t70 = t58 * t30
      t71 = t60 * t5
      t72 = t14 * t17
      t73 = t19 * t10
      t28 = t73 * (t10 * (t42 * t61 * t43 * (t4 * t5 + t70) + t13 * t9 *
     # (t49 * t5 * t31 + t70) + t11 * (t30 * (-t67 * t9 + t66 * (t37 - t
     #9)) - t68 * t26 * t69 * t19) * t8 - t71 * t56 * t7) + t30 * (t38 *
     # (t44 + t6) * t13 + t72 * (t5 * (t29 * t6 * t27 + t11 * t28) + t43
     #) + t59 * t5 * t19 * (t21 * t3 * t8 + t25)) * t16 + t53 * t41 * t1
     #9 * t10 * t16 * t3 * t6 - t71 * t8 * (t65 + t64) * t1 - t2 * t11 *
     # t36 * t8 * (t72 + t63))
      t21 = t28 + t19 * t16 * (t30 * (t11 * (t63 * (t11 * t21 * t8 - t32
     #) - t46) + t19 * (t42 * t3 * t6 * t16 * t5 + (t41 * (t11 * (t3 * t
     #22 * t8 + t25) - t24) + t42 * (t4 * t6 - t49)) * t5 * t10 + t41 * 
     #t2 * t6 * t36 * t8) + t9 * (t57 * t10 * t6 + t58 * t17) + t26 * t4
     #3 * t10 * t5 * t8 * (t52 + t51) * t33) + t5 * (-t9 * t56 * t10 * t
     #3 + t61 * t41 * t43 * t27 - t64 * t9 - t65 * t9) + t68 * t3 * ((-t
     #67 - t66) * t8 * t11 + t45) * t31)
      t2 = t60 * t62 * (t73 * t69 * t30 * t6 * t2 + t45 + t47)
      t22 = -8
      ret = t22 * (t62 * (t54 * t51 * t33 * t6 ** 2 * t4 + t60 * (t48 * 
     #t5 * t56 * t1 + t55 * t57 - t58 * t56) + t50 * t5 * t10 * (t9 * (-
     #t35 - t34) - t59) * t19) + t10 * (t8 * (t54 * t6 * (t52 * (t6 * (t
     #48 + t27) + t24) + t51 * (t50 + t49)) * t19 * t33 + t11 * (-t11 * 
     #(t47 * t23 + t46) + (t47 + t45) * t31 * t9 * t5 * t3 * t1) * t19 +
     # t10 * (t5 * (t3 * (t40 * (-t37 - t9) - t39 * (t38 + t9)) + t6 * (
     #t9 * (-t27 * t35 - t34 * t4) - t41 * t7 * t36)) + t43 * (t40 * t14
     # + t39 * t23) * t30 + t44 * t43 * t30 * (t42 + t41) * t10) * t33) 
     #+ t55 * t11 * t9 * t31)) - 12 * t2 + 4 * t21 - t61 * t5 * t18 * t1
     #2 * t36 * t16 * (t73 * (t42 * t15 + t41 * t20) + (-t13 * t15 - t17
     # * t20) * t31 * t9 * t1) + 16 * t33 * t16 * t6 * t9 * t30 * t8 * (
     #t40 + t39)

      hjetmass_box_ppmm_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function
