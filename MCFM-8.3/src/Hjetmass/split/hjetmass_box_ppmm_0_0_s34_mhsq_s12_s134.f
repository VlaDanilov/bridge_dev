
      complex*32 function hjetmass_box_ppmm_0_0_s34_mhsq_s12_s134
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i1, i2)
      t2 = za(i1, i3)
      t3 = zb(i2, i1)
      t4 = zb(i3, i1)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = za(i2, i3)
      t10 = za(i2, i4)
      t11 = t9 * t4
      t12 = t10 * t6
      t13 = t12 + t11
      t14 = t1 * t3
      t15 = t14 * t13
      t16 = zb(i3, i2)
      t17 = zb(i4, i2)
      t18 = t5 * t17
      t19 = t16 * t2 + t18
      t20 = mt ** 2
      t21 = t5 * t6
      t22 = t2 * t4
      t23 = t7 * t8
      t24 = t14 * (t23 + t21 + t22)
      t25 = 16 * t15 * t20 * t19 + 4 * t24 ** 2
      t25 = sqrt(t25)
      t24 = 2 * t24
      t26 = -t24 + t25
      t24 = -t24 - t25
      t15 = 0.1q1 / t15
      t27 = t22 + t21
      t11 = t14 * t15 * (t11 * t27 + t12 * t27 + t23 * (t12 + t11))
      t13 = t25 * t15 * t13 / 2
      t21 = 2 * t23 + 2 * t21 + 2 * t22
      t22 = -t13 + t21 - t11
      t19 = 0.1q1 / t19
      t23 = t10 * t24 * t15 / 4
      t18 = t18 / 2
      t25 = t18 * t22 * t19
      t27 = t6 * (t23 + t5) - t25
      t28 = t2 * t10
      t29 = t5 * t9
      t30 = -t28 + t29
      t31 = t24 * t15
      t32 = t31 * t22 * t19 * t3 * t30
      t11 = t13 + t21 - t11
      t13 = t10 * t26 * t15 / 4
      t18 = t18 * t11 * t19
      t6 = t6 * (t5 + t13) - t18
      t21 = t26 * t15
      t30 = t21 * t11 * t19 * t3 * t30
      t25 = t25 - t31 * t12 / 4
      t12 = t18 - t21 * t12 / 4
      t18 = 0.1q1 / t10
      t21 = 0.1q1 / t16
      t8 = 0.1q1 / t8
      t6 = 0.1q1 / t6
      t27 = 0.1q1 / t27
      t31 = 0.1q1 / t5
      t17 = 0.1q1 / t17
      t33 = t14 * t21 + t9
      t34 = t7 * t27
      t35 = t34 + t8
      t36 = t24 ** 2
      t37 = t3 ** 2
      t38 = t37 ** 2
      t39 = t3 * t37
      t40 = t26 ** 2
      t41 = t15 ** 2
      t14 = t8 * (t16 * t9 + t14)
      t42 = t3 * t6
      t43 = t16 * t8 * t31
      t44 = t7 * t6
      t45 = t3 * t27
      t46 = t24 * t22
      t47 = t32 * t27
      t48 = t30 * t6
      t49 = t36 * t22
      t50 = t40 * t11
      t51 = t17 * t18
      t52 = t21 * t10
      t53 = t52 * t3
      t54 = t29 * t18 + t2
      t55 = t17 * (t29 + t28) * t3
      t56 = t8 * t54
      t57 = (-t17 * t9 - t52) * t8
      t58 = t15 * t39
      t59 = t5 * t37 * t21
      t60 = t50 * t6
      t61 = t49 * t27
      t62 = t61 + t60
      t63 = t50 + t49
      t64 = t41 * t17
      t28 = t64 * t39 * (-t53 * t9 * (t27 * t36 + t40 * t6) + (t28 * t63
     # * t31 + t63 * t9) * t19 * t8)
      t39 = t6 + t27
      t63 = t7 * t18
      t4 = t37 * (t21 * t7 * (t51 * (t12 * t6 + t25 * t27) + t6 + t27) *
     # t5 * t3 + t63 * (t8 * (t17 * (t25 + t12) + t21 * (t23 * t4 - t22 
     #* t19 * t5 * t16 / 2 + t13 * t4 - t11 * t19 * t5 * t16 / 2)) - t17
     # * t39 * t7 * t20) + t59 * t2 * t17 * t39)
      t13 = t22 * t27
      t23 = t11 * t6
      t1 = t58 * t17 * t21 * (t3 * (t9 * t10 * (t13 * t24 * t36 + t23 * 
     #t26 * t40) * t19 * t41 + t1 * t7 * (t61 + t60) * t19 * t15) - t47 
     #* t24 - t48 * t26)
      t9 = t37 * (t3 * (t63 * (t23 + t13) * t19 * t21 * t5 ** 2 + t15 * 
     #(t46 * t35 + (t44 + t8) * t11 * t26) * t19 * t21 * t5) - t20 * t15
     # * t7 * t31 * t17 * t8 * (t26 + t24))
      ret = 8 * t8 * t3 * (t17 * (t20 * t7 ** 2 * t16 * t18 * t31 - t2 *
     # t37) - t59) + (0.5q1 / 0.8q1) * t29 * t41 * t38 * t19 * t21 * t17
     # * t62 + t28 / 4 + t1 / 8 + (0.3q1 / 0.8q1) * t64 * t52 * t2 * t38
     # * t19 * t62 - t37 * (t7 * ((t17 * (-t42 * t33 * t5 - t14) * t18 +
     # t44 + t17 * (t42 + t43) * t2) * t11 * t26 + t46 * (t17 * (-t45 * 
     #t33 * t5 - t14) * t18 + t34 + t17 * (t45 + t43) * t2)) * t19 * t15
     # + t51 * (t8 * (t32 + t30) + (t48 + t47) * t21 * t5 * t3) + t53 * 
     #(-t49 * t35 + t50 * (-t44 - t8)) * t19 * t41) / 2 + t58 * (t24 * (
     #t17 * (t45 * t5 * t54 * t21 + t56) * t19 * t22 + t57 + t21 * (t7 *
     # (-t17 * (t20 + t25) - t10) - t55) * t27) + t26 * (t17 * (t42 * t5
     # * t54 * t21 + t56) * t19 * t11 + t57 + t21 * (t7 * (-t17 * (t20 +
     # t12) - t10) - t55) * t6)) - 4 * t4 + 2 * t9

      hjetmass_box_ppmm_0_0_s34_mhsq_s12_s134 = ret/32q0/(0,1q0)
      return

      end function
