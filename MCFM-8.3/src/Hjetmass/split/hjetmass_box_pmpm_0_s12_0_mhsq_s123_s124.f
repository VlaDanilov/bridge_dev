
      complex*32 function hjetmass_box_pmpm_0_s12_0_mhsq_s123_s124
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
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = za(i3, i4)
      t5 = zb(i4, i1)
      t6 = zb(i4, i3)
      t7 = za(i2, i3)
      t8 = zb(i3, i2)
      t9 = za(i2, i4)
      t10 = zb(i4, i2)
      t11 = za(i1, i2)
      t12 = zb(i2, i1)
      t13 = t2 * t3
      t14 = t7 * t8
      t15 = t14 + t13
      t16 = t4 * t6
      t17 = t16 * t11 * t12
      t18 = t9 * t10
      t19 = t1 * t5
      t20 = t16 * (-t18 * t15 + t17 + t19 * (-t14 - t13))
      t21 = t1 * t2
      t22 = t8 * t9 + t21
      t23 = mt ** 2
      t24 = t3 * t5
      t10 = t7 * t10
      t25 = t24 + t10
      t26 = t6 ** 2
      t27 = t4 ** 2 * t26
      t28 = t27 * t25
      t29 = 16 * t16 * t23 * t22 * t28 + 4 * t20 ** 2
      t29 = sqrt(t29)
      t20 = -2 * t20
      t30 = t20 + t29
      t31 = t19 + t18
      t14 = t14 + t13
      t14 = t18 * t14 + t19 * t14 - t17
      t28 = 0.1q1 / t28
      t15 = t10 * t15 + t24 * t15
      t25 = t16 * t29 * t28 * t25 / 4
      t15 = t27 * t28 * (t18 * t15 + t19 * t15 + t17 * (-t10 - t24)) / 2
      t17 = t25 + t15 - t14
      t16 = t16 * t22
      t18 = 0.1q1 / t6
      t16 = 0.1q1 / t16
      t19 = 0.1q1 / t4
      t19 = t31 * t19 * t18 * t7
      t22 = t17 * t16 * t9
      t24 = t19 + t22
      t10 = t10 / 4
      t27 = -t10 * t30 * t28 + t8 * t24
      t32 = t30 * t28 / 4
      t24 = -t32 * t7 * t5 + t2 * t24
      t14 = -t25 + t15 - t14
      t15 = t20 - t29
      t20 = t14 * t16 * t9
      t19 = t20 + t19
      t10 = -t10 * t15 * t28 + t8 * t19
      t25 = t31 * t18 * t2
      t29 = t32 * t4 * t5 - t25
      t31 = t15 * t28 / 4
      t25 = t31 * t4 * t5 - t25
      t5 = -t31 * t7 * t5 + t2 * t19
      t8 = 0.1q1 / t8
      t19 = 0.1q1 / t10
      t11 = 0.1q1 / t11
      t31 = 0.1q1 / t27
      t32 = t9 * t12
      t33 = t15 ** 2
      t34 = t28 ** 2
      t35 = t30 ** 2
      t36 = t2 ** 2
      t37 = (t32 + t25) * t7
      t3 = (t1 * t7 + t3 * t9) * t11
      t38 = t4 * t5
      t39 = t11 * t9
      t32 = (t32 + t29) * t7
      t40 = t4 * t24
      t41 = t19 * t33
      t42 = 0.1q1 / t7
      t12 = 0.1q1 / t12
      t43 = t30 * t24
      t44 = t15 * t5
      t45 = t44 + t43
      t46 = t30 + t15
      t47 = t4 * t2
      t48 = -t23 * t31 + 1
      t49 = -t19 * t23 + 1
      t50 = t11 * t12
      t51 = t35 * t24
      t52 = t9 ** 2
      t53 = t30 * t29
      t54 = t15 * t25
      t55 = t42 * t45
      t56 = t2 * t7
      t57 = t30 * t31
      t58 = t19 * t15
      t59 = t11 * t8
      t13 = t28 * t2 * (t9 * (t58 * (-t56 + t5) + (t11 * (t6 * (t55 * t4
     # + t53 + t54) + t13 * t55) + t4 * t36 * t42 * (t14 * t15 + t17 * t
     #30) * t16) * t12 * t8 + t57 * (-t56 + t24)) + t59 * t6 * (t15 * t4
     #9 + t30 * t48) * t52 + t12 * t8 * t2 * (t47 * t46 + t53 + t54) + t
     #21 * t8 * (t50 * t45 - t56 * (t58 + t57)))
      ret = t28 * t34 * t11 * t8 * t7 * t26 * ((-t38 + t37) * t19 * t15 
     #* t33 + t31 * t30 * t35 * (-t40 + t32)) / 8 + t34 * t6 * (t39 * t7
     # * (t51 * t31 + t41 * t5) + (t4 * (t51 * (-t50 - t31) - t33 * t5 *
     # (t50 + t19)) + t11 * (t12 * (t25 * t33 + t29 * t35) + t9 * (t33 *
     # t49 + t35 * t48)) * t7) * t8 * t2) / 2 - 8 * t12 * t36 * (t1 * t3
     #6 * t8 * t18 + t59 * t23 * t52 * t42 + t9 * t2 * t18) + t34 * t8 *
     # t6 * (t41 * (t2 * (t3 * t5 + t37) + t4 * (t20 + t7) * t36 + t39 *
     # (t38 + t37) * t6) + (t2 * (t3 * t24 + t32) + t4 * (t22 + t7) * t3
     #6 + t39 * (t40 + t32) * t6) * t31 * t35) / 4 - 2 * t12 * t28 * t2 
     #* (t47 * t42 * t8 * t45 + t39 * (t46 * t8 * t2 * t23 - t43 - t44))
     # + t13 + 4 * t18 * t12 * t42 * t9 * t36 * (-t2 * t8 * (t27 + t10) 
     #+ t24 + t5)

      hjetmass_box_pmpm_0_s12_0_mhsq_s123_s124 = ret/32q0/(0,1q0)
      return

      end function
