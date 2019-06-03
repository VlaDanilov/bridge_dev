
      complex*32 function hjetmass_qqbgg_bubble_pmpm_s34
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt

      t1 = za(i2, i4)
      t2 = za(i1, i4)
      t3 = zb(i3, i1)
      t4 = zb(i3, i2)
      t5 = t2 * t3
      t6 = t1 * t4
      t7 = t6 + t5
      t8 = zb(i4, i1)
      t9 = za(i1, i3)
      t10 = za(i2, i3)
      t11 = zb(i4, i2)
      t12 = t10 * t11 + t8 * t9
      t9 = t9 * t3
      t13 = t10 * t4
      t2 = t2 * t8
      t11 = t1 * t11
      t14 = t9 + t13 - t2 - t11
      t14 = t14 ** 2
      t15 = 4 * t7 * t12 + t14
      t15 = sqrt(t15)
      t16 = -t15 + t9 + t13 - t2 - t11
      t17 = za(i1, i2)
      t18 = zb(i2, i1)
      t14 = 4 * t7 * t12 + t14
      t14 = sqrt(t14)
      t19 = -t14 + t9 + t13 - t2 - t11
      t20 = 0.1q1 / t12
      t19 = 0.1q1 / t19
      t21 = 2
      t22 = (0.1q1 / 0.2q1)
      t23 = t22 * t16 * t20
      t24 = t21 * t7
      t25 = -t24 * t19 - t23
      t15 = t15 + t9 + t13 - t2 - t11
      t26 = t20 * (t16 - t15)
      t27 = -t23 * t8 + t3
      t28 = 0.1q1 / t12
      t29 = t17 * t18
      t30 = (-t29 - t2 - t11) * t1
      t31 = t10 ** 2 * t4
      t32 = t31 * t16
      t33 = (-t29 + t9) * t28
      t34 = t28 ** 2
      t35 = t28 * t34
      t36 = t16 ** 2
      t37 = t10 * t36 * t34
      t38 = t11 + t2
      t39 = t38 * t28
      t40 = (0.3q1 / 0.4q1) * t37 * t12
      t41 = (t13 + t9) * t28
      t42 = t1 * (t29 + t9)
      t43 = (t29 - t2 - t11) * t28
      t44 = -3 * t13 * t1
      t45 = t1 * t20
      t46 = t45 * (t29 + t2 + t11)
      t38 = t38 * t34
      t47 = t1 * t7
      t48 = t7 * t28
      t49 = -t10 * t16 * t36 * t35 * t12 / 8
      t50 = -t22 * t16 * (t38 * t16 * t10 + t46) - t37 * (t29 - t9 - t13
     #) / 4 + t48 * t16 * t10 + t47 + t49
      t51 = t29 - t2 - t11
      t52 = t45 * (t29 + t9)
      t53 = t9 * t10
      t4 = t1 ** 2 * t4
      t45 = (0.3q1 / 0.2q1) * t45 * t13
      t37 = t22 * t16 * (t53 * t16 * t34 + t32 * t34 + t52) + t37 * t51 
     #/ 4 + t5 * (t16 * t10 * t28 + t1) + t4 + t49 + t45 * t16
      t49 = -t22 * t8 * t15 * t20 + t3
      t54 = t15 ** 2
      t55 = t10 * t54 * t34
      t56 = -t10 * t15 * t54 * t35 * t12 / 8
      t29 = -t22 * t15 * (t38 * t15 * t10 + t46) + t48 * t15 * t10 + t47
     # + t55 * (-t29 + t9 + t13) / 4 + t56
      t2 = t14 + t9 + t13 - t2 - t11
      t2 = 0.1q1 / t2
      t9 = t22 * t15 * t20
      t11 = -t24 * t2 - t9
      t13 = t15 * t10
      t14 = t31 * t15
      t4 = t22 * t15 * (t53 * t15 * t34 + t14 * t34 + t52) + t5 * (t13 *
     # t28 + t1) + t4 + t55 * t51 / 4 + t56 + t45 * t15
      t12 = (0.3q1 / 0.4q1) * t55 * t12
      t20 = 0.1q1 / t25
      t17 = 0.1q1 / t17
      t11 = 0.1q1 / t11
      t22 = 0.1q1 / t26
      t18 = 0.1q1 / t18
      t24 = t29 - t4
      t4 = -t29 + t4
      t25 = t50 - t37
      t26 = t22 ** 2
      t2 = t2 * t11
      t12 = t2 * t15
      t29 = t18 * t35
      t31 = t1 * t8 + t10 * t3
      t35 = t16 * t19 * t20
      t18 = t18 * t34
      ret = 64 * t29 * t26 * t17 * t7 * (t16 * (t27 * (-t50 + t37) * t19
     # * t20 ** 2 + (t25 * t8 + t27 * (-t33 * t16 * t10 + t21 * t10 * (t
     #39 * t16 - t5 - t6) - t32 * t28 - t30 + t21 * t10 * (t41 * t16 + t
     #5) + t43 * t16 * t10 + t42 - t44)) * t19 * t20) + t12 * (t24 * t8 
     #+ t49 * (t4 * t11 - t13 * t33 - t14 * t28 + t21 * t10 * (t39 * t15
     # - t5 - t6) - t30 + t21 * t10 * (t41 * t15 + t5) - t44 + t13 * t43
     # + t42))) + 24 * t18 * t22 * t17 * t7 * (-t35 * (t1 * t27 + t3 * (
     #-t23 * t10 - t1)) + t12 * (t1 * t49 + t3 * (-t9 * t10 - t1))) + 25
     #6 * t29 * t22 * t26 * t17 * t7 * (t12 * t49 * t4 + t35 * t25 * t27
     #) - 128 * t18 * t26 * t17 * t7 * (t25 * t19 * t20 * t27 + t2 * t49
     # * t24) - 4 * t29 * t22 * t17 * t7 * (-t36 * t19 * t20 * t31 + t2 
     #* t31 * t54)

      hjetmass_qqbgg_bubble_pmpm_s34 = -ret/16q0*(0,1q0)
      return

      end function
