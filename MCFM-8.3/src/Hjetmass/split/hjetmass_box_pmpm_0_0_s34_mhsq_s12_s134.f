
      complex*32 function hjetmass_box_pmpm_0_0_s34_mhsq_s12_s134
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
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i1, i2)
      t6 = za(i2, i3)
      t7 = zb(i2, i1)
      t8 = za(i2, i4)
      t9 = t6 * t2
      t10 = t8 * t4
      t11 = t10 + t9
      t12 = t5 * t7
      t13 = t12 * t11
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = zb(i3, i2)
      t17 = zb(i4, i2)
      t18 = t3 * t17
      t19 = t1 * t16 + t18
      t20 = mt ** 2
      t1 = t1 * t2
      t21 = t3 * t4
      t22 = t14 * t15
      t23 = t12 * (t22 + t1 + t21)
      t24 = 16 * t13 * t20 * t19 + 4 * t23 ** 2
      t24 = sqrt(t24)
      t13 = 0.1q1 / t13
      t9 = t10 + t9
      t25 = t12 * t13
      t9 = t25 * (t1 * t9 + t21 * t9 + t22 * t9)
      t26 = 2
      t27 = (0.1q1 / 0.2q1)
      t11 = t27 * t24 * t13 * t11
      t1 = t26 * (t22 + t1 + t21)
      t22 = t1 - t11 - t9
      t23 = t26 * t23
      t28 = -t23 - t24
      t19 = 0.1q1 / t19
      t29 = (0.1q1 / 0.4q1)
      t30 = t29 * t28 * t13 * t8
      t31 = t27 * t22 * t19 * t3
      t32 = t31 * t16 - t30 * t2
      t30 = t4 * (t3 + t30) - t31 * t17
      t1 = t1 + t11 - t9
      t9 = -t23 + t24
      t11 = t8 * t9 * t13 * t29
      t18 = t4 * (t3 + t11) - t18 * t27 * t1 * t19
      t11 = t27 * t1 * t19 * t3 * t16 - t11 * t2
      t6 = 0.1q1 / t6
      t23 = 0.1q1 / t30
      t14 = 0.1q1 / t14
      t24 = 0.1q1 / t15
      t18 = 0.1q1 / t18
      t30 = t3 * t7
      t31 = t22 * t23
      t33 = t1 * t18
      t34 = t33 + t31
      t35 = t2 * t3
      t36 = -t35 + t11
      t37 = -t35 + t32
      t38 = t8 ** 2
      t39 = t22 ** 2
      t40 = t22 * t39
      t41 = t19 ** 2
      t42 = t19 * t41
      t43 = t1 ** 2
      t44 = t1 * t43
      t45 = t2 ** 2
      t46 = t19 * (t20 * t16 * t24 - t30)
      t47 = t43 * t18
      t30 = t30 * t16 * t24
      t48 = (t22 + t1) * t14
      t49 = t24 * t5
      t50 = t2 * t24
      t51 = t44 * t18
      t52 = t39 * t23
      t53 = t41 * t6
      t54 = t5 ** 2
      t55 = t12 * t8
      t56 = t19 * t6 * t2 * (t3 * (t50 * (t52 + t47) * t19 * t7 * t54 + 
     #t55 * t14 * t24 * (t43 + t39) * t19) + t38 * t20 * (t48 * t24 + t3
     #1 + t33))
      t57 = t23 * t40 + t51
      t13 = t53 * t12 * t24 * ((t18 * t43 ** 2 + t23 * t39 ** 2) * t41 *
     # t16 ** 2 * t5 * t3 + t13 * t38 * t2 * t14 * (t28 * t39 + t43 * t9
     #))
      t39 = t53 * t50 * t54 * t7 * (t47 * t11 + t52 * t32)
      ret = t26 * t19 * (t12 * t2 * t34 * t6 * t38 + t7 * (t15 * t34 + t
     #48) * t6 * t8 * t38 + t49 * t45 * (t31 * t37 + t33 * t36) + t5 * (
     #t2 * (t47 * t46 + t31 * (t4 * t32 * t24 + t46 * t22) + t33 * t4 * 
     #t11 * t24) + t30 * t41 * (t44 + t40) * t14 + t24 * (t1 * (-t18 * t
     #20 + 1) + t22 * (-t20 * t23 + 1) + t21 * (-t33 - t31)) * t45) * t6
     # * t8) - t27 * t13 - t29 * t53 * t25 * t38 * t2 * (t52 * t28 + t47
     # * t9) + (0.5q1 / 0.2q1) * t30 * t42 * t54 * t2 * t6 * t57 + 3 * t
     #39 + (0.3q1 / 0.2q1) * t55 * t42 * t3 * t16 * t6 * t57 - t53 * t5 
     #* (t51 * t12 * t19 * t16 * t11 * t24 + t47 * t8 * (t16 * (t4 * t24
     # * t36 + t7 * t8) + t50 * (t35 - t11) * t17) + t52 * (t16 * (t7 * 
     #(t49 * t22 * t32 * t19 + t38) + t10 * t24 * t37) + t50 * t8 * (t35
     # - t32) * t17)) - 4 * t56

      hjetmass_box_pmpm_0_0_s34_mhsq_s12_s134 = ret/32q0/(0,1q0)
      return

      end function
