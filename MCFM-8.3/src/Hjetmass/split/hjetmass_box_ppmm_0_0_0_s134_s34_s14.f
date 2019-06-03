
      complex*32 function hjetmass_box_ppmm_0_0_0_s134_s34_s14
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i3, i4)
      t2 = za(i1, i4)
      t3 = zb(i4, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = mt ** 2
      t8 = t2 * t3
      t9 = t8 * t1 * t4
      t10 = t9 * (4 * t7 * t5 * t6 + t9)
      t10 = sqrt(t10)
      t11 = t10 - t9
      t12 = zb(i4, i2)
      t13 = za(i2, i3)
      t14 = zb(i2, i1)
      t15 = 0.1q1 / t4
      t16 = 0.1q1 / t1
      t4 = sqrt(-t1* t4 * t8 * (-8 * t7 * t5 * t6 - 2 * t9)) * sqrt(0
     #.2q1) * t16 * t15 / 2
      t17 = t8 + t4
      t18 = za(i2, i4)
      t19 = zb(i3, i2)
      t6 = 0.1q1 / t6
      t2 = 0.1q1 / t2
      t3 = 0.1q1 / t3
      t5 = 0.1q1 / t5
      t20 = t6 * t19
      t21 = t20 * t17
      t22 = t21 * t2
      t23 = t3 * t12
      t24 = t23 * t11
      t25 = t12 * t18 + za(i1, i2) * t14
      t26 = -t24 * t5 * t16 * t15 * t13 / 2 + t22 * t18 / 2
      t27 = -t25 + t26
      t9 = t10 + t9
      t4 = t8 - t4
      t8 = t20 * t4
      t10 = t23 * t9
      t18 = t10 * t5 * t16 * t15 * t13 / 2 + t8 * t2 * t18 / 2
      t25 = -t25 + t18
      t28 = t13 * t19
      t18 = t28 + t18
      t26 = t28 + t26
      t28 = t20 - t23
      t29 = 0.1q1 / t19
      t13 = 0.1q1 / t13
      t18 = 0.1q1 / t18
      t25 = 0.1q1 / t25
      t26 = 0.1q1 / t26
      t27 = 0.1q1 / t27
      t30 = t25 - t18
      t31 = t17 ** 2
      t32 = t15 ** 2
      t33 = t2 ** 2
      t34 = t14 ** 2
      t35 = t11 ** 2
      t36 = t4 ** 2
      t37 = t9 ** 2
      t38 = t1 ** 2
      t39 = t23 * t29
      t40 = t14 * t29
      t41 = t40 * t27
      t42 = t12 * t29
      t43 = t17 * t6
      t44 = t13 * t1
      t45 = t5 * t2
      t46 = t6 * t3
      t47 = t40 * t25
      t48 = t4 * t6
      t49 = t9 * t25
      t23 = t23 * t32 * t13
      t50 = t45 * t1
      t51 = t3 * t34
      t52 = t9 * t4
      t53 = t11 * t17
      t54 = t46 * t45
      t28 = t46 * t5 * t15 * t33 * t13 * t38 * (-t53 * t11 * t5 * t15 * 
     #t17 * t2 * t28 * t26 - t52 * t9 * t5 * t15 * t4 * t2 * t28 * t18 +
     # t54 * (t37 * t36 * t18 + t35 * t31 * t26) * t14)
      t55 = t11 * t27
      t56 = t49 - t55
      t56 = t44 * t15 * t3 * t2 * t14 * t7 * (t40 * t56 + (t42 * t56 + t
     #6 * t2 * (t53 * t27 - t49 * t4)) * t5 * t1)
      t55 = -t49 + t55
      t7 = t15 * t2 * t13 * t38 * (t55 * t5 * t12 * t1 + t14 * t55 + t54
     # * t14 * (-t52 * t18 + t53 * t26) * t7)
      ret = -t28 / 4 - t50 * (t15 * (t48 * (t30 * t2 * t13 * t19 * t38 +
     # t51 * t30) - t44 * t47 * (t10 * t16 * t15 + t8)) * t9 + t44 * (-t
     #51 * t43 * t29 * t27 + t41 * (-t24 * t16 * t15 + t21) * t15 + t22 
     #* t1 * t26 * t15) * t11 - t47 * t4 * t2 * t37 * t32 * t3 ** 2 * t1
     #2 * t13) - t50 * (t23 * (t1 * t6 * t3 * t5 * (-t26 + t27) * t33 * 
     #t31 + t27 * (-t40 * t3 + (-t39 - t6) * t5 * t1) * t2 * t17 + t41) 
     #* t35 + t49 * t3 * t13 * t4 * t14 * t1 * (t15 * (t48 * t2 - t42) +
     # t40 * t6) + t15 * t17 * (-t20 * t38 * t13 * t2 * t27 + t6 * (t26 
     #- t27) * t3 * t34 + t44 * t27 * (-t43 * t2 + t42) * t3 * t14) * t1
     #1 + t23 * (t1 * (t46 * t30 * t5 * t33 * t36 - t45 * t25 * (t39 + t
     #6) * t4) + t47) * t37) + 4 * t56 - 2 * t7

      hjetmass_box_ppmm_0_0_0_s134_s34_s14 = ret/32q0/(0,1q0)
      return

      end function
