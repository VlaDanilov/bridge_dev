
      complex*32 function hjetmass_box_pppm_0_0_0_s134_s34_s14
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
      t2 = za(i3, i4)
      t3 = zb(i4, i1)
      t4 = zb(i4, i3)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = mt ** 2
      t8 = t1 * t3
      t9 = t8 * t2 * t4
      t10 = t9 * (4 * t7 * t5 * t6 + t9)
      t10 = sqrt(t10)
      t11 = t10 - t9
      t12 = zb(i2, i1)
      t13 = za(i2, i4)
      t14 = 0.1q1 / t2
      t15 = 0.1q1 / t4
      t4 = sqrt(-t2 * t4 * t8 * (-8 * t7 * t5 * t6 - 2 * t9))*sqrt(0.
     #2q1) * t14 * t15 / 2
      t7 = zb(i3, i2)
      t16 = za(i2, i3)
      t17 = zb(i4, i2)
      t3 = 0.1q1 / t3
      t5 = 0.1q1 / t5
      t18 = 0.1q1 / t1
      t19 = 0.1q1 / t6
      t20 = t7 * t18
      t21 = t12 * za(i1, i2) + t13 * t17
      t22 = t11 * t5 * t14 * t3 * t15 * t16 * t17 / 2 - t20 * (t8 + t4) 
     #* t19 * t13 / 2
      t23 = -t22 - t21
      t9 = t10 + t9
      t4 = t9 * t5 * t14 * t3 * t15 * t16 * t17 / 2 + t20 * (t8 - t4) * 
     #t19 * t13 / 2
      t8 = -t21 + t4
      t10 = t16 * t7
      t4 = t10 + t4
      t10 = -t22 + t10
      t4 = 0.1q1 / t4
      t13 = 0.1q1 / t13
      t16 = 0.1q1 / t23
      t8 = 0.1q1 / t8
      t10 = 0.1q1 / t10
      t19 = t8 - t4
      t21 = t11 ** 2
      t22 = t11 * t21
      t23 = t9 ** 2
      t24 = t9 * t23
      t25 = (t16 - t10) * t21 + t19 * t23
      t26 = -t16 + t10
      t27 = t4 - t8
      t28 = t21 * t26 + t23 * t27
      t29 = t3 ** 2
      t30 = t15 ** 2
      t31 = t5 ** 2
      t32 = t18 * t3
      t33 = t12 * t14
      t34 = t15 * t13 * t3
      t35 = t7 * t12
      t36 = t15 * t3
      t37 = t6 * t13
      t38 = t12 * t13
      t39 = t6 * t17 ** 2
      t40 = t39 * t14 * t30
      t41 = t15 * (t14 * t17 * t15 + t13) * t7
      t42 = t7 * t13
      t39 = t39 * t29 * t18
      t43 = t3 * (t32 * t17 + t13) * t12
      t44 = t15 * t14
      ret = t36 * t31 * (t34 * t17 * t25 * t6 ** 2 + t18 * t30 * t13 * t
     #29 * t14 * t17 * (t21 ** 2 * t26 + t23 ** 2 * t27) * t31 + t35 * (
     #t28 * t15 * t14 + t32 * t28) + t36 * t17 * (t20 * t25 + t33 * t25)
     # * t6) + 2 * t36 * t5 * (t6 * (t35 * (t11 * t26 + t19 * t9) + t37 
     #* (t2 * t7 * (t10 * t11 - t4 * t9) + (-t11 * t16 + t8 * t9) * t12 
     #* t1)) + t34 * (t20 * (-t10 * t22 + t24 * t4) + t33 * (t16 * t22 -
     # t24 * t8)) * t31) - t36 * t5 * t31 * (t22 * (t32 * (t3 * ((-t33 +
     # t37) * t15 * t17 - t38 + t40) - t41) * t16 + t44 * (t15 * ((t20 -
     # t37) * t3 * t17 + t42 - t39) + t43) * t10) + t24 * (t44 * (t15 * 
     #((-t20 + t37) * t3 * t17 - t42 + t39) - t43) * t4 + t32 * (t3 * ((
     #t33 - t37) * t15 * t17 + t38 - t40) + t41) * t8)) / 2

      hjetmass_box_pppm_0_0_0_s134_s34_s14 = ret/32q0*(0,1q0)
      return

      end function
