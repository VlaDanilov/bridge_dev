
      complex*32 function hjetmass_box_pmpm_0_0_0_s134_s34_s14
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
      t2 = zb(i3, i1)
      t3 = za(i1, i4)
      t4 = za(i3, i4)
      t5 = zb(i4, i1)
      t6 = zb(i4, i3)
      t7 = za(i1, i3)
      t8 = mt ** 2
      t9 = t3 * t5
      t10 = t9 * t4 * t6
      t11 = t10 * (4 * t8 * t7 * t2 + t10)
      t11 = sqrt(t11)
      t12 = t10 + t11
      t13 = zb(i3, i2)
      t14 = 0.1q1 / t4
      t15 = 0.1q1 / t6
      t4 = sqrt(-t4 * t6* t9 * (-8 * t8 * t7 * t2 - 2 * t10)) * sqrt(
     #0.2q1) * t14 * t15 / 2
      t6 = -t4 + t9
      t16 = za(i2, i4)
      t17 = za(i2, i3)
      t18 = zb(i4, i2)
      t19 = t1 * zb(i2, i1) + t16 * t18
      t3 = 0.1q1 / t3
      t5 = 0.1q1 / t5
      t20 = 0.1q1 / t2
      t7 = 0.1q1 / t7
      t21 = t6 * t3
      t22 = t12 * t7 * t14
      t23 = t22 * t5 * t15
      t24 = t21 * t20 * t16 * t13 / 2 + t23 * t17 * t18 / 2
      t25 = t24 - t19
      t10 = t10 - t11
      t4 = t4 + t9
      t9 = t4 * t3
      t11 = t10 * t7 * t14
      t18 = t11 * t5 * t15 * t17 * t18 / 2 + t9 * t20 * t16 * t13 / 2
      t19 = t18 - t19
      t20 = t17 * t13
      t24 = t24 + t20
      t18 = t18 + t20
      t20 = t22 * t15 * t17 + t21 * t16
      t9 = t11 * t15 * t17 + t9 * t16
      t11 = 0.1q1 / t25
      t19 = 0.1q1 / t19
      t21 = 0.1q1 / t17
      t22 = 0.1q1 / t24
      t13 = 0.1q1 / t13
      t18 = 0.1q1 / t18
      t24 = t10 ** 2
      t25 = t10 * t24
      t26 = t12 ** 2
      t27 = t12 * t26
      t28 = t26 * t11
      t29 = t24 * t19
      t30 = t28 * t20 + t29 * t9
      t31 = t29 + t28
      t32 = t5 ** 2
      t33 = t16 ** 2
      t34 = t7 * t5
      t35 = t26 * t22
      t36 = t34 * t16
      t37 = t33 * t13
      t38 = t16 * t2 * t21 * t13
      t10 = t10 * t19
      t12 = t12 * t11
      t39 = (-t18 + t19) * t24
      t40 = t34 * t13
      t41 = t14 * t3
      t42 = t2 * t1 * t15
      t43 = t7 ** 2
      t43 = t7 * t43 * t15 * t14 * t5 * t32 * t3 ** 2 * t16 * t13 * (t25
     # * t4 * t18 + t27 * t6 * t22)
      ret = t34 * (-t37 * t2 * t3 * t21 * (t10 * t4 + t12 * t6) + ((t21 
     #* t30 * t3 * t1 - t14 * t30) * t13 * t5 * t2 + t41 * t33 * (t39 + 
     #t28)) * t7 * t15) + t34 * (t3 * (-t38 * t15 * t7 * t31 * t1 + t15 
     #* (t7 * (-t35 * t33 + (t25 * t9 * (t18 - t19) + (-t11 + t22) * t20
     # * t27) * t7 * t13 * t32 + t36 * (t11 * t27 + t19 * t25) * t13) + 
     #t37 * t31 * t21 * t2) * t14 + t36 * t21 * t13 * t30) + t38 * (t10 
     #* t9 + t12 * t20) + t40 * t15 * t14 * t2 ** 2 * (t26 * (t11 - t22)
     # + t39) * t17) + 2 * t40 * t2 * (t2 * (-t16 * (t10 + t12) + (-t10 
     #- t12) * t15 * t2 * t1) + t41 * t36 * t15 * (t18 * t24 + t35) * t8
     #) - 4 * t36 * t13 * t3 * t2 * t8 * (-t10 * t21 * (t42 + t16) + t12
     # * (t21 * (-t42 - t16) + t23) + t34 * t29 * t15 * t14) - t43 / 2

      hjetmass_box_pmpm_0_0_0_s134_s34_s14 = ret/32q0/(0,1q0)
      return

      end function
