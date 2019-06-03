
      complex*32 function hjetmass_box_pmpm_0_0_0_s234_s23_s34
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
      t2 = zb(i4, i3)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = mt ** 2
      t8 = t1 * t2
      t9 = t3 * t4
      t10 = 0.1q1 / t4
      t11 = 0.1q1 / t3
      t12 = (0.1q1 / 0.2q1)
      t13 = t12 *sqrt(-t9 * t8 * (-2 * t8 * t3 * t4 - 8 * t7 * t5 * t6
     #)) * sqrt(0.2q1) * t11 * t10
      t14 = -t13 + t8
      t15 = zb(i3, i1)
      t16 = zb(i2, i1)
      t17 = za(i1, i4)
      t18 = za(i1, i2)
      t19 = za(i1, i3)
      t20 = t9 * t8
      t21 = t20 * (4 * t7 * t5 * t6 + t20)
      t21 = sqrt(t21)
      t22 = t20 + t21
      t2 = 0.1q1 / t2
      t6 = 0.1q1 / t6
      t1 = 0.1q1 / t1
      t23 = 0.1q1 / t5
      t24 = t22 * t11
      t25 = t12 * (t24 * t23 * t10 * t2 * t18 * t15 + t14 * t1 * t6 * t1
     #9 * t16)
      t26 = t18 * t16
      t27 = t25 + t26
      t28 = zb(i4, i1)
      t29 = t14 * t6 * t16
      t24 = t24 * t10 * t2 * t15 + t29
      t8 = t13 + t8
      t13 = t20 - t21
      t11 = t13 * t11
      t18 = t12 * (t11 * t23 * t10 * t2 * t18 * t15 + t8 * t1 * t6 * t19
     # * t16)
      t20 = t18 + t26
      t21 = t8 * t6 * t16
      t10 = t11 * t10 * t2 * t15 + t21
      t11 = t15 * t19 + t17 * t28
      t19 = t25 - t11
      t11 = t18 - t11
      t18 = 0.1q1 / t19
      t19 = 0.1q1 / t27
      t17 = 0.1q1 / t17
      t11 = 0.1q1 / t11
      t20 = 0.1q1 / t20
      t23 = 0.1q1 / t28
      t25 = -t11 + t20
      t26 = (-t18 + t19) * t24
      t27 = t8 * t10
      t28 = t26 * t14 + t27 * t25
      t30 = t14 ** 2
      t31 = t14 * t30
      t32 = t8 ** 2
      t33 = t8 * t32
      t34 = t32 * t10
      t35 = t30 * t24
      t36 = t35 * t19 + t34 * t20
      t37 = t6 ** 2
      t38 = t15 ** 2
      t39 = t5 ** 2
      t40 = t31 * t19
      t41 = t33 * t20
      t32 = t32 * t11
      t42 = t30 * t18
      t43 = t38 * t17
      t44 = t14 * (t18 - t19) + t8 * (t11 - t20)
      t45 = t8 * t20
      t46 = t14 * t19
      t21 = t21 + t10
      t47 = t46 * t24
      t14 = t14 * t18
      t48 = t16 * t1 * (t46 + t45)
      t49 = t23 * t15
      t9 = t49 * t6 * t17 * (t7 * (t5 * (t1 * (t27 * t20 + t47) + t15 * 
     #t44) + t48 * t39) + t9 * t6 * t15 * (t42 + t32))
      ret = -t12 * t49 * t1 ** 2 * t37 * t2 * t17 * (t34 * t13 * t11 + t
     #35 * t22 * t18) + t1 * t6 * (t17 * t39 * t4 * t28 + t23 * (t28 * t
     #38 + (t6 * (t36 * t1 * t16 * t5 - t15 * t36) + t16 * (t1 * (t33 * 
     #t10 * t25 + t26 * t31) + t15 * (-t41 - t40) + t1 * (t41 + t40) * t
     #16 * t5) * t37) * t17 * t4) * t3 + t43 * t6 * t23 * (t32 * t13 + t
     #42 * t22) * t2) + 2 * t6 * (t17 * (t39 * (t15 * (t11 * t8 + t14) +
     # t48 * t5) + t49 * t1 * (t5 * (-t30 * t16 * t19 * t6 - t45 * t21 -
     # t47) + t42 * t6 * (t29 + t24) + t32 * t6 * t21) * t3) * t4 + t43 
     #* t5 * t1 * t23 * (t45 * t13 + t46 * t22) * t2 + t49 * (t44 * t38 
     #* t3 + t7 * t1 * t5 * t17 * (t27 * t11 + t14 * t24))) - 4 * t9

      hjetmass_box_pmpm_0_0_0_s234_s23_s34 = ret/32q0/(0,1q0)
      return

      end function
