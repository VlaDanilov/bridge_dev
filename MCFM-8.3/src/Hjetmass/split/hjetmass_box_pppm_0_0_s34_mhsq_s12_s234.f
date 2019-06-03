
      complex*32 function hjetmass_box_pppm_0_0_s34_mhsq_s12_s234
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt

      t1 = za(i2, i3)
      t2 = zb(i3, i2)
      t3 = za(i2, i4)
      t4 = zb(i4, i2)
      t5 = za(i1, i3)
      t6 = zb(i2, i1)
      t7 = za(i1, i4)
      t8 = t5 * t2
      t9 = t7 * t4
      t10 = t9 + t8
      t11 = za(i1, i2) * t6
      t12 = t11 * t10
      t13 = za(i3, i4)
      t14 = zb(i4, i3)
      t15 = zb(i3, i1)
      t16 = t1 * t15
      t17 = t3 * zb(i4, i1) + t16
      t18 = t1 * t2
      t4 = t3 * t4
      t19 = t13 * t14
      t20 = t11 * (t19 + t18 + t4)
      t21 = 16 * mt ** 2 * t17 * t12 + 4 * t20 ** 2
      t21 = sqrt(t21)
      t12 = 0.1q1 / t12
      t22 = t4 + t18
      t8 = t11 * t12 * (t8 * t22 + t19 * (t9 + t8) + t9 * t22)
      t9 = 2
      t11 = (0.1q1 / 0.2q1)
      t4 = t9 * (t19 + t18 + t4)
      t10 = t11 * t21 * t12 * t10
      t18 = -t10 + t4 - t8
      t19 = t9 * t20
      t20 = -t19 - t21
      t17 = 0.1q1 / t17
      t22 = (0.1q1 / 0.4q1)
      t16 = t16 * t11
      t23 = t2 * (t5 * t20 * t12 * t22 + t1) - t16 * t18 * t17
      t24 = t11 * t18 * t17 * t3 * t15 - t20 * t22 * t12 * t7 * t2
      t19 = -t19 + t21
      t4 = t10 + t4 - t8
      t8 = t22 * t19 * t12
      t10 = t11 * t4 * t17 * t3 * t15 - t8 * t7 * t2
      t1 = t2 * (t8 * t5 + t1) - t16 * t4 * t17
      t1 = 0.1q1 / t1
      t5 = 0.1q1 / t23
      t8 = 0.1q1 / t14
      t14 = 0.1q1 / t7
      t13 = 0.1q1 / t13
      t16 = t3 ** 2
      t21 = t10 ** 2
      t23 = t2 ** 2
      t25 = t23 * t16
      t26 = t25 + t21
      t27 = t24 ** 2
      t25 = (t25 + t27) * t5
      t28 = t19 * t4 * t1
      t4 = t4 * t1
      t29 = t6 * t13
      t30 = t14 * t15
      t31 = t17 * t8
      t32 = 0.1q1 / t3
      t33 = t24 * t32 - t2
      t34 = t10 * t32 - t2
      t35 = t6 ** 2
      t36 = t24 * t33
      t37 = t10 * t34
      t38 = t5 * t20
      t39 = t1 + t5
      t40 = t38 * t18
      t41 = t5 * t24
      t42 = t10 * t1
      t43 = t24 * t14
      t44 = t2 * t3
      t45 = t8 * t6
      t46 = t29 * t3
      ret = -t11 * t31 * t13 * t12 * t15 * t7 * t35 * (t28 * (-t44 + t10
     #) + t40 * (-t44 + t24)) + t22 * t13 * t8 * t12 ** 2 * t7 ** 2 * t3
     #5 * t2 * (t19 ** 2 * t1 * t34 + t33 * t5 * t20 ** 2) - t9 * t31 * 
     #t6 * (t30 * (t25 * t18 + t4 * t26) + t29 * (t25 * t20 * t18 + t28 
     #* t26) * t12) - 8 * t45 * t2 * (t21 * t14 * t1 + t42 * (t46 + t15)
     # + t41 * (t46 + t15 + t43)) + t8 * t12 * t6 * (t38 * (t2 * (t33 * 
     #t15 * t7 + t36) + t36 * t29 * t7) + (t2 * (t34 * t15 * t7 + t37) +
     # t37 * t29 * t7) * t1 * t19) + 4 * t45 * (t23 * (t3 * (t14 * (t42 
     #+ t41) + t15 * t39) + t29 * t39 * t16) + t10 * t21 * t14 * t32 * t
     #1 + t44 * t17 * (t30 * (t4 * t10 + t41 * t18) + t29 * (t28 * t10 +
     # t40 * t24) * t12) + (t15 * t32 + t29) * t1 * t21 + t5 * t27 * (t3
     #2 * (t43 + t15) + t29))

      hjetmass_box_pppm_0_0_s34_mhsq_s12_s234 = 
     & ret/32q0/(0,1q0)
      return

      end function
