
      complex*32 function hjetmass_box_pppm_0_0_s23_mhsq_s14_s234
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
      t2 = za(i2, i3)
      t3 = zb(i3, i2)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = zb(i2, i1)
      t10 = zb(i3, i1)
      t11 = t7 * t10
      t12 = t5 * t9 + t11
      t13 = mt ** 2
      t14 = za(i1, i2)
      t15 = za(i1, i3)
      t16 = t14 * t6
      t17 = t15 * t8
      t18 = t17 + t16
      t19 = t1 * t4
      t20 = t19 * t18
      t21 = t2 * t3
      t22 = t5 * t6
      t23 = t7 * t8
      t24 = t19 * (t23 + t21 + t22)
      t25 = 16 * t13 * t12 * t20 + 4 * t24 ** 2
      t25 = sqrt(t25)
      t24 = 2 * t24
      t26 = -t24 + t25
      t20 = 0.1q1 / t20
      t27 = t23 + t22
      t28 = t19 * t20
      t16 = t28 * (t16 * t27 + t21 * (t17 + t16) + t17 * t27)
      t17 = (0.1q1 / 0.2q1)
      t21 = 2 * t23 + 2 * t21 + 2 * t22
      t18 = t17 * t25 * t20 * t18
      t22 = t21 + t18 - t16
      t12 = 0.1q1 / t12
      t11 = t11 * t17
      t23 = t8 * (t15 * t26 * t20 / 4 + t7) - t11 * t22 * t12
      t24 = -t24 - t25
      t16 = t21 - t18 - t16
      t18 = t15 * t24
      t11 = t8 * (t18 * t20 / 4 + t7) - t11 * t16 * t12
      t11 = 0.1q1 / t11
      t2 = 0.1q1 / t2
      t21 = 0.1q1 / t23
      t23 = 0.1q1 / t14
      t25 = 0.1q1 / t15
      t27 = t1 ** 2
      t29 = t10 ** 2
      t30 = t26 ** 2
      t31 = t6 * t10
      t32 = t8 * t9
      t33 = t9 * (t19 + t13)
      t34 = t19 * t23
      t35 = t30 * t21
      t36 = (t26 + t24) * t2
      t37 = t36 * t3 * t4
      t38 = t24 * t11
      t39 = t26 * t21
      t40 = t3 * t11
      t41 = t20 * t1
      t42 = t4 ** 2
      t43 = t24 ** 2
      t44 = t8 * t23
      t45 = t3 * t7
      t46 = t6 * t29
      t47 = t13 * t9
      t48 = t45 * t4
      t49 = t36 * t25
      t5 = t41 * (t1 * (t39 * t23 * (t25 * (-t47 * t10 * t22 * t12 + t48
     # * t8) - t46) + t38 * (t25 * (t45 * t44 * t4 + (-t13 * t16 * t12 *
     # t23 - t6) * t10 * t9) - t46 * t23) + t4 * ((-t43 - t30) * t25 * t
     #6 - t44 * (t43 + t30)) * t20 * t2) + t9 * (t3 * (-t38 * t10 + t39 
     #* (t13 * t25 - t10)) - t49 * t14 * t9) + t29 * t23 * (t2 * (-t26 -
     # t24) - t39 * t3) * t15 + t38 * t27 * t3 * t42 * t23 * t25) + t41 
     #* (t25 * (-t35 * t27 * t4 * t6 * t8 * t23 * t20 + t38 * (t33 * t3 
     #+ t1 * (t13 * (t32 + t31) + t19 * (t8 * (-t6 * t24 * t20 + t9) - t
     #31)) * t23) + t37 * t7 + t39 * t1 * (t33 * t23 * t8 + t4 * (t34 + 
     #t9) * t3 + t31 * (t23 * (-t19 + t13) - t9))) - t37 * t23 * t5 - t4
     #0 * t18 * t29 * t23)
      t7 = t43 * t11
      t8 = t7 + t35
      t14 = t20 ** 2
      t15 = t44 * t1
      t18 = t10 * t23
      t24 = t32 * t8
      t29 = t9 ** 2
      t33 = t18 * t9
      t31 = (t32 - t31) * t9
      t32 = t3 * t23 * t42
      ret = t17 * t20 * t12 * t10 * t27 * (t7 * t28 * t16 * t6 * t23 * t
     #25 + t38 * (t25 * (-t48 * t23 + t29) + t33) * t16 + t39 * t22 * (t
     #25 * t29 + t33 + t23 * t25 * (t41 * t6 * t26 - t45) * t4)) - (0.3D
     #1 / 0.4q1) * t44 * t14 * t27 * t3 * t4 * (t7 + t35) + t25 * t14 * 
     #t12 * t27 * (t30 * (-t34 * t9 * t10 * t21 * t12 * t22 ** 2 + (-t32
     # * t1 + t31) * t21 * t22) + t7 * t16 * (t1 * (-t33 * t16 * t12 * t
     #4 - t32) + t31)) / 8 + t5 - t14 * t1 * (t6 * (-t3 * t10 * t8 + (t8
     # * t4 * t3 + t24) * t25 * t1) - t1 * t10 * t25 * t8 * t6 ** 2 + t1
     #9 * (t22 * t30 * ((-t15 * t21 - t2) * t25 * t9 - t18 * (t21 * t3 +
     # t2)) + ((-t15 * t11 - t2) * t25 * t9 - t18 * (t40 + t2)) * t16 * 
     #t43) * t12 + t24 * t3) / 4 + 2 * t41 * (t10 * (-t36 * t9 + (t3 * (
     #t39 + t38) + t36) * t23 * t13) + t49 * t47)

      hjetmass_box_pppm_0_0_s23_mhsq_s14_s234 = ret/32q0/(0,1q0)
      return

      end function
