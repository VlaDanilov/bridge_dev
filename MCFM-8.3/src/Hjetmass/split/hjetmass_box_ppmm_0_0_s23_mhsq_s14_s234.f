
      complex*32 function hjetmass_box_ppmm_0_0_s23_mhsq_s14_s234
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
      t2 = zb(i4, i1)
      t3 = za(i1, i4)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = za(i2, i4)
      t7 = zb(i4, i2)
      t8 = zb(i4, i3)
      t9 = zb(i2, i1)
      t10 = zb(i3, i1)
      t11 = t1 * t10
      t12 = t6 * t9 + t11
      t13 = mt ** 2
      t14 = za(i1, i3)
      t15 = za(i1, i2) * t7
      t16 = t14 * t8
      t17 = t16 + t15
      t18 = t3 * t2
      t19 = t18 * t17
      t20 = t4 * t5
      t21 = t6 * t7
      t22 = t1 * t8
      t23 = t18 * (t22 + t20 + t21)
      t24 = 16 * t13 * t12 * t19 + 4 * t23 ** 2
      t24 = sqrt(t24)
      t23 = 2 * t23
      t25 = -t24 - t23
      t23 = t24 - t23
      t19 = 0.1q1 / t19
      t26 = t22 + t21
      t15 = t18 * t19 * (t15 * t26 + t16 * t26 + t20 * (t16 + t15))
      t16 = (0.1q1 / 0.2q1)
      t20 = 2 * t22 + 2 * t20 + 2 * t21
      t17 = t16 * t24 * t19 * t17
      t21 = -t15 + t20 - t17
      t12 = 0.1q1 / t12
      t22 = (0.1q1 / 0.4q1)
      t24 = t22 * t25 * t19 * t14
      t26 = t16 * t21 * t12 * t1
      t27 = -t24 * t7 + t26 * t9
      t24 = t8 * (t24 + t1) - t26 * t10
      t15 = -t15 + t20 + t17
      t17 = t14 * t23 * t19 * t22
      t11 = t8 * (t1 + t17) - t11 * t16 * t15 * t12
      t17 = t16 * t15 * t12 * t1 * t9 - t17 * t7
      t11 = 0.1q1 / t11
      t20 = 0.1q1 / t24
      t4 = 0.1q1 / t4
      t24 = 0.1q1 / t5
      t6 = 0.1q1 / t6
      t26 = t23 ** 2
      t28 = t23 * t26
      t29 = t25 ** 2
      t30 = t25 * t29
      t31 = t26 * t11
      t32 = t29 * t20
      t33 = t32 + t31
      t34 = t1 * t7
      t35 = -t34 + t27
      t36 = -t34 + t17
      t37 = t29 + t26
      t38 = t1 ** 2
      t39 = t7 * t10
      t40 = t37 * t4
      t41 = t19 * t7
      t42 = t6 * t19 ** 2
      t43 = t4 * t24
      t18 = t42 * t18 * t1
      t44 = 0.1q1 / t7
      t45 = t25 * t20
      t46 = t23 * t11
      t47 = t46 + t45
      t48 = -t46 - t45
      t49 = -t25 - t23
      t50 = t9 ** 2
      t25 = t27 * t25
      t51 = t45 * t27
      t46 = t46 * t17
      t52 = t9 * t3
      t7 = t19 * (t38 * (t44 * (t25 * t4 + (t11 * t5 + t4) * t17 * t23) 
     #+ t52 * t48 + t41 * t3 * t33) * t6 * t2 + t50 * t24 * t3 * (t34 * 
     #t47 - t51)) + t19 * (t2 * (t4 * t49 + t48 * t5) * t6 * t1 * t38 + 
     #t52 * (t24 * (t10 * (-t46 - t51) + t49 * t9) + ((t43 + t11) * t17 
     #* t23 + t25 * (t43 + t20)) * t44 * t2) * t6 * t1 + (t45 * t2 * t5 
     #* t27 * t44 + t9 * t47 * t13 + (t9 * t10 * t47 + t40 * t19 * t2) *
     # t24 * t7 * t3) * t6 * t38 - t46 * t50 * t3 * t24)
      t10 = 0.1q1 / t3
      t10 = t1 * t10
      t13 = t13 * t38 * t9
      t19 = t11 + t20
      ret = t16 * t43 * t18 * (t17 * t26 + t27 * t29) - t22 * t42 * (t1 
     #* (t9 * (t14 * t37 + (-t15 * t26 - t21 * t29) * t12 * t1) * t4 * t
     #2 + t31 * (t36 * t9 * t8 + t39 * (t34 - t17)) + t32 * (t35 * t9 * 
     #t8 + t39 * (t34 - t27))) * t24 * t3 + t41 * t2 * (t28 * t11 * t36 
     #+ t35 * t20 * t30) * t24 * t3 ** 2 + t14 * t2 * t38 * (t33 * t5 + 
     #t40)) - t18 * (-t12 * t1 * t9 * (t31 * t15 + t32 * t21) + t43 * t4
     #1 * (t30 + t28) * t14) / 8 + 8 * t13 * t4 * t6 * t44 * (t24 * t9 +
     # t10) + (0.3q1 / 0.4q1) * t18 * (t31 * t17 + t32 * t27) - t7 + 4 *
     # t13 * t6 * t44 * (t10 * t19 * t5 + t19 * t9)

      hjetmass_box_ppmm_0_0_s23_mhsq_s14_s234 = ret/32q0/(0,1q0)
      return

      end function
