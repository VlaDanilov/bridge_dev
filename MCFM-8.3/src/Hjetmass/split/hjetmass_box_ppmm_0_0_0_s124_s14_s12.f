
      complex*32 function hjetmass_box_ppmm_0_0_0_s124_s14_s12
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
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = mt ** 2
      t8 = t1 * t2
      t9 = t3 * t4
      t10 = 0.1q1 / t4
      t11 = 0.1q1 / t3
      t12 = sqrt(-t9 *t8 * (-2 * t8 * t3 * t4 - 8 * t7 * t5 * t6)) * 
     #sqrt(0.2q1) * t11 * t10 / 2
      t13 = -t12 + t8
      t14 = za(i3, i4)
      t15 = zb(i3, i2)
      t16 = t9 * t8
      t17 = t16 * (4 * t7 * t5 * t6 + t16)
      t17 = sqrt(t17)
      t18 = t17 + t16
      t19 = zb(i3, i1)
      t20 = za(i1, i3)
      t21 = zb(i4, i3)
      t6 = 0.1q1 / t6
      t22 = 0.1q1 / t2
      t1 = 0.1q1 / t1
      t5 = 0.1q1 / t5
      t23 = t13 * t1
      t24 = t18 * t11 * t5
      t25 = t24 * t22 * t10 * t14 * t19 / 2 + t23 * t6 * t20 * t21 / 2
      t26 = t14 * t21
      t27 = t25 + t26
      t23 = t24 * t10 * t14 + t23 * t20
      t24 = za(i2, i3)
      t28 = t15 * t24 + t19 * t20
      t25 = -t28 + t25
      t8 = t12 + t8
      t12 = t17 - t16
      t16 = t8 * t1 * t20
      t10 = t12 * t11 * t5 * t10 * t14
      t11 = -t16 + t10
      t10 = t10 * t22 * t19 / 2 - t16 * t6 * t21 / 2
      t16 = -t28 - t10
      t10 = -t10 + t26
      t10 = 0.1q1 / t10
      t17 = 0.1q1 / t24
      t16 = 0.1q1 / t16
      t19 = 0.1q1 / t27
      t21 = 0.1q1 / t25
      t15 = 0.1q1 / t15
      t22 = t16 - t10
      t24 = t2 ** 2
      t25 = t1 ** 2
      t26 = t13 ** 2
      t27 = t13 * t26
      t28 = t14 ** 2
      t29 = t20 ** 2
      t30 = t8 ** 2
      t31 = t8 * t30
      t32 = t30 * t16
      t33 = t32 * t11
      t34 = t26 * t21
      t35 = t34 * t23
      t26 = t26 * t19
      t30 = t30 * t10
      t36 = t8 * t10
      t37 = t36 * t11
      t38 = t13 * t19
      t39 = t38 * t23
      t40 = (-t19 + t21) * t13
      t41 = t8 * t16
      t42 = t13 * t21
      t43 = t1 * t2
      t44 = (-t16 + t10) * t8
      t45 = t44 * t11 + t40 * t23
      t10 = t10 * t31 + t19 * t27
      t19 = t19 - t21
      t46 = t2 * t20
      t47 = t46 * t6
      t48 = t6 ** 2
      t9 = t43 * t15 * t17 * (t7 * (t28 * (-t36 - t38) + t6 * (t46 * (t1
     #3 * t19 + t44) + t37 - t39) * t14) - t43 * t9 * t48 * t29 * (t34 +
     # t32))
      ret = 4 * t9 - 2 * t43 * (t6 * (t4 * ((t14 * (t39 - t37) - t6 * (t
     #16 * t31 + t21 * t27) * t25 * t29 + (t14 * (t30 + t26) + t6 * (-t3
     #5 + t33)) * t1 * t20) * t15 * t3 + t20 * t28 * (t22 * t8 + t40)) +
     # t28 * t15 * (t36 * t12 - t38 * t18) * t5 + t7 * t14 * t15 * (t41 
     #* t11 - t42 * t23)) * t17 + t3 * t14 * t15 * (t36 + t38) * t2 + t3
     # * t20 * t6 * t15 * (t42 + t41) * t24) + t1 * t6 * (t17 * (t28 * t
     #45 + t3 * ((t30 * t11 * (t47 - t14) + t26 * (-t47 + t14) * t23) * 
     #t15 * t1 + t20 * (t10 * t14 + t6 * (t31 * t11 * t22 + t19 * t23 * 
     #t27) - t47 * t10) * t15 * t25)) * t4 + t43 * t14 * t20 * t6 * t17 
     #* t15 * (-t32 * t12 + t34 * t18) * t5 + t15 * t3 * t24 * t45) - t2
     #5 * t48 * t5 * t14 * t17 * t15 * (t33 * t12 + t35 * t18) / 2

      hjetmass_box_ppmm_0_0_0_s124_s14_s12 = ret/32q0/(0,1q0)
      return

      end function
