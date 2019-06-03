
      complex*32 function hjetmass_box_ppmm_0_0_s14_mhsq_s23_s124
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
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = za(i2, i4)
      t5 = zb(i4, i1)
      t6 = zb(i4, i2)
      t7 = za(i1, i2)
      t8 = za(i2, i3)
      t9 = zb(i3, i2)
      t10 = t7 * t2
      t11 = t3 * t5
      t12 = t4 * t6
      t13 = t8 * t9
      t14 = t13 * (t12 + t10 + t11)
      t15 = za(i1, i3)
      t16 = t1 * t6
      t17 = t15 * t2
      t18 = t17 + t16
      t19 = t13 * t18
      t20 = zb(i3, i1)
      t21 = t20 * t7 + t4 * zb(i4, i3)
      t22 = mt ** 2
      t23 = 16 * t19 * t22 * t21 + 4 * t14 ** 2
      t23 = sqrt(t23)
      t14 = 2 * t14
      t24 = -t23 - t14
      t19 = 0.1q1 / t19
      t25 = t12 + t10
      t26 = t13 * t19
      t16 = t26 * (t16 * t25 + t11 * (t17 + t16) + t17 * t25)
      t25 = (0.1q1 / 0.2q1)
      t18 = t25 * t23 * t19 * t18
      t10 = 2 * t12 + 2 * t10 + 2 * t11
      t11 = -t18 + t10 - t16
      t12 = 0.1q1 / t21
      t21 = (0.1q1 / 0.4q1)
      t27 = t20 * t11 * t12 * t25
      t17 = t17 * t21
      t28 = t7 * (t2 - t27) + t17 * t24 * t19
      t14 = t23 - t14
      t10 = t18 + t10 - t16
      t16 = t20 * t10 * t12 * t25
      t17 = t7 * (t2 - t16) + t17 * t14 * t19
      t18 = t21 * t24 * t19 * t1 * t2 - t27 * t4
      t16 = t21 * t14 * t19 * t1 * t2 - t16 * t4
      t19 = 0.1q1 / t28
      t5 = 0.1q1 / t5
      t6 = 0.1q1 / t6
      t23 = 0.1q1 / t3
      t17 = 0.1q1 / t17
      t27 = t2 * t4
      t28 = t27 + t16
      t29 = t27 + t18
      t30 = -t27 - t16
      t31 = t10 ** 2
      t32 = t10 * t31
      t33 = t11 ** 2
      t34 = t11 * t33
      t35 = t33 + t31
      t36 = t12 ** 2
      t37 = t4 * t15
      t38 = t31 * t17
      t39 = t33 * t19
      t40 = t35 * t5
      t41 = t6 * t36
      t42 = t41 * t2
      t43 = t2 ** 2
      t26 = t41 * t26 * t1 * t43
      t44 = 0.1q1 / t4
      t45 = t10 * t17
      t46 = t11 * t19
      t47 = t46 + t45
      t48 = -t18 * t44 - t2
      t49 = t1 ** 2
      t50 = t9 ** 2
      t51 = t10 * t16
      t52 = t5 * t8
      t53 = t46 * t18
      t54 = t45 * t16
      t41 = t41 * t23 * t8 * t4
      t55 = t43 * t8 * t6
      t43 = t43 * t6
      t6 = t12 * (t49 * ((-t11 - t10) * t6 * t2 + t53 + t54) * t23 * t9 
     #+ t43 * (-t8 * t47 * t3 * t2 + t47 * t1 * t22 - t51 * t52 * t44) -
     # t41 * t34 * t18 * t50 * t19) + t12 * (t55 * t5 * (-t10 * t2 + t11
     # * t48) + t2 * (t4 * t23 * t47 * t49 + t6 * (t8 * (-t44 * (t11 * t
     #18 + t51) * t23 * t5 - t45 * (t16 * t44 + t2) + t46 * t48) + t23 *
     # (t45 * t28 + t46 * t29) * t15) * t1 - t52 * t6 * t23 * (t16 * t31
     # + t18 * t33) * t12) * t9 + t41 * (t30 * t17 * t32 - t27 * t34 * t
     #19) * t50 - t55 * t44 * (t54 + t53) * t3)
      t10 = 0.1q1 / t9
      t10 = t10 * t2
      ret = -t21 * t26 * (t38 * t14 + t39 * t24) - t25 * t26 * t23 * t5 
     #* (t14 * t31 + t24 * t33) + 8 * t43 * t22 * t1 * t44 * t5 * (t1 * 
     #t23 + t10) - 3 * t42 * t13 * (t38 * t16 + t39 * t18) - t42 * (t20 
     #* t8 * t2 * (t3 * (t39 + t38) + t40) + ((t1 * t35 + (-t34 - t32) *
     # t12 * t4) * t5 * t20 * t8 + t38 * (t30 * t7 * t1 + t37 * t28) + t
     #39 * (t37 * t29 + (-t27 - t18) * t7 * t1)) * t23 * t9) + 2 * t6 - 
     #4 * t43 * (t22 * (-t10 * t3 * (t19 + t17) * t44 * t1 + (-t19 - t17
     #) * t44 * t49) + t13 * t36 * t4 * (t40 * t23 + t38 + t39))

      hjetmass_box_ppmm_0_0_s14_mhsq_s23_s124 = ret/32q0/(0,1q0)
      return

      end function
