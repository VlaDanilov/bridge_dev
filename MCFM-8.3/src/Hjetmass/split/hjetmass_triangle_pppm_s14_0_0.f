
      complex*32 function hjetmass_triangle_pppm_s14_0_0
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
      t2 = za(i1, i2)
      t3 = zb(i4, i2)
      t4 = za(i1, i3)
      t5 = za(i3, i4)
      t6 = zb(i2, i1)
      t7 = zb(i4, i1)
      t8 = zb(i4, i3)
      t9 = zb(i3, i1)
      t10 = t9 * t3
      t11 = t6 * t8
      t12 = t10 - t11
      t13 = za(i2, i3)
      t14 = za(i2, i4)
      t15 = t14 * t3
      t16 = t5 * t8
      t17 = t15 + t16
      t18 = zb(i3, i2)
      t19 = t18 * t7
      t20 = t10 - t11 - t19
      t21 = 3 * t19 + t12
      t22 = t1 * t7
      t23 = t22 * t20
      t24 = t13 * t3
      t25 = 4
      t26 = t4 ** 2
      t27 = t4 * t26
      t28 = t16 * t25
      t29 = t13 * t18
      t30 = t15 + t29
      t31 = t19 * t25
      t32 = t15 * t20
      t33 = t30 + 3 * t16
      t34 = t11 + t19
      t35 = t3 ** 2
      t36 = t14 ** 2
      t37 = t8 ** 2
      t38 = t5 ** 2 * t37
      t39 = t4 * t8
      t40 = t4 * t13
      t41 = t16 * t20
      t42 = t15 * t31 - t41
      t43 = 7 * t19
      t44 = t17 ** 2
      t45 = t2 ** 2
      t46 = t2 * t3
      t32 = t45 * t35 * (t25 * t7 * t8 * (t1 * t13 * t5 * t35 * t30 + t2
     #7 * (t14 * t6 - t5 * t9) * t12) + t40 * (t3 * (t34 * t36 * t35 - t
     #36 * t9 * t3 * t35 + t29 * (t28 * t19 - t32) - t23 * t33 + 3 * t38
     # * t12 + t16 * t19 * (t15 * t25 + 9 * t16)) + t39 * (-t5 * (-t11 *
     # (9 * t19 + t11) + t10 * (t11 + t31)) - 2 * t6 * (t22 * (t10 - t11
     # + t19) + t32))))
      t36 = t39 + t46
      t47 = 0.1q1 / t5
      t36 = 0.1q1 / t36
      t48 = 0.1q1 / t13
      t8 = 0.1q1 / t8
      t49 = 0.1q1 / t2
      t50 = 0.1q1 / t14
      t51 = t36 ** 2
      t52 = 0.1q1 / t3 ** 2
      t53 = 0.1q1 / t4 ** 2
      t1 = t25 * t1 ** 2 * (t46 * t39 * (t25 * t7 * t3 * (t24 * t1 * t30
     # * t33 + t27 * t14 * t9 * t12) - t4 * (t3 * (t25 * t7 * t17 * t44 
     #+ (t15 * (t16 * (13 * t19 - t12) + t15 * (t12 + t43)) - 3 * t38 * 
     #t20 + t29 * t21 * (t15 - t16) + t22 * ((2 * t16 + t15) * t20 + t29
     # * t25 * (-2 * t11 + t10))) * t13) + t40 * (t3 * (t14 * (t10 * (t1
     #1 - t31) + t11 * (-t11 + t43)) - t16 * t9 * t21) + t22 * (3 * t35 
     #* t9 ** 2 + 3 * t11 * t34 + t10 * (-6 * t11 + t19))))) + t2 * t45 
     #* t35 * (-t28 * t26 * t6 * t7 * t12 - t24 * (t4 * t6 * t17 * t20 -
     # t3 * t5 * (t16 * t21 - t23))) + t32 - t26 * t13 * t37 * (-t22 * t
     #25 * t35 * (t13 ** 2 * t18 ** 2 - t38) - t4 * (t3 * (t41 * t17 - t
     #29 * (t22 * (t34 - 5 * t10) + t42)) + t4 * (-t10 * t42 + t23 * (t1
     #0 + t11))))) * t49 * t53 * t48 * t50 * t47 * t52 * t8
      ret = -t1 * t36 * t51

      hjetmass_triangle_pppm_s14_0_0 = ret/32q0/(0,1q0)
      return

      end function
