
      complex*32 function hjetmass_box_onemass_pppp
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
      t2 = za(i2, i3)
      t3 = zb(i2, i1)
      t4 = zb(i3, i2)
      t5 = za(i1, i3)
      t6 = zb(i3, i1)
      t7 = mt ** 2
      t8 = t2 * t4
      t9 = t8 * t1 * t3
      t10 = t9 * (4 * t7 * t5 * t6 + t9)
      t10 = sqrt(t10)
      t11 = -t10 + t9
      t12 = t1 * t3
      t13 = 0.1q1 / t3
      t1 = 0.1q1 / t1
      t14 = (0.1q1 / 0.2q1)
      t15 = t14 * sqrt(-t12 * t8 * (-8 * t7 * t5 * t6 - 2 * t9))*sqrt
     #(0.2q1) * t1 * t13
      t16 = t8 + t15
      t17 = zb(i4, i1)
      t18 = za(i2, i4)
      t19 = za(i1, i4)
      t20 = zb(i4, i2)
      t21 = 0.1q1 / t6
      t22 = 0.1q1 / t4
      t2 = 0.1q1 / t2
      t5 = 0.1q1 / t5
      t23 = t19 * t17
      t24 = t14 * (t11 * t1 * t5 * t13 * t22 * t19 * t20 + t16 * t2 * t2
     #1 * t18 * t17)
      t25 = t23 + t24
      t9 = t10 + t9
      t10 = t8 - t15
      t15 = zb(i4, i3)
      t26 = t18 * t20 + za(i3, i4) * t15
      t19 = t14 * (t9 * t1 * t5 * t13 * t22 * t19 * t20 + t10 * t2 * t21
     # * t18 * t17)
      t21 = -t26 + t19
      t24 = -t26 + t24
      t19 = t23 + t19
      t23 = 0.1q1 / t25
      t24 = 0.1q1 / t24
      t19 = 0.1q1 / t19
      t18 = 0.1q1 / t18
      t21 = 0.1q1 / t21
      t25 = t11 * t24
      t26 = t9 * t21
      t27 = t26 + t25
      t28 = t11 * t23
      t29 = t9 * t19
      t30 = t29 + t28
      t31 = -t19 + t21
      t32 = -t23 + t24
      t33 = t31 * t9
      t34 = t32 * t11
      t35 = t13 * t17
      t36 = t15 * t3
      t37 = t36 * t21
      t38 = t6 * t20
      t39 = t18 * t1
      t40 = t18 * t5
      t41 = t12 * t18
      t42 = t7 * t22 * t2
      t43 = -1 + t42
      t44 = t6 * t18
      t45 = t1 * t15
      t46 = t4 * t17
      t47 = t36 * t18
      t48 = t9 ** 2
      t49 = t11 ** 2
      t50 = (t23 - t24) * t49
      t51 = (t19 - t21) * t48
      t52 = t51 + t50
      t53 = -t19 + t21
      t11 = (-t23 + t24) * t16 * t11
      t54 = t9 * t10
      t55 = t25 * t16
      t11 = t5 * t2 * t22 * (t55 * t47 + (t11 * t20 + t54 * (t19 * (-t41
     # - t20) + t20 * t21)) * t2 * t17 + t52 * t1 * t5 * t13 * t20 ** 2)
     # + t5 * (t2 * (t22 * (t18 * (t50 * t20 * t5 + t51 * t20 * t5 + t10
     # * (t38 * t53 + t37) * t9 + t11 * t38) + t45 * t20 * (t54 * t53 + 
     #t11)) - t18 * (t28 * t16 + t54 * t19) * t17) + t39 * t15 * (t54 * 
     #t21 + t55) - t28 * t41 * t22 * t16 * t2 ** 2 * t17 + t39 * t5 * t2
     #0 * t52 * t13)
      t50 = t51 * t10 + t50 * t16
      t51 = t5 ** 2
      t52 = t18 * t50
      t53 = t24 + t21
      t54 = -t23 - t19
      t55 = t7 * t2
      t56 = t4 * t54
      t3 = t15 * t7 * (t3 * ((-t23 + t24 + t21) * t2 * t17 + t44 * (t24 
     #+ t21)) + t46 * (-t23 + t24 - t19 + t21) * t1) + t44 * (t17 * (t56
     # * t7 + t12 * (t55 * t54 - t56)) + t8 * (t53 * t1 * t7 - t3 * t53)
     # * t15) + t42 * t39 * (t31 * t48 + t32 * t49) * t51 * t13 * t20 + 
     #t36 * t17 * (t4 * (t23 - t24 + t19 - t21) - t55 * t19)
      ret = t14 * t51 * t2 * t22 * t20 * (t52 * t2 * t22 + (t50 * t2 * t
     #22 * t20 + t52) * t1 * t13) - 2 * t5 * (t20 * ((t9 * (t19 * (-t7 *
     # t13 * t1 + 1) - t21) - t34) * t2 * t17 + t45 * (-t34 - t33) + t44
     # * (-t28 * t43 - t33)) + t46 * t18 * t30 + t47 * (t25 * t43 - t26)
     #) - 2 * t40 * (-t8 * t27 * t1 * t15 - t25 * t6 * t20) - 2 * t41 * 
     #t2 * t17 * (t30 * t5 + t6 * (-t10 * t19 - t16 * t23)) - 2 * t5 * (
     #t2 * (t18 * (t17 * t27 + t9 * (-t38 * t19 + t37) * t22) + t20 * ((
     #t34 + t33) * t22 * t15 + t35 * (t26 + t34)) * t1) + t39 * (t27 * t
     #13 * t20 * t6 - t15 * t30) - t35 * t39 * t30 * t4) * t7 - 2 * t36 
     #* t16 * (t32 * t2 * t17 + t44 * t24) - 2 * t36 * (t31 * t2 * t17 +
     # t44 * t21) * t10 - t11 - 4 * t3 + 8 * t40 * t7 * (-t45 * (t26 + t
     #25) + (t29 + t28) * t2 * t17)

      hjetmass_box_onemass_pppp = ret/32q0/(0,1q0)
      return

      end function
