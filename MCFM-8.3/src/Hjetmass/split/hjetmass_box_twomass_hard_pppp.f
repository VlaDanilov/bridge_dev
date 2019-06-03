
      complex*32 function hjetmass_box_twomass_hard_pppp
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
      t5 = za(i1, i2)
      t6 = za(i1, i3)
      t7 = zb(i2, i1)
      t8 = t6 * t2
      t9 = za(i1, i4) * t4
      t10 = t9 + t8
      t11 = t5 * t7
      t12 = t11 * t10
      t13 = za(i3, i4)
      t14 = zb(i4, i3)
      t15 = zb(i3, i1)
      t16 = zb(i4, i1)
      t17 = t1 * t15
      t18 = t3 * t16
      t19 = t18 + t17
      t20 = mt ** 2
      t21 = t1 * t2
      t22 = t3 * t4
      t23 = t13 * t14
      t24 = t11 * (t23 + t21 + t22)
      t25 = 16 * t20 * t19 * t12 + 4 * t24 ** 2
      t25 = sqrt(t25)
      t12 = 0.1q1 / t12
      t26 = t22 + t21
      t8 = t11 * t12 * (t8 * t26 + t9 * t26 + t23 * (t9 + t8))
      t9 = 2
      t26 = (0.1q1 / 0.2q1)
      t10 = t26 * t25 * t12 * t10
      t23 = t9 * (t23 + t21 + t22)
      t27 = -t8 + t10 + t23
      t24 = t9 * t24
      t28 = t25 - t24
      t19 = 0.1q1 / t19
      t29 = (0.1q1 / 0.4q1)
      t30 = t17 * t26
      t31 = t2 * (t6 * t28 * t12 * t29 + t1) - t30 * t27 * t19
      t24 = -t25 - t24
      t8 = -t8 - t10 + t23
      t6 = t2 * (t6 * t24 * t12 * t29 + t1) - t30 * t8 * t19
      t10 = 0.1q1 / t1
      t23 = 0.1q1 / t3
      t6 = 0.1q1 / t6
      t13 = 0.1q1 / t13
      t25 = 0.1q1 / t31
      t30 = t16 * t10
      t31 = t15 * t23
      t17 = t17 * t23
      t32 = (t17 + t16) * t2
      t18 = (t18 * t10 + t15) * t4
      t33 = t11 + t21
      t34 = (t11 * (t31 + t30) + t32 + t18) * t13
      t35 = t15 * (t23 * (t4 * t5 * t15 + t14 * t33) + t4 * (t30 * t5 + 
     #t14))
      t36 = 0.1q1 / t5
      t37 = t6 + t25
      t38 = t25 * t28
      t39 = t6 * t24
      t40 = (t28 + t24) * t13
      t41 = t4 * t10
      t42 = t2 * t14
      t17 = t17 * t37
      t43 = t23 * (t11 * t10 + t2) + t41
      t44 = t24 ** 2
      t45 = t25 * t28 ** 2 * t43
      t43 = t43 * t6 * t44
      t11 = t11 * t12 ** 2 * t4
      t46 = t39 + t38
      t47 = t7 * t14
      t48 = t2 * t16
      t49 = -t48 - t47
      t50 = t28 + t24
      t51 = t2 * t4
      t52 = t50 * t13
      t53 = t10 * t49 - t31 * t2
      t5 = t12 * ((t38 * t53 + t39 * t53) * t4 * t20 + t42 * t23 * t7 * 
     #(t38 * t33 + t39 * t33)) + t12 * t7 * (t2 * (t21 * t50 * t13 * t23
     # + t4 * t14 * t46) + t5 * (t10 * (t4 * (t48 * t46 + t52 * t7) + (t
     #51 * t44 * t12 * t6 + t39 * t49 + t38 * (t2 * (t4 * t28 * t12 - t1
     #6) - t47)) * t23 * t20) + t2 * t7 * t23 * t13 * t50 + t51 * t31 * 
     #t46) + t52 * t4 ** 2 * t10 * t3)
      t33 = (t25 * t27 + t6 * t8) * t14
      t44 = (t27 + t8) * t13
      ret = -t26 * t12 * t19 * t7 * (t8 * t24 * (t35 * t6 + t34) + (t35 
     #* t25 + t34) * t28 * t27) - t29 * t11 * t2 * (t43 + t45) + t9 * t7
     # * (t44 * t3 * t10 * t19 * t16 ** 2 + t19 * t15 ** 2 * (t44 + t33)
     # * t23 * t1 + t52 * t51 * t12 + t33 * t15 * t19 * t16) + 8 * t13 *
     # (t7 * (t15 * (t23 * (-t21 + t20) - t4) + t16 * (t10 * (-t22 + t20
     #) - t2)) + t20 * (t18 + t32) * t36) + t11 * t19 * t15 * (t45 * t27
     # + t43 * t8) / 8 - 4 * t7 * (t16 * (t42 * t37 + (-t27 - t8) * t13 
     #* t19 * t15) + t20 * (t12 * (t40 * t41 + (t14 * (t39 + t38) + t40)
     # * t23 * t2) - t31 * t14 * t37 - t30 * t14 * t37) + t42 * t17) - 4
     # * t42 * t20 * (-t16 * t37 - t17) * t36 + t5


      hjetmass_box_twomass_hard_pppp = ret/32q0/(0,1q0)
      return

      end function
