
      complex*32 function hjetmass_box_ppmm_0_0_s34_mhsq_s12_s234
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
      t2 = zb(i2, i1)
      t3 = zb(i3, i2)
      t4 = za(i2, i4)
      t5 = zb(i4, i2)
      t6 = za(i3, i4)
      t7 = zb(i4, i3)
      t8 = zb(i3, i1)
      t9 = zb(i4, i1)
      t10 = t1 * t8
      t11 = t4 * t9 + t10
      t12 = mt ** 2
      t13 = za(i1, i3)
      t14 = za(i1, i4)
      t15 = t13 * t3
      t16 = t14 * t5
      t17 = t16 + t15
      t18 = za(i1, i2) * t2
      t19 = t18 * t17
      t20 = t1 * t3
      t21 = t4 * t5
      t22 = t6 * t7
      t23 = t18 * (t22 + t20 + t21)
      t24 = 16 * t12 * t11 * t19 + 4 * t23 ** 2
      t24 = sqrt(t24)
      t23 = 2 * t23
      t25 = -t23 + t24
      t19 = 0.1q1 / t19
      t26 = t21 + t20
      t15 = t18 * t19 * (t15 * t26 + t16 * t26 + t22 * (t16 + t15))
      t26 = (0.1q1 / 0.2q1)
      t20 = 2 * t22 + 2 * t20 + 2 * t21
      t17 = t26 * t24 * t19 * t17
      t21 = t20 - t15 + t17
      t11 = 0.1q1 / t11
      t22 = (0.1q1 / 0.4q1)
      t27 = t10 * t26
      t28 = t3 * (t13 * t25 * t19 * t22 + t1) - t27 * t21 * t11
      t15 = t20 - t15 - t17
      t17 = -t23 - t24
      t20 = t3 * (t13 * t17 * t19 * t22 + t1) - t27 * t15 * t11
      t23 = t17 * t22 * t19
      t24 = t15 * t26 * t11 * t4
      t27 = -t23 * t14 * t3 + t24 * t8
      t29 = t25 * t22 * t19
      t30 = t26 * t21 * t11 * t4
      t31 = -t29 * t14 * t3 + t30 * t8
      t32 = t13 * t4
      t33 = t1 * t14 - t32
      t7 = 0.1q1 / t7
      t34 = 0.1q1 / t3
      t35 = 0.1q1 / t14
      t20 = 0.1q1 / t20
      t28 = 0.1q1 / t28
      t5 = 0.1q1 / t5
      t36 = t15 ** 2
      t37 = t21 ** 2
      t38 = t36 * t17
      t39 = t25 * t37 + t38
      t40 = t11 ** 2
      t41 = t2 ** 2
      t42 = t2 * t4
      t43 = t28 * t37
      t44 = t7 * t4
      t45 = t40 * t41
      t46 = t43 * t25
      t38 = t38 * t20
      t47 = t19 * t6
      t48 = t6 ** 2
      t49 = t5 * t35
      t50 = 0.1q1 / t4
      t51 = t28 * t6 + t7
      t52 = t20 * t6 + t7
      t53 = t7 * t35 * t6
      t54 = t48 * t28
      t55 = t2 * t7
      t56 = t3 * t6
      t57 = t21 * t37 * t28
      t32 = t55 * t11 * t5 * (t32 * t35 + t1)
      t58 = t48 * t20
      t59 = t11 * t41
      t3 = t59 * (t19 * (t25 * (t21 * (t5 * (t51 * t14 * t2 + t56 * t51)
     # * t50 * t1 - t54 + (t3 * (-t54 * t35 - t53) - t55) * t5 * t13) + 
     #t57 * t45 * t1 * t4 * t34 * t5 + t32 * t37) + t17 * t15 * (t32 * t
     #15 - t58 + (t3 * (-t58 * t35 - t53) - t55) * t5 * t13 + t45 * t1 *
     # t4 * t36 * t34 * t20 * t5 + t5 * (t52 * t14 * t2 + t56 * t52) * t
     #50 * t1)) + t42 * t49 * t10 * t6 * t34 * t40 * (t20 * t15 * t36 + 
     #t57))
      t10 = t21 + t15
      t13 = t34 * t10
      t14 = t20 + t28
      t32 = t5 * t50
      t14 = t32 * t48 * t2 * t12 * (t56 * t14 * t35 + t14 * t2)
      t32 = t32 * t55 * t12 * t6 * (t56 * t35 + t2)
      ret = -t22 * t49 * t11 * t41 * (-t7 * (t15 ** 2 * t17 * t19 * t11 
     #* t2 * t33 + t21 ** 2 * t25 * t19 * t11 * t2 * t33) + t18 * (t47 *
     # t7 * t39 * t11 + (t38 + t46) * t19 * t11 * t48)) + t26 * t3 + (0.
     #3q1 / 0.4q1) * t47 * t2 * t41 * t40 * t1 * t5 * (t38 + t46) + 8 * 
     #t32 - t45 * (t34 * (t2 * (t43 * (t6 * (t31 * t35 + t8) + t42) + (t
     #6 * (t27 * t35 + t8) + t42) * t20 * t36) * t5 * t1 + t44 * (-t19 *
     # t2 * t39 + (-t37 - t36) * t35 * t8 * t6)) + t44 * t6 * t35 * t5 *
     # (t37 + t36) * t9) + 2 * t59 * (t7 * (t2 * (-t13 * t4 + (-t21 - t1
     #5) * t5 * t1) + t6 * (t35 * (t34 * (-t27 * t15 - t21 * t31) + t5 *
     # (t15 * (-t23 * t16 + t24 * t9) + t21 * (-t29 * t16 + t30 * t9))) 
     #+ t9 * t5 * t10 - t13 * t8)) + t49 * t48 * (t15 * t20 + t21 * t28)
     # * t12) + 4 * t14

      hjetmass_box_ppmm_0_0_s34_mhsq_s12_s234 = ret/32q0/(0,1q0)
      return

      end function
