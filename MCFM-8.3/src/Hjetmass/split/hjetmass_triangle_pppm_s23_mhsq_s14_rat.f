
      complex*32 function hjetmass_triangle_pppm_s23_mhsq_s14_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

          parameter (cg = 1q0)

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = za(i3, i4)
      t8 = zb(i4, i3)
      t9 = za(i1, i4)
      t10 = za(i2, i3)
      t11 = zb(i3, i2)
      t12 = zb(i4, i1)
      t13 = t3 * t4
      t14 = t5 * t6
      t15 = t7 * t8
      t16 = t1 * t2
      t17 = t13 + t14 + t15 + t16
      t17 = -4 * t9 * t10 * t11 * t12 + t17 ** 2
      t17 = cg * sqrt(t17) + t13 + t14 + t15 + t16
      t18 = t2 * t5 + t4 * t7
      t19 = (0.1q1 / 0.2q1)
      t20 = t19 * t17
      t21 = t9 * t10
      t22 = t21 * t2 + t20 * t7
      t23 = -t11 * t22
      t24 = t17 ** 2
      t25 = t24 ** 2
      t26 = t21 * t11 * t12
      t27 = -t24 / 4 + t26
      t13 = t13 + t16
      t16 = -t19 * t17 * t13 + t26
      t28 = t20 * t5 - t21 * t4
      t29 = t11 * t28
      t30 = t5 * t11 * t12
      t31 = t20 * t4 - t30
      t32 = t9 * t31
      t33 = t9 * t12
      t34 = t33 * t20
      t13 = -t33 * t13 + t34
      t35 = t7 * t11 * t12
      t36 = t20 * t2 + t35
      t37 = -t16
      t14 = t14 + t15
      t15 = t1 * t6 + t3 * t8
      t28 = t12 * t28
      t31 = t10 * t31
      t38 = -t10 * t36
      t39 = t9 * (t1 * t11 * t12 + t20 * t8)
      t40 = t2 * t3 + t6 * t7
      t41 = t1 * t4 + t5 * t8
      t42 = t9 * (-t3 * t11 * t12 + t20 * t6)
      t43 = t33 * t18
      t3 = t11 * (t20 * t3 - t21 * t6)
      t6 = t17 * t18
      t44 = t17 / 4
      t45 = t17 * t9
      t26 = t19 * t17 * t14 - t26
      t2 = t45 * (t35 * t19 + t44 * t2)
      t22 = t12 * t22
      t35 = t10 * t11
      t1 = -t11 * (t20 * t1 + t21 * t8)
      t8 = 0.1q1 / t12
      t11 = 0.1q1 / t13
      t12 = 0.1q1 / t15
      t15 = -0.1q1 / t27
      t20 = 0.1q1 / t38
      t21 = 0.1q1 / t28
      t10 = 0.1q1 / t10
      t5 = 0.1q1 / t5
      t16 = 0.1q1 / t16
      t28 = 0.1q1 / t9
      t38 = 0.1q1 / t39
      t7 = 0.1q1 / t7
      t46 = t23 * t7
      t47 = t29 * t5
      t48 = t47 + t46
      t49 = t32 * t5
      t9 = t9 * t36 * t7
      t36 = -t9 - t49
      t50 = t13 ** 2
      t51 = t28 ** 2
      t52 = t28 * t51
      t15 = t15 ** 2
      t53 = t8 ** 2
      t54 = t8 * t53
      t27 = 0.1q1 / t27 ** 2
      t55 = t12 ** 2
      t56 = t18 * t24
      t22 = 0.1q1 / t22
      t57 = 0.1q1 / t17
      t58 = t43 * t22 * t5
      t59 = t12 * t28 * t8
      t60 = 0.1q1 / t1
      t61 = 0.1q1 / t3
      t62 = t2 * t7
      t4 = t45 * (-t30 * t19 + t44 * t4) * t5 + t62
      t30 = t43 ** 2
      t44 = t37 ** 2
      t45 = t18 ** 2
      t63 = t7 * t5 * t15
      t2 = (t63 * t43 * ((t44 * t43 * t35 * t40 * t60 * t22 ** 2 + t23 *
     # t44 * t60 * t22) * t12 * t32 + t21 * t11 * t43 * (-t17 * t29 * t4
     #3 + (-t43 * t35 * t41 * t21 - t29) * t61 * t2 * t37)) + (t4 * t12 
     #* t37 * t18 - t17 * t30 * t11 * (t47 + t46)) * t10 * t27) * t28 * 
     #t8 - t62 * t45 * t37 * t5 * t22 * t15 * t12 - t50 * t26 * t55 * t2
     #7 * t10 * t4 * t52 * t54
      t4 = t28 * t8
      t22 = 0.1q1 / t42
      t31 = 0.1q1 / t31
      t23 = t23 * t13
      t35 = t11 ** 2
      t1 = t4 * (t63 * t17 * t30 ** 2 * t1 * t35 * t21 + (t17 * t43 * t3
     #0 * t35 * (t1 * t5 + t3 * t7) + (t4 * t56 * t36 + t48 * t6) * t12 
     #* t13) * t10 * t27)
      ret = t19 * t4 * t63 * t17 * t25 * t18 * t45 * t16 * t31 * (t23 * 
     #t41 * t31 * t22 + t18 * t16 * t39) - 32 * t59 * t37 * (t43 * (t48 
     #* t10 * t27 - t58 * t46 * t15) + t15 * (t47 * t10 + t46 * (-t58 + 
     #t10)) * t57 * t12 * (-t33 * t14 + t34) * t37) + 16 * t2 + 2 * t63 
     #* t31 * t16 * t45 * t17 * t24 * (t18 * t29 + t4 * (t23 * t22 - t6)
     # * t32) + 8 * t1 - 4 * t56 * (t18 * (-t49 * t43 * t7 * t11 * t21 *
     # t15 + (t11 * t36 + t16 * t48) * t10 * t27) + t9 * t13 * t12 * t15
     # * t20 * t5 * (-t29 * t13 * t38 + t6) * t51 * t53 + t9 * t50 * t26
     # * t52 * t5 * t54 * t15 * t20 * t55 + t18 * (-t46 * t13 * t12 * t1
     #5 * t20 * t5 + t36 * t16 * t10 * t27 * t17) * t28 * t8) + t4 * t45
     # * t25 * (t59 * t47 * t50 * t40 * t7 * t15 * t20 ** 2 * t38 + (-t3
     #9 * t5 - t42 * t7) * t16 ** 2 * t10 * t27 * t18)

      hjetmass_triangle_pppm_s23_mhsq_s14_rat = ret/32q0/(0,1q0)
      return

      end function
