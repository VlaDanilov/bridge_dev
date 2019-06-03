
      complex*32 function hjetmass_triangle_pmpm_s23_mhsq_s14_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret

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
      t13 = t1 * t2
      t14 = t3 * t4
      t15 = t5 * t6
      t16 = t7 * t8
      t17 = t14 + t15 + t16 + t13
      t17 = -4 * t9 * t10 * t11 * t12 + t17 ** 2
      t15 = cg * sqrt(t17) + t13 + t14 + t15 + t16
      t16 = t2 * t5 + t4 * t7
      t17 = t15 / 2
      t18 = t9 * t10
      t19 = t17 * t7 + t18 * t2
      t13 = t14 + t13
      t14 = t9 * t12
      t20 = -t14 * t13 + t14 * t17
      t21 = t1 * t11 * t12 + t17 * t8
      t22 = t9 * t21
      t23 = t1 * t4 + t5 * t8
      t24 = (0.1q1 / 0.4q1)
      t25 = t15 ** 2
      t26 = t18 * t11 * t12
      t27 = t24 * t25 - t26
      t13 = -t15 * t13 / 2 + t26
      t26 = t9 * (-t3 * t11 * t12 + t17 * t6)
      t24 = t15 * t24
      t28 = t15 * t9
      t29 = t28 * (t24 * t4 - t5 * t11 * t12 / 2)
      t30 = -t5 * t11 * t12 + t17 * t4
      t31 = t10 * t30
      t7 = t7 * t11 * t12
      t32 = -t10 * (t17 * t2 + t7)
      t33 = t1 * t6 + t3 * t8
      t34 = t18 * t4
      t35 = t17 * t5 - t34
      t36 = t11 * t35
      t35 = t12 * t35
      t37 = t15 * t16
      t38 = t37 * t14
      t1 = t17 * t1 + t18 * t8
      t8 = -t1 * t11
      t3 = t11 * (t17 * t3 - t18 * t6)
      t6 = t14 * t16
      t14 = -t13
      t17 = t9 * t30
      t18 = t10 * t11 * t23
      t30 = t12 * t19
      t33 = 0.1q1 / t33
      t9 = 0.1q1 / t9
      t39 = 0.1q1 / t32
      t26 = 0.1q1 / t26
      t40 = 0.1q1 / t12
      t41 = 0.1q1 / t13
      t42 = 0.1q1 / t11
      t43 = 0.1q1 / t10
      t27 = 0.1q1 / t27 ** 2
      t44 = t20 * t23
      t45 = t43 * t27 * t42 * t40 * t9
      t46 = t45 * t16
      t3 = 0.1q1 / t3
      t47 = 0.1q1 / t20
      t48 = 0.1q1 / t30
      t49 = t31 ** 2
      t4 = t4 * t5
      t50 = t14 * t33
      t51 = t50 * t27
      t52 = t31 * t17
      t42 = t43 * t42 * t40 * t9
      t43 = t6 ** 2
      t53 = t15 * t35
      t54 = t15 * t36
      t30 = t42 * (t20 * (t54 * t30 * t33 * t9 * t39 ** 2 * t40 * t27 * 
     #t49 + t27 * (t54 * t17 * t26 * t41 + (-t53 * t36 * t9 * t40 + t18 
     #* t37) * t39 * t33) * t31) + t3 * (t53 * t8 ** 2 * t43 * t27 * t47
     # ** 2 + (t18 * (t38 * t14 * t17 + t15 * t8 * t43) - t36 * t35 * t3
     #8 * t8) * t27 * t47) - t4 * t37 * t22 * t41 * t26)
      t37 = 0.1q1 / t15
      ret = -2 * t46 * t15 * t25 * (t26 * (t16 * t22 ** 2 * t31 * t41 **
     # 2 - t44 * t36 * t41) - t11 * t19 * t20 * t22 * t23 * t41 * t26 **
     # 2 + t44 * t31 * t9 * t40 * t39 * t33) + 16 * t42 * ((t51 * t32 * 
     #t17 * t48 ** 2 - t51 * t36 * t48) * t35 ** 2 + t20 * (t29 * t33 * 
     #t9 * t39 * t40 * t27 * t49 + t4 * t33 * t9 * t39 * t40 * t31) + t5
     #1 * (t18 * t6 - t52) * t48 * t35 + t20 ** 2 * t29 * t31 * t10 * t2
     #1 * t9 ** 2 * t40 ** 2 * t27 * t39 * t33 ** 2 + t47 * t3 * (t36 * 
     #t15 * t12 * (t24 * t5 - t34 / 2) * t14 * t17 * t27 + (t14 * t18 * 
     #t28 * (t24 * t2 + t7 / 2) * t3 * t27 + t4) * t6 * t8)) - 8 * t30 +
     # 32 * t42 * t35 * t37 * t48 * t33 * (-t36 * t14 ** 2 * t12 * t1 * 
     #t27 * t33 + t4 * t13) + 4 * t46 * t25 * (t52 * t22 * t41 * t26 + t
     #50 * t35 * t23 * t48) - t45 * t25 ** 2 * t16 ** 2 * t22 * t23 * t4
     #1 * t26

      hjetmass_triangle_pmpm_s23_mhsq_s14_rat = ret/32q0/(0,1q0)
      return

      end function
