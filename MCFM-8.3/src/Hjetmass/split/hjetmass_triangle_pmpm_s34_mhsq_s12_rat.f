
      complex*32 function hjetmass_triangle_pmpm_s34_mhsq_s12_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret

          parameter (cg = 1q0)

      t1 = za(i1, i3)
      t2 = zb(i3, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i1)
      t7 = za(i2, i4)
      t8 = zb(i4, i2)
      t9 = t1 * t2
      t10 = t3 * t4
      t11 = t5 * t6
      t12 = t7 * t8
      t13 = t9 + t10 + t11 + t12
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t13 = -4 * t14 * t15 * t16 * t17 + t13 ** 2
      t10 = cg * sqrt(t13) + t10 + t11 + t12 + t9
      t12 = (0.1q1 / 0.2q1)
      t13 = t12 * t10
      t18 = t14 * t15
      t19 = t18 * t2
      t20 = -t13 * t7 + t19
      t21 = t16 * t20
      t9 = t9 + t11
      t11 = t18 * t16 * t17
      t22 = t12 * t10 * t9 - t11
      t20 = -t17 * t20
      t23 = t13 * t3 + t18 * t6
      t24 = -t16 * t23
      t25 = t1 * t4 + t5 * t8
      t26 = (0.1q1 / 0.4q1)
      t27 = t10 ** 2
      t28 = t26 * t27 - t11
      t29 = t2 * t3 + t6 * t7
      t30 = t14 * t16
      t31 = t30 * t29
      t32 = t13 * t5 + t18 * t4
      t33 = t17 * t32
      t30 = t30 * t13 - t30 * t9
      t18 = t17 * (-t13 * t1 + t18 * t8)
      t34 = t2 * t5 + t4 * t7
      t35 = t10 * t15
      t36 = t35 * t17 * t34
      t37 = t7 * t16 * t17
      t38 = -t13 * t2 + t37
      t39 = t15 * t38
      t3 = t3 * t16 * t17 + t13 * t6
      t6 = t15 * t3
      t40 = t15 * t17 * t34
      t4 = t5 * t16 * t17 + t13 * t4
      t5 = -t14 * t4
      t1 = t14 * (t1 * t16 * t17 - t13 * t8)
      t8 = -t22
      t38 = -t14 * t38
      t41 = t10 * t26
      t35 = t35 * (t37 * t12 - t41 * t2)
      t37 = 0.1q1 / t17
      t42 = 0.1q1 / t24
      t1 = 0.1q1 / t1
      t43 = 0.1q1 / t16
      t44 = 0.1q1 / t15
      t25 = 0.1q1 / t25
      t45 = 0.1q1 / t6
      t46 = 0.1q1 / t8
      t47 = 0.1q1 / t14
      t48 = 0.1q1 / t30
      t18 = 0.1q1 / t18
      t28 = 0.1q1 / t28 ** 2
      t49 = t43 * t25
      t50 = t49 * t45 * t47
      t51 = t29 * t5
      t44 = t47 * t44
      t52 = t44 * t28 * t43 * t37
      t2 = t2 * t7
      t53 = t44 * t43 * t37
      t4 = t53 * (t22 * t28 * (t20 * t21 ** 2 + t31 * (t21 * t40 - t34 *
     # t35)) * t25 * t42 + t21 * t34 * t31 * (t26 * t27 * t9 - t11 * t13
     #) * t6 * t42 ** 2 * t25 * t28 + t2 * t50 * t30 * t39 + t30 ** 2 * 
     #t38 * t15 * t4 * t35 * t47 ** 2 * t43 ** 2 * t25 ** 2 * t28 * t45 
     #+ t31 * t33 * t48 * t18 * (-t21 * t10 * t17 * (-t19 * t12 + t41 * 
     #t7) * t28 + t2))
      t6 = 0.1q1 / t10
      ret = -8 * t53 * ((t10 * t21 * t33 ** 2 * t28 * t18 * t48 ** 2 + t
     #33 * t36 * t28 * t18 * t48) * t31 ** 2 + t30 * (t50 * t10 * t38 * 
     #t28 * t39 ** 2 + t36 * t29 * t28 * t25 * t45 * t39) + (t22 * t33 *
     # t36 * t14 * t3 * t28 * t18 ** 2 - t22 * t36 * t38 * t28 * t18) * 
     #t48 * t31 - t2 * t51 * t10 * t1 * t46) + 16 * t4 - 2 * t52 * t29 *
     # t10 * t27 * (t1 * (t29 * t39 * t5 ** 2 * t46 ** 2 + t20 * t30 * t
     #34 * t46) + t50 * t34 * t30 * t39 + t30 * t34 * t5 * t17 * t23 * t
     #1 ** 2 * t46) - 32 * t49 * t44 * t21 * t6 * t42 * t37 * (t16 * t32
     # * t22 ** 2 * t20 * t25 * t28 - t2 * t8) + 4 * t52 * t27 * (-t51 *
     # t38 * t39 * t1 * t46 + t29 * (t38 * t1 * t46 + t50 * (-t24 * t39 
     #* t45 + t21)) * t40 * t30 + (-t20 * t31 * t48 * t18 + t21 * t29 * 
     #t42 * t25) * t34 * t22) - t52 * t27 ** 2 * t34 * t29 ** 2 * t5 * t
     #1 * t46

      hjetmass_triangle_pmpm_s34_mhsq_s12_rat = ret/32q0/(0,1q0)
      return

      end function
