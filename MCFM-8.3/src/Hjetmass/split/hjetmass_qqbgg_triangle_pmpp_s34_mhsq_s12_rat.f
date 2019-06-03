
      complex*32 function hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          real*16 cg

          cg = 1q0

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
      t13 = t10 + t11 + t12 + t9
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t18 = t14 * t16
      t13 = -0.4q1 * t18 * t15 * t17 + t13 ** 2
      t13 = cg * sqrt(t13) + t10 + t11 + t12 + t9
      t19 = t3 * t2
      t20 = t7 * t6
      t21 = t20 + t19
      t22 = 0.1q1 / 0.2q1
      t23 = t22 * t13
      t24 = -t14 * (t5 * t16 * t17 + t23 * t4)
      t25 = t14 * (t3 * t16 * t17 + t23 * t6)
      t26 = 0.1q1 / 0.4q1
      t27 = t13 ** 2
      t28 = t18 * t15 * t17
      t29 = t26 * t27 - t28
      t30 = t11 + t9
      t31 = t15 * t2
      t32 = t30 * t16
      t33 = t27 * t14
      t19 = t20 + t19
      t20 = t5 ** 2 * t6 ** 2
      t34 = t1 ** 2 * t2 ** 2
      t35 = t4 * t19 * t1
      t19 = t8 * t19 * t5
      t36 = t18 * t13
      t37 = t17 ** 2
      t30 = t30 * t14
      t38 = t30 * t17
      t39 = t1 * t7
      t40 = t11 * t9 * t13
      t41 = t22 * t36 * ((t18 + t9 + t12) * t17 * t15 + t20 + t34 + t35 
     #+ t19)
      t42 = -t26 * t33 * (t31 * t8 + t32) + t18 * ((-t39 * t37 - t38) * 
     #t16 * t15 + t40) + t41
      t43 = t18 * t21
      t44 = t14 * t15
      t45 = t17 * (t23 * t5 + t44 * t4)
      t46 = t44 * t6
      t47 = -t17 * (t23 * t3 + t46)
      t48 = -t11 - t9
      t5 = t3 * t5
      t49 = t27 * t16
      t50 = t15 ** 2
      t32 = t32 * t17 * t15
      t4 = t4 * t6
      t19 = t22 * t36 * ((t18 + t10 + t11) * t17 * t15 + t20 + t34 + t35
     # + t19)
      t20 = t26 * t49 * (t14 * t48 + t5 * t17) + t18 * (t14 * (t4 * t17 
     #* t50 - t32) + t40) + t19
      t30 = -t26 * t49 * (t39 * t17 + t30) + t41 + t18 * (t14 * (-t2 * t
     #8 * t17 * t50 - t32) + t40)
      t4 = t26 * t33 * (t4 * t15 + t16 * t48) + t19 + t18 * ((t5 * t37 -
     # t38) * t16 * t15 + t40)
      t2 = t14 * (-t7 * t16 * t17 + t23 * t2)
      t5 = t11 + t9
      t9 = t22 * t13 * t5 - t28
      t19 = t23 * t18
      t5 = -t18 * t5 + t19
      t31 = t31 * t14
      t32 = t17 * (t23 * t7 - t31)
      t26 = t13 * t26
      t33 = t13 * t17
      t10 = t10 + t12
      t6 = t1 * t6 + t3 * t8
      t34 = t15 * t17
      t11 = t11 + t12
      t12 = 0.1q1 / t7
      t35 = 0.1q1 / t16
      t4 = 0.1q1 / t4
      t36 = 0.1q1 / t3
      t37 = 0.1q1 / t42
      t38 = 0.1q1 / t14
      t39 = t4 - t37
      t40 = 0.1q1 / t29 ** 2
      t41 = t2 * t47
      t40 = t38 * t36 * t35 * t40 * t12
      t42 = t40 * t5
      t48 = t37 ** 2
      t49 = t21 ** 2
      t50 = t4 ** 2
      t30 = 0.1q1 / t30
      t20 = 0.1q1 / t20
      t51 = -t30 + t20
      t52 = t20 ** 2
      t53 = t30 ** 2
      t18 = t47 * (-t18 * t10 + t19)
      t15 = 0.1q1 / t15
      t19 = -0.1q1 / t29
      t19 = t19 ** 2
      t10 = t25 * (t22 * t13 * t10 - t28)
      ret = -0.8q1 * t40 * t43 * t13 * (t9 * (t51 * t32 * t25 + t41 * t5
     #1) + (-t20 * t47 + t30 * t47 + t9 * (t25 * (-t34 * t11 + t34 * t23
     #) + t18) * t53 + t9 * (t2 * t34 * t6 - t18) * t52) * t45 * t43 + t
     #9 * t17 * (-t23 * t1 + t44 * t8) * (-t53 + t52) * t45 * t43 ** 2) 
     #- 0.64q2 * t43 * t38 * t15 * t35 * t19 * (t33 * (-t31 * t22 + t26 
     #* t7) * t12 - t33 * (t46 * t22 + t26 * t3) * t36) + 0.2q1 * t40 * 
     #t24 * t49 * t13 * t27 * (-t25 * t37 + t25 * t4 - t5 * (t47 * (t22 
     #* t13 * t11 - t28) + t10) * t48 + t10 * t5 * t50) + 0.16q2 * t21 *
     # t27 * t38 * t15 * t35 * t19 * (t12 * t2 + t25 * t36) + 0.4q1 * t4
     #2 * t21 * t27 * (t39 * t32 * t25 + t41 * t39) - t42 * t24 * t49 * 
     #t27 ** 2 * (-t6 * t32 * t50 + (-t50 + t48) * t14 * (t1 * t16 * t17
     # - t23 * t8) * t21)

      hjetmass_qqbgg_triangle_pmpp_s34_mhsq_s12_rat = ret/32q0*(0,1q0)
      return

      end function
