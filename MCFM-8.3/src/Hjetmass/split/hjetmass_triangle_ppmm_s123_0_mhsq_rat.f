
      complex*32 function hjetmass_triangle_ppmm_s123_0_mhsq_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = za(i3, i4)
      t2 = za(i1, i4)
      t3 = zb(i2, i1)
      t4 = zb(i3, i2)
      t5 = -t1 * t4 + t2 * t3
      t6 = zb(i4, i1)
      t7 = za(i1, i3)
      t8 = za(i1, i2)
      t9 = zb(i3, i1)
      t10 = za(i2, i3)
      t11 = t8 * t3
      t12 = t7 * t9
      t13 = t10 * t4 + t11 + t12
      t14 = zb(i4, i2)
      t15 = za(i2, i4)
      t16 = zb(i4, i3)
      t17 = t2 * t6
      t18 = t15 * t14
      t19 = t1 * t16
      t20 = t19 + t17 + t18
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t17 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t17 = 0.1q1 / t17
      t18 = -2 * t2 * t13
      t19 = t18 * t14 * t17 + t4 * t7
      t11 = t18 * t6 * t17 + t11 + t12
      t12 = t15 * t19
      t20 = t6 * (t11 * t2 + t12)
      t21 = t15 * t4 + t2 * t9
      t22 = t1 * t9 + t15 * t3
      t23 = -2 * t1 * t13 * t6 * t17 - t10 * t3
      t9 = (-2 * t15 * t13 * t6 * t17 + t10 * t9) * t14
      t24 = t2 * (t11 * t6 + t9)
      t18 = t18 * t16 * t17 - t4 * t8
      t9 = t2 * (t16 * t23 + t9)
      t12 = t6 * (t1 * t18 + t12)
      t20 = 0.1q1 / t20
      t18 = 0.1q1 / t18
      t25 = 0.1q1 / t14
      t12 = 0.1q1 / t12
      t8 = 0.1q1 / t8
      t16 = 0.1q1 / t16
      t9 = 0.1q1 / t9
      t4 = 0.1q1 / t4
      t26 = 0.1q1 / t7
      t27 = 0.1q1 / t2
      t10 = 0.1q1 / t10
      t24 = 0.1q1 / t24
      t28 = 0.1q1 / t6
      t29 = 0.1q1 / t11
      t30 = t1 * t6
      t31 = t30 * t20
      t32 = t31 - t18
      t33 = t23 ** 2
      t34 = t23 * t33
      t35 = t14 * t26
      t36 = t16 * t9
      t20 = t1 * t20 * t27
      t12 = t18 * t12
      t27 = t30 * t26
      t37 = t10 * t4
      t38 = t14 ** 2
      t39 = t16 ** 2
      t40 = t24 ** 2
      t41 = t2 ** 2
      t15 = t5 * t15
      t42 = t6 * t10
      t43 = t29 * t22
      t44 = t38 * t26
      t45 = t44 * t10
      t46 = t14 * t22 * t8
      t47 = t8 * t1
      t48 = t17 * t4
      t49 = t2 * t23
      t50 = t49 * t24
      t51 = -t50 + t16
      t52 = t19 * t6 * t25
      t11 = t52 + t11
      t53 = t23 * t24
      t54 = t37 * t17 * t1
      t11 = t54 * (t5 * (t8 * (-t11 * t18 + t31 * t11) + t27 * t19 * (t1
     #2 * t11 + t20 * (t52 * t29 + 1))) - t49 * t44 * (t53 * t28 * t29 +
     # t36) * t22 + t46 * t51)
      t13 = 0.1q1 / t13
      t18 = t35 * t5
      t22 = t48 * t47 * t3 * t5 * t13 * (t42 * t7 * t25 + 1)
      ret = 8 * t37 * t3 * t1 * (-t16 * t8 + t2 * (t23 * (t24 * t8 + t36
     # * t35) + t35 * t24 * t28 * t29 * t33) + t27 * (t20 * t29 + t12) *
     # t25 * t19 ** 2 + t8 * t32 * t25 * t19) - 32 * t48 * ((t45 * t43 *
     # t40 * t34 + t46 * t40 * (t42 + t35) * t33) * t2 * t41 + t41 * (t4
     #5 * t24 * t29 * (-t15 * t24 + t43 * t28) * t34 + t14 * (-t42 * t15
     # * t8 * t40 + t35 * (-t15 * t8 * t40 + t42 * (-t1 * t21 - t15) * t
     #9 ** 2 * t16)) * t33) - t27 * t2 * t21 * t23 * t38 * t10 * t39 * t
     #9 + t47 * (t42 * t21 * t14 * t39 + t5 * t19 * t26 * t32 + t21 * t3
     #8 * t26 * t39)) + 48 * t37 * t30 * t5 * t17 * t8 * (-t50 + t16) - 
     #16 * t11 + 64 * t47 * t4 * (t7 * t3 ** 2 * t10 * t13 * t25 + t18 *
     # t17 * t51) - 80 * t18 * t54 * t49 * (t53 * t29 + t36 * t6) - 128 
     #* t22

      hjetmass_triangle_ppmm_s123_0_mhsq_rat = ret/32q0/(0,1q0)
      return

      end function
