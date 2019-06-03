
      complex*32 function hjetmass_triangle_ppmm_s234_0_mhsq_rat
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
      t2 = zb(i2, i1)
      t3 = za(i1, i2)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = zb(i4, i3)
      t7 = za(i1, i3)
      t8 = zb(i3, i1)
      t9 = za(i2, i4)
      t10 = zb(i4, i1)
      t11 = t10 * t9 + t4 * t8
      t12 = zb(i4, i2)
      t13 = t4 * t5
      t14 = t9 * t12
      t15 = t1 * t6 + t13 + t14
      t16 = za(i1, i4)
      t17 = t3 * t2
      t18 = t7 * t8
      t19 = t16 * t10
      t20 = t19 + t17 + t18
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t18 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t18 = 0.1q1 / t18
      t19 = -2 * t7 * t2 * t15 * t18 + t1 * t12
      t13 = -2 * t17 * t15 * t18 + t13 + t14
      t14 = -2 * t16 * t2 * t15 * t18 - t1 * t5
      t20 = t8 * t19
      t21 = t2 * t13
      t22 = t3 * (t21 + t20)
      t23 = t10 * t14
      t24 = t3 * (t23 + t20)
      t25 = -2 * t3
      t8 = t7 * (t25 * t8 * t15 * t18 + t6 * t9)
      t9 = t2 * (t13 * t3 + t8)
      t25 = t25 * t10 * t15 * t18 - t4 * t6
      t8 = t2 * (t16 * t25 + t8)
      t26 = 0.1q1 / t12
      t4 = 0.1q1 / t4
      t22 = 0.1q1 / t22
      t27 = 0.1q1 / t7
      t8 = 0.1q1 / t8
      t9 = 0.1q1 / t9
      t24 = 0.1q1 / t24
      t28 = 0.1q1 / t3
      t29 = 0.1q1 / t14
      t5 = 0.1q1 / t5
      t30 = 0.1q1 / t16
      t10 = 0.1q1 / t10
      t25 = 0.1q1 / t25
      t31 = 0.1q1 / t13
      t6 = 0.1q1 / t6
      t32 = t19 ** 2
      t33 = t13 ** 2
      t34 = t7 * t26
      t35 = t22 * t31
      t36 = t17 * t24
      t37 = t17 * t26
      t38 = t28 * t6
      t39 = t4 * t5
      t40 = t39 * t2
      t41 = t9 ** 2
      t42 = t7 ** 2
      t43 = t2 ** 2
      t44 = t31 ** 2
      t45 = t28 ** 2
      t46 = t20 * t8 ** 2
      t47 = t2 * t8
      t48 = t47 * t6
      t49 = t26 * t32
      t50 = t31 * t6
      t51 = t16 * t5
      t52 = t18 * t4 * t11 * t2
      t15 = 0.1q1 / t15
      t53 = t28 * t31
      t54 = t34 * t11 * t19 * t18
      t3 = t43 * t3
      t11 = t40 * t18 * t11
      t3 = t11 * (t31 * (t3 * t49 * t24 * t10 + t36 * t19 * t14 * t6 - t
     #38 * t7 * t14) + t19 * t14 * (t37 * t19 * t22 - t6) * t44 + t42 * 
     #t2 * t14 * t28 * t26 * t9 + t47 * t7 * t14 * (t34 * t13 * t30 * t2
     #5 + t6) + t32 * t31 * (t3 * t19 * t26 * t24 * t29 * t10 - t50 + t1
     #7 * (t35 * t19 * t26 + t24 * t6)) * t27 * t16)
      t14 = t54 * t39 * t43 * (t13 * t25 * t8 + t16 * t28 * t9)
      ret = -8 * t40 * t1 * (t2 * (t13 * (-t34 * t9 * t28 - t6 * t8) - t
     #34 * t8 * t30 * t25 * t33) + t6 * (-t36 + t31) * t27 * t19 - t37 *
     # (t2 * t24 * t29 * t10 + t35) * t27 * t32 + t38) - 32 * t52 * (-t4
     #9 * t6 * t44 + t26 * (t2 * (t13 * (t16 * (-t20 * t28 * t5 * t41 - 
     #t45 * t5 * t9) - t46 * t6) - t46 * t5 * t25 * t33) + t45 * t6 - t4
     #3 * t16 * t33 * t41 * t28 * t5) * t42 + t50 * t37 * t32 * t24 + t2
     #3 * t8 * t13 * t7 * (-t48 * t16 * t5 + t34 * (-t48 + (-t30 * t25 *
     #* 2 - t47 * t25) * t5 * t13)) + t51 * t6 * (-t21 * t46 + t45) * t7
     #) + 64 * t4 * t6 * t2 * (t1 ** 2 * t12 * t27 * t5 * t15 + t54 * (t
     #53 - t47)) - 16 * t3 + 48 * t11 * t16 * t19 * t6 * (t53 - t47) - 1
     #28 * t52 * t50 * t1 * t19 * t15 * (t51 * t12 * t27 + 1) - 80 * t14

      hjetmass_triangle_ppmm_s234_0_mhsq_rat = ret/32q0/(0,1q0)
      return

      end function
