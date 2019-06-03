
      complex*32 function hjetmass_triangle_ppmm_s14_s124_0_rat
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
      t4 = za(i2, i4)
      t5 = zb(i4, i2)
      t6 = zb(i4, i3)
      t7 = za(i1, i3)
      t8 = zb(i3, i1)
      t9 = za(i2, i3)
      t10 = zb(i3, i2)
      t11 = t7 * t8
      t12 = t9 * t10
      t13 = t11 + t12
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t13 = cg * sqrt(t13 ** 2) + t11 + t12
      t14 = (0.1q1 / 0.2q1)
      t15 = t14 * t13
      t16 = t7 * t2
      t17 = t3 * (-t15 * t5 + t16 * t6)
      t18 = za(i1, i4)
      t19 = 0.1q1 / t13
      t20 = 2 * t3
      t21 = t20 * t1
      t22 = t2 * (t21 * t10 * t19 + t18)
      t23 = t1 * t6
      t24 = t3 * t8
      t25 = t24 * t16 * (t23 + t11 + t12)
      t26 = -t14 * t24 * t13 * (t1 * t5 + t16) + t25
      t27 = t2 * (t24 * t1 - t15 * t4)
      t28 = t18 * zb(i4, i1) + t4 * t5
      t24 = t14 * t16 * t13 * (t4 * t6 + t24) - t25
      t25 = -t27
      t29 = -t2 * (t3 * t1 * t10 + t15 * t18)
      t30 = t2 * (t21 * t8 * t19 - t4)
      t21 = -t21 * t2 * t6 * t19 + t28
      t31 = t3 * t2
      t16 = -t12 * t20 * t16 * t8 * t19 + t7 * (-t30 * t6 + t8 * t31 * (
     #-2 * t11 * t19 + 1))
      t15 = t31 * (-t11 + t15)
      t20 = -t29
      t31 = 0.1q1 / t3
      t6 = 0.1q1 / t6
      t32 = -0.1q1 / t24
      t33 = 0.1q1 / t2
      t4 = 0.1q1 / t4
      t34 = 0.1q1 / t26
      t17 = 0.1q1 / t17
      t5 = 0.1q1 / t5
      t35 = 0.1q1 / t8
      t18 = 0.1q1 / t18
      t36 = t2 ** 2
      t37 = t31 ** 2
      t38 = t7 * t25 * t32
      t39 = t30 * t35
      t40 = t4 * t5
      t26 = -0.1q1 / t26
      t41 = 0.1q1 / t30
      t24 = 0.1q1 / t24
      t42 = 0.1q1 / t9
      t43 = t8 * t4
      t44 = t10 * t18
      t45 = t44 + t43
      t46 = t1 * t8
      t26 = t46 * t26 + t17
      t47 = t7 ** 2
      t48 = t7 * t47
      t49 = t24 ** 2
      t50 = t25 ** 2
      t51 = t6 ** 2
      t52 = t32 ** 2
      t53 = t12 * t31 * t19
      t54 = t53 * t47
      t24 = t7 * t24
      t55 = 0.1q1 / t27
      t55 = t55 * t6
      t16 = 0.1q1 / t16 ** 2
      t8 = t40 * t1 * (t47 * t8 * t26 * t19 * t36 + t31 * (-t7 * t15 * t
     #26 + t12 * (-t24 * t27 - t6)) * t19 * t2 + t44 * t31 * (t6 * (-t21
     # * t6 - t39) + t47 * t30 ** 2 * t21 * t16))
      t16 = t46 * t40 * t9 * t2 * t20 * t19 * t31 * (t24 + t55)
      ret = 8 * t40 * (t1 * (t7 * t17 * t36 + t31 * (t38 + t6) * t2) - t
     #11 * t1 ** 2 * t36 * t34 - t39 * t13 * t22 * t37 * t18 * t33 * t6)
     # + 32 * t5 * ((-t24 * t20 * t19 * t37 * t33 * t4 * t35 * t18 + t53
     # * t2 * t48 * t45 * t49) * t27 ** 2 + t1 * (-t44 * t47 * t26 * t19
     # * t36 + (t9 * (t21 * t31 * t51 * t41 * t18 * t10 ** 2 + t43 * t21
     # * t31 * t51 * t41 * t10) + t1 * t7 * t15 * t4 * t18 * t42 * t26) 
     #* t19 * t2 + t22 * t31 * t18 * t4 * t6) + t50 * (-t44 * t1 * t7 * 
     #t19 * t31 * t4 * t35 * t32 - t54 * t29 * t4 * t18 * t52) + t54 * t
     #2 * (t45 * t20 * t9 - t1 * (t14 * t13 * t28 - t23 * t3 * t2) * t45
     #) * t49 * t27 - t44 * t48 * t25 * t50 * t19 * t31 * t4 * t52) - 64
     # * t18 * t5 * t31 * t19 * t1 * (t38 * t29 * t4 + t12 * (-t55 - t24
     #) * t20 * t2) - 16 * t8 + 48 * t16

      hjetmass_triangle_ppmm_s14_s124_0_rat = ret/32q0/(0,1q0)
      return

      end function
