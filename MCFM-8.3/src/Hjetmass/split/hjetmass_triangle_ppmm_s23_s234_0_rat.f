
      complex*32 function hjetmass_triangle_ppmm_s23_s234_0_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = za(i1, i4)
      t2 = zb(i2, i1)
      t3 = za(i2, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i1)
      t7 = za(i1, i3)
      t8 = t7 * t4
      t9 = t1 * t6
      t10 = t8 + t9
      if ( real(t10) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t8 = cg * sqrt(t10 ** 2) + t8 + t9
      t9 = za(i1, i2)
      t10 = zb(i3, i2)
      t11 = zb(i4, i3)
      t12 = za(i3, i4)
      t13 = 0.1q1 / t8
      t14 = t12 * (-2 * t1 * t2 * t11 * t13 - t10)
      t15 = zb(i4, i2)
      t16 = t3 * t10 + t15 * t5
      t17 = -2 * t9 * t12
      t18 = t17 * t2 * t11 * t13 + t16
      t19 = t11 * (t17 * t6 * t13 - t3)
      t17 = (t1 * t19 + t18 * t9 + t7 * t11 * (t17 * t4 * t13 + t5)) * t
     #2
      t20 = (0.1q1 / 0.2q1)
      t21 = t20 * t8
      t22 = t7 * t2
      t23 = t22 * t11
      t24 = t12 * (t21 * t15 - t23)
      t25 = t1 * t2
      t26 = -t12 * (t21 * t10 + t25 * t11)
      t27 = t9 * t12
      t8 = -t27 * t2 * t11 + t20 * t8 * t16
      t16 = t12 * (-2 * t23 * t13 + t15)
      t20 = t2 * t18
      t23 = t4 * t16
      t28 = t6 * t14
      t29 = (t28 + t20 + t23) * t9
      t21 = -t11 * (t21 * t3 + t27 * t6)
      t27 = 0.1q1 / t18
      t19 = 0.1q1 / t19
      t15 = 0.1q1 / t15
      t30 = 0.1q1 / t6
      t31 = 0.1q1 / t9
      t32 = 0.1q1 / t5
      t21 = 0.1q1 / t21
      t29 = 0.1q1 / t29
      t11 = 0.1q1 / t11
      t17 = 0.1q1 / t17
      t10 = 0.1q1 / t10
      t33 = t18 * t19
      t34 = t33 * t10
      t35 = t34 + t11
      t36 = t2 ** 2
      t37 = t18 ** 2
      t38 = t17 ** 2
      t39 = 0.1q1 / t8 ** 2
      t40 = t7 * t6
      t41 = t14 * t11
      t42 = t2 * t14
      t3 = t13 * t32 * (t3 * t4 + t5 * t6)
      t4 = 0.1q1 / t1
      t4 = (-t9 * t2 * t29 + t27) * t4
      t5 = t3 * t10 * t2 * (-t31 * (t7 * t26 * t15 * t21 + t41 * t1 * t2
     #7) + t42 * (t33 * t7 * t15 + t1 * t11) * t17)
      ret = -32 * t3 * (t10 * t11 * t2 * (t31 ** 2 + t20 * t38 * (-t28 -
     # t20 - t23)) * t1 ** 2 - t42 * t10 * (t40 * t37 * t15 * t19 ** 2 *
     # t17 + t41 * t27 ** 2) + t15 * t10 * (t40 * t8 * t31 * t21 ** 2 - 
     #t36 * t24 * t39 * t30) * t26 + t25 * t15 * (t2 * (t34 * t16 * t17 
     #+ t18 * t7 * (-t11 * (t28 + t23) + t34 * (-t28 - t23)) * t38) + t3
     #1 * (-t24 * t21 * t10 + t35 * t31 * t7) - t7 * t37 * t35 * t38 * t
     #36) + t36 * t14 * t27 * t10 * (t2 * t16 * t15 * t30 + t41) * t29 *
     # t9) + 8 * t11 * t32 * t15 * t12 * t2 * (t4 * t14 - t20 * t17 + t3
     #1) - 48 * t41 * t3 * t22 * t15 * (t17 * t2 - t27 * t31) + 16 * t3 
     #* t11 * t15 * t2 * (t16 * (t27 * (-t42 * t9 * t29 + t1 * t31) - t2
     #5 * t17) + t24 * t26 * t39 + t4 * t14 ** 2 * t27 * t7) - 64 * t5

      hjetmass_triangle_ppmm_s23_s234_0_rat = ret/32q0/(0,1q0)
      return

      end function
