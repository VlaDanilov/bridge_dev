
      complex*32 function hjetmass_triangle_pppm_s124_0_mhsq_rat
     &     (i1,i2,i3,i4,za,zb)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          real*16 cg

      t1 = zb(i3, i1)
      t2 = zb(i3, i2)
      t3 = za(i1, i4)
      t4 = za(i2, i4)
      t5 = t1 * t3 + t2 * t4
      t6 = za(i1, i3)
      t7 = za(i2, i3)
      t8 = za(i3, i4)
      t9 = zb(i4, i3)
      t10 = t6 * t1
      t11 = t7 * t2
      t12 = t8 * t9
      t13 = t10 + t11 + t12
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i1, i2)
      t13 = zb(i4, i1)
      t14 = zb(i2, i1)
      t15 = zb(i4, i2)
      t16 = t12 * t14
      t17 = t3 * t13
      t18 = t15 * t4 + t16 + t17
      t11 = 0.1q1 / t11
      t10 = -2 * t10 * t18 * t11 + t16 + t17
      t16 = 2 * t6 * t18
      t17 = t16 * t9 * t11 - t12 * t15
      t15 = -t16 * t2 * t11 + t15 * t3
      t16 = 2 * t8 * t1 * t18 * t11 - t14 * t4
      t19 = -2 * t7 * t1 * t18 * t11 + t13 * t4
      t20 = t1 * t10
      t21 = t6 * (t19 * t2 + t20)
      t22 = t6 * t10
      t23 = t1 * (t15 * t7 + t22)
      t24 = t12 * t2 - t3 * t9
      t25 = t1 * t12 + t4 * t9
      t22 = t1 * (-t17 * t8 + t22)
      t20 = t6 * (-t16 * t9 + t20)
      t26 = 0.1q1 / t9
      t23 = 0.1q1 / t23
      t19 = 0.1q1 / t19
      t12 = 0.1q1 / t12
      t13 = 0.1q1 / t13
      t20 = 0.1q1 / t20
      t7 = 0.1q1 / t7
      t22 = 0.1q1 / t22
      t21 = 0.1q1 / t21
      t3 = 0.1q1 / t3
      t27 = 0.1q1 / t4
      t28 = 0.1q1 / t15
      t17 = 0.1q1 / t17
      t29 = t2 * t3
      t30 = t1 * t27 + t29
      t31 = t9 * t15 * t17
      t32 = t31 + t2
      t33 = t1 ** 2
      t34 = t1 * t33
      t35 = t8 ** 2
      t36 = t8 * t35
      t37 = t22 ** 2
      t38 = t23 ** 2
      t39 = t29 * t25
      t40 = t38 * t12
      t41 = t1 * t2
      t42 = t5 * t1
      t43 = t3 * t12
      t11 = t13 * t11
      t13 = 0.1q1 / t18
      t4 = t11 * t14 * t5 * t12 * t13 * (t29 * t4 + t1)
      ret = 32 * t11 * (t6 * (-t33 * t5 * t19 * t21 * t3 * t27 * t16 ** 
     #2 + t42 * (t41 * t26 * t20 * t3 * t27 + t30 * t21 * t12) * t16) + 
     #t43 * (-t33 * t2 * t35 * t15 * t25 * t38 + t34 * t35 * t15 * t24 *
     # t38 - t41 * t5 * t8 * t23 - t5 * t17 * t32) * t10 - t43 * t42 * t
     #15 * (t8 * t1 * t23 + t17) + t27 * (t1 * (t5 * t22 * t3 * t17 * t3
     #2 * t8 - t29 * t5 * t23 * t7 * t28 * t35 - t39 * t23 * t7 ** 2 * t
     #28 * t36) + t33 * (t35 * (-t31 * t5 * t37 * t3 - t40 * t2 * t25) -
     # t39 * t38 * t7 * t36) + t34 * (t24 * t38 * t7 * t3 * t36 + t24 * 
     #(-t15 * t37 * t3 * t17 + t40) * t35) - t5 * t9 * t12 * t17 ** 2) *
     # t10 ** 2 + t42 * t12 * t30 * t26) - 64 * t11 * t42 * t27 * t10 * 
     #(t1 * (t8 * (-t15 * t22 * t3 * t17 + t12 * t23) + t23 * t7 * t3 * 
     #t35) + t12 * t17) - 128 * t4

      hjetmass_triangle_pppm_s124_0_mhsq_rat = ret/32q0/(0,1q0)
      return

      end function
