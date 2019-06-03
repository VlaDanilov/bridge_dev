
      complex*32 function hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_rat
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

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i3)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = t3 * t2
      t10 = t4 * t5
      t11 = t6 * t7
      t12 = t1 * t8
      t13 = t9 + t10 + t11 + t12
      t14 = za(i1, i2)
      t15 = za(i3, i4)
      t16 = zb(i2, i1)
      t17 = zb(i4, i3)
      t18 = t14 * t16
      t13 = -0.4q1 * t18 * t15 * t17 + t13 ** 2
      t13 = cg * sqrt(t13) + t10 + t11 + t12 + t9
      t19 = 0.1q1 / 0.2q1
      t20 = t19 * t13
      t21 = t6 * t16 * t17 + t20 * t5
      t22 = -t14 * t21
      t23 = t1 * t16 * t17 - t20 * t2
      t24 = t15 * t23
      t25 = t11 + t9
      t26 = t13 ** 2
      t27 = t25 * t16
      t28 = t15 * t5
      t29 = t26 * t14
      t30 = t4 * t2
      t31 = t1 * t7
      t32 = t31 + t30
      t33 = t3 ** 2 * t2 ** 2
      t34 = t6 ** 2 * t7 ** 2
      t35 = t5 * t32 * t3
      t32 = t8 * t32 * t6
      t36 = t18 * t13
      t37 = t17 ** 2
      t25 = t25 * t14
      t38 = t25 * t17
      t4 = t4 * t6
      t39 = t11 * t9 * t13
      t40 = 0.1q1 / 0.4q1
      t41 = t19 * t36 * ((t11 + t18 + t10) * t17 * t15 + t33 + t34 + t35
     # + t32)
      t42 = -t40 * t29 * (-t28 * t7 + t27) - t18 * ((-t4 * t37 + t38) * 
     #t16 * t15 - t39) + t41
      t43 = t15 * t2
      t3 = t1 * t3
      t32 = t19 * t36 * ((t12 + t18 + t9) * t17 * t15 + t33 + t34 + t35 
     #+ t32)
      t29 = -t40 * t29 * (t43 * t8 + t27) - t18 * ((t3 * t37 + t38) * t1
     #6 * t15 - t39) + t32
      t33 = -t20 * t1 + t43 * t14
      t34 = t16 * t33
      t28 = t28 * t14 + t20 * t6
      t35 = t17 * t28
      t36 = t26 * t16
      t37 = t15 ** 2
      t27 = t27 * t17 * t15
      t4 = t40 * t36 * (t14 * (-t11 - t9) + t4 * t17) + t41 + t18 * (t14
     # * (t5 * t7 * t17 * t37 - t27) + t39)
      t3 = -t40 * t36 * (t3 * t17 + t25) + t32 - t18 * (t14 * (t2 * t8 *
     # t17 * t37 + t27) - t39)
      t5 = t1 * t5 + t2 * t6
      t6 = t31 + t30
      t7 = t18 * t6
      t8 = t11 + t9
      t25 = t18 * t15 * t17
      t27 = t19 * t13 * t8 - t25
      t30 = t26 * t40 - t25
      t8 = t20 * t18 - t18 * t8
      t18 = -t17 * t33
      t9 = t9 + t10
      t10 = t15 * t17
      t20 = t10 * t20
      t11 = t11 + t12
      t12 = t10 * t5
      t23 = -t14 * t23
      t31 = 0.1q1 / t42
      t3 = 0.1q1 / t3
      t29 = 0.1q1 / t29
      t17 = 0.1q1 / t17
      t14 = 0.1q1 / t14
      t32 = 0.1q1 / t16
      t4 = 0.1q1 / t4
      t33 = 0.1q1 / t15
      t36 = -t3 + t4
      t37 = t34 * t35
      t38 = t33 * t32 * t14 * t17
      t39 = -t31 + t29
      t30 = 0.1q1 / t30 ** 2
      t40 = t27 * (t3 ** 2 - t4 ** 2)
      t41 = t31 ** 2
      t42 = -t29 ** 2 + t41
      t14 = t32 * t14 * t30
      ret = 0.16q2 * t38 * t2 * t1 * (t37 * t36 + (t31 - t29) * t24 * t2
     #2) - 0.8q1 * t38 * t30 * t13 * (t40 * t7 * t35 ** 2 * t34 ** 2 + t
     #8 * t23 * t24 * t12 * t39 + t37 * (t12 * (t3 * (-t27 * (t19 * t13 
     #* t11 - t25) * t3 - 0.1q1) + t4) + t40 * t18 * t16 * t28) * t7) - 
     #0.2q1 * t14 * t6 * t5 * t13 * t26 * (t5 * t8 * t39 + (t29 * (-t8 *
     # (-t10 * t11 + t20) * t29 - 0.1q1) + t31) * t33 * t17 * t24 * t22)
     # - 0.4q1 * t26 * t14 * (t5 ** 2 * t7 * t27 * t36 + (t42 * t8 * t6 
     #* t24 ** 2 * t22 ** 2 + t27 * t5 * t34 * (-t18 * t3 + t4 * (-t35 *
     # t7 * (-t10 * t9 + t20) * t4 + t18)) + (-t42 * t15 * t21 * t23 - t
     #12 * (t19 * t13 * t9 - t25) * t41) * t8 * t6 * t24 * t22) * t33 * 
     #t17)

      hjetmass_qqbgg_triangle_pmpm_s34_mhsq_s12_rat = ret/32q0*(0,1q0)
      return

      end function
