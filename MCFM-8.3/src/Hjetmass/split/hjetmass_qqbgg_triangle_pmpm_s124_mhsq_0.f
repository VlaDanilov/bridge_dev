
      complex*32 function hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i1, i2)
      t4 = zb(i3, i2)
      t5 = za(i1, i4)
      t6 = zb(i4, i3)
      t7 = t3 * t4 - t5 * t6
      t8 = za(i1, i3)
      t9 = zb(i2, i1)
      t10 = zb(i4, i1)
      t11 = zb(i4, i2)
      t12 = t3 * t9
      t13 = t5 * t10
      t14 = t1 * t11 + t12 + t13
      t15 = za(i2, i3)
      t16 = za(i3, i4)
      t17 = t8 * t2
      t18 = t15 * t4
      t19 = t16 * t6
      t20 = t17 + t18 + t19
      if ( real(t20) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t19 = cg * sqrt(t20 ** 2) + t17 + t18 + t19
      t20 = 0.1q1 / t19
      t13 = -2 * t17 * t14 * t20 + t12 + t13
      t21 = 2 * t8
      t22 = -t21 * t4 * t14 * t20 + t11 * t5
      t23 = t8 * t13
      t24 = t2 * (t15 * t22 + t23)
      t11 = t21 * t14 * t6 * t20 - t11 * t3
      t21 = t1 * t6
      t25 = t3 * t2 + t21
      t26 = -2 * t15 * t2 * t14 * t20 + t1 * t10
      t27 = 2 * t16 * t2 * t14 * t20 - t1 * t9
      t28 = t2 * t13
      t29 = t8 * (t26 * t4 + t28)
      t5 = t5 * t2
      t30 = t1 * t4
      t31 = t30 + t5
      t32 = 0.1q1 / t6
      t10 = 0.1q1 / t10
      t33 = 0.1q1 / t9
      t11 = 0.1q1 / t11
      t3 = 0.1q1 / t3
      t34 = 0.1q1 / t15
      t19 = 0.1q1 / t19
      t27 = 0.1q1 / t27
      t35 = 0.1q1 / t16
      t36 = 0.1q1 / t14
      t24 = 0.1q1 / t24
      t29 = 0.1q1 / t29
      t37 = t2 * t24
      t38 = t2 ** 2
      t39 = t2 * t38
      t40 = t13 ** 2
      t41 = t24 ** 2
      t42 = t24 * t41
      t43 = mt ** 2
      t44 = t13 * (t11 * t35 + t37)
      t45 = t2 * t32
      t46 = t21 * t3
      t47 = t3 * t33
      t48 = t47 * t10
      t49 = t1 * t38
      t50 = t33 * t11
      t51 = t3 * t1
      t52 = t8 * t29
      t53 = t38 * t33 * t10
      t54 = t19 ** 2
      t55 = t15 ** 2
      t56 = t11 ** 2
      t57 = t14 ** 2
      t58 = t6 ** 2 * t40 * t56
      t59 = t39 * t16
      t60 = t59 * t4 * t40 * t10
      t61 = t7 * (t14 * (-t60 * t19 * t42 * t55 + t59 * t13 * t41 * (-t1
     #7 * t13 * t24 + 1) * t19 * t10 * t15) + t57 * (t59 * t48 * t4 ** 2
     # * t40 * t54 * t42 * t15 * t55 + t48 * (t11 * (t58 + t38) + t16 * 
     #(t2 * t38 ** 2 * t8 ** 2 * t40 * t42 + t24 * t39)) * t54 * t15)) +
     # t47 * t43 * t1 ** 2 * t2 * t34 * t36
      t62 = t2 * t22
      t63 = t4 * t13
      t64 = t63 - t62
      t23 = t23 * t39
      t24 = t38 * t24
      t39 = t63 * t38
      t32 = t48 * t1 * ((t15 * (t23 * t64 * t41 + t24 * (-t63 + t62)) + 
     #t38 + t58 + t39 * t64 * t41 * t55) * t19 * t14 + t43 * t2 * (-t2 *
     # (t27 * t32 + t52) * t34 * t26 - t44))
      t56 = t6 * t56
      t4 = t48 * t20 * t40 * t15 * t43 * (t56 * t31 * t35 + t16 * t38 * 
     #t41 * (-t2 * t7 + t25 * t4))
      t3 = t3 * t2 * ((t55 * (-t63 * t53 * t16 * t54 * t41 + t60 * t8 * 
     #t33 * t54 * t42) + t10 * t33 * (-t59 * t8 * t41 + t56) * t54 * t13
     # * t15) * t57 * t7 + t43 * t1 * t25 * t20 * t36 * (-t16 * t34 * t3
     #3 + t10))
      t6 = t51 * t19 * t14 * t7 * (t39 * t55 * t10 * t41 - t50 * (t13 * 
     #t6 * t11 + t2) + (t23 * t41 - t24) * t10 * t15)
      t8 = 4
      ret = 16 * t2 * (-t12 * t15 * t16 * t38 * t40 * t7 * t10 * t42 + t
     #48 * (t44 * t31 * t15 + t45 * t25 * (-t16 * t26 * t34 * t27 + 1) +
     # t17 * t26 * (-t16 * t25 * t34 + t31) * t29) * t20 * t43 - t21 * t
     #48 * t14 * t19 * t13 * t11 + t1 * (t10 * (t46 + t2) + t47 * (t30 +
     # t5)) * t36) + t8 * (t18 * t1 * t38 * t40 * t10 * t41 - t52 * t48 
     #* t38 * t31 * t26 - t53 * (t45 + t51) + t51 * (t11 * (-t1 * t35 + 
     #t2 * t33) + t37 * (t15 * t2 * t10 - t1)) * t7 + t28 * t10 * (t15 *
     # (-t49 * t22 * t41 - t37 * t47 * t31) + t50 * (t46 + t2))) - 64 * 
     #t61 - 8 * t32 - 48 * t28 * t48 * t43 * t25 * t20 * (t37 * t16 + t1
     #1) + 32 * t4 - 128 * t3 - 24 * t6 + 12 * t49 * t15 * t9 * t13 * t7
     # * t10 * t41

      hjetmass_qqbgg_triangle_pmpm_s124_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function
