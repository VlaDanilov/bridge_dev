
      complex*32 function hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0
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

      t1 = za(i1, i2)
      t2 = zb(i2, i1)
      t3 = za(i1, i3)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i3, i2)
      t7 = t1 * t2
      t8 = t3 * t4
      t9 = t5 * t6 + t7 + t8
      t10 = za(i1, i4)
      t11 = zb(i4, i1)
      t12 = za(i2, i4)
      t13 = zb(i4, i2)
      t14 = za(i3, i4)
      t15 = zb(i4, i3)
      t16 = t10 * t11
      t17 = t12 * t13
      t18 = t14 * t15
      t19 = t16 + t17 + t18
      if ( real(t19) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t18 = cg * sqrt(t19 ** 2) + t16 + t17 + t18
      t19 = t1 * t13 + t15 * t3
      t20 = 0.1q1 / t18
      t8 = -2 * t16 * t9 * t20 + t7 + t8
      t21 = -2 * t10 * t9
      t22 = t10 * t8
      t23 = t11 * (t12 * (t21 * t13 * t20 + t3 * t6) + t22)
      t6 = t21 * t15 * t20 - t1 * t6
      t21 = t5 * t4
      t24 = -2 * t12 * t9 * t11 * t20 + t21
      t25 = t1 * t11 - t15 * t5
      t26 = t10 * (t11 * t8 + t13 * t24)
      t27 = -2 * t14 * t9 * t11 * t20 - t2 * t5
      t26 = 0.1q1 / t26
      t28 = 0.1q1 / t1
      t2 = 0.1q1 / t2
      t23 = 0.1q1 / t23
      t29 = 0.1q1 / t6
      t27 = 0.1q1 / t27
      t30 = 0.1q1 / t5
      t14 = 0.1q1 / t14
      t31 = t4 * t19
      t32 = t11 * t6
      t33 = t15 * t8
      t34 = t33 - t32 - t31
      t35 = t11 * t30
      t36 = t15 * t28
      t37 = -t35 + t36
      t38 = t23 ** 2
      t39 = t23 * t38
      t40 = t11 ** 2
      t41 = t11 * t40
      t42 = t15 ** 2
      t43 = t12 * t19
      t44 = t33 * t29
      t45 = t43 * t28
      t46 = t40 * t12
      t47 = t4 * t15
      t48 = t11 * t2
      t49 = 0.1q1 / t9
      t18 = 0.1q1 / t18
      t50 = t12 ** 2
      t51 = mt ** 2
      t52 = t30 * t2
      t53 = t52 * t43 * t4
      t54 = t19 * t6
      t55 = t2 * t49
      t56 = t11 * t23
      t57 = t40 * t6
      t3 = t11 * (t9 * (t11 * (-t53 * t23 * t18 + t8 * (t4 * t13 * t2 + 
     #t15) * t30 * t18 * t38 * t19 * t50) + t40 * (t54 * t38 * t18 * t30
     # * t50 + t22 * t53 * t38 * t18)) + t55 * t47) + t52 * (t25 * (t11 
     #* (t32 * t23 - t14) + t33 * (-t14 * t29 + t56)) * t12 + t8 * (t57 
     #* t38 * (t11 * t19 - t13 * t25) + t15 * (t11 * t3 + t13 * t5) * t2
     #9 * t14 ** 2) * t50 + t11 * t25 * (t10 * t15 * t26 - t27) * t24) *
     # t20 * t28 * t51
      t5 = t32 + t31
      t5 = t43 * t5 * t30 - t5 * t6
      t52 = -t43 * t30 + t6
      t53 = t29 ** 2
      t58 = t18 * t2
      t59 = t32 * t8 * t38
      t60 = t30 * t18
      t10 = t60 * t23 * t19 * t50 * t40 * t9 * (t59 * (t17 + t16) + t58 
     #* (t10 * (t57 * t23 + t56 * t33) - t15 - t10 ** 2 * t41 * t8 * t6 
     #* t38 + t17 * t23 * (t33 + t32) - t59 * t50 * t13 ** 2) * t28 * t9
     #)
      t1 = -t60 * t44 * t48 * t9 * t12 * t14 + t1 * t4 * (t43 * t8 * t38
     # + t55) * t30 * t40 + t21 * t55 * t42 * t28 - t7 * t54 * t50 * t41
     # * t8 * t30 * t39
      t7 = -4
      ret = t7 * (-t48 * t37 * t27 * t24 + t48 * (-t46 * t6 * t30 + t32 
     #* (t36 * t12 - t4) + t47 * (t45 + t8)) * t23 + t46 * (t43 * t34 * 
     #t30 - t34 * t6) * t38 + t2 * (-t12 * t42 * t8 * t28 * t29 + t44 * 
     #(t35 * t12 - t4) + t4 * t11 * (-t43 * t29 * t30 + 1)) * t14 + t16 
     #* t15 * t24 * t2 * t37 * t26) - 32 * t3 - 8 * t58 * t12 * t9 * (t3
     #0 * (-t11 * (t31 * t29 + t11) + t31 * t33 * t53 - t42 * t8 ** 2 * 
     #t53 + t45 * t29 * t15 * (t29 * (-t33 + t31) + t11)) * t14 + t56 * 
     #t36 * (-t33 + t31 + t32) + t28 * t40 * (t16 * t5 + t17 * t5 + t33 
     #* (t16 * t52 + t17 * t52)) * t38) + 64 * t10 + 16 * t1 - 128 * t2 
     #* (t22 * t54 * t9 ** 2 * t18 ** 2 * t12 * t50 * t40 ** 2 * t13 * t
     #28 * t30 * t39 + (-t36 + t35) * t20 * t49 * t25 * t4 * t51)

      hjetmass_qqbgg_triangle_pmpp_s123_mhsq_0 = ret/32q0/(0,1q0)
      return

      end function
