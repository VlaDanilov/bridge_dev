
      complex*32 function hjetmass_qqbgg_triangle_pmmp_s12_s123_0
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

      t1 = za(i2, i3)
      t2 = zb(i2, i1)
      t3 = za(i1, i2)
      t4 = za(i3, i4)
      t5 = za(i2, i4)
      t6 = zb(i3, i1)
      t7 = za(i1, i3)
      t8 = zb(i3, i2)
      t9 = t7 * t6
      t10 = t1 * t8
      t11 = t9 + t10
      if ( real(t11) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t11 = cg * sqrt(t11 ** 2) + t10 + t9
      t12 = 0.1q1 / t11
      t13 = 2 * t3
      t14 = t2 * (t13 * t4 * t6 * t12 - t5)
      t15 = za(i1, i4)
      t16 = zb(i4, i1)
      t17 = zb(i4, i2)
      t18 = zb(i4, i3)
      t19 = t3 * t2
      t20 = t19 * (-2 * t9 * t12 + 1)
      t13 = t7 * (-t14 * t18 + t20 * t6) - t10 * t9 * t13 * t2 * t12
      t21 = t1 * t17 + t16 * t7
      t22 = 0.1q1 / t18
      t14 = 0.1q1 / t14
      t11 = 0.1q1 / t11
      t23 = 0.1q1 / t4
      t13 = 0.1q1 / t13
      t24 = t5 * t17
      t25 = t15 * t16
      t26 = t18 ** 2
      t27 = t1 ** 2
      t28 = t1 * t27
      t29 = t13 ** 2
      t30 = t13 * t29
      t31 = t14 ** 2
      t32 = t2 ** 2
      t33 = t7 ** 2
      t34 = t33 ** 2
      t35 = t7 * t33
      t36 = t11 ** 2
      t37 = mt ** 2
      t38 = t12 * (t25 + t24)
      t39 = t2 * t18
      t40 = t33 * t18 * t29
      t41 = t6 * t2
      t42 = t41 * t27
      t43 = t19 + t24 + t25
      t44 = t6 ** 2
      t45 = t19 * t11
      t46 = t43 * t30
      t47 = t10 * t23
      t48 = t12 * t11
      t49 = t48 * t44 * t3
      t50 = 0.1q1 / t3
      t51 = t6 * t17
      t52 = t39 + t51
      t53 = t52 * t12
      t54 = t16 * t18
      t55 = t37 * t16 * t50
      t56 = t23 * t16
      t41 = t27 * (t33 * (-t41 * t16 * t11 * t23 * t13 + t41 * (-t54 * t
     #45 + (t16 * (t53 * t15 - t45 * t8) + t53 * t24 + t53 * t19) * t23 
     #* t1) * t29) + t55 * t23 * t14 * t22 + t23 * (-t1 * t2 * t11 * t52
     # + t55) * t13 * t7 + t56 * t44 * t2 * (t19 * (-t11 + t12) + t38) *
     # t29 * t35)
      t52 = t16 ** 2
      t16 = t8 * t16
      t55 = t3 ** 2
      t57 = t37 * t21 * t12 * t50 * t23
      t53 = t53 * t3
      t5 = t6 * t27 * (t14 * (-t57 * t22 + (-t4 * t6 * t12 * t31 + t11 *
     # t14) * t32 * t1) + t35 * (t48 * t6 * t3 * t32 * ((-t39 - t51 - t1
     #6) * t23 * t1 - t54) * t29 + t23 * t12 * t18 * t6 * t32 * t1 * (-t
     #5 ** 2 * t17 ** 2 - t52 * t15 ** 2 - t32 * t55) * t30) - t56 * t49
     # * t32 * t34 * t29 - t57 * t7 * t13 + t11 * t32 * t1 * (t18 * (t23
     # * t43 - t53) - t53 * t47) * t29 * t33)
      t9 = t12 * t30 * t36 * t18 * t35 * t44 * t55 * t32 ** 2 * t28 * (t
     #9 * t18 + t10 * (t9 * t23 + t18))
      ret = -64 * t42 * (t1 * (t39 * t25 * t24 * t6 * t35 * t12 * t23 * 
     #t30 + t3 * (-t31 * t36 + (t6 * t36 * t23 * t29 + t38 * t6 * t23 * 
     #t30) * t18 * t35 + t33 * t26 * t36 * t29) * t32) + t37 * t6 * t21 
     #* t12 ** 2 * (-t22 * t31 + t40) + t40 * t3 * t32 * t36 * t27 * t8 
     #* t23) + 128 * t49 * t2 * t32 * t28 * (t35 * (-t45 * t4 * t30 * t1
     #8 * t26 + t46 * t26 + t47 * (t19 * (-t10 * t11 + 1) + t24 + t25) *
     # t30 * t18) + t46 * t23 * t6 * t18 * t34 - t45 * t7 * t34 * t44 * 
     #t18 * t23 * t30 - t4 * t14 * t31 * (t45 + 1)) + 16 * t41 - 4 * t1 
     #* (-t14 * (t2 * t27 * t50 * t23 + t22 * t52) + (-t39 * t27 * t50 *
     # t23 - t52 + t56 * (t20 * t50 + t22 * (-t16 + t51)) * t1) * t13 * 
     #t7) + 32 * t5 - 256 * t9 + 8 * t56 * t42 * t13 * t33 * (t13 * t43 
     #- t12)

      hjetmass_qqbgg_triangle_pmmp_s12_s123_0 = ret/32q0/(0,1q0)
      return

      end function
