
      complex*32 function hjetmass_qqbgg_triangle_pmpm_s12_s123_0
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
      t2 = zb(i3, i2)
      t3 = za(i3, i4)
      t4 = za(i1, i2)
      t5 = zb(i4, i2)
      t6 = za(i1, i3)
      t7 = zb(i2, i1)
      t8 = zb(i4, i3)
      t9 = zb(i3, i1)
      t10 = za(i2, i3)
      t11 = t6 * t9
      t12 = t10 * t2
      t13 = t12 + t11
      if ( real(t13) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t13 = cg * sqrt(t13 ** 2) + t11 + t12
      t14 = 0.1q1 / t13
      t15 = 2 * t6
      t16 = t15 * t7
      t17 = t4 * (t16 * t8 * t14 - t5)
      t18 = t4 * t7
      t15 = t18 * (-t15 * t9 * t14 + 1)
      t16 = t9 * (t15 * t6 - t17 * t3) - t12 * t16 * t4 * t9 * t14
      t19 = za(i1, i4)
      t20 = t1 * t2
      t21 = t19 * t9
      t22 = zb(i4, i1)
      t17 = 0.1q1 / t17
      t16 = 0.1q1 / t16
      t23 = 0.1q1 / t8
      t13 = 0.1q1 / t13
      t24 = t12 + t11
      t25 = t16 ** 2
      t26 = t16 * t25
      t27 = t17 ** 2
      t28 = t3 ** 2
      t29 = t2 ** 2
      t30 = t4 ** 2
      t31 = t4 * t30
      t32 = t9 ** 2
      t33 = t9 * t32
      t34 = t19 * t22
      t35 = t1 * t5
      t36 = t18 * t13
      t37 = t23 * t15
      t38 = t37 * t33
      t39 = t36 * t33
      t40 = t33 * t15
      t41 = 0.1q1 / t3
      t42 = 0.1q1 / t7
      t43 = t34 + t35
      t44 = mt ** 2
      t45 = t3 * t7
      t46 = t3 * t15 * t16
      t47 = t32 * t16
      t48 = t41 * t17
      t49 = t12 * t23 + t3
      t50 = t1 ** 2
      t51 = t10 * t19
      t52 = t1 * t6
      t53 = t15 * t42
      t54 = t11 * t23
      t55 = t9 * t3
      t56 = t9 * t16
      t57 = t1 * t9
      t5 = t2 * (t45 * t38 * t12 * t31 * t26 + t47 * t13 * (t1 * t29 * t
     #7 * t10 ** 2 * t16 * t23 + t12 * t16 * (-t55 * t37 + (t54 + t3) * 
     #t7 * t1) + t55 * ((-t54 - t3) * t16 * t15 + t23)) * t30 + t57 * t4
     #2 * t23 * (t56 + t48) * t44 + t10 * (t32 * (t20 * t15 * t13 * t49 
     #* t25 - t20 * t23 * t13 * t16) + t33 * (-t19 * t23 * t13 * t16 + t
     #13 * t15 * (t19 * t3 + (t52 + t51) * t23 * t2) * t25 + t37 * t2 * 
     #t3 * t42 * (t19 ** 2 * t22 ** 2 + t5 ** 2 * t50) * t26) + t37 * t6
     # * t32 ** 2 * t19 * t25 * t13 + t53 * t2 * t8 * t17 * t27) * t4)
      t7 = t53 + t4
      t19 = t43 * t42
      t1 = t47 * t23 * t2 * (t16 * (t21 * t15 * (t19 + t4) + t20 * (t4 *
     # (t18 + t15) + t34 * t7 + t35 * t7)) * t10 + t46 * t9 * t4 * (-t19
     # - t4) + t11 * t4 * t1 * t14)
      t7 = t2 * (t17 * (t4 * t32 * t42 * t23 + t41 * t50) + t56 * (t3 * 
     #t4 * t32 * t42 * t23 + t50 + t57 * (t41 * (t52 - t51) + t53) * t23
     #))
      ret = -64 * t13 * t10 * t30 * t29 * (t27 * (-t36 * t9 + (-t36 - 1)
     # * t17 * t15 * t8) + t28 * (t39 * t25 + t40 * (t18 + t35 + t34) * 
     #t26) + t3 * (t39 * t23 * t24 * t25 + t38 * (t34 * t24 + t35 * t24 
     #+ t18 * (-t6 ** 2 * t32 * t13 + t11 + t12 * (-t12 * t13 + 1))) * t
     #26) - t36 * t3 * t28 * t33 * t15 * t8 * t26) - 32 * t12 * (t2 * (-
     #t35 * t34 * t38 * t3 * t26 * t42 * t4 + t9 * (-t13 * t27 + (-t13 *
     # t43 * t23 * t25 - t37 * t43 * t26) * t32 * t3) * t30 - t45 * t33 
     #* t25 * t23 * t13 * t31) + t42 * t14 * (t20 + t21) * t44 * (t15 * 
     #t41 * t27 + t47 * (-t46 + t23) + t48 * t9 * t23)) + 128 * t45 * t4
     #0 * t13 ** 2 * t26 * t10 * t31 * t29 * (t11 * t49 + t12 * t3) + 16
     # * t5 - 8 * t1 - 4 * t7

      hjetmass_qqbgg_triangle_pmpm_s12_s123_0 = ret/32q0/(0,1q0)
      return

      end function
