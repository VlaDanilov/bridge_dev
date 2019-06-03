
      double complex function hjetmass_triangle_pmpm_s23_s234_0_dp 
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision mt
          double precision cg

      t1 = za(i2, i3)
      t2 = zb(i3, i2)
      t3 = za(i2, i4)
      t4 = zb(i4, i2)
      t5 = za(i3, i4)
      t6 = zb(i4, i3)
      t7 = t3 * t4
      t8 = t5 * t6
      t9 = t8 + t7
      if ( dreal(t9) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t9 ** 2) + t7 + t8
      t10 = za(i1, i4)
      t11 = zb(i2, i1)
      t12 = zb(i3, i1)
      t13 = t3 * t11
      t14 = t5 * t12
      t15 = t14 + t13
      t16 = zb(i4, i1)
      t17 = 0.1D1 / t9
      t18 = t2 * t16
      t19 = -2 * t18
      t20 = t1 * (t19 * t3 * t17 + t12)
      t21 = t1 * (t19 * t5 * t17 - t11)
      t22 = za(i1, i2)
      t23 = za(i1, i3)
      t24 = t22 * t11
      t25 = t23 * t12
      t19 = t19 * t10 * t1 * t17 + t24 + t25
      t26 = (t16 * t19 + t20 * t4 + t21 * t6) * t10
      t9 = 0.1D1 / t9
      t19 = 0.1D1 / t19
      t27 = 0.1D1 / t2
      t28 = 0.1D1 / t10
      t29 = 0.1D1 / t1
      t26 = 0.1D1 / t26
      t21 = 0.1D1 / t21
      t30 = 0.1D1 / t16
      t31 = t6 * t11
      t32 = -t18 + t31
      t33 = t14 + t13
      t34 = t8 + t7
      t35 = t26 ** 2
      t36 = t26 * t35
      t37 = t19 ** 2
      t38 = t20 ** 2
      t39 = t10 ** 2
      t40 = mt ** 2
      t41 = t5 ** 2
      t42 = t12 ** 2
      t43 = t1 * t15
      t44 = t43 * t12 * t6
      t45 = t16 * t26
      t46 = t18 * t3
      t47 = t9 * t20
      t48 = t16 * t36
      t49 = t20 * t15
      t50 = t19 * t28
      t51 = t50 * t30
      t52 = t3 * t12
      t53 = t15 * t21
      t54 = t1 * t2
      t55 = t54 + t24 + t25
      t56 = t5 * t21
      t57 = t56 * t9
      t58 = t26 * t9
      t59 = t53 * t20
      t17 = t59 * (t58 * t3 * (t33 * t6 - t46) + t49 * t6 * (t57 * t37 +
     # t48 * (t24 * (-t25 * t27 * t29 - 1) - t25) * t39 + t57 * t55 * t3
     #5 * t16 * t10) + t29 * (t3 * (t3 * (t50 - t45) + t56 * (-t50 + t45
     #) * t20) + (t26 * (-t45 * t39 * t20 + t3) + t20 * t37 * t30 - t51 
     #* t3) * t27 * t15 * t6) * t17 * t40)
      t50 = t15 ** 2
      t60 = t16 ** 2
      t61 = t20 * t10 * t39
      t62 = t54 * t21
      t63 = -t54 - t24 - t25
      t64 = t3 ** 2
      t65 = t54 * t41
      t66 = t10 * t20
      t67 = t64 * t28
      t68 = t21 ** 2
      t54 = t47 * t21 * t50 * t6 * (t19 * (t19 * (-t66 * t19 + t3) - t65
     # * t47 * t28 * t68) + t61 * t63 * t36 * t60 + t47 * t18 * t1 * t41
     # * t68 * t26 + t54 * t47 * t39 ** 2 * t16 * t60 * t36) + t21 * t9 
     #* t50 * t6 * (t16 * (t54 * t58 * t64 + t66 * (t62 * t47 * t41 * t6
     # + t54 * t56 * t7 * t47 + t3 * t55) * t35 + t39 * t38 * (t65 * t6 
     #** 2 * t9 + t7 * (t54 * (t7 * t9 - 1) - t24 - t25) + t8 * t63) * t
     #36) + t54 * t19 * t9 * (-t67 - t38 * t19 * (t19 * t10 + t56)) + t5
     #4 * t57 * t39 * t60 * t38 * t35)
      t24 = t25 + t24
      t1 = t59 * t35 * t10 * (t20 * (t6 * ((t13 * t24 + t14 * t24) * t29
     # * t27 + t13 + t14) + t3 * (-t24 * t29 - t2) * t16) + t12 * t6 * (
     #t24 * t27 + t1) * t15)
      t2 = t20 * (t26 * (t3 * (t21 * (t3 * t15 * t16 + t14 * t20) - t52)
     # * t29 + t53 * t12 * (t10 * t12 + (t12 * t4 - t31) * t30 * t3) * t
     #27) - t53 * (t42 * t27 * t30 + t67 * t29) * t19)
      ret = 16 * t53 * (t26 * (t39 * (t43 * t18 * t6 * t38 * t35 + t45 *
     # t9 * t20 * (t20 * t32 * t3 + t6 * (t20 * t5 + t43) * t12)) + t9 *
     # (t5 * (t32 * t3 * t21 * t38 + t44 * t21 * t20) - t44 * t3 + t41 *
     # t12 * t38 * t6 * t21) + t47 * (t20 * (t6 * (t7 * t33 + t8 * t33) 
     #- t46 * t34) + t44 * t34) * t26 * t10) + t20 * (t6 * (-t49 * t19 *
     # t37 * t10 + t48 * (t11 ** 2 * t22 ** 2 + t42 * t23 ** 2) * t20 * 
     #t15 * t39) + t52 * t40 * (-t51 + t26)) * t29 * t27) - 32 * t17 + 1
     #28 * t62 * t9 ** 2 * t20 * t50 * t6 * (t16 * (t7 * t8 * t39 * t20 
     #* t36 - t3 * t10 * t34 * t35 - t56 * t3 * t26) + t60 * (-t3 * t39 
     #* t35 + t61 * t34 * t36) + t19 * t3 * (t56 * t28 + t19)) + 64 * t5
     #4 - 8 * t1 - 4 * t2

      hjetmass_triangle_pmpm_s23_s234_0_dp = ret/32d0/(0,1d0)
      return

      end function
