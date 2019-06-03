
      double complex function hjetmass_triangle_pmpm_s14_s134_0_dp 
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

      t1 = zb(i4, i1)
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = za(i3, i4)
      t5 = zb(i4, i3)
      t6 = t2 * t3
      t7 = t4 * t5
      t8 = t6 + t7
      if ( dreal(t8) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t8 = cg * cdsqrt(t8 ** 2) + t6 + t7
      t9 = za(i2, i4)
      t10 = za(i1, i4)
      t11 = 0.1D1 / t8
      t12 = 2 * t6 * t11
      t13 = t10 * t1
      t14 = t13 * (1 - t12)
      t15 = zb(i3, i2)
      t16 = za(i2, i3)
      t17 = zb(i4, i2)
      t11 = t10 * (-2 * t2 * t15 * t1 * t11 + t17)
      t12 = t3 * (t11 * t16 + t14 * t2) - t7 * t12 * t13
      t13 = za(i1, i2)
      t18 = zb(i2, i1)
      t19 = t13 * t3
      t20 = t9 * t5
      t11 = 0.1D1 / t11
      t21 = 0.1D1 / t8
      t12 = 0.1D1 / t12
      t2 = 0.1D1 / t2
      t22 = 0.1D1 / t15
      t23 = t11 ** 2
      t24 = t10 ** 2
      t25 = t12 ** 2
      t26 = t12 * t25
      t27 = t16 ** 2
      t28 = t3 ** 2
      t29 = t3 * t28
      t30 = t2 * t25
      t31 = t7 * t14
      t32 = t30 * t7
      t33 = t32 * t22 * t16
      t34 = t22 * t25
      t35 = t2 * t23
      t36 = t7 * t21 ** 2
      t37 = t36 * t14
      t38 = 0.1D1 / t10
      t1 = 0.1D1 / t1
      t39 = 0.1D1 / t16
      t40 = t10 * t16
      t41 = t4 * t13
      t42 = t41 - t40
      t43 = t3 * t16
      t44 = t43 * t12
      t45 = t44 - t11
      t46 = t2 ** 2
      t47 = mt ** 2
      t48 = t14 ** 2
      t49 = t8 ** 2
      t50 = t1 ** 2
      t51 = t7 * t9
      t52 = t51 * t10 * t22
      t53 = t1 * t14
      t54 = t22 * t2
      t55 = t12 * t3
      t56 = t11 * t39
      t36 = t36 * t50 * t49 * t48 * t2 * t46 * t22
      t57 = t9 * t17
      t58 = t13 * t18
      t59 = (t58 + t57) * t1
      t60 = t59 + t10
      t61 = t9 ** 2
      t62 = t1 * t22
      t63 = t21 * t8
      t64 = t11 * t15
      t65 = t7 * t2
      t66 = -t59 - t10
      t67 = t53 * t38 + 1
      t68 = t9 * t28
      t46 = t55 * t22 * t14 * (t68 * t1 + t21 * t1 * (t53 * ((t20 + t19)
     # * t38 * t4 - t43) + t51) * t46 * t49 + t55 * (t4 * (t20 * (t1 * (
     #t57 * t67 + t58 * t67 + t14) + t10) + t19 * t53 * (t59 * t38 + 1))
     # + t53 * t43 * t66) * t2 * t8)
      t51 = t28 * t12
      t59 = t28 * t22
      t67 = t59 * t12
      t69 = t53 * t2
      t19 = t69 * (t1 * t4 * t47 * (t56 * t3 * t22 + t51 * (t14 * t16 * 
     #t12 - t22) - t14 * t23 * t39) * t38 * (t20 + t19) + t67 * (t21 * (
     #t4 * (-t20 - t19) + t43 * t10) + t43 * t31 * (t58 * t10 + t57 * (t
     #58 * t1 + t10)) * t25) * t8 + t62 * t28 * t45 * t47 + t37 * t10 * 
     #(t28 * (-t30 * t27 - t33) + t35 - t34 * t16 * t29) * t49)
      t20 = t2 * t48 * t29
      t4 = t63 * t7 * t10 * (t16 * (t30 * t22 * t14 * t29 * t60 + t48 * 
     #t22 * t29 * (t3 * (t1 * (-t58 - t57) + t10 * (t6 * t21 - 1)) + t4 
     #** 2 * t5 ** 2 * t10 * t21 * t2 + t65 * t66) * t26 + t51 * t54 * t
     #21 * (t10 * t3 + t69 * t8)) + t20 * t66 * t26 * t27 + t20 * t10 * 
     #t21 * t16 * t27 * t15 * t26 + t2 * t11 * (t53 * (t11 * (-t64 * t14
     # + t3) - t63 * t54 * t3) - t21 * (t48 * t15 * t23 + t59) * t10))
      t5 = t53 * t8 * (t54 * t28 * t1 * t45 + (t55 * t2 - t56 * t2 + t67
     # * t39) * t38 * t61 + t68 * t54 * t12 * (-t41 * t39 + t53) * t38)
      t6 = -8
      ret = 64 * t37 * t24 * t8 * t3 * (t28 * (t27 * (-t31 * t2 * t26 + 
     #t30) + t33) + t29 * (t16 * (-t31 * t22 * t26 + t34) - t14 * t26 * 
     #t27) - t35) + t6 * (t65 * t1 * t48 * t8 * ((-t34 * t10 * t21 + t62
     # * (t13 ** 2 * t18 ** 2 + t17 ** 2 * t61) * t26) * t16 * t29 + t21
     # * (t54 * t60 * t8 + t9) * t25 * t16 * t28 + t1 * t23 * (t63 * t2 
     #- t64)) + t8 * (t28 * (-t52 * t21 * t2 * t12 + t32 * t21 * t14 * t
     #9 * (t7 * (t53 + t10) * t22 + t40)) + t29 * (t54 * t7 * t24 * t16 
     #* t48 * t26 + t21 * t14 * (t53 * (t16 * t2 * t42 + t7 * (t41 * t2 
     #+ t9) * t22) + t52) * t25) - t36 * t11 + t54 * t50 * t14 * t3 * (t
     #9 * (-t56 + t55) + t53 * t45 * t2) * t38 * t47 + t36 * t44 + t34 *
     # t48 * t21 * t1 * t42 * t28 ** 2)) + 4 * t46 - 16 * t19 - 32 * t4 
     #+ 2 * t5

      hjetmass_triangle_pmpm_s14_s134_0_dp = ret/32d0/(0,1d0)
      return

      end function
