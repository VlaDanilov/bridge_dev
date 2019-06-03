
      double complex function hjetmass_triangle_pmpm_s14_s124_0_dp 
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

      t1 = za(i1, i2)
      t2 = za(i1, i4)
      t3 = za(i2, i4)
      t4 = zb(i2, i1)
      t5 = zb(i4, i1)
      t6 = zb(i4, i2)
      t7 = t1 * t4
      t8 = t3 * t6
      t9 = t7 + t8
      if ( dreal(t9) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t9 ** 2) + t7 + t8
      t10 = za(i2, i3)
      t11 = 0.1D1 / t9
      t12 = 2 * t7 * t11
      t13 = t2 * t5
      t14 = t13 * (1 - t12)
      t15 = za(i3, i4)
      t16 = t5 * (2 * t2 * t10 * t4 * t11 + t15)
      t17 = zb(i3, i2)
      t12 = -t8 * t13 * t12 + t1 * (t14 * t4 - t16 * t17)
      t18 = zb(i3, i1)
      t19 = zb(i4, i3)
      t20 = t5 * t17
      t21 = t2 * (-2 * t20 * t3 * t11 + t18)
      t22 = za(i1, i3)
      t23 = 0.1D1 / t17
      t12 = 0.1D1 / t12
      t16 = 0.1D1 / t16
      t9 = 0.1D1 / t9
      t14 = 0.1D1 / t14
      t24 = 0.1D1 / t10
      t25 = t4 * t19
      t26 = t20 - t25
      t27 = t5 ** 2
      t28 = t27 ** 2
      t29 = t5 * t27
      t30 = t1 ** 2
      t31 = t30 ** 2
      t32 = t1 * t30
      t33 = t12 ** 2
      t34 = t12 * t33
      t35 = mt ** 2
      t36 = t3 ** 2
      t37 = t36 ** 2
      t38 = t3 * t36
      t39 = t11 ** 2
      t40 = t2 ** 2
      t41 = t2 * t40
      t42 = t4 ** 2
      t43 = t4 * t42
      t44 = t9 ** 2
      t45 = t16 ** 2
      t46 = t3 * t17
      t47 = t3 * t19
      t48 = t9 * t24
      t49 = t48 * t11
      t50 = t39 * t17
      t51 = t24 * t3
      t52 = t17 * t18
      t53 = t9 * t39
      t54 = t53 * t5
      t55 = t54 * t42
      t56 = t55 * t40
      t57 = t11 * t26
      t58 = t14 * t27
      t59 = t58 * t4
      t60 = t18 ** 2
      t61 = t10 * t16 * t45
      t62 = t32 * t34
      t63 = t35 * t39
      t64 = t63 * t24
      t3 = t59 * t37 * (t40 * (t53 * t8 * t30 * t4 * t17 * t33 * t24 * t
     #27 - t55 * t30 * t6 * t33 * t24 * (t1 * t18 + t47)) + t64 * t16 + 
     #t13 * (t16 * t44 * t24 + (t62 * t17 * t24 * (t15 ** 2 * t19 ** 2 +
     # t60 * t22 ** 2) + t61) * t39 * t42)) + t59 * t38 * (t30 * (t50 * 
     #t4 * (-t35 * (t3 * t13 * (-2 * t8 * t11 + 1) * t24 + t21) + t40 * 
     #(-t47 * t4 * t9 * t5 + t46 * t9 * t27)) * t33 + t49 * t4 * t2 * (t
     #13 * t3 * t11 * t14 * t26 - t18) * t12) + t32 * (t46 * t41 * t42 *
     # t29 * t39 * t24 * t34 + t56 * (t51 * t26 - t52) * t33 - t56 * t18
     # * t14 * t24 * t12) + t35 * t4 * t21 * t39 * t23 * t45 + t51 * (t2
     # * (t20 * t44 + t57 * t9) + t50 * t35) * t12 * t1 - t54 * t40 * t3
     #1 * t43 * t18 * t24 * t33)
      t26 = 0.1D1 / t2
      t15 = t15 * t19
      t19 = t18 * t22
      t27 = t19 + t15
      t21 = t21 * t26 - t18
      t26 = t13 * t4
      t39 = t11 * t27
      t46 = t13 * t24
      t35 = t35 * t24
      t47 = t14 * t11
      t53 = t17 ** 2
      t54 = t14 ** 2
      t55 = t8 * t24
      t56 = t48 * t32
      t65 = t9 * t16
      t66 = t24 * t14
      t67 = t47 * t9
      t68 = t67 * t28 * t42 * t37 * t40
      t49 = t49 * t13
      t10 = t68 * (t49 * t1 * t31 * t43 * t17 * t34 - t45 * t9 + (t17 * 
     #(-t39 * t55 * t34 + t48 * t33) + t13 * t10 * t34 * t11 * t9 * t17 
     #* t53) * t4 * t32 + (-t46 * t11 * t34 + t49 * t14 * t33) * t17 * t
     #42 * t31) + t68 * (-t39 * t31 * t17 * t34 * t24 * t42 + (-t62 * t5
     #3 * t27 + t61 + t13 * (t17 * (t56 * t8 * t14 * t33 + t55 * t62 * (
     #t8 * t9 - 1) + t56 * t54 * t12) + t53 * (t32 * t9 * t14 * t33 - t6
     #2) + t65 * (-t1 * t16 * t14 + t30 * t54 * t24 + t10 * t45))) * t11
     # * t4 + t9 * t1 * (t1 * (t17 * (t66 * t12 + t55 * t33) + t33 * t53
     #) + t66 * t16))
      t31 = t58 * t7
      t34 = t64 * t14
      t8 = t62 * t50 * t14 * t44 * t5 * t28 * t43 * t37 * t41 * (t8 * t1
     #7 + t7 * (t55 + t17))
      t28 = t2 * t60
      t14 = t11 * t5 * t4 * t36 * (t14 * t16 * (t5 * t36 * t24 + t28 * t
     #23) + (t66 * t20 * t36 + t28 * t14 + t51 * t18 * (t2 * t23 * t14 *
     # (t18 * t6 - t25) + 1)) * t12 * t1)
      ret = 64 * t3 - 32 * t47 * t5 * t4 * t38 * (t30 * (t59 * t48 * t40
     # * t18 * t12 + t26 * (t52 * t13 * t9 + t51 * (t18 * (t13 * t6 * t9
     # + t57 * t22) + t57 * t13 + t15 * t57)) * t33) + t46 * t18 * t42 *
     # (-t39 + t13 * (-t11 + t9)) * t33 * t32 + t35 * t23 * t16 * t21 + 
     #t35 * t21 * t12 * t1) + 256 * t10 - 16 * t31 * t24 * t12 * t18 * t
     #38 * t2 * (t7 * t11 * (-t12 * (t13 + t15 + t19) + t11) + t9) + 128
     # * t47 * t29 * t42 * t37 * t2 * (t16 * (t1 * (-t67 * t26 * t16 + t
     #34) - t65) + t17 * (t30 * t24 * (-t31 * t40 * t11 * t9 - t27 * t9 
     #+ t63 * t7 + t13 * t9 * (-t47 * t7 * t27 - 1)) * t33 + t62 * t24 *
     # t11 * t4 * (t13 * t27 + t15 * t19) + t34 * t30 * t12)) + 512 * t8
     # - 8 * t14

      hjetmass_triangle_pmpm_s14_s124_0_dp = ret/32d0/(0,1d0)
      return

      end function
