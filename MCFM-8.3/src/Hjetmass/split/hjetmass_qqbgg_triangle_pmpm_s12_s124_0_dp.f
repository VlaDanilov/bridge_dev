
      double complex function hjetmass_qqbgg_triangle_pmpm_s12_s124_0_dp
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
      t2 = zb(i2, i1)
      t3 = za(i1, i4)
      t4 = zb(i4, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = t3 * t4
      t8 = t5 * t6
      t9 = t7 + t8
      if ( dreal(t9) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t9 ** 2) + t7 + t8
      t10 = zb(i4, i3)
      t11 = za(i3, i4)
      t12 = 0.1D1 / t9
      t13 = -2 * t7 * t12
      t14 = t1 * t2
      t15 = t14 * (1 + t13)
      t16 = za(i2, i3)
      t17 = t2 * (-2 * t1 * t11 * t4 * t12 - t16)
      t13 = t8 * t14 * t13 + t3 * (t10 * t17 + t15 * t4)
      t18 = zb(i3, i1)
      t19 = zb(i3, i2)
      t20 = za(i1, i3)
      t21 = t18 * t3 + t19 * t5
      t22 = 0.1D1 / t11
      t17 = 0.1D1 / t17
      t9 = 0.1D1 / t9
      t13 = 0.1D1 / t13
      t23 = t16 * t19
      t24 = t18 * t20
      t25 = t14 + t23 + t24
      t26 = t13 ** 2
      t27 = t13 * t26
      t28 = t4 ** 2
      t29 = t17 ** 2
      t30 = t10 ** 2
      t31 = t5 ** 2
      t32 = t5 * t31
      t33 = t2 ** 2
      t34 = t3 ** 2
      t35 = t34 ** 2
      t36 = t3 * t34
      t37 = t8 * t22
      t38 = t14 * t9
      t39 = t25 * t22
      t40 = t9 * t12
      t41 = 0.1D1 / t10
      t42 = t24 + t23
      t43 = t9 ** 2
      t44 = mt ** 2
      t45 = t1 * t33
      t46 = t45 * t5 * t43
      t47 = t46 * t26
      t48 = t12 * t22
      t49 = t48 * t5 * t4
      t50 = t44 * t4 * t21 * t12 ** 2
      t51 = t31 * t4
      t52 = t51 * t2
      t53 = 0.1D1 / t1
      t54 = t2 * t10
      t55 = t4 * t19
      t56 = t55 - t54
      t57 = t18 ** 2
      t58 = t1 ** 2
      t59 = t6 * t18
      t60 = t10 * t18
      t61 = t45 * t40
      t21 = t48 * t44 * t21 * t53
      t48 = t9 * t5
      t62 = t3 * t13
      t16 = t51 * (t17 * (-t21 * t41 + (t11 * t4 * t12 * t29 + t17 * t9)
     # * t5 * t33) + t36 * (t49 * t10 * t33 * (-t16 ** 2 * t19 ** 2 - t5
     #7 * t20 ** 2 - t33 * t58) * t27 + t61 * t4 * (t60 + (t59 - t54 + t
     #55) * t22 * t5) * t26) + t48 * t33 * (t10 * (t1 * t56 * t12 + t39)
     # + t37 * t1 * t12 * t56) * t26 * t34 + t62 * t21 + t61 * t35 * t18
     # * t28 * t22 * t26)
      t19 = -t55 + t54
      t21 = t4 * t2
      t51 = t62 * t22
      t54 = t17 * t41
      t61 = t22 * t18
      t19 = t31 * (t34 * (t21 * t18 * t22 * t9 * t13 + t21 * (t60 * t38 
     #+ (t18 * (t12 * t19 * t20 + t38 * t6) + t14 * t12 * t19 + t23 * t1
     #2 * t19) * t22 * t5) * t26) + t51 * (-t44 * t18 * t53 + t48 * t2 *
     # t56) + t54 * t44 * t18 * t53 * t22 + t61 * t28 * t2 * (-t12 * t42
     # + t14 * (-t12 + t9)) * t26 * t36)
      ret = 128 * t40 * t32 * t28 * t2 * t33 * t1 * (t36 * (-t14 * t11 *
     # t9 * t27 * t10 * t30 + t37 * (t14 * (-t8 * t9 + 1) + t23 + t24) *
     # t27 * t10 + t25 * t27 * t30) + t39 * t4 * t27 * t10 * t35 + t11 *
     # t17 * t29 * (t38 + 1) - t38 * t3 * t35 * t28 * t10 * t22 * t27) -
     # 64 * t52 * (t29 * (t50 * t41 - t46) + t34 * (t47 * t30 + (t45 * t
     #31 * t6 * t22 * t43 - t50) * t26 * t10) + (t47 * t4 * t22 + t49 * 
     #t2 * (t14 * t42 + t24 * t23) * t27) * t10 * t36) + 32 * t16 - 256 
     #* t27 * t43 * t12 * t10 * t32 * t28 * t36 * t33 ** 2 * t58 * (t8 *
     # t10 + t7 * (t37 + t10)) + 16 * t19 + 8 * t61 * t52 * t13 * t34 * 
     #(-t13 * t25 + t12) - 4 * t5 * (t57 * (t54 - t62) + t51 * t18 * (-t
     #15 * t53 + t41 * (-t59 + t55)) * t5 + t2 * t53 * (-t62 * t10 + t17
     #) * t22 * t31)

      hjetmass_qqbgg_triangle_pmpm_s12_s124_0_dp = -ret/32d0/(0,1d0)
      return

      end function
