
      double complex function hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_dp
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
      t2 = zb(i4, i3)
      t3 = za(i1, i2)
      t4 = zb(i3, i1)
      t5 = za(i2, i4)
      t6 = t3 * t4
      t7 = t2 * t5 + t6
      t8 = za(i1, i3)
      t9 = za(i2, i3)
      t10 = zb(i3, i2)
      t11 = za(i3, i4)
      t12 = t8 * t4
      t13 = t9 * t10
      t14 = t11 * t2
      t15 = t14 + t12 + t13
      if ( dreal(t15) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t14 = cg * cdsqrt(t15 ** 2) + t12 + t13 + t14
      t15 = zb(i2, i1)
      t16 = za(i1, i4)
      t17 = zb(i4, i2)
      t18 = t3 * t15
      t19 = t16 * t1
      t20 = t17 * t5 + t18 + t19
      t21 = 0.1D1 / t14
      t19 = -2 * t12 * t20 * t21 + t18 + t19
      t22 = 2 * t8 * t20
      t23 = t22 * t2 * t21 - t17 * t3
      t24 = t10 * t3 - t16 * t2
      t25 = t8 * t19
      t17 = t4 * (t9 * (-t22 * t10 * t21 + t16 * t17) + t25)
      t22 = 2 * t9 * t20
      t26 = t22 * t2 * t21 + t1 * t3
      t22 = -t22 * t4 * t21 + t1 * t5
      t27 = 2 * t11 * t4 * t20 * t21 - t15 * t5
      t28 = t8 * (t22 * t10 + t19 * t4)
      t15 = 0.1D1 / t15
      t14 = 0.1D1 / t14
      t29 = 0.1D1 / t5
      t30 = 0.1D1 / t3
      t31 = 0.1D1 / t16
      t17 = 0.1D1 / t17
      t32 = 0.1D1 / t20
      t33 = t4 * t29
      t34 = t2 * t30
      t35 = t34 + t33
      t36 = t17 ** 2
      t37 = t17 * t36
      t38 = t9 ** 2
      t39 = t4 ** 2
      t40 = t4 * t39
      t41 = mt ** 2
      t42 = t1 * t15
      t11 = 0.1D1 / t11
      t27 = 0.1D1 / t27
      t28 = 0.1D1 / t28
      t43 = 0.1D1 / t23
      t44 = t1 * t24
      t45 = t2 * t19
      t46 = t4 * t23
      t47 = -t46 - t44 - t45
      t48 = t2 ** 2
      t49 = t9 * t24
      t50 = t49 * t29
      t33 = t33 * t9
      t51 = t45 * t43
      t52 = t39 * t9
      t53 = t4 * t15
      t54 = t12 * t2
      t55 = (t3 * t31 + t42) * t32
      t56 = t4 * t17
      t57 = t39 * t24
      t16 = t4 * (t39 * (t20 * t23 * t24 * t14 * t29 * t36 * t38 + t50 *
     # t25 * t42 * t20 * t14 * t36) + t4 * (-t50 * t42 * t20 * t14 * t17
     # + t19 * (t42 * t10 - t2) * t36 * t29 * t14 * t24 * t20 * t38) - t
     #55 * t2) + t29 * t15 * (t9 * (t4 * (-t56 * t24 * t26 + t11 * t7) -
     # t45 * t7 * (t11 * t43 + t56)) - t4 * t7 * (t8 * t2 * t28 + t27) *
     # t22 + t19 * (t57 * t36 * (t10 * t26 + t46) - t2 * (t10 * t5 + t16
     # * t4) * t11 ** 2 * t43) * t38) * t21 * t30 * t41
      t26 = t46 + t44
      t26 = -t23 * t26 + t50 * t26
      t58 = t50 - t23
      t59 = t43 ** 2
      t60 = t14 * t15 * t20
      t26 = t60 * t9 * (t29 * (t4 * (t44 * t43 + t4) + t34 * t49 * t43 *
     # (t43 * (t45 + t44) + t4) + t48 * t19 ** 2 * t59 + t44 * t45 * t59
     #) * t11 + t56 * t34 * (-t46 - t44 - t45) + t30 * t39 * (t12 * t26 
     #+ t13 * t26 + t45 * (t12 * t58 + t13 * t58)) * t36)
      t8 = t57 * t17 * t29 * t14 * t20 * t38 * (t46 * t19 * t36 * (t13 +
     # t12) + t60 * (-t46 * t38 * t10 ** 2 * t19 * t36 - t40 * t8 ** 2 *
     # t19 * t23 * t36 + t39 * t8 * t23 * t17 - t54 * t19 * t17 + t2 + t
     #13 * t17 * (t46 - t45)) * t30)
      t12 = -16
      t13 = 4
      t44 = 128
      ret = t44 * (t25 * t20 ** 2 * t14 ** 2 * t9 * t38 * t39 ** 2 * t10
     # * t24 * t23 * t30 * t29 * t15 * t37 + (t31 * (t29 * t6 + t2) + t4
     #2 * t35) * t32 * t21 * t7 * t41) + t13 * (t54 * t22 * t15 * t35 * 
     #t28 + t53 * t35 * t27 * t22 + t52 * (-t23 * t47 + t47 * t50) * t36
     # + t53 * (-t52 * t23 * t29 + t46 * (-t34 * t9 - t1) - t1 * t2 * (t
     #30 * t49 + t19)) * t17 + t15 * (t48 * t9 * t19 * t30 * t43 + t51 *
     # (t33 - t1) + t1 * t4 * (t43 * t50 - 1)) * t11) + 32 * t16 + t12 *
     # (-t18 * t38 * t40 * t19 * t24 * t23 * t29 * t37 + t5 * t48 * t32 
     #* (t30 * t42 + t31) + t3 * (t1 * t19 * t36 * t49 + t55) * t29 * t3
     #9 - t51 * t33 * t20 * t14 * t11 * t15) + 8 * t26 - 64 * t8

      hjetmass_qqbgg_triangle_pmpp_s124_mhsq_0_dp = ret/32d0/(0,1d0)
      return

      end function
