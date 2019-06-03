
      double complex function hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_dp
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
      if ( dreal(t19) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t18 = cg * cdsqrt(t19 ** 2) + t16 + t17 + t18
      t19 = 0.1D1 / t18
      t8 = -2 * t16 * t9 * t19 + t7 + t8
      t20 = t1 * t13 + t15 * t3
      t21 = -2 * t10 * t9
      t22 = t21 * t15 * t19 - t1 * t6
      t21 = t21 * t13 * t19 + t3 * t6
      t23 = t10 * t8
      t24 = t11 * (t12 * t21 + t23)
      t3 = t3 * t11
      t25 = t5 * t13
      t26 = t25 + t3
      t27 = t1 * t11
      t28 = t5 * t15
      t29 = t27 - t28
      t30 = -2 * t12 * t9 * t11 * t19 + t4 * t5
      t31 = t11 * t8
      t32 = t10 * (t13 * t30 + t31)
      t33 = -2 * t14 * t9 * t11 * t19 - t2 * t5
      t1 = 0.1D1 / t1
      t22 = 0.1D1 / t22
      t4 = 0.1D1 / t4
      t18 = 0.1D1 / t18
      t34 = 0.1D1 / t2
      t24 = 0.1D1 / t24
      t35 = t11 ** 2
      t36 = t11 * t35
      t37 = t12 ** 2
      t38 = t24 ** 2
      t39 = t24 * t38
      t23 = t23 * t36 * t38
      t40 = t35 * t24
      t41 = t8 * t15 * t22
      t42 = t13 * t8
      t43 = t42 * t37 * t35
      t44 = t34 * t22
      t45 = t1 * t5
      t46 = t45 * t18
      t47 = 0.1D1 / t12
      t48 = 0.1D1 / t14
      t33 = 0.1D1 / t33
      t32 = 0.1D1 / t32
      t49 = 0.1D1 / t15
      t50 = t11 * t21
      t51 = t42 - t50
      t32 = t10 * t32
      t33 = t33 * t49 - t32
      t52 = t11 * t24
      t53 = t8 ** 2
      t54 = t22 ** 2
      t55 = mt ** 2
      t56 = t15 ** 2 * t53 * t54
      t57 = t8 * (t22 * t48 - t52)
      t58 = t11 * t55
      t59 = t45 * t34
      t60 = t14 * t35
      t61 = t60 * t38
      t54 = t15 * t54
      t62 = t34 * t19 * t4 * t1
      t6 = 0.1D1 / t6
      t63 = 0.1D1 / t9
      t3 = t25 + t3
      t15 = t45 * t15
      t25 = t32 * t26 * t30
      t32 = t14 * t47
      t64 = t9 ** 2
      t65 = t2 * t4
      t66 = t34 * t18 ** 2 * t4 * t1 * t20
      t19 = t11 * ((-t6 * (t65 + t32) + t45 * (-t32 * t34 - t4)) * t63 *
     # t19 * t29 * t55 + t66 * t8 * t12 * t64 * (t61 * (t16 * (t42 * t12
     # * t24 - 1) - t17) + t54))
      t67 = t5 * t35
      t68 = t35 * t34 * t4
      t1 = t45 * (t22 * (-t11 * t34 - t48 * t5) + t52 * (-t12 * t11 * t4
     # + t5)) * t20 + t68 * (-t11 * t49 + t45) - t17 * t67 * t53 * t4 * 
     #t38 + t31 * t4 * (t12 * (t52 * t1 * t34 * t26 + t67 * t21 * t38) +
     # t44 * (t15 - t11)) + t68 * t25 * t1
      t17 = t36 * t14
      t10 = t64 * (t66 * (t14 * (t10 ** 2 * t11 * t35 ** 2 * t53 * t39 +
     # t24 * t36) - t22 * (t56 + t35)) * t12 + t17 * t66 * t12 * t37 * t
     #53 * t13 ** 2 * t39) + t9 * (t17 * t8 * t38 * (-t16 * t8 * t24 + 1
     #) * t18 * t4 * t20 * t12 - t17 * t18 * t37 * t53 * t13 * t20 * t4 
     #* t39) - t58 * t5 * t47 * t63 * (t59 + t6)
      ret = 24 * t46 * t20 * t9 * (t43 * t4 * t38 + t44 * (t41 - t11) + 
     #(-t40 + t23) * t4 * t12) + 8 * t59 * t4 * ((t12 * (t23 * t51 + t40
     # * (-t42 + t50)) + t35 + t56 + t43 * t51 * t38) * t18 * t9 + t58 *
     # (t11 * t33 * t47 * t30 + t57)) - 32 * t62 * t53 * t12 * t55 * (t5
     #4 * t26 * t48 + t61 * (t11 * t20 - t13 * t29)) + 16 * t11 * (t63 *
     # (-t3 * t6 + t4 * (t5 * (t15 - t11) + (t28 - t27) * t6 * t2) - t59
     # * t3) + t62 * (t11 * (t32 * t29 * t30 * t33 + t29 * t49 - t25) + 
     #t57 * t26 * t12) * t55 - t7 * t60 * t12 * t53 * t20 * t4 * t39 - t
     #41 * t46 * t9 * t34 * t4) - 128 * t19 + 4 * t1 - 64 * t10 + 48 * t
     #62 * t58 * t8 * t29 * (-t52 * t14 + t22) - 12 * t65 * t67 * t12 * 
     #t8 * t20 * t38

      hjetmass_qqbgg_triangle_pmmp_s123_mhsq_0_dp = ret/32d0/(0,1d0)
      return

      end function
