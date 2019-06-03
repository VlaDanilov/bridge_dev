
      double complex function hjetmass_triangle_pmpm_s34_s234_0_dp 
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

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = za(i3, i4)
      t4 = zb(i2, i1)
      t5 = zb(i4, i3)
      t6 = za(i2, i3)
      t7 = zb(i3, i2)
      t8 = zb(i4, i2)
      t9 = t6 * t7
      t10 = t1 * t8
      t11 = t10 + t9
      if ( dreal(t11) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t11 ** 2) + t10 + t9
      t12 = 0.1D1 / t11
      t13 = t4 * t5
      t14 = 2 * t13
      t15 = t3 * (t14 * t1 * t12 - t2)
      t16 = za(i1, i2)
      t17 = za(i1, i3)
      t18 = za(i1, i4)
      t19 = zb(i4, i1)
      t20 = t17 * t2
      t21 = t18 * t19
      t22 = -t14 * t16 * t3 * t12 + t20 + t21
      t14 = t3 * (t14 * t6 * t12 + t19)
      t23 = (t14 * t7 + t15 * t8 - t22 * t4) * t16
      t24 = t6 * t2
      t25 = t1 * t19
      t26 = t25 + t24
      t27 = t3 * t5
      t28 = 0.1D1 / t4
      t23 = 0.1D1 / t23
      t22 = 0.1D1 / t22
      t29 = 0.1D1 / t5
      t30 = 0.1D1 / t3
      t31 = 0.1D1 / t16
      t14 = 0.1D1 / t14
      t11 = 0.1D1 / t11
      t32 = -t27 - t20 - t21
      t33 = t21 + t20
      t34 = t25 + t24
      t35 = t23 ** 2
      t36 = t23 * t35
      t37 = mt ** 2
      t38 = t16 ** 2
      t39 = t16 * t38
      t40 = t22 ** 2
      t41 = t1 * t31
      t42 = t11 * t14
      t43 = t13 * t1
      t44 = t11 * t23
      t45 = t26 * t7
      t46 = t14 * t26
      t47 = t46 * t15
      t48 = t10 + t9
      t49 = t4 ** 2
      t50 = t26 ** 2
      t51 = t1 * t6 * t14
      t52 = t27 * t11 ** 2
      t53 = t52 * t14
      t54 = t14 ** 2
      t55 = t15 ** 2
      t56 = t1 ** 2
      t57 = t6 ** 2
      t58 = t37 * t12 ** 2
      t59 = t10 * t11
      t60 = t6 * t55
      t61 = t57 * t55
      t62 = t56 * t31
      t63 = t16 * t15
      t64 = t15 * t6
      t65 = t14 * t50 * t7
      t66 = t63 * t44
      t52 = t63 * t65 * t35 * t4 * (t9 * (-t58 * t1 - t66 * t33) + t27 *
     # t11 * (t64 * (t42 * (t16 * t4 + t9) + t7 * t16 * (t9 * t11 - 1) *
     # t23) + t1) + t1 * t11 * t33) + t65 * (t22 * (t58 * t41 * (-t64 * 
     #t14 + t1) + t11 * t15 * t22 * (t63 * t22 + t1) + t52 * (-t60 * t22
     # * t14 + t55 * t16 * t40 + t61 * t31 * t54 + t62)) + t4 * (t23 * (
     #t56 * (t58 + t52) + t52 * t57 * t55 * t54 - t58 * t51 * t15) + t59
     # * t38 * t55 * (t27 * (t59 - 1) - t20 - t21) * t36 + t60 * t53 * t
     #10 * t16 * t35) + t52 * t38 ** 2 * t4 * t49 * t55 * t36 + t11 * t3
     #9 * t55 * t32 * t36 * t49)
      t54 = t2 ** 2
      t56 = t1 * t4
      t57 = t1 * t2
      t58 = t7 * t19
      t8 = t15 * (t23 * (t1 * (t14 * (-t64 * t2 + t56 * t26) + t57) * t3
     #0 + t46 * t2 * (t16 * t2 + (t2 * t8 - t58) * t28 * t1) * t29) + t4
     #6 * (t54 * t28 * t29 + t62 * t30) * t22)
      t58 = -t13 + t58
      t59 = t45 * t2
      t60 = t3 * t26
      t2 = t46 * (t23 * (t11 * (t6 * (t58 * t1 * t14 * t55 - t59 * t3 * 
     #t14 * t15) + t57 * t45 * t3 + t61 * t2 * t7 * t14) + t38 * (-t45 *
     # t13 * t3 * t55 * t35 + t44 * t15 * t4 * (t15 * t58 * t1 + t7 * (-
     #t60 + t64) * t2)) + t66 * (t7 * (t15 * (t24 * t48 + t25 * t48) - t
     #60 * t2 * t48) - t43 * t15 * t48)) + t15 * ((-t26 * t22 * t40 * t1
     #6 - t4 * t36 * (t18 ** 2 * t19 ** 2 + t54 * t17 ** 2) * t26 * t38)
     # * t15 * t7 + t57 * t37 * (t28 * t22 * t31 + t23)) * t30 * t29)
      ret = 32 * t47 * (t45 * (t21 * t20 * t4 * t15 * t38 * t36 + (t23 *
     # (t4 * (-t15 * t23 * t38 + t1 * t27 * (-2 * t10 * t12 + 1) * t23 *
     # t16) + t1) + t28 * (t15 * t40 + t41 * t22)) * t12 * t37) * t30 * 
     #t29 + t44 * t1 * (t34 * t7 - t43) + t7 * (t4 * (t42 * t6 * t32 * t
     #35 * t16 + t33 * t36 * t38) - t6 * t40 * t14 * t11) * t26 * t15) +
     # 128 * t53 * t50 * t15 * t7 * (t4 * (t9 * t10 * t15 * t38 * t36 - 
     #t1 * t16 * t48 * t35 - t51 * t23) + t49 * (t15 * t39 * t48 * t36 -
     # t1 * t38 * t35) + t1 * t22 * (-t6 * t31 * t14 + t22)) + 64 * t52 
     #+ 4 * t8 - 16 * t2 - 8 * t47 * t35 * t16 * (t15 * (t7 * ((-t20 * t
     #34 - t21 * t34) * t30 * t29 - t24 - t25) + t56 * (t30 * t33 + t5))
     # + t59 * (t29 * t33 + t3))

      hjetmass_triangle_pmpm_s34_s234_0_dp = ret/32d0/(0,1d0)
      return

      end function
