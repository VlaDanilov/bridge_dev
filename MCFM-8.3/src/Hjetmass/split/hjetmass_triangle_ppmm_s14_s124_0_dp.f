
      double complex function hjetmass_triangle_ppmm_s14_s124_0_dp 
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

      t1 = za(i1, i4)
      t2 = zb(i4, i1)
      t3 = za(i1, i2)
      t4 = zb(i2, i1)
      t5 = za(i2, i4)
      t6 = zb(i4, i2)
      t7 = t3 * t4
      t8 = t5 * t6
      t9 = t7 + t8
      if ( dreal(t9) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t9 = cg * cdsqrt(t9 ** 2) + t7 + t8
      t10 = za(i2, i3)
      t11 = za(i3, i4)
      t12 = zb(i3, i2)
      t13 = 0.1D1 / t9
      t14 = 2 * t7 * t13
      t15 = t1 * t2
      t16 = t15 * (1 - t14)
      t17 = zb(i3, i1)
      t18 = zb(i4, i3)
      t19 = t1 * (2 * t3 * t12 * t2 * t13 + t18)
      t14 = t8 * t15 * t14
      t20 = t4 * (-t10 * t19 + t16 * t3) - t14
      t21 = za(i1, i3)
      t22 = t2 * (2 * t1 * t10 * t4 * t13 + t11)
      t14 = -t14 + t3 * (-t12 * t22 + t16 * t4)
      t23 = t11 * t6
      t24 = t21 * t4 + t23
      t25 = 0.1D1 / t3
      t19 = 0.1D1 / t19
      t26 = 0.1D1 / t10
      t14 = 0.1D1 / t14
      t9 = 0.1D1 / t9
      t27 = 0.1D1 / t1
      t28 = 0.1D1 / t16
      t20 = 0.1D1 / t20
      t29 = 0.1D1 / t17
      t30 = 0.1D1 / t2
      t31 = t11 ** 2
      t32 = t12 ** 2
      t33 = t19 ** 2
      t34 = t4 ** 2
      t35 = t34 ** 2
      t36 = t4 * t35
      t37 = t4 * t34
      t38 = mt ** 2
      t39 = t5 ** 2
      t40 = t20 ** 2
      t41 = t20 * t40
      t42 = t1 ** 2
      t43 = t16 ** 2
      t44 = t6 ** 2
      t45 = t4 * t29
      t46 = t45 * t12 * t33
      t47 = t32 * t16 * t29 * t19 * t33
      t48 = t10 * t37
      t26 = t19 * t26
      t49 = t4 * t20
      t50 = t10 * t12
      t51 = t3 * t14
      t52 = t15 * t29
      t53 = t30 * t25
      t54 = t53 * t16
      t55 = 0.1D1 / t12
      t56 = 0.1D1 / t22
      t57 = t11 * t18
      t15 = (t57 + t15) * t29
      t58 = t15 + t21
      t59 = t10 ** 2
      t60 = t13 ** 2
      t61 = t9 ** 2
      t62 = t25 * t16
      t63 = t62 * t12
      t3 = t5 * t1 * (-t48 * t52 * t61 * t5 * t39 * t43 * t6 * t44 * t25
     # * t41 - t52 * t38 * t10 * t36 * t13 * t60 * t28 * (t55 * t56 + t5
     #1) + t48 * t41 * t25 * t43 * t58 * t9 * t44 * t39 + t8 * (t16 * (t
     #37 * (t63 * t58 * t41 * t59 + t25 * (-t15 - t21) * t40 * t10) + t1
     #6 * t58 * t41 * t10 * t35 - t47 * t25) * t9 + t52 * (t10 * (-t3 * 
     #t36 * t43 * t41 - t37 * t25 * t20) - t37 * t10 * t59 * t32 * t43 *
     # t25 * t41 - t25 * t19 * (t32 * t43 * t33 + t34)) * t61))
      t15 = t1 * t29
      t32 = t38 * t34
      t36 = t1 * t9
      t64 = t26 * t38
      t65 = t49 * t10
      t66 = t11 * t2
      t1 = t5 * (t5 * (t25 * (t10 * (t32 * t12 * t43 * t24 * t27 * t13 *
     # t29 * t30 * t40 + t43 * t6 * t37 * (-t1 * t21 + t57 * (-t21 * t30
     # - t15)) * t41) - t29 * t19 * (t64 * t12 * t43 * t24 * t13 * t27 *
     # t30 + t36 * t34 * t6)) - t45 * t38 * t16 * t6 * t13 * t30 * (t65 
     #+ t19) * t25 ** 2) + t15 * t60 * t28 * t35 * t38 * (-t55 * (t66 * 
     #t56 + 1) + t51 * (-t66 - t22)))
      t22 = t25 * (t8 + t50) + t4
      t35 = t8 * t13
      t4 = t5 * t4 * (t11 * (t16 * (-t45 * t25 * (t27 * t38 * t30 + t8 *
     # t9) * t20 - t35 * (t52 + t21) * t40 * t34 - t64 * t53 * t27 * t29
     #) + t45 * t8 * t40 * t9 * t22 * t43) + t15 * t37 * t13 * t55 * t28
     # - t35 * t31 * t34 * t16 * t18 * t29 * t40)
      t15 = t43 * t20
      t35 = t12 * t25
      t8 = t45 * t61 * t16 * t6 * t39 * t2 * t42 * (t34 * (t59 * (t63 * 
     #t8 * t41 - t35 * t40) - t8 * t25 * t40 * t10) + t37 * (t10 * (t8 *
     # t16 * t41 - t40) + t12 * t16 * t41 * t59) + t35 * t33)
      t12 = t62 * t23 * t40 * t39 * t34 * t58
      t32 = t32 * t54 * t5 * t13 * t29 * (t65 + t19)
      ret = -16 * t39 * (t54 * (t6 * (-t48 * t16 * t41 * (t31 * t18 ** 2
     # * t29 + t17 * t21 ** 2) - t46 - t47) - t45 * (t49 + t26) * t13 * 
     #t27 * t24 * t38) + t52 * t13 * t37 * (t5 * t44 * t11 * t16 * t9 * 
     #t40 + t51 * t34 * t21 * t28 * t13 + t23 * (t9 * (t16 * (t7 + t50) 
     #* t40 - t20) + t7 * t14 * t28 * t13)) - t42 * t10 * t37 * t2 * t43
     # * t6 * t25 * t29 * t41) - 64 * t3 - 32 * t1 - 24 * t23 * t52 * t2
     #0 * t9 * t39 * t34 * (t22 * t20 * t16 - t25) - 8 * t4 - 4 * t25 * 
     #t5 * (t6 * (t31 * (t49 * (-t15 * t45 * t5 * t18 * t27 * t30 + 1) +
     # t26) + t20 * t5 * t34 * (-t29 * (t15 + 1) + t16 * (-t16 * t21 * t
     #20 + t29) * t30 * t27) * t11) + t34 * t16 * t29 * t30 * (t49 * t21
     # * t5 * t27 + t19)) + 128 * t8 + 12 * t12 + 80 * t32 + 96 * t36 * 
     #t62 * t46 * t39 * t6

      hjetmass_triangle_ppmm_s14_s124_0_dp = ret/32d0/(0,1d0)
      return

      end function
