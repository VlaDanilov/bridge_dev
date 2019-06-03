
      double complex function hjetmass_triangle_pmpm_s123_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i2, i4)
      t2 = zb(i3, i1)
      t3 = zb(i2, i1)
      t4 = za(i3, i4)
      t5 = t1 * t3 + t2 * t4
      t6 = za(i1, i4)
      t7 = zb(i4, i1)
      t8 = zb(i4, i2)
      t9 = zb(i4, i3)
      t10 = t6 * t7
      t11 = t1 * t8
      t12 = t4 * t9
      t13 = t10 + t11 + t12
      if ( dreal(t13) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i2, i3)
      t13 = zb(i3, i2)
      t14 = za(i1, i2)
      t15 = za(i1, i3)
      t16 = t14 * t3
      t17 = t15 * t2
      t18 = t12 * t13 + t16 + t17
      t19 = t1 * t13 + t2 * t6
      t11 = 0.1D1 / t11
      t10 = -2 * t10 * t18 * t11 + t16 + t17
      t16 = -2 * t6 * t18
      t17 = t16 * t9 * t11 - t13 * t14
      t15 = t1 * (t16 * t8 * t11 + t13 * t15)
      t16 = t7 * (t17 * t4 + t15)
      t20 = -2 * t1 * t18 * t7 * t11 + t12 * t2
      t8 = t20 * t8
      t21 = t6 * (t10 * t7 + t8)
      t22 = -2 * t4 * t18 * t7 * t11 - t12 * t3
      t8 = t6 * (t22 * t9 + t8)
      t15 = t7 * (t10 * t6 + t15)
      t8 = 0.1D1 / t8
      t23 = 0.1D1 / t3
      t12 = 0.1D1 / t12
      t24 = 0.1D1 / t9
      t14 = 0.1D1 / t14
      t25 = 0.1D1 / t6
      t22 = 0.1D1 / t22
      t15 = 0.1D1 / t15
      t16 = 0.1D1 / t16
      t26 = 0.1D1 / t10
      t27 = 0.1D1 / t7
      t21 = 0.1D1 / t21
      t28 = 0.1D1 / t17
      t29 = 0.1D1 / t4
      t30 = 0.1D1 / t13
      t31 = t20 ** 2
      t32 = t20 * t31
      t33 = t2 * t23
      t29 = t1 * t10 * t29
      t34 = t9 * t26
      t35 = t7 * t8
      t36 = t35 * t23
      t37 = t33 * t30
      t38 = t14 * t12
      t18 = 0.1D1 / t18
      t39 = t2 ** 2
      t40 = t14 * t27
      t41 = t12 * t24
      t42 = t30 * t11
      t43 = t9 ** 2
      t44 = t7 ** 2
      t45 = t1 * t17 * t25
      t46 = t5 * t43 * t26
      t47 = t19 * t44 * t8 * t23 * t24
      t48 = t5 * t9
      t11 = t38 * t11
      t49 = t11 * t1
      t50 = t49 * t6
      t51 = t21 ** 2
      t52 = t6 ** 2
      t53 = t8 ** 2
      t3 = t1 * (-t13 * t4 + t3 * t6)
      t13 = t4 * t19
      t54 = t11 * t37
      t3 = (-t54 * t48 * t7 * t51 * t31 - t46 * t42 * t38 * t51 * t32) *
     # t6 * t52 + t52 * (t11 * (t43 * t26 * t21 * t30 * (-t5 * t26 * t27
     # + t3 * t21) + (-t13 - t3) * t22 * t23 * t53 * t44) * t32 + t54 * 
     #t9 * t7 * (-t13 * t53 + t3 * (t51 - t53)) * t31) + t1 * t2 * t39 *
     # t23 * t30 * t18 * (t41 + t40) - t47 * t11 * t6 * t4 * t32 * t22 *
     #* 2
      t4 = t50 * t37 * t20 * (t19 * t7 * t21 + t48 * t8)
      ret = 8 * t38 * t2 * t1 * (t1 * (t7 * (t1 * t24 * t25 * t15 * t26 
     #* t30 * t17 ** 2 + t33 * t24 * t15 * t30 * t17) + t10 * t23 * t16 
     #* (t2 * t30 + t29 * t28)) + t6 * (t31 * (t34 * t27 * t21 * t30 + t
     #36 * t24 * t22) + t37 * (t8 + t21) * t20)) + 64 * t42 * t18 * t23 
     #* t39 * t1 * (t19 * (-t41 * t7 + t14) + t5 * (-t40 * t9 + t12)) - 
     #16 * t49 * (t1 * (t30 * (t33 * t5 * t9 * t10 * t16 + (t45 * t15 + 
     #t33 * (t15 - t16) * t10) * t19 * t7 + t24 * t17 * t15 * (t45 * t26
     # + t33) * t19 * t44) + t29 * t5 * t23 * t16 * (t9 * t10 * t28 + t7
     #)) + t6 * (t31 * (-t46 * t27 * t21 * t30 - t47 * t22) - t37 * (t35
     # * t19 + t48 * t21) * t20)) - 80 * t50 * t31 * (t34 * t19 * t30 * 
     #t21 + t36 * t5 * t22) + 32 * t3 - 48 * t4

      hjetmass_triangle_pmpm_s123_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
