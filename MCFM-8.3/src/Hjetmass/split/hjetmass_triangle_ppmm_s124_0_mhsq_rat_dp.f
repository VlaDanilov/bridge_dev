
      double complex function hjetmass_triangle_ppmm_s124_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg

      t1 = za(i3, i4)
      t2 = zb(i2, i1)
      t3 = za(i1, i2)
      t4 = za(i1, i4)
      t5 = zb(i4, i1)
      t6 = zb(i4, i3)
      t7 = za(i1, i3)
      t8 = za(i2, i3)
      t9 = zb(i3, i1)
      t10 = -t1 * t5 + t2 * t8
      t11 = za(i2, i4)
      t12 = zb(i4, i2)
      t13 = t3 * t2
      t14 = t4 * t5
      t15 = t11 * t12
      t16 = t15 + t13 + t14
      t17 = zb(i3, i2)
      t18 = t7 * t9
      t19 = t8 * t17
      t20 = t1 * t6
      t21 = t18 + t19 + t20
      if ( dreal(t21) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t19 = cg * cdsqrt(t21 ** 2) + t18 + t19 + t20
      t19 = 0.1D1 / t19
      t20 = 2 * t1
      t21 = t20 * t9 * t16 * t19 - t11 * t2
      t22 = t20 * t17 * t16 * t19 + t2 * t4
      t23 = -2 * t8 * t9 * t16 * t19 + t11 * t5
      t13 = -2 * t18 * t16 * t19 + t13 + t14
      t24 = t9 * t13
      t25 = t7 * (t17 * t23 + t24)
      t24 = t7 * (t21 * t6 - t24)
      t26 = 2 * t7
      t27 = -t26 * t17 * t16 * t19 + t12 * t4
      t28 = t7 * t13
      t29 = t9 * (t27 * t8 + t28)
      t26 = t26 * t16 * t6 * t19 - t12 * t3
      t28 = t9 * (t1 * t26 - t28)
      t25 = 0.1D1 / t25
      t30 = 0.1D1 / t6
      t29 = 0.1D1 / t29
      t31 = 0.1D1 / t5
      t32 = 0.1D1 / t9
      t33 = 0.1D1 / t17
      t27 = 0.1D1 / t27
      t3 = 0.1D1 / t3
      t26 = 0.1D1 / t26
      t24 = 0.1D1 / t24
      t34 = 0.1D1 / t8
      t23 = 0.1D1 / t23
      t35 = 0.1D1 / t11
      t4 = 0.1D1 / t4
      t28 = 0.1D1 / t28
      t29 = t1 * t29
      t36 = t26 * t32 + t29
      t37 = t21 ** 2
      t38 = t25 * t3
      t39 = t9 * t35
      t40 = t30 * t24
      t28 = t26 * t28
      t41 = t1 * t13 ** 2
      t42 = t4 * t31
      t43 = t17 * t4
      t44 = t43 + t39
      t45 = t25 ** 2
      t46 = t9 ** 2
      t47 = t30 ** 2
      t48 = t7 ** 2
      t6 = t1 * (-t20 * t16 * t6 * t19 + t14 + t15)
      t14 = t7 * t21
      t15 = t29 * t9
      t20 = t21 * t23
      t49 = t31 * t19
      t50 = t20 * t25
      t51 = t18 * t10
      t52 = t42 * t19 * t1
      t53 = t52 * ((t3 * (t13 * (-t15 - t26) + t30 * t9) - t28 * t41 * t
     #39) * (t1 * t12 + t2 * t7) + t13 * (t17 * t36 * t3 + t1 * (t28 * t
     #17 * t13 + t29 * (t17 * t13 * t27 + t9) * t34) * t35) * t10 + t51 
     #* (t39 * (t50 * t33 + t40) - t38) * t22)
      t16 = 0.1D1 / t16
      t54 = t14 * t25
      t55 = t52 * t17 * t10 * t3 * (t54 + t30)
      t17 = t52 * t51 * t21 * t35 * (t40 * t17 + t50)
      ret = -8 * t42 * t2 * t1 * (t3 * (-t13 * t36 + t30) + t7 * (t21 * 
     #(-t40 * t39 + t38) - t39 * t25 * t33 * t23 * t37) - t41 * (t29 * t
     #27 * t34 + t28) * t35) - 32 * t49 * (t10 * ((t46 * t45 * t23 * t35
     # * t4 * t21 * t37 - t9 * t3 * t44 * t45 * t37) * t7 * t48 + t6 * t
     #43 * t7 * t46 * t47 * t24 * t35 + t43 * t48 * t46 * t21 * t35 * (-
     #t14 + t6) * t24 ** 2 * t30 + t1 * t13 * t3 * t35 * (t15 + t26)) + 
     #t1 * t9 * t47 * t3 * t44 * (t12 * t8 + t5 * t7) + t18 * t25 * t21 
     #* t10 * (t39 * t21 * t33 * t23 ** 2 * t4 + (-t43 * t3 + t39 * (t20
     # * t4 - t3)) * t25 * t7) * t22 * t8) - 16 * t53 - 128 * t49 * t1 *
     # t2 * t10 * t3 * t16 * (t43 * t11 * t32 + 1) + 64 * t3 * t31 * t1 
     #* (t11 * t2 ** 2 * t4 * t32 * t16 + t39 * t10 * t19 * (t54 + t30))
     # + 48 * t55 - 80 * t17

      hjetmass_triangle_ppmm_s124_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
