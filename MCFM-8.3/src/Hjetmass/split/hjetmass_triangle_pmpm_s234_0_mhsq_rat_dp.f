
      double complex function hjetmass_triangle_pmpm_s234_0_mhsq_rat_dp 
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
      t3 = za(i2, i3)
      t4 = zb(i4, i1)
      t5 = t1 * t4 + t2 * t3
      t6 = za(i1, i2)
      t7 = zb(i2, i1)
      t8 = za(i1, i3)
      t9 = za(i1, i4)
      t10 = t6 * t7
      t11 = t8 * t2
      t12 = t9 * t4
      t13 = t11 + t12 + t10
      if ( dreal(t13) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t11 = cg * cdsqrt(t13 ** 2) + t10 + t11 + t12
      t12 = za(i3, i4)
      t13 = zb(i4, i3)
      t14 = zb(i3, i2)
      t15 = zb(i4, i2)
      t16 = t3 * t14
      t17 = t1 * t15
      t18 = t12 * t13
      t19 = t17 + t18 + t16
      t11 = 0.1D1 / t11
      t10 = -2 * t10 * t19 * t11 + t16 + t17
      t17 = -2 * t8
      t15 = t2 * (t17 * t7 * t19 * t11 + t12 * t15)
      t20 = t6 * (t10 * t7 + t15)
      t21 = -2 * t9
      t22 = t21 * t7 * t19 * t11 - t12 * t14
      t23 = -2 * t6
      t24 = t23 * t2 * t19 * t11 + t1 * t13
      t8 = t8 * t24
      t25 = t7 * (t10 * t6 + t8)
      t23 = t23 * t4 * t19 * t11 - t13 * t3
      t8 = t7 * (t23 * t9 + t8)
      t21 = t21 * t2 * t19 * t11 + t1 * t14
      t15 = t6 * (t22 * t4 + t15)
      t14 = 0.1D1 / t14
      t26 = 0.1D1 / t9
      t27 = 0.1D1 / t19
      t13 = 0.1D1 / t13
      t28 = 0.1D1 / t6
      t29 = 0.1D1 / t10
      t3 = 0.1D1 / t3
      t12 = 0.1D1 / t12
      t30 = t1 ** 2
      t31 = t26 * t13
      t32 = t28 * t14
      t33 = t12 * t3
      t27 = t33 * t27
      t34 = 0.1D1 / t4
      t15 = 0.1D1 / t15
      t23 = 0.1D1 / t23
      t25 = 0.1D1 / t25
      t20 = 0.1D1 / t20
      t35 = 0.1D1 / t7
      t8 = 0.1D1 / t8
      t36 = t24 ** 2
      t37 = t24 * t36
      t38 = t2 * t10
      t39 = t2 * t22
      t40 = t1 * t3
      t41 = t12 * t20
      t42 = t6 * t26
      t43 = t9 * t25
      t33 = t33 * t1
      t44 = t13 * t14
      t45 = t44 * t2
      t46 = t9 ** 2
      t47 = t25 ** 2
      t48 = t8 ** 2
      t49 = t6 ** 2
      t50 = t7 ** 2
      t4 = t4 * t21
      t16 = t2 * (t17 * t2 * t19 * t11 + t16 + t18)
      t17 = t33 * t44 * t9 * t6
      t18 = t45 * t11
      t19 = t18 * t7 * t5 * (t43 * t12 * t24 * (-t29 * (t9 * t24 * t28 +
     # t21) - t40) + t3 * (t1 * t21 * t12 + t23 * t36) * t8 * t6 + t49 *
     # t24 * t21 * t26 * t3 * t23 * t8)
      t51 = t2 ** 2
      t18 = t18 * t5 * (t3 * (-t51 * t6 * t34 * t15 * (t42 * t22 + t10) 
     #+ (t6 * (t21 * t7 * t25 + t39 * (-t15 + t20)) + t9 * (-t7 * t24 * 
     #t8 + t38 * t20) - t2 * t49 * t22 ** 2 * t26 * t15 * t29) * t12 * t
     #1) + t41 * t51 * t22 * (t6 * t22 * t29 + t9) * t35)
      t20 = -32
      ret = 64 * t27 * t11 * t5 * t2 * t30 * (-t32 * t9 + t13 + (t31 * t
     #6 - t14) * t29 * t22) - 8 * t45 * t1 * (t2 * (t42 * t3 * (t1 * t22
     # * t12 + t38 * t34) * t15 + t41 * (t40 * t10 + t39 * t35)) + t7 * 
     #(t36 * (t43 * t28 * t29 * t12 + t42 * t23 * t3 * t8) + t33 * (t25 
     #+ t8) * t24)) + t20 * (t5 * ((t44 * t46 * t47 * t11 * t29 * t12 * 
     #t37 + t17 * t47 * t11 * t36) * t7 * t50 + t50 * (t44 * (t16 * t46 
     #* t47 * t29 * t12 + (-t16 - t4) * t48 * t3 * t23 * t49) * t11 * t3
     #6 + t32 * t46 * t25 * t11 * t13 * t29 ** 2 * t12 * t37 + t17 * (t1
     #6 * (-t48 + t47) - t4 * t48) * t11 * t24) - t4 * t31 * t49 * t36 *
     # t14 * t23 ** 2 * t11 * t3 * t8 * t7) - t27 * t1 * t30 * t2 * (t32
     # + t31)) - 48 * t19 + 16 * t18

      hjetmass_triangle_pmpm_s234_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
