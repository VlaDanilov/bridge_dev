
      double complex function hjetmass_triangle_pmpm_s134_0_mhsq_rat_dp 
     &     (i1,i2,i3,i4,za,zb)
          implicit double complex (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          include 'zprods_decl.f'
          double complex ret
          double precision cg


      t1 = za(i1, i2)
      t2 = za(i2, i4)
      t3 = zb(i2, i1)
      t4 = zb(i3, i1)
      t5 = za(i2, i3)
      t6 = zb(i4, i1)
      t7 = t2 * t6 + t4 * t5
      t8 = za(i1, i4)
      t9 = zb(i3, i2)
      t10 = za(i1, i3)
      t11 = za(i3, i4)
      t12 = zb(i4, i3)
      t13 = t10 * t4
      t14 = t8 * t6
      t15 = t11 * t12
      t16 = t13 + t14 + t15
      t17 = zb(i4, i2)
      t18 = t1 * t3
      t19 = t5 * t9
      t20 = t2 * t17
      t21 = t18 + t19 + t20
      if ( dreal(t21) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t19 = cg * cdsqrt(t21 ** 2) + t18 + t19 + t20
      t19 = 0.1D1 / t19
      t20 = 2 * t2
      t21 = -t20 * t9 * t16 * t19 + t4 * t8
      t13 = -2 * t18 * t16 * t19 + t13 + t14
      t22 = t20 * t3 * t16 * t19 - t11 * t4
      t23 = t17 * t22
      t24 = t1 * (t13 * t3 - t23)
      t25 = 2 * t5 * t3 * t16 * t19 + t11 * t6
      t23 = t1 * (t25 * t9 + t23)
      t26 = t1 * t4 + t2 * t12
      t27 = 2 * t1
      t28 = t27 * t9 * t16 * t19 + t12 * t8
      t10 = t2 * (t27 * t17 * t16 * t19 - t10 * t12)
      t27 = t3 * (t1 * t13 - t10)
      t10 = t3 * (t28 * t5 + t10)
      t24 = 0.1D1 / t24
      t29 = 0.1D1 / t5
      t8 = 0.1D1 / t8
      t30 = 0.1D1 / t28
      t12 = 0.1D1 / t12
      t27 = 0.1D1 / t27
      t6 = 0.1D1 / t6
      t31 = 0.1D1 / t9
      t10 = 0.1D1 / t10
      t11 = 0.1D1 / t11
      t23 = 0.1D1 / t23
      t32 = 0.1D1 / t13
      t33 = 0.1D1 / t1
      t34 = 0.1D1 / t3
      t25 = 0.1D1 / t25
      t35 = t9 ** 2
      t36 = t22 ** 2
      t37 = t22 * t36
      t38 = t3 ** 2
      t39 = t28 * t27
      t40 = t4 * t6
      t41 = t39 * t31
      t42 = t2 * t7
      t43 = t4 * t12
      t44 = t43 * t7
      t45 = t3 * t22
      t46 = t11 * t8
      t47 = t46 * t19
      t48 = t47 * t2
      t49 = 0.1D1 / t16
      t50 = t23 ** 2
      t51 = t24 ** 2
      t52 = t1 ** 2
      t53 = t4 ** 2
      t5 = t5 * t21
      t14 = t2 * (-t20 * t17 * t16 * t19 + t14 + t15)
      t15 = t47 * t35
      t16 = t47 * t43
      t17 = t11 * t31
      t20 = t8 * t34
      t5 = t5 * t47 * t1 * t38 * t7 * t36 * t31 * t6 * t25 ** 2 * t23 + 
     #(t16 * t3 * t9 * t51 * t6 * t36 - t15 * t51 * t12 * t32 * t37) * t
     #7 * t1 * t52 + (-t15 * t24 * t12 * t32 ** 2 * t34 * t37 + t16 * t4
     #5 * t6 * t9 * (t14 * (t50 - t51) + t5 * t50) + t47 * (t14 * t35 * 
     #t51 * t12 * t32 + (t5 + t14) * t25 * t50 * t6 * t38) * t36) * t7 *
     # t52 - t2 * t4 * t53 * t6 * t12 * t49 * (t17 + t20)
      t4 = t46 * t4 * t2 * (t1 * (t36 * (-t9 * t24 * t12 * t32 * t34 - t
     #3 * t6 * t31 * t23 * t25) + t43 * t6 * (t24 - t23) * t22) + t2 * (
     #t3 * (-t2 * t12 * t27 * t31 * t32 * t33 * t28 ** 2 + t41 * t43 * t
     #6) + t13 * t6 * t10 * (-t2 * t13 * t29 * t30 + t43)))
      ret = 16 * t48 * (t1 * (t7 * t35 * t24 * t12 * t32 * t34 * t36 - t
     #44 * t9 * t24 * t6 * t22) + t2 * (t12 * (t40 * t7 * t9 * t13 * t10
     # + (t39 * t2 * t33 + t40 * (-t27 - t10) * t13) * t26 * t3 + t41 * 
     #(-t2 * t28 * t32 * t33 + t40) * t26 * t38) + t42 * t13 * t6 * t10 
     #* (-t9 * t13 * t30 + t3) * t29) + t18 * t7 * t6 * (t45 * t31 * t25
     # + t43) * t23 * t21) - 32 * t5 + 48 * t48 * t44 * t1 * t6 * (t3 * 
     #t21 * t24 - t9 * t22 * t23) - 80 * t47 * t42 * t1 * t22 * (t9 * t2
     #1 * t12 * t32 * t24 + t45 * t6 * t25 * t23) + 64 * t6 * t12 * t19 
     #* t49 * t53 * t2 * (t26 * (-t17 * t3 + t8) + t7 * (-t20 * t9 + t11
     #)) - 8 * t4

      hjetmass_triangle_pmpm_s134_0_mhsq_rat_dp = ret/32d0/(0,1d0)
      return

      end function
