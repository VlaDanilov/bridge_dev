
      double complex function hjetmass_triangle_pmpm_s12_s123_0_dp 
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
      t2 = za(i1, i3)
      t3 = zb(i3, i1)
      t4 = za(i2, i3)
      t5 = zb(i3, i2)
      t6 = t2 * t3
      t7 = t4 * t5
      t8 = t7 + t6
      if ( dreal(t8) > 0d0) then; cg = 1d0; else; cg = -1d0; end if
      t8 = cg * cdsqrt(t8 ** 2) + t6 + t7
      t9 = zb(i2, i1)
      t10 = 0.1D1 / t8
      t11 = 2 * t6 * t10
      t12 = t1 * t9
      t13 = t12 * (1 - t11)
      t14 = zb(i4, i2)
      t15 = zb(i4, i3)
      t10 = t1 * (2 * t2 * t9 * t15 * t10 - t14)
      t16 = za(i3, i4)
      t11 = t3 * (-t10 * t16 + t13 * t2) - t7 * t12 * t11
      t12 = za(i1, i4)
      t17 = za(i2, i4)
      t18 = t12 * t3
      t19 = t17 * t5
      t20 = zb(i4, i1)
      t9 = 0.1D1 / t9
      t21 = 0.1D1 / t15
      t22 = 0.1D1 / t8
      t11 = 0.1D1 / t11
      t23 = 0.1D1 / t16
      t10 = 0.1D1 / t10
      t2 = 0.1D1 / t2
      t24 = 0.1D1 / t1
      t25 = t12 * t20
      t26 = t14 * t17
      t27 = t26 + t25
      t28 = t9 ** 2
      t29 = t2 ** 2
      t30 = t2 * t29
      t31 = mt ** 2
      t32 = t3 ** 2
      t33 = t3 * t32
      t34 = t11 ** 2
      t35 = t11 * t34
      t36 = t4 ** 2
      t37 = t16 ** 2
      t38 = t22 ** 2
      t39 = t8 ** 2
      t40 = t17 ** 2
      t41 = t7 * t21
      t42 = t1 * t16
      t43 = t41 * t13
      t44 = t9 * t13
      t45 = t44 * t2
      t46 = t7 * t17
      t47 = t4 * t12
      t48 = t3 * t16
      t49 = t48 * t7
      t50 = t21 * t34
      t51 = t2 * t34
      t52 = t13 * t8
      t53 = t13 ** 2
      t54 = t10 ** 2
      t55 = t1 ** 2
      t56 = t22 * t8
      t57 = t56 * t2
      t58 = t10 * t15
      t59 = t32 * t21
      t60 = t59 * t2
      t61 = t2 * t54
      t12 = t8 * (t33 * (t41 * t55 * t53 * t16 * t2 * t35 + t22 * t13 * 
     #t4 * (t45 * t16 * t12 + t19 * (t44 + t1) * t21) * t34) + t60 * (t3
     #1 * (t16 * t9 * t28 * t2 * t24 * t53 + t17 * t28 * t24 * t13) - t4
     #6 * t1 * t22) * t11 + t61 * t7 * t53 * t28 * (t58 + t57)) + t52 * 
     #(t33 * (t45 * t22 * (t36 * t5 * t12 * t21 - t1 * t37 - t42 * t41) 
     #* t34 + t43 * t16 * t28 * t2 * (t12 ** 2 * t20 ** 2 + t14 ** 2 * t
     #40) * t35) + t49 * t28 * t38 * t39 * t13 * t30 * t21 * t11 + t28 *
     # t3 * t21 * t10 * t2 * (t17 * t23 + t45) * t24 * t31 + t43 * t28 *
     # t38 * t39 * t30 * t10 + t50 * t44 * t22 * (t47 - t42) * t32 ** 2 
     #+ t51 * t7 * t22 * (t13 * (t9 * (t16 * t17 + t21 * (t42 * t8 * t2 
     #+ t46)) + t8 * t16 * t21 * t2 * t27 * t28) + t1 * t17 * (t41 + t16
     #)) * t32)
      t14 = t27 * t9
      t20 = t14 + t1
      t27 = -t14 - t1
      t28 = t1 * t22
      t30 = t2 * t10
      t42 = t48 * t11 + t10
      t41 = t51 * t41 * t16
      t62 = t3 * t21
      t63 = t59 * t11
      t38 = t7 * t38
      t31 = t45 * (t59 * t9 * t42 * t31 + t63 * (t22 * (t4 * (-t19 - t18
     #) + t48 * t1) + t49 * (t25 * t1 + t26 * (t25 * t9 + t1)) * t34 * t
     #13) * t8 + t38 * t13 * t1 * (t32 * (-t51 * t37 - t41) + t61 - t50 
     #* t16 * t33) * t39 + t9 * t4 * t31 * (t23 * (-t62 * t10 - t13 * t5
     #4) + t32 * t11 * (t13 * t16 * t11 - t21)) * t24 * (t19 + t18))
      t34 = t44 * t24 + 1
      t49 = t11 * t3
      t23 = t44 * t8 * (t60 * t9 * t42 + t60 * t11 * (-t47 * t23 + t44) 
     #* t24 * t17 + (t49 * t2 + t30 * t23 + t63 * t23) * t24 * t40)
      t38 = t38 * t52 * t3 * t55 * (t32 * (t37 * (t7 * t13 * t2 * t35 - 
     #t51) - t41) + t33 * (t16 * (t43 * t35 - t50) + t13 * t35 * t37) + 
     #t61)
      ret = -8 * t12 + 32 * t56 * t7 * t1 * (t16 * (-t60 * t22 * (t1 * t
     #3 + t45 * t8) * t11 + t21 * t53 * t33 * (t3 * (t1 * (-t6 * t22 + 1
     #) + t14) - t1 * t36 * t5 ** 2 * t22 * t2 + t7 * t20 * t2) * t35 + 
     #t50 * t2 * t13 * t33 * t27) - t28 * t16 * t37 * t33 * t53 * t15 * 
     #t2 * t35 + t30 * (-t28 * (t53 * t15 * t54 + t59) + t44 * (t10 * (-
     #t58 * t13 - t3) - t57 * t3 * t21)) + t2 * t53 * t33 * t20 * t35 * 
     #t37) - 16 * t31 + 4 * t62 * t11 * t13 * (t17 * t32 * t9 + t22 * t9
     # * (t44 * ((t19 + t18) * t24 * t4 - t48) + t46) * t29 * t39 + t49 
     #* (t4 * (t19 * (t9 * (t25 * t34 + t26 * t34 + t13) + t1) + t18 * t
     #44 * (t14 * t24 + 1)) + t44 * t48 * t27) * t2 * t8) - 64 * t38 + 2
     # * t23

      hjetmass_triangle_pmpm_s12_s123_0_dp = ret/32d0/(0,1d0)
      return

      end function
