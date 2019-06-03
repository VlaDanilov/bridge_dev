
      complex*32 function hjetmass_triangle_ppmm_s23_s234_0
     &     (i1,i2,i3,i4,za,zb,mt)
          implicit complex*32 (t)
          integer i1,i2,i3,i4
          include 'types.f'
          include 'mxpart.f'
          include 'constants.f'
          complex*32 za(mxpart,mxpart), zb(mxpart,mxpart)
          
          complex*32 ret
          double precision mt
          real*16 cg

      t1 = za(i2, i3)
      t2 = zb(i3, i2)
      t3 = za(i2, i4)
      t4 = zb(i4, i2)
      t5 = za(i3, i4)
      t6 = zb(i4, i3)
      t7 = t5 * t6
      t8 = t3 * t4
      t9 = t7 + t8
      if ( real(t9) > 0q0) then; cg = 1q0; else; cg = -1q0; end if
      t9 = cg * sqrt(t9 ** 2) + t7 + t8
      t10 = za(i1, i4)
      t11 = zb(i2, i1)
      t12 = zb(i3, i1)
      t13 = t3 * t11
      t14 = t5 * t12
      t15 = t14 + t13
      t16 = zb(i4, i1)
      t17 = 0.1q1 / t9
      t18 = t2 * t16
      t19 = -2 * t18
      t20 = t1 * (t19 * t5 * t17 - t11)
      t21 = za(i1, i2)
      t22 = za(i1, i3)
      t23 = t21 * t11
      t24 = t22 * t12
      t25 = t19 * t10 * t1 * t17 + t23 + t24
      t19 = t1 * (t19 * t3 * t17 + t12)
      t26 = (t16 * t25 + t19 * t4 + t20 * t6) * t10
      t25 = 0.1q1 / t25
      t9 = 0.1q1 / t9
      t6 = 0.1q1 / t6
      t26 = 0.1q1 / t26
      t27 = t1 * t2
      t28 = t27 + t23 + t24
      t29 = t8 * t6
      t30 = t29 + t5
      t31 = t15 ** 2
      t32 = t25 ** 2
      t33 = t10 ** 2
      t34 = t4 ** 2
      t35 = t26 ** 2
      t36 = t26 * t35
      t37 = t16 ** 2
      t38 = t27 * t5
      t39 = t20 * t33
      t40 = t20 * t6
      t41 = t27 * t9
      t42 = 0.1q1 / t1
      t43 = 0.1q1 / t10
      t44 = 0.1q1 / t2
      t45 = 0.1q1 / t16
      t46 = t16 * t26
      t47 = t25 * t43
      t48 = t47 - t46
      t49 = mt ** 2
      t50 = t15 * t4
      t51 = t9 * t5
      t52 = t50 * t6
      t53 = t16 * t10
      t54 = t1 * t15
      t55 = t20 * (t14 + t13) - t54 * t11
      t56 = t11 ** 2
      t57 = t3 * t20
      t58 = t4 * t12
      t12 = t50 * (t26 * (t33 * (t50 * t40 * t18 * t1 * t35 + t46 * t6 *
     # t9 * (t4 * (-t54 + t57) * t11 + t20 * (t58 + t18) * t5)) + t51 * 
     #t6 * (t5 * (-t58 - t18) - t8 * t11) + t9 * (t4 * (t29 * t55 + t5 *
     # t55) + t18 * t20 * t5 * t30) * t26 * t10) + t6 * (t4 * (-t15 * t2
     #0 * t25 * t32 * t10 + t39 * t16 * t36 * (t12 ** 2 * t22 ** 2 + t21
     # ** 2 * t56) * t15) + t49 * t5 * t11 * (-t47 * t45 + t26)) * t44 *
     # t42)
      t21 = t24 + t23
      t2 = t52 * t35 * t10 * (t4 * (t20 * ((t13 * t21 + t14 * t21) * t44
     # * t42 + t13 + t14) + t11 * (-t21 * t44 - t1) * t15) + t20 * t5 * 
     #(t21 * t42 + t2) * t16)
      t4 = t4 * (t26 * (t11 * t15 * (-t10 * t11 * t6 + (t58 * t6 - t11) 
     #* t45 * t5) * t44 + t5 * t6 * (t11 * (t19 * t5 - t57) - t5 * t15 *
     # t16) * t42) + t15 * t6 * (t5 ** 2 * t42 * t43 + t56 * t44 * t45) 
     #* t25)
      ret = -64 * t9 * t31 * t34 * (t16 * (t38 * t10 * t9 * t30 * t35 + 
     #t39 * (t5 * (t27 * (-t7 * t9 + 1) + t23 + t24) - t27 * t3 ** 2 * t
     #34 * t9 * t6 + t29 * t28) * t36) + t37 * (t40 * t10 * t33 * t28 * 
     #t36 + t38 * t33 * t9 * t6 * t35) + t6 * t32 * (-t38 * t9 + (t41 + 
     #1) * t25 * t20 * t10) - t41 * t40 * t33 ** 2 * t16 * t37 * t36) + 
     #32 * t52 * (t42 * (t5 * (t29 * t48 + t48 * t5) + t50 * (t26 * (t46
     # * t39 - t5) - t20 * t32 * t45 + t47 * t5 * t45) * t44) * t17 * t4
     #9 + t50 * (t51 * t28 * t35 * t16 * t10 + t51 * t32 + t39 * (t23 * 
     #(t24 * t42 * t44 + 1) + t24) * t36 * t16)) + 128 * t39 * t18 * t36
     # * t9 ** 2 * t31 * t34 * t1 * (t53 * t5 + t8 * (t53 * t6 + t5)) + 
     #16 * t12 - 8 * t2 - 4 * t4

      hjetmass_triangle_ppmm_s23_s234_0 = ret/32q0/(0,1q0)
      return

      end function
