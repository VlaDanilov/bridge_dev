
      t1 = 0.1D1 / sman
      t2 = mH ** 2
      t3 = dilogc(-t1 * t2 + 1)
      t4 = dilogc(-t2 / tman + 1)
      t5 = 0.1D1 / uman
      t2 = dilogc(-t2 * t5 + 1)
      t6 = cdlogwrap(uman * t1)
      t5 = cdlogwrap(tman * t5)
      t1 = cdlogwrap(tman * t1)
      t7 = uman ** 2
      t8 = sman ** 2
      t9 = tman ** 2
      t10 = t9 + t7 + t8
      t10 = LogMums * t10 + LogMumt * t10 + LogMumu * t10
      t1 = t1 ** 2
      t6 = t6 ** 2
      t11 = pi ** 2
      t5 = t5 ** 2
      t12 = -t5 - t1 - t6 + t11
      t13 = sman * tman
      t14 = (sman + tman) * uman + t13
      t15 = (-sman - tman) * uman - t13
      t16 = (sman + tman) * uman + t13
      t17 = t3 + t4 + t2
      t18 = t9 + t7 + t8
      t19 = LogMums * t14 + LogMumt * t14 + LogMumu * t14
      t20 = -2
      t21 = 4
      t22 = (0.2D1 / 0.3D1)
      t23 = (0.4D1 / 0.3D1)
      t24 = (0.8D1 / 0.3D1)
      c1mt0 = (0.11D2 / 0.9D1) * cA * t10 - cA * (t12 * t7 + t12 * t8 + 
     #t12 * t9) / 3 + (0.128D3 / 0.9D1) * cA * t14 + t21 * (Kg * t16 + c
     #F * t15) + t23 * ((t17 * t7 + t17 * t8 + t17 * t9) * cA + ((-t8 - 
     #t9) * tr + (-uman * tr + sman + tman) * uman + t13) * LogMuMtop) -
     # (0.20D2 / 0.9D1) * tr * t10 - (0.140D3 / 0.9D1) * tr * t14 + t24 
     #* (cA * (t16 * t2 + t16 * t3 + t16 * t4) + t15 * LogMuMtop * tr) +
     # t22 * ((t1 * t16 + t11 * t15 + t16 * t5 + t16 * t6) * cA + t18 * 
     #LogMuMtop) + t20 * (-Kg * t18 + cF * t18) + 7 * t18 * cA - (0.20D2
     # / 0.3D1) * t18 * tr + (0.22D2 / 0.9D1) * cA * t19 - (0.40D2 / 0.9
     #D1) * tr * t19
